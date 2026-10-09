#!/usr/bin/env python3
"""Build a cells × genes UMI counts h5ad from featureCounts-annotated BAMs.

Walks each `*.bam.featureCounts.bam` in `--anno-dir`, filters to primary-mate
Assigned reads, dedupes per (cell, gene, UMI), and writes a sparse AnnData.

The UMI is parsed from the read name suffix (`..._UMI_<umi>`) that
encode_umi_in_readname.py placed there before alignment.

Writes:
    X                    cells × genes UMI counts (int32 CSR)
    layers['read_counts'] cells × genes pre-dedup read counts
    obs                  per-cell QC: n_reads_mapped, n_reads_assigned, n_umis,
                         n_genes_detected, reads_per_umi
    var                  gene_id, Chr, Start, End, Strand, Length (from counts_raw.txt)
    uns['pipeline_params']  record of relaxed STAR / featureCounts params
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pysam
from scipy.sparse import csr_matrix


UMI_PAT = re.compile(r'_UMI_([ACGTN]+)$')


# -------- Gene intervaltree built from counts_raw.txt (for ambiguity resolution) ----------

def _build_gene_arrays(counts_raw: str, extension: int) -> dict:
    """chrom -> (starts, ends, ids) sorted by start; coordinates are 1-based
    gene body spans with `extension` bp added to both sides. Each gene's
    overall body is min(Start) .. max(End) from the semicolon-joined
    transcript positions in counts_raw.txt."""
    df = pd.read_csv(counts_raw, sep='\t', comment='#')
    df = df.rename(columns={'Geneid': 'gene_id'})
    rows_per_chrom: dict = {}
    for _, row in df.iterrows():
        starts = [int(s) for s in str(row['Start']).split(';')]
        ends   = [int(s) for s in str(row['End']).split(';')]
        chrom  = str(row['Chr']).split(';')[0]
        gs = max(1, min(starts) - extension)
        ge = max(ends) + 1 + extension
        if ge > gs:
            rows_per_chrom.setdefault(chrom, []).append((gs, ge, row['gene_id']))
    out: dict = {}
    for chrom, rows in rows_per_chrom.items():
        rows.sort(key=lambda r: r[0])
        s = np.array([r[0] for r in rows], dtype=np.int64)
        e = np.array([r[1] for r in rows], dtype=np.int64)
        ids = np.array([r[2] for r in rows], dtype=object)
        out[chrom] = (s, e, ids)
    return out


def _resolve_ambiguous(gene_arrays: dict,
                       chrom: str,
                       pos_lo: int, pos_hi: int) -> str | None:
    """Pick the gene with the longest overlap with [pos_lo, pos_hi); tiebreak
    alphabetical gene_id. Returns None if no overlap (shouldn't happen for
    Ambiguity reads but is a safe fallback)."""
    arr = gene_arrays.get(chrom)
    if arr is None:
        return None
    starts, ends, ids = arr
    hi = np.searchsorted(starts, pos_hi, side='right')
    if hi == 0:
        return None
    cand = np.arange(hi)
    mask = ends[cand] > pos_lo
    cand = cand[mask]
    if cand.size == 0:
        return None
    # Compute overlap for each candidate
    best_overlap = -1
    best_gene = None
    for i in cand:
        ovl = min(pos_hi, int(ends[i])) - max(pos_lo, int(starts[i]))
        gene = str(ids[i])
        if ovl > best_overlap or (ovl == best_overlap and best_gene is not None and gene < best_gene):
            best_overlap = ovl
            best_gene = gene
    return best_gene


def _barcode_from_fcbam(p: Path) -> str:
    name = p.name
    if name.endswith('.bam.featureCounts.bam'):
        return name[:-len('.bam.featureCounts.bam')]
    return name.removesuffix('.bam').removesuffix('.featureCounts')


def _parse_umi(query_name: str) -> str | None:
    m = UMI_PAT.search(query_name)
    return m.group(1) if m else None


def collect_umis(anno_dir: Path, allow_no_umi: bool,
                 gene_arrays: dict | None = None) -> tuple[dict, pd.DataFrame]:
    """umi_sets[barcode][gene_id] = set of UMI tokens; also per-cell QC.

    When allow_no_umi is True, reads without a `_UMI_<umi>` suffix are kept
    using a coordinate-based surrogate token `POS_<r1_start>_<r2_start>` so
    exact-position PCR duplicates still collapse to one count.

    When gene_arrays is provided, reads with XS=Unassigned_Ambiguity are
    rescued by picking the gene with the longest overlap with R1 (tiebreak
    alphabetical gene_id). This matches the featureCounts --largestOverlap
    behavior but also resolves the ties that fc leaves as Ambiguity.
    """
    bams = sorted(anno_dir.glob('*.bam.featureCounts.bam'))
    if not bams:
        raise RuntimeError(f'No annotated BAMs under {anno_dir}')
    umi_sets: dict[str, dict[str, set]] = {}
    qc_rows = []
    for i, bam_path in enumerate(bams):
        bc = _barcode_from_fcbam(bam_path)
        if bc in umi_sets:
            continue
        sets: dict[str, set] = defaultdict(set)
        n_reads_mapped = 0
        n_pairs = 0
        n_assigned_pairs = 0
        n_umi_tagged_pairs = 0
        n_coord_surrogate = 0
        n_ambig_resolved = 0
        try:
            bam = pysam.AlignmentFile(str(bam_path), 'rb')
        except Exception as e:
            print(f'[skip] {bc}: open failed: {e}', file=sys.stderr)
            continue
        try:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped:
                    continue
                n_reads_mapped += 1
                # Count only primary R1 per pair to avoid double-counting
                if read.is_secondary or read.is_supplementary:
                    continue
                if not read.is_read1:
                    continue
                n_pairs += 1
                try:
                    xs = read.get_tag('XS')
                except KeyError:
                    continue
                gene = None
                if xs == 'Assigned':
                    try:
                        gene = read.get_tag('XT')
                    except KeyError:
                        continue
                elif xs == 'Unassigned_Ambiguity' and gene_arrays is not None:
                    # Resolve the tie ourselves using the same gene bodies +
                    # extensions that featureCounts was using.
                    gene = _resolve_ambiguous(
                        gene_arrays,
                        read.reference_name,
                        read.reference_start + 1,     # 1-based inclusive
                        read.reference_end,           # 1-based inclusive -> half-open
                    )
                    if gene is None:
                        continue
                    n_ambig_resolved += 1
                else:
                    continue
                umi = _parse_umi(read.query_name)
                if umi is not None:
                    n_umi_tagged_pairs += 1
                elif allow_no_umi:
                    # coordinate-based surrogate: exact-position PCR dups collapse
                    umi = f'POS_{read.reference_start}_{read.next_reference_start}'
                    n_coord_surrogate += 1
                else:
                    continue
                n_assigned_pairs += 1
                sets[gene].add(umi)
        finally:
            bam.close()
        umi_sets[bc] = sets
        n_umis_total = sum(len(s) for s in sets.values())
        qc_rows.append({
            'barcode': bc,
            'n_reads_mapped': n_reads_mapped,
            'n_pairs': n_pairs,
            'n_pairs_assigned': n_assigned_pairs,
            'n_pairs_with_umi': n_umi_tagged_pairs,
            'n_pairs_coord_surrogate': n_coord_surrogate,
            'n_pairs_ambig_resolved': n_ambig_resolved,
            'n_umis': n_umis_total,
            'n_genes_detected': len(sets),
            'reads_per_umi': (n_assigned_pairs / n_umis_total) if n_umis_total else np.nan,
        })
        if (i + 1) % 50 == 0:
            print(f'  [{i+1}/{len(bams)}] {bc}: pairs={n_pairs} assigned={n_assigned_pairs} '
                  f'umis={n_umis_total}', flush=True)
    return umi_sets, pd.DataFrame(qc_rows)


def build_matrix(umi_sets: dict[str, dict[str, set]],
                 gene_order: list[str]) -> tuple[csr_matrix, list[str]]:
    barcodes = sorted(umi_sets.keys())
    gene_ix = {g: i for i, g in enumerate(gene_order)}
    rows, cols, data = [], [], []
    for ri, bc in enumerate(barcodes):
        for gene, umis in umi_sets[bc].items():
            gi = gene_ix.get(gene)
            if gi is None: continue
            rows.append(ri)
            cols.append(gi)
            data.append(len(umis))
    X = csr_matrix(
        (np.asarray(data, dtype=np.int32),
         (np.asarray(rows, dtype=np.int32),
          np.asarray(cols, dtype=np.int32))),
        shape=(len(barcodes), len(gene_order)),
    )
    return X, barcodes


def load_gene_meta(counts_raw: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(counts_raw, sep='\t', comment='#')
    bam_cols = df.columns[6:]
    barcodes = [os.path.basename(c).removesuffix('.bam') for c in bam_cols]
    gene_meta = df[['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']].rename(
        columns={'Geneid': 'gene_id'}).copy()
    counts_df = df[bam_cols].copy()
    counts_df.columns = barcodes
    counts_df.index = gene_meta['gene_id'].values
    return gene_meta, counts_df


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--anno-dir', required=True,
                    help='Directory of featureCounts annotated BAMs')
    ap.add_argument('--counts-raw', required=True,
                    help='featureCounts counts_raw.txt (gene meta + read_counts layer)')
    ap.add_argument('--out-h5ad', required=True)
    ap.add_argument('--run-min', type=int, default=15,
                    help='--run-min used by the primer-trim extraction (for uns)')
    ap.add_argument('--allow-no-umi', action='store_true',
                    help='Keep reads without a `_UMI_<umi>` suffix; use '
                         'coordinate-based surrogate (POS_<r1>_<r2>) for dedup. '
                         'Default: skip no-UMI reads (strict primer-trimmed mode).')
    ap.add_argument('--resolve-ambiguous', action='store_true',
                    help='Rescue Unassigned_Ambiguity reads by picking the '
                         'gene with the longest R1 overlap (tiebreak = '
                         'alphabetical gene_id).')
    ap.add_argument('--gene-extension', type=int, default=500,
                    help='Bp to extend each gene body on both sides when '
                         'resolving ambiguity (match featureCounts '
                         '--readExtension5/3). Default 500.')
    args = ap.parse_args()

    anno_dir = Path(args.anno_dir)
    gene_meta, read_counts = load_gene_meta(Path(args.counts_raw))
    gene_order = list(gene_meta['gene_id'].astype(str))

    gene_arrays = None
    if args.resolve_ambiguous:
        print(f'Building gene-extension-{args.gene_extension}bp intervals for '
              f'ambiguity resolution...', flush=True)
        gene_arrays = _build_gene_arrays(args.counts_raw, args.gene_extension)
        print(f'  {sum(len(a[0]) for a in gene_arrays.values())} genes across '
              f'{len(gene_arrays)} contigs', flush=True)

    print(f'Walking {len(list(anno_dir.glob("*.bam.featureCounts.bam")))} '
          f'annotated BAMs (allow_no_umi={args.allow_no_umi}, '
          f'resolve_ambiguous={args.resolve_ambiguous})...', flush=True)
    umi_sets, qc = collect_umis(anno_dir, allow_no_umi=args.allow_no_umi,
                                gene_arrays=gene_arrays)

    X_umi, barcodes = build_matrix(umi_sets, gene_order)

    # Pre-dedup read counts from counts_raw.txt, cells × genes
    read_counts = read_counts.reindex(columns=barcodes, fill_value=0)
    read_counts_T = read_counts.loc[gene_order].to_numpy(dtype=np.int32).T
    X_reads = csr_matrix(read_counts_T)

    qc = qc.set_index('barcode').reindex(barcodes)

    # featureCounts emits semicolon-joined positions when a meta-feature
    # (-g gene_id) aggregates multiple transcripts; the resulting mixed types
    # trip h5py's string writer. Coerce the string-ish var columns to str so
    # the h5ad write succeeds.
    var_df = gene_meta.set_index('gene_id').copy()
    for col in ('Chr', 'Start', 'End', 'Strand'):
        if col in var_df.columns:
            var_df[col] = var_df[col].astype(str)
    if 'Length' in var_df.columns:
        var_df['Length'] = pd.to_numeric(var_df['Length'], errors='coerce').fillna(0).astype(np.int32)

    adata = ad.AnnData(
        X=X_umi,
        obs=qc,
        var=var_df,
        layers={'read_counts': X_reads},
    )
    adata.obs_names = pd.Index(barcodes, name='barcode')
    adata.uns['pipeline_params'] = {
        'primer_trim_run_min': int(args.run_min),
        'umi_len': 10,
        'allow_no_umi': bool(args.allow_no_umi),
        'resolve_ambiguous': bool(args.resolve_ambiguous),
        'gene_extension_bp': int(args.gene_extension) if args.resolve_ambiguous else 0,
        'star_outFilterMatchNminOverLread': 0.3,
        'star_outFilterScoreMinOverLread': 0.3,
        'star_outFilterMultimapNmax': 1,
        'featurecounts_feature_type': 'transcript',
        'featurecounts_fracOverlap': 0.1,
        'featurecounts_countReadPairs': True,
    }

    out = Path(args.out_h5ad)
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    print(f'\nWrote {out}')
    print(f'  shape:            {adata.n_obs} cells × {adata.n_vars} genes')
    print(f'  X nnz (UMI):      {adata.X.nnz}  '
          f'read_counts nnz:  {adata.layers["read_counts"].nnz}')
    print(f'  Median UMIs/cell: '
          f'{np.median(np.asarray(adata.X.sum(axis=1)).ravel()):.0f}')
    print(f'  Median genes/cell:'
          f' {np.median(np.asarray((adata.X > 0).sum(axis=1)).ravel()):.0f}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
