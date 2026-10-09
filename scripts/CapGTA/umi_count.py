#!/usr/bin/env python3
"""Walk featureCounts-annotated per-cell BAMs, dedup UMIs, build a UMI counts h5ad.

Each read in a `featurecounts_annotated/<barcode>.bam.featureCounts.bam` carries:
    UB:Z:<umi>                (added by umi_tag_filter.py)
    XS:Z:<status>             (featureCounts: 'Assigned' or 'Unassigned_*')
    XT:Z:<gene_id>            (featureCounts: present when XS='Assigned')
    XN:i:<n_overlapping>      (featureCounts)

For each (barcode, gene_id) we count the number of distinct UMIs → the UMI count.
Reads lacking XT (unassigned) or UB (shouldn't happen post-filter) are skipped.

The output h5ad is cells × genes with:
    X                             sparse CSR, UMI counts (int32)
    layers['read_counts']         pre-dedup read counts from featureCounts counts_raw.txt
    obs                           per-cell QC (n_reads_rna, n_umis, reads_per_umi, ...)
    var                           gene_id (plus length, chr if available from counts_raw.txt)
    uns['filter_params']          {'run_min', 'umi_len', 'max_mm'}
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pysam
from scipy.sparse import csr_matrix


def _barcode_from_fcbam(p: Path) -> str:
    name = p.name
    # featureCounts output: <input>.featureCounts.bam
    if name.endswith('.bam.featureCounts.bam'):
        return name[:-len('.bam.featureCounts.bam')]
    return name.removesuffix('.bam').removesuffix('.featureCounts')


def collect_umis(anno_dir: Path) -> tuple[dict, pd.DataFrame]:
    """Walk per-cell annotated BAMs. Returns:
        umi_sets[barcode][gene_id] = set of UMI strings
        qc_rows: list of per-cell QC dicts
    """
    bams = sorted(anno_dir.glob('*.bam.featureCounts.bam'))
    if not bams:
        raise RuntimeError(f'No annotated BAMs under {anno_dir}')
    umi_sets: dict[str, dict[str, set]] = {}
    qc_rows = []
    for i, bam_path in enumerate(bams):
        bc = _barcode_from_fcbam(bam_path)
        if bc in umi_sets:
            # duplicated barcode (shouldn't happen) — skip
            continue
        sets: dict[str, set] = defaultdict(set)
        n_reads = 0
        n_assigned = 0
        n_umi_tagged = 0
        n_assigned_umi = 0
        try:
            bam = pysam.AlignmentFile(str(bam_path), 'rb')
        except Exception as e:
            print(f'[skip] {bc}: open failed: {e}', file=sys.stderr)
            continue
        try:
            for read in bam.fetch(until_eof=True):
                n_reads += 1
                try:
                    umi = read.get_tag('UB')
                    n_umi_tagged += 1
                except KeyError:
                    continue
                try:
                    gene = read.get_tag('XT')
                except KeyError:
                    continue
                n_assigned += 1
                n_assigned_umi += 1
                sets[gene].add(umi)
        finally:
            bam.close()
        umi_sets[bc] = sets
        n_umis_total = sum(len(s) for s in sets.values())
        qc_rows.append({
            'barcode': bc,
            'n_reads_rna': n_reads,
            'n_umi_tagged_reads': n_umi_tagged,
            'n_gene_assigned_reads': n_assigned,
            'n_umis': n_umis_total,
            'n_genes_detected': len(sets),
            'reads_per_umi': (n_assigned_umi / n_umis_total) if n_umis_total else np.nan,
        })
        if (i + 1) % 100 == 0:
            print(f'  [{i+1}/{len(bams)}] {bc}: reads={n_reads} '
                  f'assigned={n_assigned} umis={n_umis_total}', flush=True)
    return umi_sets, pd.DataFrame(qc_rows)


def build_umi_matrix(umi_sets: dict[str, dict[str, set]],
                     gene_order: list[str]) -> tuple[csr_matrix, list[str]]:
    barcodes = sorted(umi_sets.keys())
    gene_ix = {g: i for i, g in enumerate(gene_order)}
    rows, cols, data = [], [], []
    for ri, bc in enumerate(barcodes):
        for gene, umis in umi_sets[bc].items():
            gi = gene_ix.get(gene)
            if gi is None:
                continue
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


def load_read_counts(counts_raw: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse the featureCounts raw matrix.

    Returns (gene_meta, counts_df) where counts_df is gene × cell read counts.
    """
    df = pd.read_csv(counts_raw, sep='\t', comment='#')
    bam_cols = df.columns[6:]
    barcodes = [os.path.basename(c).removesuffix('.bam') for c in bam_cols]
    gene_meta = df[['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']].copy()
    gene_meta = gene_meta.rename(columns={'Geneid': 'gene_id'})
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
                    help='featureCounts raw output (counts_raw.txt) for gene metadata and read counts layer')
    ap.add_argument('--filter-stats', default=None,
                    help='Optional filter_stats.tsv from umi_tag_filter_array.sh (adds obs cols)')
    ap.add_argument('--out-h5ad', required=True)
    ap.add_argument('--run-min', type=int, default=6)
    args = ap.parse_args()

    anno_dir = Path(args.anno_dir)
    gene_meta, read_counts = load_read_counts(Path(args.counts_raw))
    gene_order = list(gene_meta['gene_id'].astype(str))

    print(f'Collecting UMIs from {anno_dir}...', flush=True)
    umi_sets, qc = collect_umis(anno_dir)
    print(f'Collected UMIs from {len(umi_sets)} cells', flush=True)

    X_umi, barcodes = build_umi_matrix(umi_sets, gene_order)

    # Align read_counts (gene × cell) to the barcode order used in X_umi
    rc_cols = [c for c in barcodes if c in read_counts.columns]
    missing = set(barcodes) - set(read_counts.columns)
    if missing:
        print(f'  warning: {len(missing)} barcodes in annotated BAMs not in counts_raw.txt',
              file=sys.stderr)
    read_counts = read_counts.reindex(columns=barcodes, fill_value=0)
    # cells × genes
    read_counts_T = read_counts.loc[gene_order].to_numpy(dtype=np.int32).T
    X_reads = csr_matrix(read_counts_T)

    # obs
    qc = qc.set_index('barcode').reindex(barcodes)
    if args.filter_stats and Path(args.filter_stats).exists():
        fs = pd.read_csv(args.filter_stats, sep='\t').drop_duplicates('barcode').set_index('barcode')
        qc = qc.join(fs, how='left', rsuffix='_filter')

    adata = ad.AnnData(
        X=X_umi,
        obs=qc,
        var=gene_meta.set_index('gene_id'),
        layers={'read_counts': X_reads},
    )
    adata.obs_names = pd.Index(barcodes, name='barcode')
    adata.uns['filter_params'] = {
        'run_min': int(args.run_min),
        'umi_len': 9,
        'max_mm': 1,
    }

    out = Path(args.out_h5ad)
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    print(f'\nWrote {out}')
    print(f'  shape:   {adata.n_obs} cells × {adata.n_vars} genes')
    print(f'  X nnz (UMI): {adata.X.nnz}  '
          f'read_counts nnz: {adata.layers["read_counts"].nnz}')
    print(f'  Median UMIs/cell: {np.median(np.asarray(adata.X.sum(axis=1)).ravel()):.0f}')
    print(f'  Median genes/cell: {np.median(np.asarray((adata.X > 0).sum(axis=1)).ravel()):.0f}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
