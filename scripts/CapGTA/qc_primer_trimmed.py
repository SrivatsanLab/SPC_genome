#!/usr/bin/env python3
"""QC + 5'/3' end position plots for the primer-trimmed pipeline output.

Reports:
  - per-cell UMIs, genes, pairs, and alignment-failure counts (saved to CSV + printed summary)
  - a plot of where R2 reads land relative to the 5' end of their assigned gene
  - a plot of where R1 (cDNA remainder after primer trim) lands relative to
    the 3' end of its assigned gene — the 3' end is the polyA site, so a
    primer-captured R1 should align with its "primer-adjacent" base (5' of the
    sequenced read, mapping to read.reference_end on + strand genes) sitting
    at or just upstream of the gene 3' end

For both plots we walk each annotated BAM, filter to XS=Assigned, pull gene_id
from XT, and compute the signed offset in gene (5'→3') coordinates.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam


def _barcode_from_fcbam(p: Path) -> str:
    name = p.name
    if name.endswith('.bam.featureCounts.bam'):
        return name[:-len('.bam.featureCounts.bam')]
    return name.removesuffix('.bam').removesuffix('.featureCounts')


def _build_gene_lookup(gene_meta: pd.DataFrame) -> dict:
    """gene_id -> (gene_start, gene_end, strand)  (1-based inclusive)."""
    def _parse_first(v):
        s = str(v)
        if ';' in s:
            s = s.split(';')[0]
        return s
    gm = gene_meta.copy()
    gm['first_start'] = gm['Start'].apply(lambda s: int(_parse_first(s)))
    gm['first_end']   = gm['End'].apply(lambda s: int(_parse_first(s)))
    gm['first_strand'] = gm['Strand'].apply(lambda s: _parse_first(s))
    return {
        row['gene_id']: (row['first_start'], row['first_end'], row['first_strand'])
        for _, row in gm.iterrows()
    }


def collect_offsets(anno_dir: Path,
                    gene_meta: pd.DataFrame,
                    max_offset: int,
                    sample_cells: int | None) -> tuple[np.ndarray, np.ndarray]:
    """Walk each annotated BAM once; return two arrays:
       r2_5p_offsets: R2 start offset from gene 5' end (TSS)
       r1_3p_offsets: R1 "primer-adjacent" end offset from gene 3' end (polyA)

    Gene coordinate convention: offset = signed distance in gene 5'→3' direction.
    Negative = inside the gene (upstream of 3' end OR downstream of 5' end).
    Positive = past the gene's 3' end (polyA site) OR upstream of its 5' end.

    For R2 vs 5' end:
        + strand: off = read_start - gene_start
        - strand: off = gene_end   - read_end
    For R1 vs 3' end (the base closest to the RT primer polyT):
        + strand (R1 cDNA aligns to - strand genome, is_reverse=True):
            primer-adjacent base = read.reference_end  → off = read_end - gene_end
        - strand (R1 cDNA aligns to + strand genome, is_reverse=False):
            primer-adjacent base = read.reference_start → off = gene_start - read_start
    """
    gene_lookup = _build_gene_lookup(gene_meta)

    r2_offsets = []
    r1_offsets = []
    bams = sorted(anno_dir.glob('*.bam.featureCounts.bam'))
    if sample_cells is not None:
        bams = bams[:sample_cells]
    for i, bam_path in enumerate(bams):
        try:
            bam = pysam.AlignmentFile(str(bam_path), 'rb')
        except Exception:
            continue
        try:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped: continue
                if read.is_secondary or read.is_supplementary: continue
                try:
                    if read.get_tag('XS') != 'Assigned': continue
                    gene_id = read.get_tag('XT')
                except KeyError:
                    continue
                meta = gene_lookup.get(gene_id)
                if meta is None: continue
                gene_start, gene_end, strand = meta
                read_start_1 = read.reference_start + 1   # 1-based
                read_end_1   = read.reference_end         # 1-based inclusive

                if read.is_read2:
                    if strand == '+':
                        off = read_start_1 - gene_start
                    else:
                        off = gene_end - read_end_1
                    if abs(off) <= max_offset:
                        r2_offsets.append(off)
                elif read.is_read1:
                    if strand == '+':
                        off = read_end_1 - gene_end
                    else:
                        off = gene_start - read_start_1
                    if abs(off) <= max_offset:
                        r1_offsets.append(off)
        finally:
            bam.close()
        if (i + 1) % 100 == 0:
            print(f'  [{i+1}/{len(bams)}] cumulative offsets: R2_5p={len(r2_offsets)}  R1_3p={len(r1_offsets)}',
                  flush=True)
    return np.asarray(r2_offsets, dtype=np.int32), np.asarray(r1_offsets, dtype=np.int32)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--counts-raw', required=True,
                    help='featureCounts counts_raw.txt (for gene_id Start/End/Strand)')
    ap.add_argument('--anno-dir', required=True,
                    help='featurecounts_annotated/ containing <cell>.bam.featureCounts.bam')
    ap.add_argument('--star-logs-dir', required=True,
                    help='Dir of per-cell STAR Log.final.out — used for unaligned counts')
    ap.add_argument('--extract-stats', default=None,
                    help='all_cells_stats.tsv from the primer extraction (optional)')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--max-offset', type=int, default=3000)
    ap.add_argument('--sample-cells-for-plot', type=int, default=None,
                    help='Sample first N cells for the 5\' plot (default: all)')
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # ---- Per-cell UMI/gene counts (from h5ad) ------------------------------
    a = ad.read_h5ad(args.h5ad)
    umis_per_cell  = np.asarray(a.X.sum(axis=1)).ravel()
    genes_per_cell = np.asarray((a.X > 0).sum(axis=1)).ravel()
    reads_per_cell = np.asarray(a.layers['read_counts'].sum(axis=1)).ravel() \
                        if 'read_counts' in a.layers else np.zeros_like(umis_per_cell)

    pc = pd.DataFrame({
        'barcode':         a.obs_names.values,
        'n_pairs':         a.obs.get('n_pairs',           pd.Series(np.nan, index=a.obs_names)).values,
        'n_pairs_assigned':a.obs.get('n_pairs_assigned',  pd.Series(np.nan, index=a.obs_names)).values,
        'n_reads_mapped':  a.obs.get('n_reads_mapped',    pd.Series(np.nan, index=a.obs_names)).values,
        'n_reads_count_matrix': reads_per_cell,
        'n_umis':          umis_per_cell,
        'n_genes_detected':genes_per_cell,
    })

    # ---- STAR per-cell input/mapped/unmapped from Log.final.out -----------
    star_dir = Path(args.star_logs_dir)
    star_rows = []
    for lf in sorted(star_dir.glob('*_Log.final.out')):
        cell = lf.name.replace('_Log.final.out', '')
        d = {'barcode': cell}
        for line in lf.read_text().splitlines():
            if '|' not in line: continue
            k, v = [s.strip() for s in line.split('|', 1)]
            try:
                v = float(v.rstrip('%'))
                if v.is_integer(): v = int(v)
            except ValueError:
                pass
            d[k] = v
        star_rows.append(d)
    star_df = pd.DataFrame(star_rows)

    if len(star_df):
        _cols_wanted = {
            'Number of input reads':                       'star_input_reads',
            'Uniquely mapped reads number':                'star_unique_mapped',
            'Number of reads mapped to multiple loci':     'star_multi_mapped',
            'Number of reads mapped to too many loci':     'star_too_many_loci',
            'Number of reads unmapped: too many mismatches': 'star_unmapped_mm',
            'Number of reads unmapped: too short':         'star_unmapped_too_short',
            'Number of reads unmapped: other':             'star_unmapped_other',
        }
        _cols_present = {k: v for k, v in _cols_wanted.items() if k in star_df.columns}
        star_slim = star_df[['barcode'] + list(_cols_present.keys())].rename(columns=_cols_present)
        for c in _cols_present.values():
            star_slim[c] = pd.to_numeric(star_slim[c], errors='coerce')
        pc = pc.merge(star_slim, on='barcode', how='left')
        pc['star_unmapped_total'] = pc[['star_unmapped_mm','star_unmapped_too_short','star_unmapped_other']].sum(axis=1)
    else:
        print('[warn] no STAR Log.final.out files found', flush=True)

    # ---- Join extraction stats --------------------------------------------
    if args.extract_stats and Path(args.extract_stats).exists():
        es = pd.read_csv(args.extract_stats, sep='\t')
        pc = pc.merge(
            es[['barcode','total_pairs','pairs_no_polyT','pairs_bad_umi',
                'pairs_dropped_short_after_trim','pairs_kept']],
            on='barcode', how='left'
        )

    pc.to_csv(out / 'per_cell_qc.csv', index=False)

    # ---- Summary stats printed to stdout + file ---------------------------
    def _desc(col):
        v = pd.to_numeric(pc[col], errors='coerce').dropna()
        if len(v) == 0:
            return 'n=0'
        return (f'n={len(v)}  '
                f'median={np.median(v):.0f}  mean={v.mean():.1f}  '
                f'p25={np.percentile(v,25):.0f}  p75={np.percentile(v,75):.0f}  '
                f'max={int(v.max())}')

    lines = [f'cells: {len(pc)}', '']
    for col in ['n_umis', 'n_genes_detected', 'n_reads_count_matrix',
                'n_pairs', 'n_pairs_assigned',
                'star_input_reads', 'star_unique_mapped',
                'star_too_many_loci', 'star_unmapped_too_short', 'star_unmapped_total',
                'total_pairs', 'pairs_kept']:
        if col in pc.columns:
            lines.append(f'{col:>30s}: {_desc(col)}')
    txt = '\n'.join(lines)
    print(txt)
    (out / 'summary.txt').write_text(txt + '\n')

    # ---- 5' + 3' end offset plots ------------------------------------------
    print('\nCollecting R2-vs-5\' and R1-vs-3\' offsets...', flush=True)
    gene_meta = pd.read_csv(args.counts_raw, sep='\t', comment='#')
    gene_meta = gene_meta.rename(columns={'Geneid': 'gene_id'})
    r2_off, r1_off = collect_offsets(Path(args.anno_dir), gene_meta,
                                     max_offset=args.max_offset,
                                     sample_cells=args.sample_cells_for_plot)
    print(f'Collected {len(r2_off)} R2 offsets, {len(r1_off)} R1 offsets')

    np.save(out / 'r2_offsets_5prime.npy', r2_off)
    np.save(out / 'r1_offsets_3prime.npy', r1_off)

    def _plot_offset(ax_full, ax_zoom, offsets, label, xlabel_full,
                     xlabel_zoom, zoom_range, color, anchor_label):
        bins_full = np.linspace(-args.max_offset, args.max_offset, 121)
        ax_full.hist(offsets, bins=bins_full, color=color, edgecolor='white')
        ax_full.axvline(0, color='k', ls='--', alpha=0.5, label=anchor_label)
        ax_full.set_xlabel(xlabel_full)
        ax_full.set_ylabel('# reads (all cells pooled)')
        ax_full.set_title(f'{label}  (n={len(offsets)})')
        ax_full.legend(loc='upper right')

        bins_zoom = np.linspace(zoom_range[0], zoom_range[1], 111)
        ax_zoom.hist(offsets, bins=bins_zoom, color=color, edgecolor='white')
        ax_zoom.axvline(0, color='k', ls='--', alpha=0.5)
        ax_zoom.set_xlabel(xlabel_zoom)
        ax_zoom.set_ylabel('# reads')
        ax_zoom.set_title(f'Zoom: {zoom_range[0]} to {zoom_range[1]} bp')

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    _plot_offset(
        axes[0, 0], axes[0, 1], r2_off,
        label="R2 start vs assigned gene 5' end",
        xlabel_full="R2 start offset from gene 5' end (bp; + = inside gene downstream of TSS)",
        xlabel_zoom="R2 start offset from gene 5' end (bp)",
        zoom_range=(-200, 2000), color='#2a7',
        anchor_label="gene 5' end (TSS)",
    )
    _plot_offset(
        axes[1, 0], axes[1, 1], r1_off,
        label="R1 primer-adjacent end vs assigned gene 3' end",
        xlabel_full="R1 offset from gene 3' end (bp; + = past polyA site, - = inside gene upstream of 3' end)",
        xlabel_zoom="R1 offset from gene 3' end (bp)",
        zoom_range=(-2000, 500), color='#a27',
        anchor_label="gene 3' end (polyA site)",
    )
    plt.tight_layout()
    plt.savefig(out / 'offset_from_gene_ends.png', dpi=150, bbox_inches='tight')
    print(f'Wrote {out}/offset_from_gene_ends.png')

    # Position-density summaries
    for label, offsets, bins_desc, outfile in [
        ("R2 start vs 5' end",  r2_off, [
            ('>5000 bp upstream',          lambda o: o < -5000),
            ('-5000..-500 bp',             lambda o: (o >= -5000) & (o < -500)),
            ('-500..0 bp just upstream',   lambda o: (o >= -500) & (o < 0)),
            ('0..500 bp downstream',       lambda o: (o >= 0) & (o < 500)),
            ('500..2000 bp',               lambda o: (o >= 500) & (o < 2000)),
            ('2000..5000 bp',              lambda o: (o >= 2000) & (o < 5000)),
            ('>5000 bp downstream',        lambda o: o >= 5000),
        ], 'r2_offset_bins.tsv'),
        ("R1 primer-end vs 3' end", r1_off, [
            ('>2000 bp upstream of 3\' end',  lambda o: o < -2000),
            ('-2000..-500 bp',                lambda o: (o >= -2000) & (o < -500)),
            ('-500..0 bp (just inside)',      lambda o: (o >= -500) & (o < 0)),
            ('0..200 bp past 3\' end',        lambda o: (o >= 0) & (o < 200)),
            ('200..1000 bp past',             lambda o: (o >= 200) & (o < 1000)),
            ('>1000 bp past',                 lambda o: o >= 1000),
        ], 'r1_offset_bins.tsv'),
    ]:
        print(f'\n{label} distribution (n={len(offsets)}):')
        with open(out / outfile, 'w') as fh:
            fh.write('bin\tn\tfrac\n')
            for bin_label, pred in bins_desc:
                n = int(pred(offsets).sum()) if len(offsets) else 0
                frac = n / max(len(offsets), 1)
                print(f'  {bin_label:<30s}: {n:>8d}  ({frac:.2%})')
                fh.write(f'{bin_label}\t{n}\t{frac:.4f}\n')

    print(f'\nAll outputs in {out}/')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
