#!/usr/bin/env python3
"""Reparse the per-pair offsets using NEAREST transcript 3' end instead of the
gene-body max end. Tests whether the "primer-side inside gene" signal is from
alternative polyA sites in different transcript isoforms.

Reads the existing r2a_mate_offsets.csv (which has gene_id + gene_strand),
then recomputes r1_off_3p / r2_off_3p vs the nearest per-transcript 3' end
of the assigned gene. Compares the two distributions.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_transcript_3p_per_gene(gtf: str) -> dict:
    """gene_id -> list of (chrom, strand, transcript_3p_pos)"""
    out: dict = {}
    attr_pat_re_gene = __import__('re').compile(r'gene_id "([^"]+)"')
    with open(gtf) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9: continue
            if parts[2] != 'transcript': continue
            chrom = parts[0]
            start, end = int(parts[3]), int(parts[4])
            strand = parts[6]
            m = attr_pat_re_gene.search(parts[8])
            if not m: continue
            gene = m.group(1)
            tx_3p = end if strand == '+' else start
            out.setdefault(gene, []).append((chrom, strand, tx_3p))
    # Convert to sorted np arrays per gene for fast nearest search
    out2 = {}
    for g, lst in out.items():
        if not lst: continue
        strand = lst[0][1]
        positions = np.array(sorted(set(t[2] for t in lst)), dtype=np.int64)
        out2[g] = (strand, positions)
    return out2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--offsets-csv', required=True)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--max-offset', type=int, default=3000)
    args = ap.parse_args()

    print('Loading transcript 3\' ends per gene from GTF...', flush=True)
    tx_lookup = load_transcript_3p_per_gene(args.gtf)
    print(f'  {len(tx_lookup):,} genes', flush=True)

    # Multi-isoform stats
    n_sites_per_gene = np.array([len(v[1]) for v in tx_lookup.values()])
    print(f'\n=== Alternative polyA sites per gene (unique transcript 3\' ends) ===')
    print(f'  genes with 1 site:  {int((n_sites_per_gene == 1).sum()):,}  ({(n_sites_per_gene == 1).mean():.1%})')
    print(f'  genes with 2 sites: {int((n_sites_per_gene == 2).sum()):,}')
    print(f'  genes with 3+ sites:{int((n_sites_per_gene >= 3).sum()):,}')
    print(f'  mean sites/gene = {n_sites_per_gene.mean():.2f}  median = {np.median(n_sites_per_gene):.0f}  max = {int(n_sites_per_gene.max())}')
    multi = [v[1] for v in tx_lookup.values() if len(v[1]) >= 2]
    if multi:
        spans = np.array([int(arr.max() - arr.min()) for arr in multi])
        print(f'  for multi-isoform genes, span between min and max 3\' end: '
              f'median={np.median(spans):.0f}  mean={spans.mean():.0f}  max={int(spans.max())}')

    df = pd.read_csv(args.offsets_csv)

    # Recompute offset from NEAREST transcript 3' end per read
    def nearest_tx_offset(row, pos_col: str) -> int:
        meta = tx_lookup.get(row['gene_id'])
        if meta is None:
            return 1_000_000  # drop
        strand, positions = meta
        pos = row[pos_col]
        if strand == '+':
            # offset = read_end - tx_3p (same convention as before)
            diffs = pos - positions
        else:
            # offset = tx_3p - read_start (same convention: + = past polyA site)
            diffs = positions - pos
        # Pick signed diff with smallest absolute value
        idx = int(np.argmin(np.abs(diffs)))
        return int(diffs[idx])

    print('Recomputing R1 end offset vs nearest transcript 3\' end...', flush=True)
    df['r1_off_nearest_3p'] = df.apply(lambda r: nearest_tx_offset(r, 'r1_end'), axis=1) \
        if False else [None] * len(df)
    # Faster: vectorize per gene_id (group by gene then apply)
    # Group by gene_id for efficiency
    out_r1 = np.empty(len(df), dtype=np.int64)
    out_r2 = np.empty(len(df), dtype=np.int64)
    gene_ids = df['gene_id'].to_numpy()
    r1_end = df.r1_end.to_numpy() if 'r1_end' in df.columns else None
    r2_end = df.r2_end.to_numpy() if 'r2_end' in df.columns else None
    r1_start = df.r1_start.to_numpy() if 'r1_start' in df.columns else None
    r2_start = df.r2_start.to_numpy() if 'r2_start' in df.columns else None

    # We don't have r1_end/r2_end explicitly in the CSV; recover them from the
    # existing r1_off_3p = r1_end - gene_end (if +) or gene_start - r1_start (if -).
    # Actually we only need the signed offset vs nearest, so re-derive the
    # absolute positions.
    # For + strand gene with meta (gs, ge, '+'): r1_off_3p = r1_end - ge → r1_end = r1_off_3p + ge
    # For - strand gene with meta (gs, ge, '-'): r1_off_3p = gs - r1_start → r1_start = gs - r1_off_3p
    # We don't carry gs/ge in the CSV. So we need a gene-coords lookup too.
    # Fall back to loading counts_raw.
    import os
    counts_raw = os.environ.get('COUNTS_RAW')
    if counts_raw is None:
        print('ERROR: set COUNTS_RAW env var to the counts_raw.txt path', flush=True)
        return 2
    cr = pd.read_csv(counts_raw, sep='\t', comment='#').rename(columns={'Geneid': 'gene_id'})
    gene_coords = {}
    for _, row in cr.iterrows():
        starts = [int(s) for s in str(row['Start']).split(';')]
        ends   = [int(s) for s in str(row['End']).split(';')]
        strand = str(row['Strand']).split(';')[0]
        gene_coords[row['gene_id']] = (min(starts), max(ends), strand)

    out_r1_nearest = np.empty(len(df), dtype=np.int64)
    out_r2_nearest = np.empty(len(df), dtype=np.int64)
    out_r1_nearest.fill(10_000_000)
    out_r2_nearest.fill(10_000_000)

    for i, row in enumerate(df.itertuples(index=False)):
        gc = gene_coords.get(row.gene_id)
        tx = tx_lookup.get(row.gene_id)
        if gc is None or tx is None: continue
        gs, ge, strand = gc
        tx_strand, positions = tx
        # Recover absolute positions
        if strand == '+':
            r1_abs = int(row.r1_off_3p) + ge  # read end
            r2_abs = int(row.r2_off_3p) + ge
            d1 = r1_abs - positions
            d2 = r2_abs - positions
        else:
            r1_abs = gs - int(row.r1_off_3p)  # read start
            r2_abs = gs - int(row.r2_off_3p)
            d1 = positions - r1_abs
            d2 = positions - r2_abs
        if d1.size:
            out_r1_nearest[i] = int(d1[np.argmin(np.abs(d1))])
        if d2.size:
            out_r2_nearest[i] = int(d2[np.argmin(np.abs(d2))])
        if (i + 1) % 500_000 == 0:
            print(f'  [{i+1:,}/{len(df):,}]', flush=True)

    df['r1_off_nearest_3p'] = out_r1_nearest
    df['r2_off_nearest_3p'] = out_r2_nearest
    df.to_csv(Path(args.out_dir) / 'r2a_mate_offsets_with_tx_nearest.csv', index=False)

    # Plot comparison: gene-body 3' end vs nearest transcript 3' end for each XP
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    M = args.max_offset
    XP_SIDES = {
        'R1T': ('r1_off_3p', 'r1_off_nearest_3p', 'R1'),
        'R1A': ('r1_off_3p', 'r1_off_nearest_3p', 'R1'),
        'R2A': ('r2_off_3p', 'r2_off_nearest_3p', 'R2'),
        'R2T': ('r2_off_3p', 'r2_off_nearest_3p', 'R2'),
    }
    xp_colors = {'R1T': '#2a7', 'R1A': '#7a2', 'R2A': '#a27', 'R2T': '#27a'}

    for xp in ['R1T', 'R1A', 'R2A', 'R2T']:
        sub = df[df.xp == xp]
        if not len(sub): continue
        prim_body, prim_nearest, prim_lbl = XP_SIDES[xp]
        body = sub[prim_body].to_numpy()
        near = sub[prim_nearest].to_numpy()
        body_in = body[np.abs(body) <= M]
        near_in = near[np.abs(near) <= M]
        fig, axes = plt.subplots(2, 2, figsize=(13, 10))
        bins_full = np.linspace(-M, M, 121)
        bins_zoom = np.linspace(-2000, 500, 111)
        axes[0, 0].hist(body_in, bins=bins_full, color='#888', edgecolor='white')
        axes[0, 0].axvline(0, color='k', ls='--', alpha=0.5, label='gene-body 3\' end')
        axes[0, 0].set_title(f'XP={xp}  {prim_lbl} vs gene-body 3\' end  (n={len(body_in):,})')
        axes[0, 0].set_xlabel(f'{prim_lbl} offset from gene-body 3\' end (bp)')
        axes[0, 0].legend()
        axes[0, 1].hist(body_in[(body_in >= -2000) & (body_in <= 500)],
                        bins=bins_zoom, color='#888', edgecolor='white')
        axes[0, 1].axvline(0, color='k', ls='--', alpha=0.5)
        axes[0, 1].set_title(f'XP={xp}  gene-body 3\' end — Zoom')
        axes[0, 1].set_xlabel(f'{prim_lbl} offset (bp)')
        c = xp_colors[xp]
        axes[1, 0].hist(near_in, bins=bins_full, color=c, edgecolor='white')
        axes[1, 0].axvline(0, color='k', ls='--', alpha=0.5, label='nearest transcript 3\' end')
        axes[1, 0].set_title(f'XP={xp}  {prim_lbl} vs NEAREST transcript 3\' end  (n={len(near_in):,})')
        axes[1, 0].set_xlabel(f'{prim_lbl} offset from nearest tx 3\' end (bp)')
        axes[1, 0].legend()
        axes[1, 1].hist(near_in[(near_in >= -2000) & (near_in <= 500)],
                        bins=bins_zoom, color=c, edgecolor='white')
        axes[1, 1].axvline(0, color='k', ls='--', alpha=0.5)
        axes[1, 1].set_title(f'XP={xp}  nearest transcript 3\' end — Zoom')
        axes[1, 1].set_xlabel(f'{prim_lbl} offset (bp)')
        plt.tight_layout()
        fp = out / f'primer_3p_body_vs_nearest_{xp}.png'
        plt.savefig(fp, dpi=140, bbox_inches='tight')
        plt.close(fig)
        print(f'Wrote {fp}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
