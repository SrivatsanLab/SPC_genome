#!/usr/bin/env python3
"""Classify where featureCounts-unassigned R1 reads fall relative to gene bodies.

For each sampled cell BAM, walk primary R1 reads. For reads with XS tag:
  Unassigned_NoFeatures      → compute signed distance to nearest gene
                               (- = before gene start on reference,
                                + = after gene end, inside gene = overlap)
  Unassigned_Overlapping_Length → compute actual overlap fraction with best-overlapping gene
  Unassigned_Ambiguity       → count how many genes the read overlaps

Writes TSVs + a 3-panel plot showing (a) NoFeatures distance distribution,
(b) Overlapping_Length fraction (how close to 0.1 cutoff), (c) Ambiguity gene-count.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam


def build_gene_arrays(counts_raw: str) -> dict:
    """chrom -> (starts_sorted, ends_sorted_by_start, gene_ids).

    Sorted by start. ends[i] corresponds to the same gene as starts[i].
    Supports fast O(log N) nearest-gene queries via np.searchsorted.
    """
    df = pd.read_csv(counts_raw, sep='\t', comment='#')
    df = df.rename(columns={'Geneid': 'gene_id'})
    out: dict = {}
    rows_per_chrom: dict = defaultdict(list)
    for _, row in df.iterrows():
        starts = [int(s) for s in str(row['Start']).split(';')]
        ends   = [int(s) for s in str(row['End']).split(';')]
        chrom  = str(row['Chr']).split(';')[0]
        gs = min(starts)
        ge = max(ends) + 1  # half-open
        if ge > gs:
            rows_per_chrom[chrom].append((gs, ge, row['gene_id']))
    for chrom, rows in rows_per_chrom.items():
        rows.sort(key=lambda r: r[0])
        starts = np.array([r[0] for r in rows], dtype=np.int64)
        ends   = np.array([r[1] for r in rows], dtype=np.int64)
        ids    = np.array([r[2] for r in rows], dtype=object)
        out[chrom] = (starts, ends, ids)
    return out


def _overlapping_indices(starts, ends, pos_lo, pos_hi) -> np.ndarray:
    """Return indices of genes whose [start, end) overlaps read [pos_lo, pos_hi)."""
    # Candidates: starts < pos_hi AND ends > pos_lo
    # starts is sorted. First candidate: all starts < pos_hi.
    hi = np.searchsorted(starts, pos_hi, side='right')
    cand = np.arange(hi)
    if cand.size == 0:
        return cand
    mask = ends[cand] > pos_lo
    return cand[mask]


def nearest_distance(gene_arr: tuple, pos_lo: int, pos_hi: int) -> int:
    starts, ends, _ = gene_arr
    if _overlapping_indices(starts, ends, pos_lo, pos_hi).size:
        return 0
    # Nearest gene whose end <= pos_lo (upstream of read)
    # Nearest gene whose start >= pos_hi (downstream of read)
    left_dist = None
    right_dist = None
    # Left: largest end <= pos_lo. ends is NOT sorted, so scan ends in starts-sorted order
    # Use starts to find candidates ending before pos_lo
    idx_left = np.searchsorted(starts, pos_lo, side='right')
    if idx_left > 0:
        subset_ends = ends[:idx_left]
        max_end = int(subset_ends.max())
        if max_end <= pos_lo:
            left_dist = pos_lo - max_end + 1
    idx_right = np.searchsorted(starts, pos_hi, side='left')
    if idx_right < len(starts):
        min_start = int(starts[idx_right])
        if min_start >= pos_hi:
            right_dist = min_start - pos_hi + 1
    if left_dist is None and right_dist is None:
        return 10_000_000
    if left_dist is None:   return -right_dist   # read is upstream of a gene
    if right_dist is None:  return  left_dist    # read is downstream of a gene
    return left_dist if left_dist <= right_dist else -right_dist


def best_overlap_frac(gene_arr: tuple, pos_lo: int, pos_hi: int) -> float:
    starts, ends, _ = gene_arr
    read_len = pos_hi - pos_lo
    if read_len <= 0: return 0.0
    idx = _overlapping_indices(starts, ends, pos_lo, pos_hi)
    if idx.size == 0: return 0.0
    best = 0.0
    for i in idx:
        ovl = max(0, min(pos_hi, int(ends[i])) - max(pos_lo, int(starts[i])))
        f = ovl / read_len
        if f > best: best = f
    return best


def count_overlaps(gene_arr: tuple, pos_lo: int, pos_hi: int) -> int:
    starts, ends, _ = gene_arr
    return int(_overlapping_indices(starts, ends, pos_lo, pos_hi).size)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--anno-dir', required=True)
    ap.add_argument('--counts-raw', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--n-cells', type=int, default=100)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    print('Building gene coordinate arrays...', flush=True)
    trees = build_gene_arrays(args.counts_raw)
    print(f'  {sum(len(t[0]) for t in trees.values())} genes across {len(trees)} contigs', flush=True)

    bams = sorted(Path(args.anno_dir).glob('*.bam.featureCounts.bam'))[:args.n_cells]
    print(f'Scanning {len(bams)} cells...', flush=True)

    nofeat_dists = []    # distance from nearest gene for NoFeatures reads
    overlap_fracs = []   # overlap frac for Overlapping_Length reads
    amb_counts = []      # # overlapping genes for Ambiguity reads

    for i, bam_path in enumerate(bams):
        try:
            bam = pysam.AlignmentFile(str(bam_path), 'rb')
        except Exception:
            continue
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped: continue
            if read.is_secondary or read.is_supplementary: continue
            if not read.is_read1: continue
            try:
                xs = read.get_tag('XS')
            except KeyError:
                continue
            chrom = read.reference_name
            tree = trees.get(chrom)
            if tree is None: continue
            # pysam is 0-based; our gene coords are 1-based (featureCounts). Convert.
            pos_lo = read.reference_start + 1
            pos_hi = read.reference_end         # 1-based inclusive -> use as half-open
            if xs == 'Unassigned_NoFeatures':
                d = nearest_distance(tree, pos_lo, pos_hi)
                nofeat_dists.append(d)
            elif xs == 'Unassigned_Overlapping_Length':
                overlap_fracs.append(best_overlap_frac(tree, pos_lo, pos_hi))
            elif xs == 'Unassigned_Ambiguity':
                amb_counts.append(count_overlaps(tree, pos_lo, pos_hi))
        bam.close()
        if (i + 1) % 25 == 0:
            print(f'  [{i+1}/{len(bams)}] no_feat={len(nofeat_dists)} ovl={len(overlap_fracs)} amb={len(amb_counts)}',
                  flush=True)

    nofeat = np.asarray(nofeat_dists, dtype=np.int64)
    ovl = np.asarray(overlap_fracs, dtype=np.float32)
    amb = np.asarray(amb_counts, dtype=np.int32)

    np.save(out / 'nofeat_distances.npy', nofeat)
    np.save(out / 'overlap_fracs.npy', ovl)
    np.save(out / 'ambig_counts.npy', amb)

    # ----- summaries -----
    def bin_report(name, arr, bins):
        print(f'\n=== {name} (n={len(arr)}) ===')
        for lo, hi, label in bins:
            if hi is None:
                mask = (arr >= lo)
            elif lo is None:
                mask = (arr < hi)
            else:
                mask = (arr >= lo) & (arr < hi)
            n = int(mask.sum())
            pct = n/max(len(arr),1)*100
            print(f'  {label:<35s}: {n:>8d}  ({pct:.2f}%)')

    bin_report('NoFeatures  — signed distance to nearest gene (bp)', nofeat, [
        (None,     -10000, '< -10 kb upstream'),
        (-10000,   -1000,  '-10 .. -1 kb upstream'),
        (-1000,    -200,   '-1000 .. -200 bp upstream'),
        (-200,     0,      '-200 .. 0 bp (just upstream)'),
        (0,        1,      '= 0 (inside a gene already!?)'),  # should be ~0
        (1,        200,    '1 .. 200 bp downstream'),
        (200,      1000,   '200 .. 1000 bp downstream'),
        (1000,     10000,  '1 .. 10 kb downstream'),
        (10000,    None,   '> 10 kb downstream'),
    ])

    bin_report('Overlapping_Length — best overlap frac (cutoff was 0.10)', ovl, [
        (None, 0.01, '< 1%'),
        (0.01, 0.02, '1-2%'),
        (0.02, 0.05, '2-5%'),
        (0.05, 0.08, '5-8%'),
        (0.08, 0.10, '8-10% (just below cutoff)'),
    ])

    bin_report('Ambiguity — # overlapping gene bodies', amb, [
        (0, 1,  '0'),
        (1, 2,  '1'),
        (2, 3,  '2'),
        (3, 4,  '3'),
        (4, None, '≥4'),
    ])

    # ----- plots -----
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    # (1) NoFeatures distance distribution — zoom to +/- 2 kb
    clipped = np.clip(nofeat, -2000, 2000)
    axes[0].hist(clipped, bins=81, color='#2a7', edgecolor='white')
    axes[0].axvline(0, color='k', ls='--', alpha=0.5, label='gene edge')
    axes[0].set_xlabel('signed distance to nearest gene body (bp)\n(- = read is upstream of a gene; + = downstream)')
    axes[0].set_ylabel('# NoFeatures R1 reads')
    axes[0].set_title(f'NoFeatures distances (n={len(nofeat)})')
    axes[0].legend()
    # (2) Overlap frac for Overlapping_Length
    axes[1].hist(ovl, bins=np.linspace(0, 0.1, 21), color='#a27', edgecolor='white')
    axes[1].axvline(0.1, color='k', ls='--', alpha=0.5, label='current cutoff (10%)')
    axes[1].set_xlabel('best gene overlap fraction')
    axes[1].set_ylabel('# Overlapping_Length R1 reads')
    axes[1].set_title(f'Reads just under 10% overlap cutoff (n={len(ovl)})')
    axes[1].legend()
    # (3) Ambiguity gene counts
    if len(amb):
        bins = np.arange(2, max(5, int(amb.max())+2)) - 0.5
        axes[2].hist(amb, bins=bins, color='#27a', edgecolor='white')
    axes[2].set_xlabel('# gene bodies overlapped')
    axes[2].set_ylabel('# Ambiguity R1 reads')
    axes[2].set_title(f'Ambiguity: multi-gene overlap (n={len(amb)})')

    plt.tight_layout()
    plt.savefig(out / 'unassigned_breakdown.png', dpi=150, bbox_inches='tight')
    print(f'\nWrote plot: {out}/unassigned_breakdown.png')
    print(f'Done. All outputs in {out}/')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
