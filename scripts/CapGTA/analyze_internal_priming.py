#!/usr/bin/env python3
"""Test the internal-priming hypothesis.

For each XP's primer-side reads, split by distance from the nearest transcript 3' end:
  - "near_3p":  |offset from nearest tx 3' end| <= near_window (default 200 bp)
  - "internal": offset < -far_threshold (default -500 bp inside the gene)

For each group, compute the fraction of reads whose primer-adjacent base lies
within check_window bp of a genomic A{min_run+,} or T{min_run+,} tract on the
+ strand. If internal priming is real, "internal" reads should overlap genomic
A/T tracts far more often than a random background does; "near_3p" reads should
also overlap tracts (the mRNA polyA tail aligns to a genomic A-rich region).

Also computes an INTERGENIC background by random-sampling positions.
"""
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


def find_tracts(fasta_path: str, min_run: int) -> dict:
    """chrom -> ({'A': (starts, ends), 'T': (starts, ends)}) on + strand, 1-based inclusive."""
    pat_A = re.compile(rb'[Aa]{' + str(min_run).encode() + rb',}')
    pat_T = re.compile(rb'[Tt]{' + str(min_run).encode() + rb',}')
    out = {}
    cur = None; buf = []
    def flush():
        if cur is None: return
        seq = b''.join(buf)
        a = [(m.start() + 1, m.end()) for m in pat_A.finditer(seq)]
        t = [(m.start() + 1, m.end()) for m in pat_T.finditer(seq)]
        out[cur] = {
            'A': (np.array([x[0] for x in a], dtype=np.int64),
                  np.array([x[1] for x in a], dtype=np.int64)),
            'T': (np.array([x[0] for x in t], dtype=np.int64),
                  np.array([x[1] for x in t], dtype=np.int64)),
            'len': len(seq),
        }
    with open(fasta_path, 'rb') as f:
        for line in f:
            if line.startswith(b'>'):
                flush()
                cur = line[1:].strip().split()[0].decode()
                buf = []
            else:
                buf.append(line.strip())
        flush()
    return out


def within_tract(starts: np.ndarray, ends: np.ndarray, pos: int, window: int) -> bool:
    """True if any tract overlaps [pos - window, pos + window]."""
    if starts.size == 0: return False
    lo = pos - window
    hi = pos + window
    idx = np.searchsorted(starts, hi + 1, side='right')
    if idx == 0: return False
    return bool((ends[:idx] >= lo).any())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--offsets-csv', required=True,
                    help='CSV with xp, gene_id, gene_strand, r1_off_nearest_3p, r2_off_nearest_3p + absolute positions')
    ap.add_argument('--counts-raw', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--near-window', type=int, default=200,
                    help='|offset| <= this bp = "near 3-prime end" (default 200)')
    ap.add_argument('--far-threshold', type=int, default=500,
                    help='offset < -(this) bp (upstream) = "internal" (default 500)')
    ap.add_argument('--check-window', type=int, default=30,
                    help='Primer base within this bp of a genomic A/T tract (default 30)')
    ap.add_argument('--min-run', type=int, default=10,
                    help='Minimum genomic tract length (default 10)')
    ap.add_argument('--n-bg-samples', type=int, default=200_000,
                    help='Random genomic positions for the background (default 200k)')
    args = ap.parse_args()

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    print(f'Finding genomic A/T tracts ≥ {args.min_run} bp...', flush=True)
    tracts = find_tracts(args.fasta, args.min_run)
    tot_bp = sum(v['len'] for v in tracts.values())
    tot_A = sum(len(v['A'][0]) for v in tracts.values())
    tot_T = sum(len(v['T'][0]) for v in tracts.values())
    print(f'  {tot_A:,} polyA tracts and {tot_T:,} polyT tracts across {tot_bp:,} bp', flush=True)

    print('Loading gene coords...', flush=True)
    cr = pd.read_csv(args.counts_raw, sep='\t', comment='#', low_memory=False).rename(columns={'Geneid': 'gene_id'})
    gene_meta = {}
    for _, row in cr.iterrows():
        starts = [int(s) for s in str(row['Start']).split(';')]
        ends   = [int(s) for s in str(row['End']).split(';')]
        strand = str(row['Strand']).split(';')[0]
        chrom  = str(row['Chr']).split(';')[0]
        gene_meta[row['gene_id']] = (chrom, min(starts), max(ends), strand)

    print('Loading offsets CSV...', flush=True)
    df = pd.read_csv(args.offsets_csv, low_memory=False)
    print(f'  {len(df):,} rows', flush=True)

    # The CSV from plot_tx_nearest_3p.py has r1/r2_off_nearest_3p computed.
    # Recover the primer-adjacent absolute position: for +strand gene, read_end;
    # for -strand gene, read_start. We rebuild from off_3p + gene_body_3p coords.
    # For primer side (R1 for R1T/R1A, R2 for R2A/R2T):
    XP_PRIMER = {
        'R1T': 'r1',  'R1A': 'r1',
        'R2A': 'r2',  'R2T': 'r2',
    }

    def primer_abs_pos(row):
        xp = row['xp']
        side = XP_PRIMER.get(xp)
        if side is None: return None
        meta = gene_meta.get(row['gene_id'])
        if meta is None: return None
        chrom, gs, ge, strand = meta
        if strand == '+':
            # r{side}_off_3p = r{side}_end - ge → r{side}_end = off + ge
            end = int(row[f'{side}_off_3p']) + ge
            return chrom, end
        else:
            # r{side}_off_3p = gs - r{side}_start → r{side}_start = gs - off
            start = gs - int(row[f'{side}_off_3p'])
            return chrom, start

    # Compute primer-position and group label per row
    print('Classifying rows into near/internal/other...', flush=True)
    results = []
    rng = np.random.default_rng(0)
    xp_group_counts = defaultdict(lambda: {'near_3p': 0, 'internal': 0, 'other': 0})
    xp_group_hits = defaultdict(lambda: {'near_3p': 0, 'internal': 0, 'other': 0})

    # Need near-offset = nearest-tx-3' offset for the primer side
    side_offs = {'R1T': 'r1_off_nearest_3p', 'R1A': 'r1_off_nearest_3p',
                 'R2A': 'r2_off_nearest_3p', 'R2T': 'r2_off_nearest_3p'}

    for i, row in enumerate(df.itertuples(index=False)):
        xp = row.xp
        off_col = side_offs.get(xp)
        if off_col is None: continue
        try:
            off = int(getattr(row, off_col))
        except Exception:
            continue
        if abs(off) >= 1_000_000: continue  # sentinel: no transcript lookup
        if abs(off) <= args.near_window:
            grp = 'near_3p'
        elif off < -args.far_threshold:
            grp = 'internal'
        else:
            grp = 'other'
        xp_group_counts[xp][grp] += 1

        # Only check tract overlap for the two interesting groups
        if grp == 'other':
            continue
        meta = gene_meta.get(row.gene_id)
        if meta is None: continue
        chrom, gs, ge, strand = meta
        if strand == '+':
            pos = int(getattr(row, 'r1_off_3p' if XP_PRIMER[xp] == 'r1' else 'r2_off_3p')) + ge
        else:
            pos = gs - int(getattr(row, 'r1_off_3p' if XP_PRIMER[xp] == 'r1' else 'r2_off_3p'))
        tract = tracts.get(chrom)
        if tract is None: continue
        hit = (within_tract(tract['A'][0], tract['A'][1], pos, args.check_window) or
               within_tract(tract['T'][0], tract['T'][1], pos, args.check_window))
        if hit:
            xp_group_hits[xp][grp] += 1
        if (i + 1) % 500_000 == 0:
            print(f'  [{i+1:,}/{len(df):,}]', flush=True)

    # Background: random genomic positions
    print(f'Background: sampling {args.n_bg_samples:,} random genomic positions...', flush=True)
    bg_hits = 0
    chrom_list = list(tracts.keys())
    chrom_weights = np.array([tracts[c]['len'] for c in chrom_list], dtype=np.float64)
    chrom_weights /= chrom_weights.sum()
    for _ in range(args.n_bg_samples):
        ci = int(rng.choice(len(chrom_list), p=chrom_weights))
        chrom = chrom_list[ci]
        pos = int(rng.integers(1, tracts[chrom]['len'] + 1))
        if (within_tract(tracts[chrom]['A'][0], tracts[chrom]['A'][1], pos, args.check_window) or
            within_tract(tracts[chrom]['T'][0], tracts[chrom]['T'][1], pos, args.check_window)):
            bg_hits += 1
    bg_rate = bg_hits / args.n_bg_samples

    # Report
    print(f'\n=== Fraction of primer-side reads near a genomic A/T tract (≥{args.min_run} bp, window ±{args.check_window} bp) ===')
    print(f'Background rate (random genomic positions): {bg_rate:.3%}\n')
    rows = []
    for xp in ('R1T', 'R1A', 'R2A', 'R2T'):
        n_near = xp_group_counts[xp]['near_3p']
        n_int  = xp_group_counts[xp]['internal']
        h_near = xp_group_hits[xp]['near_3p']
        h_int  = xp_group_hits[xp]['internal']
        f_near = h_near / max(n_near, 1)
        f_int  = h_int / max(n_int, 1)
        print(f'{xp}:')
        print(f'  near 3\' end (|off|≤{args.near_window}bp): {n_near:>8,} reads, '
              f'{h_near:>8,} hit tract ({f_near:.2%})  — enrichment vs bg: {f_near/bg_rate:.2f}×')
        print(f'  internal (off<-{args.far_threshold}bp):     {n_int:>8,} reads, '
              f'{h_int:>8,} hit tract ({f_int:.2%})  — enrichment vs bg: {f_int/bg_rate:.2f}×')
        rows.append({
            'xp': xp,
            'near3p_n': n_near, 'near3p_hits': h_near, 'near3p_frac': f_near,
            'internal_n': n_int, 'internal_hits': h_int, 'internal_frac': f_int,
            'bg_frac': bg_rate,
            'near3p_enrichment_vs_bg': f_near / bg_rate,
            'internal_enrichment_vs_bg': f_int / bg_rate,
        })
    pd.DataFrame(rows).to_csv(out / 'internal_priming_tract_overlap.tsv', sep='\t', index=False)
    print(f'\nWrote {out}/internal_priming_tract_overlap.tsv')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
