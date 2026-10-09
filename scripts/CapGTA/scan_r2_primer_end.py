#!/usr/bin/env python3
"""Scan demuxed R2 ends for the RT primer signature.

For a polyT RT primer-captured molecule with primer at the opposite end from
the SPC barcodes (the normal structure), R2 reads from the barcode-side end
through the cDNA and may reach the primer at the END of R2.

What we expect to find at the end of R2 (reading 5'→3' on R2's own strand):
   Case A (R2 reads mRNA-sense strand): ...[cDNA]...[polyA_tail_as_polyA][rc_UMI][rc_adapter]
      → polyA run NEAR the end, then an "rc_UMI" (10bp random), then rev-comp of adapter
   Case B (R2 reads mRNA-antisense strand): ...[cDNA_antisense]...[polyT][UMI][adapter]
      → polyT run near the end, then UMI, then adapter

For both cases we scan the last ~50 bp of R2 and look for an A{15,} or T{15,}
run and an adjacent 10 bp UMI window.
"""
from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path


PAT_T = re.compile(r'T{15,}')
PAT_A = re.compile(r'A{15,}')


def umi_window_ok(window: str) -> bool:
    if len(window) != 10: return False
    if any(b not in 'ACGT' for b in window): return False
    if window.count('A') >= 8 or window.count('T') >= 8: return False
    return True


def iter_fastq(path: str, limit=None):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        n = 0
        while True:
            h = fh.readline()
            if not h: break
            s = fh.readline().strip()
            fh.readline()
            q = fh.readline().strip()
            yield h.strip(), s, q
            n += 1
            if limit and n >= limit: break


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--demux-dir', required=True)
    ap.add_argument('--n-cells', type=int, default=100)
    args = ap.parse_args()

    demux = Path(args.demux_dir)
    pairs = []
    for r2 in sorted(demux.rglob('*_R2.fastq.gz'))[:args.n_cells]:
        r1 = Path(str(r2).replace('_R2.', '_R1.'))
        if r1.exists():
            pairs.append((r1, r2))
    print(f'Scanning {len(pairs)} cells...', flush=True)

    totals = {
        'pairs_scanned':            0,
        'r2_short_lt_40bp':         0,   # too short to look at both ends
        # Scan the LAST 60 bp of R2
        'r2_tail_polyT_run':        0,
        'r2_tail_polyT_umi_ok':     0,
        'r2_tail_polyA_run':        0,
        'r2_tail_polyA_umi_ok':     0,
        # Where in R2 the tail run starts (distance from read end)
        'r2_tail_polyT_dist_0_15':   0,
        'r2_tail_polyT_dist_15_30':  0,
        'r2_tail_polyT_dist_30_60':  0,
        'r2_tail_polyA_dist_0_15':   0,
        'r2_tail_polyA_dist_15_30':  0,
        'r2_tail_polyA_dist_30_60':  0,
    }
    tail_examples = {'polyT': [], 'polyA': []}

    for r1, r2 in pairs:
        for h2, s2, q2 in iter_fastq(str(r2)):
            totals['pairs_scanned'] += 1
            L = len(s2)
            if L < 40:
                totals['r2_short_lt_40bp'] += 1
                continue
            # Look at the LAST 60 bp of R2 for polyT or polyA
            tail = s2[-60:] if L >= 60 else s2
            offset_in_tail_to_read = L - len(tail)   # add this to get read-level index

            # Scan for polyT run in the tail (R2 reading into primer polyT)
            mT = PAT_T.search(tail)
            if mT:
                run_abs_start = offset_in_tail_to_read + mT.start()
                dist_from_end = L - run_abs_start - (mT.end() - mT.start())
                totals['r2_tail_polyT_run'] += 1
                if dist_from_end < 15:
                    totals['r2_tail_polyT_dist_0_15'] += 1
                elif dist_from_end < 30:
                    totals['r2_tail_polyT_dist_15_30'] += 1
                else:
                    totals['r2_tail_polyT_dist_30_60'] += 1
                # UMI should be 10bp AFTER polyT in reads of Case B
                umi_start_abs = run_abs_start + (mT.end() - mT.start())
                if umi_start_abs + 10 <= L:
                    umi = s2[umi_start_abs : umi_start_abs + 10]
                    if umi_window_ok(umi):
                        totals['r2_tail_polyT_umi_ok'] += 1
                        if len(tail_examples['polyT']) < 10:
                            tail_examples['polyT'].append(
                                (h2.split()[0], s2[max(0, run_abs_start-5):], L, run_abs_start))

            # Scan for polyA run in the tail (R2 reading into mRNA polyA tail)
            mA = PAT_A.search(tail)
            if mA:
                run_abs_start = offset_in_tail_to_read + mA.start()
                dist_from_end = L - run_abs_start - (mA.end() - mA.start())
                totals['r2_tail_polyA_run'] += 1
                if dist_from_end < 15:
                    totals['r2_tail_polyA_dist_0_15'] += 1
                elif dist_from_end < 30:
                    totals['r2_tail_polyA_dist_15_30'] += 1
                else:
                    totals['r2_tail_polyA_dist_30_60'] += 1
                umi_start_abs = run_abs_start + (mA.end() - mA.start())
                if umi_start_abs + 10 <= L:
                    umi = s2[umi_start_abs : umi_start_abs + 10]
                    if umi_window_ok(umi):
                        totals['r2_tail_polyA_umi_ok'] += 1
                        if len(tail_examples['polyA']) < 10:
                            tail_examples['polyA'].append(
                                (h2.split()[0], s2[max(0, run_abs_start-5):], L, run_abs_start))

    n = totals['pairs_scanned']
    print(f'\n=== Scanned {n:,} pairs across {len(pairs)} cells ===\n')
    for k, v in totals.items():
        if k == 'pairs_scanned': continue
        pct = v / max(n, 1) * 100
        print(f'  {k:<32s}: {v:>10,}  ({pct:.4f}%)')

    print(f'\n--- 10 R2-tail polyT hits (polyT in last 60bp + UMI window after) ---')
    for name, seq, L, start in tail_examples['polyT']:
        print(f'  {name}  R2_len={L}  polyT_at={start}')
        print(f'    tail: ...{seq}')
    print(f'\n--- 10 R2-tail polyA hits (polyA in last 60bp + UMI window after) ---')
    for name, seq, L, start in tail_examples['polyA']:
        print(f'  {name}  R2_len={L}  polyA_at={start}')
        print(f'    tail: ...{seq}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
