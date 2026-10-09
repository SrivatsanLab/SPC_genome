#!/usr/bin/env python3
"""Scan demuxed per-cell R2 FASTQs for primer signatures near the read start.

The spc-demultiplex pipeline strips the 45 bp barcode block from R2, so
position 0 of the demuxed R2 insert is where the molecule's cDNA begins
(from the barcode-side end of the fragment).

If the library is non-directional, some molecules will have the RT primer
tail at the R2-side of the fragment, meaning:
  (a) R2 reads from first-strand cDNA → starts with [UMI][polyT][cDNA]
      (same structure as the primer-side R1 we already extract)
  (b) R2 reads from second-strand cDNA → starts with [polyA][rc_UMI][rc_adapter_seq]
      (read goes INTO the primer tail from the barcode end — only if fragment
       is short enough that R2 reaches the primer end)

We check both cases: polyT at the SAME expected offset (7-12) and polyA at
the SAME offset. For each hit we also check the UMI window.
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
            fh.readline()  # +
            q = fh.readline().strip()
            yield h.strip(), s, q
            n += 1
            if limit and n >= limit: break


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--demux-dir', required=True,
                    help='Directory of per-cell demuxed FASTQs (spc-demultiplex output)')
    ap.add_argument('--n-cells', type=int, default=100)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    demux = Path(args.demux_dir)
    # Find per-cell R2 FASTQs; also open R1 to compare orientation
    r2_files = []
    for root in sorted(demux.rglob('*_R2.fastq.gz'))[:args.n_cells]:
        r1 = Path(str(root).replace('_R2.', '_R1.'))
        if r1.exists():
            r2_files.append((r1, root))
    print(f'Found {len(r2_files)} R1/R2 pairs; scanning...', flush=True)

    totals = {
        'pairs_scanned': 0,
        'r1_T_pos7_12': 0,
        'r1_T_pos7_12_umi_ok': 0,
        'r2_T_pos7_12': 0,
        'r2_T_pos7_12_umi_ok': 0,
        'r2_A_pos7_12': 0,
        'r2_A_pos7_12_umi_ok': 0,
        'r2_T_pos0_6':  0,     # polyT very close to R2 start (no UMI slot)
        'r2_A_pos0_6':  0,
    }
    examples = {'r2_polyT': [], 'r2_polyA': []}

    for r1, r2 in r2_files:
        g1 = iter_fastq(str(r1))
        g2 = iter_fastq(str(r2))
        for (h1, s1, q1), (h2, s2, q2) in zip(g1, g2):
            totals['pairs_scanned'] += 1

            # ---- R1 primer check (same rule as extraction: polyT ≥15 at pos 7-12) ----
            m = PAT_T.search(s1)
            if m and 7 <= m.start() <= 12:
                totals['r1_T_pos7_12'] += 1
                if m.start() >= 10:
                    umi = s1[m.start()-10 : m.start()]
                    if umi_window_ok(umi):
                        totals['r1_T_pos7_12_umi_ok'] += 1

            # ---- R2 polyT check (primer tail orientation same as R1) ----
            m = PAT_T.search(s2)
            if m and 7 <= m.start() <= 12:
                totals['r2_T_pos7_12'] += 1
                if m.start() >= 10:
                    umi = s2[m.start()-10 : m.start()]
                    if umi_window_ok(umi):
                        totals['r2_T_pos7_12_umi_ok'] += 1
                        if len(examples['r2_polyT']) < 10:
                            examples['r2_polyT'].append((h2.split()[0], s2[:60]))
            if m and 0 <= m.start() <= 6:
                totals['r2_T_pos0_6'] += 1

            # ---- R2 polyA check (primer tail reversed; "primer→cDNA" from barcode side) ----
            m = PAT_A.search(s2)
            if m and 7 <= m.start() <= 12:
                totals['r2_A_pos7_12'] += 1
                if m.start() >= 10:
                    umi = s2[m.start()-10 : m.start()]
                    if umi_window_ok(umi):
                        totals['r2_A_pos7_12_umi_ok'] += 1
                        if len(examples['r2_polyA']) < 10:
                            examples['r2_polyA'].append((h2.split()[0], s2[:60]))
            if m and 0 <= m.start() <= 6:
                totals['r2_A_pos0_6'] += 1

    # ---- Print summary ----
    n = totals['pairs_scanned']
    print(f'\n=== Scanned {n:,} pairs across {len(r2_files)} cells ===\n')
    for k, v in totals.items():
        if k == 'pairs_scanned': continue
        pct = v / max(n, 1) * 100
        print(f'  {k:<28s}: {v:>10,}  ({pct:.3f}%)')

    print(f'\n--- 10 example R2 primer-polyT hits (polyT at 7-12 + good UMI) ---')
    for name, seq in examples['r2_polyT']:
        print(f'  {name}: {seq}')
    print(f'\n--- 10 example R2 primer-polyA hits (polyA at 7-12 + good UMI) ---')
    for name, seq in examples['r2_polyA']:
        print(f'  {name}: {seq}')

    with open(out / 'r2_primer_scan_summary.tsv', 'w') as fh:
        for k, v in totals.items():
            fh.write(f'{k}\t{v}\n')
    print(f'\nWrote {out}/r2_primer_scan_summary.tsv')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
