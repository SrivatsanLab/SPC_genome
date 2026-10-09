#!/usr/bin/env python3
"""Filter a single-cell BAM to primer-signature RNA reads and tag them with UMIs.

For each primary non-duplicate mapped read, we look at the soft-clipped
sequence adjacent to the alignment. The RT primer tail leaves a polyT at the
inner edge of the left SC (and a polyA, with reverse-complemented UMI, at the
inner edge of the right SC). We require:

  - polyA/T run length >= L (default 6)                                 [--run-min]
  - a 9 bp window on the far side of the run that looks like a UMI     [DNA bases,
    (not another A/T stretch; all ACGT)                                 not >= 7 A's or T's]

The UMI window is 9 bp (not 10) because the base adjacent to the primer polyT is
ambiguous: a UMI base of T there would be absorbed into the polyT-run detector,
and a UMI base of A there is tolerated as a 1 bp mismatch inside the run. Taking
the 9 bp just outside the ambiguous position gives a window that is unaffected
by which base happened to be at the polyT boundary. All reads from the same
molecule share that boundary base, so dedup is unaffected.

If both conditions pass, the read is written to the output BAM with:
    UB:Z:<umi>       # 9 bp UMI (reverse-complemented for right-side hits so
                     # left and right UMIs share a common orientation)
    XP:A:L or R      # which end carried the primer signature
    XR:i:<run_len>   # observed polyA/T run length

Reads failing either filter are dropped. Per-cell summary stats are written to a
TSV sidecar (unless --no-stats).
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pysam


UMI_LEN = 9
MAX_MM = 1
_RC = str.maketrans('ACGTN', 'TGCAN')


def _revcomp(s: str) -> str:
    return s.translate(_RC)[::-1]


def _longest_run_from_edge(seq: str, base: str, from_end: str) -> int:
    """Length of the longest base-run anchored at one edge of seq.

    from_end='right' -> anchor at seq[-1] (left SC: polyT at alignment side).
    from_end='left'  -> anchor at seq[0]  (right SC: polyA at alignment side).
    Allows up to MAX_MM non-matching bases interior to the run.
    """
    if not seq:
        return 0
    s = seq[::-1] if from_end == 'right' else seq
    if s[0] != base:
        return 0
    run = 0
    mm = 0
    for ch in s:
        if ch == base:
            run += 1
        else:
            if mm >= MAX_MM:
                break
            mm += 1
            run += 1
    while run > 0 and s[run - 1] != base:
        run -= 1
    return run


def _umi_window_ok(window: str) -> bool:
    if len(window) != UMI_LEN:
        return False
    if any(b not in 'ACGT' for b in window):
        return False
    # Reject windows that are themselves ~all A or ~all T (would be false UMIs
    # arising from the polyA/T region bleeding past our detector)
    if window.count('T') >= 7 or window.count('A') >= 7:
        return False
    return True


def classify_read(read, run_min: int, require_umi: bool = True):
    """Return (umi, side, run_len) if read passes, else (None, None, 0).

    Prefers the left SC when both ends qualify. With require_umi=False, a read
    that has a polyA/T run >= run_min at the inner SC edge passes even if the
    adjacent 9 bp don't look like a valid UMI (or the SC isn't long enough to
    carry a UMI). In that case the returned 'umi' is None but 'side' and
    'run_len' are set.
    """
    if (read.is_secondary or read.is_supplementary or read.is_duplicate
            or read.is_unmapped):
        return None, None, 0
    cig = read.cigartuples
    seq = read.query_sequence
    if not cig or seq is None:
        return None, None, 0

    if cig[0][0] == 4:
        sc_len = cig[0][1]
        sc = seq[:sc_len]
        run_T = _longest_run_from_edge(sc, 'T', from_end='right')
        if run_T >= run_min:
            umi_end = sc_len - run_T
            if umi_end >= UMI_LEN:
                umi = sc[umi_end - UMI_LEN:umi_end]
                if _umi_window_ok(umi):
                    return umi, 'L', run_T
            if not require_umi:
                return None, 'L', run_T

    if cig[-1][0] == 4:
        sc_len = cig[-1][1]
        sc = seq[-sc_len:]
        run_A = _longest_run_from_edge(sc, 'A', from_end='left')
        if run_A >= run_min:
            umi_start = run_A
            if sc_len - umi_start >= UMI_LEN:
                umi = sc[umi_start:umi_start + UMI_LEN]
                if _umi_window_ok(umi):
                    return _revcomp(umi), 'R', run_A
            if not require_umi:
                return None, 'R', run_A

    return None, None, 0


def filter_bam(src_bam: str, out_bam: str, run_min: int,
               require_umi: bool = True) -> dict:
    stats = dict(
        barcode=Path(src_bam).name.removesuffix('.bam'),
        n_reads=0, n_primary=0, n_mapped=0,
        n_pass=0, n_pass_left=0, n_pass_right=0,
        n_pass_with_umi=0,
        sum_run_len=0,
    )
    bam = pysam.AlignmentFile(src_bam, 'rb')
    out = pysam.AlignmentFile(out_bam, 'wb', template=bam)
    try:
        for read in bam.fetch(until_eof=True):
            stats['n_reads'] += 1
            if (read.is_secondary or read.is_supplementary or read.is_duplicate):
                continue
            stats['n_primary'] += 1
            if read.is_unmapped:
                continue
            stats['n_mapped'] += 1
            umi, side, run_len = classify_read(read, run_min, require_umi=require_umi)
            if side is None:
                continue
            if umi is not None:
                read.set_tag('UB', umi, value_type='Z')
                stats['n_pass_with_umi'] += 1
            read.set_tag('XP', side, value_type='A')
            read.set_tag('XR', run_len, value_type='i')
            out.write(read)
            stats['n_pass'] += 1
            stats['sum_run_len'] += run_len
            if side == 'L':
                stats['n_pass_left'] += 1
            else:
                stats['n_pass_right'] += 1
    finally:
        out.close()
        bam.close()
    pysam.index(out_bam)
    stats['mean_run_len'] = (stats['sum_run_len'] / stats['n_pass']) if stats['n_pass'] else 0.0
    return stats


def _append_stats(stats_tsv: str, stats: dict) -> None:
    cols = ['barcode', 'n_reads', 'n_primary', 'n_mapped',
            'n_pass', 'n_pass_left', 'n_pass_right', 'n_pass_with_umi',
            'mean_run_len']
    import csv
    is_new = not os.path.exists(stats_tsv)
    with open(stats_tsv, 'a', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t')
        if is_new:
            w.writeheader()
        w.writerow({k: stats[k] for k in cols})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--src-bam', required=True)
    ap.add_argument('--out-bam', required=True)
    ap.add_argument('--run-min', type=int, default=6,
                    help='Minimum polyA/T run length to accept (default 6)')
    ap.add_argument('--require-umi', dest='require_umi', action='store_true', default=True,
                    help='Only keep reads where a 9 bp UMI window is recoverable (default)')
    ap.add_argument('--no-require-umi', dest='require_umi', action='store_false',
                    help='Keep polyA/T-passing reads even if 9 bp UMI window is not recoverable; those reads have no UB tag')
    ap.add_argument('--stats-tsv', default=None,
                    help='Append per-cell stats as a row to this TSV (header auto)')
    args = ap.parse_args()

    Path(args.out_bam).parent.mkdir(parents=True, exist_ok=True)
    stats = filter_bam(args.src_bam, args.out_bam, args.run_min,
                       require_umi=args.require_umi)
    print(f'{stats["barcode"]}: '
          f'mapped={stats["n_mapped"]} '
          f'pass={stats["n_pass"]} (L={stats["n_pass_left"]} R={stats["n_pass_right"]}'
          f' umi={stats["n_pass_with_umi"]}) '
          f'mean_run_len={stats["mean_run_len"]:.2f}', flush=True)
    if args.stats_tsv:
        _append_stats(args.stats_tsv, stats)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
