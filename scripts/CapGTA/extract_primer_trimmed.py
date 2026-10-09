#!/usr/bin/env python3
"""Extract primer-captured read pairs from a single cell's R1/R2 FASTQs.

The RT primer leaves a [10 bp UMI][polyT] signature at R1's 5' end on read
pairs that captured the original mRNA 3' end. These pairs fail paired-end
alignment because STAR either rejects them (too-short alignable region) or
force-aligns the polyT into genomic T-runs (with no masking in the reference,
there are plenty). Pre-trimming the primer out of R1 and storing the UMI as a
header tag lets both mates align normally and preserves UMI for dedup.

For each pair we scan R1:

  - polyT run anchored at the LEFT edge, starting at some position in
    [0, max_start] (default [0, 20]). Length must be >= run_min (default 15).
    Up to max_mm (default 1) HIGH-quality interior mismatches are allowed.
    Non-T bases with quality <= qual_threshold (default 20) are treated as
    ambiguous: they count as part of the run and don't consume the mm budget.
    Trailing true (high-Q) mismatches are trimmed off.

  - UMI window of 9 or 10 bp immediately before the polyT run:
      * 10 bp if polyT starts at position >=10
      * 9 bp if polyT starts at position 9
    (A 9 bp window keeps reads that lost one UMI base to the primer boundary.)
    Must be all ACGT, and must not be itself AT-rich (>=7/9 or >=8/10 A or T),
    which would indicate the UMI window has been contaminated by polyT bleed.

If both pass, the pair is kept with:
  - R1 trimmed to remove [UMI][polyT] (just the cDNA tail)
  - R2 unchanged
  - UB:Z:<umi> appended to both R1 and R2 headers (existing CB:Z: tag preserved)
  - Both R1 cDNA tail and R2 must be >= min_cdna_len (default 20) after trim
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# PolyT detection with quality-aware mismatch budget
# ---------------------------------------------------------------------------

def find_homopolymer_left_edge(seq: str, qual: str, base: str,
                               run_min: int = 15, max_start: int = 20,
                               max_mm: int = 1, qual_threshold: int = 20):
    """Find longest homopolymer run of `base` ('T' or 'A') starting in [0, max_start].

    Quality-aware semantics (same as before, parameterised on `base`):
      - Non-`base` bases with qual <= qual_threshold are ambiguous → counted
        in the run as if they were `base`, and do NOT consume the mm budget.
      - Non-`base` bases with qual > qual_threshold are TRUE mismatches:
        counted in run length, consuming budget. The (max_mm+1)th breaks it.
      - Trailing true mismatches are stripped from the run length.

    Returns (start, run_length) or (None, 0).
    """
    if not seq or len(seq) != len(qual):
        return (None, 0)

    best_start = None
    best_len = 0

    for start in range(min(max_start + 1, len(seq))):
        q0 = ord(qual[start]) - 33
        if seq[start] != base and q0 > qual_threshold:
            continue

        mm = 0
        last_contrib_idx = -1
        for i in range(start, len(seq)):
            ch = seq[i]
            q = ord(qual[i]) - 33
            if ch == base:
                last_contrib_idx = i
            elif q <= qual_threshold:
                last_contrib_idx = i
            else:
                mm += 1
                if mm > max_mm:
                    break

        if last_contrib_idx < start:
            run_len = 0
        else:
            run_len = last_contrib_idx - start + 1

        if run_len > best_len:
            best_len = run_len
            best_start = start

    if best_len >= run_min:
        return (best_start, best_len)
    return (None, 0)


# Keep old name as a thin wrapper for backwards compatibility (and the self-test).
def find_polyT_left_edge(seq: str, qual: str, run_min: int = 15,
                         max_start: int = 20, max_mm: int = 1,
                         qual_threshold: int = 20):
    return find_homopolymer_left_edge(seq, qual, 'T', run_min=run_min,
                                      max_start=max_start, max_mm=max_mm,
                                      qual_threshold=qual_threshold)


_RC_TABLE = str.maketrans('ACGTN', 'TGCAN')


def revcomp(s: str) -> str:
    return s.translate(_RC_TABLE)[::-1]


def umi_window_ok(window: str) -> bool:
    """Window must be 9 or 10 bp, all ACGT, and not itself AT-rich.

    AT-rich threshold scales with length: >=7 of 9 or >=8 of 10 A or T is rejected.
    """
    L = len(window)
    if L not in (9, 10):
        return False
    if any(b not in 'ACGT' for b in window):
        return False
    at_thresh = 7 if L == 9 else 8
    if window.count('T') >= at_thresh or window.count('A') >= at_thresh:
        return False
    return True


def extract_umi(seq: str, polyT_start: int) -> str:
    """UMI = 9 or 10 bp immediately before polyT (whichever fits)."""
    if polyT_start >= 10:
        return seq[polyT_start - 10:polyT_start]
    elif polyT_start == 9:
        return seq[0:9]
    return ''  # polyT_start < 9: no viable UMI window


# ---------------------------------------------------------------------------
# FASTQ IO
# ---------------------------------------------------------------------------

class FastqCorrupt(Exception):
    """Raised when an input FASTQ record fails structural validation."""


def _validate_record(lines, label, rec_num):
    """Raise FastqCorrupt if the 4-line block is not a valid FASTQ record."""
    h, s, p, q = lines
    if not h.startswith('@'):
        raise FastqCorrupt(f'{label} record {rec_num}: header does not start with @ '
                           f'(got {h[:40]!r})')
    if not p.startswith('+'):
        raise FastqCorrupt(f'{label} record {rec_num}: separator line not + '
                           f'(got {p[:40]!r})')
    if len(s) != len(q):
        raise FastqCorrupt(f'{label} record {rec_num}: '
                           f'seq_len={len(s)} != qual_len={len(q)}')


def fastq_pair_iter(r1_path: str, r2_path: str):
    """Yield (r1_header, r1_seq, r1_qual, r2_header, r2_seq, r2_qual) tuples.

    Validates each 4-line record on both sides (header starts with @, separator
    starts with +, len(seq) == len(qual)). Raises FastqCorrupt the moment a bad
    record is seen, so callers can abort before writing out garbage. A truncated
    tail (one side ends with only partial lines of the final record) is also
    treated as corruption.
    """
    rec_num = 0
    with gzip.open(r1_path, 'rt') as f1, gzip.open(r2_path, 'rt') as f2:
        while True:
            r1_lines = [f1.readline() for _ in range(4)]
            r2_lines = [f2.readline() for _ in range(4)]
            # Clean EOF on both sides.
            if not r1_lines[0] and not r2_lines[0]:
                break
            rec_num += 1
            # Partial record on either side → corruption.
            if any(l == '' for l in r1_lines):
                raise FastqCorrupt(f'R1 record {rec_num}: truncated (partial record at EOF)')
            if any(l == '' for l in r2_lines):
                raise FastqCorrupt(f'R2 record {rec_num}: truncated (partial record at EOF)')
            r1_stripped = [l.rstrip('\n') for l in r1_lines]
            r2_stripped = [l.rstrip('\n') for l in r2_lines]
            _validate_record(r1_stripped, 'R1', rec_num)
            _validate_record(r2_stripped, 'R2', rec_num)
            yield (r1_stripped[0], r1_stripped[1], r1_stripped[3],
                   r2_stripped[0], r2_stripped[1], r2_stripped[3])


def append_umi_tag(header: str, umi: str) -> str:
    """Append UB:Z:<umi> to a FASTQ header line (preserving any existing tags)."""
    return f'{header} UB:Z:{umi}'


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

_CASES_BIDIR = [
    # (code, which read, which base)
    ('R1T', 'R1', 'T'),
    ('R2T', 'R2', 'T'),
    ('R1A', 'R1', 'A'),
    ('R2A', 'R2', 'A'),
]
_CASES_R1T_ONLY = [('R1T', 'R1', 'T')]


def _best_signature(r1s, r1q, r2s, r2q, cases, run_min, max_mm, qual_threshold, max_start):
    """Try each candidate case; return the best (longest run, earliest tie-breaker) or None."""
    best = None  # (run_len, case_rank, code, which, base, start)
    for rank, (code, which, base) in enumerate(cases):
        seq, qual = (r1s, r1q) if which == 'R1' else (r2s, r2q)
        start, run_len = find_homopolymer_left_edge(
            seq, qual, base, run_min=run_min, max_start=max_start,
            max_mm=max_mm, qual_threshold=qual_threshold)
        if start is None:
            continue
        key = (run_len, -rank, code, which, base, start)
        if best is None or key > best:
            best = key
    return best


def process(r1_path: str, r2_path: str, out_r1_path: str, out_r2_path: str,
            run_min: int, max_mm: int, qual_threshold: int, min_cdna_len: int,
            max_start: int = 20, require_umi: bool = True,
            any_orientation: bool = False) -> dict:
    """Scan pairs for a left-edge polyT/polyA signature on R1 or R2.

    Modes:
      - any_orientation=False (default): only R1 5' polyT (original behaviour).
      - any_orientation=True: try R1-polyT, R2-polyT, R1-polyA, R2-polyA; pick
        the one with the longest run (R1-polyT > R2-polyT > R1-polyA > R2-polyA
        on ties). The read bearing the signature is trimmed; the other mate
        is written unchanged. UMI from polyA cases is reverse-complemented so
        all UMIs share a canonical orientation.

    When require_umi=True, pairs with a signature but a failing UMI window are
    DROPPED. When False, those pairs are KEPT with no UB:Z: tag.
    """
    cases = _CASES_BIDIR if any_orientation else _CASES_R1T_ONLY

    n_total = 0
    n_no_sig = 0
    n_bad_umi = 0
    n_too_short = 0
    n_kept_with_umi = 0
    n_kept_no_umi = 0
    run_len_dist: dict = {}
    start_dist: dict = {}
    umi_len_dist: dict = {}
    case_dist: dict = {}

    Path(out_r1_path).parent.mkdir(parents=True, exist_ok=True)

    corrupt_at = None
    corrupt_msg = None
    with gzip.open(out_r1_path, 'wt', compresslevel=5) as o1, \
         gzip.open(out_r2_path, 'wt', compresslevel=5) as o2:
      try:
        for r1h, r1s, r1q, r2h, r2s, r2q in fastq_pair_iter(r1_path, r2_path):
            n_total += 1
            best = _best_signature(r1s, r1q, r2s, r2q, cases,
                                   run_min, max_mm, qual_threshold, max_start)
            if best is None:
                n_no_sig += 1
                continue
            run_len, _rank_neg, code, which, base, start = best

            trim_seq, trim_qual = (r1s, r1q) if which == 'R1' else (r2s, r2q)
            umi_raw = extract_umi(trim_seq, start)
            umi_valid = umi_window_ok(umi_raw)
            if not umi_valid and require_umi:
                n_bad_umi += 1
                continue

            # Trim the read that carries the signature; the other is kept intact.
            trim_end = start + run_len
            trimmed_s = trim_seq[trim_end:]
            trimmed_q = trim_qual[trim_end:]
            if which == 'R1':
                r1_out_s, r1_out_q = trimmed_s, trimmed_q
                r2_out_s, r2_out_q = r2s, r2q
            else:
                r1_out_s, r1_out_q = r1s, r1q
                r2_out_s, r2_out_q = trimmed_s, trimmed_q

            if len(r1_out_s) < min_cdna_len or len(r2_out_s) < min_cdna_len:
                n_too_short += 1
                continue

            # Canonicalise UMI across polyT/polyA cases: polyA → RC.
            umi_canon = revcomp(umi_raw) if (base == 'A' and umi_raw) else umi_raw
            # Build headers: always include an XP:Z:<code> tag; UMI only when valid.
            def _tag_header(h):
                parts = [h]
                if umi_valid:
                    parts.append(f'UB:Z:{umi_canon}')
                parts.append(f'XP:Z:{code}')
                return ' '.join(parts)
            new_r1h = _tag_header(r1h)
            new_r2h = _tag_header(r2h)

            o1.write(f'{new_r1h}\n{r1_out_s}\n+\n{r1_out_q}\n')
            o2.write(f'{new_r2h}\n{r2_out_s}\n+\n{r2_out_q}\n')

            if umi_valid:
                n_kept_with_umi += 1
                umi_len_dist[len(umi_canon)] = umi_len_dist.get(len(umi_canon), 0) + 1
            else:
                n_kept_no_umi += 1
            run_len_dist[run_len] = run_len_dist.get(run_len, 0) + 1
            start_dist[start] = start_dist.get(start, 0) + 1
            case_dist[code] = case_dist.get(code, 0) + 1
      except FastqCorrupt as e:
        corrupt_at = n_total + 1  # record that failed validation
        corrupt_msg = str(e)
        # Writers already have the records before the corruption flushed via with.

    n_kept = n_kept_with_umi + n_kept_no_umi
    return dict(
        r1_path=r1_path,
        r2_path=r2_path,
        run_min=run_min,
        max_mm=max_mm,
        qual_threshold=qual_threshold,
        min_cdna_len=min_cdna_len,
        max_start=max_start,
        require_umi=require_umi,
        any_orientation=any_orientation,
        total_pairs=n_total,
        pairs_no_signature=n_no_sig,
        pairs_bad_umi=n_bad_umi,
        pairs_dropped_short_after_trim=n_too_short,
        pairs_kept=n_kept,
        pairs_kept_with_umi=n_kept_with_umi,
        pairs_kept_no_umi=n_kept_no_umi,
        frac_kept=(n_kept / n_total) if n_total else 0.0,
        case_dist={k: case_dist.get(k, 0) for k in ('R1T', 'R2T', 'R1A', 'R2A')},
        run_len_dist={str(k): v for k, v in sorted(run_len_dist.items())},
        start_dist={str(k): v for k, v in sorted(start_dist.items())},
        umi_len_dist={str(k): v for k, v in sorted(umi_len_dist.items())},
        input_corrupt=corrupt_at is not None,
        corrupt_at_record=corrupt_at,
        corrupt_msg=corrupt_msg,
    )


# ---------------------------------------------------------------------------
# Sanity self-test (--self-test)
# ---------------------------------------------------------------------------

def _self_test() -> int:
    """Hand-constructed cases that must pass with default thresholds."""
    def hq(seq):
        return 'J' * len(seq)       # Q41 for every base — all mismatches are "true"
    # Use UMI stems that contain NO T, so no spurious polyT runs start within the UMI.
    UMI10 = 'AGCCAGCGAC'   # 10 bp, no T, balanced
    UMI9  = 'AGCCAGCGA'    # 9 bp, no T

    # case 1: canonical UMI (10 bp, no T) + 20 T + cDNA
    s = UMI10 + 'T' * 20 + 'ACGACGACGA'
    r = find_polyT_left_edge(s, hq(s))
    assert r == (10, 20), f'case 1 got {r}'
    # case 2: polyT starts at pos 9 (UMI = 9 bp, no T)
    s = UMI9 + 'T' * 20 + 'ACGACG'
    r = find_polyT_left_edge(s, hq(s))
    assert r == (9, 20), f'case 2 got {r}'
    # case 3: 1 interior high-Q mismatch → still accepted (run length includes the mm)
    s = UMI10 + 'T' * 10 + 'G' + 'T' * 9
    r = find_polyT_left_edge(s, hq(s))
    assert r == (10, 20), f'case 3 got {r}'
    # case 4: 2 interior high-Q mismatches → run truncates at the 2nd (use run_min=10 to accept)
    s = UMI10 + 'TTTTT' + 'G' + 'TTTT' + 'G' + 'TTT'
    r = find_polyT_left_edge(s, hq(s), run_min=10)
    # From start=10: TTTTT (last=14) + G (mm=1) + TTTT (last=19)
    #   next G (mm=2) → break. run_len = 19-10+1 = 10.
    assert r == (10, 10), f'case 4 got {r}'
    # case 5: low-Q non-T is FREE
    s = UMI10 + 'TTTTTTTTTT' + 'C' + 'TTTTTTTTT'
    q = 'J' * 20 + '!' + 'J' * (len(s) - 21)   # pos 20 is Q0 (ambiguous)
    r = find_polyT_left_edge(s, q)
    assert r == (10, 20), f'case 5 got {r}'
    # case 6: polyT too short → reject
    s = UMI10 + 'T' * 14 + 'ACGA'
    assert find_polyT_left_edge(s, hq(s), run_min=15) == (None, 0)
    # case 7: no polyT at all
    s = UMI10 + 'ACGACGACGACGACGACGACGACGA'
    assert find_polyT_left_edge(s, hq(s)) == (None, 0)
    # case 8: polyT starts after position 20 → not anchored to left edge
    s = 'A' * 25 + 'T' * 20 + 'AC'
    assert find_polyT_left_edge(s, hq(s)) == (None, 0)
    # case 9: UMI window check
    assert umi_window_ok('AGCCTGTTAC') is True
    assert umi_window_ok('AGCCTGTTA')  is True
    assert umi_window_ok('AAAAAAAAAC') is False     # 9A/10 → reject
    assert umi_window_ok('TTTTTTTTT')  is False     # 9T/9 → reject
    assert umi_window_ok('AGCCTGTTACN') is False    # 11 bp → reject
    assert umi_window_ok('AGCCTGN')    is False     # 7 bp → reject
    assert umi_window_ok('NGCCTGTTAC') is False     # has N → reject
    # case 10: UMI extraction
    assert extract_umi('AGCCTGTTAC' + 'T' * 20, 10) == 'AGCCTGTTAC'
    assert extract_umi('AGCCTGTTA'  + 'T' * 20, 9)  == 'AGCCTGTTA'
    assert extract_umi('AGCCTGTT'   + 'T' * 20, 8)  == ''  # polyT_start < 9 → no window
    # case 11: homopolymer polyA path is symmetric
    # UMI stem with no A (and no T) to avoid confounding runs starting inside the stem.
    UMI10_noA = 'GCCTGCTCGC'
    s = UMI10_noA + 'A' * 20 + 'CGCCGCCGCC'   # pretend reverse-orientation read
    r = find_homopolymer_left_edge(s, hq(s), 'A')
    assert r == (10, 20), f'case 11 got {r}'
    # case 12: reverse-complement helper
    assert revcomp('ACGT') == 'ACGT'
    assert revcomp('AAAA') == 'TTTT'
    assert revcomp('AGCCTGTTAC') == 'GTAACAGGCT'
    print('self-test: all cases passed')
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--self-test', action='store_true',
                    help='Run built-in sanity tests and exit')
    ap.add_argument('--r1', help='Input R1 fastq.gz')
    ap.add_argument('--r2', help='Input R2 fastq.gz')
    ap.add_argument('--out-r1', help='Output primer-trimmed R1 fastq.gz')
    ap.add_argument('--out-r2', help='Output R2 fastq.gz (UMI-tagged header)')
    ap.add_argument('--out-stats', help='Per-cell stats JSON')
    ap.add_argument('--run-min', type=int, default=15,
                    help='Minimum polyT run length to accept (default 15)')
    ap.add_argument('--max-mm', type=int, default=1,
                    help='Max HIGH-Q interior mismatches in polyT run (default 1)')
    ap.add_argument('--qual-threshold', type=int, default=20,
                    help='Non-T bases with Phred <= this are treated as T (default 20)')
    ap.add_argument('--max-start', type=int, default=20,
                    help='PolyT must start at position <= this in R1 (default 20)')
    ap.add_argument('--min-cdna-len', type=int, default=20,
                    help='After trim, both R1 cDNA and R2 must be >= this bp (default 20)')
    ap.add_argument('--require-umi', dest='require_umi', action='store_true', default=True,
                    help='Drop pairs whose UMI window is missing/invalid (default)')
    ap.add_argument('--no-require-umi', dest='require_umi', action='store_false',
                    help='Keep polyT pairs even when the UMI window is missing or invalid; '
                         'those pairs get no UB:Z: tag')
    ap.add_argument('--any-orientation', dest='any_orientation', action='store_true',
                    help='Scan R1 AND R2, polyT AND polyA. Pick the longest matching '
                         'homopolymer across the four cases; trim the read bearing it; '
                         'reverse-complement UMI for polyA cases to canonicalise.')
    args = ap.parse_args()

    if args.self_test:
        return _self_test()

    required = ['r1', 'r2', 'out_r1', 'out_r2', 'out_stats']
    missing = [f'--{k.replace("_", "-")}' for k in required if getattr(args, k) is None]
    if missing:
        ap.error(f'required args missing: {" ".join(missing)}')

    stats = process(args.r1, args.r2, args.out_r1, args.out_r2,
                    run_min=args.run_min, max_mm=args.max_mm,
                    qual_threshold=args.qual_threshold,
                    min_cdna_len=args.min_cdna_len, max_start=args.max_start,
                    require_umi=args.require_umi,
                    any_orientation=args.any_orientation)
    Path(args.out_stats).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_stats, 'w') as fh:
        json.dump(stats, fh, indent=2)
    tag = ''
    if stats.get('input_corrupt'):
        tag = f'  [INPUT CORRUPT @ rec {stats["corrupt_at_record"]}: {stats["corrupt_msg"]}]'
    print(f'{Path(args.r1).name}: kept {stats["pairs_kept"]}/{stats["total_pairs"]} '
          f'({100*stats["frac_kept"]:.2f}%){tag}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
