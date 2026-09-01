#!/usr/bin/env python3
"""Filter a BAM to reads that are EITHER spliced (CIGAR contains 'N') OR contain a run of >=10 A or 10 T.

Union of the two RNA-signal filters used in CapGTA_dev diagnostics —
lets us build a transcript-count matrix that includes both spliced and
polyA/T-carrying reads as RNA evidence in one pass over the BAM.

Usage:
    filter_spliced_polyat_reads.py <in.bam> <out.bam> [--min-run 10]

Also indexes the output.
"""

import argparse
import re
import sys
from pathlib import Path

import pysam


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('in_bam',  type=Path)
    ap.add_argument('out_bam', type=Path)
    ap.add_argument('--min-run', type=int, default=10)
    args = ap.parse_args()

    polyat = re.compile(rf'A{{{args.min_run},}}|T{{{args.min_run},}}')

    args.out_bam.parent.mkdir(parents=True, exist_ok=True)

    n_in = n_spliced = n_polyat = n_out = 0
    with pysam.AlignmentFile(str(args.in_bam), 'rb') as ib, \
         pysam.AlignmentFile(str(args.out_bam), 'wb', template=ib) as ob:
        for read in ib:
            n_in += 1
            cigar = read.cigarstring
            seq   = read.query_sequence
            has_junc  = bool(cigar) and ('N' in cigar)
            has_polyat = bool(seq) and bool(polyat.search(seq))
            if has_junc or has_polyat:
                ob.write(read)
                n_out += 1
                if has_junc:   n_spliced += 1
                if has_polyat: n_polyat  += 1

    pysam.index(str(args.out_bam))
    print(f'{args.in_bam.name}: {n_in} in -> {n_out} out '
          f'(spliced={n_spliced}, polyat={n_polyat}, union rate={n_out/max(n_in,1):.4f})')
    return 0


if __name__ == '__main__':
    sys.exit(main())
