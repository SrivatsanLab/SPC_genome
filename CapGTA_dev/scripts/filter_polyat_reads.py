#!/usr/bin/env python3
"""Filter a BAM to reads containing a run of >= 10 A or 10 T in the sequence.

polyT-primed cDNA in this assay → RNA-derived fragments carry the polyA/T
signature (either in aligned or soft-clipped bases). We treat any read
with >= 10 consecutive A or T anywhere in its sequence as a candidate
RNA read.

Usage:
    filter_polyat_reads.py <in.bam> <out.bam> [--min-run 10]

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
    ap.add_argument('in_bam', type=Path)
    ap.add_argument('out_bam', type=Path)
    ap.add_argument('--min-run', type=int, default=10)
    args = ap.parse_args()

    pattern = re.compile(rf'A{{{args.min_run},}}|T{{{args.min_run},}}')

    args.out_bam.parent.mkdir(parents=True, exist_ok=True)

    n_in = n_out = 0
    with pysam.AlignmentFile(str(args.in_bam), 'rb') as ib, \
         pysam.AlignmentFile(str(args.out_bam), 'wb', template=ib) as ob:
        for read in ib:
            n_in += 1
            seq = read.query_sequence
            if seq and pattern.search(seq):
                ob.write(read)
                n_out += 1

    pysam.index(str(args.out_bam))
    print(f'{args.in_bam.name}: {n_in} in -> {n_out} out ({n_out/max(n_in,1):.4f})')
    return 0


if __name__ == '__main__':
    sys.exit(main())
