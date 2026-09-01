#!/usr/bin/env python3
"""Filter a BAM to reads containing a C. elegans SL1 or SL2 trans-splice sequence.

SL trans-splicing appends a spliced-leader sequence to the 5' end of
~70-80% of C. elegans mRNAs. Reads that captured the trans-splice
junction have SL sequence at the 5' end (or its reverse complement at
the 3' end, for reverse-mapped reads).

Detection uses 14-nt seeds from the 3' end of each SL (the portion most
often captured in short reads, adjacent to the mRNA), matched anywhere
in the read sequence. Exact-match only — 14 nt is specific enough
(chance rate ~1 per 3e8 bp, negligible next to 100k reads/cell).

    SL1 (22 nt):  GGTTTAATTACCCAAGTTTGAG   → seed: TACCCAAGTTTGAG
    SL2 (22 nt):  GGTTTTAACCCAGTTACTCAAG   → seed: CCCAGTTACTCAAG

Also checks reverse complements to catch minus-strand alignments.

Usage:
    filter_sl_reads.py <in.bam> <out.bam>
"""

import argparse
import sys
from pathlib import Path

import pysam


SL_SEEDS = (
    'TACCCAAGTTTGAG',   # SL1 3' seed
    'CTCAAACTTGGGTA',   # SL1 3' seed, revcomp
    'CCCAGTTACTCAAG',   # SL2 3' seed
    'CTTGAGTAACTGGG',   # SL2 3' seed, revcomp
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('in_bam',  type=Path)
    ap.add_argument('out_bam', type=Path)
    args = ap.parse_args()

    args.out_bam.parent.mkdir(parents=True, exist_ok=True)

    n_in = n_out = 0
    with pysam.AlignmentFile(str(args.in_bam), 'rb') as ib, \
         pysam.AlignmentFile(str(args.out_bam), 'wb', template=ib) as ob:
        for read in ib:
            n_in += 1
            seq = read.query_sequence
            if seq and any(s in seq for s in SL_SEEDS):
                ob.write(read)
                n_out += 1

    pysam.index(str(args.out_bam))
    print(f'{args.in_bam.name}: {n_in} in -> {n_out} out ({n_out/max(n_in,1):.4f})')
    return 0


if __name__ == '__main__':
    sys.exit(main())
