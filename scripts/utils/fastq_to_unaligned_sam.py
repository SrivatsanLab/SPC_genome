#!/usr/bin/env python3
"""
Convert paired-end FASTQ files to unaligned SAM format with CB:Z (cell barcode) tags.
Extracts barcode from read headers (added by atrandi_demux.py) and adds as CB:Z tag.
"""

import gzip
import sys
import argparse
from pathlib import Path

def fastq_to_unaligned_sam(read1_path, read2_path, output_sam):
    """
    Convert FASTQ files to unaligned SAM with CB:Z tags.

    Args:
        read1_path: Path to R1 FASTQ file
        read2_path: Path to R2 FASTQ file
        output_sam: Path to output SAM file
    """

    # Determine if files are gzipped
    r1_open = gzip.open if read1_path.endswith('.gz') else open
    r2_open = gzip.open if read2_path.endswith('.gz') else open

    with r1_open(read1_path, 'rt') as r1_file, \
         r2_open(read2_path, 'rt') as r2_file, \
         open(output_sam, 'w') as sam_file:

        # Write SAM header
        sam_file.write("@HD\tVN:1.5\tSO:unsorted\n")

        while True:
            # Read 4 lines from each FASTQ (1 read)
            r1_header = r1_file.readline().strip()
            r1_seq = r1_file.readline().strip()
            r1_plus = r1_file.readline().strip()
            r1_qual = r1_file.readline().strip()

            r2_header = r2_file.readline().strip()
            r2_seq = r2_file.readline().strip()
            r2_plus = r2_file.readline().strip()
            r2_qual = r2_file.readline().strip()

            # Check for end of file
            if not r1_header or not r2_header:
                break

            # Extract read name and barcode from header
            # Header format from atrandi_demux.py: @READNAME:BARCODE
            r1_parts = r1_header[1:].split(':')
            read_name = r1_parts[0] if len(r1_parts) > 0 else r1_header[1:]
            barcode = r1_parts[-1] if len(r1_parts) > 1 else ""

            # Write R1 as first in pair (FLAG=77: paired, unmapped, mate unmapped, first in pair)
            sam_file.write(f"{read_name}\t77\t*\t0\t0\t*\t*\t0\t0\t{r1_seq}\t{r1_qual}\tCB:Z:{barcode}\n")

            # Write R2 as second in pair (FLAG=141: paired, unmapped, mate unmapped, second in pair)
            sam_file.write(f"{read_name}\t141\t*\t0\t0\t*\t*\t0\t0\t{r2_seq}\t{r2_qual}\tCB:Z:{barcode}\n")

def main():
    parser = argparse.ArgumentParser(description='Convert FASTQ to unaligned SAM with CB:Z tags')
    parser.add_argument('read1', help='Read 1 FASTQ file')
    parser.add_argument('read2', help='Read 2 FASTQ file')
    parser.add_argument('output', help='Output SAM file')

    args = parser.parse_args()

    fastq_to_unaligned_sam(args.read1, args.read2, args.output)
    print(f"Converted {args.read1} and {args.read2} to {args.output}")

if __name__ == '__main__':
    main()
