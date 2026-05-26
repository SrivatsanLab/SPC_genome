#!/usr/bin/env python3
"""
Convert paired-end FASTQ files to unaligned SAM format with CB:Z (cell barcode) tags.
Extracts barcode from read headers (added by atrandi_demux.py) and adds as CB:Z tag.
"""

import gzip
import sys
import argparse
from pathlib import Path

def fastq_to_unaligned_sam(read1_path, read2_path, output_sam, single_end=False):
    """
    Convert FASTQ files to unaligned SAM with CB:Z tags.

    Args:
        read1_path: Path to R1 FASTQ file
        read2_path: Path to R2 FASTQ file (ignored if single_end=True)
        output_sam: Path to output SAM file
        single_end: If True, only emit R1 as an unpaired unmapped record (FLAG=4)
    """

    # Determine if files are gzipped
    r1_open = gzip.open if read1_path.endswith('.gz') else open

    if single_end:
        with r1_open(read1_path, 'rt') as r1_file, \
             open(output_sam, 'w') as sam_file:

            sam_file.write("@HD\tVN:1.5\tSO:unsorted\n")

            while True:
                r1_header = r1_file.readline().strip()
                r1_seq = r1_file.readline().strip()
                r1_file.readline()  # plus line
                r1_qual = r1_file.readline().strip()

                if not r1_header:
                    break

                r1_parts = r1_header[1:].split(':')
                read_name = r1_parts[0] if len(r1_parts) > 0 else r1_header[1:]
                barcode = r1_parts[-1] if len(r1_parts) > 1 else ""

                # FLAG=4: unmapped, unpaired
                sam_file.write(f"{read_name}\t4\t*\t0\t0\t*\t*\t0\t0\t{r1_seq}\t{r1_qual}\tCB:Z:{barcode}\n")
        return

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
    parser.add_argument('positional', nargs='+',
                        help='Paired-end: read1 read2 output. Single-end (with --single-end): read1 output')
    parser.add_argument('--single-end', action='store_true',
                        help='Single-end mode: emit only R1 as unpaired unmapped records')

    args = parser.parse_args()

    if args.single_end:
        if len(args.positional) != 2:
            parser.error('--single-end requires exactly: read1 output')
        read1, output = args.positional
        fastq_to_unaligned_sam(read1, None, output, single_end=True)
        print(f"Converted {read1} to {output} (single-end)")
    else:
        if len(args.positional) != 3:
            parser.error('paired-end requires exactly: read1 read2 output')
        read1, read2, output = args.positional
        fastq_to_unaligned_sam(read1, read2, output)
        print(f"Converted {read1} and {read2} to {output}")

if __name__ == '__main__':
    main()
