#!/usr/bin/env python3
"""
Combine DNA and RNA read counts for CapGTA data.

Reads two CSV files (readcounts_dna.csv and readcounts_rna.csv) and combines them
into a single readcounts.csv file with columns: barcode, dna_reads, rna_reads, total_reads

Usage:
    python combine_readcounts_gta.py <dna_readcounts.csv> <rna_readcounts.csv> <output.csv>
"""

import sys
import csv
import argparse


def combine_readcounts(dna_file, rna_file, output_file):
    """
    Combine DNA and RNA read counts.

    Args:
        dna_file: Path to DNA readcounts CSV (barcode, read_count)
        rna_file: Path to RNA readcounts CSV (barcode, read_count)
        output_file: Path to output combined CSV
    """
    # Read DNA read counts
    dna_counts = {}
    with open(dna_file, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)  # Skip header
        for row in reader:
            if len(row) >= 2:
                barcode = row[0]
                count = int(row[1])
                dna_counts[barcode] = count

    # Read RNA read counts
    rna_counts = {}
    with open(rna_file, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)  # Skip header
        for row in reader:
            if len(row) >= 2:
                barcode = row[0]
                count = int(row[1])
                rna_counts[barcode] = count

    # Combine counts - include all barcodes from both files
    all_barcodes = set(dna_counts.keys()) | set(rna_counts.keys())

    # Write combined output
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        # Include 'read_count' column for compatibility with detect_cells.py
        writer.writerow(['barcode', 'dna_reads', 'rna_reads', 'total_reads', 'read_count'])

        # Sort barcodes for consistent output
        for barcode in sorted(all_barcodes):
            dna_reads = dna_counts.get(barcode, 0)
            rna_reads = rna_counts.get(barcode, 0)
            total_reads = dna_reads + rna_reads

            # read_count is same as total_reads (for compatibility with detect_cells.py)
            writer.writerow([barcode, dna_reads, rna_reads, total_reads, total_reads])

    print(f"Combined {len(all_barcodes)} barcodes", file=sys.stderr)
    print(f"  DNA barcodes: {len(dna_counts)}", file=sys.stderr)
    print(f"  RNA barcodes: {len(rna_counts)}", file=sys.stderr)
    print(f"Output written to: {output_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Combine DNA and RNA read counts for CapGTA data'
    )
    parser.add_argument('dna_file', help='DNA readcounts CSV file')
    parser.add_argument('rna_file', help='RNA readcounts CSV file')
    parser.add_argument('output_file', help='Output combined CSV file')

    args = parser.parse_args()

    combine_readcounts(args.dna_file, args.rna_file, args.output_file)


if __name__ == '__main__':
    main()
