#!/usr/bin/env python3
"""
Compile exonic enrichment metrics from individual cell files into a single CSV.

Reads exonic enrichment values from QC metrics directory and compiles them into
a single CSV file with columns:
  - barcode
  - dna_enrichment_pre (pre-rescue)
  - rna_enrichment_pre (pre-rescue)
  - dna_enrichment_post (post-rescue, if available)
  - rna_enrichment_post (post-rescue, if available)

Usage:
    python compile_exonic_enrichment.py <qc_metrics_dir> <output.csv>
"""

import os
import sys
import csv
import glob
import argparse


def extract_barcode_from_filename(filename):
    """
    Extract barcode from filename.

    Example: AAAA-BBBB-CCCC-DDDD_exonic_enrichment_dna.txt -> AAAA-BBBB-CCCC-DDDD
    """
    basename = os.path.basename(filename)
    # Remove _exonic_enrichment_*.txt suffix
    barcode = basename.split('_exonic_enrichment_')[0]
    return barcode


def read_enrichment_value(filepath):
    """
    Read single enrichment value from file.

    File should contain a single floating-point number.
    """
    try:
        with open(filepath, 'r') as f:
            value = float(f.read().strip())
        return value
    except (ValueError, IOError) as e:
        print(f"Warning: Could not read {filepath}: {e}", file=sys.stderr)
        return None


def compile_exonic_enrichment(qc_dir, output_file):
    """
    Compile exonic enrichment metrics from QC directory.

    Args:
        qc_dir: Directory containing *_exonic_enrichment_*.txt files
        output_file: Output CSV path
    """
    # Find all enrichment files
    dna_pre_files = glob.glob(os.path.join(qc_dir, '*_exonic_enrichment_dna.txt'))
    rna_pre_files = glob.glob(os.path.join(qc_dir, '*_exonic_enrichment_rna.txt'))
    dna_post_files = glob.glob(os.path.join(qc_dir, '*_exonic_enrichment_dna_postrescue.txt'))
    rna_post_files = glob.glob(os.path.join(qc_dir, '*_exonic_enrichment_rna_postrescue.txt'))

    print(f"Found {len(dna_pre_files)} DNA pre-rescue enrichment files", file=sys.stderr)
    print(f"Found {len(rna_pre_files)} RNA pre-rescue enrichment files", file=sys.stderr)
    print(f"Found {len(dna_post_files)} DNA post-rescue enrichment files", file=sys.stderr)
    print(f"Found {len(rna_post_files)} RNA post-rescue enrichment files", file=sys.stderr)

    # Build dict of enrichment values by barcode
    enrichment_data = {}

    # Process pre-rescue DNA
    for filepath in dna_pre_files:
        barcode = extract_barcode_from_filename(filepath)
        value = read_enrichment_value(filepath)
        if barcode not in enrichment_data:
            enrichment_data[barcode] = {}
        enrichment_data[barcode]['dna_pre'] = value

    # Process pre-rescue RNA
    for filepath in rna_pre_files:
        barcode = extract_barcode_from_filename(filepath)
        value = read_enrichment_value(filepath)
        if barcode not in enrichment_data:
            enrichment_data[barcode] = {}
        enrichment_data[barcode]['rna_pre'] = value

    # Process post-rescue DNA (if available)
    for filepath in dna_post_files:
        barcode = extract_barcode_from_filename(filepath)
        value = read_enrichment_value(filepath)
        if barcode not in enrichment_data:
            enrichment_data[barcode] = {}
        enrichment_data[barcode]['dna_post'] = value

    # Process post-rescue RNA (if available)
    for filepath in rna_post_files:
        barcode = extract_barcode_from_filename(filepath)
        value = read_enrichment_value(filepath)
        if barcode not in enrichment_data:
            enrichment_data[barcode] = {}
        enrichment_data[barcode]['rna_post'] = value

    # Write compiled output
    has_post_rescue = len(dna_post_files) > 0 or len(rna_post_files) > 0

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)

        # Write header
        if has_post_rescue:
            header = ['barcode', 'dna_enrichment_pre', 'rna_enrichment_pre',
                      'dna_enrichment_post', 'rna_enrichment_post']
        else:
            header = ['barcode', 'dna_enrichment_pre', 'rna_enrichment_pre']

        writer.writerow(header)

        # Write data rows (sorted by barcode)
        for barcode in sorted(enrichment_data.keys()):
            data = enrichment_data[barcode]

            if has_post_rescue:
                row = [
                    barcode,
                    data.get('dna_pre', ''),
                    data.get('rna_pre', ''),
                    data.get('dna_post', ''),
                    data.get('rna_post', '')
                ]
            else:
                row = [
                    barcode,
                    data.get('dna_pre', ''),
                    data.get('rna_pre', '')
                ]

            writer.writerow(row)

    print(f"Compiled {len(enrichment_data)} cells", file=sys.stderr)
    print(f"Output written to: {output_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Compile exonic enrichment metrics from individual cell files'
    )
    parser.add_argument('qc_dir', help='QC metrics directory containing enrichment files')
    parser.add_argument('output_file', help='Output CSV file')

    args = parser.parse_args()

    if not os.path.isdir(args.qc_dir):
        print(f"Error: QC directory not found: {args.qc_dir}", file=sys.stderr)
        sys.exit(1)

    compile_exonic_enrichment(args.qc_dir, args.output_file)


if __name__ == '__main__':
    main()
