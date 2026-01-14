#!/usr/bin/env python
"""
Compile QC metrics from Picard output files into summary tables

Usage:
    python scripts/compile_qc_metrics.py <qc_metrics_dir> <output_csv>

Example:
    python scripts/compile_qc_metrics.py data/benchmarking/qc_metrics/PTA PTA_qc_summary.csv
"""

import sys
import os
import pandas as pd
from pathlib import Path
import re

def parse_picard_metrics(filepath):
    """Parse a Picard metrics file and return metrics as a dictionary"""
    metrics = {}

    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Find the metrics section (starts after a line with column headers)
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith('## METRICS CLASS'):
            header_idx = i + 1
            break

    if header_idx is None:
        return metrics

    # Parse header and data
    header = lines[header_idx].strip().split('\t')
    data = lines[header_idx + 1].strip().split('\t')

    # Create dictionary from header and data
    for h, d in zip(header, data):
        try:
            # Try to convert to numeric
            metrics[h] = float(d) if '.' in d else int(d)
        except ValueError:
            metrics[h] = d

    return metrics

def compile_qc_metrics(qc_dir):
    """Compile all QC metrics from a directory into a single DataFrame"""
    qc_path = Path(qc_dir)

    if not qc_path.exists():
        print(f"ERROR: Directory not found: {qc_dir}")
        sys.exit(1)

    # Find all sample names (from alignment metrics files)
    alignment_files = list(qc_path.glob("*_alignment_metrics.txt"))

    if not alignment_files:
        print(f"ERROR: No alignment metrics files found in {qc_dir}")
        sys.exit(1)

    print(f"Found {len(alignment_files)} samples in {qc_dir}")

    results = []

    for alignment_file in sorted(alignment_files):
        # Extract sample name
        sample = alignment_file.stem.replace('_alignment_metrics', '')

        # Parse all metrics files for this sample
        sample_metrics = {'sample': sample}

        # Alignment metrics
        align_file = qc_path / f"{sample}_alignment_metrics.txt"
        if align_file.exists():
            align_data = parse_picard_metrics(align_file)
            # Extract key alignment metrics
            if 'PCT_PF_READS_ALIGNED' in align_data:
                sample_metrics['pct_reads_aligned'] = align_data['PCT_PF_READS_ALIGNED']
            if 'MEAN_READ_LENGTH' in align_data:
                sample_metrics['mean_read_length'] = align_data['MEAN_READ_LENGTH']
            if 'PF_READS_ALIGNED' in align_data:
                sample_metrics['reads_aligned'] = align_data['PF_READS_ALIGNED']
            if 'PF_READS' in align_data:
                sample_metrics['total_reads'] = align_data['PF_READS']

        # GC metrics
        gc_summary_file = qc_path / f"{sample}_gc_summary.txt"
        if gc_summary_file.exists():
            gc_data = parse_picard_metrics(gc_summary_file)
            # Extract GC content
            if 'MEAN_GC' in gc_data:
                sample_metrics['mean_gc'] = gc_data['MEAN_GC']
            if 'AT_DROPOUT' in gc_data:
                sample_metrics['at_dropout'] = gc_data['AT_DROPOUT']
            if 'GC_DROPOUT' in gc_data:
                sample_metrics['gc_dropout'] = gc_data['GC_DROPOUT']

        # Duplicate metrics
        dup_file = qc_path / f"{sample}_duplicate_metrics.txt"
        if dup_file.exists():
            dup_data = parse_picard_metrics(dup_file)
            if 'PERCENT_DUPLICATION' in dup_data:
                sample_metrics['pct_duplicates'] = dup_data['PERCENT_DUPLICATION']
            if 'READ_PAIRS_EXAMINED' in dup_data:
                sample_metrics['read_pairs_examined'] = dup_data['READ_PAIRS_EXAMINED']
            if 'READ_PAIR_DUPLICATES' in dup_data:
                sample_metrics['read_pair_duplicates'] = dup_data['READ_PAIR_DUPLICATES']

        # WGS metrics
        wgs_file = qc_path / f"{sample}_wgs_metrics.txt"
        if wgs_file.exists():
            wgs_data = parse_picard_metrics(wgs_file)
            if 'MEAN_COVERAGE' in wgs_data:
                sample_metrics['mean_coverage'] = wgs_data['MEAN_COVERAGE']
            if 'MEDIAN_COVERAGE' in wgs_data:
                sample_metrics['median_coverage'] = wgs_data['MEDIAN_COVERAGE']
            if 'PCT_1X' in wgs_data:
                sample_metrics['pct_1x'] = wgs_data['PCT_1X']
            if 'PCT_10X' in wgs_data:
                sample_metrics['pct_10x'] = wgs_data['PCT_10X']

        results.append(sample_metrics)

    # Convert to DataFrame
    df = pd.DataFrame(results)

    return df

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python compile_qc_metrics.py <qc_metrics_dir> <output_csv>")
        print("")
        print("Example:")
        print("  python scripts/compile_qc_metrics.py data/benchmarking/qc_metrics/PTA PTA_qc_summary.csv")
        sys.exit(1)

    qc_dir = sys.argv[1]
    output_csv = sys.argv[2]

    print(f"Compiling QC metrics from: {qc_dir}")

    df = compile_qc_metrics(qc_dir)

    # Save to CSV
    df.to_csv(output_csv, index=False)

    print(f"Saved summary to: {output_csv}")
    print(f"Total samples: {len(df)}")
    print("")
    print("Summary statistics:")
    print(df.describe())
