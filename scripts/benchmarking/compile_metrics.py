#!/usr/bin/env python3
"""
Metrics Compiler for Downsampling Analysis

Purpose: Parse all Picard WGS metrics files into single CSV
Input:   Directory containing *_wgs_metrics.txt files
Output:  Compiled CSV with all metrics and metadata
"""

import argparse
import pandas as pd
import sys
from pathlib import Path
import re


def parse_picard_metrics(metrics_file):
    """
    Parse a Picard metrics file.

    Returns a dictionary with metric values.
    Picard files have:
    - Comment lines starting with #
    - A header line with column names
    - A data line with values
    """
    with open(metrics_file, 'r') as f:
        lines = f.readlines()

    # Find the header line (first non-comment line)
    header_idx = None
    for i, line in enumerate(lines):
        if not line.startswith('#') and line.strip():
            header_idx = i
            break

    if header_idx is None:
        raise ValueError(f"No header found in {metrics_file}")

    # Header is at header_idx, data is at header_idx + 1
    header = lines[header_idx].strip().split('\t')

    if header_idx + 1 >= len(lines):
        raise ValueError(f"No data line found in {metrics_file}")

    data = lines[header_idx + 1].strip().split('\t')

    if len(header) != len(data):
        raise ValueError(
            f"Header/data mismatch in {metrics_file}: "
            f"{len(header)} columns in header, {len(data)} in data"
        )

    # Create dictionary
    metrics = dict(zip(header, data))

    return metrics


def extract_metadata_from_filename(metrics_path):
    """
    Extract dataset, cell_id, and depth_label from metrics filename.

    Expected format: {dataset}/{cell_id}_{depth_label}_wgs_metrics.txt
    E.g., HSC4/AAAA-BBBB-CCCC-DDDD_5M_wgs_metrics.txt
    """
    # Get parent directory (dataset)
    dataset = metrics_path.parent.name

    # Parse filename: {cell_id}_{depth_label}_wgs_metrics.txt
    filename = metrics_path.name

    # Remove _wgs_metrics.txt suffix
    if not filename.endswith('_wgs_metrics.txt'):
        raise ValueError(f"Unexpected filename format: {filename}")

    base = filename[:-len('_wgs_metrics.txt')]

    # Split on last underscore to separate cell_id and depth_label
    # This assumes depth_label doesn't contain underscores (e.g., "5M", "10M")
    parts = base.rsplit('_', 1)

    if len(parts) != 2:
        raise ValueError(f"Could not parse cell_id and depth from: {filename}")

    cell_id, depth_label = parts

    # Convert depth_label to numeric target_depth
    depth_match = re.match(r'(\d+)([KM])', depth_label)
    if depth_match:
        num, unit = depth_match.groups()
        multiplier = 1_000 if unit == 'K' else 1_000_000
        target_depth = int(num) * multiplier
    else:
        # Assume it's already a number
        try:
            target_depth = int(depth_label)
        except ValueError:
            target_depth = None

    return {
        'dataset': dataset,
        'cell_id': cell_id,
        'depth_label': depth_label,
        'target_depth': target_depth
    }


def main():
    parser = argparse.ArgumentParser(
        description="Compile Picard WGS metrics from downsampling analysis"
    )
    parser.add_argument(
        "--metrics-dir",
        required=True,
        help="Directory containing metrics files (e.g., data/benchmarking/qc_metrics)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV file (e.g., results/benchmarking/compiled_metrics.csv)"
    )

    args = parser.parse_args()

    metrics_dir = Path(args.metrics_dir)
    output_file = Path(args.output)

    if not metrics_dir.exists():
        print(f"Error: Metrics directory not found: {metrics_dir}", file=sys.stderr)
        sys.exit(1)

    print("=" * 70)
    print("Metrics Compiler")
    print("=" * 70)
    print(f"Metrics directory: {metrics_dir}")
    print(f"Output file: {output_file}")
    print()

    # Find all metrics files
    metrics_files = list(metrics_dir.glob("**/*_wgs_metrics.txt"))

    if not metrics_files:
        print(f"Error: No metrics files found in {metrics_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(metrics_files)} metrics files")
    print()

    # Parse all metrics files
    all_metrics = []
    failed_count = 0

    for i, metrics_path in enumerate(sorted(metrics_files), 1):
        if i % 50 == 0:
            print(f"Processing file {i}/{len(metrics_files)}...")

        try:
            # Extract metadata from filename
            metadata = extract_metadata_from_filename(metrics_path)

            # Parse metrics file
            metrics = parse_picard_metrics(metrics_path)

            # Combine metadata and metrics
            combined = {**metadata, **metrics}
            all_metrics.append(combined)

        except Exception as e:
            print(f"Warning: Failed to parse {metrics_path}: {e}", file=sys.stderr)
            failed_count += 1
            continue

    if not all_metrics:
        print("Error: No metrics were successfully parsed", file=sys.stderr)
        sys.exit(1)

    # Create DataFrame
    df = pd.DataFrame(all_metrics)

    # Reorder columns: metadata first, then metrics
    metadata_cols = ['dataset', 'cell_id', 'depth_label', 'target_depth']
    metric_cols = [col for col in df.columns if col not in metadata_cols]
    df = df[metadata_cols + metric_cols]

    # Convert numeric columns
    numeric_cols = [
        'target_depth', 'GENOME_TERRITORY', 'MEAN_COVERAGE', 'SD_COVERAGE',
        'MEDIAN_COVERAGE', 'MAD_COVERAGE', 'PCT_EXC_ADAPTER', 'PCT_EXC_MAPQ',
        'PCT_EXC_DUPE', 'PCT_EXC_UNPAIRED', 'PCT_EXC_BASEQ', 'PCT_EXC_OVERLAP',
        'PCT_EXC_CAPPED', 'PCT_EXC_TOTAL', 'PCT_1X', 'PCT_5X', 'PCT_10X',
        'PCT_15X', 'PCT_20X', 'PCT_25X', 'PCT_30X', 'PCT_40X', 'PCT_50X',
        'PCT_60X', 'PCT_70X', 'PCT_80X', 'PCT_90X', 'PCT_100X'
    ]

    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # Sort by dataset, cell_id, target_depth
    df = df.sort_values(['dataset', 'cell_id', 'target_depth'])

    # Write output
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, index=False)

    # Summary
    print("=" * 70)
    print("Compilation complete!")
    print("=" * 70)
    print(f"Successfully parsed: {len(all_metrics)} files")
    print(f"Failed to parse: {failed_count} files")
    print()
    print("Metrics per dataset:")
    print(df['dataset'].value_counts())
    print()
    print("Metrics per depth:")
    print(df['depth_label'].value_counts().sort_index())
    print()
    print(f"Output written to: {output_file}")
    print(f"Total rows: {len(df)}")
    print(f"Total columns: {len(df.columns)}")
    print()

    # Show column names
    print("Available metrics:")
    for col in df.columns:
        print(f"  - {col}")
    print()


if __name__ == "__main__":
    main()
