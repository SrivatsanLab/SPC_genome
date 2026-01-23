#!/usr/bin/env python3
"""
Parse Picard WGS metrics to extract coverage statistics for normalization.

Extracts:
- Mean coverage
- Callable megabases at different thresholds (5x, 10x, 20x, 30x)
"""

import pandas as pd
import argparse
from pathlib import Path


def parse_wgs_metrics(metrics_file):
    """Parse Picard CollectWgsMetrics output file.

    Parameters
    ----------
    metrics_file : str
        Path to Picard WGS metrics file

    Returns
    -------
    dict
        Dictionary with coverage statistics
    """
    with open(metrics_file, 'r') as f:
        lines = f.readlines()

    # Find the metrics section (starts after ## METRICS CLASS line)
    metrics_start = None
    for i, line in enumerate(lines):
        if line.startswith('## METRICS CLASS'):
            metrics_start = i + 1
            break

    if metrics_start is None:
        raise ValueError(f"Could not find metrics section in {metrics_file}")

    # Read header and values
    header = lines[metrics_start].strip().split('\t')
    values = lines[metrics_start + 1].strip().split('\t')

    # Create dict
    metrics = dict(zip(header, values))

    # Extract key metrics
    result = {
        'mean_coverage': float(metrics['MEAN_COVERAGE']),
        'median_coverage': float(metrics['MEDIAN_COVERAGE']),
        'pct_5x': float(metrics.get('PCT_5X', 0.0)),  # Use get() in case not present
        'pct_10x': float(metrics['PCT_10X']),
        'pct_20x': float(metrics['PCT_20X']),
        'pct_30x': float(metrics['PCT_30X']),
        'genome_territory': int(metrics['GENOME_TERRITORY']),
        'sd_coverage': float(metrics['SD_COVERAGE'])
    }

    # Calculate callable bases at different thresholds
    result['callable_bases_5x'] = int(result['genome_territory'] * result['pct_5x'])
    result['callable_bases_10x'] = int(result['genome_territory'] * result['pct_10x'])
    result['callable_bases_20x'] = int(result['genome_territory'] * result['pct_20x'])
    result['callable_bases_30x'] = int(result['genome_territory'] * result['pct_30x'])

    # Calculate callable megabases
    result['callable_Mb_5x'] = result['callable_bases_5x'] / 1e6
    result['callable_Mb_10x'] = result['callable_bases_10x'] / 1e6
    result['callable_Mb_20x'] = result['callable_bases_20x'] / 1e6
    result['callable_Mb_30x'] = result['callable_bases_30x'] / 1e6

    return result


def main():
    parser = argparse.ArgumentParser(
        description='Parse Picard WGS metrics for coverage normalization'
    )
    parser.add_argument(
        'metrics_dir',
        help='Directory containing WGS metrics files (*_wgs_metrics.txt)'
    )
    parser.add_argument(
        '-o', '--output',
        default='coverage_summary.tsv',
        help='Output TSV file (default: coverage_summary.tsv)'
    )

    args = parser.parse_args()

    # Find all WGS metrics files
    metrics_dir = Path(args.metrics_dir)
    metrics_files = sorted(metrics_dir.glob('*_wgs_metrics.txt'))

    if not metrics_files:
        print(f"No WGS metrics files found in {metrics_dir}")
        return

    print(f"Found {len(metrics_files)} WGS metrics files")

    # Parse all files
    results = []
    for metrics_file in metrics_files:
        sample_name = metrics_file.stem.replace('_wgs_metrics', '')
        print(f"Processing {sample_name}...")

        try:
            metrics = parse_wgs_metrics(metrics_file)
            metrics['sample'] = sample_name
            results.append(metrics)
        except Exception as e:
            print(f"Error parsing {metrics_file}: {e}")
            continue

    # Create DataFrame
    df = pd.DataFrame(results)

    # Reorder columns
    cols = ['sample', 'mean_coverage', 'median_coverage', 'sd_coverage',
            'pct_5x', 'pct_10x', 'pct_20x', 'pct_30x',
            'callable_Mb_5x', 'callable_Mb_10x', 'callable_Mb_20x', 'callable_Mb_30x',
            'callable_bases_5x', 'callable_bases_10x', 'callable_bases_20x', 'callable_bases_30x',
            'genome_territory']
    df = df[cols]

    # Save to file
    df.to_csv(args.output, sep='\t', index=False, float_format='%.2f')
    print(f"\nSaved coverage summary to {args.output}")

    # Print summary
    print("\n" + "="*80)
    print("COVERAGE SUMMARY")
    print("="*80)
    print(df[['sample', 'mean_coverage', 'callable_Mb_5x', 'callable_Mb_10x']].to_string(index=False))
    print("\nFor mutation rate normalization:")
    print("  - Mutations per Mb: n_mutations / mean_coverage")
    print("  - Mutations per callable Mb (5x): n_mutations / callable_Mb_5x")
    print("  - Mutations per callable Mb (10x): n_mutations / callable_Mb_10x")


if __name__ == '__main__':
    main()
