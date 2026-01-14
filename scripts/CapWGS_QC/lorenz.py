#!/usr/bin/env python3
"""
Generate Lorenz curve from bigwig file for coverage uniformity analysis.

This script reads coverage from a bigwig file, computes the Lorenz curve,
calculates the Gini coefficient, and outputs results to CSV.
"""

import argparse
import numpy as np
import pandas as pd
import pyBigWig
import sys
from pathlib import Path

def create_uniform_bins(chrom_lengths, binsize):
    """Create uniform genomic bins."""
    uniform_bins = []
    for chrom, length in chrom_lengths.items():
        for start in range(0, length, binsize):
            end = min(start + binsize, length)
            uniform_bins.append((chrom, start, end))
    return uniform_bins

def map_coverage_to_bins(bigwig_file, uniform_bins):
    """Map coverage values from bigwig to uniform bins."""
    coverage = []
    bw = pyBigWig.open(bigwig_file)
    for chrom, start, end in uniform_bins:
        try:
            values = bw.values(chrom, start, end, numpy=False)
            interval_value = values[0] if values else 0
            coverage.append(interval_value if interval_value is not None else 0)
        except RuntimeError:
            # Handle chromosomes not in bigwig
            coverage.append(0)
    bw.close()
    return np.array(coverage)

def compute_lorenz_curve(coverage, n_points=100):
    """
    Compute Lorenz curve from coverage data.

    Returns:
        x: x-axis values (cumulative fraction of genome)
        y: y-axis values (cumulative fraction of reads)
        gini: Gini coefficient
    """
    # Compute proportion of total coverage
    total_coverage = coverage.sum()
    if total_coverage == 0:
        raise ValueError("Total coverage is zero - cannot compute Lorenz curve")

    prop = coverage / total_coverage

    # Filter out zeros and sort
    prop = prop[prop > 0]
    prop_sorted = np.sort(prop)

    # Compute cumulative sum
    y = np.cumsum(prop_sorted)

    # Create x-axis (cumulative fraction of bins)
    x = np.arange(len(prop)) / len(prop)

    # Prepend zeros for proper Lorenz curve
    x = np.concatenate(([0], x))
    y = np.concatenate(([0], y))

    # Compute Gini coefficient
    auc = np.trapz(y, x)
    gini = 1 - 2 * auc

    # Interpolate to common scale for consistent output
    common_x = np.linspace(0, 1, n_points)
    common_y = np.interp(common_x, x, y)

    return common_x, common_y, gini

def main():
    parser = argparse.ArgumentParser(
        description='Compute Lorenz curve from bigwig file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate Lorenz curve with default 1kb bins
  %(prog)s sample.bw -o sample_lorenz.csv

  # Use 10kb bins and save to specific directory
  %(prog)s sample.bw -b 10000 -o results/sample_lorenz.csv

  # Change number of output points
  %(prog)s sample.bw -o sample_lorenz.csv -n 200
"""
    )

    parser.add_argument('bigwig', type=str,
                        help='Input bigwig file')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Output CSV file for Lorenz curve')
    parser.add_argument('-b', '--binsize', type=int, default=1000,
                        help='Bin size for coverage calculation (default: 1000)')
    parser.add_argument('-n', '--npoints', type=int, default=100,
                        help='Number of points in output Lorenz curve (default: 100)')
    parser.add_argument('--gini-output', type=str,
                        help='Optional file to write Gini coefficient (default: print to stdout)')

    args = parser.parse_args()

    # Validate inputs
    if not Path(args.bigwig).exists():
        print(f"ERROR: Bigwig file not found: {args.bigwig}", file=sys.stderr)
        sys.exit(1)

    # Create output directory if needed
    output_dir = Path(args.output).parent
    if output_dir != Path('.') and not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading bigwig: {args.bigwig}")

    # Get chromosome lengths from bigwig
    bw = pyBigWig.open(args.bigwig)
    chrom_lengths = bw.chroms()
    bw.close()

    # Create bins and extract coverage
    print(f"Creating {args.binsize}bp bins...")
    bins = create_uniform_bins(chrom_lengths, args.binsize)
    print(f"Extracting coverage for {len(bins):,} bins...")
    coverage = map_coverage_to_bins(args.bigwig, bins)

    # Compute Lorenz curve
    print("Computing Lorenz curve...")
    x, y, gini = compute_lorenz_curve(coverage, n_points=args.npoints)

    # Save Lorenz curve
    lorenz_df = pd.DataFrame({
        'cumulative_fraction_genome': x,
        'cumulative_fraction_reads': y
    })
    lorenz_df.to_csv(args.output, index=False)
    print(f"Lorenz curve saved to: {args.output}")

    # Output Gini coefficient
    gini_message = f"Gini coefficient: {gini:.6f}"
    print(gini_message)

    if args.gini_output:
        with open(args.gini_output, 'w') as f:
            f.write(f"{gini}\n")
        print(f"Gini coefficient saved to: {args.gini_output}")

    print("Done!")

if __name__ == '__main__':
    main()
