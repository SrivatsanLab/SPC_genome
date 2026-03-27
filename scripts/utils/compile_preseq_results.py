#!/usr/bin/env python3
"""
Compile preseq results from multiple samples into a single CSV file.
Extracts both complexity curves (observed) and extrapolation predictions.
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path
import argparse

def parse_c_curve(filepath):
    """Parse preseq c_curve output to get observed complexity."""
    try:
        df = pd.read_csv(filepath, sep='\t')
        if len(df) == 0:
            return None
        # Get the maximum observed values
        max_row = df.iloc[-1]
        return {
            'observed_total_reads': max_row['total_reads'],
            'observed_distinct_reads': max_row['distinct_reads'],
        }
    except Exception as e:
        print(f"Error parsing {filepath}: {e}", file=sys.stderr)
        return None

def parse_lc_extrap(filepath, target_depths=[1e6, 5e6, 10e6, 20e6, 50e6, 100e6]):
    """Parse preseq lc_extrap output to get predicted complexity at target depths."""
    try:
        df = pd.read_csv(filepath, sep='\t')
        if len(df) == 0:
            return None

        result = {}
        for depth in target_depths:
            # Find the closest extrapolation point
            closest_idx = (df['TOTAL_READS'] - depth).abs().idxmin()
            closest_row = df.loc[closest_idx]

            depth_str = f"{int(depth/1e6)}M" if depth >= 1e6 else f"{int(depth/1e3)}K"
            result[f'extrap_{depth_str}_reads'] = closest_row['TOTAL_READS']
            result[f'extrap_{depth_str}_distinct'] = closest_row['EXPECTED_DISTINCT']
            result[f'extrap_{depth_str}_lower_ci'] = closest_row['LOWER_0.95CI']
            result[f'extrap_{depth_str}_upper_ci'] = closest_row['UPPER_0.95CI']

        return result
    except Exception as e:
        print(f"Error parsing {filepath}: {e}", file=sys.stderr)
        return None

def compile_preseq_results(preseq_dir, output_file):
    """Compile all preseq results into a single CSV."""

    results = []

    # Iterate through experiment directories
    for exp_dir in Path(preseq_dir).iterdir():
        if not exp_dir.is_dir():
            continue

        experiment = exp_dir.name
        print(f"Processing experiment: {experiment}")

        # Find all c_curve files
        c_curve_files = list(exp_dir.glob("*_c_curve.txt"))

        for c_curve_file in c_curve_files:
            # Extract sample name from filename
            sample_name = c_curve_file.stem.replace("_c_curve", "")

            # Corresponding lc_extrap file
            lc_extrap_file = exp_dir / f"{sample_name}_lc_extrap.txt"

            if not lc_extrap_file.exists():
                print(f"Warning: No lc_extrap file for {sample_name}", file=sys.stderr)
                continue

            # Parse both files
            c_curve_data = parse_c_curve(c_curve_file)
            lc_extrap_data = parse_lc_extrap(lc_extrap_file)

            if c_curve_data is None or lc_extrap_data is None:
                print(f"Warning: Could not parse data for {sample_name}", file=sys.stderr)
                continue

            # Combine results
            row = {
                'sample': sample_name,
                'experiment': experiment,
                **c_curve_data,
                **lc_extrap_data
            }

            results.append(row)

        print(f"  Found {len(c_curve_files)} samples")

    # Create DataFrame
    df = pd.DataFrame(results)

    # Calculate some additional metrics
    if len(df) > 0:
        df['observed_fraction_unique'] = df['observed_distinct_reads'] / df['observed_total_reads']

        # For each extrapolation depth, calculate fraction unique
        for col in df.columns:
            if col.startswith('extrap_') and col.endswith('_distinct'):
                depth_label = col.replace('extrap_', '').replace('_distinct', '')
                reads_col = f'extrap_{depth_label}_reads'
                if reads_col in df.columns:
                    df[f'extrap_{depth_label}_fraction_unique'] = df[col] / df[reads_col]

    # Save to CSV
    df.to_csv(output_file, index=False)
    print(f"\nSaved {len(df)} samples to {output_file}")

    return df

def main():
    parser = argparse.ArgumentParser(description='Compile preseq results into a single CSV')
    parser.add_argument('--preseq_dir', type=str, default='results/preseq_analysis',
                      help='Directory containing preseq results (default: results/preseq_analysis)')
    parser.add_argument('--output', type=str, default='results/preseq_analysis/compiled_preseq_results.csv',
                      help='Output CSV file (default: results/preseq_analysis/compiled_preseq_results.csv)')

    args = parser.parse_args()

    df = compile_preseq_results(args.preseq_dir, args.output)

    # Print summary statistics
    print("\n=== Summary ===")
    print(f"Total samples: {len(df)}")
    print(f"\nSamples per experiment:")
    print(df['experiment'].value_counts())
    print(f"\nObserved library complexity:")
    print(df.groupby('experiment')['observed_fraction_unique'].describe())

if __name__ == '__main__':
    main()
