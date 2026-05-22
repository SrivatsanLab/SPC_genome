#!/usr/bin/env python3
"""
Compile per-cell Lorenz curves and Gini coefficients into a single long-format CSV.

Scans data/benchmarking/lorenz/{dataset}/ for *_lorenz.csv files,
pairs each with its *_gini.txt, and outputs a combined CSV.

Usage:
    python scripts/benchmarking/compile_lorenz.py \
        --lorenz-dir data/benchmarking/lorenz/public_WGA \
        --output results/benchmarking/public_WGS_lorenz_series.csv
"""

import argparse
import pandas as pd
from pathlib import Path
import re
import sys


def parse_filename(path):
    """Extract cell_id and depth_label from filename like SRR5365341_10M_lorenz.csv."""
    stem = path.stem.replace('_lorenz', '')
    match = re.match(r'^(.+?)_(\d+M)$', stem)
    if not match:
        return None, None
    return match.group(1), match.group(2)


def main():
    parser = argparse.ArgumentParser(description='Compile Lorenz curves into long-format CSV')
    parser.add_argument('--lorenz-dir', required=True,
                        help='Directory containing *_lorenz.csv and *_gini.txt files')
    parser.add_argument('--output', required=True,
                        help='Output CSV path')
    args = parser.parse_args()

    lorenz_dir = Path(args.lorenz_dir)
    dataset = lorenz_dir.name

    lorenz_files = sorted(lorenz_dir.glob('*_lorenz.csv'))
    if not lorenz_files:
        print(f"Error: No *_lorenz.csv files found in {lorenz_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(lorenz_files)} Lorenz curve files in {lorenz_dir}")

    records = []
    missing_gini = []

    for lf in lorenz_files:
        cell_id, depth_label = parse_filename(lf)
        if cell_id is None:
            print(f"Warning: Could not parse filename {lf.name}, skipping")
            continue

        gini_file = lf.parent / lf.name.replace('_lorenz.csv', '_gini.txt')
        if not gini_file.exists():
            missing_gini.append(lf.name)
            gini = float('nan')
        else:
            gini = float(gini_file.read_text().strip())

        df = pd.read_csv(lf)
        df['dataset'] = dataset
        df['cell_id'] = cell_id
        df['depth_label'] = depth_label
        df['gini'] = gini

        records.append(df)

    if missing_gini:
        print(f"Warning: Missing Gini files for {len(missing_gini)} cells: {missing_gini[:5]}...")

    compiled = pd.concat(records, ignore_index=True)
    compiled = compiled[['dataset', 'cell_id', 'depth_label', 'gini',
                          'cumulative_fraction_genome', 'cumulative_fraction_reads']]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    compiled.to_csv(output_path, index=False)

    print(f"Compiled {len(compiled)} rows ({len(lorenz_files)} cells × 100 points) → {output_path}")


if __name__ == '__main__':
    main()
