#!/usr/bin/env python3
"""
Compile Lorenz curves from multiple CSV files in a directory.

Usage:
    python compile_lorenz_curves.py <input_dir> <output_file> [--pattern PATTERN]
    
Example:
    python compile_lorenz_curves.py data/HSC_enzyme_coverage/sc_outputs results/lorenz_curves.csv
    python compile_lorenz_curves.py data/coverage_benchmarking results/public_lorenz.csv --pattern "*_lorenz.csv"
"""

import sys
import os
import argparse
import pandas as pd
import glob
from pathlib import Path
from typing import Dict, List


def find_lorenz_files(input_dir: str, pattern: str = "*_lorenz.csv") -> List[str]:
    """
    Find all Lorenz curve files in the input directory.
    
    Args:
        input_dir: Directory to search for Lorenz curve files
        pattern: Glob pattern for matching files (default: "*_lorenz.csv")
    
    Returns:
        Sorted list of file paths
    """
    search_path = os.path.join(input_dir, pattern)
    lorenz_files = sorted(glob.glob(search_path))
    return lorenz_files


def read_lorenz_file(filepath: str) -> pd.Series:
    """
    Read a single Lorenz curve file and extract the relevant data.
    
    Args:
        filepath: Path to the Lorenz curve CSV file
    
    Returns:
        Pandas Series with cumulative_fraction_genome as index and 
        cumulative_fraction_reads as values
    
    Raises:
        ValueError: If expected columns are not found
    """
    df = pd.read_csv(filepath)
    
    # Check for expected columns
    if 'cumulative_fraction_genome' not in df.columns or \
       'cumulative_fraction_reads' not in df.columns:
        raise ValueError(
            f"Unexpected column names in {filepath}: {df.columns.tolist()}"
        )
    
    return df.set_index('cumulative_fraction_genome')['cumulative_fraction_reads']


def compile_lorenz_curves(input_dir: str, output_file: str, 
                          pattern: str = "*_lorenz.csv",
                          barcode_suffix: str = "_lorenz.csv") -> pd.DataFrame:
    """
    Compile multiple Lorenz curve files into a single dataframe.
    
    Args:
        input_dir: Directory containing Lorenz curve files
        output_file: Path to save the compiled CSV
        pattern: Glob pattern for matching files
        barcode_suffix: Suffix to remove from filenames to get barcode names
    
    Returns:
        Combined DataFrame with all Lorenz curves
    
    Raises:
        FileNotFoundError: If no Lorenz curve files are found
        ValueError: If no valid Lorenz curves could be read
    """
    print(f"Input directory: {input_dir}")
    print(f"Output file: {output_file}")
    print(f"Search pattern: {pattern}")
    print()
    
    # Find all Lorenz curve files
    lorenz_files = find_lorenz_files(input_dir, pattern)
    
    if len(lorenz_files) == 0:
        raise FileNotFoundError(f"No Lorenz curve files found in {input_dir}")
    
    print(f"Found {len(lorenz_files)} Lorenz curve files")
    
    # Read all files and combine
    dfs = {}
    skipped = []
    
    for filepath in lorenz_files:
        barcode = os.path.basename(filepath).replace(barcode_suffix, '')
        try:
            dfs[barcode] = read_lorenz_file(filepath)
            print(f"  ✓ {barcode}")
        except ValueError as e:
            print(f"  ✗ {barcode}: {e}")
            skipped.append(filepath)
    
    if not dfs:
        raise ValueError(
            f"No valid Lorenz curves found. Skipped {len(skipped)} files."
        )
    
    # Combine into single dataframe
    combined_df = pd.DataFrame(dfs)
    combined_df.index.name = 'cumulative_fraction_genome'
    
    # Create output directory if it doesn't exist
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save to file
    combined_df.to_csv(output_file)
    
    print()
    print(f"✓ Successfully compiled {len(dfs)} Lorenz curves")
    print(f"  Output: {output_file}")
    print(f"  Shape: {combined_df.shape}")
    
    if skipped:
        print(f"  Warning: Skipped {len(skipped)} files with errors")
    
    return combined_df


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Compile Lorenz curves from multiple CSV files into a single dataframe.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s data/HSC_enzyme_coverage/sc_outputs results/hsc_lorenz.csv
  %(prog)s data/coverage_benchmarking results/public_lorenz.csv --pattern "*_lorenz.csv"
        """
    )
    
    parser.add_argument(
        'input_dir',
        help='Directory containing Lorenz curve CSV files'
    )
    
    parser.add_argument(
        'output_file',
        help='Path to save the compiled Lorenz curves CSV'
    )
    
    parser.add_argument(
        '--pattern',
        default='*_lorenz.csv',
        help='Glob pattern for matching Lorenz curve files (default: *_lorenz.csv)'
    )
    
    parser.add_argument(
        '--barcode-suffix',
        default='_lorenz.csv',
        help='Suffix to remove from filenames to extract barcode names (default: _lorenz.csv)'
    )
    
    args = parser.parse_args()
    
    try:
        compile_lorenz_curves(
            args.input_dir,
            args.output_file,
            pattern=args.pattern,
            barcode_suffix=args.barcode_suffix
        )
        return 0
    except (FileNotFoundError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"UNEXPECTED ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
