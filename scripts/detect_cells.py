#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.optimize import minimize

import argparse
import sys

def ordmag_algorithm(df, expected_cells=None):
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    # Function to compute OrdMag(x)
    def ordmag(x):
        # Calculate the 99th percentile of the top x readcounts
        m = df['read_count'].iloc[:int(x)].quantile(0.99)
        # Count barcodes with readcounts >= m / 10
        return (df['read_count'] >= m / 10).sum()

    if expected_cells is None:
        # Loss function for optimization
        def loss_function(x):
            ordmag_x = ordmag(x)
            return ((ordmag_x - x) ** 2) / x

        # Minimize the loss function to find the optimal number of barcodes
        result = minimize(loss_function, x0=len(df)//2, bounds=[(1, len(df))])
        expected_cells = int(result.x[0])

    # Calculate the threshold using the expected or estimated number of cells
    m = df['read_count'].iloc[:expected_cells].quantile(0.99)
    threshold = m / 10

    # Identify barcodes that are considered cells (readcounts >= threshold)
    cell_barcodes = df[df['read_count'] >= threshold]['barcode'].values

    return cell_barcodes, expected_cells

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Knee plot and cell detection using barcode read counts.')
    parser.add_argument('-i', '--input', type=str, help='Path to the input CSV file with barcode read counts. If not specified, reads from stdin.')
    parser.add_argument('-o', '--output', type=str, help='Output file for cell barcodes. If not specified, prints to stdout.')
    parser.add_argument('--plot', type=str, help='File path to save the knee plot image.')

    args = parser.parse_args()

    # Read input file or standard input
    if args.input:
        read_counts = pd.read_csv(args.input)
    else:
        read_counts = pd.read_csv(sys.stdin)

    # Ensure data is sorted by read_count in descending order
    read_counts = read_counts.sort_values(by='read_count', ascending=False)

    cell_barcodes, optimal_cells = ordmag_algorithm(read_counts)

    if args.plot:
        plt.figure(figsize=(8, 6))
        plt.plot(read_counts['barcode'].to_numpy(), read_counts['read_count'].to_numpy(), linestyle='-', color='blue')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(1,len(read_counts['barcode']))
        plt.axvline(x=optimal_cells, color='red', label=f'predicted cells: {optimal_cells}', linestyle='--', linewidth=1)
        plt.xlabel('Barcode Rank', fontdict={'weight':'bold'})
        plt.ylabel('Read Count', fontdict={'weight':'bold'})
        plt.legend()
        plt.savefig(args.plot, dpi=300, bbox_inches='tight')

    # Output cell barcodes to file or stdout
    if args.output:
        np.savetxt(args.output, cell_barcodes, fmt='%s')
        # print(f"Cell barcodes saved to {args.output}")
    else:
        print("\n".join(cell_barcodes))
        
if __name__ == "__main__":
    main()