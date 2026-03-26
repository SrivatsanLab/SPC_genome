#!/usr/bin/env python3
"""
Generate Lander-Waterman coverage plot from compiled QC metrics

This script creates a scatter plot of reads per reference genome megabase (X-axis)
vs proportion of reference genome covered (Y-axis) with the theoretical
Lander-Waterman curve overlay.

Usage:
    python plot_lander_waterman.py <qc_csv> <reference_fasta> <output_plot> [--experiment-col COLUMN]

Arguments:
    qc_csv: Path to compiled QC metrics CSV (from compile_qc_metrics.py)
    reference_fasta: Path to reference genome FASTA file
    output_plot: Path for output plot (e.g., lander_waterman.png)
    --experiment-col: Optional column name for grouping/coloring points (default: None)
"""

import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def get_genome_length_from_fasta(fasta_path):
    """
    Calculate total genome length from FASTA index (.fai file)

    Args:
        fasta_path: Path to FASTA file (will look for .fai index)

    Returns:
        Total genome length in base pairs
    """
    fai_path = Path(str(fasta_path) + '.fai')

    if not fai_path.exists():
        raise FileNotFoundError(f"FASTA index not found: {fai_path}\nRun: samtools faidx {fasta_path}")

    total_length = 0
    with open(fai_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                total_length += int(fields[1])  # Second column is sequence length

    return total_length

def plot_lander_waterman(qc_csv, genome_length, output_plot, experiment_col=None):
    """
    Generate Lander-Waterman coverage plot

    Args:
        qc_csv: Path to QC metrics CSV
        genome_length: Genome length in bp
        output_plot: Output plot path
        experiment_col: Optional column for hue grouping
    """
    # Read QC data
    df = pd.read_csv(qc_csv)

    # Calculate reads per megabase
    df['reads_per_mb'] = (df['total_reads'] / genome_length) * 1e6

    # Generate theoretical Lander-Waterman curve
    # L-W: fraction covered = 1 - exp(-coverage)
    # coverage = (reads * read_length) / genome_length
    # For reads_per_mb: coverage = reads_per_mb / 1e6 * genome_length / genome_length
    #                             = reads_per_mb / 1e6
    x_theory = np.linspace(0, df['reads_per_mb'].max() * 1.05, 1000)
    y_theory = 1 - np.exp(-x_theory / 1e6)  # x is reads/Mb, coverage = reads/Mb / 1e6

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))

    # Plot data points
    if experiment_col and experiment_col in df.columns:
        # Use seaborn for colored scatter by experiment
        sns.scatterplot(
            data=df,
            x='reads_per_mb',
            y='pct_1x',
            hue=experiment_col,
            edgecolor='black',
            s=100,
            ax=ax
        )
    else:
        # Simple scatter plot
        ax.scatter(
            df['reads_per_mb'],
            df['pct_1x'],
            s=100,
            edgecolor='black',
            alpha=0.7
        )

    # Plot theoretical curve
    ax.plot(
        x_theory,
        y_theory,
        color='black',
        linestyle='--',
        linewidth=2,
        label='Lander-Waterman',
        zorder=0
    )

    # Set log scales
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Labels and title
    ax.set_xlabel('Reads per Reference Genome Megabase', fontsize=14)
    ax.set_ylabel('Proportion of Reference Genome Covered (≥1X)', fontsize=14)
    ax.set_title('Coverage vs Sequencing Depth with Lander-Waterman Curve', fontsize=16)

    # Add grid
    ax.grid(True, which='both', alpha=0.3)

    # Legend
    ax.legend(fontsize=12)

    # Tight layout
    plt.tight_layout()

    # Save figure
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    print(f"Lander-Waterman plot saved to: {output_plot}")

    # Print summary statistics
    print(f"\nSummary statistics:")
    print(f"  Total cells: {len(df)}")
    print(f"  Reads per Mb range: {df['reads_per_mb'].min():.2f} - {df['reads_per_mb'].max():.2f}")
    print(f"  Coverage range: {df['pct_1x'].min():.4f} - {df['pct_1x'].max():.4f}")
    print(f"  Mean coverage: {df['pct_1x'].mean():.4f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate Lander-Waterman coverage plot from QC metrics'
    )
    parser.add_argument(
        'qc_csv',
        help='Path to compiled QC metrics CSV'
    )
    parser.add_argument(
        'reference_fasta',
        help='Path to reference genome FASTA file'
    )
    parser.add_argument(
        'output_plot',
        help='Output plot path (e.g., lander_waterman.png)'
    )
    parser.add_argument(
        '--experiment-col',
        default=None,
        help='Column name for grouping/coloring points (optional)'
    )

    args = parser.parse_args()

    # Validate inputs
    if not Path(args.qc_csv).exists():
        print(f"ERROR: QC CSV file not found: {args.qc_csv}")
        sys.exit(1)

    if not Path(args.reference_fasta).exists():
        print(f"ERROR: Reference FASTA file not found: {args.reference_fasta}")
        sys.exit(1)

    # Get genome length from FASTA index
    print(f"Reading genome length from: {args.reference_fasta}.fai")
    genome_length = get_genome_length_from_fasta(args.reference_fasta)
    print(f"Total genome length: {genome_length:,} bp ({genome_length/1e6:.2f} Mb)")

    # Generate plot
    plot_lander_waterman(
        args.qc_csv,
        genome_length,
        args.output_plot,
        args.experiment_col
    )
