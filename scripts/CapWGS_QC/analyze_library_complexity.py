#!/usr/bin/env python3
"""
Analyze library complexity by comparing:
1. Preseq predictions from shallow NextSeq runs
2. Lander-Waterman theoretical predictions
3. Observed coverage from deep NovaSeq runs

This script investigates why shallow sequencing appears to show higher library
quality than is realized in deeper sequencing.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def lander_waterman(reads_per_mb, genome_length_mb, read_length):
    """
    Calculate expected fraction of genome covered using Lander-Waterman equation.

    coverage = (reads_per_mb * read_length) / 1e6
    fraction_covered = 1 - exp(-coverage)
    """
    coverage = (reads_per_mb * read_length) / 1e6
    return 1 - np.exp(-coverage)

def load_and_merge_data():
    """Load coverage QC data and preseq results, merge by sample."""

    # Load coverage data
    cap_data = pd.read_csv('results/benchmarking_coverage/ALL_SPC_experiments_qc.csv', index_col=0)
    cap_data['reads_per_mb'] = (cap_data['total_reads'] / cap_data['genome_length_reference']) * 1e6

    # Load preseq data
    preseq = pd.read_csv('results/preseq_analysis/compiled_preseq_results.csv')

    # Fix naming inconsistency
    cap_data['experiment'] = cap_data['experiment'].replace('worm_unsyc_ns', 'worm_unsync_ns')

    return cap_data, preseq

def create_comparison_df(cap_data, preseq):
    """
    Create comparison DataFrame with:
    - NextSeq shallow data + preseq predictions
    - NovaSeq deep data (observed)
    """

    # Define experiment pairs (NextSeq -> NovaSeq)
    pairs = {
        'HSC_ns': ['HSC_nova', 'HSC2_nova_big'],
        'worm_unsync_ns': [],  # No corresponding NovaSeq data
        'worm_L4_ns': ['worm_L4_nova']
    }

    comparison_data = []

    for ns_exp, nova_exps in pairs.items():
        # Get NextSeq data
        ns_data = cap_data[cap_data['experiment'] == ns_exp].copy()
        ns_preseq = preseq[preseq['experiment'] == ns_exp].copy()

        # Merge on sample name
        ns_data['sample'] = [idx.split('_' + ns_exp)[0] for idx in ns_data.index]
        merged_ns = ns_data.merge(ns_preseq, on='sample', how='inner', suffixes=('_obs', '_preseq'))

        for _, row in merged_ns.iterrows():
            # NextSeq observed data point
            comparison_data.append({
                'sample': row['sample'],
                'experiment_pair': f"{ns_exp} -> {','.join(nova_exps) if nova_exps else 'None'}",
                'sequencer': 'NextSeq',
                'experiment': ns_exp,
                'reads_per_mb': row['reads_per_mb'],
                'pct_1x': row['pct_1x'],
                'total_reads': row['total_reads'],
                'mean_read_length': row['mean_read_length'],
                'genome_length': row['genome_length_reference'],
                'pct_duplicates': row['pct_duplicates'],
                'observed_fraction_unique': row['observed_fraction_unique'],
                'data_type': 'Observed'
            })

            # Preseq predictions at various depths
            genome_mb = row['genome_length_reference'] / 1e6
            for depth_label in ['10M', '20M', '50M', '100M']:
                total_reads = row[f'extrap_{depth_label}_reads']
                distinct_reads = row[f'extrap_{depth_label}_distinct']

                # Calculate expected coverage from these distinct reads
                # Assume each distinct read contributes to coverage
                reads_per_mb_pred = total_reads / genome_mb

                # Estimate fraction covered using distinct/total as proxy for library quality
                # Then apply Lander-Waterman with effective coverage
                frac_unique = distinct_reads / total_reads if total_reads > 0 else 0
                effective_reads = total_reads * frac_unique
                effective_reads_per_mb = effective_reads / genome_mb
                pct_1x_pred = lander_waterman(effective_reads_per_mb, genome_mb, row['mean_read_length'])

                comparison_data.append({
                    'sample': row['sample'],
                    'experiment_pair': f"{ns_exp} -> {','.join(nova_exps) if nova_exps else 'None'}",
                    'sequencer': f'Preseq_{depth_label}',
                    'experiment': ns_exp,
                    'reads_per_mb': reads_per_mb_pred,
                    'pct_1x': pct_1x_pred,
                    'total_reads': total_reads,
                    'mean_read_length': row['mean_read_length'],
                    'genome_length': row['genome_length_reference'],
                    'pct_duplicates': 1 - frac_unique,
                    'observed_fraction_unique': frac_unique,
                    'data_type': f'Preseq {depth_label}'
                })

        # Get corresponding NovaSeq data
        for nova_exp in nova_exps:
            nova_data = cap_data[cap_data['experiment'] == nova_exp].copy()
            nova_data['sample'] = [idx.split('_' + nova_exp)[0] for idx in nova_data.index]

            # Only include samples that were in NextSeq run
            nova_data = nova_data[nova_data['sample'].isin(merged_ns['sample'])]

            for _, row in nova_data.iterrows():
                comparison_data.append({
                    'sample': row['sample'],
                    'experiment_pair': f"{ns_exp} -> {','.join(nova_exps)}",
                    'sequencer': 'NovaSeq',
                    'experiment': nova_exp,
                    'reads_per_mb': row['reads_per_mb'],
                    'pct_1x': row['pct_1x'],
                    'total_reads': row['total_reads'],
                    'mean_read_length': row['mean_read_length'],
                    'genome_length': row['genome_length_reference'],
                    'pct_duplicates': row['pct_duplicates'],
                    'observed_fraction_unique': np.nan,  # Not directly available
                    'data_type': 'Observed'
                })

    return pd.DataFrame(comparison_data)

def plot_comparison(comparison_df, output_dir='results/preseq_analysis'):
    """Create comparison plots."""

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Color palette
    colors = {
        'NextSeq': '#438CFD',
        'Preseq_10M': '#FF8C00',
        'Preseq_20M': '#FF1B5E',
        'Preseq_50M': '#5941A9',
        'Preseq_100M': '#E980FC',
        'NovaSeq': '#80F15E'
    }

    # Plot for each experiment pair
    for pair in comparison_df['experiment_pair'].unique():
        pair_data = comparison_df[comparison_df['experiment_pair'] == pair]

        if len(pair_data) == 0:
            continue

        # Get genome info for Lander-Waterman curve
        genome_length = pair_data['genome_length'].iloc[0]
        read_length = pair_data['mean_read_length'].mean()
        genome_mb = genome_length / 1e6

        # Create Lander-Waterman theoretical curve
        x_theory = np.logspace(0, np.log10(pair_data['reads_per_mb'].max() * 1.2), 1000)
        y_theory = lander_waterman(x_theory, genome_mb, read_length)

        # Create plot
        fig, ax = plt.subplots(1, 1, figsize=(14, 8))

        # Plot Lander-Waterman
        ax.plot(x_theory, y_theory, 'k--', linewidth=2, label='Lander-Waterman', zorder=1)

        # Plot observed NextSeq data
        ns_data = pair_data[pair_data['sequencer'] == 'NextSeq']
        if len(ns_data) > 0:
            ax.scatter(ns_data['reads_per_mb'], ns_data['pct_1x'],
                      c=colors['NextSeq'], s=50, alpha=0.7,
                      edgecolors='black', linewidth=0.5,
                      label='NextSeq (observed)', zorder=3)

        # Plot preseq predictions
        for depth in ['10M', '20M', '50M', '100M']:
            preseq_data = pair_data[pair_data['sequencer'] == f'Preseq_{depth}']
            if len(preseq_data) > 0:
                ax.scatter(preseq_data['reads_per_mb'], preseq_data['pct_1x'],
                          c=colors[f'Preseq_{depth}'], s=30, alpha=0.5,
                          edgecolors='black', linewidth=0.5,
                          label=f'Preseq @ {depth}', zorder=2)

        # Plot observed NovaSeq data
        nova_data = pair_data[pair_data['sequencer'] == 'NovaSeq']
        if len(nova_data) > 0:
            ax.scatter(nova_data['reads_per_mb'], nova_data['pct_1x'],
                      c=colors['NovaSeq'], s=100, alpha=0.8,
                      edgecolors='black', linewidth=1,
                      marker='s', label='NovaSeq (observed)', zorder=4)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Reads per Megabase', fontsize=14, fontweight='bold')
        ax.set_ylabel('Fraction of Genome Covered (≥1x)', fontsize=14, fontweight='bold')
        ax.set_title(f'Library Complexity: {pair}', fontsize=16, fontweight='bold')
        ax.legend(loc='lower right', fontsize=10)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        # Save
        safe_name = pair.replace(' ', '_').replace('->', 'to').replace(',', '_')
        plt.savefig(f'{output_dir}/{safe_name}_complexity_comparison.png', dpi=300)
        plt.close()

        print(f"Saved plot: {output_dir}/{safe_name}_complexity_comparison.png")

def print_summary_statistics(comparison_df):
    """Print summary statistics comparing predictions vs observations."""

    print("\n" + "="*80)
    print("LIBRARY COMPLEXITY ANALYSIS SUMMARY")
    print("="*80)

    for pair in comparison_df['experiment_pair'].unique():
        pair_data = comparison_df[comparison_df['experiment_pair'] == pair]

        print(f"\n{pair}")
        print("-" * 80)

        # NextSeq stats
        ns_data = pair_data[pair_data['sequencer'] == 'NextSeq']
        if len(ns_data) > 0:
            print(f"\nNextSeq (shallow sequencing):")
            print(f"  N samples: {len(ns_data)}")
            print(f"  Reads/Mb: {ns_data['reads_per_mb'].mean():.1f} ± {ns_data['reads_per_mb'].std():.1f}")
            print(f"  Coverage (≥1x): {ns_data['pct_1x'].mean():.3f} ± {ns_data['pct_1x'].std():.3f}")
            print(f"  Duplicates: {ns_data['pct_duplicates'].mean():.3f} ± {ns_data['pct_duplicates'].std():.3f}")

        # Preseq predictions
        print(f"\nPreseq predictions from NextSeq:")
        for depth in ['10M', '20M', '50M', '100M']:
            preseq_data = pair_data[pair_data['sequencer'] == f'Preseq_{depth}']
            if len(preseq_data) > 0:
                print(f"  @ {depth} reads:")
                print(f"    Predicted reads/Mb: {preseq_data['reads_per_mb'].mean():.1f}")
                print(f"    Predicted coverage: {preseq_data['pct_1x'].mean():.3f} ± {preseq_data['pct_1x'].std():.3f}")
                print(f"    Predicted duplicates: {preseq_data['pct_duplicates'].mean():.3f} ± {preseq_data['pct_duplicates'].std():.3f}")

        # NovaSeq stats
        nova_data = pair_data[pair_data['sequencer'] == 'NovaSeq']
        if len(nova_data) > 0:
            print(f"\nNovaSeq (deep sequencing, OBSERVED):")
            print(f"  N samples: {len(nova_data)}")
            print(f"  Reads/Mb: {nova_data['reads_per_mb'].mean():.1f} ± {nova_data['reads_per_mb'].std():.1f}")
            print(f"  Coverage (≥1x): {nova_data['pct_1x'].mean():.3f} ± {nova_data['pct_1x'].std():.3f}")
            print(f"  Duplicates: {nova_data['pct_duplicates'].mean():.3f} ± {nova_data['pct_duplicates'].std():.3f}")

            # Compare to preseq prediction at similar depth
            nova_avg_reads_mb = nova_data['reads_per_mb'].mean()

            # Find closest preseq prediction
            for depth in ['10M', '20M', '50M', '100M']:
                preseq_data = pair_data[pair_data['sequencer'] == f'Preseq_{depth}']
                if len(preseq_data) > 0:
                    preseq_avg_reads_mb = preseq_data['reads_per_mb'].mean()

                    if abs(np.log10(preseq_avg_reads_mb) - np.log10(nova_avg_reads_mb)) < 0.5:
                        print(f"\n  Comparison to Preseq @ {depth} (similar depth):")
                        print(f"    Predicted coverage: {preseq_data['pct_1x'].mean():.3f}")
                        print(f"    Observed coverage: {nova_data['pct_1x'].mean():.3f}")
                        print(f"    Discrepancy: {(preseq_data['pct_1x'].mean() - nova_data['pct_1x'].mean()):.3f}")
                        print(f"    Predicted duplicates: {preseq_data['pct_duplicates'].mean():.3f}")
                        print(f"    Observed duplicates: {nova_data['pct_duplicates'].mean():.3f}")
                        print(f"    Duplicate discrepancy: {(nova_data['pct_duplicates'].mean() - preseq_data['pct_duplicates'].mean()):.3f}")

def main():
    """Main analysis function."""

    print("Loading data...")
    cap_data, preseq = load_and_merge_data()

    print("Creating comparison dataset...")
    comparison_df = create_comparison_df(cap_data, preseq)

    # Save comparison data
    output_file = 'results/preseq_analysis/library_complexity_comparison.csv'
    comparison_df.to_csv(output_file, index=False)
    print(f"Saved comparison data: {output_file}")

    print("\nGenerating plots...")
    plot_comparison(comparison_df)

    print_summary_statistics(comparison_df)

    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)

if __name__ == '__main__':
    main()
