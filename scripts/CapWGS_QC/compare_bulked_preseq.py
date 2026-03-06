#!/usr/bin/env python3
"""
Compare preseq results from NextSeq vs NovaSeq bulked alignments.

This analysis addresses the question: Why do shallow NextSeq runs appear to show
better library quality than is observed in deeper NovaSeq sequencing?

We compare:
1. Preseq predictions from NextSeq bulked BAMs
2. Preseq predictions from NovaSeq bulked BAMs
3. Observed coverage metrics from both sequencers
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys

def parse_preseq_files(preseq_dir):
    """Parse all preseq output files in a directory."""

    preseq_dir = Path(preseq_dir)

    if not preseq_dir.exists():
        print(f"Error: Directory {preseq_dir} does not exist")
        sys.exit(1)

    results = []

    # Find all c_curve files
    c_curve_files = list(preseq_dir.glob("*_c_curve.txt"))

    if len(c_curve_files) == 0:
        print(f"Error: No c_curve files found in {preseq_dir}")
        sys.exit(1)

    for c_curve_file in c_curve_files:
        sample_name = c_curve_file.stem.replace("_c_curve", "")
        lc_extrap_file = preseq_dir / f"{sample_name}_lc_extrap.txt"

        if not lc_extrap_file.exists():
            print(f"Warning: No lc_extrap file for {sample_name}")
            continue

        # Parse c_curve
        try:
            c_curve_df = pd.read_csv(c_curve_file, sep='\t')
            if len(c_curve_df) > 0:
                max_row = c_curve_df.iloc[-1]
                observed_total = max_row['total_reads']
                observed_distinct = max_row['distinct_reads']
            else:
                observed_total = 0
                observed_distinct = 0
        except Exception as e:
            print(f"Warning: Error parsing {c_curve_file}: {e}")
            observed_total = 0
            observed_distinct = 0

        # Parse lc_extrap
        try:
            lc_extrap_df = pd.read_csv(lc_extrap_file, sep='\t')

            # Get predictions at specific depths
            extrap_data = {}
            for depth in [1e6, 5e6, 10e6, 20e6, 50e6, 100e6, 200e6]:
                closest_idx = (lc_extrap_df['TOTAL_READS'] - depth).abs().idxmin()
                row = lc_extrap_df.loc[closest_idx]

                depth_label = f"{int(depth/1e6)}M"
                extrap_data[f'{depth_label}_reads'] = row['TOTAL_READS']
                extrap_data[f'{depth_label}_distinct'] = row['EXPECTED_DISTINCT']
                extrap_data[f'{depth_label}_lower_ci'] = row['LOWER_0.95CI']
                extrap_data[f'{depth_label}_upper_ci'] = row['UPPER_0.95CI']
                extrap_data[f'{depth_label}_frac_unique'] = row['EXPECTED_DISTINCT'] / row['TOTAL_READS']
        except Exception as e:
            print(f"Warning: Error parsing {lc_extrap_file}: {e}")
            extrap_data = {}

        result = {
            'sample': sample_name,
            'observed_total_reads': observed_total,
            'observed_distinct_reads': observed_distinct,
            'observed_frac_unique': observed_distinct / observed_total if observed_total > 0 else np.nan,
            **extrap_data
        }

        results.append(result)

    return pd.DataFrame(results)

def annotate_samples(df):
    """Add metadata about sample type and sequencer."""

    def classify_sample(name):
        if 'Enzyme_amp_HSC_UDI_5' in name:
            return 'HSC_Enzyme', 'NextSeq'
        elif 'Cell_Div_HSC_UDI_6' in name:
            return 'HSC_CellDiv', 'NextSeq'
        elif 'L4_worm_sci_5_ns' in name:
            return 'Worm_L4', 'NextSeq'
        elif 'worm_CapGTA_unsync_pilot_rerun' in name:
            return 'Worm_unsync', 'NextSeq'
        elif 'HSC_enzyme_coverage' in name:
            return 'HSC_Enzyme', 'NovaSeq'
        elif 'HSC_CellGrowth_coverage' in name:
            return 'HSC_CellDiv', 'NovaSeq'
        elif 'worm_CapGTA_UDI_5' in name:
            return 'Worm_L4', 'NovaSeq'
        else:
            return 'Unknown', 'Unknown'

    df[['experiment', 'sequencer']] = df['sample'].apply(
        lambda x: pd.Series(classify_sample(x))
    )

    return df

def plot_complexity_curves(df, output_dir='results/preseq_analysis/bulked_comparison'):
    """Create comparison plots."""

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Colors
    colors = {
        'NextSeq': '#438CFD',
        'NovaSeq': '#80F15E'
    }

    experiments = df['experiment'].unique()
    experiments = [e for e in experiments if e != 'Unknown']

    for exp in experiments:
        exp_data = df[df['experiment'] == exp]

        fig, axes = plt.subplots(1, 2, figsize=(16, 6))

        # Plot 1: Fraction unique vs total reads
        ax = axes[0]

        # Extrapolation curves
        for _, row in exp_data.iterrows():
            seq = row['sequencer']

            # Build curve from extrapolation data
            depths = []
            frac_unique = []

            for depth_label in ['1M', '5M', '10M', '20M', '50M', '100M', '200M']:
                if f'{depth_label}_reads' in row and not pd.isna(row[f'{depth_label}_reads']):
                    depths.append(row[f'{depth_label}_reads'])
                    frac_unique.append(row[f'{depth_label}_frac_unique'])

            if len(depths) > 0:
                ax.plot(depths, frac_unique,
                       color=colors[seq],
                       linewidth=2,
                       label=f'{seq}' if row.name == exp_data[exp_data['sequencer']==seq].index[0] else '',
                       alpha=0.8)

                # Add observed point
                if row['observed_total_reads'] > 0:
                    ax.scatter([row['observed_total_reads']],
                             [row['observed_frac_unique']],
                             color=colors[seq],
                             s=150,
                             edgecolors='black',
                             linewidth=2,
                             zorder=5)

        ax.set_xlabel('Total Reads', fontsize=12, fontweight='bold')
        ax.set_ylabel('Fraction of Reads that are Unique', fontsize=12, fontweight='bold')
        ax.set_xscale('log')
        ax.set_title(f'{exp}: Library Complexity', fontsize=14, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        # Plot 2: Distinct reads vs total reads
        ax = axes[1]

        for _, row in exp_data.iterrows():
            seq = row['sequencer']

            # Build curve
            total_reads = []
            distinct_reads = []

            for depth_label in ['1M', '5M', '10M', '20M', '50M', '100M', '200M']:
                if f'{depth_label}_reads' in row and not pd.isna(row[f'{depth_label}_reads']):
                    total_reads.append(row[f'{depth_label}_reads'])
                    distinct_reads.append(row[f'{depth_label}_distinct'])

            if len(total_reads) > 0:
                ax.plot(total_reads, distinct_reads,
                       color=colors[seq],
                       linewidth=2,
                       label=f'{seq}' if row.name == exp_data[exp_data['sequencer']==seq].index[0] else '',
                       alpha=0.8)

                # Add observed point
                if row['observed_total_reads'] > 0:
                    ax.scatter([row['observed_total_reads']],
                             [row['observed_distinct_reads']],
                             color=colors[seq],
                             s=150,
                             edgecolors='black',
                             linewidth=2,
                             zorder=5)

        # Add diagonal (perfect uniqueness)
        max_val = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([0, max_val], [0, max_val], 'k--', linewidth=1, alpha=0.5, label='All reads unique')

        ax.set_xlabel('Total Reads', fontsize=12, fontweight='bold')
        ax.set_ylabel('Distinct Reads', fontsize=12, fontweight='bold')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(f'{exp}: Distinct vs Total Reads', fontsize=14, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(f'{output_dir}/{exp}_preseq_comparison.png', dpi=300)
        plt.close()

        print(f"Saved: {output_dir}/{exp}_preseq_comparison.png")

def print_comparison_summary(df):
    """Print detailed comparison statistics."""

    print("\n" + "="*80)
    print("PRESEQ COMPARISON: NextSeq vs NovaSeq Bulked Alignments")
    print("="*80)

    experiments = df['experiment'].unique()
    experiments = [e for e in experiments if e != 'Unknown']

    for exp in experiments:
        exp_data = df[df['experiment'] == exp]

        print(f"\n{exp}")
        print("-" * 80)

        for seq in ['NextSeq', 'NovaSeq']:
            seq_data = exp_data[exp_data['sequencer'] == seq]

            if len(seq_data) == 0:
                continue

            row = seq_data.iloc[0]

            print(f"\n{seq}:")
            print(f"  Observed:")
            print(f"    Total reads: {row['observed_total_reads']:,.0f}")
            print(f"    Distinct reads: {row['observed_distinct_reads']:,.0f}")
            print(f"    Fraction unique: {row['observed_frac_unique']:.4f}")
            print(f"    Duplicates: {1 - row['observed_frac_unique']:.4f}")

            print(f"\n  Predictions:")
            for depth_label in ['10M', '50M', '100M', '200M']:
                if f'{depth_label}_reads' in row and not pd.isna(row[f'{depth_label}_reads']):
                    print(f"    @ {depth_label} reads:")
                    print(f"      Distinct: {row[f'{depth_label}_distinct']:,.0f}")
                    print(f"      Fraction unique: {row[f'{depth_label}_frac_unique']:.4f}")
                    print(f"      Duplicates: {1 - row[f'{depth_label}_frac_unique']:.4f}")

        # Direct comparison
        ns_data = exp_data[exp_data['sequencer'] == 'NextSeq']
        nova_data = exp_data[exp_data['sequencer'] == 'NovaSeq']

        if len(ns_data) > 0 and len(nova_data) > 0:
            print(f"\n  Comparison at similar depths:")

            for depth_label in ['50M', '100M']:
                ns_row = ns_data.iloc[0]
                nova_row = nova_data.iloc[0]

                if f'{depth_label}_frac_unique' in ns_row and f'{depth_label}_frac_unique' in nova_row:
                    if not pd.isna(ns_row[f'{depth_label}_frac_unique']) and not pd.isna(nova_row[f'{depth_label}_frac_unique']):
                        ns_frac = ns_row[f'{depth_label}_frac_unique']
                        nova_frac = nova_row[f'{depth_label}_frac_unique']

                        print(f"\n    @ {depth_label}:")
                        print(f"      NextSeq predicted fraction unique: {ns_frac:.4f}")
                        print(f"      NovaSeq predicted fraction unique: {nova_frac:.4f}")
                        print(f"      Difference: {ns_frac - nova_frac:.4f}")
                        print(f"      Fold overestimate: {ns_frac / nova_frac if nova_frac > 0 else np.inf:.2f}x")

def main():
    """Main analysis."""

    preseq_dir = 'results/preseq_analysis/bulked_comparison'

    print(f"Loading preseq results from: {preseq_dir}")
    df = parse_preseq_files(preseq_dir)

    print(f"Found {len(df)} samples")

    # Annotate
    df = annotate_samples(df)

    # Save
    output_file = f'{preseq_dir}/compiled_bulked_preseq.csv'
    df.to_csv(output_file, index=False)
    print(f"Saved compiled results to: {output_file}")

    # Plot
    print("\nGenerating comparison plots...")
    plot_complexity_curves(df, preseq_dir)

    # Print summary
    print_comparison_summary(df)

    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)

if __name__ == '__main__':
    main()
