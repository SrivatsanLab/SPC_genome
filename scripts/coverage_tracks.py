#!/usr/bin/env python3
"""
Generate coverage track plots from bigwig files

Supports:
- Individual contig or whole genome plotting
- Customizable colors and line styles
- Dark mode for presentations
- Flexible output paths
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pyBigWig

mpl.rcParams['agg.path.chunksize'] = 10000


def parse_args():
    parser = argparse.ArgumentParser(description='Generate coverage track plots from bigwig files')

    # Required arguments
    parser.add_argument('bigwig', help='Path to bigwig file')
    parser.add_argument('-o', '--output', required=True, help='Output file path (PNG)')

    # Plotting options
    parser.add_argument('-c', '--contig', default=None,
                        help='Plot specific contig (e.g., chr1). If not specified, plots all contigs')
    parser.add_argument('--color', default='#4A90E2',
                        help='Line color (default: #4A90E2 - blue)')
    parser.add_argument('--linewidth', type=float, default=1.0,
                        help='Line width (default: 1.0)')
    parser.add_argument('--linestyle', default='-',
                        choices=['-', '--', '-.', ':'],
                        help='Line style (default: -)')

    # Display options
    parser.add_argument('--dark-mode', action='store_true',
                        help='Use dark mode (white lines/labels on dark background)')
    parser.add_argument('--log-scale', action='store_true', default=True,
                        help='Use log scale for y-axis (default: True)')
    parser.add_argument('--no-log-scale', dest='log_scale', action='store_false',
                        help='Disable log scale for y-axis')
    parser.add_argument('--ylim', nargs=2, type=float, default=[0, 15],
                        help='Y-axis limits (default: 0 15)')

    # Whole genome plot options
    parser.add_argument('--gap-size', type=int, default=10000000,
                        help='Gap size between contigs in whole genome plot (default: 10000000)')
    parser.add_argument('--bin-size', type=int, default=1000,
                        help='Bin size for coverage averaging (default: 1000)')
    parser.add_argument('--show-contig-labels', action='store_true', default=True,
                        help='Show contig labels (default: True)')
    parser.add_argument('--no-contig-labels', dest='show_contig_labels', action='store_false',
                        help='Hide contig labels')
    parser.add_argument('--show-dividers', action='store_true', default=True,
                        help='Show vertical dividers between contigs (default: True)')
    parser.add_argument('--no-dividers', dest='show_dividers', action='store_false',
                        help='Hide vertical dividers between contigs')

    # Figure options
    parser.add_argument('--width', type=float, default=16,
                        help='Figure width in inches (default: 16)')
    parser.add_argument('--height', type=float, default=2,
                        help='Figure height in inches (default: 2)')
    parser.add_argument('--dpi', type=int, default=600,
                        help='Figure DPI (default: 600)')

    return parser.parse_args()


def setup_dark_mode():
    """Configure matplotlib for dark mode"""
    plt.style.use('dark_background')
    return {
        'text_color': 'white',
        'divider_color': '#555555',
        'spine_color': 'white'
    }


def setup_light_mode():
    """Configure matplotlib for light mode"""
    sns.set_style("white")
    return {
        'text_color': 'black',
        'divider_color': 'gray',
        'spine_color': 'black'
    }


def get_contig_order(bw):
    """
    Determine contig order from bigwig file
    Prioritizes standard chromosomes in numerical order
    """
    chrom_sizes = bw.chroms()
    contigs = list(chrom_sizes.keys())

    # Try to extract numeric chromosomes and sort them
    numeric_chroms = []
    other_chroms = []

    for contig in contigs:
        # Handle various naming conventions
        contig_base = contig.split('_')[0]  # Remove _hg38, _mm10 suffixes

        if contig_base.startswith('chr'):
            chrom_num = contig_base[3:]
            try:
                # Try to convert to int (for chr1, chr2, etc.)
                numeric_chroms.append((int(chrom_num), contig))
            except ValueError:
                # Handle chrX, chrY, chrM
                if chrom_num in ['X', 'Y', 'M', 'MT']:
                    order = {'X': 23, 'Y': 24, 'M': 25, 'MT': 25}
                    numeric_chroms.append((order.get(chrom_num, 99), contig))
                else:
                    other_chroms.append(contig)
        else:
            other_chroms.append(contig)

    # Sort numeric chromosomes and combine with others
    numeric_chroms.sort(key=lambda x: x[0])
    ordered_contigs = [c[1] for c in numeric_chroms] + sorted(other_chroms)

    return ordered_contigs


def plot_single_contig(bw, contig, args, style_config):
    """Plot coverage for a single contig"""
    chrom_len = bw.chroms()[contig]
    coverage = bw.values(contig, 0, chrom_len)

    # Bin the coverage
    binned_positions = np.arange(0, chrom_len, args.bin_size)
    binned_coverage = [
        np.nanmean(coverage[i:i + args.bin_size])
        for i in range(0, len(coverage), args.bin_size)
    ]

    # If there's an extra bin for the remainder, extend binned_positions
    if len(binned_coverage) > len(binned_positions):
        binned_positions = np.append(binned_positions, chrom_len)

    # Create plot
    fig, ax = plt.subplots(1, 1, figsize=(args.width, args.height), dpi=args.dpi)

    ax.plot(binned_positions[:len(binned_coverage)], binned_coverage,
            color=args.color, linewidth=args.linewidth, linestyle=args.linestyle)

    if args.log_scale:
        ax.set_yscale('log')

    ax.set_ylim(args.ylim)
    ax.set_xlabel(contig, fontsize=12, color=style_config['text_color'])
    ax.set_ylabel('Coverage', fontsize=12, color=style_config['text_color'])

    # Style the axes
    ax.tick_params(colors=style_config['text_color'])
    for spine in ax.spines.values():
        spine.set_edgecolor(style_config['spine_color'])

    sns.despine()

    return fig


def plot_all_contigs(bw, args, style_config):
    """Plot coverage for all contigs with gaps between them"""
    chrom_sizes = bw.chroms()
    chrom_order = get_contig_order(bw)

    binned_positions = []
    binned_coverage = []
    chrom_offsets = {}
    curr_offset = 0

    print(f"Plotting {len(chrom_order)} contigs:")

    for contig in chrom_order:
        chrom_len = chrom_sizes[contig]
        print(f"  {contig}: {chrom_len:,} bp")

        coverage = bw.values(contig, 0, chrom_len)

        # Bin the coverage
        binned_pos = np.arange(curr_offset, curr_offset + chrom_len, args.bin_size)
        binned_cov = [
            np.nanmean(coverage[i:i + args.bin_size])
            for i in range(0, len(coverage), args.bin_size)
        ]

        if len(binned_cov) > len(binned_pos):
            binned_pos = np.append(binned_pos, curr_offset + chrom_len)

        chrom_offsets[contig] = (curr_offset, curr_offset + chrom_len)
        curr_offset += chrom_len + args.gap_size

        binned_positions.extend(binned_pos)
        binned_coverage.extend(binned_cov)

    # Create plot
    fig, ax = plt.subplots(1, 1, figsize=(args.width, args.height), dpi=args.dpi)

    ax.plot(binned_positions[:len(binned_coverage)], binned_coverage,
            color=args.color, linewidth=args.linewidth, linestyle=args.linestyle)

    if args.log_scale:
        ax.set_yscale('log')

    ax.set_ylim(args.ylim)
    ax.set_ylabel('Coverage', fontsize=12, color=style_config['text_color'])

    # Add contig dividers
    if args.show_dividers:
        for contig, (start, end) in chrom_offsets.items():
            ax.axvline(x=end + args.gap_size/2,
                      color=style_config['divider_color'],
                      linestyle='--', linewidth=0.5)

    # Add contig labels
    if args.show_contig_labels:
        tick_positions = [(start + end) / 2 for (start, end) in chrom_offsets.values()]
        # Clean up contig names (remove _hg38, _mm10 suffixes)
        tick_labels = [contig.split('_')[0] for contig in chrom_order]
        ax.set_xticks(tick_positions, labels=tick_labels, rotation=45,
                     color=style_config['text_color'])
    else:
        ax.set_xticks([])

    # Style the axes
    ax.tick_params(colors=style_config['text_color'])
    for spine in ax.spines.values():
        spine.set_edgecolor(style_config['spine_color'])

    sns.despine()

    return fig


def main():
    args = parse_args()

    # Setup style
    if args.dark_mode:
        style_config = setup_dark_mode()
    else:
        style_config = setup_light_mode()

    # Open bigwig file
    print(f"Opening bigwig: {args.bigwig}")
    bw = pyBigWig.open(args.bigwig)

    # Generate plot
    if args.contig:
        print(f"Plotting single contig: {args.contig}")
        fig = plot_single_contig(bw, args.contig, args, style_config)
    else:
        print("Plotting all contigs")
        fig = plot_all_contigs(bw, args, style_config)

    bw.close()

    # Save figure
    print(f"Saving to: {args.output}")
    fig.savefig(args.output, bbox_inches='tight',
                facecolor=fig.get_facecolor() if args.dark_mode else 'white')

    print("Done!")


if __name__ == '__main__':
    main()
