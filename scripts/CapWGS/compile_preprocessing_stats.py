#!/usr/bin/env python3
"""
Aggregate per-chunk preprocessing JSON files into summary statistics.

Reads chunk_*.json files from a preprocessing_stats directory and produces:
  - preprocessing_summary.json (overall stats)
  - adapter_histogram_r1.csv (combined adapter position histogram for R1)
  - adapter_histogram_r2.csv (combined adapter position histogram for R2)

Usage:
    python compile_preprocessing_stats.py <stats_dir> <output_dir>
"""

import argparse
import glob
import json
import os
import csv
from collections import defaultdict


def load_chunk_stats(stats_dir):
    """Load all chunk_*.json files from stats_dir."""
    pattern = os.path.join(stats_dir, 'chunk_*.json')
    files = sorted(glob.glob(pattern))
    chunks = []
    for f in files:
        with open(f) as fh:
            chunks.append(json.load(fh))
    return chunks


def aggregate_adapter_histogram(chunks, read_key):
    """Sum adapter position histograms across chunks for a given read (r1 or r2).

    Returns dict mapping length -> total count.
    """
    combined = defaultdict(int)
    for chunk in chunks:
        trimming = chunk.get('trimming', {})
        read_data = trimming.get(read_key, {})
        for length, count in read_data.get('adapter_position_histogram', []):
            combined[length] += count
    return dict(sorted(combined.items()))


def compile_summary(chunks):
    """Build summary dict from all chunk stats."""
    # Barcode mapping - weighted average
    total_input_reads = 0
    total_mapped_reads = 0
    for c in chunks:
        demux = c.get('demux', {})
        tr = demux.get('total_reads', 0)
        mr = demux.get('mapped_reads', 0)
        total_input_reads += tr
        total_mapped_reads += mr

    barcode_mapping_rate = total_mapped_reads / total_input_reads if total_input_reads > 0 else 0

    # Trimming stats - weighted averages
    total_trimmed_reads_r1 = 0
    total_with_adapters_r1 = 0
    total_trimmed_reads_r2 = 0
    total_with_adapters_r2 = 0
    total_bp_processed_r1 = 0
    total_bp_quality_trimmed_r1 = 0
    total_bp_written_r1 = 0
    total_bp_processed_r2 = 0
    total_bp_quality_trimmed_r2 = 0
    total_bp_written_r2 = 0
    total_pairs_removed = 0
    total_pairs_analyzed = 0

    for c in chunks:
        trimming = c.get('trimming', {})
        r1 = trimming.get('r1', {})
        r2 = trimming.get('r2', {})

        tr1 = r1.get('total_reads', 0)
        total_trimmed_reads_r1 += tr1
        total_with_adapters_r1 += r1.get('reads_with_adapters', 0)
        total_bp_processed_r1 += r1.get('total_bp_processed', 0)
        total_bp_quality_trimmed_r1 += r1.get('bp_quality_trimmed', 0)
        total_bp_written_r1 += r1.get('bp_written', 0)

        tr2 = r2.get('total_reads', 0)
        total_trimmed_reads_r2 += tr2
        total_with_adapters_r2 += r2.get('reads_with_adapters', 0)
        total_bp_processed_r2 += r2.get('total_bp_processed', 0)
        total_bp_quality_trimmed_r2 += r2.get('bp_quality_trimmed', 0)
        total_bp_written_r2 += r2.get('bp_written', 0)

        total_pairs_removed += trimming.get('pairs_removed', 0)
        total_pairs_analyzed += tr1  # R1 count == pairs count

    # Reference genome and BWA version (from first chunk)
    reference_genome = ''
    bwa_version = ''
    insert_size_mean = None
    insert_size_median = None
    if chunks:
        aln = chunks[0].get('alignment', {})
        reference_genome = aln.get('reference_index_path', '')
        bwa_version = aln.get('bwa_version', '')

    # Aggregate alignment stats across chunks
    total_reads_processed = 0
    for c in chunks:
        aln = c.get('alignment', {})
        total_reads_processed += aln.get('reads_processed', 0)

    # Aggregate insert size stats (average across chunks)
    insert_means = []
    insert_medians = []
    for c in chunks:
        aln = c.get('alignment', {})
        if 'insert_size_mean' in aln:
            insert_means.append(aln['insert_size_mean'])
        if 'insert_size_median' in aln:
            insert_medians.append(aln['insert_size_median'])
    if insert_means:
        insert_size_mean = sum(insert_means) / len(insert_means)
    if insert_medians:
        insert_size_median = sum(insert_medians) / len(insert_medians)

    summary = {
        'n_chunks': len(chunks),
        'demux': {
            'total_input_reads': total_input_reads,
            'total_mapped_reads': total_mapped_reads,
            'barcode_mapping_rate': round(barcode_mapping_rate, 6),
        },
        'trimming': {
            'r1': {
                'total_reads': total_trimmed_reads_r1,
                'reads_with_adapters': total_with_adapters_r1,
                'pct_reads_with_adapters': round(
                    100 * total_with_adapters_r1 / total_trimmed_reads_r1, 2
                ) if total_trimmed_reads_r1 > 0 else 0,
                'total_bp_processed': total_bp_processed_r1,
                'bp_quality_trimmed': total_bp_quality_trimmed_r1,
                'bp_written': total_bp_written_r1,
            },
            'r2': {
                'total_reads': total_trimmed_reads_r2,
                'reads_with_adapters': total_with_adapters_r2,
                'pct_reads_with_adapters': round(
                    100 * total_with_adapters_r2 / total_trimmed_reads_r2, 2
                ) if total_trimmed_reads_r2 > 0 else 0,
                'total_bp_processed': total_bp_processed_r2,
                'bp_quality_trimmed': total_bp_quality_trimmed_r2,
                'bp_written': total_bp_written_r2,
            },
            'pairs_removed': total_pairs_removed,
            'pct_pairs_removed': round(
                100 * total_pairs_removed / total_pairs_analyzed, 2
            ) if total_pairs_analyzed > 0 else 0,
        },
        'alignment': {
            'reference_genome': reference_genome,
            'bwa_version': bwa_version,
            'total_reads_processed': total_reads_processed,
            'insert_size_mean': round(insert_size_mean, 2) if insert_size_mean else None,
            'insert_size_median': round(insert_size_median, 1) if insert_size_median else None,
        },
    }
    return summary


def write_histogram_csv(histogram, output_path):
    """Write adapter position histogram as CSV."""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['length', 'count'])
        for length in sorted(histogram.keys()):
            writer.writerow([length, histogram[length]])


def main():
    parser = argparse.ArgumentParser(
        description='Compile per-chunk preprocessing stats into summary files.'
    )
    parser.add_argument('stats_dir', help='Directory containing chunk_*.json files')
    parser.add_argument('output_dir', help='Directory to write summary files')
    args = parser.parse_args()

    chunks = load_chunk_stats(args.stats_dir)
    if not chunks:
        print(f"Warning: No chunk stats found in {args.stats_dir}")
        return

    print(f"Loaded {len(chunks)} chunk stats files")

    os.makedirs(args.output_dir, exist_ok=True)

    # Compile and write summary
    summary = compile_summary(chunks)
    summary_path = os.path.join(args.output_dir, 'preprocessing_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Written: {summary_path}")

    # Write adapter histograms
    for read_key, suffix in [('r1', 'r1'), ('r2', 'r2')]:
        hist = aggregate_adapter_histogram(chunks, read_key)
        if hist:
            hist_path = os.path.join(args.output_dir, f'adapter_histogram_{suffix}.csv')
            write_histogram_csv(hist, hist_path)
            print(f"Written: {hist_path}")

    print("Preprocessing stats compilation complete.")


if __name__ == '__main__':
    main()
