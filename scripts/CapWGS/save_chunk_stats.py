#!/usr/bin/env python3
"""
Save per-chunk preprocessing stats as structured JSON.

Called at the end of PP_array.sh after demux, trimming, and alignment.
Reads captured log files from tee and writes a single JSON file.

Usage:
    python save_chunk_stats.py \
        --chunk-id 00 \
        --sample-name C_elegans \
        --reference-genome /path/to/genome \
        --demux-log /tmp/demux_stats_chunk_00.txt \
        --trimming-log /tmp/trimming_stats_chunk_00.txt \
        --bwa-log /tmp/bwa_stats_chunk_00.txt \
        --output results/preprocessing_stats/chunk_00.json
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from preprocessing_parsers import parse_barcode_mapping, parse_trimming_sections, parse_pairs_removed, parse_bwa_stats


def read_file_safe(filepath):
    """Read file contents, returning empty string if file doesn't exist."""
    if not filepath or not os.path.exists(filepath):
        return ''
    with open(filepath) as f:
        return f.read()


def main():
    parser = argparse.ArgumentParser(description='Save per-chunk preprocessing stats as JSON.')
    parser.add_argument('--chunk-id', required=True, help='Chunk identifier')
    parser.add_argument('--sample-name', required=True, help='Sample name')
    parser.add_argument('--reference-genome', required=True, help='Reference genome path')
    parser.add_argument('--demux-log', required=True, help='Path to captured demux output')
    parser.add_argument('--trimming-log', required=True, help='Path to captured Trim Galore output')
    parser.add_argument('--bwa-log', required=True, help='Path to captured BWA stderr')
    parser.add_argument('--output', required=True, help='Output JSON file path')
    args = parser.parse_args()

    demux_text = read_file_safe(args.demux_log)
    trimming_text = read_file_safe(args.trimming_log)
    bwa_text = read_file_safe(args.bwa_log)

    # Parse each section
    demux = parse_barcode_mapping(demux_text)
    trimming_sections = parse_trimming_sections(trimming_text)
    pairs = parse_pairs_removed(trimming_text)
    alignment = parse_bwa_stats(bwa_text)

    trimming = {}
    if len(trimming_sections) >= 1:
        trimming['r1'] = trimming_sections[0]
    if len(trimming_sections) >= 2:
        trimming['r2'] = trimming_sections[1]
    trimming.update(pairs)

    stats = {
        'chunk_id': args.chunk_id,
        'sample_name': args.sample_name,
        'reference_genome': args.reference_genome,
        'demux': demux,
        'trimming': trimming,
        'alignment': alignment,
    }

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, 'w') as f:
        json.dump(stats, f, indent=2)


if __name__ == '__main__':
    main()
