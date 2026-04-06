#!/usr/bin/env python3
"""
Retroactively parse PP_array SLURM log files to extract preprocessing statistics.

This produces the same structured output as the pipeline's built-in stats capture,
allowing dashboard generation for runs that were completed before stats capture
was added.

Usage:
    python parse_slurm_logs.py <slurm_log_dir> <results_dir> [--job-id JOBID]

Examples:
    # Auto-detect job ID from log filenames
    python parse_slurm_logs.py SLURM_outs/array_outs results/my_sample/

    # Specify a particular job ID
    python parse_slurm_logs.py SLURM_outs/array_outs results/my_sample/ --job-id 49804952
"""

import argparse
import glob
import json
import os
import re
import sys

# Allow importing from same directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from preprocessing_parsers import parse_full_log
from compile_preprocessing_stats import load_chunk_stats, compile_summary, aggregate_adapter_histogram, write_histogram_csv


def find_pp_array_logs(slurm_dir, job_id=None):
    """Find PP_array log files, optionally filtered by job ID.

    Returns list of (task_id, filepath) tuples sorted by task_id.
    """
    if job_id:
        pattern = os.path.join(slurm_dir, f'PP_array_{job_id}_*.out')
    else:
        pattern = os.path.join(slurm_dir, 'PP_array_*_*.out')

    files = glob.glob(pattern)
    if not files:
        return []

    results = []
    for f in files:
        basename = os.path.basename(f)
        # PP_array_JOBID_TASKID.out
        m = re.match(r'PP_array_(\d+)_(\d+)\.out$', basename)
        if m:
            results.append((m.group(1), m.group(2), f))

    # If no job_id specified and multiple jobs found, pick the one with most logs
    if not job_id and results:
        from collections import Counter
        job_counts = Counter(r[0] for r in results)
        best_job = job_counts.most_common(1)[0][0]
        print(f"Auto-detected job ID: {best_job} ({job_counts[best_job]} log files)")
        results = [(jid, tid, f) for jid, tid, f in results if jid == best_job]

    return sorted(results, key=lambda x: int(x[1]))


def parse_single_log(filepath):
    """Parse one PP_array SLURM log file."""
    with open(filepath) as f:
        text = f.read()
    return parse_full_log(text)


def main():
    parser = argparse.ArgumentParser(
        description='Parse PP_array SLURM logs to extract preprocessing statistics.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('slurm_log_dir', help='Directory containing PP_array_*_*.out files')
    parser.add_argument('results_dir', help='Results directory to write stats to')
    parser.add_argument('--job-id', help='SLURM job ID to filter logs (auto-detected if omitted)')
    args = parser.parse_args()

    logs = find_pp_array_logs(args.slurm_log_dir, args.job_id)
    if not logs:
        print(f"Error: No PP_array log files found in {args.slurm_log_dir}", file=sys.stderr)
        if args.job_id:
            print(f"  (filtered for job ID: {args.job_id})", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(logs)} PP_array log files")

    # Create output directory
    stats_dir = os.path.join(args.results_dir, 'preprocessing_stats')
    os.makedirs(stats_dir, exist_ok=True)

    # Parse each log and write chunk JSON
    success = 0
    failed = 0
    for job_id, task_id, filepath in logs:
        try:
            stats = parse_single_log(filepath)
            stats['chunk_id'] = task_id
            stats['source_log'] = os.path.basename(filepath)

            output_path = os.path.join(stats_dir, f'chunk_{task_id}.json')
            with open(output_path, 'w') as f:
                json.dump(stats, f, indent=2)
            success += 1
        except Exception as e:
            print(f"Warning: Failed to parse {filepath}: {e}", file=sys.stderr)
            failed += 1

    print(f"Parsed {success} logs successfully ({failed} failures)")

    if success == 0:
        print("Error: No logs parsed successfully", file=sys.stderr)
        sys.exit(1)

    # Compile summary
    chunks = load_chunk_stats(stats_dir)
    summary = compile_summary(chunks)

    summary_path = os.path.join(args.results_dir, 'preprocessing_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Written: {summary_path}")

    # Write adapter histograms
    for read_key, suffix in [('r1', 'r1'), ('r2', 'r2')]:
        hist = aggregate_adapter_histogram(chunks, read_key)
        if hist:
            hist_path = os.path.join(args.results_dir, f'adapter_histogram_{suffix}.csv')
            write_histogram_csv(hist, hist_path)
            print(f"Written: {hist_path}")

    print("Done! Preprocessing stats have been extracted from SLURM logs.")


if __name__ == '__main__':
    main()
