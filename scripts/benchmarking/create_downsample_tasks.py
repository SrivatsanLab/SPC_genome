#!/usr/bin/env python3
"""
Task Expander for Downsampling Analysis

Purpose: Generate downsampling task list with calculated fractions
Input:   bam_manifest.tsv (from create_bam_manifest.sh)
Output:  downsample_tasks.tsv with all (cell × depth) combinations
"""

import argparse
import pandas as pd
import sys
from pathlib import Path


def format_depth_label(depth):
    """Convert depth integer to human-readable label (e.g., 5000000 -> '5M')"""
    if depth >= 1_000_000:
        return f"{depth // 1_000_000}M"
    elif depth >= 1_000:
        return f"{depth // 1_000}K"
    else:
        return str(depth)


def main():
    parser = argparse.ArgumentParser(
        description="Generate downsampling task list from BAM manifest"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input BAM manifest TSV (from create_bam_manifest.sh)"
    )
    parser.add_argument(
        "--depths",
        required=True,
        help="Comma-separated target depths (e.g., '1000000,5000000,10000000')"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output downsample tasks TSV"
    )
    parser.add_argument(
        "--base-dir",
        default="data/benchmarking",
        help="Base directory for output files (default: data/benchmarking)"
    )

    args = parser.parse_args()

    # Parse target depths
    try:
        target_depths = [int(d.strip()) for d in args.depths.split(",")]
    except ValueError:
        print("Error: --depths must be comma-separated integers", file=sys.stderr)
        sys.exit(1)

    # Read manifest
    print("=" * 60)
    print("Downsampling Task Expander")
    print("=" * 60)
    print(f"Input manifest: {args.input}")
    print(f"Target depths: {', '.join(format_depth_label(d) for d in target_depths)}")
    print(f"Output tasks: {args.output}")
    print()

    try:
        manifest = pd.read_csv(args.input, sep="\t")
    except Exception as e:
        print(f"Error reading manifest: {e}", file=sys.stderr)
        sys.exit(1)

    # Validate required columns
    required_cols = ["dataset", "cell_id", "bam_path", "reference", "total_mapped_reads"]
    missing_cols = set(required_cols) - set(manifest.columns)
    if missing_cols:
        print(f"Error: Manifest missing columns: {missing_cols}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(manifest)} cells from manifest")
    print()

    # Generate tasks
    tasks = []
    skipped_count = 0

    for _, row in manifest.iterrows():
        dataset = row["dataset"]
        cell_id = row["cell_id"]
        input_bam = row["bam_path"]
        reference = row["reference"]
        total_reads = row["total_mapped_reads"]

        for target_depth in target_depths:
            # Calculate fraction
            fraction = target_depth / total_reads

            # Skip if insufficient reads
            if fraction > 1.0:
                skipped_count += 1
                continue

            # Generate output path
            depth_label = format_depth_label(target_depth)
            output_bam = f"{args.base_dir}/downsampled/{dataset}/{cell_id}_{depth_label}.bam"

            # Add task
            tasks.append({
                "dataset": dataset,
                "cell_id": cell_id,
                "input_bam": input_bam,
                "reference": reference,
                "target_depth": target_depth,
                "depth_label": depth_label,
                "fraction": fraction,
                "output_bam": output_bam,
            })

    # Create DataFrame and add task_id
    tasks_df = pd.DataFrame(tasks)
    tasks_df.insert(0, "task_id", range(1, len(tasks_df) + 1))

    # Write output
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    tasks_df.to_csv(args.output, sep="\t", index=False)

    # Summary
    print("=" * 60)
    print("Task generation complete!")
    print("=" * 60)
    print(f"Total tasks generated: {len(tasks_df)}")
    print(f"Skipped (insufficient reads): {skipped_count}")
    print()
    print("Tasks per dataset:")
    print(tasks_df.groupby("dataset").size())
    print()
    print("Tasks per depth:")
    print(tasks_df.groupby("depth_label").size())
    print()
    print(f"Output written to: {args.output}")
    print()

    # Preview
    print("Preview of tasks:")
    print(tasks_df.head(10).to_string(index=False))
    print()


if __name__ == "__main__":
    main()
