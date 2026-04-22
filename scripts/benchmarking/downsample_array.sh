#!/bin/bash
#SBATCH -J downsample
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 2:00:00

##########################################################################################################################
# Downsampling Array Job
#
# Purpose: Downsample BAMs to specified read depths using fixed random seed
# Input:   downsample_tasks.tsv (from create_downsample_tasks.py)
# Output:  Downsampled BAM files in data/benchmarking/downsampled/{dataset}/
#
# Usage:
#   NTASKS=$(tail -n +2 downsample_tasks.tsv | wc -l)
#   sbatch --array=1-${NTASKS} downsample_array.sh downsample_tasks.tsv
##########################################################################################################################

set -euo pipefail

TASK_FILE="$1"

# Fixed random seed for reproducible, nested downsampling
RANDOM_SEED=42

# Validate task file exists
if [ ! -f "$TASK_FILE" ]; then
    echo "Error: Task file not found: $TASK_FILE"
    exit 1
fi

# Get the task for this array index (skip header line)
TASK_LINE=$(tail -n +2 "$TASK_FILE" | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [ -z "$TASK_LINE" ]; then
    echo "Error: No task found for array task ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Parse task line (tab-delimited)
# Columns: task_id, dataset, cell_id, input_bam, reference, target_depth, depth_label, fraction, output_bam
IFS=$'\t' read -r TASK_ID DATASET CELL_ID INPUT_BAM REFERENCE TARGET_DEPTH DEPTH_LABEL FRACTION OUTPUT_BAM <<< "$TASK_LINE"

echo "==========================================="
echo "Downsampling Task ${TASK_ID}"
echo "==========================================="
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Dataset: ${DATASET}"
echo "Cell: ${CELL_ID}"
echo "Target depth: ${DEPTH_LABEL} (${TARGET_DEPTH} reads)"
echo "Fraction: ${FRACTION}"
echo "Input BAM: ${INPUT_BAM}"
echo "Output BAM: ${OUTPUT_BAM}"
echo "Random seed: ${RANDOM_SEED}"
echo ""

# Validate input BAM exists
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM not found: $INPUT_BAM"
    exit 1
fi

# Create output directory
OUTPUT_DIR=$(dirname "$OUTPUT_BAM")
mkdir -p "$OUTPUT_DIR"

# Load samtools
module load SAMtools

# Downsample
echo "Downsampling BAM..."
# Strip leading "0." from fraction if present to avoid "42.0.123" format
FRAC_STRIPPED=$(echo "$FRACTION" | sed 's/^0\.//')
samtools view -s ${RANDOM_SEED}.${FRAC_STRIPPED} -b -o "$OUTPUT_BAM" "$INPUT_BAM"

echo "✓ Downsampling complete"
echo ""

# Index
echo "Indexing downsampled BAM..."
samtools index "$OUTPUT_BAM"

echo "✓ Indexing complete"
echo ""

# Verify output
if [ ! -f "$OUTPUT_BAM" ] || [ ! -f "${OUTPUT_BAM}.bai" ]; then
    echo "Error: Output files not created successfully"
    exit 1
fi

# Count reads in output (for verification)
OUTPUT_READS=$(samtools view -c -F 260 "$OUTPUT_BAM")
echo "Verification:"
echo "  Expected reads: ${TARGET_DEPTH}"
echo "  Actual reads: ${OUTPUT_READS}"
echo "  Difference: $((OUTPUT_READS - TARGET_DEPTH))"
echo ""

echo "==========================================="
echo "✓ Task ${TASK_ID} complete!"
echo "==========================================="
echo ""
