#!/bin/bash
#SBATCH -J lorenz
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 2:00:00

##########################################################################################################################
# Lorenz Curve Array Job
#
# Purpose: Generate Lorenz curves and Gini coefficients for downsampled BAMs
# Input:   downsample_tasks.tsv (from create_downsample_tasks.py)
# Output:  Lorenz CSVs and Gini coefficients in data/benchmarking/lorenz/{dataset}/
#
# Usage:
#   NTASKS=$(tail -n +2 downsample_tasks.tsv | wc -l)
#   sbatch --array=1-${NTASKS} scripts/benchmarking/lorenz_array.sh downsample_tasks.tsv
##########################################################################################################################

set -euo pipefail

TASK_FILE="$1"
LORENZ_BASE_DIR="data/benchmarking/lorenz"

if [ ! -f "$TASK_FILE" ]; then
    echo "Error: Task file not found: $TASK_FILE"
    exit 1
fi

TASK_LINE=$(tail -n +2 "$TASK_FILE" | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [ -z "$TASK_LINE" ]; then
    echo "Error: No task found for array task ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Columns: task_id, dataset, cell_id, input_bam, reference, target_depth, depth_label, fraction, output_bam
IFS=$'\t' read -r TASK_ID DATASET CELL_ID INPUT_BAM REFERENCE TARGET_DEPTH DEPTH_LABEL FRACTION OUTPUT_BAM <<< "$TASK_LINE"

echo "==========================================="
echo "Lorenz Task ${TASK_ID}"
echo "==========================================="
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Dataset: ${DATASET}"
echo "Cell: ${CELL_ID}"
echo "Depth: ${DEPTH_LABEL}"
echo "Downsampled BAM: ${OUTPUT_BAM}"
echo ""

if [ ! -f "$OUTPUT_BAM" ]; then
    echo "Error: Downsampled BAM not found: $OUTPUT_BAM"
    exit 1
fi

LORENZ_DIR="${LORENZ_BASE_DIR}/${DATASET}"
mkdir -p "$LORENZ_DIR"

BIGWIG="${LORENZ_DIR}/${CELL_ID}_${DEPTH_LABEL}.bw"
LORENZ_CSV="${LORENZ_DIR}/${CELL_ID}_${DEPTH_LABEL}_lorenz.csv"
GINI_FILE="${LORENZ_DIR}/${CELL_ID}_${DEPTH_LABEL}_gini.txt"

##########################################################################################################################
# Step 1: BAM → bigwig
##########################################################################################################################

echo "Step 1: Generating bigwig (1kb bins)..."

module load deepTools

bamCoverage -b "$OUTPUT_BAM" -o "$BIGWIG" --binSize 1000 -of bigwig -p ${SLURM_CPUS_PER_TASK}

module unload deepTools

echo "✓ Bigwig complete"
echo ""

##########################################################################################################################
# Step 2: Lorenz curve + Gini coefficient
##########################################################################################################################

echo "Step 2: Computing Lorenz curve and Gini coefficient..."

eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

python scripts/utils/lorenz.py "$BIGWIG" -o "$LORENZ_CSV" --gini-output "$GINI_FILE"

echo "✓ Lorenz curve complete"
echo ""

##########################################################################################################################
# Cleanup
##########################################################################################################################

rm "$BIGWIG"
echo "✓ Bigwig removed"
echo ""

echo "==========================================="
echo "✓ Task ${TASK_ID} complete!"
echo "==========================================="
echo ""
