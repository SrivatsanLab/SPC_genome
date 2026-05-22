#!/bin/bash
#SBATCH -J markdup_wgs
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 6:00:00

##########################################################################################################################
# MarkDuplicates + CollectWgsMetrics Array Job
#
# Purpose: Mark duplicates on already-downsampled BAMs, then re-run WGS metrics
#          (Use when the original BAMs did not have duplicates marked upstream)
# Input:   downsample_tasks.tsv (from create_downsample_tasks.py)
# Output:  Downsampled BAMs are OVERWRITTEN with duplicate-marked versions
#          WGS metrics are OVERWRITTEN in data/benchmarking/qc_metrics/{dataset}/
#          MarkDuplicates metrics stored in data/benchmarking/qc_metrics/{dataset}/
#
# Usage:
#   NTASKS=$(tail -n +2 downsample_tasks.tsv | wc -l)
#   sbatch --array=1-${NTASKS} markdup_and_metrics_array.sh downsample_tasks.tsv
##########################################################################################################################

set -euo pipefail

TASK_FILE="$1"
QC_BASE_DIR="data/benchmarking/qc_metrics"

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
echo "MarkDup + WGS Metrics Task ${TASK_ID}"
echo "==========================================="
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Dataset: ${DATASET}"
echo "Cell: ${CELL_ID}"
echo "Depth: ${DEPTH_LABEL}"
echo "Downsampled BAM: ${OUTPUT_BAM}"
echo "Reference: ${REFERENCE}"
echo ""

if [ ! -f "$OUTPUT_BAM" ]; then
    echo "Error: Downsampled BAM not found: $OUTPUT_BAM"
    exit 1
fi

if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference not found: $REFERENCE"
    exit 1
fi

QC_DIR="${QC_BASE_DIR}/${DATASET}"
mkdir -p "$QC_DIR"

WGS_METRICS="${QC_DIR}/${CELL_ID}_${DEPTH_LABEL}_wgs_metrics.txt"
DUP_METRICS="${QC_DIR}/${CELL_ID}_${DEPTH_LABEL}_duplicate_metrics.txt"
MARKDUP_BAM="${OUTPUT_BAM%.bam}.markdup.bam"

##########################################################################################################################
# Step 1: MarkDuplicates
##########################################################################################################################

echo "Step 1: Running MarkDuplicates..."

module load picard SAMtools

java -Xmx20g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    INPUT="$OUTPUT_BAM" \
    OUTPUT="$MARKDUP_BAM" \
    METRICS_FILE="$DUP_METRICS" \
    ASSUME_SORT_ORDER=coordinate \
    CREATE_INDEX=false

echo "✓ MarkDuplicates complete"
echo ""

# Replace original BAM with marked version
mv "$MARKDUP_BAM" "$OUTPUT_BAM"
samtools index "$OUTPUT_BAM"

# Remove stale picard-style .bai sitting next to the replaced BAM, if any
if [ -f "${OUTPUT_BAM%.bam}.bai" ]; then
    rm "${OUTPUT_BAM%.bam}.bai"
fi

echo "✓ BAM replaced and re-indexed"
echo ""

##########################################################################################################################
# Step 2: CollectWgsMetrics
##########################################################################################################################

echo "Step 2: Running CollectWgsMetrics..."

module load R

java -Xmx20g -jar $EBROOTPICARD/picard.jar CollectWgsMetrics \
    I="$OUTPUT_BAM" \
    O="$WGS_METRICS" \
    R="$REFERENCE"

echo "✓ CollectWgsMetrics complete"
echo ""

if [ ! -f "$WGS_METRICS" ]; then
    echo "Error: Metrics file not created"
    exit 1
fi

echo "Metrics summary:"
grep -A 2 "^GENOME_TERRITORY" "$WGS_METRICS" || true
echo ""

echo "==========================================="
echo "✓ Task ${TASK_ID} complete!"
echo "==========================================="
echo ""
