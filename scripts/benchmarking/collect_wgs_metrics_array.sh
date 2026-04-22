#!/bin/bash
#SBATCH -J wgs_metrics
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 4:00:00

##########################################################################################################################
# CollectWgsMetrics Array Job
#
# Purpose: Collect Picard WGS metrics on downsampled BAM files
# Input:   downsample_tasks.tsv (from create_downsample_tasks.py)
# Output:  Picard WGS metrics in data/benchmarking/qc_metrics/{dataset}/
#
# Usage:
#   NTASKS=$(tail -n +2 downsample_tasks.tsv | wc -l)
#   sbatch --array=1-${NTASKS} collect_wgs_metrics_array.sh downsample_tasks.tsv
##########################################################################################################################

set -euo pipefail

TASK_FILE="$1"
QC_BASE_DIR="data/benchmarking/qc_metrics"

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
echo "CollectWgsMetrics Task ${TASK_ID}"
echo "==========================================="
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Dataset: ${DATASET}"
echo "Cell: ${CELL_ID}"
echo "Depth: ${DEPTH_LABEL}"
echo "Downsampled BAM: ${OUTPUT_BAM}"
echo "Reference: ${REFERENCE}"
echo ""

# Validate downsampled BAM exists
if [ ! -f "$OUTPUT_BAM" ]; then
    echo "Error: Downsampled BAM not found: $OUTPUT_BAM"
    echo "Make sure downsampling array job completed successfully"
    exit 1
fi

# Validate reference exists
if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference not found: $REFERENCE"
    exit 1
fi

# Create output directory
QC_DIR="${QC_BASE_DIR}/${DATASET}"
mkdir -p "$QC_DIR"

# Define output path
WGS_METRICS="${QC_DIR}/${CELL_ID}_${DEPTH_LABEL}_wgs_metrics.txt"

echo "Output metrics: ${WGS_METRICS}"
echo ""

# Load required modules
module load picard R

# Run CollectWgsMetrics
echo "Running CollectWgsMetrics..."
echo ""

java -Xmx14g -jar $EBROOTPICARD/picard.jar CollectWgsMetrics \
    I="$OUTPUT_BAM" \
    O="$WGS_METRICS" \
    R="$REFERENCE"

echo ""
echo "✓ CollectWgsMetrics complete"
echo ""

# Verify output
if [ ! -f "$WGS_METRICS" ]; then
    echo "Error: Metrics file not created"
    exit 1
fi

# Show summary from metrics file
echo "Metrics summary:"
echo ""
grep -A 2 "^GENOME_TERRITORY" "$WGS_METRICS" || true
echo ""

echo "==========================================="
echo "✓ Task ${TASK_ID} complete!"
echo "==========================================="
echo ""
