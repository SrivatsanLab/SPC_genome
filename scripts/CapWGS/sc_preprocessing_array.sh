#!/bin/bash
#SBATCH -J sc_preprocess
#SBATCH --output=SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 4:00:00

###########################################################################################################################
# Single-cell preprocessing: MarkDuplicates and BQSR
# Runs on individual single-cell BAMs in parallel (GATK mode only)
#
# This array job performs:
# 1. MarkDuplicates with Picard
# 2. Base Quality Score Recalibration (BQSR) with GATK (if known sites available)
#
# Outputs: Preprocessed BAM ready for variant calling
###########################################################################################################################

set -euo pipefail

SAMPLE_LIST="$1"
BAM_DIR="$2"
OUTPUT_DIR="$3"
REFERENCE_DIR="$4"
SCRIPTS_DIR="$5"

# Get cell barcode from array task ID
BARCODE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

# Input/output files
INPUT_BAM="${BAM_DIR}/${BARCODE}.bam"
MARKDUP_BAM="${OUTPUT_DIR}/${BARCODE}.markdup.bam"
MARKDUP_METRICS="${OUTPUT_DIR}/${BARCODE}.markdup_metrics.txt"
RECAL_TABLE="${OUTPUT_DIR}/${BARCODE}.recal_data.table"
FINAL_BAM="${OUTPUT_DIR}/${BARCODE}.preprocessed.bam"

REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"

echo "=========================================="
echo "Preprocessing cell: ${BARCODE}"
echo "=========================================="
echo "Input BAM: ${INPUT_BAM}"
echo "Output BAM: ${FINAL_BAM}"
echo ""

# Verify input BAM exists
if [ ! -f "${INPUT_BAM}" ]; then
    echo "ERROR: Input BAM not found: ${INPUT_BAM}"
    exit 1
fi

##########################################################################################################################
# Step 1: Mark Duplicates
##########################################################################################################################

echo "Step 1: Marking duplicates with Picard..."
echo "Started at: $(date)"

module load picard

java -Xmx6g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I="${INPUT_BAM}" \
    O="${MARKDUP_BAM}" \
    M="${MARKDUP_METRICS}" \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=true

echo "MarkDuplicates complete at: $(date)"
echo ""

##########################################################################################################################
# Step 2: Base Quality Score Recalibration (BQSR)
##########################################################################################################################

# Detect known sites for BQSR
source "${SCRIPTS_DIR}/scripts/utils/detect_known_sites.sh"
detect_known_sites "${REFERENCE_DIR}"

if [ -n "${KNOWN_SITES_ARGS}" ]; then
    echo "Step 2: Running BQSR with known sites..."
    echo "Started at: $(date)"

    module load GATK

    # BaseRecalibrator: Generate recalibration table
    echo "  Running BaseRecalibrator..."
    gatk BaseRecalibrator \
        -R "${REFERENCE}" \
        -I "${MARKDUP_BAM}" \
        ${KNOWN_SITES_ARGS} \
        -O "${RECAL_TABLE}"

    # ApplyBQSR: Apply recalibration to BAM
    echo "  Applying BQSR..."
    gatk ApplyBQSR \
        -R "${REFERENCE}" \
        -I "${MARKDUP_BAM}" \
        --bqsr-recal-file "${RECAL_TABLE}" \
        -O "${FINAL_BAM}"

    echo "BQSR complete at: $(date)"
    echo ""

    # Remove intermediate markdup BAM to save space
    echo "Cleaning up intermediate files..."
    rm "${MARKDUP_BAM}" "${MARKDUP_BAM}.bai"

else
    echo "Step 2: Skipping BQSR (no known sites files found)"
    echo ""

    # No BQSR available, just rename markdup BAM as final output
    mv "${MARKDUP_BAM}" "${FINAL_BAM}"
    mv "${MARKDUP_BAM}.bai" "${FINAL_BAM}.bai"
fi

##########################################################################################################################
# Summary
##########################################################################################################################

##########################################################################################################################
# Clean up: Remove raw BAM to save disk space
##########################################################################################################################

echo "Cleaning up raw BAM (no longer needed after preprocessing)..."
rm "${INPUT_BAM}" "${INPUT_BAM}.bai"

##########################################################################################################################
# Summary
##########################################################################################################################

echo "=========================================="
echo "Preprocessing complete for ${BARCODE}!"
echo "=========================================="
echo "Output files:"
echo "  Preprocessed BAM: ${FINAL_BAM}"
echo "  Preprocessed index: ${FINAL_BAM}.bai"
echo "  Duplicate metrics: ${MARKDUP_METRICS}"
if [ -n "${KNOWN_SITES_ARGS}" ]; then
    echo "  BQSR table: ${RECAL_TABLE}"
fi
echo ""
echo "Raw BAM deleted to save disk space"
echo "Finished at: $(date)"
