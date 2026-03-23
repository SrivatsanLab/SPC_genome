#!/bin/bash
#SBATCH -J sc_extract_preprocess
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 2:00:00

##########################################################################################################################
# Unified single-cell extraction and preprocessing array job
# For each cell: extracts reads from bulk BAM, then runs MarkDuplicates + BQSR
# This avoids race conditions by combining both steps into a single array job
##########################################################################################################################

set -euo pipefail

BULK_BAM="$1"
BARCODE_FILE="$2"
OUTPUT_DIR="$3"
REFERENCE_DIR="$4"
SCRIPTS_DIR="$5"

# Get the cell barcode for this array task
BARCODE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BARCODE_FILE")

echo "=========================================="
echo "Processing cell: ${BARCODE}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "=========================================="
echo ""

# Define paths
RAW_BAM="${OUTPUT_DIR}/${BARCODE}.bam"
PREPROCESSED_BAM="${OUTPUT_DIR}/${BARCODE}.preprocessed.bam"
REFERENCE="${REFERENCE_DIR}/genome.fa"

# Known sites for BQSR (if available)
KNOWN_SITES=""
if [ -f "${REFERENCE_DIR}/dbsnp.vcf.gz" ]; then
    KNOWN_SITES="--known-sites ${REFERENCE_DIR}/dbsnp.vcf.gz"
fi

##########################################################################################################################
# Step 1: Extract single-cell BAM from bulk BAM
##########################################################################################################################

echo "Step 1: Extracting reads for barcode ${BARCODE}..."
echo "Bulk BAM: ${BULK_BAM}"
echo "Output BAM: ${RAW_BAM}"
echo ""

module load SAMtools

# Extract reads with this barcode and sort
samtools view -h "${BULK_BAM}" | \
    grep -E "^@|CB:Z:${BARCODE}$" | \
    samtools sort -@ 4 -o "${RAW_BAM}"

# Index the raw BAM
samtools index "${RAW_BAM}"

echo "Extraction complete: ${RAW_BAM}"
echo ""

##########################################################################################################################
# Step 2: Run MarkDuplicates
##########################################################################################################################

echo "Step 2: Running MarkDuplicates..."

module load picard

MARKDUP_BAM="${OUTPUT_DIR}/${BARCODE}.markdup.bam"
MARKDUP_METRICS="${OUTPUT_DIR}/${BARCODE}.markdup_metrics.txt"

java -Xmx12g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    INPUT="${RAW_BAM}" \
    OUTPUT="${MARKDUP_BAM}" \
    METRICS_FILE="${MARKDUP_METRICS}" \
    ASSUME_SORT_ORDER=coordinate \
    CREATE_INDEX=true

echo "MarkDuplicates complete"
echo ""

# Clean up: Remove raw BAM to save disk space
echo "Removing raw BAM (no longer needed after MarkDuplicates)..."
rm "${RAW_BAM}" "${RAW_BAM}.bai"

##########################################################################################################################
# Step 3: Run BQSR (if known sites available)
##########################################################################################################################

if [ -n "${KNOWN_SITES}" ]; then
    echo "Step 3: Running Base Quality Score Recalibration..."

    module load GATK

    RECAL_TABLE="${OUTPUT_DIR}/${BARCODE}.recal_data.table"

    # Generate recalibration table
    gatk BaseRecalibrator \
        -R "${REFERENCE}" \
        -I "${MARKDUP_BAM}" \
        ${KNOWN_SITES} \
        -O "${RECAL_TABLE}"

    # Apply recalibration
    gatk ApplyBQSR \
        -R "${REFERENCE}" \
        -I "${MARKDUP_BAM}" \
        --bqsr-recal-file "${RECAL_TABLE}" \
        -O "${PREPROCESSED_BAM}"

    # Clean up intermediate files
    rm "${MARKDUP_BAM}" "${MARKDUP_BAM}.bai" "${RECAL_TABLE}"

    echo "BQSR complete"
else
    echo "Step 3: Skipping BQSR (no known sites available)"
    echo "Using MarkDuplicates output as final preprocessed BAM"

    # Just rename markdup BAM to preprocessed BAM
    mv "${MARKDUP_BAM}" "${PREPROCESSED_BAM}"
    mv "${MARKDUP_BAM}.bai" "${PREPROCESSED_BAM}.bai"
fi

echo ""

##########################################################################################################################
# Summary
##########################################################################################################################

echo "=========================================="
echo "Processing complete for ${BARCODE}!"
echo "=========================================="
echo ""
echo "Output: ${PREPROCESSED_BAM}"
echo ""

# Verify output exists
if [ ! -f "${PREPROCESSED_BAM}" ]; then
    echo "ERROR: Preprocessed BAM not created"
    exit 1
fi

echo "Done!"
echo ""
