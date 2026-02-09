#!/bin/bash
#SBATCH -J markdup_bqsr
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 24:00:00

##########################################################################################################################
# Mark duplicates and perform Base Quality Score Recalibration (BQSR) on bulk BAM
# This follows GATK best practices for data preprocessing before variant calling
#
# Steps:
# 1. Mark duplicates with Picard MarkDuplicates
# 2. Detect known sites files for BQSR (dbSNP, known indels, Mills indels)
# 3. Run BaseRecalibrator and ApplyBQSR (if known sites available)
#
# This script is called by CapWGS_PP.sh after bulk BAM creation and before single-cell extraction
##########################################################################################################################

set -euo pipefail

# Arguments
INPUT_BAM="$1"              # Input bulk BAM file
OUTPUT_DIR="$2"             # Output directory for processed BAM
REFERENCE_INPUT="$3"        # Reference genome (can be directory or FASTA file path)
SCRIPTS_DIR="$4"            # Scripts directory (for sourcing utilities)
OUTPUT_NAME="$5"            # Sample name

echo "=========================================="
echo "Mark Duplicates and BQSR"
echo "=========================================="
echo "Input BAM: ${INPUT_BAM}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Reference input: ${REFERENCE_INPUT}"
echo "Sample name: ${OUTPUT_NAME}"
echo "=========================================="
echo ""

# Verify input BAM exists
if [ ! -f "${INPUT_BAM}" ]; then
    echo "ERROR: Input BAM file not found: ${INPUT_BAM}"
    exit 1
fi

# Create output directories
mkdir -p "${OUTPUT_DIR}"
METRICS_DIR="${OUTPUT_DIR}/metrics"
mkdir -p "${METRICS_DIR}"

# Determine reference FASTA and directory (handle both file and directory inputs)
if [ -f "${REFERENCE_INPUT}" ]; then
    # Reference input is a FASTA file
    REFERENCE="${REFERENCE_INPUT}"
    REFERENCE_DIR=$(dirname "${REFERENCE_INPUT}")
elif [ -d "${REFERENCE_INPUT}" ]; then
    # Reference input is a directory, find FASTA file
    REFERENCE_DIR="${REFERENCE_INPUT}"

    # Try GATK bundle naming first
    if [ -f "${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta" ]; then
        REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
    else
        # Search for any .fasta or .fa file
        REFERENCE=$(find "${REFERENCE_DIR}" -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \) | head -n 1)
        if [ -z "${REFERENCE}" ] || [ ! -f "${REFERENCE}" ]; then
            echo "ERROR: Reference genome not found in ${REFERENCE_DIR}"
            exit 1
        fi
    fi
else
    echo "ERROR: Reference input is neither a file nor directory: ${REFERENCE_INPUT}"
    exit 1
fi

echo "Using reference FASTA: ${REFERENCE}"
echo "Reference directory: ${REFERENCE_DIR}"

# Define output files
MARKDUP_BAM="${OUTPUT_DIR}/${OUTPUT_NAME}.markdup.bam"
MARKDUP_METRICS="${METRICS_DIR}/${OUTPUT_NAME}.markdup_metrics.txt"
RECAL_TABLE="${METRICS_DIR}/${OUTPUT_NAME}.recal_data.table"
BQSR_BAM="${OUTPUT_DIR}/${OUTPUT_NAME}.bqsr.bam"

##########################################################################################################################
# Step 1: Mark Duplicates with Picard
##########################################################################################################################

echo ""
echo "Step 1: Marking duplicates with Picard..."
echo ""

module load picard

java -Xmx48g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I="${INPUT_BAM}" \
    O="${MARKDUP_BAM}" \
    M="${MARKDUP_METRICS}" \
    ASSUME_SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT

module unload picard

echo ""
echo "Duplicate marking complete."
echo "Metrics saved to: ${MARKDUP_METRICS}"
echo ""

# Display duplicate metrics summary
echo "Duplicate Metrics Summary:"
grep -A 2 "LIBRARY" "${MARKDUP_METRICS}" || echo "(Metrics table not found)"
echo ""

##########################################################################################################################
# Step 2: Detect known sites files for BQSR
##########################################################################################################################

echo ""
echo "Step 2: Detecting known sites files for BQSR..."
echo ""

# Source the known sites detection utility
source "${SCRIPTS_DIR}/scripts/utils/detect_known_sites.sh"

# Detect known sites (sets DBSNP, KNOWN_INDELS, MILLS_INDELS, KNOWN_SITES_ARGS)
if detect_known_sites "${REFERENCE_DIR}"; then
    echo "Known sites files found. BQSR will be performed."
    RUN_BQSR=true
else
    echo "No known sites files found. Skipping BQSR."
    echo "Using mark duplicates BAM for downstream analysis."
    RUN_BQSR=false
fi

##########################################################################################################################
# Step 3: Base Quality Score Recalibration (BQSR) - if known sites available
##########################################################################################################################

if [ "$RUN_BQSR" = true ]; then
    echo ""
    echo "Step 3: Base Quality Score Recalibration..."
    echo ""

    module load GATK

    # BaseRecalibrator - generates recalibration table
    echo "Running BaseRecalibrator..."
    echo "Known sites arguments: ${KNOWN_SITES_ARGS}"
    echo ""

    gatk --java-options "-Xmx48g" BaseRecalibrator \
        -R "${REFERENCE}" \
        -I "${MARKDUP_BAM}" \
        ${KNOWN_SITES_ARGS} \
        -O "${RECAL_TABLE}"

    echo ""
    echo "BaseRecalibrator complete."
    echo ""

    # ApplyBQSR - applies recalibration to BAM
    echo "Running ApplyBQSR..."
    echo ""

    gatk --java-options "-Xmx48g" ApplyBQSR \
        -R "${REFERENCE}" \
        -I "${MARKDUP_BAM}" \
        --bqsr-recal-file "${RECAL_TABLE}" \
        -O "${BQSR_BAM}"

    module unload GATK

    echo ""
    echo "BQSR complete."
    echo "Recalibrated BAM: ${BQSR_BAM}"
    echo ""

    # The recalibrated BAM will be used for downstream analysis
    FINAL_BAM="${BQSR_BAM}"
else
    echo ""
    echo "Step 3: Skipping BQSR (no known sites files available)."
    echo ""

    # Use mark duplicates BAM directly
    FINAL_BAM="${MARKDUP_BAM}"
fi

##########################################################################################################################
# Create symlink to final BAM with consistent name
##########################################################################################################################

# Create a symlink named {sample}.preprocessed.bam pointing to the final BAM
# This makes it easier for downstream scripts to find the correct BAM
PREPROCESSED_BAM="${OUTPUT_DIR}/${OUTPUT_NAME}.preprocessed.bam"
PREPROCESSED_BAI="${OUTPUT_DIR}/${OUTPUT_NAME}.preprocessed.bam.bai"

# Remove old symlinks if they exist
if [ -L "${PREPROCESSED_BAM}" ]; then
    rm "${PREPROCESSED_BAM}"
fi
if [ -L "${PREPROCESSED_BAI}" ]; then
    rm "${PREPROCESSED_BAI}"
fi

# Create new symlinks
ln -s "$(basename ${FINAL_BAM})" "${PREPROCESSED_BAM}"
ln -s "$(basename ${FINAL_BAM}).bai" "${PREPROCESSED_BAI}"

##########################################################################################################################
# Summary
##########################################################################################################################

echo ""
echo "=========================================="
echo "Mark Duplicates and BQSR Complete"
echo "=========================================="
echo ""
echo "Output files:"
echo "  Mark duplicates BAM: ${MARKDUP_BAM}"
echo "  Mark duplicates BAM index: ${MARKDUP_BAM}.bai"
echo "  Duplicate metrics: ${MARKDUP_METRICS}"

if [ "$RUN_BQSR" = true ]; then
    echo "  Recalibration table: ${RECAL_TABLE}"
    echo "  Recalibrated BAM: ${BQSR_BAM}"
    echo "  Recalibrated BAM index: ${BQSR_BAM}.bai"
    echo ""
    echo "Final preprocessed BAM: ${BQSR_BAM}"
else
    echo ""
    echo "Final preprocessed BAM: ${MARKDUP_BAM}"
fi

echo ""
echo "Symlink created: ${PREPROCESSED_BAM} -> $(basename ${FINAL_BAM})"
echo ""
echo "Ready for single-cell extraction."
echo ""
