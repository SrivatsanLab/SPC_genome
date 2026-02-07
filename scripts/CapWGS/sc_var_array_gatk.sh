#!/bin/bash
#SBATCH -J sc_var_gatk
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 12:00:00

##########################################################################################################################
# GATK HaplotypeCaller for single-cell variant calling (GVCF mode)
# This script is used with preprocessed single-cell BAMs (after mark duplicates and BQSR)
#
# Steps:
# 1. Add read groups to BAM (required by GATK)
# 2. Call variants with GATK HaplotypeCaller in GVCF mode
##########################################################################################################################

set -euo pipefail

BARCODE_FILE="$1"
REFERENCE_DIR="$2"
SC_OUTPUTS_DIR="$3"

# Get the barcode for this array task
barcode=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BARCODE_FILE")

if [ -z "$barcode" ]; then
    echo "ERROR: Could not read barcode for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "=========================================="
echo "GATK HaplotypeCaller for cell: ${barcode}"
echo "=========================================="

# Define input/output files
SC_BAM="${SC_OUTPUTS_DIR}/${barcode}.bam"
GATK_BAM="${SC_OUTPUTS_DIR}/${barcode}_GATK.bam"
OUTPUT_GVCF="${SC_OUTPUTS_DIR}/${barcode}.g.vcf.gz"

# Verify input BAM exists
if [ ! -f "${SC_BAM}" ]; then
    echo "ERROR: Single-cell BAM not found: ${SC_BAM}"
    exit 1
fi

# Construct reference path
REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"

# Check if reference exists, try alternative naming
if [ ! -f "${REFERENCE}" ]; then
    # Try to find any .fasta or .fa file in the reference directory
    REFERENCE=$(find "${REFERENCE_DIR}" -maxdepth 1 -name "*.fasta" -o -name "*.fa" | head -n 1)
    if [ -z "${REFERENCE}" ] || [ ! -f "${REFERENCE}" ]; then
        echo "ERROR: Reference genome not found in ${REFERENCE_DIR}"
        exit 1
    fi
fi

echo "Input BAM: ${SC_BAM}"
echo "Reference: ${REFERENCE}"
echo "Output GVCF: ${OUTPUT_GVCF}"
echo ""

##########################################################################################################################
# Step 1: Add read groups for GATK
##########################################################################################################################

echo "Step 1: Adding read groups..."

module load SAMtools picard

java -Xmx4g -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I="${SC_BAM}" \
    O="${GATK_BAM}" \
    RGID="${barcode}" \
    RGSM="${barcode}" \
    RGPL=illumina \
    RGLB=lib1 \
    RGPU=unit1

# Index the BAM
samtools index "${GATK_BAM}"

module unload SAMtools picard

echo "Read groups added."
echo ""

##########################################################################################################################
# Step 2: Call variants with GATK HaplotypeCaller (GVCF mode)
##########################################################################################################################

echo "Step 2: Calling variants with GATK HaplotypeCaller..."

module load GATK

gatk --java-options "-Xmx6g" HaplotypeCaller \
    -R "${REFERENCE}" \
    -I "${GATK_BAM}" \
    -O "${OUTPUT_GVCF}" \
    -ERC GVCF

module unload GATK

echo ""
echo "Variant calling complete."
echo "Output GVCF: ${OUTPUT_GVCF}"
echo ""
