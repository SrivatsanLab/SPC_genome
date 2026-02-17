#!/bin/bash
#SBATCH -J sc_var_bcf
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 6:00:00

##########################################################################################################################
# BCFtools variant calling for single-cell BAMs
# This script calls variants on single-cell BAMs using bcftools mpileup and call
#
# Steps:
# 1. Run bcftools mpileup to generate pileup
# 2. Run bcftools call to call variants
# 3. Index the output VCF
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
echo "BCFtools variant calling for cell: ${barcode}"
echo "=========================================="

# Define input/output files
SC_BAM="${SC_OUTPUTS_DIR}/${barcode}.bam"
OUTPUT_VCF="${SC_OUTPUTS_DIR}/${barcode}.vcf.gz"

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
echo "Output VCF: ${OUTPUT_VCF}"
echo ""

##########################################################################################################################
# Run bcftools mpileup and call
##########################################################################################################################

module load BCFtools SAMtools

echo "Running bcftools mpileup and call..."

bcftools mpileup \
    --threads 2 \
    -f "${REFERENCE}" \
    -a FORMAT/AD,FORMAT/DP \
    -O u \
    "${SC_BAM}" | \
bcftools call \
    --threads 2 \
    -m \
    -A \
    -O z \
    -o "${OUTPUT_VCF}"

# Index VCF
echo "Indexing VCF..."
bcftools index "${OUTPUT_VCF}"

module unload BCFtools SAMtools

echo ""
echo "Variant calling complete for cell: ${barcode}"
echo "Output VCF: ${OUTPUT_VCF}"
echo ""
