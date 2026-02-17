#!/bin/bash
#SBATCH -J sc_from_bam
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1

##########################################################################################################################
# Extract single-cell BAMs from bulk BAM(s)
# This is a unified utility script used by all pipelines (CapWGS, CapGTA, etc.)
#
# Usage: sc_from_bam.sh <bulk_bam> <barcode_file> <output_dir> <scripts_dir> [suffix]
#
# Arguments:
#   bulk_bam      - Bulk BAM file to extract cells from
#   barcode_file  - File with list of real cell barcodes (one per line)
#   output_dir    - Output directory for single-cell BAMs
#   scripts_dir   - Scripts directory (contains scripts/utils/)
#   suffix        - Optional suffix for output filenames (e.g., "_dna", "_rna")
#
# For each detected cell, submits an array job to extract reads with that cell's barcode
# Outputs sorted and indexed BAMs to output_dir
##########################################################################################################################

set -euo pipefail

BULK_BAM="$1"           # Bulk BAM file
BARCODE_FILE="$2"       # File with list of real cell barcodes
SC_OUTPUTS_DIR="$3"     # Output directory for single-cell BAMs
SCRIPTS_DIR="$4"        # Scripts directory
SUFFIX="${5:-}"         # Optional suffix for output filenames

echo "=========================================="
echo "Extracting Single Cells from Bulk BAM"
echo "=========================================="
echo "Bulk BAM: ${BULK_BAM}"
echo "Barcode file: ${BARCODE_FILE}"
echo "Output directory: ${SC_OUTPUTS_DIR}"
if [ -n "$SUFFIX" ]; then
    echo "Output suffix: ${SUFFIX}"
fi
echo "=========================================="
echo ""

# Verify inputs
if [ ! -f "${BULK_BAM}" ]; then
    echo "ERROR: Bulk BAM file not found: ${BULK_BAM}"
    exit 1
fi

if [ ! -f "${BARCODE_FILE}" ]; then
    echo "ERROR: Barcode file not found: ${BARCODE_FILE}"
    exit 1
fi

# Create output directory
mkdir -p "${SC_OUTPUTS_DIR}"

# Count cells
cell_count=$(wc -l < "${BARCODE_FILE}")

echo "Found ${cell_count} cells to extract"
echo ""

if [ "$cell_count" -eq 0 ]; then
    echo "ERROR: No cells to extract. Barcode file is empty."
    exit 1
fi

# Submit array job to extract all cells
echo "Submitting array job to extract ${cell_count} cells..."

if [ -n "$SUFFIX" ]; then
    array_ID=$(sbatch --parsable \
        --array=1-${cell_count} \
        "${SCRIPTS_DIR}/scripts/utils/extract_sc_from_bam_array.sh" \
        "${BULK_BAM}" \
        "${BARCODE_FILE}" \
        "${SC_OUTPUTS_DIR}" \
        "${SUFFIX}")
else
    array_ID=$(sbatch --parsable \
        --array=1-${cell_count} \
        "${SCRIPTS_DIR}/scripts/utils/extract_sc_from_bam_array.sh" \
        "${BULK_BAM}" \
        "${BARCODE_FILE}" \
        "${SC_OUTPUTS_DIR}")
fi

echo "Array job submitted: ${array_ID}"
echo "Array size: 1-${cell_count}"
echo ""
echo "Single-cell BAMs will be written to: ${SC_OUTPUTS_DIR}"
echo ""
