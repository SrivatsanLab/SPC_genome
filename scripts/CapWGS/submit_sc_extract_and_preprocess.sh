#!/bin/bash
#SBATCH -J submit_sc_extract_preprocess
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH -t 5

##########################################################################################################################
# Submit unified single-cell extraction and preprocessing array job
# Reads cell count from real_cells.txt and submits array that does both extraction and preprocessing
##########################################################################################################################

set -euo pipefail

BULK_BAM="$1"
BARCODE_FILE="$2"
OUTPUT_DIR="$3"
REFERENCE_DIR="$4"
SCRIPTS_DIR="$5"

echo "=========================================="
echo "Single-Cell Extract & Preprocess Submission"
echo "=========================================="
echo "Bulk BAM: ${BULK_BAM}"
echo "Barcode file: ${BARCODE_FILE}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Reference directory: ${REFERENCE_DIR}"
echo "=========================================="
echo ""

# Verify inputs exist
if [ ! -f "${BULK_BAM}" ]; then
    echo "ERROR: Bulk BAM not found: ${BULK_BAM}"
    exit 1
fi

if [ ! -f "${BARCODE_FILE}" ]; then
    echo "ERROR: Barcode file not found: ${BARCODE_FILE}"
    exit 1
fi

# Count cells
CELL_COUNT=$(wc -l < "${BARCODE_FILE}")
echo "Found ${CELL_COUNT} cells to process"
echo ""

if [ "$CELL_COUNT" -eq 0 ]; then
    echo "ERROR: No cells to process"
    exit 1
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Submit unified array job (extraction + preprocessing)
echo "Submitting unified extraction and preprocessing array for ${CELL_COUNT} cells..."

array_job_ID=$(sbatch --parsable \
    --array=1-${CELL_COUNT} \
    "${SCRIPTS_DIR}/scripts/CapWGS/sc_extract_and_preprocess_array.sh" \
    "${BULK_BAM}" \
    "${BARCODE_FILE}" \
    "${OUTPUT_DIR}" \
    "${REFERENCE_DIR}" \
    "${SCRIPTS_DIR}")

echo "Array job submitted: ${array_job_ID}"
echo "Array size: 1-${CELL_COUNT}"
echo ""
echo "=========================================="
echo "Submission complete!"
echo "=========================================="
echo ""
echo "Monitor progress with:"
echo "  squeue -j ${array_job_ID}"
echo ""
echo "Check logs in:"
echo "  SLURM_outs/array_outs/sc_extract_preprocess_${array_job_ID}_*.out"
echo ""
echo "Output BAMs will be in: ${OUTPUT_DIR}"
echo "=========================================="
echo ""
