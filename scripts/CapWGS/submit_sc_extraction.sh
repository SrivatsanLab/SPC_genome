#!/bin/bash
#SBATCH -J submit_sc_extract
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH -t 5

##########################################################################################################################
# Submit single-cell extraction array job
# This wrapper reads the cell count from real_cells.txt and submits the extraction array
# Runs after concatenation completes, when real_cells.txt is available
##########################################################################################################################

set -euo pipefail

BULK_BAM="$1"
BARCODE_FILE="$2"
SC_OUTPUTS_DIR="$3"
SCRIPTS_DIR="$4"

echo "=========================================="
echo "Single-Cell Extraction Submission"
echo "=========================================="
echo "Bulk BAM: ${BULK_BAM}"
echo "Barcode file: ${BARCODE_FILE}"
echo "Output directory: ${SC_OUTPUTS_DIR}"
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
echo "Found ${CELL_COUNT} cells to extract"
echo ""

if [ "$CELL_COUNT" -eq 0 ]; then
    echo "ERROR: No cells to extract"
    exit 1
fi

# Submit array job
echo "Submitting extraction array for ${CELL_COUNT} cells..."

array_job_ID=$(sbatch --parsable \
    --array=1-${CELL_COUNT} \
    "${SCRIPTS_DIR}/scripts/utils/extract_sc_from_bam_array.sh" \
    "${BULK_BAM}" \
    "${BARCODE_FILE}" \
    "${SC_OUTPUTS_DIR}")

echo "Extraction array job submitted: ${array_job_ID}"
echo "Array size: 1-${CELL_COUNT}"
echo ""
echo "Waiting for extraction array to complete..."
echo "(This wrapper will hold until all ${CELL_COUNT} extraction jobs finish)"
echo ""

# Wait for the array job to complete using srun with dependency
# This ensures downstream jobs depending on THIS wrapper will wait for extraction to finish
srun --dependency=afterok:${array_job_ID} --ntasks=1 --cpus-per-task=1 --mem=100M --time=1 /bin/true

echo ""
echo "=========================================="
echo "Extraction array complete!"
echo "=========================================="
echo ""
