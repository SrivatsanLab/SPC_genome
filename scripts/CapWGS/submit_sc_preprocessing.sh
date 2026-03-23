#!/bin/bash
#SBATCH -J submit_preprocess
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH -t 10

###########################################################################################################################
# Wrapper script that submits single-cell preprocessing array job
# This script runs AFTER single cell BAM extraction completes
#
# Submits array job to run MarkDuplicates and BQSR on each cell in parallel
###########################################################################################################################

set -euo pipefail

REAL_CELLS_FILE="$1"
RAW_BAM_DIR="$2"
PREPROCESSED_BAM_DIR="$3"
REFERENCE_DIR="$4"
SCRIPTS_DIR="$5"

echo "=========================================="
echo "Single-Cell Preprocessing Submission"
echo "=========================================="
echo "Cell list: ${REAL_CELLS_FILE}"
echo "Input BAM directory: ${RAW_BAM_DIR}"
echo "Output BAM directory: ${PREPROCESSED_BAM_DIR}"
echo "Reference directory: ${REFERENCE_DIR}"
echo "=========================================="
echo ""

# Verify real_cells.txt exists
if [ ! -f "$REAL_CELLS_FILE" ]; then
    echo "ERROR: Cell list not found at: $REAL_CELLS_FILE"
    exit 1
fi

# Count cells
cell_count=$(wc -l < "$REAL_CELLS_FILE")

echo "Found ${cell_count} cells in ${REAL_CELLS_FILE}"

if [ "$cell_count" -eq 0 ]; then
    echo "ERROR: No cells detected. Cannot submit preprocessing array."
    exit 1
fi

# Create output directory
mkdir -p "${PREPROCESSED_BAM_DIR}"

# Submit preprocessing array job
echo ""
echo "Submitting preprocessing array for ${cell_count} cells..."

preprocess_array_ID=$(sbatch --parsable \
    --array=1-${cell_count} \
    "${SCRIPTS_DIR}/scripts/CapWGS/sc_preprocessing_array.sh" \
    "${REAL_CELLS_FILE}" \
    "${RAW_BAM_DIR}" \
    "${PREPROCESSED_BAM_DIR}" \
    "${REFERENCE_DIR}" \
    "${SCRIPTS_DIR}")

echo ""
echo "=========================================="
echo "Preprocessing array submitted!"
echo "=========================================="
echo "Job ID: ${preprocess_array_ID}"
echo "Array size: 1-${cell_count}"
echo ""
echo "Waiting for preprocessing array to complete..."
echo "(This wrapper will hold until all ${cell_count} preprocessing jobs finish)"
echo ""

# Wait for the array job to complete using srun with dependency
# This ensures downstream jobs depending on THIS wrapper will wait for preprocessing to finish
srun --dependency=afterok:${preprocess_array_ID} --ntasks=1 --cpus-per-task=1 --mem=100M --time=1 /bin/true

echo ""
echo "=========================================="
echo "Preprocessing array complete!"
echo "=========================================="
echo "Output BAMs in: ${PREPROCESSED_BAM_DIR}"
echo "=========================================="
