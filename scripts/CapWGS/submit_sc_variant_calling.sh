#!/bin/bash
#SBATCH -J submit_sc_var
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH -t 10

###########################################################################################################################
# Wrapper script that submits single cell variant calling array job
# This script runs AFTER concatenation completes and reads real_cells.txt to determine array size
###########################################################################################################################

set -euo pipefail

TMP_DIR="$1"
REAL_CELLS_FILE="$2"
REFERENCE_GENOME="$3"
OUTPUT_DIR="$4"
SCRIPTS_DIR="$5"

# Verify real_cells.txt exists
if [ ! -f "$REAL_CELLS_FILE" ]; then
    echo "ERROR: real_cells.txt not found at: $REAL_CELLS_FILE"
    exit 1
fi

# Count cells
cell_count=$(wc -l < "$REAL_CELLS_FILE")

echo "Found ${cell_count} cells in ${REAL_CELLS_FILE}"

if [ "$cell_count" -eq 0 ]; then
    echo "ERROR: No cells detected. Cannot submit variant calling array."
    exit 1
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Submit single cell variant calling array
sc_var_array_ID=$(sbatch --parsable \
    --array=1-${cell_count} \
    "${SCRIPTS_DIR}/scripts/sc_var_array.sh" \
    "${TMP_DIR}" \
    "${REAL_CELLS_FILE}" \
    "${REFERENCE_GENOME}" \
    "${OUTPUT_DIR}")

echo "Single cell variant calling array submitted: ${sc_var_array_ID}"
echo "Array size: 1-${cell_count}"
