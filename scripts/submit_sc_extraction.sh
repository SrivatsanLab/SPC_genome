#!/bin/bash
#SBATCH -J submit_sc_extract
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH -t 10

###########################################################################################################################
# Wrapper script that submits single cell BAM extraction array job
# This script runs AFTER concatenation completes and reads real_cells.txt to determine array size
###########################################################################################################################

set -euo pipefail

BULK_BAM="$1"
REAL_CELLS_FILE="$2"
OUTPUT_DIR="$3"

# Verify real_cells.txt exists
if [ ! -f "$REAL_CELLS_FILE" ]; then
    echo "ERROR: real_cells.txt not found at: $REAL_CELLS_FILE"
    exit 1
fi

# Count cells
cell_count=$(wc -l < "$REAL_CELLS_FILE")

echo "Found ${cell_count} cells in ${REAL_CELLS_FILE}"

if [ "$cell_count" -eq 0 ]; then
    echo "ERROR: No cells detected. Cannot submit extraction array."
    exit 1
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Submit single cell extraction array
sc_extract_ID=$(sbatch --parsable \
    --array=1-${cell_count} \
    ./scripts/extract_sc_array.sh \
    "${BULK_BAM}" \
    "${REAL_CELLS_FILE}")

echo "Single cell extraction array submitted: ${sc_extract_ID}"
echo "Array size: 1-${cell_count}"
echo "Output directory: ${OUTPUT_DIR}"
