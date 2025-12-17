#!/bin/bash
#SBATCH -J submit_bw_lz
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH -t 10

###########################################################################################################################
# Wrapper script that submits bigwig and Lorenz curve generation array job
# This script runs AFTER single cell BAM extraction completes and reads real_cells.txt to determine array size
###########################################################################################################################

set -euo pipefail

REAL_CELLS_FILE="$1"
BAM_DIR="$2"
OUTPUT_DIR="$3"
BINSIZE="${4:-1000}"

# Verify real_cells.txt exists
if [ ! -f "$REAL_CELLS_FILE" ]; then
    echo "ERROR: real_cells.txt not found at: $REAL_CELLS_FILE"
    exit 1
fi

# Count cells
cell_count=$(wc -l < "$REAL_CELLS_FILE")

echo "Found ${cell_count} cells in ${REAL_CELLS_FILE}"

if [ "$cell_count" -eq 0 ]; then
    echo "ERROR: No cells detected. Cannot submit bigwig/lorenz array."
    exit 1
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Submit bigwig and Lorenz curve generation array
bw_lz_ID=$(sbatch --parsable \
    --array=1-${cell_count} \
    ./scripts/generate_bigwig_and_lorenz_array.sh \
    "${REAL_CELLS_FILE}" \
    "${BAM_DIR}" \
    "${OUTPUT_DIR}" \
    "${BINSIZE}")

echo "Bigwig and Lorenz curve array submitted: ${bw_lz_ID}"
echo "Array size: 1-${cell_count}"
echo "Bin size: ${BINSIZE}bp"
echo "Output directory: ${OUTPUT_DIR}"
