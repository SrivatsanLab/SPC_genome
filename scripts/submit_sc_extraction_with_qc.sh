#!/bin/bash
#SBATCH -J submit_sc_qc
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH -t 10

###########################################################################################################################
# Wrapper script that submits single cell BAM extraction and downstream QC (bigwig + Lorenz curves)
# This script runs AFTER concatenation completes and reads real_cells.txt to determine array size
###########################################################################################################################

set -euo pipefail

BULK_BAM="$1"
REAL_CELLS_FILE="$2"
OUTPUT_DIR="$3"
BINSIZE="${4:-1000}"  # Default 1kb bins

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

# Determine BAM directory and SC output directory
BAM_DIR=$(dirname "$BULK_BAM")
SC_OUTPUT_DIR="${BAM_DIR}/sc_outputs"

# Create output directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${SC_OUTPUT_DIR}"

# Submit single cell extraction array
sc_extract_ID=$(sbatch --parsable \
    --array=1-${cell_count} \
    ./scripts/extract_sc_array.sh \
    "${BULK_BAM}" \
    "${REAL_CELLS_FILE}")

echo "Single cell extraction array submitted: ${sc_extract_ID}"
echo "Array size: 1-${cell_count}"
echo "Output directory: ${SC_OUTPUT_DIR}"

# Submit bigwig and Lorenz curve generation (depends on sc extraction)
bw_lz_ID=$(sbatch --parsable \
    --dependency=afterok:${sc_extract_ID} \
    --array=1-${cell_count} \
    ./scripts/generate_bigwig_and_lorenz_array.sh \
    "${REAL_CELLS_FILE}" \
    "${SC_OUTPUT_DIR}" \
    "${SC_OUTPUT_DIR}" \
    "${BINSIZE}")

echo "Bigwig and Lorenz curve array submitted: ${bw_lz_ID}"
echo "Bin size: ${BINSIZE}bp"
