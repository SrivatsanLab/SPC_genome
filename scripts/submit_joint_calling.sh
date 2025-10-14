#!/bin/bash
#SBATCH -J submit_jc
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH -t 30

###########################################################################################################################
# Wrapper script that submits joint calling and h5ad compilation
# This script runs AFTER single-cell variant calling completes
###########################################################################################################################

set -euo pipefail

REFERENCE_GENOME="$1"
BIN_DIR="$2"
ALIGNED_DIR="$3"
RESULTS_DIR="$4"
SCRIPTS_DIR="$5"
OUTPUT_NAME="$6"

# Paths
REFERENCE="${REFERENCE_GENOME}/Homo_sapiens_assembly38.fasta"
INTERVALS_FILE="${BIN_DIR}/intervals.bed"
BARCODES_MAP="${BIN_DIR}/barcodes.map"
TMP_JC_DIR="/hpc/temp/srivatsan_s/joint_calling_${OUTPUT_NAME}"
FINAL_H5AD="${RESULTS_DIR}/${OUTPUT_NAME}.h5ad"

mkdir -p "${TMP_JC_DIR}"
mkdir -p "${RESULTS_DIR}"

# Generate intervals for parallelization (1Mb windows)
echo "Generating genomic intervals..."
bash "${SCRIPTS_DIR}/scripts/make_intervals.sh" "${REFERENCE}" 1000000 "${INTERVALS_FILE}"

interval_count=$(wc -l < "${INTERVALS_FILE}")
echo "Created ${interval_count} intervals"

# Create barcodes.map from VCF files
echo "Creating barcodes.map from VCF files..."
> "${BARCODES_MAP}"  # Clear file
for f in "${ALIGNED_DIR}"/*.g.vcf.gz; do
  if [[ -f "${f}.tbi" ]]; then
    bc=$(basename "$f" .g.vcf.gz)
    echo -e "${bc}\t${f}" >> "${BARCODES_MAP}"
  fi
done

# Verify barcodes.map has entries
if [ ! -f "$BARCODES_MAP" ] || [ ! -s "$BARCODES_MAP" ]; then
    echo "ERROR: barcodes.map not found or empty at: $BARCODES_MAP"
    exit 1
fi

cell_count=$(wc -l < "$BARCODES_MAP")
echo "Found ${cell_count} cells in ${BARCODES_MAP}"

if [ "$cell_count" -eq 0 ]; then
    echo "ERROR: No cells in barcodes.map. Cannot perform joint calling."
    exit 1
fi

# Submit joint calling array job
echo "Submitting joint calling array..."
jc_array_ID=$(sbatch --parsable \
    --array=1-${interval_count} \
    "${SCRIPTS_DIR}/scripts/joint_calling_array.sh" \
    "${TMP_JC_DIR}" \
    "${BARCODES_MAP}" \
    "${INTERVALS_FILE}" \
    "${REFERENCE}" \
    "${SCRIPTS_DIR}")

echo "Joint calling array submitted: ${jc_array_ID}"

# Submit h5ad compilation job (after joint calling completes)
echo "Submitting h5ad compilation..."
compile_ID=$(sbatch --parsable \
    --dependency=afterok:${jc_array_ID} \
    "${SCRIPTS_DIR}/scripts/compile_h5ad.sh" \
    "${TMP_JC_DIR}" \
    "${FINAL_H5AD}")

echo "H5ad compilation job submitted: ${compile_ID}"
echo "Final h5ad will be saved to: ${FINAL_H5AD}"
