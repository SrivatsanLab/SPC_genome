#!/bin/bash
#SBATCH -J submit_sc_unified
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH -t 5

##########################################################################################################################
# Submit unified single-cell processing array (extraction + preprocessing + QC all in one)
##########################################################################################################################

set -euo pipefail

BULK_BAM="$1"
BARCODE_FILE="$2"
SC_OUTPUT_DIR="$3"
QC_METRICS_DIR="$4"
REFERENCE_DIR="$5"
SCRIPTS_DIR="$6"
BINSIZE="${7:-1000}"

echo "=========================================="
echo "Single-Cell Unified Processing Submission"
echo "=========================================="
echo "Bulk BAM: ${BULK_BAM}"
echo "Barcode file: ${BARCODE_FILE}"
echo "SC output directory: ${SC_OUTPUT_DIR}"
echo "QC metrics directory: ${QC_METRICS_DIR}"
echo "Reference directory: ${REFERENCE_DIR}"
echo "Bigwig bin size: ${BINSIZE}bp"
echo "=========================================="
echo ""

# Verify inputs
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

# Create output directories
mkdir -p "${SC_OUTPUT_DIR}" "${QC_METRICS_DIR}"

# Submit unified array job
echo "Submitting unified processing array for ${CELL_COUNT} cells..."
echo "(Each job does: extraction → preprocessing → QC)"
echo ""

array_job_ID=$(sbatch --parsable \
    --array=1-${CELL_COUNT} \
    "${SCRIPTS_DIR}/scripts/CapWGS/sc_extract_preprocess_qc_array.sh" \
    "${BULK_BAM}" \
    "${BARCODE_FILE}" \
    "${SC_OUTPUT_DIR}" \
    "${QC_METRICS_DIR}" \
    "${REFERENCE_DIR}" \
    "${SCRIPTS_DIR}" \
    "${BINSIZE}")

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
echo "  SLURM_outs/array_outs/sc_extract_pp_qc_${array_job_ID}_*.out"
echo ""
echo "Outputs will be in:"
echo "  Preprocessed BAMs: ${SC_OUTPUT_DIR}/*.preprocessed.bam"
echo "  Bigwigs: ${SC_OUTPUT_DIR}/*.bw"
echo "  Lorenz curves: ${SC_OUTPUT_DIR}/*_lorenz.csv"
echo "  QC metrics: ${QC_METRICS_DIR}/"
echo "=========================================="
echo ""
