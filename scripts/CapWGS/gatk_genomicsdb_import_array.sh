#!/bin/bash
#SBATCH -J genomicsdb_interval
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 4:00:00

##########################################################################################################################
# Per-interval GenomicsDB import for CapWGS (SLURM array job)
# Uses balanced genomic intervals (~100 intervals) for optimal parallelization
# Each job creates a GenomicsDB for one genomic interval
##########################################################################################################################

set -euo pipefail

# Arguments
GVCF_DIR="$1"           # Directory containing single-cell GVCFs
GENOMICSDB_DIR="$2"     # Base directory for GenomicsDB workspaces
SAMPLE_MAP="$3"         # Sample map file
REFERENCE_DIR="$4"      # Reference directory
INTERVAL_LIST="$5"      # File with one interval file path per line

# Get the interval file for this array task
INTERVAL_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INTERVAL_LIST")

# Extract interval name from filename (e.g., "0042-scattered" from "0042-scattered.interval_list")
INTERVAL_NAME=$(basename "$INTERVAL_FILE" .interval_list)

echo "============================================"
echo "CapWGS GenomicsDB Import for interval: ${INTERVAL_NAME}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "============================================"
echo ""

# Define paths
REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
GENOMICS_DB="${GENOMICSDB_DIR}/${INTERVAL_NAME}"
TMP_DIR="${GENOMICSDB_DIR}/tmp_${INTERVAL_NAME}"

mkdir -p "${TMP_DIR}"

# Verify sample map exists
if [ ! -f "${SAMPLE_MAP}" ]; then
    echo "ERROR: Sample map not found at ${SAMPLE_MAP}"
    exit 1
fi

n_samples=$(wc -l < "${SAMPLE_MAP}")
echo "Sample map: ${SAMPLE_MAP}"
echo "Cells: ${n_samples}"
echo "Interval: ${INTERVAL_NAME}"
echo "Interval file: ${INTERVAL_FILE}"
echo ""

# Verify interval file exists
if [ ! -f "${INTERVAL_FILE}" ]; then
    echo "ERROR: Interval file not found at ${INTERVAL_FILE}"
    exit 1
fi

# Remove incomplete GenomicsDB if it exists for this interval
if [ -d "${GENOMICS_DB}" ]; then
    if [ ! -f "${GENOMICS_DB}/callset.json" ]; then
        echo "Removing incomplete GenomicsDB for ${INTERVAL_NAME}..."
        chmod -R +w "${GENOMICS_DB}" 2>/dev/null || true
        find "${GENOMICS_DB}" -type f -delete 2>/dev/null || true
        find "${GENOMICS_DB}" -depth -type d -exec rmdir {} \; 2>/dev/null || true
        rm -rf "${GENOMICS_DB}" 2>/dev/null || true
        if [ -d "${GENOMICS_DB}" ]; then
            echo "Could not fully remove GenomicsDB, renaming it..."
            mv "${GENOMICS_DB}" "${GENOMICS_DB}_old_$(date +%s)" || true
        fi
    else
        echo "GenomicsDB already exists for ${INTERVAL_NAME}. Using existing database."
        echo "If you want to recreate it, remove: ${GENOMICS_DB}"
        exit 0
    fi
fi

echo "Starting GenomicsDB import for ${INTERVAL_NAME}..."
echo "Started at: $(date)"
echo ""

module load GATK

# Import GVCFs for this genomic interval only
gatk --java-options "-Xmx56g -Xms56g" GenomicsDBImport \
    --sample-name-map "${SAMPLE_MAP}" \
    --genomicsdb-workspace-path "${GENOMICS_DB}" \
    -R "${REFERENCE}" \
    -L "${INTERVAL_FILE}" \
    --tmp-dir "${TMP_DIR}" \
    --reader-threads 8

module unload GATK

# Clean up temp directory
rm -rf "${TMP_DIR}"

echo ""
echo "GenomicsDB import complete for ${INTERVAL_NAME} at: $(date)"
echo ""
echo "============================================"
echo "GenomicsDB creation successful!"
echo "============================================"
echo ""
echo "Output: ${GENOMICS_DB}"
echo ""
