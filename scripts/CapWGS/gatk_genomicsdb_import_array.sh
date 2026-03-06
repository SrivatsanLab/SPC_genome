#!/bin/bash
#SBATCH -J genomicsdb_chr
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 4:00:00

##########################################################################################################################
# Per-chromosome GenomicsDB import for CapWGS (SLURM array job)
# Much faster than importing all chromosomes at once
# Each job creates a GenomicsDB for one chromosome only
##########################################################################################################################

set -euo pipefail

# Arguments
GVCF_DIR="$1"           # Directory containing single-cell GVCFs
GENOMICSDB_DIR="$2"     # Base directory for GenomicsDB workspaces
SAMPLE_MAP="$3"         # Sample map file
REFERENCE_DIR="$4"      # Reference directory
CHROMOSOME_FILE="$5"    # File with one chromosome per line

# Get the chromosome for this array task
CHROMOSOME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHROMOSOME_FILE")

echo "============================================"
echo "CapWGS GenomicsDB Import for: ${CHROMOSOME}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "============================================"
echo ""

# Define paths
REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
GENOMICS_DB="${GENOMICSDB_DIR}/${CHROMOSOME}"
TMP_DIR="${GENOMICSDB_DIR}/tmp_${CHROMOSOME}"

mkdir -p "${TMP_DIR}"

# Verify sample map exists
if [ ! -f "${SAMPLE_MAP}" ]; then
    echo "ERROR: Sample map not found at ${SAMPLE_MAP}"
    exit 1
fi

n_samples=$(wc -l < "${SAMPLE_MAP}")
echo "Sample map: ${SAMPLE_MAP}"
echo "Cells: ${n_samples}"
echo "Chromosome: ${CHROMOSOME}"
echo ""

# Remove incomplete GenomicsDB if it exists for this chromosome
if [ -d "${GENOMICS_DB}" ]; then
    if [ ! -f "${GENOMICS_DB}/callset.json" ]; then
        echo "Removing incomplete GenomicsDB for ${CHROMOSOME}..."
        chmod -R +w "${GENOMICS_DB}" 2>/dev/null || true
        find "${GENOMICS_DB}" -type f -delete 2>/dev/null || true
        find "${GENOMICS_DB}" -depth -type d -exec rmdir {} \; 2>/dev/null || true
        rm -rf "${GENOMICS_DB}" 2>/dev/null || true
        if [ -d "${GENOMICS_DB}" ]; then
            echo "Could not fully remove GenomicsDB, renaming it..."
            mv "${GENOMICS_DB}" "${GENOMICS_DB}_old_$(date +%s)" || true
        fi
    else
        echo "GenomicsDB already exists for ${CHROMOSOME}. Using existing database."
        echo "If you want to recreate it, remove: ${GENOMICS_DB}"
        exit 0
    fi
fi

echo "Starting GenomicsDB import for ${CHROMOSOME}..."
echo "Started at: $(date)"
echo ""

module load GATK

# Import GVCFs for this chromosome only
gatk --java-options "-Xmx56g -Xms56g" GenomicsDBImport \
    --sample-name-map "${SAMPLE_MAP}" \
    --genomicsdb-workspace-path "${GENOMICS_DB}" \
    -R "${REFERENCE}" \
    -L "${CHROMOSOME}" \
    --tmp-dir "${TMP_DIR}" \
    --reader-threads 8

module unload GATK

# Clean up temp directory
rm -rf "${TMP_DIR}"

echo ""
echo "GenomicsDB import complete for ${CHROMOSOME} at: $(date)"
echo ""
echo "============================================"
echo "GenomicsDB creation successful!"
echo "============================================"
echo ""
echo "Output: ${GENOMICS_DB}"
echo ""
