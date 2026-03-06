#!/bin/bash

##########################################################################################################################
# Master script for parallelized joint calling by chromosome
# This script:
# 1. Creates/validates GenomicsDB if needed
# 2. Generates chromosome list
# 3. Submits array job for per-chromosome joint calling
# 4. Submits merge job to combine per-chromosome VCFs
##########################################################################################################################

set -euo pipefail

# Configuration
COHORT_NAME="K562_mut_accumulation"
OUTPUT_DIR="data/${COHORT_NAME}"
REFERENCE_DIR="/shared/biodata/reference/GATK/hg38"

# Script paths
SCRIPT_DIR="scripts/bulk"
ARRAY_SCRIPT="${SCRIPT_DIR}/gatk_joint_calling_parallel_array.sh"
MERGE_SCRIPT="${SCRIPT_DIR}/gatk_joint_calling_parallel_merge.sh"

# Define paths
REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
INTERVALS="${REFERENCE_DIR}/wgs_calling_regions.hg38.interval_list"
GVCF_DIR="${OUTPUT_DIR}/gvcfs"
GENOMICS_DB="${OUTPUT_DIR}/vcfs/genomicsdb_${COHORT_NAME}"
SAMPLE_MAP="${OUTPUT_DIR}/vcfs/${COHORT_NAME}_sample_map.txt"
CHROMOSOME_FILE="${OUTPUT_DIR}/vcfs/chromosome_list.txt"
PER_CHR_DIR="${OUTPUT_DIR}/vcfs/per_chromosome"
FINAL_VCF="${OUTPUT_DIR}/vcfs/${COHORT_NAME}_joint.vcf.gz"
TMP_DIR="${OUTPUT_DIR}/vcfs/tmp"

echo "============================================"
echo "Parallelized Joint Calling for: ${COHORT_NAME}"
echo "============================================"
echo ""
echo "Configuration:"
echo "  GenomicsDB: ${GENOMICS_DB}"
echo "  Per-chromosome VCFs: ${PER_CHR_DIR}"
echo "  Final VCF: ${FINAL_VCF}"
echo ""

mkdir -p "${PER_CHR_DIR}"
mkdir -p "${TMP_DIR}"

##########################################################################################################################
# Step 1: Create/validate GenomicsDB
##########################################################################################################################

echo "Step 1: Checking GenomicsDB..."

if [ ! -f "${GENOMICS_DB}/callset.json" ]; then
    echo "GenomicsDB not found or incomplete. Creating GenomicsDB..."
    echo "Started at: $(date)"
    echo ""

    # Remove incomplete GenomicsDB if it exists
    if [ -d "${GENOMICS_DB}" ]; then
        echo "Removing incomplete GenomicsDB..."
        rm -rf "${GENOMICS_DB}"
    fi

    module load GATK

    # Import all GVCFs into GenomicsDB
    gatk --java-options "-Xmx56g -Xms56g" GenomicsDBImport \
        --sample-name-map "${SAMPLE_MAP}" \
        --genomicsdb-workspace-path "${GENOMICS_DB}" \
        -R "${REFERENCE}" \
        -L "${INTERVALS}" \
        --tmp-dir "${TMP_DIR}" \
        --reader-threads 8

    module unload GATK

    echo ""
    echo "GenomicsDB creation complete at: $(date)"
    echo ""
else
    echo "GenomicsDB already exists and is valid."
    echo ""
fi

##########################################################################################################################
# Step 2: Generate chromosome list from GATK interval list
##########################################################################################################################

echo "Step 2: Generating chromosome list from interval list..."

# Extract unique chromosome names from GATK wgs_calling_regions.hg38.interval_list
# This ensures we only call variants on regions recommended by GATK
grep -v "^@" "${INTERVALS}" | cut -f1 | sort -u > "${CHROMOSOME_FILE}"

N_CHROMOSOMES=$(wc -l < "${CHROMOSOME_FILE}")
echo "Extracted ${N_CHROMOSOMES} chromosomes from ${INTERVALS}"
echo ""

##########################################################################################################################
# Step 3: Submit array job for per-chromosome joint calling
##########################################################################################################################

echo "Step 3: Submitting array job for per-chromosome joint calling..."

ARRAY_JOB_ID=$(sbatch --parsable \
    --array=1-${N_CHROMOSOMES} \
    "${ARRAY_SCRIPT}" \
    "${GENOMICS_DB}" \
    "${PER_CHR_DIR}" \
    "${REFERENCE_DIR}" \
    "${CHROMOSOME_FILE}")

echo "Array job submitted: ${ARRAY_JOB_ID}"
echo "  Processing ${N_CHROMOSOMES} chromosomes in parallel"
echo ""

##########################################################################################################################
# Step 4: Submit merge job (with dependency on array completion)
##########################################################################################################################

echo "Step 4: Submitting merge job..."

MERGE_JOB_ID=$(sbatch --parsable \
    --dependency=afterok:${ARRAY_JOB_ID} \
    "${MERGE_SCRIPT}" \
    "${PER_CHR_DIR}" \
    "${FINAL_VCF}" \
    "${CHROMOSOME_FILE}")

echo "Merge job submitted: ${MERGE_JOB_ID}"
echo "  (Will run after array job ${ARRAY_JOB_ID} completes)"
echo ""

##########################################################################################################################
# Summary
##########################################################################################################################

echo "============================================"
echo "Submission complete!"
echo "============================================"
echo ""
echo "Monitor jobs with:"
echo "  squeue -u \$USER"
echo "  sacct -j ${ARRAY_JOB_ID}"
echo "  sacct -j ${MERGE_JOB_ID}"
echo ""
echo "Check logs in:"
echo "  SLURM_outs/array_outs/joint_call_chr_${ARRAY_JOB_ID}_*.out"
echo "  SLURM_outs/merge_joint_vcfs_${MERGE_JOB_ID}.out"
echo ""
echo "Per-chromosome VCFs will be in:"
echo "  ${PER_CHR_DIR}/"
echo ""
echo "Final merged VCF will be:"
echo "  ${FINAL_VCF}"
echo ""
