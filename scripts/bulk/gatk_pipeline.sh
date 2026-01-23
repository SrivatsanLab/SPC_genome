#!/bin/bash

##########################################################################################################################
# Master submission script for bulk alignment and joint variant calling pipeline
# This script:
# 1. Creates necessary directories
# 2. Generates sample list
# 3. Submits array job for per-sample processing
# 4. (Optional) Submits joint calling job with dependency on array completion
##########################################################################################################################

set -euo pipefail

# Configuration
FASTQ_DIR="/fh/fast/srivatsan_s/pub/projects/00_genome_transcriptome_coassay/res/260111_VH00738_362_AACLC53HV"
OUTPUT_DIR="data/bulk_spectra"
REFERENCE_DIR="/shared/biodata/reference/GATK/hg38"
TMP_BASE_DIR="/hpc/temp/srivatsan_s/SPC_genome_preprocessing"
COHORT_NAME="AAVS_PolE_clones"

# Script paths
SCRIPT_DIR="scripts/bulk"
GENERATE_SAMPLES="${SCRIPT_DIR}/generate_sample_list.sh"
ARRAY_SCRIPT="${SCRIPT_DIR}/gatk_align_call_array.sh"
JOINT_SCRIPT="${SCRIPT_DIR}/gatk_joint_calling.sh"

# Output files
SAMPLE_LIST="${OUTPUT_DIR}/sample_list.txt"
JOINT_VCF="${OUTPUT_DIR}/vcfs/${COHORT_NAME}_joint.vcf.gz"

echo "============================================"
echo "Bulk Alignment and Variant Calling Pipeline"
echo "============================================"
echo ""
echo "Configuration:"
echo "  FASTQ directory: ${FASTQ_DIR}"
echo "  Output directory: ${OUTPUT_DIR}"
echo "  Reference: ${REFERENCE_DIR}"
echo "  Temp directory: ${TMP_BASE_DIR}"
echo "  Cohort name: ${COHORT_NAME}"
echo ""

##########################################################################################################################
# Step 1: Create directories
##########################################################################################################################

echo "Step 1: Creating directories..."

# Create temp directory for processing
mkdir -p "${TMP_BASE_DIR}"
echo "  Created: ${TMP_BASE_DIR}"

# Create output directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/bams"
mkdir -p "${OUTPUT_DIR}/gvcfs"
mkdir -p "${OUTPUT_DIR}/vcfs"
mkdir -p "${OUTPUT_DIR}/metrics"
mkdir -p "SLURM_outs/array_outs"
echo "  Created output directories in ${OUTPUT_DIR}"
echo ""

##########################################################################################################################
# Step 2: Generate sample list
##########################################################################################################################

echo "Step 2: Generating sample list..."

bash "${GENERATE_SAMPLES}" "${FASTQ_DIR}" "${SAMPLE_LIST}" "*Clone*"

N_SAMPLES=$(wc -l < "${SAMPLE_LIST}")
echo ""
echo "Found ${N_SAMPLES} samples"
echo ""

if [ "${N_SAMPLES}" -eq 0 ]; then
    echo "ERROR: No samples found!"
    exit 1
fi

##########################################################################################################################
# Step 3: Submit array job for per-sample processing
##########################################################################################################################

echo "Step 3: Submitting array job for per-sample processing..."

ARRAY_JOB_ID=$(sbatch --parsable \
    --array=1-${N_SAMPLES} \
    "${ARRAY_SCRIPT}" \
    "${SAMPLE_LIST}" \
    "${FASTQ_DIR}" \
    "${OUTPUT_DIR}" \
    "${REFERENCE_DIR}" \
    "${TMP_BASE_DIR}")

echo "Array job submitted: ${ARRAY_JOB_ID}"
echo ""

##########################################################################################################################
# Step 4: Submit joint calling job (with dependency on array completion)
##########################################################################################################################

echo "Step 4: Do you want to submit the joint calling job now? (y/n)"
echo "  (Joint calling will wait for all array jobs to complete)"
echo ""

read -r -p "Submit joint calling? [y/n]: " response

if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
    echo ""
    echo "Submitting joint calling job with dependency on array job ${ARRAY_JOB_ID}..."

    JOINT_JOB_ID=$(sbatch --parsable \
        --dependency=afterok:${ARRAY_JOB_ID} \
        "${JOINT_SCRIPT}" \
        "${OUTPUT_DIR}/gvcfs" \
        "${JOINT_VCF}" \
        "${REFERENCE_DIR}" \
        "${COHORT_NAME}")

    echo "Joint calling job submitted: ${JOINT_JOB_ID}"
    echo "  (Will run after array job ${ARRAY_JOB_ID} completes)"
    echo ""
else
    echo ""
    echo "Joint calling job NOT submitted."
    echo "You can submit it later with:"
    echo ""
    echo "  sbatch ${JOINT_SCRIPT} \\"
    echo "    ${OUTPUT_DIR}/gvcfs \\"
    echo "    ${JOINT_VCF} \\"
    echo "    ${REFERENCE_DIR} \\"
    echo "    ${COHORT_NAME}"
    echo ""
fi

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
echo ""
echo "Check logs in:"
echo "  SLURM_outs/array_outs/"
echo ""
echo "Final outputs will be in:"
echo "  BAMs: ${OUTPUT_DIR}/bams/"
echo "  GVCFs: ${OUTPUT_DIR}/gvcfs/"
echo "  Joint VCF: ${JOINT_VCF}"
echo ""
