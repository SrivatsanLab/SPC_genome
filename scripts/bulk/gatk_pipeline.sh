#!/bin/bash

##########################################################################################################################
# Master submission script for bulk alignment and joint variant calling pipeline
# This script:
# 1. Creates necessary directories
# 2. Generates sample list
# 3. Submits array job for per-sample processing
# 4. Submits joint calling job with dependency on array completion
#
# Usage:
#   bash gatk_pipeline.sh                         # Default: parallel joint calling
#   bash gatk_pipeline.sh --no-joint-calling      # Skip joint calling
#   bash gatk_pipeline.sh --sequential-joint-calling  # Use sequential method
#
# Options:
#   (none)                       Default: Submit parallel joint calling (by chromosome)
#   --no-joint-calling           Skip joint calling entirely
#   --sequential-joint-calling   Use sequential joint calling (slower, ~3-5 days)
#
# Default: Submits parallel joint calling (recommended, ~12-24 hours)
##########################################################################################################################

set -euo pipefail

# Parse command line arguments
SUBMIT_JOINT_CALLING=true
PARALLEL_METHOD=true

for arg in "$@"; do
    case $arg in
        --no-joint-calling)
            SUBMIT_JOINT_CALLING=false
            shift
            ;;
        --sequential-joint-calling)
            PARALLEL_METHOD=false
            shift
            ;;
        *)
            # Unknown option
            ;;
    esac
done

# Configuration
FASTQ_DIR="/fh/fast/srivatsan_s/SR/ngs/illumina/sanjay/20260130_LH00740_0176_A23F3HLLT4/Unaligned/MutationAccumulationHuman"
OUTPUT_DIR="data/K562_mut_accumulation"
REFERENCE_DIR="/shared/biodata/reference/GATK/hg38"
TMP_BASE_DIR="/hpc/temp/srivatsan_s/SPC_genome_preprocessing"
COHORT_NAME="K562_mut_accumulation"

# Script paths
SCRIPT_DIR="scripts/bulk"
GENERATE_SAMPLES="${SCRIPT_DIR}/generate_sample_list.sh"
ARRAY_SCRIPT="${SCRIPT_DIR}/gatk_align_call_array.sh"
JOINT_SCRIPT="${SCRIPT_DIR}/gatk_joint_calling.sh"
GENOMICSDB_SCRIPT="${SCRIPT_DIR}/gatk_genomicsdb_import.sh"
PARALLEL_ARRAY_SCRIPT="${SCRIPT_DIR}/gatk_joint_calling_parallel_array.sh"
PARALLEL_MERGE_SCRIPT="${SCRIPT_DIR}/gatk_joint_calling_parallel_merge.sh"

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
if [[ "$SUBMIT_JOINT_CALLING" == true ]]; then
    if [[ "$PARALLEL_METHOD" == true ]]; then
        echo "Joint calling: PARALLEL (by chromosome, ~12-24 hours)"
    else
        echo "Joint calling: SEQUENTIAL (~3-5 days)"
    fi
else
    echo "Joint calling: SKIPPED"
fi
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
    "${TMP_BASE_DIR}" \
    ".")

echo "Array job submitted: ${ARRAY_JOB_ID}"
echo ""

##########################################################################################################################
# Step 4: Submit joint calling job (with dependency on array completion)
##########################################################################################################################

if [[ "$SUBMIT_JOINT_CALLING" == true ]]; then
    echo "Step 4: Submitting joint calling job..."
    echo "  (Joint calling will wait for all array jobs to complete)"
    echo ""

    if [[ "$PARALLEL_METHOD" == true ]]; then
        echo ""
        echo "Method: PARALLEL joint calling (interval-based)"
        echo ""

        # Define paths for parallel pipeline
        TMP_DIR="${TMP_BASE_DIR}/${COHORT_NAME}_joint_calling"
        GENOMICS_DB="${OUTPUT_DIR}/vcfs/genomicsdb_${COHORT_NAME}"
        INTERVAL_LIST="${TMP_DIR}/interval_list.txt"
        PER_INTERVAL_DIR="${TMP_DIR}/per_interval_vcfs"
        SCATTER_COUNT=100  # Number of intervals to generate

        mkdir -p "${OUTPUT_DIR}/vcfs"
        mkdir -p "${TMP_DIR}"
        mkdir -p "${PER_INTERVAL_DIR}"

        # Stage 0: Generate balanced genomic intervals
        echo "Stage 0/4: Submitting interval generation job..."
        INTERVAL_GEN_SCRIPT="${SCRIPT_DIR}/generate_intervals.sh"

        # Check if CapWGS generate_intervals.sh exists, use it for bulk too
        if [ ! -f "${INTERVAL_GEN_SCRIPT}" ]; then
            INTERVAL_GEN_SCRIPT="scripts/CapWGS/generate_intervals.sh"
        fi

        INTERVAL_GEN_JOB_ID=$(sbatch --parsable \
            --dependency=afterok:${ARRAY_JOB_ID} \
            "${INTERVAL_GEN_SCRIPT}" \
            "${REFERENCE_DIR}" \
            "${TMP_DIR}" \
            "${SCATTER_COUNT}")
        echo "  Interval generation job: ${INTERVAL_GEN_JOB_ID}"
        echo "  Creating ~${SCATTER_COUNT} balanced intervals"
        echo "  (Will run after per-sample array job ${ARRAY_JOB_ID} completes)"
        echo ""

        # Create sample map (needed for GenomicsDB import)
        SAMPLE_MAP="${OUTPUT_DIR}/vcfs/${COHORT_NAME}_sample_map.txt"
        > "${SAMPLE_MAP}"
        for gvcf in "${OUTPUT_DIR}/gvcfs"/*.g.vcf.gz; do
            if [ -f "$gvcf" ]; then
                sample_name=$(basename "$gvcf" .g.vcf.gz)
                echo -e "${sample_name}\t${gvcf}" >> "${SAMPLE_MAP}"
            fi
        done

        # Stage 1: GenomicsDB import array (depends on interval generation)
        echo "Stage 1/4: Submitting GenomicsDB import array job..."
        GENOMICSDB_ARRAY_SCRIPT="${SCRIPT_DIR}/gatk_genomicsdb_import_array.sh"

        # Check if bulk-specific script exists, otherwise use CapWGS version
        if [ ! -f "${GENOMICSDB_ARRAY_SCRIPT}" ]; then
            GENOMICSDB_ARRAY_SCRIPT="scripts/CapWGS/gatk_genomicsdb_import_array.sh"
        fi

        GENOMICSDB_JOB_ID=$(sbatch --parsable \
            --array=1-${SCATTER_COUNT} \
            --dependency=afterok:${INTERVAL_GEN_JOB_ID} \
            "${GENOMICSDB_ARRAY_SCRIPT}" \
            "${OUTPUT_DIR}/gvcfs" \
            "${GENOMICS_DB}" \
            "${SAMPLE_MAP}" \
            "${REFERENCE_DIR}" \
            "${INTERVAL_LIST}")
        echo "  GenomicsDB import array: ${GENOMICSDB_JOB_ID}"
        echo "  Importing ${SCATTER_COUNT} intervals in parallel"
        echo "  (Will run after interval generation ${INTERVAL_GEN_JOB_ID} completes)"
        echo ""

        # Stage 2: Per-interval joint calling
        echo "Stage 2/4: Submitting per-interval joint calling array job..."
        PARALLEL_ARRAY_JOB_ID=$(sbatch --parsable \
            --array=1-${SCATTER_COUNT} \
            --dependency=afterok:${GENOMICSDB_JOB_ID} \
            "${PARALLEL_ARRAY_SCRIPT}" \
            "${GENOMICS_DB}" \
            "${PER_INTERVAL_DIR}" \
            "${REFERENCE_DIR}" \
            "${INTERVAL_LIST}")
        echo "  Per-interval array job: ${PARALLEL_ARRAY_JOB_ID}"
        echo "  Processing ${SCATTER_COUNT} intervals in parallel"
        echo "  (Will run after GenomicsDB import ${GENOMICSDB_JOB_ID} completes)"
        echo ""

        # Stage 3: Merge per-interval VCFs
        echo "Stage 3/4: Submitting merge job..."
        MERGE_JOB_ID=$(sbatch --parsable \
            --dependency=afterok:${PARALLEL_ARRAY_JOB_ID} \
            "${PARALLEL_MERGE_SCRIPT}" \
            "${PER_INTERVAL_DIR}" \
            "${JOINT_VCF}" \
            "${INTERVAL_LIST}")
        echo "  Merge job: ${MERGE_JOB_ID}"
        echo "  (Will run after all ${SCATTER_COUNT} interval jobs complete)"
        echo ""

        echo "Parallel pipeline submitted successfully!"
        echo ""
        echo "Job chain:"
        echo "  0. Interval generation: ${INTERVAL_GEN_JOB_ID} (~${SCATTER_COUNT} balanced intervals)"
        echo "  1. Per-sample processing: ${ARRAY_JOB_ID}"
        echo "  2. GenomicsDB import: ${GENOMICSDB_JOB_ID} (${SCATTER_COUNT} intervals)"
        echo "  3. Per-interval calling: ${PARALLEL_ARRAY_JOB_ID} (${SCATTER_COUNT} jobs)"
        echo "  4. Merge VCFs: ${MERGE_JOB_ID}"
        echo ""

    else
        echo ""
        echo "Method: SEQUENTIAL joint calling"
        echo ""

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
    fi
else
    echo "Step 4: Skipping joint calling (--no-joint-calling flag set)"
    echo ""
    echo "You can submit it later with:"
    echo ""
    echo "PARALLEL (recommended):"
    echo "  1. sbatch --dependency=afterok:${ARRAY_JOB_ID} ${GENOMICSDB_SCRIPT}"
    echo "  2. Then run: bash scripts/bulk/submit_parallel_joint_calling.sh <GENOMICSDB_JOB_ID>"
    echo ""
    echo "SEQUENTIAL:"
    echo "  sbatch --dependency=afterok:${ARRAY_JOB_ID} ${JOINT_SCRIPT} \\"
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
echo ""
echo "Check logs in:"
echo "  SLURM_outs/array_outs/ (per-sample processing)"
if [[ "$SUBMIT_JOINT_CALLING" == true ]] && [[ "$PARALLEL_METHOD" == true ]]; then
    echo "  SLURM_outs/generate_intervals_${INTERVAL_GEN_JOB_ID}.out (interval generation)"
    echo "  SLURM_outs/array_outs/genomicsdb_interval_${GENOMICSDB_JOB_ID}_*.out (GenomicsDB import)"
    echo "  SLURM_outs/array_outs/jc_interval_${PARALLEL_ARRAY_JOB_ID}_*.out (per-interval calling)"
    echo "  SLURM_outs/merge_joint_vcfs_${MERGE_JOB_ID}.out (merge)"
fi
echo ""
echo "Final outputs will be in:"
echo "  BAMs: ${OUTPUT_DIR}/bams/"
echo "  GVCFs: ${OUTPUT_DIR}/gvcfs/"
if [[ "$SUBMIT_JOINT_CALLING" == true ]]; then
    echo "  Joint VCF: ${JOINT_VCF}"
    if [[ "$PARALLEL_METHOD" == true ]]; then
        echo "  Per-interval VCFs (temp): ${PER_INTERVAL_DIR}/"
        echo "  Intervals (temp): ${TMP_DIR}/intervals/"
    fi
fi
echo ""
