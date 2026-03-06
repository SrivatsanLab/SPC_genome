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
        echo "Method: PARALLEL joint calling (by chromosome)"
        echo ""

        # Define paths for parallel pipeline
        GENOMICS_DB="${OUTPUT_DIR}/vcfs/genomicsdb_${COHORT_NAME}"
        CHROMOSOME_FILE="${OUTPUT_DIR}/vcfs/chromosome_list.txt"
        PER_CHR_DIR="${OUTPUT_DIR}/vcfs/per_chromosome"

        mkdir -p "${OUTPUT_DIR}/vcfs"

        # Create chromosome list for parallel processing
        cat > "${CHROMOSOME_FILE}" << 'EOF'
chr1
chr2
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
chr19
chr20
chr21
chr22
chrX
chrY
chrM
EOF

        N_CHROMOSOMES=$(wc -l < "${CHROMOSOME_FILE}")
        mkdir -p "${PER_CHR_DIR}"

        # Stage 1: GenomicsDB import
        echo "Stage 1/3: Submitting GenomicsDB import job..."
        GENOMICSDB_JOB_ID=$(sbatch --parsable \
            --dependency=afterok:${ARRAY_JOB_ID} \
            "${GENOMICSDB_SCRIPT}")
        echo "  GenomicsDB import job: ${GENOMICSDB_JOB_ID}"
        echo "  (Will run after per-sample array job ${ARRAY_JOB_ID} completes)"
        echo ""

        # Stage 2: Per-chromosome joint calling
        echo "Stage 2/3: Submitting per-chromosome joint calling array job..."
        PARALLEL_ARRAY_JOB_ID=$(sbatch --parsable \
            --array=1-${N_CHROMOSOMES} \
            --dependency=afterok:${GENOMICSDB_JOB_ID} \
            "${PARALLEL_ARRAY_SCRIPT}" \
            "${GENOMICS_DB}" \
            "${PER_CHR_DIR}" \
            "${REFERENCE_DIR}" \
            "${CHROMOSOME_FILE}")
        echo "  Per-chromosome array job: ${PARALLEL_ARRAY_JOB_ID}"
        echo "  Processing ${N_CHROMOSOMES} chromosomes in parallel"
        echo "  (Will run after GenomicsDB import ${GENOMICSDB_JOB_ID} completes)"
        echo ""

        # Stage 3: Merge per-chromosome VCFs
        echo "Stage 3/3: Submitting merge job..."
        MERGE_JOB_ID=$(sbatch --parsable \
            --dependency=afterok:${PARALLEL_ARRAY_JOB_ID} \
            "${PARALLEL_MERGE_SCRIPT}" \
            "${PER_CHR_DIR}" \
            "${JOINT_VCF}" \
            "${CHROMOSOME_FILE}")
        echo "  Merge job: ${MERGE_JOB_ID}"
        echo "  (Will run after all ${N_CHROMOSOMES} chromosome jobs complete)"
        echo ""

        echo "Parallel pipeline submitted successfully!"
        echo ""
        echo "Job chain:"
        echo "  1. Per-sample processing: ${ARRAY_JOB_ID}"
        echo "  2. GenomicsDB import: ${GENOMICSDB_JOB_ID}"
        echo "  3. Per-chromosome calling: ${PARALLEL_ARRAY_JOB_ID} (${N_CHROMOSOMES} jobs)"
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
    echo "  SLURM_outs/genomicsdb_import_${GENOMICSDB_JOB_ID}.out (GenomicsDB import)"
    echo "  SLURM_outs/array_outs/joint_call_chr_${PARALLEL_ARRAY_JOB_ID}_*.out (per-chromosome calling)"
    echo "  SLURM_outs/merge_joint_vcfs_${MERGE_JOB_ID}.out (merge)"
fi
echo ""
echo "Final outputs will be in:"
echo "  BAMs: ${OUTPUT_DIR}/bams/"
echo "  GVCFs: ${OUTPUT_DIR}/gvcfs/"
if [[ "$SUBMIT_JOINT_CALLING" == true ]]; then
    echo "  Joint VCF: ${JOINT_VCF}"
    if [[ "$PARALLEL_METHOD" == true ]]; then
        echo "  Per-chromosome VCFs: ${PER_CHR_DIR}/"
    fi
fi
echo ""
