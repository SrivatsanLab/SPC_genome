#!/bin/bash

##########################################################################################################################
# BCFtools variant calling pipeline for bulk sequencing data
# This script submits the complete pipeline:
# 1. Alignment array job (BWA-MEM)
# 2. Variant calling array job (BCFtools mpileup/call)
# 3. Joint calling to merge all samples into multi-sample VCF (BCFtools merge)
##########################################################################################################################

set -euo pipefail

# Function to display help message
show_help() {
    cat << EOF
Usage: $0 [OPTIONS]

This script submits the complete BCFtools variant calling pipeline for bulk sequencing data.

Required arguments:
  -s, --sample-list FILE      Sample list file (one sample name per line)
  -f, --fastq-dir DIR         Directory containing input FASTQ files
  -o, --output-dir DIR        Output directory for BAMs, VCFs, and results
  -n, --cohort-name NAME      Name for this cohort (used in output filenames)

Optional arguments:
  -g, --reference FILE        Reference genome FASTA file for variant calling
                              (default: /shared/biodata/reference/GATK/hg38/Homo_sapiens_assembly38.fasta)
  -r, --reference-bwa FILE    BWA index prefix for alignment
                              (default: /shared/biodata/reference/GATK/hg38/BWAIndex/Homo_sapiens_assembly38.fasta.64)
  -h, --help                  Show this help message and exit

Default configuration:
  Sample list:    data/bulk_spectra_pilot/BCF/sample_list.txt
  FASTQ dir:      /fh/fast/srivatsan_s/pub/projects/00_genome_transcriptome_coassay/res/260111_VH00738_362_AACLC53HV
  Output dir:     data/bulk_spectra_pilot/BCF
  Cohort name:    AAVS_PolE_bcf
  Reference:      /shared/biodata/reference/GATK/hg38/Homo_sapiens_assembly38.fasta
  BWA index:      /shared/biodata/reference/GATK/hg38/BWAIndex/Homo_sapiens_assembly38.fasta.64

Example:
  # Use default configuration
  $0

  # Custom configuration
  $0 -s my_samples.txt -f /path/to/fastqs -o output -n my_cohort

EOF
}

##########################################################################################################################
# Parse command-line arguments
##########################################################################################################################

# Default values
SAMPLE_LIST="data/bulk_spectra_pilot/BCF/sample_list.txt"
FASTQ_DIR="/fh/fast/srivatsan_s/pub/projects/00_genome_transcriptome_coassay/res/260111_VH00738_362_AACLC53HV"
OUTPUT_DIR="data/bulk_spectra_pilot/BCF"
REFERENCE_FA="/shared/biodata/reference/GATK/hg38/Homo_sapiens_assembly38.fasta"
REFERENCE_BWA="/shared/biodata/reference/GATK/hg38/BWAIndex/Homo_sapiens_assembly38.fasta.64"
COHORT_NAME="AAVS_PolE_bcf"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--sample-list)
            SAMPLE_LIST="$2"
            shift 2
            ;;
        -f|--fastq-dir)
            FASTQ_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -n|--cohort-name)
            COHORT_NAME="$2"
            shift 2
            ;;
        -g|--reference)
            REFERENCE_FA="$2"
            shift 2
            ;;
        -r|--reference-bwa)
            REFERENCE_BWA="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

##########################################################################################################################
# Validate inputs
##########################################################################################################################

if [ ! -f "${SAMPLE_LIST}" ]; then
    echo "ERROR: Sample list file not found: ${SAMPLE_LIST}"
    exit 1
fi

if [ ! -d "${FASTQ_DIR}" ]; then
    echo "ERROR: FASTQ directory not found: ${FASTQ_DIR}"
    exit 1
fi

if [ ! -f "${REFERENCE_FA}" ]; then
    echo "ERROR: Reference FASTA not found: ${REFERENCE_FA}"
    exit 1
fi

if [ ! -f "${REFERENCE_BWA}.amb" ]; then
    echo "ERROR: BWA index not found: ${REFERENCE_BWA}"
    echo "Expected files like ${REFERENCE_BWA}.amb, ${REFERENCE_BWA}.bwt, etc."
    exit 1
fi

##########################################################################################################################
# Configuration
##########################################################################################################################

SCRIPT_DIR="scripts/bulk"
ALIGN_SCRIPT="${SCRIPT_DIR}/bcftools_align_array.sh"
VCALL_SCRIPT="${SCRIPT_DIR}/bcftools_call_array.sh"
JOINT_SCRIPT="${SCRIPT_DIR}/bcftools_merge_and_normalize.sh"

# Output paths
BAM_DIR="${OUTPUT_DIR}/bams"
VCF_DIR="${OUTPUT_DIR}/vcfs"
JOINT_VCF="${OUTPUT_DIR}/${COHORT_NAME}_joint.vcf.gz"

# Count samples
N_SAMPLES=$(wc -l < "${SAMPLE_LIST}")

# Create directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${BAM_DIR}"
mkdir -p "${VCF_DIR}"
mkdir -p SLURM_outs

##########################################################################################################################
# Display configuration
##########################################################################################################################

echo "============================================"
echo "BCFtools Variant Calling Pipeline"
echo "============================================"
echo ""
echo "Configuration:"
echo "  Sample list:       ${SAMPLE_LIST}"
echo "  Number of samples: ${N_SAMPLES}"
echo "  FASTQ directory:   ${FASTQ_DIR}"
echo "  Output directory:  ${OUTPUT_DIR}"
echo "  Reference FASTA:   ${REFERENCE_FA}"
echo "  BWA index:         ${REFERENCE_BWA}"
echo "  Cohort name:       ${COHORT_NAME}"
echo ""
echo "Output files:"
echo "  BAM files:       ${BAM_DIR}/"
echo "  VCF files:       ${VCF_DIR}/"
echo "  Joint VCF:       ${JOINT_VCF}"
echo ""

##########################################################################################################################
# Step 1: Submit alignment array job
##########################################################################################################################

echo "============================================"
echo "Step 1: Submitting alignment array job..."
echo "============================================"
echo ""

ALIGN_JOB_ID=$(sbatch --parsable \
    --array=1-${N_SAMPLES} \
    "${ALIGN_SCRIPT}" \
    "${SAMPLE_LIST}" \
    "${FASTQ_DIR}" \
    "${REFERENCE_BWA}" \
    "${BAM_DIR}" \
    "scripts")

echo "Alignment array job submitted: ${ALIGN_JOB_ID}"
echo ""

##########################################################################################################################
# Step 2: Submit variant calling array job (depends on alignment completion)
##########################################################################################################################

echo "============================================"
echo "Step 2: Submitting variant calling array job..."
echo "============================================"
echo ""

VCALL_JOB_ID=$(sbatch --parsable \
    --dependency=afterok:${ALIGN_JOB_ID} \
    --array=1-${N_SAMPLES} \
    "${VCALL_SCRIPT}" \
    "${SAMPLE_LIST}" \
    "${BAM_DIR}" \
    "${VCF_DIR}" \
    "${REFERENCE_FA}")

echo "Variant calling array job submitted: ${VCALL_JOB_ID}"
echo "  (Will run after alignment job ${ALIGN_JOB_ID} completes)"
echo ""

##########################################################################################################################
# Step 3: Submit merge and normalize job (depends on variant calling completion)
##########################################################################################################################

echo "============================================"
echo "Step 3: Submitting merge and normalize job..."
echo "============================================"
echo ""

JOINT_JOB_ID=$(sbatch --parsable \
    --dependency=afterok:${VCALL_JOB_ID} \
    "${JOINT_SCRIPT}" \
    "${VCF_DIR}" \
    "${REFERENCE_FA}" \
    "${JOINT_VCF}" \
    "${COHORT_NAME}")

echo "Joint calling job submitted: ${JOINT_JOB_ID}"
echo "  (Will run after variant calling job ${VCALL_JOB_ID} completes)"
echo ""

##########################################################################################################################
# Summary
##########################################################################################################################

echo "============================================"
echo "Pipeline submission complete!"
echo "============================================"
echo ""
echo "Job dependencies:"
echo "  1. Alignment:       ${ALIGN_JOB_ID} (${N_SAMPLES} tasks)"
echo "  2. Variant calling: ${VCALL_JOB_ID} (${N_SAMPLES} tasks, depends on ${ALIGN_JOB_ID})"
echo "  3. Joint calling:   ${JOINT_JOB_ID} (depends on ${VCALL_JOB_ID})"
echo ""
echo "Monitor progress with:"
echo "  squeue -u \$USER"
echo "  sacct -j ${ALIGN_JOB_ID},${VCALL_JOB_ID},${JOINT_JOB_ID}"
echo ""
echo "Expected outputs:"
echo "  ${N_SAMPLES} BAM files in ${BAM_DIR}/"
echo "  ${N_SAMPLES} VCF files in ${VCF_DIR}/"
echo "  1 merged VCF in ${JOINT_VCF}"
echo ""
echo "Done!"
