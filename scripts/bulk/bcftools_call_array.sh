#!/bin/bash
#SBATCH --job-name=bulk_vcall_array
#SBATCH --output=SLURM_outs/bulk_vcall_%A_%a.out
#SBATCH --array=1-17
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 4:00:00

######################################################################################################
# Array job to call variants on bulk samples using bcftools_call_single.sh
# Reads sample names from a sample list file
######################################################################################################

set -euo pipefail

# Configuration (passed from bcftools_pipeline.sh or use defaults)
SAMPLE_LIST="${1:-data/bulk_spectra_pilot/BCF/sample_list.txt}"
BAM_DIR="${2:-data/bulk_spectra_pilot/BCF/bams}"
VCF_DIR="${3:-data/bulk_spectra_pilot/BCF/vcfs}"
REFERENCE="${4:-/shared/biodata/reference/GATK/hg38/Homo_sapiens_assembly38.fasta}"

# Create output directories
mkdir -p SLURM_outs/
mkdir -p "${VCF_DIR}"

# Get sample name from array index
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

echo "========================================"
echo "Processing sample ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_MAX}: ${SAMPLE}"
echo "========================================"

# Construct file paths
BAM_FILE="${BAM_DIR}/${SAMPLE}.bam"
OUTPUT_VCF="${VCF_DIR}/${SAMPLE}.vcf.gz"

# Verify BAM exists
if [ ! -f "${BAM_FILE}" ]; then
    echo "Error: BAM file not found: ${BAM_FILE}"
    exit 1
fi

# Run variant calling script
bash scripts/bulk/bcftools_call_single.sh \
    "${BAM_FILE}" \
    "${REFERENCE}" \
    "${OUTPUT_VCF}"

echo "Completed variant calling for ${SAMPLE}"
