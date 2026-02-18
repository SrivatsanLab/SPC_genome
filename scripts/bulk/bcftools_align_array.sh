#!/bin/bash
#SBATCH --job-name=bulk_align_array
#SBATCH --output=SLURM_outs/bulk_align_%A_%a.out
#SBATCH --array=1-17
#SBATCH -c 16
#SBATCH -t 5-00:00:00

######################################################################################################
# Array job to align bulk samples using bcftools_align_single.sh
# Reads sample names from a sample list file
######################################################################################################

set -euo pipefail

# Configuration (passed from bcftools_pipeline.sh or use defaults)
SAMPLE_LIST="${1:-data/bulk_spectra_pilot/BCF/sample_list.txt}"
FASTQ_DIR="${2:-/fh/fast/srivatsan_s/pub/projects/00_genome_transcriptome_coassay/res/260111_VH00738_362_AACLC53HV}"
REFERENCE="${3:-/shared/biodata/reference/GATK/hg38/BWAIndex/Homo_sapiens_assembly38.fasta.64}"
OUTPUT_DIR="${4:-data/bulk_spectra_pilot/BCF/bams}"

# Create output directories
mkdir -p SLURM_outs/
mkdir -p "${OUTPUT_DIR}"

# Get sample name from array index
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

echo "========================================"
echo "Processing sample ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_TASK_MAX}: ${SAMPLE}"
echo "========================================"

# Construct FASTQ file paths
# The FASTQ files follow the pattern: SAMPLE_S##_R1_001.fastq.gz
R1_FILE="${FASTQ_DIR}/${SAMPLE}_S*_R1_001.fastq.gz"
R2_FILE="${FASTQ_DIR}/${SAMPLE}_S*_R2_001.fastq.gz"

# Verify files exist (wildcards need to expand)
R1_EXPANDED=$(ls ${R1_FILE} 2>/dev/null || echo "")
R2_EXPANDED=$(ls ${R2_FILE} 2>/dev/null || echo "")

if [ -z "${R1_EXPANDED}" ] || [ -z "${R2_EXPANDED}" ]; then
    echo "Error: FASTQ files not found for sample ${SAMPLE}"
    echo "  Expected R1: ${R1_FILE}"
    echo "  Expected R2: ${R2_FILE}"
    exit 1
fi

echo "FASTQ files found:"
echo "  R1: ${R1_EXPANDED}"
echo "  R2: ${R2_EXPANDED}"

# Run alignment script (use path relative to submission directory)
SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]}")"
if [[ "$SCRIPT_DIR" == "/var/tmp"* ]]; then
    # SLURM copied script to temp dir, use relative path from project root
    SCRIPT_DIR="scripts/bulk"
fi

bash "${SCRIPT_DIR}/bcftools_align_single.sh" \
    -o "${SAMPLE}" \
    -1 "${R1_EXPANDED}" \
    -2 "${R2_EXPANDED}" \
    -g "${REFERENCE}" \
    -O "${OUTPUT_DIR}"

echo "Completed sample ${SAMPLE}"
