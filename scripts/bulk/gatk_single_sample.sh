#!/bin/bash
#SBATCH -J gatk_single_sample
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 3-00:00:00

##########################################################################################################################
# Single-sample bulk WGS variant calling following GATK best practices
#
# This script processes a single bulk WGS sample through the complete GATK workflow:
# 1. Merge fastq files (if multiple lanes)
# 2. Trim adapters with trim_galore
# 3. Align to reference with BWA-MEM
# 4. Mark duplicates with Picard MarkDuplicates
# 5. Perform Base Quality Score Recalibration (BQSR)
# 6. Call variants with GATK HaplotypeCaller (produces final VCF, not GVCF)
#
# Usage:
#   sbatch scripts/bulk/gatk_single_sample.sh \
#       -n SAMPLE_NAME \
#       -1 "path/to/*_R1*.fastq.gz" \
#       -2 "path/to/*_R2*.fastq.gz" \
#       -r /path/to/reference/dir \
#       -o /path/to/output/dir \
#       [-t /path/to/temp/dir] \
#       [-s /path/to/scripts/root]
#
# Arguments:
#   -n    Sample name (used for output file naming)
#   -1    Read 1 fastq file(s) - can be wildcard pattern in quotes
#   -2    Read 2 fastq file(s) - can be wildcard pattern in quotes
#   -r    Reference directory (contains fasta, BWAIndex/, bundle/)
#   -o    Output directory (will create bams/, vcfs/, metrics/ subdirs)
#   -t    Temp directory (optional, default: /hpc/temp/srivatsan_s)
#   -s    Scripts root directory (optional, default: .)
##########################################################################################################################

set -euo pipefail

# Parse command line arguments
while getopts "n:1:2:r:o:t:s:" opt; do
    case $opt in
        n) SAMPLE_NAME="$OPTARG" ;;
        1) R1_PATTERN="$OPTARG" ;;
        2) R2_PATTERN="$OPTARG" ;;
        r) REFERENCE_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) TMP_BASE_DIR="$OPTARG" ;;
        s) SCRIPTS_ROOT="$OPTARG" ;;
        *) echo "Usage: $0 -n SAMPLE -1 R1 -2 R2 -r REF_DIR -o OUT_DIR [-t TMP] [-s SCRIPTS]"; exit 1 ;;
    esac
done

# Check required arguments
if [ -z "${SAMPLE_NAME:-}" ] || [ -z "${R1_PATTERN:-}" ] || [ -z "${R2_PATTERN:-}" ] || \
   [ -z "${REFERENCE_DIR:-}" ] || [ -z "${OUTPUT_DIR:-}" ]; then
    echo "ERROR: Missing required arguments"
    echo "Usage: $0 -n SAMPLE -1 R1 -2 R2 -r REF_DIR -o OUT_DIR [-t TMP] [-s SCRIPTS]"
    exit 1
fi

# Set defaults for optional arguments
TMP_BASE_DIR="${TMP_BASE_DIR:-/hpc/temp/srivatsan_s}"
SCRIPTS_ROOT="${SCRIPTS_ROOT:-.}"

# Activate conda environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate spc_genome

# Set up paths
REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
BWA_INDEX="${REFERENCE_DIR}/BWAIndex/Homo_sapiens_assembly38.fasta.64"

# Create output directories
mkdir -p "${OUTPUT_DIR}/bams"
mkdir -p "${OUTPUT_DIR}/vcfs"
mkdir -p "${OUTPUT_DIR}/metrics"

# Create temp directory
SAMPLE_TMP="${TMP_BASE_DIR}/${SAMPLE_NAME}_$$"
mkdir -p "${SAMPLE_TMP}"

echo "============================================"
echo "GATK Single-Sample Variant Calling Pipeline"
echo "============================================"
echo "Sample: ${SAMPLE_NAME}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Temp directory: ${SAMPLE_TMP}"
echo "Reference: ${REFERENCE}"
echo "BWA index: ${BWA_INDEX}"
echo "============================================"

# Detect known sites files for BQSR
source "${SCRIPTS_ROOT}/scripts/utils/detect_known_sites.sh"
detect_known_sites "${REFERENCE_DIR}"

# Output files
MERGED_R1="${SAMPLE_TMP}/${SAMPLE_NAME}_R1.fastq.gz"
MERGED_R2="${SAMPLE_TMP}/${SAMPLE_NAME}_R2.fastq.gz"
TRIMMED_R1="${SAMPLE_TMP}/${SAMPLE_NAME}_R1_val_1.fq.gz"
TRIMMED_R2="${SAMPLE_TMP}/${SAMPLE_NAME}_R2_val_2.fq.gz"
SAM_FILE="${SAMPLE_TMP}/${SAMPLE_NAME}.sam"
RAW_BAM="${SAMPLE_TMP}/${SAMPLE_NAME}.raw.bam"
SORTED_BAM="${SAMPLE_TMP}/${SAMPLE_NAME}.sorted.bam"
MARKDUP_BAM="${OUTPUT_DIR}/bams/${SAMPLE_NAME}.markdup.bam"
MARKDUP_METRICS="${OUTPUT_DIR}/metrics/${SAMPLE_NAME}.markdup_metrics.txt"
RECAL_TABLE="${OUTPUT_DIR}/metrics/${SAMPLE_NAME}.recal_data.table"
BQSR_BAM="${OUTPUT_DIR}/bams/${SAMPLE_NAME}.bqsr.bam"
OUTPUT_VCF="${OUTPUT_DIR}/vcfs/${SAMPLE_NAME}.vcf.gz"

##########################################################################################################################
# Step 1: Merge fastq files (if multiple lanes)
##########################################################################################################################

echo ""
echo "Step 1: Processing fastq files..."

# Check if pattern matches multiple files
R1_FILES=(${R1_PATTERN})
R2_FILES=(${R2_PATTERN})

if [ ${#R1_FILES[@]} -gt 1 ]; then
    echo "Multiple R1 files found, concatenating..."
    cat ${R1_PATTERN} > "${MERGED_R1}"
    cat ${R2_PATTERN} > "${MERGED_R2}"
    echo "  Merged R1: ${MERGED_R1}"
    echo "  Merged R2: ${MERGED_R2}"
else
    echo "Single pair of fastq files found, creating symlinks..."
    ln -s "$(readlink -f ${R1_PATTERN})" "${MERGED_R1}"
    ln -s "$(readlink -f ${R2_PATTERN})" "${MERGED_R2}"
    echo "  R1: ${MERGED_R1}"
    echo "  R2: ${MERGED_R2}"
fi

##########################################################################################################################
# Step 2: Trim adapters
##########################################################################################################################

echo ""
echo "Step 2: Trimming adapters with Trim Galore..."

module load cutadapt/4.1-GCCcore-11.2.0 Trim_Galore/0.6.7-GCCcore-11.2.0

trim_galore --paired --cores 4 -o "${SAMPLE_TMP}" --illumina --gzip "${MERGED_R1}" "${MERGED_R2}"

module unload cutadapt/4.1-GCCcore-11.2.0 Trim_Galore/0.6.7-GCCcore-11.2.0

# Clean up merged files
rm "${MERGED_R1}" "${MERGED_R2}"

echo "Trimmed files:"
echo "  R1: ${TRIMMED_R1}"
echo "  R2: ${TRIMMED_R2}"

##########################################################################################################################
# Step 3: Align with BWA-MEM
##########################################################################################################################

echo ""
echo "Step 3: Aligning with BWA-MEM..."

module load BWA SAMtools

# Add read group for GATK compatibility
RG="@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:ILLUMINA\tLB:${SAMPLE_NAME}"
bwa mem -t 8 -R "${RG}" -o "${SAM_FILE}" "${BWA_INDEX}" "${TRIMMED_R1}" "${TRIMMED_R2}"

# Clean up trimmed files
rm "${TRIMMED_R1}" "${TRIMMED_R2}"

echo "Alignment complete: ${SAM_FILE}"

##########################################################################################################################
# Step 4: Convert to BAM and sort
##########################################################################################################################

echo ""
echo "Step 4: Converting to BAM and sorting..."

samtools view -@ 8 -bS "${SAM_FILE}" | samtools sort -@ 8 -o "${SORTED_BAM}"

# Clean up SAM file
rm "${SAM_FILE}"

echo "Sorted BAM created: ${SORTED_BAM}"

##########################################################################################################################
# Step 5: Mark duplicates with Picard
##########################################################################################################################

echo ""
echo "Step 5: Marking duplicates with Picard..."

module load picard

java -Xmx48g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I="${SORTED_BAM}" \
    O="${MARKDUP_BAM}" \
    M="${MARKDUP_METRICS}" \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=true

# Clean up sorted BAM
rm "${SORTED_BAM}"

echo "Marked duplicates BAM: ${MARKDUP_BAM}"
echo "Metrics: ${MARKDUP_METRICS}"

##########################################################################################################################
# Step 6: Base Quality Score Recalibration (BQSR)
##########################################################################################################################

if [ -n "${KNOWN_SITES_ARGS}" ]; then
    echo ""
    echo "Step 6: Running BQSR with known sites..."

    module load GATK

    echo "  Running BaseRecalibrator..."
    gatk BaseRecalibrator \
        -R "${REFERENCE}" \
        -I "${MARKDUP_BAM}" \
        ${KNOWN_SITES_ARGS} \
        -O "${RECAL_TABLE}"

    echo "  Applying BQSR..."
    gatk ApplyBQSR \
        -R "${REFERENCE}" \
        -I "${MARKDUP_BAM}" \
        --bqsr-recal-file "${RECAL_TABLE}" \
        -O "${BQSR_BAM}"

    # Update BAM for variant calling
    FINAL_BAM="${BQSR_BAM}"

    echo "BQSR complete: ${BQSR_BAM}"
else
    echo ""
    echo "Step 6: Skipping BQSR (no known sites files found)"
    FINAL_BAM="${MARKDUP_BAM}"
fi

##########################################################################################################################
# Step 7: Call variants with GATK HaplotypeCaller (standard mode, produces final VCF)
##########################################################################################################################

echo ""
echo "Step 7: Calling variants with GATK HaplotypeCaller..."

module load GATK

# Call variants directly without GVCF mode for single samples
gatk HaplotypeCaller \
    -R "${REFERENCE}" \
    -I "${FINAL_BAM}" \
    -O "${OUTPUT_VCF}"

echo "VCF created: ${OUTPUT_VCF}"

##########################################################################################################################
# Step 8: Cleanup temp directory
##########################################################################################################################

echo ""
echo "Step 8: Cleaning up temp directory..."

rm -rf "${SAMPLE_TMP}"

echo ""
echo "============================================"
echo "Single-sample variant calling complete!"
echo "============================================"
echo "Outputs:"
echo "  BAM: ${FINAL_BAM}"
echo "  VCF: ${OUTPUT_VCF}"
echo "  Metrics: ${MARKDUP_METRICS}"
if [ -n "${KNOWN_SITES_ARGS}" ]; then
    echo "  BQSR table: ${RECAL_TABLE}"
fi
echo "============================================"
