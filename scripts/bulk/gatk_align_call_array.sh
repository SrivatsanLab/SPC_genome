#!/bin/bash
#SBATCH -J bulk_align_var
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 24:00:00

##########################################################################################################################
# Bulk DNA alignment and GVCF calling following GATK best practices
# This script:
# 1. Merges all passage fastqs for a given clone
# 2. Trims adapters with trim_galore
# 3. Aligns to GRCh38 with BWA-MEM
# 4. Marks duplicates with Picard MarkDuplicates
# 5. Performs Base Quality Score Recalibration (BQSR)
# 6. Calls GVCFs with GATK HaplotypeCaller (for joint calling)
##########################################################################################################################

set -euo pipefail

# Arguments
SAMPLE_LIST="$1"        # File with one sample name per line (e.g., AAVS_Clone_4)
FASTQ_DIR="$2"          # Directory containing input fastq files
OUTPUT_DIR="$3"         # Output directory for BAMs and VCFs
REFERENCE_DIR="$4"      # Reference directory (e.g., /shared/biodata/reference/GATK/hg38)
TMP_BASE_DIR="${5:-}"   # Optional: base temp directory (default: /hpc/temp/srivatsan_s)

# Get the sample for this array task
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

# Set default temp directory if not provided
if [ -z "$TMP_BASE_DIR" ]; then
    TMP_BASE_DIR="/hpc/temp/srivatsan_s"
fi

# Create sample-specific temp directory
SAMPLE_TMP="${TMP_BASE_DIR}/${SAMPLE}"
mkdir -p "${SAMPLE_TMP}"

echo "============================================"
echo "Processing sample: ${SAMPLE}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Temp directory: ${SAMPLE_TMP}"
echo "============================================"

# Create output directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/bams"
mkdir -p "${OUTPUT_DIR}/gvcfs"
mkdir -p "${OUTPUT_DIR}/metrics"
mkdir -p SLURM_outs/array_outs

# Define paths
REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
BWA_INDEX="${REFERENCE_DIR}/BWAIndex/Homo_sapiens_assembly38.fasta.64"
DBSNP="${REFERENCE_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf"
KNOWN_INDELS="${REFERENCE_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz"
MILLS_INDELS="${REFERENCE_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

# Output files
MERGED_R1="${SAMPLE_TMP}/${SAMPLE}_R1.fastq.gz"
MERGED_R2="${SAMPLE_TMP}/${SAMPLE}_R2.fastq.gz"
TRIMMED_R1="${SAMPLE_TMP}/${SAMPLE}_R1_val_1.fq.gz"
TRIMMED_R2="${SAMPLE_TMP}/${SAMPLE}_R2_val_2.fq.gz"
SAM_FILE="${SAMPLE_TMP}/${SAMPLE}.sam"
RAW_BAM="${SAMPLE_TMP}/${SAMPLE}.raw.bam"
SORTED_BAM="${SAMPLE_TMP}/${SAMPLE}.sorted.bam"
MARKDUP_BAM="${OUTPUT_DIR}/bams/${SAMPLE}.markdup.bam"
MARKDUP_METRICS="${OUTPUT_DIR}/metrics/${SAMPLE}.markdup_metrics.txt"
RECAL_TABLE="${OUTPUT_DIR}/metrics/${SAMPLE}.recal_data.table"
BQSR_BAM="${OUTPUT_DIR}/bams/${SAMPLE}.bqsr.bam"
OUTPUT_GVCF="${OUTPUT_DIR}/gvcfs/${SAMPLE}.g.vcf.gz"

##########################################################################################################################
# Step 1: Merge fastq files for all passages of this clone
##########################################################################################################################

echo ""
echo "Step 1: Merging fastq files..."
echo "Sample pattern: ${SAMPLE}_P*"

# Find all R1 and R2 files for this sample (all passages)
# Use ls instead of find for better compatibility with symlinks
R1_FILES=$(ls "${FASTQ_DIR}"/${SAMPLE}_P*_R1_001.fastq.gz 2>/dev/null | sort)
R2_FILES=$(ls "${FASTQ_DIR}"/${SAMPLE}_P*_R2_001.fastq.gz 2>/dev/null | sort)

# Check if files were found
if [ -z "$R1_FILES" ]; then
    echo "ERROR: No R1 fastq files found for sample ${SAMPLE}"
    exit 1
fi

# Count files
N_R1=$(echo "$R1_FILES" | wc -l)
N_R2=$(echo "$R2_FILES" | wc -l)

echo "Found ${N_R1} R1 files and ${N_R2} R2 files"
echo "R1 files:"
echo "$R1_FILES"
echo ""
echo "R2 files:"
echo "$R2_FILES"

# Merge R1 files
echo "Merging R1 files..."
cat $R1_FILES > "${MERGED_R1}"

# Merge R2 files
echo "Merging R2 files..."
cat $R2_FILES > "${MERGED_R2}"

echo "Merged files created:"
echo "  R1: ${MERGED_R1}"
echo "  R2: ${MERGED_R2}"

##########################################################################################################################
# Step 2: Adapter trimming with trim_galore
##########################################################################################################################

echo ""
echo "Step 2: Trimming adapters with trim_galore..."

module load cutadapt/4.1-GCCcore-11.2.0 Trim_Galore/0.6.7-GCCcore-11.2.0

trim_galore --paired --cores 4 -o "${SAMPLE_TMP}" --illumina --gzip "${MERGED_R1}" "${MERGED_R2}"

module unload cutadapt/4.1-GCCcore-11.2.0 Trim_Galore/0.6.7-GCCcore-11.2.0

# Clean up merged files
rm "${MERGED_R1}" "${MERGED_R2}"

echo "Trimming complete."
echo "  Trimmed R1: ${TRIMMED_R1}"
echo "  Trimmed R2: ${TRIMMED_R2}"

##########################################################################################################################
# Step 3: Alignment with BWA-MEM
##########################################################################################################################

echo ""
echo "Step 3: Aligning with BWA-MEM..."

module load BWA SAMtools

# Add read group information during alignment
# RGID: Read Group ID (sample name)
# RGSM: Sample name (used by GATK)
# RGPL: Platform (Illumina)
# RGLB: Library
# RGPU: Platform unit (flowcell-barcode.lane)

bwa mem -t 8 \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
    "${BWA_INDEX}" \
    "${TRIMMED_R1}" \
    "${TRIMMED_R2}" \
    > "${SAM_FILE}"

echo "Alignment complete. Converting to BAM..."

# Convert SAM to BAM
samtools view -@ 8 -bS "${SAM_FILE}" -o "${RAW_BAM}"

# Sort BAM
echo "Sorting BAM..."
samtools sort -@ 8 -o "${SORTED_BAM}" "${RAW_BAM}"

# Index sorted BAM
echo "Indexing sorted BAM..."
samtools index "${SORTED_BAM}"

# Clean up intermediate files
rm "${SAM_FILE}" "${RAW_BAM}" "${TRIMMED_R1}" "${TRIMMED_R2}"

echo "Alignment statistics:"
samtools flagstat "${SORTED_BAM}"

module unload BWA SAMtools

##########################################################################################################################
# Step 4: Mark Duplicates with Picard
##########################################################################################################################

echo ""
echo "Step 4: Marking duplicates with Picard..."

module load picard

java -Xmx32g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I="${SORTED_BAM}" \
    O="${MARKDUP_BAM}" \
    M="${MARKDUP_METRICS}" \
    ASSUME_SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT

module unload picard

# Clean up sorted BAM
rm "${SORTED_BAM}" "${SORTED_BAM}.bai"

echo "Duplicate marking complete."
echo "Metrics saved to: ${MARKDUP_METRICS}"

# Display duplicate metrics summary
echo ""
echo "Duplicate Metrics Summary:"
grep -A 2 "LIBRARY" "${MARKDUP_METRICS}"

##########################################################################################################################
# Step 5: Base Quality Score Recalibration (BQSR)
##########################################################################################################################

echo ""
echo "Step 5: Base Quality Score Recalibration..."

module load GATK

# BaseRecalibrator - generates recalibration table
echo "Running BaseRecalibrator..."

# Build known sites arguments
KNOWN_SITES_ARGS=""
if [ -f "${DBSNP}" ]; then
    KNOWN_SITES_ARGS="${KNOWN_SITES_ARGS} --known-sites ${DBSNP}"
fi
if [ -f "${KNOWN_INDELS}" ]; then
    KNOWN_SITES_ARGS="${KNOWN_SITES_ARGS} --known-sites ${KNOWN_INDELS}"
fi
if [ -f "${MILLS_INDELS}" ]; then
    KNOWN_SITES_ARGS="${KNOWN_SITES_ARGS} --known-sites ${MILLS_INDELS}"
fi

# If no known sites files exist, skip BQSR and use markdup BAM directly
if [ -z "$KNOWN_SITES_ARGS" ]; then
    echo "WARNING: No known sites files found. Skipping BQSR step."
    echo "Using mark duplicates BAM for variant calling."
    BQSR_BAM="${MARKDUP_BAM}"
else
    gatk --java-options "-Xmx32g" BaseRecalibrator \
        -R "${REFERENCE}" \
        -I "${MARKDUP_BAM}" \
        ${KNOWN_SITES_ARGS} \
        -O "${RECAL_TABLE}"

    # ApplyBQSR - applies recalibration to BAM
    echo "Running ApplyBQSR..."
    gatk --java-options "-Xmx32g" ApplyBQSR \
        -R "${REFERENCE}" \
        -I "${MARKDUP_BAM}" \
        --bqsr-recal-file "${RECAL_TABLE}" \
        -O "${BQSR_BAM}"

    echo "BQSR complete."
fi

##########################################################################################################################
# Step 6: GVCF Calling with HaplotypeCaller
##########################################################################################################################

echo ""
echo "Step 6: Calling GVCFs with GATK HaplotypeCaller..."

# Call variants in GVCF mode for joint calling
gatk --java-options "-Xmx32g" HaplotypeCaller \
    -R "${REFERENCE}" \
    -I "${BQSR_BAM}" \
    -O "${OUTPUT_GVCF}" \
    -ERC GVCF

module unload GATK

# Clean up temporary directory
echo ""
echo "Cleaning up temporary files in ${SAMPLE_TMP}..."
rm -rf "${SAMPLE_TMP}"

##########################################################################################################################
# Summary
##########################################################################################################################

echo ""
echo "============================================"
echo "Pipeline complete for sample: ${SAMPLE}"
echo "============================================"
echo ""
echo "Output files:"
echo "  BAM (BQSR): ${BQSR_BAM}"
echo "  BAM index: ${BQSR_BAM}.bai"
echo "  GVCF: ${OUTPUT_GVCF}"
echo "  GVCF index: ${OUTPUT_GVCF}.tbi"
echo "  Duplicate metrics: ${MARKDUP_METRICS}"
if [ -f "${RECAL_TABLE}" ]; then
    echo "  Recalibration table: ${RECAL_TABLE}"
fi
echo ""
echo "Done! Ready for joint calling."
