#!/bin/bash
#SBATCH -J pp_bulk_call
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 1-12:00:00

##########################################################################################################################
# Bulk truth-set preprocessing + single-sample variant calling.
#
# For benchmarking, public_WGA bulks were never marked for duplicates upstream, so the
# pipeline here is: AddOrReplaceReadGroups -> MarkDuplicates -> BaseRecalibrator -> ApplyBQSR
# -> HaplotypeCaller (single-sample, EMIT_VARIANTS_ONLY).
#
# Input TSV (with header) columns:
#   donor_id  sample_id  bam_path  reference_dir  output_dir
#
# Output: <output_dir>/<donor_id>/<sample_id>.bqsr.bam, .vcf.gz
##########################################################################################################################

set -euo pipefail

TASK_FILE="$1"
SCRIPTS_DIR="${2:-$(pwd)}"

if [ ! -f "$TASK_FILE" ]; then
    echo "ERROR: task file not found: $TASK_FILE"
    exit 1
fi

TASK_LINE=$(tail -n +2 "$TASK_FILE" | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ -z "$TASK_LINE" ]; then
    echo "ERROR: no task for array index ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

IFS=$'\t' read -r DONOR_ID SAMPLE_ID BAM_PATH REFERENCE_DIR OUTPUT_DIR <<< "$TASK_LINE"

REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
DONOR_OUT="${OUTPUT_DIR}/${DONOR_ID}"
RG_BAM="${DONOR_OUT}/${SAMPLE_ID}.rg.bam"
MARKDUP_BAM="${DONOR_OUT}/${SAMPLE_ID}.markdup.bam"
DUP_METRICS="${DONOR_OUT}/${SAMPLE_ID}.duplicate_metrics.txt"
RECAL_TABLE="${DONOR_OUT}/${SAMPLE_ID}.recal.table"
BQSR_BAM="${DONOR_OUT}/${SAMPLE_ID}.bqsr.bam"
OUTPUT_VCF="${DONOR_OUT}/${SAMPLE_ID}.vcf.gz"

mkdir -p "${DONOR_OUT}"

echo "=========================================="
echo "Bulk preprocess+call: donor=${DONOR_ID} sample=${SAMPLE_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Input BAM:     ${BAM_PATH}"
echo "Reference:     ${REFERENCE}"
echo "Output VCF:    ${OUTPUT_VCF}"
echo "=========================================="
echo ""

if [ ! -f "$BAM_PATH" ]; then
    echo "ERROR: input BAM not found: $BAM_PATH"
    exit 1
fi

if [ -f "${OUTPUT_VCF}" ] && [ -f "${OUTPUT_VCF}.tbi" ]; then
    echo "VCF already exists; skipping."
    exit 0
fi

##########################################################################################################################
# Detect known sites
##########################################################################################################################

source "${SCRIPTS_DIR}/scripts/utils/detect_known_sites.sh"
if detect_known_sites "${REFERENCE_DIR}"; then
    echo "BQSR known sites: ${KNOWN_SITES_ARGS}"
else
    echo "ERROR: no known-sites files found in ${REFERENCE_DIR}/bundle"
    exit 1
fi
echo ""

##########################################################################################################################
# Step 1: AddReadGroups
##########################################################################################################################

echo "Step 1: AddOrReplaceReadGroups..."
module load picard SAMtools

java -Xmx20g -jar "$EBROOTPICARD/picard.jar" AddOrReplaceReadGroups \
    I="${BAM_PATH}" \
    O="${RG_BAM}" \
    RGID="${SAMPLE_ID}" \
    RGSM="${SAMPLE_ID}" \
    RGPL=illumina \
    RGLB="lib_${SAMPLE_ID}" \
    RGPU="unit_${SAMPLE_ID}" \
    SORT_ORDER=coordinate

samtools index "${RG_BAM}"

##########################################################################################################################
# Step 2: MarkDuplicates
##########################################################################################################################

echo "Step 2: MarkDuplicates..."
java -Xmx20g -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
    INPUT="${RG_BAM}" \
    OUTPUT="${MARKDUP_BAM}" \
    METRICS_FILE="${DUP_METRICS}" \
    ASSUME_SORT_ORDER=coordinate \
    CREATE_INDEX=true

rm -f "${RG_BAM}" "${RG_BAM}.bai" "${RG_BAM%.bam}.bai"
module unload picard

##########################################################################################################################
# Step 3: BQSR
##########################################################################################################################

echo "Step 3: BaseRecalibrator..."
module load GATK

gatk BaseRecalibrator \
    -R "${REFERENCE}" \
    -I "${MARKDUP_BAM}" \
    ${KNOWN_SITES_ARGS} \
    -O "${RECAL_TABLE}"

echo "Step 4: ApplyBQSR..."
gatk ApplyBQSR \
    -R "${REFERENCE}" \
    -I "${MARKDUP_BAM}" \
    --bqsr-recal-file "${RECAL_TABLE}" \
    -O "${BQSR_BAM}"

# Drop intermediate markdup BAM (keep .bqsr as the analysis-ready bulk BAM)
rm -f "${MARKDUP_BAM}" "${MARKDUP_BAM%.bam}.bai" "${MARKDUP_BAM}.bai" "${RECAL_TABLE}"

##########################################################################################################################
# Step 5: HaplotypeCaller (single-sample, variants-only)
##########################################################################################################################

echo "Step 5: HaplotypeCaller (single-sample)..."
gatk --java-options "-Xmx20g" HaplotypeCaller \
    -R "${REFERENCE}" \
    -I "${BQSR_BAM}" \
    -O "${OUTPUT_VCF}"

module unload GATK SAMtools

echo ""
echo "=========================================="
echo "Done: ${SAMPLE_ID} -> ${OUTPUT_VCF}"
echo "=========================================="
