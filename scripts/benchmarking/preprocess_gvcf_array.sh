#!/bin/bash
#SBATCH -J pp_gvcf
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 12:00:00

##########################################################################################################################
# Per-cell preprocessing + GVCF generation for benchmarking samples
#
# Adapted from scripts/CapWGS/sc_extract_preprocess_qc_array.sh and sc_var_array_gatk.sh.
# Differences vs. the CapWGS scripts:
#   - The input BAM is already a per-cell BAM (e.g., 100M downsampled public_WGA), so the
#     SAMtools extraction step is removed.
#   - MarkDuplicates is assumed to have run upstream (markdup_and_metrics_array.sh on the
#     downsampled public BAMs; HSC4 BAMs were marked during sc_extract_preprocess_qc_array).
#   - GVCF output is placed in a per-donor directory so the joint-calling step can pick
#     up only that donor's cells.
#
# Pipeline per task: AddOrReplaceReadGroups -> BaseRecalibrator -> ApplyBQSR
#                    -> HaplotypeCaller (GVCF mode)
#
# Input: a TSV (with header) of per-cell tasks with columns:
#   donor_id  sample_id  bam_path  reference_dir  output_dir
# One SLURM array task per data row.
##########################################################################################################################

set -euo pipefail

TASK_FILE="$1"
SCRIPTS_DIR="${2:-$(pwd)}"

if [ ! -f "$TASK_FILE" ]; then
    echo "ERROR: task file not found: $TASK_FILE"
    exit 1
fi

# Skip header, get this array task's row
TASK_LINE=$(tail -n +2 "$TASK_FILE" | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ -z "$TASK_LINE" ]; then
    echo "ERROR: no task for array index ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

IFS=$'\t' read -r DONOR_ID SAMPLE_ID BAM_PATH REFERENCE_DIR OUTPUT_DIR <<< "$TASK_LINE"

REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
DONOR_OUT="${OUTPUT_DIR}/${DONOR_ID}"
RG_BAM="${DONOR_OUT}/${SAMPLE_ID}.rg.bam"
RECAL_TABLE="${DONOR_OUT}/${SAMPLE_ID}.recal.table"
PREPROCESSED_BAM="${DONOR_OUT}/${SAMPLE_ID}.preprocessed.bam"
OUTPUT_GVCF="${DONOR_OUT}/${SAMPLE_ID}.g.vcf.gz"

mkdir -p "${DONOR_OUT}"

echo "=========================================="
echo "Preprocess + GVCF: donor=${DONOR_ID}  sample=${SAMPLE_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Input BAM   : ${BAM_PATH}"
echo "Reference   : ${REFERENCE}"
echo "Output GVCF : ${OUTPUT_GVCF}"
echo "=========================================="
echo ""

if [ ! -f "$BAM_PATH" ]; then
    echo "ERROR: input BAM not found: $BAM_PATH"
    exit 1
fi

# Skip if final GVCF already exists (idempotent reruns)
if [ -f "${OUTPUT_GVCF}" ] && [ -f "${OUTPUT_GVCF}.tbi" ]; then
    echo "GVCF already exists; skipping."
    exit 0
fi

##########################################################################################################################
# Detect known sites for BQSR
##########################################################################################################################

source "${SCRIPTS_DIR}/scripts/utils/detect_known_sites.sh"

if detect_known_sites "${REFERENCE_DIR}"; then
    echo "BQSR known sites: ${KNOWN_SITES_ARGS}"
else
    KNOWN_SITES_ARGS=""
    echo "ERROR: no known-sites files found in ${REFERENCE_DIR}/bundle"
    echo "BQSR is required for this pipeline; aborting."
    exit 1
fi
echo ""

##########################################################################################################################
# Step 1: Add/replace read groups (HaplotypeCaller and BQSR require an @RG)
##########################################################################################################################

echo "Step 1: AddOrReplaceReadGroups..."
module load picard SAMtools

java -Xmx12g -jar "$EBROOTPICARD/picard.jar" AddOrReplaceReadGroups \
    I="${BAM_PATH}" \
    O="${RG_BAM}" \
    RGID="${SAMPLE_ID}" \
    RGSM="${SAMPLE_ID}" \
    RGPL=illumina \
    RGLB="lib_${SAMPLE_ID}" \
    RGPU="unit_${SAMPLE_ID}" \
    SORT_ORDER=coordinate

samtools index "${RG_BAM}"

module unload picard

##########################################################################################################################
# Step 2: BQSR
##########################################################################################################################

echo "Step 2: BaseRecalibrator..."
module load GATK

gatk BaseRecalibrator \
    -R "${REFERENCE}" \
    -I "${RG_BAM}" \
    ${KNOWN_SITES_ARGS} \
    -O "${RECAL_TABLE}"

echo "Step 3: ApplyBQSR..."
gatk ApplyBQSR \
    -R "${REFERENCE}" \
    -I "${RG_BAM}" \
    --bqsr-recal-file "${RECAL_TABLE}" \
    -O "${PREPROCESSED_BAM}"

# Drop the intermediate read-group BAM
rm -f "${RG_BAM}" "${RG_BAM}.bai" "${RG_BAM%.bam}.bai" "${RECAL_TABLE}"

##########################################################################################################################
# Step 4: HaplotypeCaller in GVCF mode
##########################################################################################################################

echo "Step 4: HaplotypeCaller (GVCF mode)..."

gatk --java-options "-Xmx12g" HaplotypeCaller \
    -R "${REFERENCE}" \
    -I "${PREPROCESSED_BAM}" \
    -O "${OUTPUT_GVCF}" \
    -ERC GVCF

module unload GATK SAMtools

echo ""
echo "=========================================="
echo "Done: ${SAMPLE_ID} -> ${OUTPUT_GVCF}"
echo "=========================================="
