#!/bin/bash
#SBATCH -J sc_extract_pp_qc
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 12:00:00

##########################################################################################################################
# Unified single-cell processing array job
# For each cell: extraction → MarkDuplicates → BQSR → Bigwig → Lorenz → Picard QC metrics
# This eliminates all race conditions by doing everything in a single job per cell
##########################################################################################################################

set -euo pipefail

BULK_BAM="$1"
BARCODE_FILE="$2"
SC_OUTPUT_DIR="$3"
QC_METRICS_DIR="$4"
REFERENCE_INPUT="$5"
SCRIPTS_DIR="$6"
BINSIZE="${7:-1000}"  # Bigwig bin size (default: 1000)

# Determine if input is a directory or file, and find the FASTA
if [ -f "${REFERENCE_INPUT}" ]; then
    # Direct path to FASTA file
    REFERENCE="${REFERENCE_INPUT}"
elif [ -d "${REFERENCE_INPUT}" ]; then
    # Directory - search for FASTA file in standard locations
    if [ -f "${REFERENCE_INPUT}/BWAIndex/genome.fa" ]; then
        REFERENCE="${REFERENCE_INPUT}/BWAIndex/genome.fa"
    elif [ -f "${REFERENCE_INPUT}/genome.fa" ]; then
        REFERENCE="${REFERENCE_INPUT}/genome.fa"
    elif compgen -G "${REFERENCE_INPUT}/BWAIndex/"*.fa > /dev/null 2>&1; then
        REFERENCE=$(ls "${REFERENCE_INPUT}/BWAIndex/"*.fa 2>/dev/null | head -1)
    elif compgen -G "${REFERENCE_INPUT}/BWAIndex/"*.fna > /dev/null 2>&1; then
        REFERENCE=$(ls "${REFERENCE_INPUT}/BWAIndex/"*.fna 2>/dev/null | head -1)
    elif compgen -G "${REFERENCE_INPUT}/"*.fa > /dev/null 2>&1; then
        REFERENCE=$(ls "${REFERENCE_INPUT}/"*.fa 2>/dev/null | head -1)
    elif compgen -G "${REFERENCE_INPUT}/"*.fna > /dev/null 2>&1; then
        REFERENCE=$(ls "${REFERENCE_INPUT}/"*.fna 2>/dev/null | head -1)
    else
        echo "Error: Could not find FASTA file in ${REFERENCE_INPUT}"
        exit 1
    fi
else
    echo "Error: Reference input is neither file nor directory: ${REFERENCE_INPUT}"
    exit 1
fi

# Get the cell barcode for this array task
BARCODE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BARCODE_FILE")

echo "=========================================="
echo "Processing cell: ${BARCODE}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Reference: ${REFERENCE}"
echo "=========================================="
echo ""

# Define paths
RAW_BAM="${SC_OUTPUT_DIR}/${BARCODE}.bam"
PREPROCESSED_BAM="${SC_OUTPUT_DIR}/${BARCODE}.preprocessed.bam"
BIGWIG="${SC_OUTPUT_DIR}/${BARCODE}.bw"

# QC output files (all in qc_metrics directory)
ALIGNMENT_METRICS="${QC_METRICS_DIR}/${BARCODE}_alignment_metrics.txt"
GC_METRICS="${QC_METRICS_DIR}/${BARCODE}_gc_metrics.txt"
DUPLICATE_METRICS="${QC_METRICS_DIR}/${BARCODE}_duplicate_metrics.txt"
WGS_METRICS="${QC_METRICS_DIR}/${BARCODE}_wgs_metrics.txt"
LORENZ_CSV="${QC_METRICS_DIR}/${BARCODE}_lorenz.csv"
GINI_TXT="${QC_METRICS_DIR}/${BARCODE}_gini.txt"

# Known sites for BQSR (if available)
# Get directory containing reference FASTA for checking dbsnp
REFERENCE_DIR=$(dirname "${REFERENCE}")
KNOWN_SITES=""
if [ -f "${REFERENCE_DIR}/dbsnp.vcf.gz" ]; then
    KNOWN_SITES="--known-sites ${REFERENCE_DIR}/dbsnp.vcf.gz"
fi

mkdir -p "${SC_OUTPUT_DIR}" "${QC_METRICS_DIR}"

##########################################################################################################################
# Step 1: Extract single-cell BAM from bulk BAM
##########################################################################################################################

echo "Step 1: Extracting reads for barcode ${BARCODE}..."

module load SAMtools

# Extract reads for this barcode, filtering out secondary/supplementary alignments to avoid duplicates
# -F 0x900 filters out secondary (0x100) and supplementary (0x800) alignments
samtools view -h -F 0x900 "${BULK_BAM}" | \
    grep -E "^@|CB:Z:${BARCODE}" | \
    samtools sort -@ 4 -o "${RAW_BAM}"

samtools index "${RAW_BAM}"

echo "✓ Extraction complete"
echo ""

##########################################################################################################################
# Step 2: MarkDuplicates
##########################################################################################################################

echo "Step 2: Running MarkDuplicates..."

module load picard

MARKDUP_BAM="${SC_OUTPUT_DIR}/${BARCODE}.markdup.bam"

java -Xmx12g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    INPUT="${RAW_BAM}" \
    OUTPUT="${MARKDUP_BAM}" \
    METRICS_FILE="${DUPLICATE_METRICS}" \
    ASSUME_SORT_ORDER=coordinate \
    CREATE_INDEX=true

echo "✓ MarkDuplicates complete"

# Clean up raw BAM
rm "${RAW_BAM}" "${RAW_BAM}.bai"

##########################################################################################################################
# Step 3: BQSR (if known sites available)
##########################################################################################################################

if [ -n "${KNOWN_SITES}" ]; then
    echo "Step 3: Running BQSR..."

    module load GATK

    RECAL_TABLE="${SC_OUTPUT_DIR}/${BARCODE}.recal_data.table"

    gatk BaseRecalibrator \
        -R "${REFERENCE}" \
        -I "${MARKDUP_BAM}" \
        ${KNOWN_SITES} \
        -O "${RECAL_TABLE}"

    gatk ApplyBQSR \
        -R "${REFERENCE}" \
        -I "${MARKDUP_BAM}" \
        --bqsr-recal-file "${RECAL_TABLE}" \
        -O "${PREPROCESSED_BAM}"

    rm "${MARKDUP_BAM}" "${RECAL_TABLE}"

    # Handle index file naming
    if [ -f "${MARKDUP_BAM%.bam}.bai" ]; then
        rm "${MARKDUP_BAM%.bam}.bai"
    elif [ -f "${MARKDUP_BAM}.bai" ]; then
        rm "${MARKDUP_BAM}.bai"
    fi

    echo "✓ BQSR complete"
else
    echo "Step 3: Skipping BQSR (no known sites)"

    mv "${MARKDUP_BAM}" "${PREPROCESSED_BAM}"

    if [ -f "${MARKDUP_BAM%.bam}.bai" ]; then
        mv "${MARKDUP_BAM%.bam}.bai" "${PREPROCESSED_BAM}.bai"
    elif [ -f "${MARKDUP_BAM}.bai" ]; then
        mv "${MARKDUP_BAM}.bai" "${PREPROCESSED_BAM}.bai"
    fi
fi

echo ""

##########################################################################################################################
# Step 4: Generate Bigwig
##########################################################################################################################

echo "Step 4: Generating bigwig (${BINSIZE}bp bins)..."

module load deepTools

bamCoverage -b "${PREPROCESSED_BAM}" -o "${BIGWIG}" \
    --binSize ${BINSIZE} -of bigwig -p 4

module unload deepTools

echo "✓ Bigwig complete"
echo ""

##########################################################################################################################
# Step 5: Generate Lorenz curve
##########################################################################################################################

echo "Step 5: Generating Lorenz curve and Gini coefficient..."

eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

python "${SCRIPTS_DIR}/scripts/CapWGS_QC/lorenz.py" \
    "${BIGWIG}" -o "${LORENZ_CSV}" --gini-output "${GINI_TXT}"

echo "✓ Lorenz curve and Gini coefficient complete"
echo ""

##########################################################################################################################
# Step 6: Collect QC metrics
##########################################################################################################################

echo "Step 6: Collecting benchmarking QC metrics..."

module load picard GATK R

PICARD="java -jar ${EBROOTPICARD}/picard.jar"

# Alignment summary metrics
${PICARD} CollectAlignmentSummaryMetrics \
    I="${PREPROCESSED_BAM}" \
    O="${ALIGNMENT_METRICS}" \
    R="${REFERENCE}"

# GC bias metrics
${PICARD} CollectGcBiasMetrics \
    I="${PREPROCESSED_BAM}" \
    O="${GC_METRICS}" \
    CHART="${QC_METRICS_DIR}/${BARCODE}_gc_bias.pdf" \
    S="${QC_METRICS_DIR}/${BARCODE}_gc_summary.txt" \
    R="${REFERENCE}"

# WGS metrics
${PICARD} CollectWgsMetrics \
    I="${PREPROCESSED_BAM}" \
    O="${WGS_METRICS}" \
    R="${REFERENCE}"

echo "✓ QC metrics complete"
echo ""

##########################################################################################################################
# Summary
##########################################################################################################################

echo "=========================================="
echo "✓ All processing complete for ${BARCODE}!"
echo "=========================================="
echo ""
echo "Outputs:"
echo "  Preprocessed BAM: ${PREPROCESSED_BAM}"
echo "  Bigwig: ${BIGWIG}"
echo "  Lorenz curve: ${LORENZ_CSV}"
echo "  QC metrics: ${QC_METRICS_DIR}/${BARCODE}_*.txt"
echo ""

# Verify key outputs exist
if [ ! -f "${PREPROCESSED_BAM}" ] || [ ! -f "${BIGWIG}" ] || [ ! -f "${LORENZ_CSV}" ]; then
    echo "ERROR: Some outputs missing!"
    exit 1
fi

echo "Done!"
echo ""
