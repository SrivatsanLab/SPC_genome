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
REFERENCE_DIR="$5"
SCRIPTS_DIR="$6"
BINSIZE="${7:-1000}"  # Bigwig bin size (default: 1000)

# Get the cell barcode for this array task
BARCODE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BARCODE_FILE")

echo "=========================================="
echo "Processing cell: ${BARCODE}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "=========================================="
echo ""

# Define paths
RAW_BAM="${SC_OUTPUT_DIR}/${BARCODE}.bam"
PREPROCESSED_BAM="${SC_OUTPUT_DIR}/${BARCODE}.preprocessed.bam"
BIGWIG="${SC_OUTPUT_DIR}/${BARCODE}.bw"
LORENZ_CSV="${SC_OUTPUT_DIR}/${BARCODE}_lorenz.csv"
REFERENCE="${REFERENCE_DIR}/genome.fa"

# QC output files
ALIGNMENT_METRICS="${QC_METRICS_DIR}/${BARCODE}_alignment_metrics.txt"
GC_METRICS="${QC_METRICS_DIR}/${BARCODE}_gc_metrics.txt"
DUPLICATE_METRICS="${QC_METRICS_DIR}/${BARCODE}_duplicate_metrics.txt"
WGS_METRICS="${QC_METRICS_DIR}/${BARCODE}_wgs_metrics.txt"

# Known sites for BQSR (if available)
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

samtools view -h "${BULK_BAM}" | \
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

echo "Step 5: Generating Lorenz curve..."

eval "$(micromamba shell hook --shell bash)"
micromamba activate spc_genome

python "${SCRIPTS_DIR}/scripts/CapWGS_QC/lorenz_curve.py" \
    "${PREPROCESSED_BAM}" "${LORENZ_CSV}"

echo "✓ Lorenz curve complete"
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
