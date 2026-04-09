#!/bin/bash
#SBATCH -J sc_extract_pp_qc_gta
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 12:00:00

##########################################################################################################################
# Unified single-cell processing array job for CapGTA
# For each cell: extract DNA & RNA BAMs → MarkDuplicates (DNA, QC only) → Bigwig → Lorenz → Picard QC → Exonic enrichment
# This eliminates all race conditions by doing everything in a single job per cell
##########################################################################################################################

set -euo pipefail

BULK_DNA_BAM="$1"
BULK_RNA_BAM="$2"
BARCODE_FILE="$3"
SC_OUTPUT_DIR="$4"
QC_METRICS_DIR="$5"
REFERENCE_INPUT="$6"
GTF_FILE="$7"
SCRIPTS_DIR="$8"
BINSIZE="${9:-1000}"          # Bigwig bin size (default: 1000)

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
echo "GTF: ${GTF_FILE}"
echo "=========================================="
echo ""

# Define paths
DNA_BAM="${SC_OUTPUT_DIR}/${BARCODE}_dna.bam"
RNA_BAM="${SC_OUTPUT_DIR}/${BARCODE}_rna.bam"
DNA_BIGWIG="${SC_OUTPUT_DIR}/${BARCODE}_dna.bw"

# QC output files (all in qc_metrics directory)
ALIGNMENT_METRICS="${QC_METRICS_DIR}/${BARCODE}_alignment_metrics.txt"
GC_METRICS="${QC_METRICS_DIR}/${BARCODE}_gc_metrics.txt"
DUPLICATE_METRICS="${QC_METRICS_DIR}/${BARCODE}_duplication_metrics.txt"
WGS_METRICS="${QC_METRICS_DIR}/${BARCODE}_wgs_metrics.txt"
LORENZ_CSV="${QC_METRICS_DIR}/${BARCODE}_lorenz.csv"
GINI_TXT="${QC_METRICS_DIR}/${BARCODE}_gini.txt"
EXONIC_ENRICHMENT_DNA="${QC_METRICS_DIR}/${BARCODE}_exonic_enrichment_dna.txt"
EXONIC_ENRICHMENT_RNA="${QC_METRICS_DIR}/${BARCODE}_exonic_enrichment_rna.txt"

mkdir -p "${SC_OUTPUT_DIR}" "${QC_METRICS_DIR}"

##########################################################################################################################
# Step 1: Extract single-cell DNA and RNA BAMs from bulk BAMs
##########################################################################################################################

echo "Step 1: Extracting DNA and RNA reads for barcode ${BARCODE}..."

module load SAMtools

# Extract DNA reads for this barcode, filtering out secondary/supplementary alignments to avoid duplicates
# -F 0x900 filters out secondary (0x100) and supplementary (0x800) alignments
samtools view -h -F 0x900 "${BULK_DNA_BAM}" | \
    grep -E "^@|CB:Z:${BARCODE}" | \
    samtools sort -@ 4 -o "${DNA_BAM}"

samtools index "${DNA_BAM}"

# Extract RNA reads for this barcode
samtools view -h -F 0x900 "${BULK_RNA_BAM}" | \
    grep -E "^@|CB:Z:${BARCODE}" | \
    samtools sort -@ 4 -o "${RNA_BAM}"

samtools index "${RNA_BAM}"

echo "✓ Extraction complete"
echo "  DNA BAM: ${DNA_BAM}"
echo "  RNA BAM: ${RNA_BAM}"
echo ""

##########################################################################################################################
# Step 2: MarkDuplicates on DNA BAM (QC only, not removed)
##########################################################################################################################

echo "Step 2: Running MarkDuplicates on DNA BAM (QC only)..."

module load picard

MARKDUP_BAM="${SC_OUTPUT_DIR}/${BARCODE}_dna.markdup.bam"

java -Xmx12g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    INPUT="${DNA_BAM}" \
    OUTPUT="${MARKDUP_BAM}" \
    METRICS_FILE="${DUPLICATE_METRICS}" \
    ASSUME_SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    REMOVE_DUPLICATES=false

echo "✓ MarkDuplicates complete"
echo ""

# Note: Keep original DNA BAM for variant calling, use markdup BAM for QC metrics only
# Clean up markdup BAM after QC metrics are collected

##########################################################################################################################
# Step 3: Generate DNA Bigwig
##########################################################################################################################

echo "Step 3: Generating DNA bigwig (${BINSIZE}bp bins)..."

module load deepTools

bamCoverage -b "${DNA_BAM}" -o "${DNA_BIGWIG}" \
    --binSize ${BINSIZE} -of bigwig -p 4

module unload deepTools

echo "✓ DNA Bigwig complete"
echo ""

##########################################################################################################################
# Step 4: Generate Lorenz curve and Gini coefficient from DNA coverage
##########################################################################################################################

echo "Step 4: Generating Lorenz curve and Gini coefficient..."

eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

python "${SCRIPTS_DIR}/scripts/utils/lorenz.py" \
    "${DNA_BIGWIG}" -o "${LORENZ_CSV}" --gini-output "${GINI_TXT}"

echo "✓ Lorenz curve and Gini coefficient complete"
echo ""

##########################################################################################################################
# Step 5: Collect Picard QC metrics from DNA BAM
##########################################################################################################################

echo "Step 5: Collecting Picard QC metrics from DNA BAM..."

module load picard GATK R

PICARD="java -jar ${EBROOTPICARD}/picard.jar"

# Alignment summary metrics (use markdup BAM for accurate duplication stats)
${PICARD} CollectAlignmentSummaryMetrics \
    I="${MARKDUP_BAM}" \
    O="${ALIGNMENT_METRICS}" \
    R="${REFERENCE}"

# GC bias metrics
${PICARD} CollectGcBiasMetrics \
    I="${MARKDUP_BAM}" \
    O="${GC_METRICS}" \
    CHART="${QC_METRICS_DIR}/${BARCODE}_gc_bias.pdf" \
    S="${QC_METRICS_DIR}/${BARCODE}_gc_summary.txt" \
    R="${REFERENCE}"

# WGS metrics
${PICARD} CollectWgsMetrics \
    I="${MARKDUP_BAM}" \
    O="${WGS_METRICS}" \
    R="${REFERENCE}"

echo "✓ QC metrics complete"
echo ""

# Clean up markdup BAM after QC metrics collected
rm "${MARKDUP_BAM}"
if [ -f "${MARKDUP_BAM%.bam}.bai" ]; then
    rm "${MARKDUP_BAM%.bam}.bai"
elif [ -f "${MARKDUP_BAM}.bai" ]; then
    rm "${MARKDUP_BAM}.bai"
fi

##########################################################################################################################
# Step 6: Calculate exonic enrichment for DNA and RNA BAMs
##########################################################################################################################

echo "Step 6: Calculating exonic enrichment..."

# Calculate DNA exonic enrichment
python "${SCRIPTS_DIR}/scripts/CapGTA/calculate_exonic_enrichment.py" \
    "${DNA_BAM}" \
    "${GTF_FILE}" \
    "${REFERENCE}" \
    > "${EXONIC_ENRICHMENT_DNA}"

# Calculate RNA exonic enrichment
python "${SCRIPTS_DIR}/scripts/CapGTA/calculate_exonic_enrichment.py" \
    "${RNA_BAM}" \
    "${GTF_FILE}" \
    "${REFERENCE}" \
    > "${EXONIC_ENRICHMENT_RNA}"

echo "✓ Exonic enrichment calculation complete"
echo "  DNA enrichment: ${EXONIC_ENRICHMENT_DNA}"
echo "  RNA enrichment: ${EXONIC_ENRICHMENT_RNA}"
echo ""

##########################################################################################################################
# Summary
##########################################################################################################################

echo "=========================================="
echo "✓ All processing complete for ${BARCODE}!"
echo "=========================================="
echo ""
echo "Outputs:"
echo "  DNA BAM: ${DNA_BAM}"
echo "  RNA BAM: ${RNA_BAM}"
echo "  DNA Bigwig: ${DNA_BIGWIG}"
echo "  Lorenz curve: ${LORENZ_CSV}"
echo "  Gini coefficient: ${GINI_TXT}"
echo "  QC metrics: ${QC_METRICS_DIR}/${BARCODE}_*.txt"
echo "  Exonic enrichment: ${QC_METRICS_DIR}/${BARCODE}_exonic_enrichment_*.txt"
echo ""

# Verify key outputs exist
if [ ! -f "${DNA_BAM}" ] || [ ! -f "${RNA_BAM}" ] || [ ! -f "${DNA_BIGWIG}" ] || \
   [ ! -f "${LORENZ_CSV}" ] || [ ! -f "${EXONIC_ENRICHMENT_DNA}" ] || [ ! -f "${EXONIC_ENRICHMENT_RNA}" ]; then
    echo "ERROR: Some outputs missing!"
    exit 1
fi

echo "Done!"
echo ""
