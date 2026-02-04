#!/bin/bash
#SBATCH -J bench_qc
#SBATCH --output=SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 4:00:00

###########################################################################################################################
# Collect QC metrics for benchmarking BAM files using Picard and GATK
# Metrics: Alignment summary, GC bias, Duplicates, WGS metrics
###########################################################################################################################

set -euo pipefail

# Arguments
SAMPLE_LIST="$1"      # File containing sample names (one per line)
OUTPUT_DIR="$2"       # Output directory for QC metrics
REFERENCE="$3"        # Reference genome fasta
BAM_DIR="$4"          # Directory containing BAM files

# Get sample name from array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

# Input BAM file
BAM_FILE="${BAM_DIR}/${SAMPLE}.bam"

# Output files
ALIGNMENT_METRICS="${OUTPUT_DIR}/${SAMPLE}_alignment_metrics.txt"
GC_METRICS="${OUTPUT_DIR}/${SAMPLE}_gc_metrics.txt"
GC_SUMMARY="${OUTPUT_DIR}/${SAMPLE}_gc_summary.txt"
DUPLICATE_METRICS="${OUTPUT_DIR}/${SAMPLE}_duplicate_metrics.txt"
WGS_METRICS="${OUTPUT_DIR}/${SAMPLE}_wgs_metrics.txt"

echo "Processing sample: ${SAMPLE}"
echo "BAM file: ${BAM_FILE}"

# Load required modules
module load picard
module load GATK
module load R

# Set up Picard command
PICARD="java -jar ${EBROOTPICARD}/picard.jar"

# 1. Collect Alignment Summary Metrics
# Provides: Fraction of reads aligned, MAPQ scores, mismatch rates
echo "Collecting alignment summary metrics..."
${PICARD} CollectAlignmentSummaryMetrics \
    I="${BAM_FILE}" \
    O="${ALIGNMENT_METRICS}" \
    R="${REFERENCE}"

# 2. Collect GC Bias Metrics
# Provides: GC content distribution
echo "Collecting GC bias metrics..."
${PICARD} CollectGcBiasMetrics \
    I="${BAM_FILE}" \
    O="${GC_METRICS}" \
    CHART="${OUTPUT_DIR}/${SAMPLE}_gc_bias.pdf" \
    S="${GC_SUMMARY}" \
    R="${REFERENCE}"

# 3. Mark Duplicates (metrics only, don't write output BAM)
# Provides: Percent duplicates
echo "Collecting duplicate metrics..."
${PICARD} MarkDuplicates \
    I="${BAM_FILE}" \
    O=/dev/null \
    M="${DUPLICATE_METRICS}" \
    ASSUME_SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=LENIENT

# 4. Collect WGS Metrics (optional but useful for coverage statistics)
# Provides: Mean coverage, coverage distribution
echo "Collecting WGS metrics..."
${PICARD} CollectWgsMetrics \
    I="${BAM_FILE}" \
    O="${WGS_METRICS}" \
    R="${REFERENCE}" \
    VALIDATION_STRINGENCY=LENIENT

echo "QC metrics collection complete for ${SAMPLE}"
