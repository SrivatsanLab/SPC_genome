#!/bin/bash

###########################################################################################################################
# Submit QC analysis for any set of BAM files
# Usage: bash scripts/submit_benchmarking_qc.sh <sample_list> <output_dir> <bam_dir> [reference]
#
# Arguments:
#   sample_list: Path to file containing sample names (one per line, without .bam extension)
#   output_dir:  Directory where QC metrics will be saved
#   bam_dir:     Directory containing BAM files
#   reference:   Reference genome fasta (optional, defaults to hg38)
#
# Example:
#   bash scripts/submit_benchmarking_qc.sh bin/benchmarking/PTA_names.txt data/benchmarking/qc_metrics/PTA data/benchmarking/aligned/PTA
#   bash scripts/submit_benchmarking_qc.sh bin/benchmarking/LIANTI_names.txt data/benchmarking/qc_metrics/LIANTI data/benchmarking/aligned/LIANTI
###########################################################################################################################

set -euo pipefail

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 <sample_list> <output_dir> <bam_dir> [reference]"
    echo ""
    echo "Arguments:"
    echo "  sample_list: Path to file containing sample names (one per line)"
    echo "  output_dir:  Directory where QC metrics will be saved"
    echo "  bam_dir:     Directory containing BAM files"
    echo "  reference:   Reference genome fasta (optional, defaults to hg38)"
    exit 1
fi

SAMPLE_LIST="$1"
OUTPUT_DIR="$2"
BAM_DIR="$3"
REFERENCE="${4:-/shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/version0.5.x/genome.fa}"

# Validate input file exists
if [ ! -f "${SAMPLE_LIST}" ]; then
    echo "ERROR: Sample list file not found: ${SAMPLE_LIST}"
    exit 1
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Count samples
SAMPLE_COUNT=$(wc -l < "${SAMPLE_LIST}")

echo "=========================================="
echo "Submitting QC analysis job array"
echo "=========================================="
echo "Sample list: ${SAMPLE_LIST}"
echo "Sample count: ${SAMPLE_COUNT}"
echo "BAM directory: ${BAM_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Reference: ${REFERENCE}"
echo ""

# Submit array job
JOB_ID=$(sbatch --parsable \
    --array=1-${SAMPLE_COUNT} \
    scripts/CapWGS_QC/benchmarking_qc_array.sh \
    "${SAMPLE_LIST}" \
    "${OUTPUT_DIR}" \
    "${REFERENCE}" \
    "${BAM_DIR}")

echo "Job submitted with ID: ${JOB_ID}"
echo ""
echo "Monitor progress with:"
echo "  squeue -j ${JOB_ID}"
echo ""
echo "Results will be saved to: ${OUTPUT_DIR}"
