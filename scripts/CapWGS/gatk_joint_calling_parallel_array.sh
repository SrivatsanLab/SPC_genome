#!/bin/bash
#SBATCH -J jc_interval
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 2-00:00:00

##########################################################################################################################
# Per-interval joint genotyping using existing GenomicsDB for CapWGS pipeline
# This script is run as a SLURM array job, with one job per genomic interval
# Uses balanced intervals (~100 intervals) for optimal parallelization
##########################################################################################################################

set -euo pipefail

# Arguments
GENOMICSDB_DIR="$1"     # Base GenomicsDB directory (contains per-interval databases)
OUTPUT_DIR="$2"         # Output directory for per-interval VCFs
REFERENCE_DIR="$3"      # Reference directory (e.g., /shared/biodata/reference/GATK/hg38)
INTERVAL_LIST="$4"      # File with one interval file path per line

# Get the interval file for this array task
INTERVAL_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INTERVAL_LIST")

# Extract interval name from filename (e.g., "0042-scattered" from "0042-scattered.interval_list")
INTERVAL_NAME=$(basename "$INTERVAL_FILE" .interval_list)

echo "============================================"
echo "CapWGS joint calling interval: ${INTERVAL_NAME}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "============================================"
echo ""

# Define paths
REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
GENOMICS_DB="${GENOMICSDB_DIR}/${INTERVAL_NAME}"  # Per-interval GenomicsDB
OUTPUT_VCF="${OUTPUT_DIR}/${INTERVAL_NAME}.vcf.gz"
TMP_DIR="${OUTPUT_DIR}/tmp_${INTERVAL_NAME}"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TMP_DIR}"

# Verify interval file exists
if [ ! -f "${INTERVAL_FILE}" ]; then
    echo "ERROR: Interval file not found at ${INTERVAL_FILE}"
    exit 1
fi

# Verify GenomicsDB exists for this interval
if [ ! -f "${GENOMICS_DB}/callset.json" ]; then
    echo "ERROR: GenomicsDB not found or incomplete for ${INTERVAL_NAME} at ${GENOMICS_DB}"
    exit 1
fi

##########################################################################################################################
# Joint genotyping with GenotypeGVCFs for this interval
##########################################################################################################################

echo "Running GenotypeGVCFs for ${INTERVAL_NAME}..."
echo "Interval file: ${INTERVAL_FILE}"
echo "Started at: $(date)"
echo ""

module load GATK

# Perform joint genotyping for this genomic interval only
gatk --java-options "-Xmx28g" GenotypeGVCFs \
    -R "${REFERENCE}" \
    -V gendb://"${GENOMICS_DB}" \
    -L "${INTERVAL_FILE}" \
    -O "${OUTPUT_VCF}" \
    --tmp-dir "${TMP_DIR}"

echo ""
echo "GenotypeGVCFs completed for ${INTERVAL_NAME} at: $(date)"
echo ""

# Clean up temp directory
rm -rf "${TMP_DIR}"

##########################################################################################################################
# Generate statistics
##########################################################################################################################

echo "Generating variant statistics for ${INTERVAL_NAME}..."

module load BCFtools/1.18-GCC-12.2.0

bcftools stats "${OUTPUT_VCF}" > "${OUTPUT_VCF%.vcf.gz}_stats.txt"

echo ""
echo "============================================"
echo "Joint calling complete for ${INTERVAL_NAME}!"
echo "============================================"
echo ""
echo "Output files:"
echo "  VCF: ${OUTPUT_VCF}"
echo "  VCF index: ${OUTPUT_VCF}.tbi"
echo "  Statistics: ${OUTPUT_VCF%.vcf.gz}_stats.txt"
echo ""

# Variant summary
echo "=== Variant Summary for ${INTERVAL_NAME} ==="
bcftools view -H "${OUTPUT_VCF}" | wc -l | awk '{print "Total variant sites: " $1}'
bcftools view -H -v snps "${OUTPUT_VCF}" | wc -l | awk '{print "SNP sites: " $1}'
bcftools view -H -v indels "${OUTPUT_VCF}" | wc -l | awk '{print "Indel sites: " $1}'

module unload BCFtools GATK

echo ""
echo "Done!"
