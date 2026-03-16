#!/bin/bash
#SBATCH -J jc_chr
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 2-00:00:00

##########################################################################################################################
# Per-chromosome joint genotyping using existing GenomicsDB for CapWGS pipeline
# This script is run as a SLURM array job, with one job per chromosome
# Much faster than processing all chromosomes sequentially
##########################################################################################################################

set -euo pipefail

# Arguments
GENOMICSDB_DIR="$1"     # Base GenomicsDB directory (contains per-chromosome databases)
OUTPUT_DIR="$2"         # Output directory for per-chromosome VCFs
REFERENCE_DIR="$3"      # Reference directory (e.g., /shared/biodata/reference/GATK/hg38)
CHROMOSOME_FILE="$4"    # File with one chromosome per line

# Get the chromosome for this array task
CHROMOSOME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHROMOSOME_FILE")

echo "============================================"
echo "CapWGS joint calling chromosome: ${CHROMOSOME}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "============================================"
echo ""

# Define paths
REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
GENOMICS_DB="${GENOMICSDB_DIR}/${CHROMOSOME}"  # Per-chromosome GenomicsDB
OUTPUT_VCF="${OUTPUT_DIR}/${CHROMOSOME}.vcf.gz"
TMP_DIR="${OUTPUT_DIR}/tmp_${CHROMOSOME}"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TMP_DIR}"

# Verify GenomicsDB exists for this chromosome
if [ ! -f "${GENOMICS_DB}/callset.json" ]; then
    echo "ERROR: GenomicsDB not found or incomplete for ${CHROMOSOME} at ${GENOMICS_DB}"
    exit 1
fi

##########################################################################################################################
# Joint genotyping with GenotypeGVCFs for this chromosome
##########################################################################################################################

echo "Running GenotypeGVCFs for ${CHROMOSOME}..."
echo "Started at: $(date)"
echo ""

module load GATK

# Perform joint genotyping for this chromosome only
gatk --java-options "-Xmx28g" GenotypeGVCFs \
    -R "${REFERENCE}" \
    -V gendb://"${GENOMICS_DB}" \
    -L "${CHROMOSOME}" \
    -O "${OUTPUT_VCF}" \
    --tmp-dir "${TMP_DIR}"

echo ""
echo "GenotypeGVCFs completed for ${CHROMOSOME} at: $(date)"
echo ""

# Clean up temp directory
rm -rf "${TMP_DIR}"

##########################################################################################################################
# Generate statistics
##########################################################################################################################

echo "Generating variant statistics for ${CHROMOSOME}..."

module load BCFtools/1.18-GCC-12.2.0

bcftools stats "${OUTPUT_VCF}" > "${OUTPUT_VCF%.vcf.gz}_stats.txt"

echo ""
echo "============================================"
echo "Joint calling complete for ${CHROMOSOME}!"
echo "============================================"
echo ""
echo "Output files:"
echo "  VCF: ${OUTPUT_VCF}"
echo "  VCF index: ${OUTPUT_VCF}.tbi"
echo "  Statistics: ${OUTPUT_VCF%.vcf.gz}_stats.txt"
echo ""

# Variant summary
echo "=== Variant Summary for ${CHROMOSOME} ==="
bcftools view -H "${OUTPUT_VCF}" | wc -l | awk '{print "Total variant sites: " $1}'
bcftools view -H -v snps "${OUTPUT_VCF}" | wc -l | awk '{print "SNP sites: " $1}'
bcftools view -H -v indels "${OUTPUT_VCF}" | wc -l | awk '{print "Indel sites: " $1}'

module unload BCFtools GATK

echo ""
echo "Done!"
