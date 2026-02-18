#!/bin/bash
#SBATCH -J bulk_vcall
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 4:00:00

###########################################################################################################################
# Bulk variant calling on DNA BAM using BCFtools
# Identifies germline variants from bulk DNA sequencing
# NOTE: This script does NOT normalize variants - normalization should be done on the merged VCF
###########################################################################################################################

set -euo pipefail

BAM_FILE="$1"
REFERENCE_FA="$2"
OUTPUT_VCF="$3"

# Verify inputs
if [ ! -f "${BAM_FILE}" ]; then
    echo "Error: BAM file not found: ${BAM_FILE}"
    exit 1
fi

if [ ! -f "${REFERENCE_FA}" ]; then
    echo "Error: Reference FASTA not found: ${REFERENCE_FA}"
    exit 1
fi

# Load BCFtools and SAMtools
module load BCFtools SAMtools

echo "Starting bulk variant calling..."
echo "BAM file: ${BAM_FILE}"
echo "Reference: ${REFERENCE_FA}"
echo "Output VCF: ${OUTPUT_VCF}"
echo ""

# Check BAM file
echo "Checking BAM file..."
samtools flagstat "${BAM_FILE}"
echo ""

# Step 1: Generate mpileup and call variants
echo "Running bcftools mpileup and call..."
echo "Note: Normalization will be performed on the merged multi-sample VCF"

bcftools mpileup \
    --threads 8 \
    -f "${REFERENCE_FA}" \
    -a FORMAT/AD,FORMAT/DP \
    -O u \
    "${BAM_FILE}" | \
bcftools call \
    --threads 8 \
    -m \
    --gvcf 1 \
    -O z \
    -o "${OUTPUT_VCF}"

# Step 2: Index VCF
echo "Indexing VCF..."
bcftools index "${OUTPUT_VCF}"

# Generate variant statistics
echo "Generating variant statistics..."
bcftools stats "${OUTPUT_VCF}" > "${OUTPUT_VCF%.vcf.gz}_stats.txt"

echo ""
echo "Variant calling complete!"
echo "Output VCF: ${OUTPUT_VCF}"
echo "Stats: ${OUTPUT_VCF%.vcf.gz}_stats.txt"

# Quick summary
echo ""
echo "=== Variant Summary ==="
bcftools view -H "${OUTPUT_VCF}" | wc -l | awk '{print "Total variants: " $1}'
bcftools view -H -v snps "${OUTPUT_VCF}" | wc -l | awk '{print "SNPs: " $1}'
bcftools view -H -v indels "${OUTPUT_VCF}" | wc -l | awk '{print "Indels: " $1}'

module unload BCFtools SAMtools
