#!/bin/bash
#SBATCH -J sc_vcall_bcf
#SBATCH -o SLURM_outs/%x_%A_%a.out
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 2:00:00

###########################################################################################################################
# Single-cell variant calling on DNA BAMs using BCFtools
# Uses same parameters as bulk_variant_calling.sh
###########################################################################################################################

set -euo pipefail

barcode_file="$1"
sc_bam_dir="$2"
reference_fa="$3"
output_dir="$4"

# Get barcode for this array task
barcode=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$barcode_file")

# Input/output files
BAM_FILE="${sc_bam_dir}/${barcode}_dna.bam"
OUTPUT_VCF="${output_dir}/${barcode}.g.vcf.gz"

# Verify inputs
if [ ! -f "${BAM_FILE}" ]; then
    echo "Error: BAM file not found: ${BAM_FILE}"
    exit 1
fi

if [ ! -f "${reference_fa}" ]; then
    echo "Error: Reference FASTA not found: ${reference_fa}"
    exit 1
fi

# Load BCFtools and SAMtools
module load BCFtools SAMtools

echo "Starting single-cell variant calling for: ${barcode}"
echo "BAM file: ${BAM_FILE}"
echo "Reference: ${reference_fa}"
echo "Output VCF: ${OUTPUT_VCF}"
echo ""

# Check BAM file
echo "Checking BAM file..."
samtools flagstat "${BAM_FILE}"
echo ""

# Step 1: Generate mpileup and call variants
echo "Running bcftools mpileup and call..."

TMP_VCF="${OUTPUT_VCF%.g.vcf.gz}_raw.g.vcf.gz"

bcftools mpileup \
    --threads 2 \
    -f "${reference_fa}" \
    -a FORMAT/AD,FORMAT/DP \
    -O u \
    "${BAM_FILE}" | \
bcftools call \
    --threads 2 \
    -m \
    -A \
    --gvcf 1 \
    -O z \
    -o "${TMP_VCF}"

# Step 2: Normalize variants - split multiallelic sites and left-align indels
echo "Normalizing variants (splitting multiallelic sites, left-aligning indels)..."

bcftools norm \
    --threads 2 \
    -m -both \
    -f "${reference_fa}" \
    -O z \
    -o "${OUTPUT_VCF}" \
    "${TMP_VCF}"

# Clean up temporary file
rm "${TMP_VCF}"

echo "Indexing normalized VCF..."
bcftools index "${OUTPUT_VCF}"

# Quick summary
echo ""
echo "=== Variant Summary for ${barcode} ==="
bcftools view -H "${OUTPUT_VCF}" | wc -l | awk '{print "Total variants: " $1}'
bcftools view -H -v snps "${OUTPUT_VCF}" | wc -l | awk '{print "SNPs: " $1}'
bcftools view -H -v indels "${OUTPUT_VCF}" | wc -l | awk '{print "Indels: " $1}'

echo ""
echo "Variant calling complete for ${barcode}!"

module unload BCFtools SAMtools
