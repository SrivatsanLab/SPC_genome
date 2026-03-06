#!/bin/bash
#SBATCH -J merge_jc_vcfs
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00

##########################################################################################################################
# Merge per-chromosome VCFs into a single joint-called VCF for CapWGS pipeline
# This runs after all per-chromosome joint calling jobs complete
##########################################################################################################################

set -euo pipefail

# Arguments
PER_CHR_DIR="$1"        # Directory containing per-chromosome VCFs
OUTPUT_VCF="$2"         # Final merged VCF output
CHROMOSOME_FILE="$3"    # File with one chromosome per line

echo "============================================"
echo "Merging CapWGS per-chromosome VCFs"
echo "============================================"
echo ""
echo "Input directory: ${PER_CHR_DIR}"
echo "Output VCF: ${OUTPUT_VCF}"
echo ""

# Check that all per-chromosome VCFs exist
echo "Checking for per-chromosome VCFs..."
MISSING=0
while read -r chr; do
    if [ ! -f "${PER_CHR_DIR}/${chr}.vcf.gz" ]; then
        echo "ERROR: Missing VCF for ${chr}"
        MISSING=1
    fi
done < "${CHROMOSOME_FILE}"

if [ ${MISSING} -eq 1 ]; then
    echo "ERROR: Not all per-chromosome VCFs are present. Cannot merge."
    exit 1
fi

echo "All per-chromosome VCFs found."
echo ""

##########################################################################################################################
# Merge VCFs with bcftools concat
##########################################################################################################################

echo "Merging VCFs with bcftools concat..."
echo "Started at: $(date)"
echo ""

module load BCFtools

# Create list of VCF files in chromosome order
VCF_LIST="${PER_CHR_DIR}/vcf_list.txt"
> "${VCF_LIST}"

while read -r chr; do
    echo "${PER_CHR_DIR}/${chr}.vcf.gz" >> "${VCF_LIST}"
done < "${CHROMOSOME_FILE}"

# Concatenate VCFs (they're already sorted by chromosome)
bcftools concat \
    --file-list "${VCF_LIST}" \
    --output-type z \
    --output "${OUTPUT_VCF}" \
    --threads 15

# Index the merged VCF
echo "Indexing merged VCF..."
bcftools index --tbi "${OUTPUT_VCF}"

echo ""
echo "Merge complete at: $(date)"
echo ""

##########################################################################################################################
# Generate statistics
##########################################################################################################################

echo "Generating variant statistics..."

bcftools stats "${OUTPUT_VCF}" > "${OUTPUT_VCF%.vcf.gz}_stats.txt"

echo ""
echo "============================================"
echo "Merge complete!"
echo "============================================"
echo ""
echo "Output files:"
echo "  Merged VCF: ${OUTPUT_VCF}"
echo "  VCF index: ${OUTPUT_VCF}.tbi"
echo "  Statistics: ${OUTPUT_VCF%.vcf.gz}_stats.txt"
echo ""

# Variant summary
echo "=== Variant Summary ==="
n_samples=$(bcftools query -l "${OUTPUT_VCF}" | wc -l)
echo "Total cells: ${n_samples}"
bcftools view -H "${OUTPUT_VCF}" | wc -l | awk '{print "Total variant sites: " $1}'
bcftools view -H -v snps "${OUTPUT_VCF}" | wc -l | awk '{print "SNP sites: " $1}'
bcftools view -H -v indels "${OUTPUT_VCF}" | wc -l | awk '{print "Indel sites: " $1}'

module unload BCFtools

echo ""
echo "Done!"
echo ""
echo "Per-chromosome VCFs can be deleted if desired:"
echo "  rm -rf ${PER_CHR_DIR}"
