#!/bin/bash
#SBATCH -J bcf_joint_call
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 4:00:00

##########################################################################################################################
# Joint calling of bulk samples using BCFtools merge
# This script merges individual VCFs from BCFtools variant calling into a single multi-sample VCF
##########################################################################################################################

set -euo pipefail

# Arguments
VCF_DIR="$1"           # Directory containing individual VCF files
OUTPUT_VCF="$2"        # Output merged VCF file
COHORT_NAME="$3"       # Name for this cohort (e.g., "AAVS_PolE_clones_bcf")

echo "============================================"
echo "BCFtools joint calling for cohort: ${COHORT_NAME}"
echo "============================================"
echo ""

# Define paths
OUTPUT_DIR=$(dirname "${OUTPUT_VCF}")
VCF_LIST="${OUTPUT_DIR}/${COHORT_NAME}_vcf_list.txt"

mkdir -p "${OUTPUT_DIR}"

##########################################################################################################################
# Step 1: Create list of VCF files to merge
##########################################################################################################################

echo "Step 1: Finding VCF files to merge..."

# Find all VCF files and create list
> "${VCF_LIST}"

for vcf in "${VCF_DIR}"/*.vcf.gz; do
    if [ -f "$vcf" ]; then
        # Verify index exists
        if [ ! -f "${vcf}.csi" ] && [ ! -f "${vcf}.tbi" ]; then
            echo "WARNING: Index not found for ${vcf}, creating index..."
            module load BCFtools
            bcftools index "${vcf}"
        fi
        echo "${vcf}" >> "${VCF_LIST}"
    fi
done

# Check if any VCFs were found
n_samples=$(wc -l < "${VCF_LIST}")
if [ "$n_samples" -eq 0 ]; then
    echo "ERROR: No VCF files found in ${VCF_DIR}"
    exit 1
fi

echo "Found ${n_samples} VCF files for merging:"
cat "${VCF_LIST}"
echo ""

##########################################################################################################################
# Step 2: Merge VCFs using BCFtools merge
##########################################################################################################################

echo "Step 2: Merging VCFs with BCFtools..."

module load BCFtools

# Merge all VCFs into a single multi-sample VCF
# Options:
#   -m none: do not merge records (keep sites separate if they differ across samples)
#   -O z: output compressed VCF
#   --threads: use multiple threads for compression
bcftools merge \
    --threads 8 \
    -m none \
    -O z \
    -o "${OUTPUT_VCF}" \
    --file-list "${VCF_LIST}"

echo "VCF merge complete."
echo ""

##########################################################################################################################
# Step 3: Index merged VCF
##########################################################################################################################

echo "Step 3: Indexing merged VCF..."

bcftools index --threads 8 "${OUTPUT_VCF}"

echo "Indexing complete."
echo ""

##########################################################################################################################
# Step 4: Generate statistics
##########################################################################################################################

echo "Step 4: Generating variant statistics..."

bcftools stats "${OUTPUT_VCF}" > "${OUTPUT_VCF%.vcf.gz}_stats.txt"

echo ""
echo "============================================"
echo "BCFtools joint calling complete!"
echo "============================================"
echo ""
echo "Output files:"
echo "  Merged VCF: ${OUTPUT_VCF}"
echo "  VCF index: ${OUTPUT_VCF}.csi"
echo "  VCF list: ${VCF_LIST}"
echo "  Statistics: ${OUTPUT_VCF%.vcf.gz}_stats.txt"
echo ""

# Variant summary
echo "=== Variant Summary ==="
echo "Total samples: ${n_samples}"
bcftools view -H "${OUTPUT_VCF}" | wc -l | awk '{print "Total variant sites: " $1}'
bcftools view -H -v snps "${OUTPUT_VCF}" | wc -l | awk '{print "SNP sites: " $1}'
bcftools view -H -v indels "${OUTPUT_VCF}" | wc -l | awk '{print "Indel sites: " $1}'

# Show sample names in merged VCF
echo ""
echo "=== Sample Names in Merged VCF ==="
bcftools query -l "${OUTPUT_VCF}"

module unload BCFtools

echo ""
echo "Done!"
