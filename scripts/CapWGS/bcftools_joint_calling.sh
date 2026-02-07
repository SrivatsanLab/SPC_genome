#!/bin/bash
#SBATCH -J bcf_joint
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 24:00:00

##########################################################################################################################
# BCFtools joint calling for single-cell data
# This script:
# 1. Merges all single-cell VCFs into a multi-sample VCF
# 2. Normalizes the merged VCF to split multiallelic sites and left-align indels
# 3. Filters out reference-only sites
##########################################################################################################################

set -euo pipefail

# Arguments
VCF_DIR="$1"           # Directory containing single-cell VCF files
REFERENCE_DIR="$2"     # Reference genome directory
OUTPUT_VCF="$3"        # Output merged and normalized VCF file
SAMPLE_NAME="$4"       # Sample name (e.g., "HSC_pilot")

echo "============================================"
echo "BCFtools joint calling for: ${SAMPLE_NAME}"
echo "============================================"
echo ""

# Construct reference path
REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"

# Check if reference exists, try alternative naming
if [ ! -f "${REFERENCE}" ]; then
    # Try to find any .fasta or .fa file in the reference directory
    REFERENCE=$(find "${REFERENCE_DIR}" -maxdepth 1 -name "*.fasta" -o -name "*.fa" | head -n 1)
    if [ -z "${REFERENCE}" ] || [ ! -f "${REFERENCE}" ]; then
        echo "ERROR: Reference genome not found in ${REFERENCE_DIR}"
        exit 1
    fi
fi

echo "Reference: ${REFERENCE}"
echo "VCF directory: ${VCF_DIR}"
echo "Output VCF: ${OUTPUT_VCF}"
echo ""

# Define intermediate files
OUTPUT_DIR=$(dirname "${OUTPUT_VCF}")
VCF_LIST="${OUTPUT_DIR}/${SAMPLE_NAME}_vcf_list.txt"
RAW_MERGED="${OUTPUT_VCF%.vcf.gz}_raw_merged.vcf.gz"

mkdir -p "${OUTPUT_DIR}"

##########################################################################################################################
# Step 1: Create list of VCF files to merge
##########################################################################################################################

echo "Step 1: Finding VCF files to merge..."

# Find all VCF files (excluding GVCF files if present)
find "${VCF_DIR}" -name "*.vcf.gz" ! -name "*.g.vcf.gz" > "${VCF_LIST}"

# Check if any VCFs were found
n_cells=$(wc -l < "${VCF_LIST}")
if [ "$n_cells" -eq 0 ]; then
    echo "ERROR: No VCF files found in ${VCF_DIR}"
    echo "Note: Looking for *.vcf.gz files (excluding *.g.vcf.gz files)"
    exit 1
fi

echo "Found ${n_cells} VCF files for merging"
echo ""

# Ensure all VCFs are indexed
echo "Checking VCF indices..."
module load BCFtools

while IFS= read -r vcf; do
    if [ ! -f "${vcf}.csi" ] && [ ! -f "${vcf}.tbi" ]; then
        echo "Creating index for $(basename ${vcf})..."
        bcftools index "${vcf}" || echo "Warning: Failed to index ${vcf}"
    fi
done < "${VCF_LIST}"

echo "Index check complete."
echo ""

##########################################################################################################################
# Step 2: Merge VCFs using BCFtools merge
##########################################################################################################################

if [ -f "${RAW_MERGED}" ]; then
    echo "Step 2: Raw merged VCF already exists, skipping merge step..."
    echo "Using existing: ${RAW_MERGED}"
    echo ""
else
    echo "Step 2: Merging ${n_cells} VCFs with BCFtools..."

    # Merge all VCFs into a single multi-sample VCF
    # Options:
    #   -m none: do not merge records (keep sites separate if they differ across samples)
    #   -O z: output compressed VCF
    #   --threads: use multiple threads for compression
    bcftools merge \
        --threads 15 \
        -m none \
        -O z \
        -o "${RAW_MERGED}" \
        --file-list "${VCF_LIST}"

    echo "VCF merge complete."
    echo ""
fi

##########################################################################################################################
# Step 3: Normalize merged VCF (split multiallelic sites, left-align indels) and filter
##########################################################################################################################

echo "Step 3: Normalizing merged VCF and filtering out reference-only sites..."

# Normalize and filter out sites with ALT = "."
# Pipe: normalize -> filter out ALT="." -> compress
bcftools norm \
    --threads 8 \
    -m -both \
    -f "${REFERENCE}" \
    -O u \
    "${RAW_MERGED}" | \
bcftools view \
    --threads 7 \
    -e 'ALT="."' \
    -O z \
    -o "${OUTPUT_VCF}"

echo "Normalization and filtering complete."
echo ""

##########################################################################################################################
# Step 4: Index final VCF
##########################################################################################################################

echo "Step 4: Indexing final VCF..."

bcftools index --threads 15 "${OUTPUT_VCF}"

echo "Indexing complete."
echo ""

##########################################################################################################################
# Step 5: Generate statistics
##########################################################################################################################

echo "Step 5: Generating variant statistics..."

bcftools stats "${OUTPUT_VCF}" > "${OUTPUT_VCF%.vcf.gz}_stats.txt"

echo ""
echo "============================================"
echo "BCFtools joint calling complete!"
echo "============================================"
echo ""
echo "Output files:"
echo "  Merged & normalized VCF: ${OUTPUT_VCF}"
echo "  VCF index: ${OUTPUT_VCF}.csi"
echo "  VCF list: ${VCF_LIST}"
echo "  Statistics: ${OUTPUT_VCF%.vcf.gz}_stats.txt"
echo ""

# Variant summary
echo "=== Variant Summary ==="
echo "Total cells: ${n_cells}"
bcftools view -H "${OUTPUT_VCF}" | wc -l | awk '{print "Total variant sites: " $1}'
bcftools view -H -v snps "${OUTPUT_VCF}" | wc -l | awk '{print "SNP sites: " $1}'
bcftools view -H -v indels "${OUTPUT_VCF}" | wc -l | awk '{print "Indel sites: " $1}'

# Show first few cell names in merged VCF
echo ""
echo "=== Sample Names in Merged VCF (first 10) ==="
bcftools query -l "${OUTPUT_VCF}" | head -n 10
echo "..."
echo "(Total: ${n_cells} cells)"

module unload BCFtools

echo ""
echo "Done!"
