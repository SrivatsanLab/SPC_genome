#!/bin/bash
#SBATCH -J bcf_merge_norm
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 32
#SBATCH -t 24:00:00

##########################################################################################################################
# Merge individual VCFs and normalize the merged result
# This script:
# 1. Merges all individual VCFs (mix of raw and normalized) into multi-sample VCF
# 2. Normalizes the merged VCF to split multiallelic sites and left-align indels
##########################################################################################################################

set -euo pipefail

# Arguments
VCF_DIR="$1"           # Directory containing individual VCF files
REFERENCE="$2"         # Reference FASTA file
OUTPUT_VCF="$3"        # Output merged and normalized VCF file
COHORT_NAME="$4"       # Name for this cohort (e.g., "AAVS_PolE_bcf")

echo "============================================"
echo "BCFtools merge and normalize for cohort: ${COHORT_NAME}"
echo "============================================"
echo ""

# Define paths
OUTPUT_DIR=$(dirname "${OUTPUT_VCF}")
VCF_LIST="${OUTPUT_DIR}/${COHORT_NAME}_vcf_list.txt"
RAW_MERGED="${OUTPUT_VCF%.vcf.gz}_raw_merged.vcf.gz"

mkdir -p "${OUTPUT_DIR}"

##########################################################################################################################
# Step 1: Create list of VCF files to merge (prefer raw VCFs)
##########################################################################################################################

echo "Step 1: Finding VCF files to merge..."

# Create list, preferring _raw.vcf.gz files where available
> "${VCF_LIST}"

for vcf in "${VCF_DIR}"/*.vcf.gz; do
    if [ -f "$vcf" ]; then
        basename=$(basename "$vcf" .vcf.gz)

        # Skip if this is a raw file (we'll add it below)
        if [[ "$basename" == *"_raw" ]]; then
            continue
        fi

        # Check if raw version exists
        raw_vcf="${VCF_DIR}/${basename}_raw.vcf.gz"
        if [ -f "$raw_vcf" ]; then
            # Use raw version
            vcf_to_use="$raw_vcf"
        else
            # Use normalized version
            vcf_to_use="$vcf"
        fi

        # Verify index exists
        if [ ! -f "${vcf_to_use}.csi" ] && [ ! -f "${vcf_to_use}.tbi" ]; then
            echo "Creating index for $(basename ${vcf_to_use})..."
            module load BCFtools
            bcftools index "${vcf_to_use}" || echo "Warning: Failed to index ${vcf_to_use}, will try without index"
        fi

        echo "${vcf_to_use}" >> "${VCF_LIST}"
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
# Step 2: Merge VCFs using BCFtools merge (skip if already exists)
##########################################################################################################################

module load BCFtools

if [ -f "${RAW_MERGED}" ]; then
    echo "Step 2: Raw merged VCF already exists, skipping merge step..."
    echo "Using existing: ${RAW_MERGED}"
    echo ""
else
    echo "Step 2: Merging VCFs with BCFtools..."

    # Merge all VCFs into a single multi-sample VCF
    # Options:
    #   -m none: do not merge records (keep sites separate if they differ across samples)
    #   -O z: output compressed VCF
    #   --threads: use multiple threads for compression
    bcftools merge \
        --threads 31 \
        -m none \
        -O z \
        -o "${RAW_MERGED}" \
        --file-list "${VCF_LIST}"

    echo "VCF merge complete. Output: ${RAW_MERGED}"
    echo ""
fi

##########################################################################################################################
# Step 3: Normalize merged VCF (split multiallelic sites, left-align indels) and filter
##########################################################################################################################

echo "Step 3: Normalizing merged VCF and filtering out reference-only sites..."
echo "Using 32 threads"

# Normalize and filter out sites with ALT = "."
# Pipe: normalize -> filter out ALT="." -> compress
bcftools norm \
    --threads 16 \
    -m -both \
    -f "${REFERENCE}" \
    -O u \
    "${RAW_MERGED}" | \
bcftools view \
    --threads 15 \
    -e 'ALT="."' \
    -O z \
    -o "${OUTPUT_VCF}"

echo "Normalization and filtering complete."
echo ""

##########################################################################################################################
# Step 4: Clean up and index
##########################################################################################################################

echo "Step 4: Cleaning up and indexing..."

# Remove raw merged VCF (keep it to avoid regenerating if job needs to be rerun)
# rm "${RAW_MERGED}"

# Index normalized VCF
bcftools index --threads 31 "${OUTPUT_VCF}"

echo "Indexing complete."
echo ""

##########################################################################################################################
# Step 5: Generate statistics
##########################################################################################################################

echo "Step 5: Generating variant statistics..."

bcftools stats "${OUTPUT_VCF}" > "${OUTPUT_VCF%.vcf.gz}_stats.txt"

echo ""
echo "============================================"
echo "BCFtools merge and normalize complete!"
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
