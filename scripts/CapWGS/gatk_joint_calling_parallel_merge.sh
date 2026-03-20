#!/bin/bash
#SBATCH -J merge_jc_vcfs
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 12:00:00

##########################################################################################################################
# Merge per-interval VCFs into a single joint-called VCF for CapWGS pipeline
# This runs after all per-interval joint calling jobs complete
# Handles both chromosome-based (legacy) and interval-based (new) parallelization
##########################################################################################################################

set -euo pipefail

# Arguments
PER_INTERVAL_DIR="$1"   # Directory containing per-interval VCFs
OUTPUT_VCF="$2"         # Final merged VCF output
INTERVAL_LIST="$3"      # File with one interval file path per line (or chromosome names for legacy mode)

echo "============================================"
echo "Merging CapWGS per-interval VCFs"
echo "============================================"
echo ""
echo "Input directory: ${PER_INTERVAL_DIR}"
echo "Output VCF: ${OUTPUT_VCF}"
echo ""

# Check that all per-interval VCFs exist
echo "Checking for per-interval VCFs..."
MISSING=0
while read -r interval_path; do
    # Extract interval name from file path (e.g., "0042-scattered" from "/path/0042-scattered.interval_list")
    # Or just use the name directly if it's a chromosome name (legacy mode)
    if [[ "$interval_path" == *.interval_list ]]; then
        # New interval-based mode
        interval_name=$(basename "$interval_path" .interval_list)
    else
        # Legacy chromosome-based mode
        interval_name="$interval_path"
    fi

    if [ ! -f "${PER_INTERVAL_DIR}/${interval_name}.vcf.gz" ]; then
        echo "ERROR: Missing VCF for ${interval_name}"
        MISSING=1
    fi
done < "${INTERVAL_LIST}"

if [ ${MISSING} -eq 1 ]; then
    echo "ERROR: Not all per-interval VCFs are present. Cannot merge."
    exit 1
fi

echo "All per-interval VCFs found."
echo ""

##########################################################################################################################
# Merge VCFs with bcftools concat
##########################################################################################################################

echo "Merging VCFs with bcftools concat..."
echo "Started at: $(date)"
echo ""

module load BCFtools

# Create list of VCF files in genomic order
VCF_LIST="${PER_INTERVAL_DIR}/vcf_list.txt"
> "${VCF_LIST}"

while read -r interval_path; do
    # Extract interval name from file path (or use chromosome name for legacy mode)
    if [[ "$interval_path" == *.interval_list ]]; then
        # New interval-based mode
        interval_name=$(basename "$interval_path" .interval_list)
    else
        # Legacy chromosome-based mode
        interval_name="$interval_path"
    fi

    echo "${PER_INTERVAL_DIR}/${interval_name}.vcf.gz" >> "${VCF_LIST}"
done < "${INTERVAL_LIST}"

# Concatenate VCFs (they're already sorted by genomic position)
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
