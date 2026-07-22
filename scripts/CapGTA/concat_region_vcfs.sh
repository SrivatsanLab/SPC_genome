#!/bin/bash
#SBATCH -J concat_vcfs
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -p campus-new
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 4:00:00

###########################################################################################################################
# Concatenate per-region joint VCFs (from joint_variant_calling_array.sh) into a single joint VCF.
#
# Usage:
#   sbatch scripts/CapGTA/concat_region_vcfs.sh <per_region_vcf_dir> <output.vcf.gz>
#
# The per-region directory is expected to contain one or more `region_*.vcf.gz` files.
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <per_region_vcf_dir> <output.vcf.gz>" >&2
    exit 1
fi

VCF_DIR="$1"
OUTPUT_VCF="$2"

if [ ! -d "${VCF_DIR}" ]; then
    echo "Error: per-region VCF dir not found: ${VCF_DIR}" >&2
    exit 1
fi

mkdir -p "$(dirname "${OUTPUT_VCF}")"

module load BCFtools

FILE_LIST="${OUTPUT_VCF%.vcf.gz}_file_list.txt"

# Sort lexically — region_000001.vcf.gz, region_000002.vcf.gz, ...
find "${VCF_DIR}" -maxdepth 1 -name 'region_*.vcf.gz' -print | sort > "${FILE_LIST}"

N_VCFS=$(wc -l < "${FILE_LIST}")
if [ "${N_VCFS}" -eq 0 ]; then
    echo "Error: no region_*.vcf.gz files found in ${VCF_DIR}" >&2
    exit 1
fi

echo "=========================================="
echo "Concatenating ${N_VCFS} per-region VCFs"
echo "=========================================="
echo "Input dir:  ${VCF_DIR}"
echo "Output:     ${OUTPUT_VCF}"
echo "Started:    $(date -Iseconds)"

# -a: allow overlaps and re-sort (safe if regions are adjacent, no overlap by construction)
# --naive is faster but assumes identical headers with no overlap; we use -a to be robust.
bcftools concat \
    --threads 4 \
    -a \
    -O z \
    -o "${OUTPUT_VCF}" \
    --file-list "${FILE_LIST}"

bcftools index "${OUTPUT_VCF}"

rm -f "${FILE_LIST}"

echo "Finished:   $(date -Iseconds)"
echo ""
echo "=== Joint VCF Summary ==="
N_SAMPLES=$(bcftools query -l "${OUTPUT_VCF}" | wc -l)
N_VARIANTS=$(bcftools view -H "${OUTPUT_VCF}" | wc -l)
N_SNPS=$(bcftools view -H -v snps "${OUTPUT_VCF}" | wc -l)
N_INDELS=$(bcftools view -H -v indels "${OUTPUT_VCF}" | wc -l)
echo "Samples:  ${N_SAMPLES}"
echo "Variants: ${N_VARIANTS}"
echo "  SNPs:   ${N_SNPS}"
echo "  Indels: ${N_INDELS}"
echo ""
echo "Output: ${OUTPUT_VCF}"
