#!/bin/bash
#SBATCH -J concat_pb
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -p campus-new
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 4:00:00

###########################################################################################################################
# Concat per-region joint VCFs into one worm6 pseudobulk joint VCF.
#
# Usage:
#   sbatch scripts/CapGTA/pseudobulk_gatk/concat_and_final.sh <output_dir> <regions_file>
#
# Writes:
#   <output_dir>/joint.vcf.gz[.csi]
#   <output_dir>/joint_vcf_list.txt
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: sbatch $0 <output_dir> <regions_file>" >&2
    exit 1
fi

OUT="$1"
REGIONS_FILE="$2"

[ -f "${REGIONS_FILE}" ] || { echo "Error: regions_file not found: ${REGIONS_FILE}" >&2; exit 1; }
[ -d "${OUT}/joint" ]     || { echo "Error: joint region dir not found: ${OUT}/joint" >&2; exit 1; }

# Build the concat list in the same order as regions_file (matches the array task ordering).
VCF_LIST="${OUT}/joint_vcf_list.txt"
: > "${VCF_LIST}"
while read -r REGION; do
    [ -z "${REGION}" ] && continue
    REGION_TAG="${REGION//:/_}"
    REGION_TAG="${REGION_TAG//-/_}"
    VCF="${OUT}/joint/${REGION_TAG}.vcf.gz"
    if [ ! -f "${VCF}" ]; then
        echo "Error: missing joint VCF for region ${REGION}: ${VCF}" >&2
        exit 1
    fi
    echo "${VCF}" >> "${VCF_LIST}"
done < "${REGIONS_FILE}"

N=$(wc -l < "${VCF_LIST}")
FINAL_VCF="${OUT}/joint.vcf.gz"

echo "=========================================="
echo "bcftools concat  n_regions=${N}"
echo "started=$(date -Iseconds)"
echo "=========================================="

module load BCFtools

bcftools concat -a --threads 4 -O z -f "${VCF_LIST}" -o "${FINAL_VCF}"
# Both indices: .csi for bcftools' preferred format, .tbi because GATK requires it.
bcftools index --threads 4 -c "${FINAL_VCF}"
bcftools index --threads 4 -t "${FINAL_VCF}"

N_VAR=$(bcftools view -H "${FINAL_VCF}" | wc -l)
N_SAMP=$(bcftools query -l "${FINAL_VCF}" | wc -l)

echo ""
echo "=========================================="
echo "Done  finished=$(date -Iseconds)"
echo "Final VCF: ${FINAL_VCF}"
echo "Variants:  ${N_VAR}"
echo "Samples:   ${N_SAMP}"
echo "=========================================="
