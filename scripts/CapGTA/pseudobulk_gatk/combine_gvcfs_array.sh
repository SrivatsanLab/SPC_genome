#!/bin/bash
#SBATCH -J comb_gt
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -p campus-new
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 4:00:00

###########################################################################################################################
# CombineGVCFs across worms + GenotypeGVCFs, one SLURM array task per region.
#
# Joint genotyping across worms is important: it forces every worm to be genotyped at every site,
# giving a complete sites x worms matrix rather than a ragged union of per-worm calls.
#
# CombineGVCFs is used (not GenomicsDBImport) because worm6 has ~17 worms — well below the
# ~30-sample threshold where GenomicsDBImport wins.
#
# Usage:
#   sbatch --array=1-$(wc -l < <regions_file>) \
#       scripts/CapGTA/pseudobulk_gatk/combine_gvcfs_array.sh <output_dir> <reference_fa> <regions_file>
#
# Writes:
#   <output_dir>/combined_gvcf/<region_tag>.g.vcf.gz
#   <output_dir>/joint/<region_tag>.vcf.gz
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: sbatch --array=1-N $0 <output_dir> <reference_fa> <regions_file>" >&2
    exit 1
fi

OUT="$1"
REF="$2"
REGIONS_FILE="$3"

for f in "${REF}" "${REF}.fai" "${REGIONS_FILE}" "${OUT}/worm_lists/passing_worms.txt"; do
    [ -f "${f}" ] || { echo "Error: required input not found: ${f}" >&2; exit 1; }
done

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID not set — submit as an array." >&2
    exit 1
fi

REGION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${REGIONS_FILE}")
if [ -z "${REGION}" ]; then
    echo "Error: no region at line ${SLURM_ARRAY_TASK_ID} of ${REGIONS_FILE}" >&2
    exit 1
fi

REGION_TAG="${REGION//:/_}"
REGION_TAG="${REGION_TAG//-/_}"

COMB_DIR="${OUT}/combined_gvcf"
JOINT_DIR="${OUT}/joint"
mkdir -p "${COMB_DIR}" "${JOINT_DIR}"

COMB_GVCF="${COMB_DIR}/${REGION_TAG}.g.vcf.gz"
JOINT_VCF="${JOINT_DIR}/${REGION_TAG}.vcf.gz"

# Assemble -V args from all worm GVCFs for this region.
V_ARGS=()
while read -r WORM; do
    [ -z "${WORM}" ] && continue
    GVCF="${OUT}/gvcf/${WORM}/${REGION_TAG}.g.vcf.gz"
    if [ ! -f "${GVCF}" ]; then
        echo "Error: missing per-worm GVCF for ${WORM} at ${REGION}: ${GVCF}" >&2
        exit 1
    fi
    V_ARGS+=("-V" "${GVCF}")
done < "${OUT}/worm_lists/passing_worms.txt"

echo "=========================================="
echo "CombineGVCFs + GenotypeGVCFs   region=${REGION}"
echo "task_id=${SLURM_ARRAY_TASK_ID}  n_gvcfs=$(( ${#V_ARGS[@]} / 2 ))"
echo "started=$(date -Iseconds)"
echo "=========================================="

module load GATK

gatk --java-options "-Xmx28g" CombineGVCFs \
    -R "${REF}" \
    "${V_ARGS[@]}" \
    -O "${COMB_GVCF}"

gatk --java-options "-Xmx28g" GenotypeGVCFs \
    -R "${REF}" \
    -V "${COMB_GVCF}" \
    -O "${JOINT_VCF}" \
    --sample-ploidy 2

echo ""
echo "=========================================="
echo "Done: ${REGION}  finished=$(date -Iseconds)"
echo "Combined GVCF: ${COMB_GVCF}"
echo "Joint VCF:     ${JOINT_VCF}"
echo "=========================================="
