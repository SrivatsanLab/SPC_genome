#!/bin/bash
#SBATCH -J hc_pb
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -p campus-new
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 8:00:00

###########################################################################################################################
# Per-worm × per-region GATK HaplotypeCaller in GVCF mode.
#
# Task decoding (0-based indexing internally):
#   t          = SLURM_ARRAY_TASK_ID - 1
#   n_regions  = wc -l < regions.txt
#   worm_idx   = t / n_regions      -> WORM   = passing_worms[worm_idx]
#   region_idx = t % n_regions      -> REGION = regions[region_idx]
#
# Contamination handling:
#   HC's --contamination-fraction-to-filter is set to the worm's med_alpha from worm_summary.csv.
#   This is a coarse GATK-side filter; the real contamination correction is downstream in
#   regenotype_from_ad.py (see plan §5).
#
# Usage:
#   N_WORMS=$(wc -l < <output_dir>/worm_lists/passing_worms.txt)
#   N_REG=$(wc -l < <regions_file>)
#   sbatch --array=1-$((N_WORMS * N_REG)) \
#       scripts/CapGTA/pseudobulk_gatk/haplotype_caller_array.sh <output_dir> <reference_fa> <regions_file>
#
# Writes:
#   <output_dir>/gvcf/<worm>/<region_tag>.g.vcf.gz[.tbi]
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: sbatch --array=1-N $0 <output_dir> <reference_fa> <regions_file>" >&2
    exit 1
fi

OUT="$1"
REF="$2"
REGIONS_FILE="$3"

for f in "${REF}" "${REF}.fai" "${REGIONS_FILE}" \
         "${OUT}/worm_lists/passing_worms.txt" "${OUT}/worm_summary.csv"; do
    [ -f "${f}" ] || { echo "Error: required input not found: ${f}" >&2; exit 1; }
done
REF_DICT="${REF%.fa}.dict"
[ -f "${REF_DICT}" ] || REF_DICT="${REF%.fasta}.dict"
[ -f "${REF_DICT}" ] || { echo "Error: reference dict not found alongside ${REF}" >&2; exit 1; }

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID not set — submit as an array." >&2
    exit 1
fi

N_REG=$(wc -l < "${REGIONS_FILE}")
T=$((SLURM_ARRAY_TASK_ID - 1))
WORM_IDX=$((T / N_REG))
REGION_IDX=$((T % N_REG))

WORM=$(sed -n "$((WORM_IDX + 1))p" "${OUT}/worm_lists/passing_worms.txt")
REGION=$(sed -n "$((REGION_IDX + 1))p" "${REGIONS_FILE}")

if [ -z "${WORM}" ] || [ -z "${REGION}" ]; then
    echo "Error: task_id=${SLURM_ARRAY_TASK_ID} decoded to empty (worm='${WORM}', region='${REGION}')" >&2
    exit 1
fi

WORM_BAM="${OUT}/worm_bams/${WORM}.bam"
[ -f "${WORM_BAM}" ] || { echo "Error: worm BAM not found: ${WORM_BAM}" >&2; exit 1; }

# Look up the worm's median alpha_hap from worm_summary.csv (skip header row).
ALPHA=$(awk -F, -v w="${WORM}" 'NR>1 && $1==w {print $3; exit}' "${OUT}/worm_summary.csv")
if [ -z "${ALPHA}" ]; then
    echo "Error: could not find med_alpha for worm '${WORM}' in ${OUT}/worm_summary.csv" >&2
    exit 1
fi

# Sanitize region for filesystem (`:` and `-` are legal on POSIX but noisy).
REGION_TAG="${REGION//:/_}"
REGION_TAG="${REGION_TAG//-/_}"

GVCF_DIR="${OUT}/gvcf/${WORM}"
mkdir -p "${GVCF_DIR}"
OUT_GVCF="${GVCF_DIR}/${REGION_TAG}.g.vcf.gz"

echo "=========================================="
echo "HaplotypeCaller   worm=${WORM}   region=${REGION}"
echo "task_id=${SLURM_ARRAY_TASK_ID}  worm_idx=${WORM_IDX}  region_idx=${REGION_IDX}"
echo "alpha=${ALPHA}  started=$(date -Iseconds)"
echo "=========================================="

module load GATK

gatk --java-options "-Xmx14g" HaplotypeCaller \
    -R "${REF}" \
    -I "${WORM_BAM}" \
    -O "${OUT_GVCF}" \
    -L "${REGION}" \
    -ERC GVCF \
    --sample-ploidy 2 \
    --contamination-fraction-to-filter "${ALPHA}" \
    --min-base-quality-score 20 \
    --disable-read-filter MappingQualityAvailableReadFilter
# STAR sets MAPQ=255 for uniquely-mapped reads (with --outFilterMultimapNmax 1),
# but GATK treats 255 as the "MAPQ not available" sentinel and drops every such
# read via MappingQualityAvailableReadFilter — leaving HC 0 reads to call from.
# Disabling this filter is the standard fix for STAR-aligned input.

echo ""
echo "=========================================="
echo "Done: ${WORM} ${REGION}  finished=$(date -Iseconds)"
echo "Output: ${OUT_GVCF}"
echo "=========================================="
