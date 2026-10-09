#!/bin/bash
###########################################################################################################################
# Submit joint variant calling for a set of single-cell BAMs.
#
# Reads a real_cells.csv (columns: barcode, bam_path) and submits:
#   1. joint_variant_calling_array.sh — one SLURM task per genomic region (from make_regions.sh)
#   2. concat_region_vcfs.sh          — depends afterok on the array; stitches per-region VCFs
#
# Usage:
#   bash CapGTA_joint_variant_calling.sh \
#       <real_cells.csv> \
#       <reference.fa> \
#       <output_dir> \
#       [region_size_bp]      # default: 1000000
#
# Layout under <output_dir>:
#   bam_list.txt              (derived from real_cells.csv bam_path column)
#   regions.txt               (one bcftools -r region per line)
#   per_region/region_000001.vcf.gz ...
#   joint_variants.vcf.gz     (final)
###########################################################################################################################

set -euo pipefail

if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <real_cells.csv> <reference.fa> <output_dir> [region_size_bp]" >&2
    exit 1
fi

CELLS_CSV="$1"
REFERENCE_FA="$2"
OUTPUT_DIR="$3"
REGION_SIZE="${4:-1000000}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${REPO_ROOT}/scripts/CapGTA"

for f in "${CELLS_CSV}" "${REFERENCE_FA}" "${REFERENCE_FA}.fai"; do
    if [ ! -f "${f}" ]; then
        echo "Error: required input not found: ${f}" >&2
        exit 1
    fi
done

mkdir -p "${OUTPUT_DIR}" "${OUTPUT_DIR}/per_region" SLURM_outs SLURM_outs/array_outs

BAM_LIST="${OUTPUT_DIR}/bam_list.txt"
REGIONS_FILE="${OUTPUT_DIR}/regions.txt"
PER_REGION_DIR="${OUTPUT_DIR}/per_region"
FINAL_VCF="${OUTPUT_DIR}/joint_variants.vcf.gz"

# --- Build bam_list from real_cells.csv ---------------------------------------
# Expect header: barcode,bam_path
HEADER=$(head -1 "${CELLS_CSV}")
if [[ "${HEADER}" != "barcode,bam_path" ]]; then
    echo "Warning: expected header 'barcode,bam_path' in ${CELLS_CSV} (got '${HEADER}'); assuming column 2 is bam_path anyway." >&2
fi

awk -F, 'NR > 1 && $2 != "" { print $2 }' "${CELLS_CSV}" > "${BAM_LIST}"
N_BAMS=$(wc -l < "${BAM_LIST}")
if [ "${N_BAMS}" -eq 0 ]; then
    echo "Error: no BAM paths parsed from ${CELLS_CSV}" >&2
    exit 1
fi

# Sanity-check the first and last BAM exist to catch a stale CSV early.
FIRST_BAM=$(head -1 "${BAM_LIST}")
LAST_BAM=$(tail -1 "${BAM_LIST}")
for bam in "${FIRST_BAM}" "${LAST_BAM}"; do
    if [ ! -f "${bam}" ]; then
        echo "Error: BAM not found (spot-check failed): ${bam}" >&2
        exit 1
    fi
done

# --- Generate regions ---------------------------------------------------------
bash "${SCRIPTS_DIR}/make_regions.sh" "${REFERENCE_FA}" "${REGION_SIZE}" "${REGIONS_FILE}"
N_REGIONS=$(wc -l < "${REGIONS_FILE}")

echo "=========================================="
echo "Joint variant calling submission"
echo "=========================================="
echo "Cells CSV:    ${CELLS_CSV}"
echo "Cells (BAMs): ${N_BAMS}"
echo "Reference:    ${REFERENCE_FA}"
echo "Region size:  ${REGION_SIZE} bp"
echo "Regions:      ${N_REGIONS}"
echo "Output dir:   ${OUTPUT_DIR}"
echo "Final VCF:    ${FINAL_VCF}"
echo "=========================================="

# --- Submit joint-calling array ----------------------------------------------
array_job=$(sbatch --parsable \
    --array="1-${N_REGIONS}" \
    "${SCRIPTS_DIR}/joint_variant_calling_array.sh" \
    "${BAM_LIST}" \
    "${REGIONS_FILE}" \
    "${REFERENCE_FA}" \
    "${PER_REGION_DIR}")

echo "Submitted joint-calling array:  ${array_job} (1-${N_REGIONS})"

# --- Submit concat with dependency -------------------------------------------
concat_job=$(sbatch --parsable \
    --dependency="afterok:${array_job}" \
    "${SCRIPTS_DIR}/concat_region_vcfs.sh" \
    "${PER_REGION_DIR}" \
    "${FINAL_VCF}")

echo "Submitted concat job:           ${concat_job} (depends on ${array_job})"
echo ""
echo "Monitor with:  squeue -j ${array_job},${concat_job}"
