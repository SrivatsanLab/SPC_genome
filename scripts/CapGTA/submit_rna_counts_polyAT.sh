#!/bin/bash
###########################################################################################################################
# Submit the polyA/T-only (no UMI dedup) RNA counting pipeline.
#
# Companion to submit_rna_counts_umi.sh. Three stages:
#   1. umi_tag_filter_array.sh --no-require-umi:   SLURM array, one task per cell. Keep reads
#                                                   with polyA/T >= run_min at inner SC edge
#                                                   regardless of whether a UMI is recoverable.
#   2. featureCounts_assign.sh:                     featureCounts with -R BAM — same wrapper as
#                                                   the UMI pipeline.
#   3. polyAT_count.sh:                             build h5ad from the featureCounts raw matrix
#                                                   (cells × genes read counts, no dedup).
#
# Usage:
#   scripts/CapGTA/submit_rna_counts_polyAT.sh \
#       <real_cells.csv> \
#       <annotation.gtf> \
#       <output_dir> \
#       [run_min=6]
#
# Layout under <output_dir>:
#   bam_list_source.txt / bam_list_umi.txt
#   umi_bams/<barcode>.bam[.bai]
#   filter_stats.tsv
#   featurecounts_annotated/<barcode>.bam.featureCounts.bam
#   counts_raw.txt / counts_raw.txt.summary
#   rna_polyAT_matrix.h5ad                 # cells × genes read counts, filter_params in uns
###########################################################################################################################

set -euo pipefail

if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <real_cells.csv> <annotation.gtf> <output_dir> [run_min=6]" >&2
    exit 1
fi

CELLS_CSV="$1"
GTF="$2"
OUTPUT_DIR="$3"
RUN_MIN="${4:-6}"

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for f in "${CELLS_CSV}" "${GTF}"; do
    if [ ! -f "${f}" ]; then
        echo "Error: required input not found: ${f}" >&2
        exit 1
    fi
done

mkdir -p "${OUTPUT_DIR}" SLURM_outs SLURM_outs/array_outs
BAM_LIST_SRC="${OUTPUT_DIR}/bam_list_source.txt"
UMI_DIR="${OUTPUT_DIR}/umi_bams"
BAM_LIST_UMI="${OUTPUT_DIR}/bam_list_umi.txt"
STATS_TSV="${OUTPUT_DIR}/filter_stats.tsv"
ANNO_DIR="${OUTPUT_DIR}/featurecounts_annotated"
COUNTS_RAW="${OUTPUT_DIR}/counts_raw.txt"
H5AD="${OUTPUT_DIR}/rna_polyAT_matrix.h5ad"

mkdir -p "${UMI_DIR}" "${ANNO_DIR}"

HEADER=$(head -1 "${CELLS_CSV}")
if [[ "${HEADER}" != "barcode,bam_path" ]]; then
    echo "Warning: expected header 'barcode,bam_path' in ${CELLS_CSV} (got '${HEADER}'); assuming column 2 is bam_path." >&2
fi
awk -F, 'NR > 1 && $2 != "" { print $2 }' "${CELLS_CSV}" > "${BAM_LIST_SRC}"
N_BAMS=$(wc -l < "${BAM_LIST_SRC}")
if [ "${N_BAMS}" -eq 0 ]; then
    echo "Error: no BAM paths parsed from ${CELLS_CSV}" >&2
    exit 1
fi
awk -v d="${UMI_DIR}" '{
    n = split($0, parts, "/")
    print d "/" parts[n]
}' "${BAM_LIST_SRC}" > "${BAM_LIST_UMI}"

rm -f "${STATS_TSV}"

echo "=========================================="
echo "polyA/T-only (no UMI) RNA counts submission"
echo "=========================================="
echo "Cells CSV:      ${CELLS_CSV}"
echo "Cells (BAMs):   ${N_BAMS}"
echo "GTF:            ${GTF}"
echo "Output dir:     ${OUTPUT_DIR}"
echo "run_min (L):    ${RUN_MIN}"
echo "require_umi:    0  (--no-require-umi)"
echo "=========================================="

# --- 1) Filter array (polyA/T only, no UMI required) ------------------------
filter_job=$(sbatch --parsable \
    --array="1-${N_BAMS}" \
    "${SCRIPTS_DIR}/umi_tag_filter_array.sh" \
    "${BAM_LIST_SRC}" \
    "${UMI_DIR}" \
    "${RUN_MIN}" \
    0)
echo "Submitted filter array:         ${filter_job} (1-${N_BAMS})"

# --- 2) featureCounts gene assignment ----------------------------------------
fc_job=$(sbatch --parsable \
    --dependency="afterok:${filter_job}" \
    "${SCRIPTS_DIR}/featureCounts_assign.sh" \
    "${BAM_LIST_UMI}" \
    "${GTF}" \
    "${OUTPUT_DIR}")
echo "Submitted featureCounts job:    ${fc_job} (depends on ${filter_job})"

# --- 3) h5ad build (read counts from featureCounts raw matrix) ---------------
count_job=$(sbatch --parsable \
    --dependency="afterok:${fc_job}" \
    "${SCRIPTS_DIR}/polyAT_count.sh" \
    "${COUNTS_RAW}" "${STATS_TSV}" "${RUN_MIN}" "${H5AD}")
echo "Submitted polyAT_count job:     ${count_job} (depends on ${fc_job})"

echo ""
echo "Monitor with:  squeue -j ${filter_job},${fc_job},${count_job}"
