#!/bin/bash
###########################################################################################################################
# Submit RNA count matrix generation for a set of single-cell BAMs.
#
# Two-stage pipeline:
#   1. filter_spliced_reads_array.sh — SLURM array, one task per cell, writes spliced-only BAMs.
#   2. create_rna_count_matrix.sh    — depends afterok on the array; runs featureCounts on the
#                                      filtered BAMs and writes a gene_id x barcode matrix.
#
# Usage:
#   scripts/CapGTA/submit_rna_counts.sh \
#       <real_cells.csv> \
#       <annotation.gtf> \
#       <output_dir>
#
# Layout under <output_dir>:
#   bam_list_source.txt       # bam_path column from real_cells.csv
#   spliced_bams/{basename}.bam[.bai]
#   bam_list_spliced.txt      # paths of filtered BAMs (same order/basename as source)
#   rna_counts_raw.txt
#   rna_counts_matrix.csv
#   rna_counts_summary.csv
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <real_cells.csv> <annotation.gtf> <output_dir>" >&2
    exit 1
fi

CELLS_CSV="$1"
GTF="$2"
OUTPUT_DIR="$3"

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for f in "${CELLS_CSV}" "${GTF}"; do
    if [ ! -f "${f}" ]; then
        echo "Error: required input not found: ${f}" >&2
        exit 1
    fi
done

SPLICED_DIR="${OUTPUT_DIR}/spliced_bams"
BAM_LIST_SRC="${OUTPUT_DIR}/bam_list_source.txt"
BAM_LIST_SPLICED="${OUTPUT_DIR}/bam_list_spliced.txt"
OUTPUT_PREFIX="${OUTPUT_DIR}/rna_counts"

mkdir -p "${OUTPUT_DIR}" "${SPLICED_DIR}" SLURM_outs SLURM_outs/array_outs

# --- Build source bam_list from real_cells.csv --------------------------------
HEADER=$(head -1 "${CELLS_CSV}")
if [[ "${HEADER}" != "barcode,bam_path" ]]; then
    echo "Warning: expected header 'barcode,bam_path' in ${CELLS_CSV} (got '${HEADER}'); assuming column 2 is bam_path anyway." >&2
fi

awk -F, 'NR > 1 && $2 != "" { print $2 }' "${CELLS_CSV}" > "${BAM_LIST_SRC}"
N_BAMS=$(wc -l < "${BAM_LIST_SRC}")
if [ "${N_BAMS}" -eq 0 ]; then
    echo "Error: no BAM paths parsed from ${CELLS_CSV}" >&2
    exit 1
fi

# Spot-check first/last BAMs exist.
for bam in "$(head -1 "${BAM_LIST_SRC}")" "$(tail -1 "${BAM_LIST_SRC}")"; do
    if [ ! -f "${bam}" ]; then
        echo "Error: source BAM not found (spot-check): ${bam}" >&2
        exit 1
    fi
done

# --- Pre-compute spliced-BAM paths (same basenames, in SPLICED_DIR) ----------
awk -v d="${SPLICED_DIR}" '{
    n = split($0, parts, "/")
    print d "/" parts[n]
}' "${BAM_LIST_SRC}" > "${BAM_LIST_SPLICED}"

echo "=========================================="
echo "RNA counts submission"
echo "=========================================="
echo "Cells CSV:      ${CELLS_CSV}"
echo "Cells (BAMs):   ${N_BAMS}"
echo "GTF:            ${GTF}"
echo "Spliced BAMs:   ${SPLICED_DIR}"
echo "Output prefix:  ${OUTPUT_PREFIX}"
echo "=========================================="

# --- Submit spliced-filter array ---------------------------------------------
filter_job=$(sbatch --parsable \
    --array="1-${N_BAMS}" \
    "${SCRIPTS_DIR}/filter_spliced_reads_array.sh" \
    "${BAM_LIST_SRC}" \
    "${SPLICED_DIR}")

echo "Submitted spliced-filter array: ${filter_job} (1-${N_BAMS})"

# --- Submit counts job with dependency ---------------------------------------
counts_job=$(sbatch --parsable \
    --dependency="afterok:${filter_job}" \
    "${SCRIPTS_DIR}/create_rna_count_matrix.sh" \
    "${BAM_LIST_SPLICED}" \
    "${GTF}" \
    "${OUTPUT_PREFIX}")

echo "Submitted count-matrix job:     ${counts_job} (depends on ${filter_job})"
echo ""
echo "Monitor with:  squeue -j ${filter_job},${counts_job}"
