#!/bin/bash
###########################################################################################################################
# Dev variant of scripts/CapGTA/submit_rna_counts.sh — filters BAMs to the UNION of spliced and polyA/T-containing
# reads (instead of spliced-only) before counting.
#
# Three-stage submission (chained via SLURM afterok):
#   1. filter_spliced_polyat_reads_array.sh — SLURM array, one task per cell, writes filtered BAMs.
#   2. create_rna_count_matrix.sh           — depends afterok on the array; runs featureCounts on the
#                                             filtered BAMs and writes a gene_id x barcode matrix.
#   3. counts_matrix_to_h5ad.py             — depends afterok on the counts job; produces sparse .h5ad.
#
# Usage:
#   CapGTA_dev/scripts/submit_rna_counts_spliced_polyat.sh \
#       <real_cells.csv> \
#       <annotation.gtf> \
#       <output_dir>
#
# Layout under <output_dir>:
#   bam_list_source.txt              # bam_path column from real_cells.csv
#   splicedpolyat_bams/{basename}.bam[.bai]
#   bam_list_splicedpolyat.txt
#   rna_counts_raw.txt
#   rna_counts_matrix.csv
#   rna_counts_summary.csv
#   rna_counts.h5ad
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <real_cells.csv> <annotation.gtf> <output_dir>" >&2
    exit 1
fi

CELLS_CSV="$1"
GTF="$2"
OUTPUT_DIR="$3"

DEV_SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN_SCRIPTS="$(cd "${DEV_SCRIPTS}/../.." && pwd)/scripts/CapGTA"

for f in "${CELLS_CSV}" "${GTF}"; do
    if [ ! -f "${f}" ]; then
        echo "Error: required input not found: ${f}" >&2
        exit 1
    fi
done

BAM_LIST_SRC="${OUTPUT_DIR}/bam_list_source.txt"
FILTERED_DIR="${OUTPUT_DIR}/splicedpolyat_bams"
BAM_LIST_FILTERED="${OUTPUT_DIR}/bam_list_splicedpolyat.txt"
OUTPUT_PREFIX="${OUTPUT_DIR}/rna_counts"

mkdir -p "${OUTPUT_DIR}" "${FILTERED_DIR}" SLURM_outs SLURM_outs/array_outs

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

for bam in "$(head -1 "${BAM_LIST_SRC}")" "$(tail -1 "${BAM_LIST_SRC}")"; do
    if [ ! -f "${bam}" ]; then
        echo "Error: source BAM not found (spot-check): ${bam}" >&2
        exit 1
    fi
done

# Pre-compute filtered-BAM paths (same basenames, in FILTERED_DIR).
awk -v d="${FILTERED_DIR}" '{
    n = split($0, parts, "/")
    print d "/" parts[n]
}' "${BAM_LIST_SRC}" > "${BAM_LIST_FILTERED}"

echo "=========================================="
echo "RNA counts submission (spliced ∪ polyA/T)"
echo "=========================================="
echo "Cells CSV:      ${CELLS_CSV}"
echo "Cells (BAMs):   ${N_BAMS}"
echo "GTF:            ${GTF}"
echo "Output prefix:  ${OUTPUT_PREFIX}"
echo "Filtered BAMs:  ${FILTERED_DIR}"
echo "=========================================="

filter_job=$(sbatch --parsable \
    --array="1-${N_BAMS}" \
    "${DEV_SCRIPTS}/filter_spliced_polyat_reads_array.sh" \
    "${BAM_LIST_SRC}" \
    "${FILTERED_DIR}")

echo "Submitted spliced∪polyat filter array: ${filter_job} (1-${N_BAMS})"

counts_job=$(sbatch --parsable \
    --dependency="afterok:${filter_job}" \
    "${MAIN_SCRIPTS}/create_rna_count_matrix.sh" \
    "${BAM_LIST_FILTERED}" \
    "${GTF}" \
    "${OUTPUT_PREFIX}")

echo "Submitted count-matrix job:            ${counts_job} (depends on ${filter_job})"

# --- h5ad conversion with dependency ------------------------------------------
h5ad_job=$(sbatch --parsable \
    --dependency="afterok:${counts_job}" \
    -J counts_to_h5ad \
    -o SLURM_outs/%x_%j.out \
    -p campus-new -c 1 --mem=8G -t 30:00 \
    --wrap="python3 '${MAIN_SCRIPTS}/counts_matrix_to_h5ad.py' \
        '${OUTPUT_PREFIX}_matrix.csv' \
        '${OUTPUT_PREFIX}.h5ad' \
        --summary '${OUTPUT_PREFIX}_summary.csv'")

echo "Submitted h5ad conversion:             ${h5ad_job} (depends on ${counts_job})"
echo ""
echo "Monitor with:  squeue -j ${filter_job},${counts_job},${h5ad_job}"
