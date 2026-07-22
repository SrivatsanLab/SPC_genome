#!/bin/bash
###########################################################################################################################
# Submit RNA count matrix generation for a set of single-cell BAMs.
#
# By default (spliced-only mode), three-stage pipeline:
#   1. filter_spliced_reads_array.sh — SLURM array, one task per cell, writes spliced-only BAMs.
#   2. create_rna_count_matrix.sh    — depends afterok on the array; runs featureCounts on the
#                                      filtered BAMs and writes a gene_id x barcode matrix.
#   3. counts_matrix_to_h5ad.py      — depends afterok on the counts job; converts the matrix
#                                      CSV into a sparse .h5ad AnnData for scanpy.
#
# With --all-reads: the spliced-filter step is skipped and featureCounts runs directly on the
# source BAMs. Output layout is the same minus spliced_bams/ and bam_list_spliced.txt.
#
# Usage:
#   scripts/CapGTA/submit_rna_counts.sh [--all-reads] \
#       <real_cells.csv> \
#       <annotation.gtf> \
#       <output_dir>
#
# Layout under <output_dir>:
#   bam_list_source.txt              # bam_path column from real_cells.csv
#   spliced_bams/{basename}.bam[.bai]      (spliced-only mode)
#   bam_list_spliced.txt                    (spliced-only mode)
#   rna_counts_raw.txt
#   rna_counts_matrix.csv
#   rna_counts_summary.csv
#   rna_counts.h5ad                  # scanpy-ready AnnData
###########################################################################################################################

set -euo pipefail

FILTER_SPLICED=1
if [ "${1:-}" = "--all-reads" ]; then
    FILTER_SPLICED=0
    shift
fi

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 [--all-reads] <real_cells.csv> <annotation.gtf> <output_dir>" >&2
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

BAM_LIST_SRC="${OUTPUT_DIR}/bam_list_source.txt"
OUTPUT_PREFIX="${OUTPUT_DIR}/rna_counts"

mkdir -p "${OUTPUT_DIR}" SLURM_outs SLURM_outs/array_outs

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

echo "=========================================="
echo "RNA counts submission"
echo "=========================================="
echo "Cells CSV:      ${CELLS_CSV}"
echo "Cells (BAMs):   ${N_BAMS}"
echo "GTF:            ${GTF}"
echo "Output prefix:  ${OUTPUT_PREFIX}"
if [ "${FILTER_SPLICED}" -eq 1 ]; then
    echo "Mode:           spliced-only (filter → count)"
else
    echo "Mode:           all reads (count directly on source BAMs)"
fi
echo "=========================================="

if [ "${FILTER_SPLICED}" -eq 1 ]; then
    SPLICED_DIR="${OUTPUT_DIR}/spliced_bams"
    BAM_LIST_SPLICED="${OUTPUT_DIR}/bam_list_spliced.txt"
    mkdir -p "${SPLICED_DIR}"

    # Pre-compute spliced-BAM paths (same basenames, in SPLICED_DIR).
    awk -v d="${SPLICED_DIR}" '{
        n = split($0, parts, "/")
        print d "/" parts[n]
    }' "${BAM_LIST_SRC}" > "${BAM_LIST_SPLICED}"

    filter_job=$(sbatch --parsable \
        --array="1-${N_BAMS}" \
        "${SCRIPTS_DIR}/filter_spliced_reads_array.sh" \
        "${BAM_LIST_SRC}" \
        "${SPLICED_DIR}")

    echo "Submitted spliced-filter array: ${filter_job} (1-${N_BAMS})"

    counts_job=$(sbatch --parsable \
        --dependency="afterok:${filter_job}" \
        "${SCRIPTS_DIR}/create_rna_count_matrix.sh" \
        "${BAM_LIST_SPLICED}" \
        "${GTF}" \
        "${OUTPUT_PREFIX}")

    echo "Submitted count-matrix job:     ${counts_job} (depends on ${filter_job})"
    monitor_jobs="${filter_job},${counts_job}"
else
    counts_job=$(sbatch --parsable \
        "${SCRIPTS_DIR}/create_rna_count_matrix.sh" \
        "${BAM_LIST_SRC}" \
        "${GTF}" \
        "${OUTPUT_PREFIX}")

    echo "Submitted count-matrix job:     ${counts_job}"
    monitor_jobs="${counts_job}"
fi

# --- Submit h5ad conversion with dependency ----------------------------------
h5ad_job=$(sbatch --parsable \
    --dependency="afterok:${counts_job}" \
    -J counts_to_h5ad \
    -o SLURM_outs/%x_%j.out \
    -p campus-new -c 1 --mem=8G -t 30:00 \
    --wrap="python3 '${SCRIPTS_DIR}/counts_matrix_to_h5ad.py' \
        '${OUTPUT_PREFIX}_matrix.csv' \
        '${OUTPUT_PREFIX}.h5ad' \
        --summary '${OUTPUT_PREFIX}_summary.csv'")

echo "Submitted h5ad conversion:      ${h5ad_job} (depends on ${counts_job})"
echo ""
echo "Monitor with:  squeue -j ${monitor_jobs},${h5ad_job}"
