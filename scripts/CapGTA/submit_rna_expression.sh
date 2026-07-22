#!/bin/bash
###########################################################################################################################
# Submit RNA expression estimation for a set of single-cell BAMs (simple-excess model).
#
# Prerequisite: an all-reads featureCounts run has already been done for the same cells
# (rna_counts_raw.txt from `submit_rna_counts.sh --all-reads`). The R_obs matrix and per-gene
# Length column are read from that file rather than recomputed.
#
# Three stages:
#   1. make_intergenic_bed.py           — one-shot (if intergenic.bed doesn't already exist).
#   2. count_intergenic_reads_array.sh  — SLURM array, one task per cell.
#   3. assemble_expression_matrix.py    — depends afterok on the array; writes CSV + h5ad.
#
# Usage:
#   scripts/CapGTA/submit_rna_expression.sh \
#       <real_cells.csv> \
#       <reference.fa> \
#       <annotation.gtf> \
#       <featurecounts_raw.txt> \
#       <output_dir> \
#       [--flank BP]  [--min-size BP]
#
# Layout under <output_dir>:
#   bam_list.txt
#   intergenic.bed
#   per_cell/{barcode}.tsv
#   rna_expression_matrix.csv
#   rna_expression.h5ad
###########################################################################################################################

set -euo pipefail

FLANK=5000
MIN_SIZE=10000

POSITIONAL=()
while [ "$#" -gt 0 ]; do
    case "$1" in
        --flank)    FLANK="$2";    shift 2 ;;
        --min-size) MIN_SIZE="$2"; shift 2 ;;
        *)          POSITIONAL+=("$1"); shift ;;
    esac
done
set -- "${POSITIONAL[@]}"

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 [--flank BP] [--min-size BP] <real_cells.csv> <reference.fa> <annotation.gtf> <featurecounts_raw.txt> <output_dir>" >&2
    exit 1
fi

CELLS_CSV="$1"
REFERENCE_FA="$2"
GTF="$3"
FC_RAW="$4"
OUTPUT_DIR="$5"

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FAI="${REFERENCE_FA}.fai"

for f in "${CELLS_CSV}" "${REFERENCE_FA}" "${FAI}" "${GTF}" "${FC_RAW}"; do
    if [ ! -f "${f}" ]; then
        echo "Error: required input not found: ${f}" >&2
        exit 1
    fi
done

INTERGENIC_BED="${OUTPUT_DIR}/intergenic.bed"
BAM_LIST="${OUTPUT_DIR}/bam_list.txt"
PER_CELL_DIR="${OUTPUT_DIR}/per_cell"
OUTPUT_PREFIX="${OUTPUT_DIR}/rna_expression"

mkdir -p "${OUTPUT_DIR}" "${PER_CELL_DIR}" SLURM_outs SLURM_outs/array_outs

# --- Build bam_list from real_cells.csv --------------------------------------
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

for bam in "$(head -1 "${BAM_LIST}")" "$(tail -1 "${BAM_LIST}")"; do
    if [ ! -f "${bam}" ]; then
        echo "Error: BAM not found (spot-check): ${bam}" >&2
        exit 1
    fi
done

# --- Generate intergenic BED (or reuse) --------------------------------------
if [ ! -s "${INTERGENIC_BED}" ]; then
    echo "Generating intergenic BED (flank=${FLANK} min_size=${MIN_SIZE})..."
    python3 "${SCRIPTS_DIR}/make_intergenic_bed.py" \
        "${GTF}" "${FAI}" "${INTERGENIC_BED}" \
        --flank "${FLANK}" --min-size "${MIN_SIZE}"
else
    echo "Reusing existing intergenic BED: ${INTERGENIC_BED}"
fi

echo "=========================================="
echo "RNA expression submission"
echo "=========================================="
echo "Cells CSV:      ${CELLS_CSV}"
echo "Cells (BAMs):   ${N_BAMS}"
echo "Reference:      ${REFERENCE_FA}"
echo "GTF:            ${GTF}"
echo "featureCounts:  ${FC_RAW}"
echo "Intergenic BED: ${INTERGENIC_BED}"
echo "Output prefix:  ${OUTPUT_PREFIX}"
echo "=========================================="

# --- Submit per-cell intergenic counting array -------------------------------
array_job=$(sbatch --parsable \
    --array="1-${N_BAMS}" \
    "${SCRIPTS_DIR}/count_intergenic_reads_array.sh" \
    "${BAM_LIST}" \
    "${INTERGENIC_BED}" \
    "${PER_CELL_DIR}")

echo "Submitted intergenic array:     ${array_job} (1-${N_BAMS})"

# --- Submit assembler with dependency ----------------------------------------
assemble_job=$(sbatch --parsable \
    --dependency="afterok:${array_job}" \
    -J assemble_expr \
    -o SLURM_outs/%x_%j.out \
    -p campus-new -c 2 --mem=16G -t 1:00:00 \
    --wrap="python3 '${SCRIPTS_DIR}/assemble_expression_matrix.py' \
        '${FC_RAW}' \
        '${INTERGENIC_BED}' \
        '${PER_CELL_DIR}' \
        '${OUTPUT_PREFIX}'")

echo "Submitted assembler:            ${assemble_job} (depends on ${array_job})"
echo ""
echo "Monitor with:  squeue -j ${array_job},${assemble_job}"
