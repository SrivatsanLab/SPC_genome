#!/bin/bash
###########################################################################################################################
# Submit RNA expression estimation for a set of single-cell BAMs (simple-excess model,
# featureCounts-NoFeatures baseline).
#
# Prerequisite: an all-reads featureCounts run has already been done for the same cells
# (rna_counts_raw.txt from `submit_rna_counts.sh --all-reads`, alongside .summary).
# The R_obs matrix, per-gene Length, and per-cell NoFeatures fragment counts are all
# read from files that job produced — no per-cell SLURM work is required.
#
# One stage: assemble_expression_matrix.py (writes CSV + h5ad).
#
# Usage:
#   scripts/CapGTA/submit_rna_expression.sh \
#       <reference.fa> \
#       <featurecounts_raw.txt> \
#       <output_dir>
#
# Layout under <output_dir>:
#   rna_expression_matrix.csv
#   rna_expression.h5ad
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <reference.fa> <featurecounts_raw.txt> <output_dir>" >&2
    exit 1
fi

REFERENCE_FA="$1"
FC_RAW="$2"
OUTPUT_DIR="$3"

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FAI="${REFERENCE_FA}.fai"
FC_SUMMARY="${FC_RAW}.summary"

for f in "${REFERENCE_FA}" "${FAI}" "${FC_RAW}" "${FC_SUMMARY}"; do
    if [ ! -f "${f}" ]; then
        echo "Error: required input not found: ${f}" >&2
        exit 1
    fi
done

OUTPUT_PREFIX="${OUTPUT_DIR}/rna_expression"

mkdir -p "${OUTPUT_DIR}" SLURM_outs

echo "=========================================="
echo "RNA expression submission"
echo "=========================================="
echo "Reference:      ${REFERENCE_FA}"
echo "featureCounts:  ${FC_RAW}"
echo "Summary:        ${FC_SUMMARY}"
echo "Output prefix:  ${OUTPUT_PREFIX}"
echo "=========================================="

assemble_job=$(sbatch --parsable \
    -J assemble_expr \
    -o SLURM_outs/%x_%j.out \
    -p campus-new -c 2 --mem=16G -t 1:00:00 \
    --wrap="python3 '${SCRIPTS_DIR}/assemble_expression_matrix.py' \
        '${FC_RAW}' \
        '${FC_SUMMARY}' \
        '${FAI}' \
        '${OUTPUT_PREFIX}'")

echo "Submitted assembler: ${assemble_job}"
echo ""
echo "Monitor with:  squeue -j ${assemble_job}"
