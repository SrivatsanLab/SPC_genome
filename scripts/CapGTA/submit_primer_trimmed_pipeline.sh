#!/bin/bash
###########################################################################################################################
# End-to-end submission: primer-trimmed per-cell FASTQs → aligned BAMs → featureCounts → UMI counts h5ad.
#
# Three stages chained with afterok:
#   1. align_primer_trimmed_array.sh   SLURM array (N tasks each handling BATCH_SIZE cells)
#   2. featurecounts_primer_trimmed.sh single job on the full BAM list
#   3. build_umi_matrix.sh             single job: UMI dedup + h5ad write
#
# Usage:
#   scripts/CapGTA/submit_primer_trimmed_pipeline.sh \
#       <cells_list.txt> <input_root> <output_root> [batch_size=20] [run_min=15] \
#       [input_suffix=primer_trimmed] [allow_no_umi=0]
#
#   cells_list.txt  — one cell barcode per line (used by alignment array)
#   input_root      — <input_root>/<UDI>/<cell>_R{1,2}.<input_suffix>.fastq.gz
#   output_root     — new directory for pipeline outputs
#   batch_size      — cells per SLURM array task (default 20)
#   run_min         — primer-trim --run-min used upstream (stored in h5ad uns, default 15)
#   input_suffix    — file-name stem (default "primer_trimmed"; use "polyT_trimmed" for the UMI-optional variant)
#   allow_no_umi    — 0 = skip reads without UB:Z: UMI (strict); 1 = include them with coordinate dedup
#
# Output layout under <output_root>:
#   bams/<cell>.bam[.bai]
#   star_logs/<cell>_Log.final.out
#   bam_list.txt
#   counts_raw.txt + counts_raw.txt.summary
#   featurecounts_annotated/<cell>.bam.featureCounts.bam
#   umi_matrix.h5ad                      # cells × genes UMI counts
###########################################################################################################################

set -euo pipefail

if [ "$#" -lt 3 ] || [ "$#" -gt 7 ]; then
    echo "Usage: $0 <cells_list.txt> <input_root> <output_root> [batch_size=20] [run_min=15] [input_suffix=primer_trimmed] [allow_no_umi=0]" >&2
    exit 1
fi

CELLS_LIST="$1"
INPUT_ROOT="$2"
OUTPUT_ROOT="$3"
BATCH_SIZE="${4:-20}"
RUN_MIN="${5:-15}"
INPUT_SUFFIX="${6:-primer_trimmed}"
ALLOW_NO_UMI="${7:-0}"

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for f in "${CELLS_LIST}"; do
    if [ ! -f "${f}" ]; then
        echo "Error: required input not found: ${f}" >&2
        exit 1
    fi
done
if [ ! -d "${INPUT_ROOT}" ]; then
    echo "Error: input_root not found: ${INPUT_ROOT}" >&2
    exit 1
fi

GTF="/shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"
if [ ! -f "${GTF}" ]; then
    echo "Error: GTF not found: ${GTF}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_ROOT}" SLURM_outs SLURM_outs/array_outs

N_CELLS=$(wc -l < "${CELLS_LIST}")
N_TASKS=$(( (N_CELLS + BATCH_SIZE - 1) / BATCH_SIZE ))
BAMS_DIR="${OUTPUT_ROOT}/bams"
COUNTS_RAW="${OUTPUT_ROOT}/counts_raw.txt"
ANNO_DIR="${OUTPUT_ROOT}/featurecounts_annotated"
H5AD="${OUTPUT_ROOT}/umi_matrix.h5ad"

echo "=========================================="
echo "Primer-trimmed pipeline"
echo "=========================================="
echo "cells_list:    ${CELLS_LIST}  (${N_CELLS} cells)"
echo "input_root:    ${INPUT_ROOT}"
echo "output_root:   ${OUTPUT_ROOT}"
echo "batch_size:    ${BATCH_SIZE}  (${N_TASKS} array tasks)"
echo "run_min:       ${RUN_MIN}"
echo "input_suffix:  ${INPUT_SUFFIX}"
echo "allow_no_umi:  ${ALLOW_NO_UMI}"
echo "GTF:           ${GTF}"
echo "=========================================="

# --- 1) Per-cell alignment array ----------------------------------------------
align_job=$(sbatch --parsable \
    --array="1-${N_TASKS}%50" \
    "${SCRIPTS_DIR}/align_primer_trimmed_array.sh" \
    "${CELLS_LIST}" "${INPUT_ROOT}" "${OUTPUT_ROOT}" "${BATCH_SIZE}" "${INPUT_SUFFIX}")
echo "Submitted alignment array:    ${align_job} (1-${N_TASKS}%50)"

# --- 2) featureCounts on the full BAM directory --------------------------------
fc_job=$(sbatch --parsable \
    --dependency="afterok:${align_job}" \
    "${SCRIPTS_DIR}/featurecounts_primer_trimmed.sh" \
    "${BAMS_DIR}" "${GTF}" "${OUTPUT_ROOT}")
echo "Submitted featureCounts:      ${fc_job} (depends on ${align_job})"

# --- 3) UMI dedup + h5ad ------------------------------------------------------
h5ad_job=$(sbatch --parsable \
    --dependency="afterok:${fc_job}" \
    "${SCRIPTS_DIR}/build_umi_matrix.sh" \
    "${ANNO_DIR}" "${COUNTS_RAW}" "${H5AD}" "${RUN_MIN}" "${ALLOW_NO_UMI}")
echo "Submitted UMI matrix build:   ${h5ad_job} (depends on ${fc_job})"

echo ""
echo "Final UMI counts h5ad will land at:"
echo "  ${H5AD}"
echo ""
echo "Monitor with: squeue -j ${align_job},${fc_job},${h5ad_job}"
