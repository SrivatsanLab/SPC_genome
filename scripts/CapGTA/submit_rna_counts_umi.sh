#!/bin/bash
###########################################################################################################################
# Submit the UMI-based RNA counting pipeline for a set of single-cell BAMs.
#
# Three-stage pipeline (compare with submit_rna_counts.sh, which uses the spliced-only + fragment
# count path):
#   1. umi_tag_filter_array.sh        SLURM array, one task per cell. For each read with a
#                                     primer-adjacent polyA/T >= run_min bp AND a recoverable 10 bp
#                                     UMI window, write a UB:Z-tagged copy to umi_bams/.
#   2. featureCounts_assign.sh        single job, 16 threads. featureCounts with -R BAM so each
#                                     passing read is annotated with XT:Z:<gene_id>.
#   3. umi_count.py                   single job. Walk per-cell annotated BAMs, dedup by
#                                     (barcode × gene × UMI), build sparse UMI matrix h5ad.
#
# Usage:
#   scripts/CapGTA/submit_rna_counts_umi.sh \
#       <real_cells.csv> \
#       <annotation.gtf> \
#       <output_dir> \
#       [run_min=6]
#
# Layout under <output_dir>:
#   bam_list_source.txt
#   umi_bams/<barcode>.bam[.bai]
#   bam_list_umi.txt
#   filter_stats.tsv
#   featurecounts_annotated/<barcode>.bam.featureCounts.bam
#   counts_raw.txt            # featureCounts raw matrix
#   counts_raw.txt.summary    # featureCounts QC
#   rna_umi_matrix.h5ad       # cells × genes UMI counts (+ read_counts layer, obs QC)
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
H5AD="${OUTPUT_DIR}/rna_umi_matrix.h5ad"

mkdir -p "${UMI_DIR}" "${ANNO_DIR}"

# --- Source bam_list from real_cells.csv --------------------------------------
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
for bam in "$(head -1 "${BAM_LIST_SRC}")" "$(tail -1 "${BAM_LIST_SRC}")"; do
    if [ ! -f "${bam}" ]; then
        echo "Error: source BAM not found (spot-check): ${bam}" >&2
        exit 1
    fi
done

# --- Pre-compute UMI-filter output paths --------------------------------------
awk -v d="${UMI_DIR}" '{
    n = split($0, parts, "/")
    print d "/" parts[n]
}' "${BAM_LIST_SRC}" > "${BAM_LIST_UMI}"

# Fresh stats file for this run
rm -f "${STATS_TSV}"

echo "=========================================="
echo "UMI-based RNA counts submission"
echo "=========================================="
echo "Cells CSV:      ${CELLS_CSV}"
echo "Cells (BAMs):   ${N_BAMS}"
echo "GTF:            ${GTF}"
echo "Output dir:     ${OUTPUT_DIR}"
echo "run_min (L):    ${RUN_MIN}"
echo "=========================================="

# --- 1) Per-cell UMI filter array --------------------------------------------
filter_job=$(sbatch --parsable \
    --array="1-${N_BAMS}" \
    "${SCRIPTS_DIR}/umi_tag_filter_array.sh" \
    "${BAM_LIST_SRC}" \
    "${UMI_DIR}" \
    "${RUN_MIN}")
echo "Submitted UMI filter array:     ${filter_job} (1-${N_BAMS})"

# --- 2) featureCounts gene assignment ----------------------------------------
fc_job=$(sbatch --parsable \
    --dependency="afterok:${filter_job}" \
    "${SCRIPTS_DIR}/featureCounts_assign.sh" \
    "${BAM_LIST_UMI}" \
    "${GTF}" \
    "${OUTPUT_DIR}")
echo "Submitted featureCounts job:    ${fc_job} (depends on ${filter_job})"

# --- 3) UMI dedup + h5ad build -----------------------------------------------
umi_count_job=$(sbatch --parsable \
    --dependency="afterok:${fc_job}" \
    "${SCRIPTS_DIR}/umi_count.sh" \
    "${ANNO_DIR}" "${COUNTS_RAW}" "${STATS_TSV}" "${RUN_MIN}" "${H5AD}")
echo "Submitted UMI count+h5ad job:   ${umi_count_job} (depends on ${fc_job})"

echo ""
echo "Monitor with:  squeue -j ${filter_job},${fc_job},${umi_count_job}"
