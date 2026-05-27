#!/bin/bash
##########################################################################################################################
# Driver for the single-cell-fastq demux test pipeline.
# Run from the project root on the login node. Submits the per-cell array and a follow-up compile job.
##########################################################################################################################

set -euo pipefail

show_help() {
    cat <<EOF
Usage: $0 -i <input_cells_dir> -s <sample_name> -g <star_index_dir> \\
          -r <reference_fasta> -a <gtf_file> [-x <sample_suffix>] [-n <max_cells>]

Required:
  -i <dir>        Directory containing per-cell paired FASTQs named
                  {cell}_{suffix}_R{1,2}.fastq.gz
  -s <name>       Output sample name (creates data/{name}/ and results/{name}/)
  -g <dir>        STAR index directory (e.g. .../STARIndex)
  -r <fasta>      Reference FASTA whose contig names match the STAR index. A
                  sibling .dict (Picard CreateSequenceDictionary) must exist.
  -a <gtf>        GTF annotation file (for exonic enrichment)

Optional:
  -x <suffix>     Suffix between cell ID and "_R1.fastq.gz" in input filenames.
                  Default: same as -s (i.e. assumes filenames look like
                  {cell}_{sample}_R1.fastq.gz).
  -n <int>        Limit the array to the first N cells (for smoke-testing).
EOF
}

SAMPLE_SUFFIX=""
MAX_CELLS=""

while getopts ":i:s:g:r:a:x:n:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        s) SAMPLE="$OPTARG" ;;
        g) STAR_INDEX_DIR="$OPTARG" ;;
        r) REFERENCE_FASTA="$OPTARG" ;;
        a) GTF_FILE="$OPTARG" ;;
        x) SAMPLE_SUFFIX="$OPTARG" ;;
        n) MAX_CELLS="$OPTARG" ;;
        h) show_help; exit 0 ;;
        \?) echo "Invalid option: -$OPTARG" >&2; show_help; exit 1 ;;
        :)  echo "Option -$OPTARG requires an argument." >&2; show_help; exit 1 ;;
    esac
done

if [ -z "${INPUT_DIR:-}" ] || [ -z "${SAMPLE:-}" ] || [ -z "${STAR_INDEX_DIR:-}" ] || \
   [ -z "${REFERENCE_FASTA:-}" ] || [ -z "${GTF_FILE:-}" ]; then
    echo "Missing required argument." >&2
    show_help
    exit 1
fi

if [ -z "${SAMPLE_SUFFIX}" ]; then
    SAMPLE_SUFFIX="${SAMPLE}"
fi

SCRIPTS_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"   # project root
BAM_DIR="${SCRIPTS_DIR}/data/${SAMPLE}/sc_outputs"
QC_DIR="${SCRIPTS_DIR}/data/${SAMPLE}/qc_metrics"
RESULTS_DIR="${SCRIPTS_DIR}/results/${SAMPLE}"
PER_CELL_CSV_DIR="${RESULTS_DIR}/per_cell_qc"
CELLS_FILE="${RESULTS_DIR}/cells.txt"
FINAL_CSV="${RESULTS_DIR}/sc_qc_metrics.csv"

mkdir -p "${BAM_DIR}" "${QC_DIR}" "${PER_CELL_CSV_DIR}" SLURM_outs/array_outs

echo "=========================================="
echo "single-cell demux test driver"
echo "  input dir:    ${INPUT_DIR}"
echo "  sample:       ${SAMPLE}"
echo "  suffix:       ${SAMPLE_SUFFIX}"
echo "  STAR index:   ${STAR_INDEX_DIR}"
echo "  ref FASTA:    ${REFERENCE_FASTA}"
echo "  gtf:          ${GTF_FILE}"
echo "  bam dir:      ${BAM_DIR}"
echo "  results dir:  ${RESULTS_DIR}"
echo "=========================================="

# Enumerate cells from R1 filenames
echo "Enumerating cells in ${INPUT_DIR}..."
( cd "${INPUT_DIR}" && \
  ls *_"${SAMPLE_SUFFIX}"_R1.fastq.gz 2>/dev/null | \
  sed "s/_${SAMPLE_SUFFIX}_R1\\.fastq\\.gz\$//" ) > "${CELLS_FILE}"

N_CELLS=$(wc -l < "${CELLS_FILE}")
if [ "${N_CELLS}" -eq 0 ]; then
    echo "ERROR: no cells found in ${INPUT_DIR} matching *_${SAMPLE_SUFFIX}_R1.fastq.gz" >&2
    exit 1
fi
echo "Found ${N_CELLS} cells → ${CELLS_FILE}"

if [ -n "${MAX_CELLS}" ]; then
    head -"${MAX_CELLS}" "${CELLS_FILE}" > "${CELLS_FILE}.tmp" && mv "${CELLS_FILE}.tmp" "${CELLS_FILE}"
    N_CELLS=$(wc -l < "${CELLS_FILE}")
    echo "Limited to first ${N_CELLS} cells via -n"
fi

# Submit per-cell array
ARRAY_JOB=$(sbatch --parsable --array=1-${N_CELLS} \
    "${SCRIPTS_DIR}/scripts/CapGTA/test_new_demux/sc_align_qc_array.sh" \
    "${CELLS_FILE}" \
    "${INPUT_DIR}" \
    "${SAMPLE_SUFFIX}" \
    "${BAM_DIR}" \
    "${QC_DIR}" \
    "${PER_CELL_CSV_DIR}" \
    "${STAR_INDEX_DIR}" \
    "${REFERENCE_FASTA}" \
    "${GTF_FILE}" \
    "${SCRIPTS_DIR}")
echo "Submitted per-cell array job: ${ARRAY_JOB} (1-${N_CELLS})"

# Submit compile job, dependent on array (afterany so partial successes still compile)
COMPILE_JOB=$(sbatch --parsable --dependency=afterany:${ARRAY_JOB} \
    "${SCRIPTS_DIR}/scripts/CapGTA/test_new_demux/compile_qc.sh" \
    "${PER_CELL_CSV_DIR}" \
    "${FINAL_CSV}")
echo "Submitted compile job: ${COMPILE_JOB} (will run after ${ARRAY_JOB})"

echo ""
echo "Outputs:"
echo "  BAMs:     ${BAM_DIR}/{cell}.bam"
echo "  per-cell: ${PER_CELL_CSV_DIR}/{cell}.csv"
echo "  final:    ${FINAL_CSV}"
