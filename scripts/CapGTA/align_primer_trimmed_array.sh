#!/bin/bash
#SBATCH -J align_primer_trim
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -p campus-new
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 2:00:00

###########################################################################################################################
# Paired-end STAR alignment for primer-trimmed per-cell FASTQs — SLURM array, one task per BATCH of cells.
#
# Each STAR invocation uses --genomeLoad NoSharedMemory (default) so there's no cross-task contention
# when SLURM puts multiple array tasks on the same node. For a 100 Mb C. elegans index, STAR load is
# ~15-30s, so 20 cells per task costs ~5-10 min of load overhead per task.
#
# For each cell in the batch:
#   1. Encode UB:Z:<umi> from the FASTQ header into the read name as _UMI_<umi>  (keeps UMI through STAR)
#   2. STAR paired-end align with RELAXED thresholds (match/score min 0.3; unique mapping only)
#   3. Index the resulting BAM
#
# Usage:
#   sbatch --array=1-N scripts/CapGTA/align_primer_trimmed_array.sh \
#       <cells_list.txt>  <input_root>  <output_root>  [batch_size=20] [input_suffix=primer_trimmed]
###########################################################################################################################

set -euo pipefail

if [ "$#" -lt 3 ] || [ "$#" -gt 5 ]; then
    echo "Usage: sbatch --array=1-N $0 <cells_list.txt> <input_root> <output_root> [batch_size=20] [input_suffix=primer_trimmed]" >&2
    exit 1
fi

CELLS_LIST="$1"
INPUT_ROOT="$2"
OUTPUT_ROOT="$3"
BATCH_SIZE="${4:-20}"
INPUT_SUFFIX="${5:-primer_trimmed}"

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID not set — submit as an array." >&2
    exit 1
fi

STAR_INDEX="/shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/STAR2Index"
# SLURM copies the submitted script to /var/tmp/slurmd/job<N>/, so BASH_SOURCE
# points there, not at scripts/CapGTA. Hard-code the real location.
SCRIPTS_DIR="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/scripts/CapGTA"
BAMS_DIR="${OUTPUT_ROOT}/bams"
LOGS_DIR="${OUTPUT_ROOT}/star_logs"
mkdir -p "${BAMS_DIR}" "${LOGS_DIR}"

SCRATCH="${TMPDIR:-/tmp}/align_primer_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${SCRATCH}"
trap 'rm -rf "${SCRATCH}"' EXIT

START=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH_SIZE + 1 ))
END=$(( SLURM_ARRAY_TASK_ID * BATCH_SIZE ))
N_TOTAL=$(wc -l < "${CELLS_LIST}")
if [ "${END}" -gt "${N_TOTAL}" ]; then END="${N_TOTAL}"; fi

echo "=========================================="
echo "Primer-trimmed alignment array task"
echo "  task_id       = ${SLURM_ARRAY_TASK_ID}"
echo "  batch         = ${START}..${END}  (${BATCH_SIZE} per task)"
echo "  cells_list    = ${CELLS_LIST}"
echo "  input_root    = ${INPUT_ROOT}"
echo "  output_root   = ${OUTPUT_ROOT}"
echo "  input_suffix  = ${INPUT_SUFFIX}"
echo "  scratch       = ${SCRATCH}"
echo "  started       = $(date -Iseconds)"
echo "=========================================="

module load STAR/2.7.11b-GCC-13.3.0 SAMtools/1.21-GCC-13.3.0
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

echo "which STAR:     $(which STAR)"
echo "which samtools: $(which samtools)"
echo "which python:   $(which python)"
echo ""

for i in $(seq "${START}" "${END}"); do
    CELL=$(sed -n "${i}p" "${CELLS_LIST}" | awk '{print $1}')
    if [ -z "${CELL}" ]; then continue; fi
    UDI="${CELL#*_*_*_*_}"
    R1_IN="${INPUT_ROOT}/${UDI}/${CELL}_R1.${INPUT_SUFFIX}.fastq.gz"
    R2_IN="${INPUT_ROOT}/${UDI}/${CELL}_R2.${INPUT_SUFFIX}.fastq.gz"
    if [ ! -f "${R1_IN}" ] || [ ! -f "${R2_IN}" ]; then
        echo "  [skip] ${CELL}: inputs not found (${R1_IN})"
        continue
    fi

    # Encode UMI into read name
    R1_ENC="${SCRATCH}/${CELL}_R1.enc.fastq.gz"
    R2_ENC="${SCRATCH}/${CELL}_R2.enc.fastq.gz"
    python "${SCRIPTS_DIR}/encode_umi_in_readname.py" "${R1_IN}" "${R1_ENC}"
    python "${SCRIPTS_DIR}/encode_umi_in_readname.py" "${R2_IN}" "${R2_ENC}"

    OUT_BAM="${BAMS_DIR}/${CELL}.bam"
    PFX="${SCRATCH}/${CELL}_"

    # STAR with NoSharedMemory (default). Stderr goes to the SLURM output so
    # errors are visible immediately.
    STAR --runThreadN 4 \
         --genomeDir "${STAR_INDEX}" \
         --readFilesIn "${R1_ENC}" "${R2_ENC}" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${PFX}" \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes NH HI AS nM NM MD \
         --outSAMunmapped Within \
         --outSAMattrRGline ID:"${CELL}" SM:"${CELL}" PL:ILLUMINA LB:"${CELL}" \
         --outFilterMultimapNmax 1 \
         --outFilterMatchNminOverLread 0.3 \
         --outFilterScoreMinOverLread 0.3 \
         --limitBAMsortRAM 2000000000

    mv "${PFX}Aligned.sortedByCoord.out.bam" "${OUT_BAM}"
    samtools index "${OUT_BAM}"
    cp "${PFX}Log.final.out" "${LOGS_DIR}/${CELL}_Log.final.out"
    # Clean up STAR's intermediate files to save scratch
    rm -f "${R1_ENC}" "${R2_ENC}" "${PFX}"*.out "${PFX}"*.tab "${PFX}"*.progress "${PFX}"SJ.out.tab 2>/dev/null || true
    rm -rf "${PFX}"_STARtmp 2>/dev/null || true

    echo "  [${i}/${N_TOTAL}] ${CELL} -> $(basename ${OUT_BAM}) at $(date -Iseconds)"
done

echo ""
echo "Task ${SLURM_ARRAY_TASK_ID} done at $(date -Iseconds)"
