#!/bin/bash
#SBATCH -J filter_spliced
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -p campus-new
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 1:00:00

###########################################################################################################################
# Filter a per-cell BAM to spliced reads (CIGAR contains 'N') — one BAM per SLURM array task.
#
# Used upstream of RNA read counting for STAR-aligned coassay BAMs, where each cell's BAM
# contains both DNA and RNA reads. Spliced reads (N in CIGAR = intron skip) are a conservative
# marker for RNA-origin reads.
#
# Usage:
#   sbatch --array=1-$(wc -l < bam_list.txt) \
#       scripts/CapGTA/filter_spliced_reads_array.sh <bam_list.txt> <output_dir>
#
# Writes:
#   <output_dir>/<source_basename>.bam        # spliced-only reads
#   <output_dir>/<source_basename>.bam.bai
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: sbatch --array=1-N $0 <bam_list.txt> <output_dir>" >&2
    exit 1
fi

BAM_LIST="$1"
OUTPUT_DIR="$2"

if [ ! -f "${BAM_LIST}" ]; then
    echo "Error: bam_list not found: ${BAM_LIST}" >&2
    exit 1
fi

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID not set — submit as an array." >&2
    exit 1
fi

SRC_BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${BAM_LIST}")
if [ -z "${SRC_BAM}" ] || [ ! -f "${SRC_BAM}" ]; then
    echo "Error: source BAM at line ${SLURM_ARRAY_TASK_ID} missing or unreadable: '${SRC_BAM}'" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
OUT_BAM="${OUTPUT_DIR}/$(basename "${SRC_BAM}")"

module load SAMtools

# CIGAR expression filter requires samtools >= 1.15.
# `cigar =~ "N"` matches reads whose CIGAR contains an N operator (intron skip).
samtools view \
    --threads 2 \
    -b \
    -e 'cigar =~ "N"' \
    -o "${OUT_BAM}" \
    "${SRC_BAM}"

samtools index "${OUT_BAM}"

N_IN=$(samtools view -c "${SRC_BAM}")
N_OUT=$(samtools view -c "${OUT_BAM}")
echo "Filtered spliced reads: ${SRC_BAM}"
echo "  in:  ${N_IN}"
echo "  out: ${N_OUT}"
echo "  -> ${OUT_BAM}"
