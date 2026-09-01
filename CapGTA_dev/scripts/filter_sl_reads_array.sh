#!/bin/bash
#SBATCH -J filter_sl
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -p campus-new
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 1:00:00

###########################################################################################################################
# One SLURM array task per BAM: filter to reads containing an SL1 or SL2 3'-end seed.
#
# Usage:
#   sbatch --array=1-$(wc -l < bam_list.txt) \
#       CapGTA_dev/scripts/filter_sl_reads_array.sh <bam_list.txt> <output_dir>
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: sbatch --array=1-N $0 <bam_list.txt> <output_dir>" >&2
    exit 1
fi

BAM_LIST="$1"
OUTPUT_DIR="$2"

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID not set — submit as an array." >&2
    exit 1
fi

SRC_BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${BAM_LIST}")
if [ -z "${SRC_BAM}" ] || [ ! -f "${SRC_BAM}" ]; then
    echo "Error: source BAM at line ${SLURM_ARRAY_TASK_ID} missing: '${SRC_BAM}'" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
OUT_BAM="${OUTPUT_DIR}/$(basename "${SRC_BAM}")"

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
python3 "${REPO_ROOT}/CapGTA_dev/scripts/filter_sl_reads.py" "${SRC_BAM}" "${OUT_BAM}"
