#!/bin/bash
#SBATCH -J umi_filter
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -p campus-new
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 1:00:00

###########################################################################################################################
# Filter a per-cell BAM to primer-signature RNA reads and tag them with UMIs — one BAM per SLURM array task.
#
# Replaces filter_spliced_reads_array.sh in the UMI-based counting pipeline. Instead of selecting
# reads by spliced alignment (CIGAR has N), we select reads whose soft clip adjacent to the
# alignment carries the RT primer polyA/T tail + a recoverable 10 bp UMI, and tag passing reads
# with UB:Z:<umi>, XP:A:<side>, XR:i:<run_len>.
#
# Usage:
#   sbatch --array=1-$(wc -l < bam_list.txt) \
#       scripts/CapGTA/umi_tag_filter_array.sh <bam_list.txt> <output_dir> <run_min>
#
# Writes:
#   <output_dir>/<source_basename>.bam        # filtered + UB-tagged
#   <output_dir>/<source_basename>.bam.bai
#   <output_dir>/filter_stats.tsv             # appended per-task row
###########################################################################################################################

set -euo pipefail

if [ "$#" -lt 2 ] || [ "$#" -gt 4 ]; then
    echo "Usage: sbatch --array=1-N $0 <bam_list.txt> <output_dir> [run_min=6] [require_umi=1]" >&2
    echo "  require_umi: 1 -> only keep reads with recoverable UMI (default)" >&2
    echo "               0 -> keep all polyA/T-passing reads, no UB tag on UMI-less reads" >&2
    exit 1
fi

BAM_LIST="$1"
OUTPUT_DIR="$2"
RUN_MIN="${3:-6}"
REQUIRE_UMI="${4:-1}"

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
STATS_TSV="${OUTPUT_DIR}/filter_stats.tsv"

eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

# SLURM copies this script to a sandbox dir, so ${BASH_SOURCE[0]} doesn't locate
# sibling scripts. Resolve via an absolute repo-rooted path instead.
REPO="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome"
UMI_FLAG="--require-umi"
if [ "${REQUIRE_UMI}" = "0" ]; then
    UMI_FLAG="--no-require-umi"
fi
python "${REPO}/scripts/CapGTA/umi_tag_filter.py" \
    --src-bam "${SRC_BAM}" \
    --out-bam "${OUT_BAM}" \
    --run-min "${RUN_MIN}" \
    ${UMI_FLAG} \
    --stats-tsv "${STATS_TSV}"
