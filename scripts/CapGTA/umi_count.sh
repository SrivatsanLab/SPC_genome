#!/bin/bash
#SBATCH -J umi_count
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -p campus-new
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 2:00:00

###########################################################################################################################
# Walk featureCounts-annotated BAMs, dedup by (cell, gene, UMI), build h5ad.
#
# Usage:
#   sbatch scripts/CapGTA/umi_count.sh <anno_dir> <counts_raw.txt> <filter_stats.tsv|NA> <run_min> <out.h5ad>
###########################################################################################################################

set -euo pipefail

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <anno_dir> <counts_raw.txt> <filter_stats.tsv|NA> <run_min> <out.h5ad>" >&2
    exit 1
fi

ANNO_DIR="$1"
COUNTS_RAW="$2"
FILTER_STATS="$3"
RUN_MIN="$4"
OUT_H5AD="$5"

REPO="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome"

eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

ARGS=(--anno-dir "${ANNO_DIR}" --counts-raw "${COUNTS_RAW}"
      --run-min "${RUN_MIN}" --out-h5ad "${OUT_H5AD}")
if [ "${FILTER_STATS}" != "NA" ] && [ -f "${FILTER_STATS}" ]; then
    ARGS+=(--filter-stats "${FILTER_STATS}")
fi

python "${REPO}/scripts/CapGTA/umi_count.py" "${ARGS[@]}"
