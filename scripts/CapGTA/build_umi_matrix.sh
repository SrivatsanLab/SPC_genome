#!/bin/bash
#SBATCH -J build_umi_matrix
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -p campus-new
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 2:00:00

set -euo pipefail

if [ "$#" -lt 3 ] || [ "$#" -gt 6 ]; then
    echo "Usage: $0 <anno_dir> <counts_raw.txt> <out_h5ad> [run_min=15] [allow_no_umi=0] [resolve_ambiguous=0]" >&2
    exit 1
fi

ANNO_DIR="$1"
COUNTS_RAW="$2"
OUT_H5AD="$3"
RUN_MIN="${4:-15}"
ALLOW_NO_UMI="${5:-0}"
RESOLVE_AMBIG="${6:-0}"

# SLURM copies the submitted script to /var/tmp/slurmd/..., so BASH_SOURCE
# resolves there, not at scripts/CapGTA. Hard-code the real location.
SCRIPTS_DIR="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/scripts/CapGTA"

eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

ALLOW_FLAG=""
if [ "${ALLOW_NO_UMI}" = "1" ]; then ALLOW_FLAG="--allow-no-umi"; fi
RESOLVE_FLAG=""
if [ "${RESOLVE_AMBIG}" = "1" ]; then RESOLVE_FLAG="--resolve-ambiguous"; fi

python "${SCRIPTS_DIR}/build_umi_matrix.py" \
    --anno-dir   "${ANNO_DIR}" \
    --counts-raw "${COUNTS_RAW}" \
    --out-h5ad   "${OUT_H5AD}" \
    --run-min    "${RUN_MIN}" \
    ${ALLOW_FLAG} ${RESOLVE_FLAG}
