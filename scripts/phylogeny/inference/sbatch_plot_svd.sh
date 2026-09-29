#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=0:30:00
#SBATCH --output=SLURM_outs/inference/%x_%j.out
#SBATCH --error=SLURM_outs/inference/%x_%j.err

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${PANEL_TAG:?PANEL_TAG must be set}"
OUT_ROOT="${OUT_ROOT:-results/worm6_final/DNA_analysis/phylogeny/inference}"

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

PANEL_ROOT="${PANEL_ROOT:-results/worm6_final/DNA_analysis/phylogeny/panels}"
SVD_SUBDIR="${SVD_SUBDIR:-svd}"

python scripts/phylogeny/inference/plot_svd_trees.py \
    --inference-dir "${OUT_ROOT}/${PANEL_TAG}" \
    --panels-dir    "${PANEL_ROOT}/${PANEL_TAG}" \
    --svd-subdir    "$SVD_SUBDIR" \
    --out-dir       "${OUT_ROOT}/${PANEL_TAG}/figures/svd"
