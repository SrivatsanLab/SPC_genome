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
METHODS="${METHODS:-cellphy distance mp}"

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

python scripts/phylogeny/inference/plot_all_methods.py \
    --inference-dir "${OUT_ROOT}/${PANEL_TAG}" \
    --out-dir       "${OUT_ROOT}/${PANEL_TAG}/figures/trees" \
    --methods       $METHODS
