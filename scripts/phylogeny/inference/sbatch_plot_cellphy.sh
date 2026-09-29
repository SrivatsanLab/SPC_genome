#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=0:30:00
#SBATCH --output=SLURM_outs/inference/%x_%j.out
#SBATCH --error=SLURM_outs/inference/%x_%j.err
#
# Render per-worm CellPhy phylograms, one PNG per (worm x matrix) + an overview grid.
#
# Usage:
#   sbatch --job-name=plot_<tag> \
#       --export=ALL,PANEL_TAG=<tag>[,MATRICES="real covmask shuffled"] \
#       scripts/phylogeny/inference/sbatch_plot_cellphy.sh

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${PANEL_TAG:?PANEL_TAG must be set}"
OUT_ROOT="${OUT_ROOT:-results/worm6_final/DNA_analysis/phylogeny/inference}"
MATRICES="${MATRICES:-real covmask shuffled}"

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

for M in $MATRICES; do
    echo "[$(date -Iseconds)] rendering matrix=${M}"
    python scripts/phylogeny/inference/plot_cellphy_trees.py \
        --inference-dir "${OUT_ROOT}/${PANEL_TAG}" \
        --matrix        "$M" \
        --out-dir       "${OUT_ROOT}/${PANEL_TAG}/figures/trees"
done

echo "[$(date -Iseconds)] done"
