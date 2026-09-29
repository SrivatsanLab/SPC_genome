#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --output=SLURM_outs/inference/%x_%j.out
#SBATCH --error=SLURM_outs/inference/%x_%j.err
#
# D1-D3 panel diagnostics per plan §1.
#
# Usage:
#   sbatch --job-name=paneldiag_<tag> \
#       --export=ALL,PANEL_TAG=<tag> \
#       scripts/phylogeny/inference/sbatch_panel_diagnostics.sh

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${PANEL_TAG:?PANEL_TAG must be set}"
PANEL_ROOT="${PANEL_ROOT:-results/worm6_final/DNA_analysis/phylogeny/panels}"

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

python scripts/phylogeny/inference/panel_diagnostics.py \
    --panels-dir "${PANEL_ROOT}/${PANEL_TAG}" \
    --out-dir    "${PANEL_ROOT}/${PANEL_TAG}/diagnostics"
