#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=1:00:00
#SBATCH --output=SLURM_outs/inference/%x_%j.out
#SBATCH --error=SLURM_outs/inference/%x_%j.err
#
# Build one panel tag from a config.
#
# Usage:
#   sbatch --job-name=panel_tct_kept \
#       --export=ALL,CONFIG=scripts/phylogeny/configs/tct_kept.yaml,PANEL_TAG=tct_kept \
#       scripts/phylogeny/inference/sbatch_build_panel.sh

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${CONFIG:?CONFIG must be set}"
: "${PANEL_TAG:?PANEL_TAG must be set}"

echo "[$(date -Iseconds)] building panel tag=${PANEL_TAG} from ${CONFIG}"

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

python scripts/phylogeny/inference/build_panel_at_tag.py \
    --config     "$CONFIG" \
    --panel-tag  "$PANEL_TAG"

echo "[$(date -Iseconds)] done"
