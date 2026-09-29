#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=SLURM_outs/inference/%x_%j.out
#SBATCH --error=SLURM_outs/inference/%x_%j.err
#
# One-worm MP encoding + model sweep. Runs mp_encoding_model_sweep.py which
# in turn shells out to the phangorn runner for each (K, k) config.
#
# Usage:
#   sbatch --job-name=mp_sweep_<worm> \
#       --export=ALL,PANEL_TAG=<tag>,WORM=<worm>,[BOOTSTRAP=100],[OUT_ROOT=scratch/mp_sweep] \
#       scripts/phylogeny/inference/sbatch_mp_encoding_model_sweep.sh

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${PANEL_TAG:?PANEL_TAG must be set}"
: "${WORM:?WORM must be set (e.g. worm10)}"
BOOTSTRAP="${BOOTSTRAP:-100}"
OUT_ROOT="${OUT_ROOT:-scratch/mp_sweep}"
# Optional: whitespace-separated list of config labels to run (default = all).
CONFIGS="${CONFIGS:-}"

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

CFG_ARGS=()
if [[ -n "$CONFIGS" ]]; then
    CFG_ARGS=(--configs $CONFIGS)
fi

python scripts/phylogeny/inference/mp_encoding_model_sweep.py \
    --panel-tag "$PANEL_TAG" \
    --worm "$WORM" \
    --out-root "$OUT_ROOT" \
    --bootstrap "$BOOTSTRAP" \
    --threads "$SLURM_CPUS_PER_TASK" \
    "${CFG_ARGS[@]}"
