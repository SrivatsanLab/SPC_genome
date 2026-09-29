#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=1:00:00
#SBATCH --output=SLURM_outs/inference/%x_%j.out
#SBATCH --error=SLURM_outs/inference/%x_%j.err

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${CONFIG:?CONFIG must be set}"
: "${OUT_DIR:?OUT_DIR must be set}"

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

python scripts/phylogeny/inference/panel_filter_audit.py \
    --config "$CONFIG" \
    --out-dir "$OUT_DIR"
