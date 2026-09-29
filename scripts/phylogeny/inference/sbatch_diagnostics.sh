#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=1:00:00
#SBATCH --output=SLURM_outs/inference/%x_%A_%a.out
#SBATCH --error=SLURM_outs/inference/%x_%A_%a.err
#SBATCH --array=0-15
#
# Cross-tree diagnostics (max-Jaccard + permutation null) per worm.
#
# Usage:
#   sbatch --job-name=diag_<tag> \
#       --export=ALL,PANEL_TAG=<tag> \
#       scripts/phylogeny/inference/sbatch_diagnostics.sh

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${PANEL_TAG:?PANEL_TAG must be set}"
PANEL_ROOT="${PANEL_ROOT:-results/worm6_final/DNA_analysis/phylogeny/panels}"
OUT_ROOT="${OUT_ROOT:-results/worm6_final/DNA_analysis/phylogeny/inference}"
N_PERM="${N_PERM:-500}"

WORMS="${WORMS:-worm01 worm03 worm05 worm06 worm07 worm08 worm09 worm10 worm11 worm12 worm13 worm15 worm16 worm17 worm18 worm20}"
readarray -t WORM_ARR < <(printf '%s\n' $WORMS)
N=${#WORM_ARR[@]}
if (( SLURM_ARRAY_TASK_ID >= N )); then echo "task oob"; exit 0; fi
WORM="${WORM_ARR[$SLURM_ARRAY_TASK_ID]}"

PANEL_H5AD="${PANEL_ROOT}/${PANEL_TAG}/worm_${WORM}.h5ad"
INF_DIR="${OUT_ROOT}/${PANEL_TAG}"
OUT_DIR="${INF_DIR}/summary/worm_${WORM}"

if [[ ! -f "${INF_DIR}/cellphy/worm_${WORM}/real/run.raxml.bestTree" ]]; then
    echo "  SKIP: no real bestTree for $WORM"; exit 0
fi

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

python scripts/phylogeny/inference/diagnostics.py \
    --worm "$WORM" \
    --inference-dir "$INF_DIR" \
    --panel-h5ad    "$PANEL_H5AD" \
    --out-dir       "$OUT_DIR" \
    --n-perm        "$N_PERM"
