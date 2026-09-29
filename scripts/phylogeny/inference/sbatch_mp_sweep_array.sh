#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1:30:00
#SBATCH --output=SLURM_outs/inference/%x_%A_%a.out
#SBATCH --error=SLURM_outs/inference/%x_%A_%a.err
#SBATCH --array=0-15
#
# 16-task array wrapper that runs the (K, gain, loss) parsimony sweep for
# every worm of one PANEL_TAG. Each task runs the sweep script for one worm,
# iterating over all configured (K, cost-matrix) tuples serially.
#
# Total wall per task ≈ 20-30 min (11 configs × ~2 min at bootstrap=100).
# Set --time=1:30:00 to give slack (fits within the maintenance window when
# submitted 24+ h before Sunday 6 am).
#
# Usage:
#   sbatch --job-name=mpsw_<tag> \
#       --export=ALL,PANEL_TAG=<tag>,[CONFIGS="a b c"],[BOOTSTRAP=100],[OUT_ROOT=scratch/mp_sweep] \
#       scripts/phylogeny/inference/sbatch_mp_sweep_array.sh

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${PANEL_TAG:?PANEL_TAG must be set}"
BOOTSTRAP="${BOOTSTRAP:-100}"
OUT_ROOT="${OUT_ROOT:-scratch/mp_sweep}"
CONFIGS="${CONFIGS:-}"

WORMS="${WORMS:-worm01 worm03 worm05 worm06 worm07 worm08 worm09 worm10 worm11 worm12 worm13 worm15 worm16 worm17 worm18 worm20}"
readarray -t WORM_ARR < <(printf '%s\n' $WORMS)
N=${#WORM_ARR[@]}
if (( SLURM_ARRAY_TASK_ID >= N )); then exit 0; fi
WORM="${WORM_ARR[$SLURM_ARRAY_TASK_ID]}"

PANEL_H5="results/worm6_final/DNA_analysis/phylogeny/panels/${PANEL_TAG}/worm_${WORM}.h5ad"
if [[ ! -f "$PANEL_H5" ]]; then
    echo "skip: no panel h5ad ${PANEL_H5}"; exit 0
fi

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
