#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --output=SLURM_outs/inference/%x_%A_%a.out
#SBATCH --error=SLURM_outs/inference/%x_%A_%a.err
#SBATCH --array=0-15
#
# Tree-guided haplotype imputation per worm.
#
# Usage:
#   sbatch --job-name=impute_<src_tag> \
#       --export=ALL,SRC_TAG=<src>,IMPUTED_TAG=<src>__svdImp,TREE_SUBDIR=svd_mquad5 \
#       scripts/phylogeny/inference/sbatch_impute.sh

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${SRC_TAG:?SRC_TAG must be set}"
: "${IMPUTED_TAG:?IMPUTED_TAG must be set}"
TREE_SUBDIR="${TREE_SUBDIR:-svd_mquad5}"
PSEUDO_ALT="${PSEUDO_ALT:-1}"
PSEUDO_REF="${PSEUDO_REF:-0}"          # >0 turns on symmetric non-carrier imputation
TRUST_MAX_DEPTH="${TRUST_MAX_DEPTH:-}"

PANEL_ROOT="${PANEL_ROOT:-results/worm6_final/DNA_analysis/phylogeny/panels}"
OUT_ROOT="${OUT_ROOT:-results/worm6_final/DNA_analysis/phylogeny/inference}"

WORMS="${WORMS:-worm01 worm03 worm05 worm06 worm07 worm08 worm09 worm10 worm11 worm12 worm13 worm15 worm16 worm17 worm18 worm20}"
readarray -t WORM_ARR < <(printf '%s\n' $WORMS)
N=${#WORM_ARR[@]}
if (( SLURM_ARRAY_TASK_ID >= N )); then exit 0; fi
WORM="${WORM_ARR[$SLURM_ARRAY_TASK_ID]}"

SRC_H5="${PANEL_ROOT}/${SRC_TAG}/worm_${WORM}.h5ad"
TREE_JSON="${OUT_ROOT}/${SRC_TAG}/${TREE_SUBDIR}/worm_${WORM}/tree.json"
OUT_H5="${PANEL_ROOT}/${IMPUTED_TAG}/worm_${WORM}.h5ad"

[[ ! -f "$SRC_H5" ]] && { echo "no panel h5ad"; exit 0; }
[[ ! -f "$TREE_JSON" ]] && { echo "no tree.json"; exit 0; }
[[ -f "$OUT_H5" ]] && { echo "skip: exists"; exit 0; }

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

TRUST_ARG=()
[[ -n "$TRUST_MAX_DEPTH" ]] && TRUST_ARG=(--trust-max-depth "$TRUST_MAX_DEPTH")

python scripts/phylogeny/inference/impute_panel.py \
    --panel-h5ad "$SRC_H5" \
    --tree-json  "$TREE_JSON" \
    --out-h5ad   "$OUT_H5" \
    --pseudo-alt "$PSEUDO_ALT" \
    --pseudo-ref "$PSEUDO_REF" \
    "${TRUST_ARG[@]}"
