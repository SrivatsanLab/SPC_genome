#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0:45:00
#SBATCH --output=SLURM_outs/inference/%x_%A_%a.out
#SBATCH --error=SLURM_outs/inference/%x_%A_%a.err
#SBATCH --array=0-15
#
# Top-down SVD sweep: recursive K=2 bipartitioning down to 2-cell leaves.
# 16-task array (one per worm). Each task runs the full 9-config grid
# serially: ncomp ∈ {1,2,3} × core_selection ∈ {none, mquad@5, mquad@10}.
#
# Output layout:
#   ${OUT_ROOT}/${PANEL_TAG}/<config>/worm_<W>/tree.json
#
# Usage:
#   sbatch --job-name=svdtd_<tag> \
#       --export=ALL,PANEL_TAG=all_variants_relax03 \
#       scripts/phylogeny/inference/sbatch_svd_topdown.sh

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${PANEL_TAG:?PANEL_TAG must be set}"
PANEL_ROOT="${PANEL_ROOT:-results/worm6_final/DNA_analysis/phylogeny/panels}"
OUT_ROOT="${OUT_ROOT:-results/worm6_final/DNA_analysis/phylogeny/topdown}"
MIN_CLADE_SIZE="${MIN_CLADE_SIZE:-2}"
MAX_DEPTH="${MAX_DEPTH:-20}"

WORMS="${WORMS:-worm01 worm03 worm05 worm06 worm07 worm08 worm09 worm10 worm11 worm12 worm13 worm15 worm16 worm17 worm18 worm20}"
readarray -t WORM_ARR < <(printf '%s\n' $WORMS)
N=${#WORM_ARR[@]}
if (( SLURM_ARRAY_TASK_ID >= N )); then exit 0; fi
WORM="${WORM_ARR[$SLURM_ARRAY_TASK_ID]}"

PANEL_H5="${PANEL_ROOT}/${PANEL_TAG}/worm_${WORM}.h5ad"
if [[ ! -f "$PANEL_H5" ]]; then
    echo "skip: no panel h5ad ${PANEL_H5}"; exit 0
fi

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

# Config grid: (name, ncomp, core_selection, mquad_bic)
# Original: none / mquad5 / mquad10 across ncomp={1,2,3}
# Hybrid variants: mquad above hybrid_switch_cells (default 20), set-op filter below.
CONFIGS=(
    "n1_none    1 none    5"
    "n2_none    2 none    5"
    "n3_none    3 none    5"
    "n1_mquad5  1 mquad   5"
    "n2_mquad5  2 mquad   5"
    "n3_mquad5  3 mquad   5"
    "n1_mquad10 1 mquad  10"
    "n2_mquad10 2 mquad  10"
    "n3_mquad10 3 mquad  10"
    "n1_hybrid  1 hybrid  5"
    "n2_hybrid  2 hybrid  5"
    "n3_hybrid  3 hybrid  5"
)

for row in "${CONFIGS[@]}"; do
    read -r NAME NCOMP CORE_SEL MQUAD_BIC <<< "$row"
    OUT_DIR="${OUT_ROOT}/${PANEL_TAG}/${NAME}/worm_${WORM}"
    if [[ -f "${OUT_DIR}/tree.json" ]]; then
        echo "[skip] ${NAME}/${WORM}: tree.json exists"
        continue
    fi
    echo "[run]  ${NAME}/${WORM}: ncomp=${NCOMP} core=${CORE_SEL} bic=${MQUAD_BIC}"
    python scripts/phylogeny/inference/run_svd.py \
        --panel-h5ad "$PANEL_H5" --out-dir "$OUT_DIR" \
        --ncomp "$NCOMP" \
        --core-selection "$CORE_SEL" \
        --mquad-delta-bic-threshold "$MQUAD_BIC" \
        --min-clade-size "$MIN_CLADE_SIZE" \
        --max-depth "$MAX_DEPTH" \
        --mquad-nproc "$SLURM_CPUS_PER_TASK"
done
