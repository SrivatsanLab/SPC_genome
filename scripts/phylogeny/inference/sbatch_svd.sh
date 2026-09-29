#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=0:30:00
#SBATCH --output=SLURM_outs/inference/%x_%A_%a.out
#SBATCH --error=SLURM_outs/inference/%x_%A_%a.err
#SBATCH --array=0-15

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${PANEL_TAG:?PANEL_TAG must be set}"
PANEL_ROOT="${PANEL_ROOT:-results/worm6_final/DNA_analysis/phylogeny/panels}"
OUT_ROOT="${OUT_ROOT:-results/worm6_final/DNA_analysis/phylogeny/inference}"
NCOMP="${NCOMP:-3}"
CORE_SEL="${CORE_SEL:-none}"     # none | centroid | mquad
N_CORE="${N_CORE:-}"
MQUAD_BIC="${MQUAD_BIC:-5.0}"
CORE_TAG="${CORE_TAG:-}"          # optional subdir tag e.g. "centroid500"
MIN_CLADE_SIZE="${MIN_CLADE_SIZE:-10}"   # SVD stop rule: don't split below this many cells
MAX_DEPTH="${MAX_DEPTH:-6}"              # SVD stop rule: cap recursion depth

WORMS="${WORMS:-worm01 worm03 worm05 worm06 worm07 worm08 worm09 worm10 worm11 worm12 worm13 worm15 worm16 worm17 worm18 worm20}"
readarray -t WORM_ARR < <(printf '%s\n' $WORMS)
N=${#WORM_ARR[@]}
if (( SLURM_ARRAY_TASK_ID >= N )); then exit 0; fi
WORM="${WORM_ARR[$SLURM_ARRAY_TASK_ID]}"

PANEL_H5AD="${PANEL_ROOT}/${PANEL_TAG}/worm_${WORM}.h5ad"
# Route output into a core-selection-specific subdir when CORE_TAG is set.
if [[ -n "$CORE_TAG" ]]; then
    OUT_DIR="${OUT_ROOT}/${PANEL_TAG}/svd_${CORE_TAG}/worm_${WORM}"
else
    OUT_DIR="${OUT_ROOT}/${PANEL_TAG}/svd/worm_${WORM}"
fi
[[ ! -f "$PANEL_H5AD" ]] && { echo "no panel"; exit 0; }
[[ -f "${OUT_DIR}/tree.json" ]] && { echo "skip"; exit 0; }

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy
N_CORE_ARG=()
[[ -n "$N_CORE" ]] && N_CORE_ARG=(--n-core "$N_CORE")

python scripts/phylogeny/inference/run_svd.py \
    --panel-h5ad "$PANEL_H5AD" --out-dir "$OUT_DIR" --ncomp "$NCOMP" \
    --core-selection "$CORE_SEL" \
    "${N_CORE_ARG[@]}" \
    --mquad-delta-bic-threshold "$MQUAD_BIC" \
    --min-clade-size "$MIN_CLADE_SIZE" \
    --max-depth "$MAX_DEPTH" \
    --mquad-nproc "$SLURM_CPUS_PER_TASK"
