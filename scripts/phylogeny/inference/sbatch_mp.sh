#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=SLURM_outs/inference/%x_%A_%a.out
#SBATCH --error=SLURM_outs/inference/%x_%A_%a.err
#SBATCH --array=0-15
#
# Camin-Sokal MP (phangorn Fitch on presence-only encoding) + 100 bootstrap.
# Runs on 'real' matrix only.

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${PANEL_TAG:?PANEL_TAG must be set}"
PANEL_ROOT="${PANEL_ROOT:-results/worm6_final/DNA_analysis/phylogeny/panels}"
OUT_ROOT="${OUT_ROOT:-results/worm6_final/DNA_analysis/phylogeny/inference}"
BOOTSTRAP="${BOOTSTRAP:-100}"

WORMS="${WORMS:-worm01 worm03 worm05 worm06 worm07 worm08 worm09 worm10 worm11 worm12 worm13 worm15 worm16 worm17 worm18 worm20}"
readarray -t WORM_ARR < <(printf '%s\n' $WORMS)
N=${#WORM_ARR[@]}
if (( SLURM_ARRAY_TASK_ID >= N )); then exit 0; fi
WORM="${WORM_ARR[$SLURM_ARRAY_TASK_ID]}"

PANEL_H5AD="${PANEL_ROOT}/${PANEL_TAG}/worm_${WORM}.h5ad"
OUT_DIR="${OUT_ROOT}/${PANEL_TAG}/mp/worm_${WORM}"
[[ ! -f "$PANEL_H5AD" ]] && { echo "no panel"; exit 0; }
[[ -f "${OUT_DIR}/real.majority.newick" ]] && { echo "skip"; exit 0; }

# Encode matrix.
eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy
python scripts/phylogeny/inference/encode_mp_matrix.py \
    --panel-h5ad "$PANEL_H5AD" --out-dir "$OUT_DIR"

# Run phangorn.
module load fhR
Rscript scripts/phylogeny/inference/mp_run.R \
    --tsv "${OUT_DIR}/real.tsv" \
    --out-dir "$OUT_DIR" \
    --tag real \
    --bootstrap "$BOOTSTRAP" \
    --k 10 --seed 0 --threads "$SLURM_CPUS_PER_TASK"
