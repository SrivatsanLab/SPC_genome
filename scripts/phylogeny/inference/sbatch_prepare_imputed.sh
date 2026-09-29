#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=SLURM_outs/inference/%x_%A_%a.out
#SBATCH --error=SLURM_outs/inference/%x_%A_%a.err
#SBATCH --array=0-15

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${PANEL_TAG:?PANEL_TAG must be set (imputed tag)}"
PANEL_ROOT="${PANEL_ROOT:-results/worm6_final/DNA_analysis/phylogeny/panels}"
OUT_ROOT="${OUT_ROOT:-results/worm6_final/DNA_analysis/phylogeny/inference}"
JOINT_VCF="${JOINT_VCF:-results/worm6_final/DNA_analysis/joint_variants/joint_variants.vcf.gz}"
SEED="${SEED:-0}"

WORMS="${WORMS:-worm01 worm03 worm05 worm06 worm07 worm08 worm09 worm10 worm11 worm12 worm13 worm15 worm16 worm17 worm18 worm20}"
readarray -t WORM_ARR < <(printf '%s\n' $WORMS)
N=${#WORM_ARR[@]}
if (( SLURM_ARRAY_TASK_ID >= N )); then exit 0; fi
WORM="${WORM_ARR[$SLURM_ARRAY_TASK_ID]}"

PANEL_H5="${PANEL_ROOT}/${PANEL_TAG}/worm_${WORM}.h5ad"
OUT_DIR="${OUT_ROOT}/${PANEL_TAG}/vcf/worm_${WORM}"

[[ ! -f "$PANEL_H5" ]] && { echo "no panel"; exit 0; }

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

python scripts/phylogeny/inference/prepare_vcf.py \
    --panel-h5ad "$PANEL_H5" \
    --joint-vcf  "$JOINT_VCF" \
    --out-dir    "$OUT_DIR" \
    --seed       "$SEED" \
    --imputed-panel-h5ad "$PANEL_H5" \
    --skip-existing
