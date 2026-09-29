#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=1:00:00
#SBATCH --output=SLURM_outs/inference/%x_%A_%a.out
#SBATCH --error=SLURM_outs/inference/%x_%A_%a.err
#SBATCH --array=0-15
#
# Prepare per-worm real/covmask/shuffled VCFs from a panel h5ad dir.
#
# Usage:
#   sbatch --job-name=prep_v0 \
#       --export=ALL,PANEL_TAG=v0_original \
#       scripts/phylogeny/inference/sbatch_prepare.sh
#
# Env:
#   PANEL_TAG        subdir name under phylogeny/panels/ (required)
#   PANEL_ROOT       phylogeny/panels root (default: repo layout)
#   OUT_ROOT         phylogeny/inference root (default: repo layout)
#   JOINT_VCF        joint VCF path (default: repo layout)
#   WORMS            space-separated worm IDs (default: 16-worm set)
#   SEED             shuffle seed (default 0)

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${PANEL_TAG:?PANEL_TAG must be set (subdir under phylogeny/panels/)}"
PANEL_ROOT="${PANEL_ROOT:-results/worm6_final/DNA_analysis/phylogeny/panels}"
OUT_ROOT="${OUT_ROOT:-results/worm6_final/DNA_analysis/phylogeny/inference}"
JOINT_VCF="${JOINT_VCF:-results/worm6_final/DNA_analysis/joint_variants/joint_variants.vcf.gz}"
SEED="${SEED:-0}"

WORMS="${WORMS:-worm01 worm03 worm05 worm06 worm07 worm08 worm09 worm10 worm11 worm12 worm13 worm15 worm16 worm17 worm18 worm20}"
readarray -t WORM_ARR < <(printf '%s\n' $WORMS)
N=${#WORM_ARR[@]}

if (( SLURM_ARRAY_TASK_ID >= N )); then
    echo "task id $SLURM_ARRAY_TASK_ID out of range (0..$((N-1)))"; exit 0
fi
WORM="${WORM_ARR[$SLURM_ARRAY_TASK_ID]}"

PANEL_H5AD="${PANEL_ROOT}/${PANEL_TAG}/worm_${WORM}.h5ad"
OUT_DIR="${OUT_ROOT}/${PANEL_TAG}/vcf/worm_${WORM}"

echo "[$(date -Iseconds)] task=${SLURM_ARRAY_TASK_ID} worm=${WORM} panel=${PANEL_TAG}"

if [[ ! -f "$PANEL_H5AD" ]]; then
    echo "  SKIP: panel h5ad does not exist: $PANEL_H5AD"; exit 0
fi

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

python scripts/phylogeny/inference/prepare_vcf.py \
    --panel-h5ad "$PANEL_H5AD" \
    --joint-vcf  "$JOINT_VCF" \
    --out-dir    "$OUT_DIR" \
    --seed       "$SEED" \
    --skip-existing

echo "[$(date -Iseconds)] done"
