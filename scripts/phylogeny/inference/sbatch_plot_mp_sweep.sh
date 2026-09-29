#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --output=SLURM_outs/inference/%x_%A_%a.out
#SBATCH --error=SLURM_outs/inference/%x_%A_%a.err
#SBATCH --array=0-17

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

SWEEP_ROOT="${SWEEP_ROOT:-scratch/mp_sweep}"
OUT_DIR="${OUT_DIR:-results/worm6_final/DNA_analysis/phylogeny/figures/mp_sweep}"

PANELS=(
    all_variants_relax03
    all_variants_relax03__svdImp
    all_variants_relax03__svdImp_mquad10_cs3
    all_variants_relax03__svdImp_mquad10_cs5
    all_variants_relax03__svdImp_mquad5_cs5
    all_variants_relax03__svdImpSym
    all_variants_relax03__svdImpSym_mquad10_cs3
    all_variants_relax03__svdImpSym_mquad10_cs5
    all_variants_relax03__svdImpSym_mquad5_cs5
    tct_kept_relax03
    tct_kept_relax03__svdImp
    tct_kept_relax03__svdImp_mquad10_cs3
    tct_kept_relax03__svdImp_mquad10_cs5
    tct_kept_relax03__svdImp_mquad5_cs5
    tct_kept_relax03__svdImpSym
    tct_kept_relax03__svdImpSym_mquad10_cs3
    tct_kept_relax03__svdImpSym_mquad10_cs5
    tct_kept_relax03__svdImpSym_mquad5_cs5
)

N=${#PANELS[@]}
if (( SLURM_ARRAY_TASK_ID >= N )); then exit 0; fi
PANEL="${PANELS[$SLURM_ARRAY_TASK_ID]}"

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

python scripts/phylogeny/inference/plot_mp_sweep.py \
    --panel-tag "$PANEL" \
    --sweep-root "$SWEEP_ROOT" \
    --out-dir "$OUT_DIR"
