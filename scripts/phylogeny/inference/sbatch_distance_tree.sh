#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=SLURM_outs/inference/%x_%A_%a.out
#SBATCH --error=SLURM_outs/inference/%x_%A_%a.err
#SBATCH --array=0-15
#
# Distance-based NJ + UPGMA with soft-PL Hamming (default) or hard Hamming.
# Bootstrap = 100 resamples of variants with replacement.

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${PANEL_TAG:?PANEL_TAG must be set}"
METRIC="${METRIC:-soft}"      # soft | hard
N_BOOT="${N_BOOT:-100}"
PANEL_ROOT="${PANEL_ROOT:-results/worm6_final/DNA_analysis/phylogeny/panels}"
OUT_ROOT="${OUT_ROOT:-results/worm6_final/DNA_analysis/phylogeny/inference}"

WORMS="${WORMS:-worm01 worm03 worm05 worm06 worm07 worm08 worm09 worm10 worm11 worm12 worm13 worm15 worm16 worm17 worm18 worm20}"
readarray -t WORM_ARR < <(printf '%s\n' $WORMS)
N=${#WORM_ARR[@]}
if (( SLURM_ARRAY_TASK_ID >= N )); then exit 0; fi
WORM="${WORM_ARR[$SLURM_ARRAY_TASK_ID]}"

PANEL_H5AD="${PANEL_ROOT}/${PANEL_TAG}/worm_${WORM}.h5ad"
VCF="${OUT_ROOT}/${PANEL_TAG}/vcf/worm_${WORM}/real.vcf.gz"
OUT_DIR="${OUT_ROOT}/${PANEL_TAG}/distance/worm_${WORM}"
[[ ! -f "$PANEL_H5AD" ]] && { echo "no panel"; exit 0; }

VCF_ARG=()
if [[ -f "$VCF" && "$METRIC" == "soft" ]]; then
    VCF_ARG=(--vcf "$VCF")
fi

# Skip if v2 outputs exist (all pipelines × both rootings)
if [[ -f "${OUT_DIR}/${METRIC}__nj_nni__og.support" \
   && -f "${OUT_DIR}/${METRIC}__bme_nni__og.support" \
   && -f "${OUT_DIR}/${METRIC}__nj_nni__mad.support" ]]; then
    echo "skip: outputs exist"; exit 0
fi

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

# scikit-bio is not in the base env; install if missing
python - <<'PY'
try:
    import skbio  # noqa
except ImportError:
    import subprocess; subprocess.check_call(["pip","install","-q","scikit-bio","dendropy"])
PY

python scripts/phylogeny/inference/run_distance_tree.py \
    --panel-h5ad "$PANEL_H5AD" \
    "${VCF_ARG[@]}" \
    --out-dir "$OUT_DIR" \
    --metric "$METRIC" \
    --n-bootstraps "$N_BOOT" \
    --seed 0
