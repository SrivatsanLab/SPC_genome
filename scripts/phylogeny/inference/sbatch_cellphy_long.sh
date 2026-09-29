#!/bin/bash
#SBATCH --partition=campus-new
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=SLURM_outs/inference/%x_%A_%a.out
#SBATCH --error=SLURM_outs/inference/%x_%A_%a.err
#
# Long-time / larger-panel variant of sbatch_cellphy.sh. Same workflow, same
# outputs — just 12h walltime on campus-new. Use for timed-out reruns.
#
# Usage:
#   sbatch --job-name=cp_long_<tag> \
#       --export=ALL,PANEL_TAG=<tag>,WORMS="worm07 worm10",MATRICES=real \
#       --array=0-<N-1> \
#       scripts/phylogeny/inference/sbatch_cellphy_long.sh

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${PANEL_TAG:?PANEL_TAG must be set}"
OUT_ROOT="${OUT_ROOT:-results/worm6_final/DNA_analysis/phylogeny/inference}"
BOOTSTRAP="${BOOTSTRAP:-100}"

WORMS="${WORMS:-worm01 worm03 worm05 worm06 worm07 worm08 worm09 worm10 worm11 worm12 worm13 worm15 worm16 worm17 worm18 worm20}"
MATRICES="${MATRICES:-real covmask shuffled}"
readarray -t WORM_ARR < <(printf '%s\n' $WORMS)
readarray -t MAT_ARR  < <(printf '%s\n' $MATRICES)
NW=${#WORM_ARR[@]}; NM=${#MAT_ARR[@]}
N=$((NW * NM))
if (( SLURM_ARRAY_TASK_ID >= N )); then exit 0; fi
WORM="${WORM_ARR[$((SLURM_ARRAY_TASK_ID / NM))]}"
MATRIX="${MAT_ARR[$((SLURM_ARRAY_TASK_ID % NM))]}"

VCF="${OUT_ROOT}/${PANEL_TAG}/vcf/worm_${WORM}/${MATRIX}.vcf.gz"
OUT_DIR="${OUT_ROOT}/${PANEL_TAG}/cellphy/worm_${WORM}/${MATRIX}"
mkdir -p "$OUT_DIR"

if [[ ! -f "$VCF" ]]; then echo "no VCF"; exit 0; fi
if [[ -f "${OUT_DIR}/sup.raxml.supportFBP" ]]; then echo "skip: supportFBP done"; exit 0; fi

# Clean any partial run since we want a fresh start.
rm -f "${OUT_DIR}"/run.raxml.* "${OUT_DIR}"/bs.raxml.* "${OUT_DIR}"/sup.raxml.*

CELLPHY=/home/dmullane/cellphy/cellphy.sh
RAXML=/home/dmullane/cellphy/bin/raxml-ng-cellphy-linux

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

cp "$VCF" "${OUT_DIR}/input.vcf.gz"
cp "${VCF}.tbi" "${OUT_DIR}/input.vcf.gz.tbi" 2>/dev/null || true

pushd "$OUT_DIR" > /dev/null
"$CELLPHY" SEARCH -t "$SLURM_CPUS_PER_TASK" -p run -r input.vcf.gz 2>&1 | tee cellphy.log
if [[ -f run.raxml.bestModel && -f run.raxml.bestTree ]]; then
    "$RAXML" --bootstrap --msa input.vcf.gz --model run.raxml.bestModel \
        --bs-trees "$BOOTSTRAP" --threads "$SLURM_CPUS_PER_TASK" \
        --prefix bs --seed 0 2>&1 | tee bootstrap.log
    if [[ -f bs.raxml.bootstraps ]]; then
        "$RAXML" --support --tree run.raxml.bestTree --bs-trees bs.raxml.bootstraps \
            --bs-metric fbp,tbe --prefix sup 2>&1 | tee support.log
    fi
fi
popd > /dev/null
