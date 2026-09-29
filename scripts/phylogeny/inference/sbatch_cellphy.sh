#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=6:00:00
#SBATCH --output=SLURM_outs/inference/%x_%A_%a.out
#SBATCH --error=SLURM_outs/inference/%x_%A_%a.err
#SBATCH --array=0-47
#
# CellPhy SEARCH (GT16 GL mode + bootstrap 100) over worm x matrix.
# Default array covers 16 worms x 3 matrices = 48 tasks.
#
# Usage:
#   sbatch --job-name=cp_v0 \
#       --export=ALL,PANEL_TAG=v0_original \
#       scripts/phylogeny/inference/sbatch_cellphy.sh
#
# Env:
#   PANEL_TAG    subdir under phylogeny/inference/ (required)
#   OUT_ROOT     phylogeny/inference root (default: repo layout)
#   WORMS        space-separated worm IDs (default: 16-worm set)
#   MATRICES     space-separated matrix names (default: real covmask shuffled)
#   BOOTSTRAP    replicate count (default 100)

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
NW=${#WORM_ARR[@]}
NM=${#MAT_ARR[@]}
N=$((NW * NM))

if (( SLURM_ARRAY_TASK_ID >= N )); then
    echo "task id $SLURM_ARRAY_TASK_ID out of range (0..$((N-1)))"; exit 0
fi
WORM="${WORM_ARR[$((SLURM_ARRAY_TASK_ID / NM))]}"
MATRIX="${MAT_ARR[$((SLURM_ARRAY_TASK_ID % NM))]}"

VCF="${OUT_ROOT}/${PANEL_TAG}/vcf/worm_${WORM}/${MATRIX}.vcf.gz"
OUT_DIR="${OUT_ROOT}/${PANEL_TAG}/cellphy/worm_${WORM}/${MATRIX}"
mkdir -p "$OUT_DIR"

echo "[$(date -Iseconds)] task=${SLURM_ARRAY_TASK_ID} worm=${WORM} matrix=${MATRIX} cpus=${SLURM_CPUS_PER_TASK}"

if [[ ! -f "$VCF" ]]; then
    echo "  SKIP: input VCF does not exist: $VCF"; exit 0
fi
if [[ -f "${OUT_DIR}/run.raxml.support" ]] || [[ -f "${OUT_DIR}/run.raxml.bestTree" ]]; then
    echo "  SKIP: outputs already exist"; exit 0
fi

# Free the bundled raxml-ng binary to consume PL directly (GL mode is default;
# `-l` opts into ML mode, which we do NOT want).
CELLPHY=/home/dmullane/cellphy/cellphy.sh
if [[ ! -x "$CELLPHY" ]]; then
    echo "  ERROR: cellphy not found at $CELLPHY"; exit 1
fi

# bcftools is needed by cellphy.sh for VCF preprocessing.
eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

# Copy VCF into per-run dir so cellphy's aux files live alongside it.
cp "$VCF" "${OUT_DIR}/input.vcf.gz"
cp "${VCF}.tbi" "${OUT_DIR}/input.vcf.gz.tbi" 2>/dev/null || true

pushd "$OUT_DIR" > /dev/null
"$CELLPHY" SEARCH \
    -t "$SLURM_CPUS_PER_TASK" \
    -p run \
    -r \
    input.vcf.gz \
    2>&1 | tee cellphy.log

# Bootstrap and support mapping — cellphy.sh does this in FULL mode, but SEARCH
# skips it. Run bootstrap + support via raxml-ng directly against the best tree.
RAXML=/home/dmullane/cellphy/bin/raxml-ng-cellphy-linux
if [[ -f run.raxml.bestModel ]] && [[ -f run.raxml.bestTree ]]; then
    "$RAXML" --bootstrap \
        --msa input.vcf.gz \
        --model run.raxml.bestModel \
        --bs-trees "$BOOTSTRAP" \
        --threads "$SLURM_CPUS_PER_TASK" \
        --prefix bs \
        --seed 0 \
        2>&1 | tee bootstrap.log
    if [[ -f bs.raxml.bootstraps ]]; then
        "$RAXML" --support \
            --tree run.raxml.bestTree \
            --bs-trees bs.raxml.bootstraps \
            --bs-metric fbp,tbe \
            --prefix sup \
            2>&1 | tee support.log
    fi
fi
popd > /dev/null

echo "[$(date -Iseconds)] done"
