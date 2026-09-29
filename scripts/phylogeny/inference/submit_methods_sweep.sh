#!/bin/bash
# Submit SVD + MP + distance-tree runs on a panel tag. Optional --include-cellphy
# triggers prepare_vcf → cellphy → plot as well (needed for freshly-built panels).
#
# Usage:
#   submit_methods_sweep.sh <PANEL_TAG> [--after JOBID] [--include-cellphy]

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome

PANEL_TAG="$1"; shift
DEP_ARG=""
INCLUDE_CELLPHY="no"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --after) DEP_ARG="--dependency=afterok:$2"; shift 2;;
        --include-cellphy) INCLUDE_CELLPHY="yes"; shift;;
        *) echo "unknown arg $1"; exit 1;;
    esac
done

echo "# ${PANEL_TAG}"

SVD=$(sbatch --job-name="svd_${PANEL_TAG}" --export=ALL,PANEL_TAG="$PANEL_TAG" \
    $DEP_ARG --parsable scripts/phylogeny/inference/sbatch_svd.sh)
echo "svd=${SVD}"

MP=$(sbatch --job-name="mp_${PANEL_TAG}" --export=ALL,PANEL_TAG="$PANEL_TAG" \
    $DEP_ARG --parsable scripts/phylogeny/inference/sbatch_mp.sh)
echo "mp=${MP}"

# CellPhy pipeline for freshly-built panels.
if [[ "$INCLUDE_CELLPHY" == "yes" ]]; then
    PREP=$(sbatch --job-name="prep_${PANEL_TAG}" --export=ALL,PANEL_TAG="$PANEL_TAG" \
        $DEP_ARG --parsable scripts/phylogeny/inference/sbatch_prepare.sh)
    echo "prepare_vcf=${PREP}"
    CP=$(sbatch --job-name="cp_${PANEL_TAG}" --export=ALL,PANEL_TAG="$PANEL_TAG" \
        --dependency=afterok:${PREP} --parsable scripts/phylogeny/inference/sbatch_cellphy.sh)
    echo "cellphy=${CP}"
    PLOT=$(sbatch --job-name="plot_${PANEL_TAG}" --export=ALL,PANEL_TAG="$PANEL_TAG" \
        --dependency=afterany:${CP} --parsable scripts/phylogeny/inference/sbatch_plot_cellphy.sh)
    echo "plot=${PLOT}"
    PANELDIAG=$(sbatch --job-name="paneldiag_${PANEL_TAG}" --export=ALL,PANEL_TAG="$PANEL_TAG" \
        $DEP_ARG --parsable scripts/phylogeny/inference/sbatch_panel_diagnostics.sh)
    echo "panel_diagnostics=${PANELDIAG}"

    # Distance NJ/UPGMA depends on the VCF being prepared (soft-PL metric).
    DIST_DEP="--dependency=afterok:${PREP}"
else
    DIST_DEP="$DEP_ARG"
fi

DIST=$(sbatch --job-name="dist_${PANEL_TAG}" --export=ALL,PANEL_TAG="$PANEL_TAG",METRIC=soft,METHOD=both \
    $DIST_DEP --parsable scripts/phylogeny/inference/sbatch_distance_tree.sh)
echo "distance=${DIST}"
