#!/bin/bash
# Submit the full pipeline for one panel tag as a chained slurm dependency graph:
#
#   build_panel  →  prepare_vcf (16 tasks)  →  cellphy (48 tasks)  →  { plot, diagnostics }
#
# Panel diagnostics can run in parallel with build_panel's downstream since it
# only depends on build.
#
# Usage:
#   scripts/phylogeny/inference/submit_full_pipeline.sh <CONFIG> <PANEL_TAG>
#
# Prints one line per stage: STAGE=JOBID

set -euo pipefail
CONFIG="$1"
PANEL_TAG="$2"

cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome

BUILD=$(sbatch --job-name="panel_${PANEL_TAG}" \
    --export=ALL,CONFIG="${CONFIG}",PANEL_TAG="${PANEL_TAG}" \
    --parsable scripts/phylogeny/inference/sbatch_build_panel.sh)
echo "build_panel=${BUILD}"

PREP=$(sbatch --job-name="prep_${PANEL_TAG}" \
    --export=ALL,PANEL_TAG="${PANEL_TAG}" \
    --dependency=afterok:${BUILD} \
    --parsable scripts/phylogeny/inference/sbatch_prepare.sh)
echo "prepare_vcf=${PREP}"

CP=$(sbatch --job-name="cp_${PANEL_TAG}" \
    --export=ALL,PANEL_TAG="${PANEL_TAG}" \
    --dependency=afterok:${PREP} \
    --parsable scripts/phylogeny/inference/sbatch_cellphy.sh)
echo "cellphy=${CP}"

PLOT=$(sbatch --job-name="plot_${PANEL_TAG}" \
    --export=ALL,PANEL_TAG="${PANEL_TAG}" \
    --dependency=afterany:${CP} \
    --parsable scripts/phylogeny/inference/sbatch_plot_cellphy.sh)
echo "plot=${PLOT}"

DIAG=$(sbatch --job-name="diag_${PANEL_TAG}" \
    --export=ALL,PANEL_TAG="${PANEL_TAG}" \
    --dependency=afterany:${CP} \
    --parsable scripts/phylogeny/inference/sbatch_diagnostics.sh)
echo "diagnostics=${DIAG}"

PANELDIAG=$(sbatch --job-name="paneldiag_${PANEL_TAG}" \
    --export=ALL,PANEL_TAG="${PANEL_TAG}" \
    --dependency=afterok:${BUILD} \
    --parsable scripts/phylogeny/inference/sbatch_panel_diagnostics.sh)
echo "panel_diagnostics=${PANELDIAG}"
