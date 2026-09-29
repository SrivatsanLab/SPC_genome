#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=1:00:00
#SBATCH --output=SLURM_outs/inference/%x_%j.out
#SBATCH --error=SLURM_outs/inference/%x_%j.err
#
# Post-process a top-down SVD sweep: rerun summary + shared/private diagnostic,
# then plot per-config tree grids.
#
# Usage:
#   sbatch --dependency=afterany:<sweep_jobid> \
#          --job-name=svdtd_postproc \
#          --export=ALL,PANEL_TAG=all_variants_relax03 \
#          scripts/phylogeny/inference/sbatch_topdown_postproc.sh

set -euo pipefail
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome
mkdir -p SLURM_outs/inference

: "${PANEL_TAG:?PANEL_TAG must be set}"
TD_ROOT="results/worm6_final/DNA_analysis/phylogeny/topdown"
PANEL_ROOT="results/worm6_final/DNA_analysis/phylogeny/panels"
FIG_ROOT="${TD_ROOT}/${PANEL_TAG}/figures"

eval "$(micromamba shell hook -s bash)"
micromamba activate cellphy

echo "[postproc] summarize_topdown_sweep"
python scripts/phylogeny/inference/summarize_topdown_sweep.py --panel-tag "$PANEL_TAG"

echo "[postproc] shared_private diagnostic"
python scripts/phylogeny/inference/topdown_shared_private_diag.py --panel-tag "$PANEL_TAG"

echo "[postproc] tree.json -> Newick conversion"
python scripts/phylogeny/inference/topdown_tree_to_newick.py --panel-tag "$PANEL_TAG"

echo "[postproc] plotting per-config trees"
mkdir -p "$FIG_ROOT"
for cfg_dir in "${TD_ROOT}/${PANEL_TAG}"/*/; do
    cfg=$(basename "$cfg_dir")
    [[ "$cfg" == "figures" ]] && continue
    [[ "$cfg" == "_summaries" ]] && continue
    # only render configs that actually have >0 trees
    n_trees=$(find "$cfg_dir" -name tree.json 2>/dev/null | wc -l)
    if (( n_trees == 0 )); then
        echo "  [skip] $cfg: no trees"
        continue
    fi
    echo "  [plot] $cfg ($n_trees trees)"
    python scripts/phylogeny/inference/plot_svd_trees.py \
        --inference-dir "${TD_ROOT}/${PANEL_TAG}" \
        --svd-subdir "$cfg" \
        --panels-dir "${PANEL_ROOT}/${PANEL_TAG}" \
        --out-dir "${FIG_ROOT}/${cfg}" || echo "  [warn] plot failed for $cfg"
done

echo "[postproc] done. summaries: ${TD_ROOT}/_summaries/  figs: ${FIG_ROOT}/"
