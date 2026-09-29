#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=SLURM_outs/phylo/%x_%j.out
#SBATCH --error=SLURM_outs/phylo/%x_%j.err
# Usage:
#   sbatch --job-name=<run_name> scripts/phylogeny/sweeps/run_config.sh <config.yaml>
# where <config.yaml> is either absolute or relative to scripts/phylogeny/.
# The wrapper cd's to scripts/phylogeny/ first, so bare "configs/foo.yaml" works.
# mquad_nproc autodetects from $SLURM_CPUS_PER_TASK.
set -euo pipefail
CFG_IN="${1:-configs/default.yaml}"
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/scripts/phylogeny
# Absolutify to remove the cd-vs-cwd footgun.
if [[ "$CFG_IN" = /* ]]; then CFG="$CFG_IN"; else CFG="$PWD/$CFG_IN"; fi
if [[ ! -f "$CFG" ]]; then
    echo "ERROR: config not found: $CFG" >&2
    exit 2
fi
echo "[$(date -Iseconds)] config=$CFG cpus=$SLURM_CPUS_PER_TASK job=$SLURM_JOB_ID"
python -m phylo panels --config "$CFG" --force
python -m phylo trees  --config "$CFG" --force
echo "[$(date -Iseconds)] done"
