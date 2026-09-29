#!/bin/bash
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --output=SLURM_outs/phylo/%x_%j.out
#SBATCH --error=SLURM_outs/phylo/%x_%j.err
# Usage:
#   sbatch --job-name=<name> scripts/phylogeny/sweeps/run_stage.sh <stage> <config.yaml>
# where <stage> is one of {g0, panels, trees, clade_matching, all}.
set -euo pipefail
STAGE="${1:?stage required}"
CFG_IN="${2:?config required}"
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/scripts/phylogeny
if [[ "$CFG_IN" = /* ]]; then CFG="$CFG_IN"; else CFG="$PWD/$CFG_IN"; fi
if [[ ! -f "$CFG" ]]; then echo "ERROR: config not found: $CFG" >&2; exit 2; fi
echo "[$(date -Iseconds)] stage=$STAGE config=$CFG cpus=$SLURM_CPUS_PER_TASK job=$SLURM_JOB_ID"
python -m phylo "$STAGE" --config "$CFG" --force
echo "[$(date -Iseconds)] done"
