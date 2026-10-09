#!/bin/bash
#SBATCH -J primer_trim_agg
#SBATCH -o SLURM_outs/array_outs/%x_%j.out
#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 00:30:00
#
# After the primer_trim array completes, roll up <barcode>_stats.json into
# <output_root>/all_cells_stats.tsv.
#
# Usage: sbatch --dependency=afterok:<array_jobid> \
#            scripts/CapGTA/submit_primer_trim_aggregate.sh [output_root]

set -euo pipefail
OUTPUT_ROOT="${1:-$HOME/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/primer_trimmed_fastqs}"
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "${SCRIPT_DIR}/aggregate_primer_trim_stats.py" "${OUTPUT_ROOT}"
echo "aggregate done at $(date)"
