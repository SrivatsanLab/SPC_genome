#!/bin/bash
#SBATCH -J polyat_count
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -p campus-new
#SBATCH -c 16
#SBATCH --mem=8G
#SBATCH -t 4:00:00

###########################################################################################################################
# Count polyA/T reads (SEQ contains A{10,} or T{10,}) per cell in an h5ad and (optionally) inject
# the counts as a new obs column. See count_polyat_per_cell.py for the rule + safety design.
#
# Usage:
#   sbatch scripts/CapGTA/count_polyat_per_cell.sh \
#       <adata_h5ad> <out_csv> <mm_env> [--inject]
###########################################################################################################################

set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage: sbatch $0 <adata_h5ad> <out_csv> <mm_env> [--inject]" >&2
    exit 1
fi

H5AD="$1"
OUT_CSV="$2"
MMENV="$3"
shift 3
EXTRA=("$@")

PY="${SLURM_SUBMIT_DIR}/scripts/CapGTA/count_polyat_per_cell.py"

for f in "${H5AD}" "${PY}"; do
    [ -e "${f}" ] || { echo "Error: not found: ${f}" >&2; exit 1; }
done

module load SAMtools

micromamba run -n "${MMENV}" python3 "${PY}" \
    "${H5AD}" "${OUT_CSV}" \
    --threads 16 --samtools-threads 2 --min-run 10 \
    "${EXTRA[@]}"
