#!/bin/bash
#SBATCH -J compile_sc_qc_test
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 0:30:00
#SBATCH -p short

##########################################################################################################################
# Concatenate per-cell QC CSVs into a single sc_qc_metrics.csv.
##########################################################################################################################

set -euo pipefail

PER_CELL_CSV_DIR="$1"   # results/L4-4_test3/per_cell_qc
OUTPUT_CSV="$2"         # results/L4-4_test3/sc_qc_metrics.csv

mkdir -p SLURM_outs

n=$(ls "${PER_CELL_CSV_DIR}"/*.csv 2>/dev/null | wc -l)
echo "Found ${n} per-cell CSVs in ${PER_CELL_CSV_DIR}"
if [ "${n}" -eq 0 ]; then
    echo "ERROR: no per-cell CSVs found" >&2
    exit 1
fi

# Header from first file, then all data rows from all files (skip their headers).
# Use a glob -> array to avoid `ls | head -1` getting SIGPIPE under `set -o pipefail`.
shopt -s nullglob
csvs=( "${PER_CELL_CSV_DIR}"/*.csv )
shopt -u nullglob
head -1 "${csvs[0]}" > "${OUTPUT_CSV}"
for f in "${csvs[@]}"; do
    tail -n +2 "${f}" >> "${OUTPUT_CSV}"
done

rows=$(($(wc -l < "${OUTPUT_CSV}") - 1))
echo "Wrote ${OUTPUT_CSV} with ${rows} cell rows."
