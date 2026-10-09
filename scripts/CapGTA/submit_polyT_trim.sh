#!/bin/bash
#SBATCH -J polyT_trim
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -p short
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 2:00:00

###########################################################################################################################
# Per-cell SLURM array wrapper for extract_primer_trimmed.py in --no-require-umi mode.
#
# Keeps every pair with a left-edge polyT run (same quality-aware criteria as the primer_trim run)
# regardless of whether a valid UMI window exists in front of it. Pairs with a valid UMI window
# still receive the UB:Z:<umi> tag; pairs without one are kept but have no UB tag.
#
# Usage: identical to submit_primer_trim.sh — just points at a different output root.
#   sbatch --array=1-$(wc -l < primer_trim_cells.txt) \
#       scripts/CapGTA/submit_polyT_trim.sh primer_trim_cells.txt
#
# Writes:
#   <output_root>/<UDI>/<barcode>_R1.polyT_trimmed.fastq.gz
#   <output_root>/<UDI>/<barcode>_R2.polyT_trimmed.fastq.gz
#   <output_root>/<UDI>/<barcode>_stats.json
###########################################################################################################################

set -euo pipefail

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID not set — submit as an array." >&2
    exit 1
fi

CELL_LIST="${1:-primer_trim_cells.txt}"
INPUT_ROOT="${2:-/fh/working/srivatsan_s/Dustin/worm6-dubnova-100k-align-r2trim/combined_demux}"
OUTPUT_ROOT="${3:-$HOME/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/polyT_trimmed_fastqs}"
RUN_MIN="${4:-15}"
QUAL_THRESHOLD="${5:-20}"

if [ ! -f "${CELL_LIST}" ]; then
    echo "Error: cell list not found: ${CELL_LIST}" >&2
    exit 1
fi

BARCODE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CELL_LIST}")
if [ -z "${BARCODE}" ]; then
    echo "Error: no barcode at line ${SLURM_ARRAY_TASK_ID} of ${CELL_LIST}" >&2
    exit 1
fi

if [[ ! "${BARCODE}" =~ worm6_UDI[0-9]+$ ]]; then
    echo "Error: barcode '${BARCODE}' does not end in worm6_UDI<n>" >&2
    exit 1
fi
UDI="${BASH_REMATCH[0]}"

R1="${INPUT_ROOT}/${UDI}/results/${UDI}/cells/${BARCODE}_R1.fastq.gz"
R2="${INPUT_ROOT}/${UDI}/results/${UDI}/cells/${BARCODE}_R2.fastq.gz"
for f in "${R1}" "${R2}"; do
    if [ ! -f "${f}" ]; then
        echo "Error: input FASTQ missing: ${f}" >&2
        exit 1
    fi
done

OUT_DIR="${OUTPUT_ROOT}/${UDI}"
mkdir -p "${OUT_DIR}"
OUT_R1="${OUT_DIR}/${BARCODE}_R1.polyT_trimmed.fastq.gz"
OUT_R2="${OUT_DIR}/${BARCODE}_R2.polyT_trimmed.fastq.gz"
OUT_STATS="${OUT_DIR}/${BARCODE}_stats.json"

eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${SCRIPT_DIR}/extract_primer_trimmed.py"
if [ ! -f "${PY}" ]; then
    for cand in \
        "$HOME/SrivatsanLab/Dustin/SPC_genome/scripts/CapGTA/extract_primer_trimmed.py" \
        "./scripts/CapGTA/extract_primer_trimmed.py"; do
        if [ -f "${cand}" ]; then PY="${cand}"; break; fi
    done
fi

echo "=== $(date) : ${BARCODE} (UDI=${UDI}, task=${SLURM_ARRAY_TASK_ID}) [no-require-umi] ==="
python "${PY}" \
    --r1 "${R1}" --r2 "${R2}" \
    --out-r1 "${OUT_R1}" --out-r2 "${OUT_R2}" --out-stats "${OUT_STATS}" \
    --run-min "${RUN_MIN}" --qual-threshold "${QUAL_THRESHOLD}" \
    --no-require-umi
echo "=== $(date) : done ==="
