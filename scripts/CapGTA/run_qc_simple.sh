#!/bin/bash
set -euo pipefail
DIR="$1"
STATS="$2"
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter
EXTRA=""
if [ -n "$STATS" ] && [ -f "$STATS" ]; then EXTRA="--extract-stats $STATS"; fi
mkdir -p "${DIR}/qc"
python /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/scripts/CapGTA/qc_primer_trimmed.py \
    --h5ad "${DIR}/umi_matrix.h5ad" \
    --counts-raw "${DIR}/counts_raw.txt" \
    --anno-dir "${DIR}/featurecounts_annotated" \
    --star-logs-dir "${DIR}/star_logs" \
    $EXTRA \
    --out-dir "${DIR}/qc"
