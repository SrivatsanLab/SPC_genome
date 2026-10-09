#!/bin/bash
set -euo pipefail
P="$1"
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter
DIR="/home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/${P}_pipeline_relaxed_fc"
python /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/scripts/CapGTA/check_r1_strand.py \
    --anno-dir "${DIR}/featurecounts_annotated" \
    --counts-raw "${DIR}/counts_raw.txt" \
    --n-cells 200
