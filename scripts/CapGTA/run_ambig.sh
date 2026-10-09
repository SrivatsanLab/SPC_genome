#!/bin/bash
set -euo pipefail
P="$1"
EXT="$2"
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter
DIR="/home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/${P}_pipeline_relaxed_fc"
python /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/scripts/CapGTA/analyze_ambig_pairs.py \
    --anno-dir "${DIR}/featurecounts_annotated" \
    --counts-raw "${DIR}/counts_raw.txt" \
    --n-cells 300 \
    --gene-extension ${EXT}
