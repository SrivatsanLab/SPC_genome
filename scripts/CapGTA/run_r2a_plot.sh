#!/bin/bash
set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter
SIDE=/home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/polyAT_bidirectional_side
python /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/scripts/CapGTA/plot_r2a_mate_positions.py \
    --bams-dir "$SIDE/bams" \
    --fastq-root /home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/polyAT_bidirectional_fastqs \
    --cells-list "$SIDE/cells20.txt" \
    --counts-raw "$SIDE/counts_raw.txt" \
    --anno-dir "$SIDE/featurecounts_annotated" \
    --out-dir "$SIDE/qc"
