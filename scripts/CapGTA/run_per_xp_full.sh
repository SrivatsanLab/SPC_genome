#!/bin/bash
set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter
PIPE=/home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/polyAT_bidirectional_pipeline
mkdir -p "$PIPE/qc/per_xp_full"
python /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/scripts/CapGTA/plot_per_xp_full.py \
    --fastq-root /home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/polyAT_bidirectional_fastqs \
    --cells-list /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/polyAT_bidir_cells.txt \
    --counts-raw "$PIPE/counts_raw.txt" \
    --anno-dir "$PIPE/featurecounts_annotated" \
    --out-dir "$PIPE/qc/per_xp_full"
