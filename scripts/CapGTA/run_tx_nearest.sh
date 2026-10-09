#!/bin/bash
set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter
PIPE=/home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/polyAT_bidirectional_pipeline
export COUNTS_RAW="$PIPE/counts_raw.txt"
python /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/scripts/CapGTA/plot_tx_nearest_3p.py \
    --offsets-csv "$PIPE/qc/per_xp_full/r2a_mate_offsets.csv" \
    --gtf /shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf \
    --out-dir "$PIPE/qc/per_xp_full"
