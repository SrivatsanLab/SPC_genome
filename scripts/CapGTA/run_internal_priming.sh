#!/bin/bash
set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter
PIPE=/home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/polyAT_bidirectional_pipeline
python /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/scripts/CapGTA/analyze_internal_priming.py \
    --offsets-csv "$PIPE/qc/per_xp_full/r2a_mate_offsets_with_tx_nearest.csv" \
    --counts-raw "$PIPE/counts_raw.txt" \
    --fasta /shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa \
    --out-dir "$PIPE/qc/internal_priming" \
    --near-window 200 \
    --far-threshold 500 \
    --check-window 30 \
    --min-run 10
