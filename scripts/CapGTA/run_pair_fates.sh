#!/bin/bash
set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter
PIPE=/home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/polyAT_bidirectional_pipeline
python /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/scripts/CapGTA/classify_pair_fates.py \
    --fastq-root /home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/polyAT_bidirectional_fastqs \
    --cells-list /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/polyAT_bidir_cells.txt \
    --counts-raw "$PIPE/counts_raw.txt" \
    --gtf /shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf \
    --anno-dir "$PIPE/featurecounts_annotated" \
    --out-dir "$PIPE/qc/pair_fates" \
    --near-window 200
