#!/bin/bash
set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter
PIPE=/home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/polyAT_bidirectional_pipeline
python /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/scripts/CapGTA/polyAT_count.py \
    --fasta /shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa \
    --fastq-root /home/dmullane/SrivatsanLab/Dustin/SPC_genome/results/worm6_final/GEX/polyAT_bidirectional_fastqs \
    --cells-list /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/polyAT_bidir_cells.txt \
    --anno-dir "$PIPE/featurecounts_annotated" \
    --out-dir "$PIPE/qc/polyAT_tract_overlap" \
    --min-run 15
