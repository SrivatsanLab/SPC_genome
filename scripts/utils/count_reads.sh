#!/bin/bash
DATA_DIR="/fh/fast/srivatsan_s/pub/projects/00_genome_transcriptome_coassay/res/260202_VH00738_369_AAHML2JM5_L4_worm_pilot"
for sample in L4_worm_sci_5 L4_worm_sci_6 L4_worm_sci_7 L4_worm_sci_8; do
    r1=$(ls ${DATA_DIR}/${sample}_S*_R1_*.fastq.gz)
    lines=$(zcat "$r1" | wc -l)
    reads=$((lines / 4))
    echo "${sample}: ${reads} reads"
done
