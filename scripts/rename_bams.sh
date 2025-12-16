#!/bin/bash
# Move and rename sc BAMs, removing sample prefix

for bam in data/HSC_CellGrowth_coverage/aligned/HSC_CellGrowth_coverage_*.bam; do
  barcode=$(basename "$bam" .bam | sed 's/HSC_CellGrowth_coverage_//')
  mv "$bam" "data/HSC_CellGrowth_coverage/sc_outputs/${barcode}.bam"
  if [ -f "${bam}.bai" ]; then
    mv "${bam}.bai" "data/HSC_CellGrowth_coverage/sc_outputs/${barcode}.bam.bai"
  fi
done

echo "BAMs moved and renamed"
ls data/HSC_CellGrowth_coverage/sc_outputs/*.bam | wc -l
