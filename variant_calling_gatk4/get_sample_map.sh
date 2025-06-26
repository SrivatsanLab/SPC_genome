#!/bin/bash

for file in ./K562_gvcfs/*.g.vcf.gz; do
    sample=$(basename "$file" .g.vcf.gz)  # Extract sample name
    echo -e "$sample\t$file" >> sample_map_K562.txt
done
