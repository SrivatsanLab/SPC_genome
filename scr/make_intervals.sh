#!/bin/bash

FASTA_FILE=$1
WINDOW_SIZE=$2      # Size of each interval in base pairs (e.g., 1000000 for 1Mb)
OUTPUT_BED_FILE=$3

while read -r chrom len; do
    start=0
    while [ $start -lt $len ]; do
        end=$((start + WINDOW_SIZE))
        if [ $end -gt $len ]; then
            end=$len
        fi
        echo -e "$chrom\t$start\t$end" >> $OUTPUT_BED_FILE
        start=$end
    done
done < <(cut -f1,2 $FASTA_FILE.fai)