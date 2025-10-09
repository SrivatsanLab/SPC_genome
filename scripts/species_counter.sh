#!/bin/bash

module load SAMtools

# Check if input BAM file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input_bam_file>"
    exit 1
fi

bam_file="$1"  # Input BAM file
output_file="barnyard.csv"  # Output CSV file

# Temporary files
mouse_reads="mouse_reads.txt"
human_reads="human_reads.txt"

# Step 1: Extract reads aligned to mouse and human chromosomes
samtools view "$bam_file" | awk '/chr.*_mm10/ {print}' > "$mouse_reads"
samtools view "$bam_file" | awk '/chr.*_hg38/ {print}' > "$human_reads"

# Step 2: Count reads per cell barcode
count_reads() {
    local input_file=$1
    local output_file=$2

    awk '{
        for (i=12; i<=NF; i++) {
            if ($i ~ /^CB:Z:/) {
                split($i, arr, ":");
                barcode=arr[3];
                count[barcode]++;
            }
        }
    }
    END {
        for (barcode in count) {
            print barcode "," count[barcode];
        }
    }' "$input_file" | sort > "$output_file"
}

count_reads "$mouse_reads" "mouse_counts.csv"
count_reads "$human_reads" "human_counts.csv"

# Combine results
echo "barcode,mouse_count,human_count" > "$output_file"
join -t, -a 1 -a 2 -e 0 -o 0,1.2,2.2 mouse_counts.csv human_counts.csv

# Cleanup
rm "$mouse_reads" "$human_reads" mouse_counts.csv human_counts.csv