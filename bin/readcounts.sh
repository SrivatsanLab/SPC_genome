#!/bin/bash

module load parallel

# Parse command line args
# Check command-line arguments
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 <file.fastq.gz> <CPUs> [ouput_dir (default: .)]"
    exit 1
fi

fastq_file="$1"
declare -i CPU_threads="$2"
OUTPUT_DIR="${3:-.}"

# Outputs
headers="fastq_headers.txt"
output_file="${OUTPUT_DIR}/barcode_readcounts.csv"

echo "counting reads per barcode from ${fastq_file} using ${CPU_threads} CPUs, writing results to ${output_file}"

# Extract fastq headers
zgrep '^@' "${fastq_file}" > "${headers}"

# Create a temporary directory for processing chunks
temp_dir=$(mktemp -d)

declare -i total_lines=$(wc -l < "${headers}")

# Calculate lines per chunk, rounding up
lines_per_chunk=$(expr \( "$total_lines" + "$CPU_threads" - 1 \) / "$CPU_threads")
# Split the headers file into chunks with specified number of lines
split -l "$lines_per_chunk" "${headers}" "${temp_dir}/chunk_"

# Define a function to process each chunk
process_chunk() {
    chunk_file=$1
    sed -n 's/.*CB:Z:\(.*\)/\1/p' "${chunk_file}" | sort | uniq -c > "${chunk_file}_barcodes.txt"
}

export -f process_chunk

# Use GNU parallel to process each chunk in parallel
find "$temp_dir" -type f -name 'chunk_*' | parallel process_chunk

# Combine results from all chunks
cat "$temp_dir"/chunk_*_barcodes.txt | awk '{counts[$2]+=$1} END {for (barcode in counts) print barcode, counts[barcode]}' | sort > "${temp_dir}/combined_barcodes.txt"

# Create the output CSV file
echo "barcode,readcount" > "${output_file}"

# Write combined results to CSV
awk '{print $1 "," $2}' "${temp_dir}/combined_barcodes.txt" >> "${output_file}"

# Clean up temporary files
rm -rf "${temp_dir}"
rm "${headers}"