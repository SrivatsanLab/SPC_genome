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
OUTPUT_FILE="${OUTPUT_DIR}/barcodes.txt"

echo "extracting a list of unique barcodes from ${fastq_file} using ${CPU_threads} CPUs, writing results to ${OUTPUT_FILE}"

# Create a temporary directory for processing chunks
temp_dir=$(mktemp -d)	

headers="${temp_dir}/fastq_headers.txt"

# Extract fastq headers
zgrep '^@' "${fastq_file}" > "${headers}"

declare -i total_lines=$(wc -l < "${headers}")

# Calculate lines per chunk, rounding up
lines_per_chunk=$(expr \( "$total_lines" + "$CPU_threads" - 1 \) / "$CPU_threads")

# Split the headers file into chunks with specified number of lines
split -l "$lines_per_chunk" "${headers}" "${temp_dir}/chunk_"

# Define a function to process each chunk
process_chunk() {
    chunk_file=$1
    sed -n 's/.*CB:Z:\(.*\)/\1/p' "${chunk_file}" | sort | uniq > "${chunk_file}_barcodes.txt"
}

export -f process_chunk


# Use GNU parallel to process each chunk in parallel
find "$temp_dir" -type f -name 'chunk_*' | parallel process_chunk

# Combine results from all chunks
cat "$temp_dir"/chunk_*_barcodes.txt > "${OUTPUT_FILE}"

# Clean up temporary files
rm -rf "${temp_dir}"