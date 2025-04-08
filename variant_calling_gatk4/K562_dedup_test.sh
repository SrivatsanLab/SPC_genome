#!/bin/bash
#SBATCH --job-name=K562_dedup_test       # Job name
#SBATCH --output=K562_dedup_test_%j.out  # Output file
#SBATCH --error=K562_dedup_test_%j.err   # Error file
#SBATCH --cpus-per-task=16       # Number of CPUs
#SBATCH --mem=16G                # Memory allocation
#SBATCH --time=12:00:00          # Maximum runtime

# Directory containing the .sam files
input_dir="K562_split_test"
output_dir="K562_dedup_test"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# # Path to GATK
# gatk_path="/path/to/gatk"

# Iterate over each .sam file in the input directory
for sam_file in "$input_dir"/*.sam; do
    # Extract the base filename (without extension)
    base_name=$(basename "$sam_file" .sam)

    # Run GATK MarkDuplicatesSpark
    echo "Processing $sam_file..."
    python ../add_rg_header2.py "$sam_file" "${sam_file%.sam}_processed.sam"

done

for sam_file in "$input_dir"/*_processed.sam; do

    echo "Dedupping $sam_file..."
    base_name=$(basename "$sam_file" _processed.sam)
    # Define output files
    metrics_file="$output_dir/${base_name}_dedup_metrics.txt"
    output_bam="$output_dir/${base_name}_dedup_reads.bam"
    gatk MarkDuplicatesSpark \
        -I "$sam_file" \
        -M "$metrics_file" \
        -O "$output_bam"

done
echo "All files finished."
