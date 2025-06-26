#!/bin/bash
#SBATCH --job-name=K562_calling_test       # Job name
#SBATCH --output=K562_calling_test_%j.out  # Output file
#SBATCH --error=K562_calling_test_%j.err   # Error file
#SBATCH --cpus-per-task=16               # Number of CPUs
#SBATCH --mem=16G                        # Memory allocation
#SBATCH --time=12:00:00                  # Maximum runtime

# Define variables
input_dir="K562_dedup_test"
output_dir="K562_gvcfs"
reference="../hg38_new/hg38_new.fasta"  # Replace with the actual reference genome path
bed1="../BIRC6.bed"
bed2="../FANCC.bed"
bed3="../BRCA1.bed"
bed4="../AKT3_chr1.bed"
mkdir -p "$output_dir"  # Create output directory if it doesn't exist

# Loop through BAM files in the input directory
for bam_file in "$input_dir"/*.bam; do
    # Extract the base name of the BAM file (e.g., file.bam -> file)
    base_name=$(basename "$bam_file" .bam)
    
    # Define the output GVCF file name
    output_gvcf="$output_dir/${base_name}.g.vcf.gz"
    
    # Run GATK HaplotypeCaller
    gatk HaplotypeCaller \
        -R "$reference" \
        -I "$bam_file" \
        -O "$output_gvcf" \
        -ERC GVCF \
        --native-pair-hmm-threads 16 \
        -L "$bed1" \
        -L "$bed2" \
        -L "$bed3" \
        -L "$bed4"
done
