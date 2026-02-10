#!/bin/bash
#SBATCH -J extract_sc
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -t 12:00:00

##########################################################################################################################
# This job extracts reads from single cells from a bulked sam file                                                       #
##########################################################################################################################

module load SAMtools

input_file="$1"  # Can be SAM or BAM
barcode_file="$2"
output_dir="${3:-}"  # Optional: explicit output directory

# Get the barcode for the current array task
barcode=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$barcode_file")

# Determine output directory
if [ -n "${output_dir}" ]; then
    # Use explicitly provided output directory
    sc_output_dir="${output_dir}"
else
    # Default: create sc_outputs next to input file
    dir_name=$(dirname "$input_file")
    sc_output_dir="${dir_name}/sc_outputs"
fi

mkdir -p "${sc_output_dir}"

# Output directly as barcode.bam (no sample prefix)
output_bam="${sc_output_dir}/${barcode}.bam"

# samtools view -h -b -d CB:${barcode} "${input_file}" > "${output_bam}"

samtools view -h ${input_file} | grep -E "^@|CB:Z:${barcode}" | samtools view -b -o "${output_bam}"

# Index the BAM file
samtools index "${output_bam}"