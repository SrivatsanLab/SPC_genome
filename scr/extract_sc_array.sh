#!/bin/bash
#SBATCH -J extract_sc
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -p short
#SBATCH -t 3:00:00

##########################################################################################################################
# This job extracts reads from single cells from a bulked sam file                                                       #
##########################################################################################################################

module load SAMtools

input_BAM="$1"
barcode_file="$2"

# Get the barcode for the current array task
barcode=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$barcode_file")

dir_name=$(dirname "$input_BAM")
base_name=$(basename "$input_BAM" .bam)

output_bam="${dir_name}/${base_name}_${barcode}.bam"

# samtools view -h -b -d CB:${barcode} "${input_BAM}" > "${output_bam}"

samtools view -h ${input_BAM} | grep -E "^@|CB:Z:${barcode}" | samtools view -b -o "${output_bam}"