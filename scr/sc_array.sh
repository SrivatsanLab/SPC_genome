#!/bin/bash
#SBATCH -J sc_array
#SBATCH -o SLURM_outs/array_outs/%x_%j.out
#SBATCH -c 4
#SBATCH -p short

##########################################################################################################################
# This job extracts reads from single cells, performs preliminary single sample variant calling                          #
##########################################################################################################################

BAM_file=""
barcode_file=""
OUTPUT_DIR="$3"
SCRIPTS_DIR="$4"

SC_outputs="${OUTPUT_DIR}/sc_outputs"
mkdir -p "${SC_outputs}"

# Get the barcode for the current array task
barcode=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$barcode_file")
output_bam="${SC_outputs}/${barcode}.bam"

echo "Processing barcode ${barcode} in SLURM task ${SLURM_ARRAY_TASK_ID}"

# Filter BAM for the barcode and save to temp directory

## For some reason -d isn't working since switching to atrandi combinatorial indices, idk why??
# samtools view -h -b -@ 8 -d CB:${barcode} "${input_bam}" > "${output_bam}"
## Workaround:
samtools view -h ${input_bam} | grep -E "^@|CB:Z:${barcode}" | samtools view -b -o "${output_bam}"

samtools index "${output_bam}"