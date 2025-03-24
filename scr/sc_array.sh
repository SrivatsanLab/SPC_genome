#!/bin/bash
#SBATCH -J sc_array
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH -p short

##########################################################################################################################
# This job extracts reads from single cells, performs preliminary single sample variant calling                          #
##########################################################################################################################

module load SAMtools GATK

input_bam="$1"
barcode_file="$2"
OUTPUT_DIR="$3"
SCRIPTS_DIR="$4"
genome="$5"

SC_outputs="${OUTPUT_DIR}/sc_outputs"
mkdir -p "${SC_outputs}"

# Get the barcode for the current array task
barcode=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$barcode_file")
output_bam="${SC_outputs}/${barcode}.bam"
output_vcf="${SC_outputs}/${barcode}.g.vcf.gz"

echo "Processing barcode ${barcode} in SLURM task ${SLURM_ARRAY_TASK_ID}"

## For some reason -d isn't working since switching to atrandi combinatorial indices, idk why??
# samtools view -h -b -@ 8 -d CB:${barcode} "${input_bam}" > "${output_bam}"
## Workaround:
samtools view -h ${input_bam} | grep -E "^@|CB:Z:${barcode}" | samtools view -b -o "${output_bam}"

samtools index "${output_bam}"

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R "${}" \
   -I "${output_bam}" \
   -O "${output_vcf}" \
   -ERC GVCF