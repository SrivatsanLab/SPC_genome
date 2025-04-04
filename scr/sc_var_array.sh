#!/bin/bash
#SBATCH -J sc_var
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out

##########################################################################################################################
# This job performs preliminary single sample variant calling                         									 #
##########################################################################################################################

module load SAMtools GATK

TMP_dir="$1"
barcode_file="$2"
genome="$3"
SC_outputs="$4"
mkdir -p "${SC_outputs}"

# Get the barcode for the current array task
barcode=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$barcode_file")
output_bam="${SC_outputs}/${barcode}.bam"
output_vcf="${SC_outputs}/${barcode}.g.vcf.gz"

# Merge chunks into single cell bam:
echo "Merging chunks for cell: ${barcode}"

ls ${TMP_dir}/*${barcode}.bam > "${TMP_dir}/${barcode}_files.txt"

samtools merge -o "${output_bam}" -b "${TMP_dir}/${barcode}_files.txt"
samtools index "${output_bam}"

# variant Calling with GATK
echo "Variant Calling for cell: ${barcode}"

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R "${genome}" \
   -I "${output_bam}" \
   -O "${output_vcf}" \
   -ERC GVCF