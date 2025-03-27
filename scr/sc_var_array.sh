#!/bin/bash
#SBATCH -J sc_var
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 2

##########################################################################################################################
# This job performs preliminary single sample variant calling                         									 #
##########################################################################################################################

module load SAMtools GATK

TMP_dir="$1"
barcode_file="$2"
output_dir="$3"
genome="$4"

SC_outputs="${OUTPUT_dir}/sc_outputs"

mkdir -p "${SC_outputs}"

# Get the barcode for the current array task
barcode=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$barcode_file")
output_bam="${SC_outputs}/${barcode}.bam"
output_vcf="${SC_outputs}/${barcode}.g.vcf.gz"

# Merge chunks into single cell bam:
echo "Merging chunks for cell: ${barcode}"

bam_list=$("${TMP_dir}/*${barcode}*.bam")

samtools merge -o "${output_bam}" -b "${bam_list}"
samtools index "${output_bam}"

# variant Calling with GATK
echo "Variant Calling for cell: ${barcode}"

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R "${genome}" \
   -I "${output_bam}" \
   -O "${output_vcf}" \
   -ERC GVCF