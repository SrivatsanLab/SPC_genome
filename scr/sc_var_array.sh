#!/bin/bash
#SBATCH -J sc_var
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out

##########################################################################################################################
# This job performs preliminary single sample variant calling                         									 #
##########################################################################################################################

module load SAMtools picard

TMP_dir="$1"
barcode_file="$2"
genome="$3"
SC_outputs="$4"
mkdir -p "${SC_outputs}"

# Get the barcode for the current array task
barcode=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$barcode_file")
sc_bam="${SC_outputs}/${barcode}.bam"
GATK_bam="${SC_outputs}/${barcode}_GATK.bam"
output_vcf="${SC_outputs}/${barcode}.g.vcf.gz"

# Merge chunks into single cell bam:
echo "Merging chunks for cell: ${barcode}"

ls ${TMP_dir}/*${barcode}.bam > "${TMP_dir}/${barcode}_files.txt"

samtools merge -o "${sc_bam}" -b "${TMP_dir}/${barcode}_files.txt"
samtools sort -o "${sc_bam%.bam}.sorted.bam}" "${sc_bam}"
mv "${sc_bam%.bam}.sorted.bam}" "${sc_bam}"
samtools index "${sc_bam}"

# add readgroups to header for GATK
echo "adding GATK required readgroups for cell: ${barcode}"
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
      I="${sc_bam}" \
      O="${GATK_bam}" \
      SM="${barcode}" \
	  RGPL=illumina \
	  RGLB=lib1 \
	  RGPU=unit1

samtools index "${GATK_bam}"

module unload picard
module load GATK

# rm "${sc_bam}" "${sc_bam}.bai"

# variant Calling with GATK
echo "Variant Calling for cell: ${barcode}"

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R "${genome}" \
   -I "${GATK_bam}" \
   -O "${output_vcf}" \
   -ERC GVCF
   
gatk IndexFeatureFile -F "${output_vcf}"