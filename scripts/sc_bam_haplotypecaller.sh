#!/bin/bash
#SBATCH -J sc_hc
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH -t 12:00:00

##########################################################################################################################
# Variant calling on already-extracted single-cell BAM files using GATK HaplotypeCaller                                #
# Use this when you have sc BAMs extracted from bulk (not chunks to merge)                                            #
##########################################################################################################################

set -euo pipefail

SC_BAM_DIR="$1"
BARCODE_FILE="$2"
REFERENCE_GENOME="$3"
OUTPUT_DIR="$4"
SAMPLE_PREFIX="$5"

mkdir -p "${OUTPUT_DIR}"

# Get the barcode for the current array task
barcode=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BARCODE_FILE")

# Input and output paths
sc_bam="${SC_BAM_DIR}/${SAMPLE_PREFIX}_${barcode}.bam"
GATK_bam="${OUTPUT_DIR}/${barcode}_GATK.bam"
output_vcf="${OUTPUT_DIR}/${barcode}.g.vcf.gz"

# Check if input BAM exists
if [ ! -f "${sc_bam}" ]; then
    echo "ERROR: Single-cell BAM not found: ${sc_bam}"
    exit 1
fi

echo "Processing cell: ${barcode}"

# Load modules
module load SAMtools picard

# Add read groups to header for GATK
echo "Adding GATK required readgroups for cell: ${barcode}"
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
      I="${sc_bam}" \
      O="${GATK_bam}" \
      SM="${barcode}" \
      RGID="${barcode}" \
      RGPL=illumina \
      RGLB=lib1 \
      RGPU=unit1

samtools index "${GATK_bam}"

module unload SAMtools picard
module load GATK

# Variant Calling with GATK
echo "Variant calling for cell: ${barcode}"

# Construct reference path
REFERENCE="${REFERENCE_GENOME}/Homo_sapiens_assembly38.fasta"

gatk --java-options "-Xmx16g" HaplotypeCaller  \
   -R "${REFERENCE}" \
   -I "${GATK_bam}" \
   -O "${output_vcf}" \
   -ERC GVCF

echo "Completed variant calling for cell: ${barcode}"
echo "Output: ${output_vcf}"
