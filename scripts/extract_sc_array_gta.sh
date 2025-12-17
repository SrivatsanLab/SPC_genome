#!/bin/bash
#SBATCH -J extract_sc_gta
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -t 12:00:00

##########################################################################################################################
# This job extracts DNA and RNA reads from single cells from chunked SAM files                                          #
##########################################################################################################################

module load SAMtools

chunk="$1"        # Chunk ID (e.g., "00", "01")
TMP_dir="$2"      # Temporary directory containing SAM files
barcode_file="$3" # File containing barcodes (one per line)

# Get the barcode for the current array task
barcode=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$barcode_file")

# Input files
dna_sam="${TMP_dir}/${chunk}_dna.sam"
rna_sam="${TMP_dir}/${chunk}_rna.sam"

# Create sc_outputs subdirectory if it doesn't exist
sc_output_dir="${TMP_dir}/sc_outputs"
mkdir -p "${sc_output_dir}"

# Output files
dna_bam="${sc_output_dir}/${barcode}_dna.bam"
rna_bam="${sc_output_dir}/${barcode}_rna.bam"

# Extract DNA reads for this barcode
if [ -f "${dna_sam}" ]; then
    samtools view -h "${dna_sam}" | grep -E "^@|CB:Z:${barcode}" | samtools view -b -o "${dna_bam}"
else
    echo "Warning: DNA SAM file not found: ${dna_sam}"
fi

# Extract RNA reads for this barcode
if [ -f "${rna_sam}" ]; then
    samtools view -h "${rna_sam}" | grep -E "^@|CB:Z:${barcode}" | samtools view -b -o "${rna_bam}"
else
    echo "Warning: RNA SAM file not found: ${rna_sam}"
fi

echo "Extracted reads for barcode ${barcode}:"
echo "  DNA: ${dna_bam}"
echo "  RNA: ${rna_bam}"
