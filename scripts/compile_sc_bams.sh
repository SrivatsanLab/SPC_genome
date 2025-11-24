#!/bin/bash
#SBATCH -J compile_sc
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH -t 2:00:00

##########################################################################################################################
# This job compiles per-chunk single-cell BAMs into final single-cell BAMs with QC metrics
# (No variant calling - QC only)
##########################################################################################################################

module load SAMtools

TMP_dir="$1"
barcode_file="$2"
SC_outputs="$3"
mkdir -p "${SC_outputs}"

# Get the barcode for the current array task
barcode=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$barcode_file")
sc_bam="${SC_outputs}/${barcode}.bam"

# Merge chunks into single cell bam:
echo "Merging chunks for cell: ${barcode}"

ls ${TMP_dir}/*${barcode}.bam > "${TMP_dir}/${barcode}_files.txt"

samtools merge -@ 4 -o "${sc_bam}" -b "${TMP_dir}/${barcode}_files.txt"
samtools sort -@ 4 -o "${sc_bam%.bam}.sorted.bam" "${sc_bam}"
mv "${sc_bam%.bam}.sorted.bam" "${sc_bam}"
samtools index -@ 4 "${sc_bam}"

# Generate QC metrics
echo "Generating QC metrics for cell: ${barcode}"
samtools flagstat "${sc_bam}" > "${SC_outputs}/${barcode}_flagstat.txt"
samtools stats "${sc_bam}" > "${SC_outputs}/${barcode}_stats.txt"

# Clean up temporary per-chunk BAMs for this cell
rm ${TMP_dir}/*${barcode}.bam
rm "${TMP_dir}/${barcode}_files.txt"

echo "Completed processing cell: ${barcode}"
echo "  BAM: ${sc_bam}"
echo "  Flagstat: ${SC_outputs}/${barcode}_flagstat.txt"
echo "  Stats: ${SC_outputs}/${barcode}_stats.txt"
