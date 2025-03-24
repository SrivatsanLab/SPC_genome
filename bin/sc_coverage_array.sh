#!/bin/bash
#SBATCH -J coverage
#SBATCH -o array_outs/%A_%a.out
#SBATCH -p short
#SBATCH -t 02:00:00
#SBATCH -c 8

mkdir -p array_outs

module load fhPython
module load SAMtools
module load deepTools

# Input arguments
barcode_file="$1"
input_bam="$2"
temp_dir="$3"
chrom_sizes="$4"
scripts_dir="$5/scripts"

# Get the barcode for the current array task
barcode=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$barcode_file")
output_bam="${temp_dir}/${barcode}.bam"

echo "Processing barcode ${barcode} in SLURM task ${SLURM_ARRAY_TASK_ID}"

# Filter BAM for the barcode and save to temp directory

## For some reason -d isn't working since switching to atrandi combinatorial indices, idk why??
# samtools view -h -b -@ 8 -d CB:${barcode} "${input_bam}" > "${output_bam}"
## Workaround:
samtools view -h ${input_bam} | grep -E "^@|CB:Z:${barcode}" | samtools view -b -o "${output_bam}"

samtools index "${output_bam}"
echo "Filtered BAM for barcode ${barcode}"

# Generate BigWig for this barcode
bamCoverage -b "${output_bam}" -o "${temp_dir}/${barcode}.bw" --binSize 1 -of bigwig -p 8
bamCoverage -b "${output_bam}" -o "${temp_dir}/${barcode}_1kb.bw" --binSize 1000 -of bigwig -p 8
echo "created bigwigs for ${barcode}"

echo "running python script to compute coverage for ${barcode}"
python "${scripts_dir}/compute_coverage.py" "${barcode_file}" "${chrom_sizes}" "${temp_dir}"

echo "running python script to compute readcounts vs coverage in bins for lorenz curves"
python "${scripts_dir}/binned_coverage.py" "${barcode_file}" ref_genome "${temp_dir}" "${species}"

