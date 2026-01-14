#!/bin/bash
#SBATCH -J merge_vcfs
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 4:00:00

###########################################################################################################################
# Merge single-cell VCFs into multi-sample VCF for anndata creation
###########################################################################################################################

set -euo pipefail

sc_vcf_dir="$1"
output_vcf="$2"

# Load BCFtools
module load BCFtools

echo "Merging single-cell VCFs..."
echo "Input directory: ${sc_vcf_dir}"
echo "Output VCF: ${output_vcf}"
echo ""

# Get list of all VCF files
vcf_files=(${sc_vcf_dir}/*_variants.vcf.gz)

echo "Found ${#vcf_files[@]} VCF files"
echo ""

# Create temporary directory for reheadered VCFs
tmp_dir="${sc_vcf_dir}/tmp_reheadered"
mkdir -p "${tmp_dir}"

echo "Reheadering VCFs with barcode sample names..."
reheadered_vcfs=()

for vcf in "${vcf_files[@]}"; do
    # Extract barcode from filename
    filename=$(basename "$vcf")
    barcode="${filename%_variants.vcf.gz}"

    # Create reheadered VCF with barcode as sample name
    tmp_vcf="${tmp_dir}/${filename}"
    echo "${barcode}" | bcftools reheader -s /dev/stdin -o "${tmp_vcf}" "${vcf}"
    bcftools index "${tmp_vcf}"

    reheadered_vcfs+=("${tmp_vcf}")
done

echo "Reheadering complete!"
echo ""

# Create file list for bcftools merge
vcf_list="${output_vcf%.vcf.gz}_file_list.txt"
printf "%s\n" "${reheadered_vcfs[@]}" > "${vcf_list}"

echo "Merging VCFs..."
bcftools merge \
    --threads 8 \
    --merge all \
    -O z \
    -o "${output_vcf}" \
    --file-list "${vcf_list}"

# Clean up temporary reheadered VCFs
echo "Cleaning up temporary files..."
rm -rf "${tmp_dir}"

echo "Indexing merged VCF..."
bcftools index "${output_vcf}"

# Clean up file list
rm "${vcf_list}"

# Summary
echo ""
echo "=== Merged VCF Summary ==="
bcftools view -h "${output_vcf}" | grep "^#CHROM" | tr '\t' '\n' | tail -n +10 | wc -l | awk '{print "Total samples: " $1}'
bcftools view -H "${output_vcf}" | wc -l | awk '{print "Total variants: " $1}'
bcftools view -H -v snps "${output_vcf}" | wc -l | awk '{print "SNPs: " $1}'
bcftools view -H -v indels "${output_vcf}" | wc -l | awk '{print "Indels: " $1}'

echo ""
echo "VCF merging complete!"
echo "Output: ${output_vcf}"

module unload BCFtools
