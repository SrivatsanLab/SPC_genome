#!/bin/bash
#SBATCH -J create_adata
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 2:00:00

###########################################################################################################################
# Create AnnData objects for variants and RNA counts
###########################################################################################################################

set -euo pipefail

merged_vcf="$1"
rna_count_matrix="$2"
output_dir="$3"
scripts_dir="$4"

echo "Creating AnnData objects..."
echo "Merged VCF: ${merged_vcf}"
echo "RNA count matrix: ${rna_count_matrix}"
echo "Output directory: ${output_dir}"
echo ""

# Create output directory
mkdir -p "${output_dir}"

# Output files
variant_h5ad="${output_dir}/variants.h5ad"
rna_h5ad="${output_dir}/rna_counts.h5ad"

# Activate micromamba environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

# Create variant AnnData
echo "Creating variant AnnData from merged VCF..."
python "${scripts_dir}/scripts/make_variant_anndata.py" \
    "${merged_vcf}" \
    "${variant_h5ad}"

echo ""
echo "Variant AnnData created: ${variant_h5ad}"
echo ""

# Create RNA AnnData
echo "Creating RNA AnnData from count matrix..."
python "${scripts_dir}/scripts/rna_matrix_to_anndata.py" \
    "${rna_count_matrix}" \
    "${rna_h5ad}"

echo ""
echo "RNA AnnData created: ${rna_h5ad}"
echo ""

echo "================================"
echo "AnnData creation complete!"
echo "================================"
echo "Outputs:"
echo "  Variants: ${variant_h5ad}"
echo "  RNA counts: ${rna_h5ad}"
