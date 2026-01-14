#!/bin/bash
#SBATCH -J rna_counts
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 4:00:00

###########################################################################################################################
# Create RNA count matrix (genes x cells) from single-cell RNA BAMs
# Uses featureCounts to count reads per gene per cell
###########################################################################################################################

set -euo pipefail

SC_RNA_DIR="$1"       # Directory containing *_rna.bam files
GTF="$2"              # GTF annotation file
OUTPUT_PREFIX="$3"    # Output prefix (e.g., results/PolE_worm_pilot/rna_counts)

mkdir -p "$(dirname ${OUTPUT_PREFIX})"

echo "Creating RNA count matrix..."
echo "SC RNA BAMs: ${SC_RNA_DIR}"
echo "GTF: ${GTF}"
echo "Output prefix: ${OUTPUT_PREFIX}"
echo ""

# Load required modules
module load Subread  # Contains featureCounts

# Get list of all RNA BAMs
RNA_BAMS=(${SC_RNA_DIR}/*_rna.bam)
NUM_CELLS=${#RNA_BAMS[@]}

echo "Found ${NUM_CELLS} RNA BAMs"
echo ""

# Run featureCounts on all RNA BAMs
# -a: GTF file
# -o: output file
# -T: number of threads
# -t: feature type (exon)
# -g: group by gene_id
# -M: count multi-mapping reads
# --fracOverlap: minimum fraction of overlapping bases (0.5 = 50%)
# --primary: only count primary alignments
# -p: paired-end mode
# -B: only count read pairs with both ends mapped
# -C: do not count chimeric fragments

echo "Running featureCounts..."
featureCounts \
    -a "${GTF}" \
    -o "${OUTPUT_PREFIX}_raw.txt" \
    -T 16 \
    -t exon \
    -g gene_id \
    -M \
    --fracOverlap 0.5 \
    --primary \
    -p \
    -B \
    -C \
    "${RNA_BAMS[@]}"

echo ""
echo "featureCounts complete! Raw output: ${OUTPUT_PREFIX}_raw.txt"
echo ""

# Convert featureCounts output to clean matrix format
# featureCounts output has 6 metadata columns, then one column per BAM
# We'll extract just the gene IDs and counts, and clean up BAM file names

echo "Converting to clean count matrix format..."

python3 << PYTHON_SCRIPT
import pandas as pd
import os

# Read featureCounts output
input_file = "${OUTPUT_PREFIX}_raw.txt"
output_prefix = "${OUTPUT_PREFIX}"

# Read the count table (skip comment line starting with #)
counts = pd.read_csv(input_file, sep='\t', comment='#')

# Extract gene IDs and count columns
# First 6 columns are: Geneid, Chr, Start, End, Strand, Length
gene_ids = counts['Geneid']
count_cols = counts.columns[6:]  # Skip first 6 metadata columns

# Extract cell barcodes from BAM filenames
# Format: /path/to/BARCODE_rna.bam -> BARCODE
cell_barcodes = []
for col in count_cols:
    # Extract filename from path
    basename = os.path.basename(col)
    # Remove _rna.bam suffix
    barcode = basename.replace('_rna.bam', '')
    cell_barcodes.append(barcode)

# Create count matrix with clean names
count_matrix = counts[count_cols].copy()
count_matrix.columns = cell_barcodes
count_matrix.insert(0, 'gene_id', gene_ids)

# Save as CSV
output_csv = output_prefix + '_matrix.csv'
count_matrix.to_csv(output_csv, index=False)
print(f"Saved count matrix: {output_csv}")
print(f"Dimensions: {len(gene_ids)} genes x {len(cell_barcodes)} cells")

# Print summary statistics
total_counts = count_matrix[cell_barcodes].sum(axis=0)
genes_detected = (count_matrix[cell_barcodes] > 0).sum(axis=0)

print("\nPer-cell summary:")
print(f"  Mean counts per cell: {total_counts.mean():.0f}")
print(f"  Median counts per cell: {total_counts.median():.0f}")
print(f"  Mean genes detected per cell: {genes_detected.mean():.0f}")
print(f"  Median genes detected per cell: {genes_detected.median():.0f}")

# Also create a summary file
summary = pd.DataFrame({
    'cell_barcode': cell_barcodes,
    'total_rna_counts': total_counts.values,
    'genes_detected': genes_detected.values
})
summary_file = output_prefix + '_summary.csv'
summary.to_csv(summary_file, index=False)
print(f"\nSaved per-cell summary: {summary_file}")
PYTHON_SCRIPT

echo ""
echo "====================================="
echo "RNA count matrix generation complete!"
echo "====================================="
echo "Outputs:"
echo "  Raw counts: ${OUTPUT_PREFIX}_raw.txt"
echo "  Count matrix (CSV): ${OUTPUT_PREFIX}_matrix.csv"
echo "  Per-cell summary: ${OUTPUT_PREFIX}_summary.csv"
echo ""
