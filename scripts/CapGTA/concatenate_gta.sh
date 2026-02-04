#!/bin/bash
#SBATCH -J concat_gta
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 8
#SBATCH --mem=64G

###########################################################################################################################
# This job concatenates DNA and RNA SAM chunks separately, creates sorted BAMs, and detects real cells                  #
# using combined DNA+RNA read counts per barcode.                                                                        #
###########################################################################################################################

# Activate conda environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate spc_genome

OUTPUT_NAME="$1"
TMP_dir="$2"
ALIGNED_DIR="$3"  # Directory for aligned BAM/SAM files (data/[experiment]/aligned/)
BIN_DIR="$4"      # Directory for metadata files (bin/[experiment]/)
SCRIPTS_DIR="$5"

# Extract just the sample name (basename) for file naming
SAMPLE_NAME=$(basename "${OUTPUT_NAME}")

module load SAMtools

###########################################################################################################################
# Concatenate DNA SAM chunks
###########################################################################################################################

echo "Listing DNA SAM chunks..."
ls $TMP_dir/*_dna.sam > "${BIN_DIR}/sam_list_dna.txt"

DNA_BAM="${ALIGNED_DIR}/${SAMPLE_NAME}_dna.bam"

echo "Merging DNA SAM chunks directly to sorted BAM..."
samtools merge -@ 8 -b "${BIN_DIR}/sam_list_dna.txt" -O SAM - | samtools sort -@ 8 -o "${DNA_BAM}"

echo "Indexing DNA BAM file..."
samtools index -@ 8 "${DNA_BAM}"

echo "DNA BAM created: ${DNA_BAM}"

###########################################################################################################################
# Concatenate RNA SAM chunks
###########################################################################################################################

echo "Listing RNA SAM chunks..."
ls $TMP_dir/*_rna.sam > "${BIN_DIR}/sam_list_rna.txt"

RNA_BAM="${ALIGNED_DIR}/${SAMPLE_NAME}_rna.bam"

echo "Merging RNA SAM chunks directly to sorted BAM..."
samtools merge -@ 8 -b "${BIN_DIR}/sam_list_rna.txt" -O SAM - | samtools sort -@ 8 -o "${RNA_BAM}"

echo "Indexing RNA BAM file..."
samtools index -@ 8 "${RNA_BAM}"

echo "RNA BAM created: ${RNA_BAM}"

###########################################################################################################################
# Count reads per barcode from BOTH DNA and RNA BAMs (combined)
###########################################################################################################################

echo "Counting reads per barcode from DNA BAM..."
samtools view -@ 8 "${DNA_BAM}" | python "${SCRIPTS_DIR}/scripts/utils/readcounts.py" -o "${BIN_DIR}/readcounts_dna.csv"

echo "Counting reads per barcode from RNA BAM..."
samtools view -@ 8 "${RNA_BAM}" | python "${SCRIPTS_DIR}/scripts/utils/readcounts.py" -o "${BIN_DIR}/readcounts_rna.csv"

echo "Combining DNA and RNA read counts..."
python - <<EOF
import pandas as pd
import sys

# Read DNA and RNA read counts
dna = pd.read_csv("${BIN_DIR}/readcounts_dna.csv", index_col=0)
rna = pd.read_csv("${BIN_DIR}/readcounts_rna.csv", index_col=0)

# Combine counts (sum for shared barcodes, keep unique barcodes)
combined = dna.add(rna, fill_value=0).astype(int)

# Sort by read count descending
combined = combined.sort_values(by=combined.columns[0], ascending=False)

# Save combined counts
combined.to_csv("${BIN_DIR}/readcounts.csv")

print(f"Combined read counts saved to ${BIN_DIR}/readcounts.csv")
print(f"Total DNA reads: {dna.sum().values[0]}")
print(f"Total RNA reads: {rna.sum().values[0]}")
print(f"Total combined reads: {combined.sum().values[0]}")
EOF

###########################################################################################################################
# Detect real cells using combined counts
###########################################################################################################################

echo "Creating knee plot and detecting real cells using combined DNA+RNA counts..."

cat "${BIN_DIR}/readcounts.csv" | python "${SCRIPTS_DIR}/scripts/CapWGS/detect_cells.py" --plot "${BIN_DIR}/kneeplot.png" > "${BIN_DIR}/real_cells.txt"

cell_count=$(wc -l < "${BIN_DIR}/real_cells.txt")
echo "Detected ${cell_count} real cells"

echo "Concatenation and cell detection complete!"
