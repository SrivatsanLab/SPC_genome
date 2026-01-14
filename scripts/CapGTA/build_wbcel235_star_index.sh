#!/bin/bash
#SBATCH -J build_STAR_index
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 2:00:00

###########################################################################################################################
# Build STAR index for C. elegans WBcel235 using STAR 2.7.9a
# Uses existing genome FASTA and GTF from Ensembl WBcel235
###########################################################################################################################

set -euo pipefail

# Load STAR 2.7.9a
module load STAR/2.7.9a-GCC-11.2.0

# Paths
GENOME_FA="/shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa"
GTF_FILE="/shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"
OUTPUT_DIR="./data/reference/WBcel235_STAR279a"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "Building STAR index for C. elegans WBcel235..."
echo "Genome FASTA: ${GENOME_FA}"
echo "GTF file: ${GTF_FILE}"
echo "Output directory: ${OUTPUT_DIR}"
echo "STAR version: $(STAR --version)"

# Build STAR index
# C. elegans genome is ~100Mb, so using --genomeSAindexNbases 12
STAR --runMode genomeGenerate \
     --runThreadN 16 \
     --genomeDir "${OUTPUT_DIR}" \
     --genomeFastaFiles "${GENOME_FA}" \
     --sjdbGTFfile "${GTF_FILE}" \
     --sjdbOverhang 100 \
     --genomeSAindexNbases 12

echo "STAR index build complete!"
echo "Index location: ${OUTPUT_DIR}"
