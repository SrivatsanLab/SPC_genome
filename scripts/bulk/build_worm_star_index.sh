#!/bin/bash
#SBATCH -J build_STAR_worm
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 2:00:00

###########################################################################################################################
# Build STAR index for C. elegans GCA_028201515.1 (N2 strain, University of Oregon 2023)
# Downloads WormBase annotation and builds STAR index
###########################################################################################################################

set -euo pipefail

# Load STAR module
module load STAR/2.7.9a-GCC-11.2.0

# Paths
ORIGINAL_GENOME_FA="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/data/reference/worm_GCA_028201515.1/GCA_028201515.1_genomic.fna"
REFERENCE_DIR="/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/data/reference/worm_GCA_028201515.1_STAR"
RENAMED_GENOME_FA="${REFERENCE_DIR}/GCA_028201515.1_genomic_renamed.fna"
OUTPUT_DIR="${REFERENCE_DIR}/STAR_index"
GTF_FILE="${REFERENCE_DIR}/c_elegans.PRJNA13758.WS295.canonical_geneset.gtf"

# Create reference directory for STAR-compatible version
mkdir -p "${REFERENCE_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Rename chromosomes to match WormBase convention (I, II, III, etc.)
# Original uses RefSeq accessions: CP116366.1 -> I, CP116367.1 -> II, etc.
if [ ! -f "${RENAMED_GENOME_FA}" ]; then
    echo "Renaming chromosomes to match WormBase convention..."
    sed 's/>CP116366.1.*/>I Caenorhabditis elegans strain Bristol N2 chromosome I/' "${ORIGINAL_GENOME_FA}" | \
    sed 's/>CP116367.1.*/>II Caenorhabditis elegans strain Bristol N2 chromosome II/' | \
    sed 's/>CP116368.1.*/>III Caenorhabditis elegans strain Bristol N2 chromosome III/' | \
    sed 's/>CP116369.1.*/>IV Caenorhabditis elegans strain Bristol N2 chromosome IV/' | \
    sed 's/>CP116370.1.*/>V Caenorhabditis elegans strain Bristol N2 chromosome V/' | \
    sed 's/>CP116371.1.*/>X Caenorhabditis elegans strain Bristol N2 chromosome X/' \
    > "${RENAMED_GENOME_FA}"
    echo "Chromosome renaming complete: ${RENAMED_GENOME_FA}"
else
    echo "Renamed genome already exists: ${RENAMED_GENOME_FA}"
fi

# Use renamed genome for indexing
GENOME_FA="${RENAMED_GENOME_FA}"

# Use Ensembl WBcel235 GTF annotation (already available on cluster)
# Copy from shared biodata if not already present
ENSEMBL_GTF="/shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"
if [ ! -f "${GTF_FILE}" ]; then
    echo "Copying Ensembl WBcel235 GTF annotation..."
    cp "${ENSEMBL_GTF}" "${GTF_FILE}"
    echo "GTF copy complete"
else
    echo "GTF file already exists: ${GTF_FILE}"
fi

# Verify files exist
if [ ! -f "${GENOME_FA}" ]; then
    echo "ERROR: Genome FASTA not found: ${GENOME_FA}"
    exit 1
fi

if [ ! -f "${GTF_FILE}" ]; then
    echo "ERROR: GTF file not found: ${GTF_FILE}"
    exit 1
fi

echo "============================================"
echo "Building STAR index for C. elegans"
echo "============================================"
echo "Genome FASTA: ${GENOME_FA}"
echo "GTF file: ${GTF_FILE}"
echo "Output directory: ${OUTPUT_DIR}"
echo "STAR version: $(STAR --version)"
echo ""

# Build STAR index
# C. elegans genome is ~100Mb, so using --genomeSAindexNbases 12
STAR --runMode genomeGenerate \
     --runThreadN 16 \
     --genomeDir "${OUTPUT_DIR}" \
     --genomeFastaFiles "${GENOME_FA}" \
     --sjdbGTFfile "${GTF_FILE}" \
     --sjdbOverhang 100 \
     --genomeSAindexNbases 12

echo ""
echo "============================================"
echo "STAR index build complete!"
echo "============================================"
echo "Index location: ${OUTPUT_DIR}"
echo ""
echo "To use this index with STAR alignment, specify:"
echo "  --genomeDir ${OUTPUT_DIR}"
