#!/bin/bash
#SBATCH -J concat_sams
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 4

###########################################################################################################################
# This job concatenates the processed chunks, and converts to sorted and indexed bam format. Next it computes read counts #
# per barcode to distinguish real cells from empty SPCs.                                                                  #
###########################################################################################################################

# Activate conda environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate spc_genome

OUTPUT_NAME="$1"
TMP_dir="$2"
ALIGNED_DIR="$3"  # Directory for aligned BAM/SAM files (data/[experiment]/aligned/)
BIN_DIR="$4"      # Directory for metadata files (bin/[experiment]/)
SCRIPTS_DIR="$5"

ls $TMP_dir/*.sam > "${BIN_DIR}/sam_list.txt"

SAM_FILE="${ALIGNED_DIR}/${OUTPUT_NAME}.sam"
BAM_FILE="${ALIGNED_DIR}/${OUTPUT_NAME}.bam"

module load SAMtools

echo "Concatenating SAM Chunks"

samtools merge -@ 4 -o "${SAM_FILE}" -b "${BIN_DIR}/sam_list.txt"

echo "Converting SAM to sorted BAM and indexing..."
samtools view -@ 4 -bS "${SAM_FILE}" | samtools sort -@ 4 -o "${BAM_FILE}" && samtools index -@ 4 "${BAM_FILE}"

echo "Creating knee plot and detecting real cells..."

samtools view -@ 4 "${BAM_FILE}" | python "${SCRIPTS_DIR}/scripts/readcounts.py" -o "${BIN_DIR}/readcounts.csv"

cat "${BIN_DIR}/readcounts.csv" | python "${SCRIPTS_DIR}/scripts/detect_cells.py" --plot "${BIN_DIR}/kneeplot.png" > "${BIN_DIR}/real_cells.txt"
