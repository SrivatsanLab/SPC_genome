#!/bin/bash
#SBATCH -J concat_sams
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 4

###########################################################################################################################
# This job concatenates the processed chunks, and converts to sorted and indexed bam format. Next it computes read counts #
# per barcode to distinguish real cells from empty SPCs.                                                                  #
###########################################################################################################################

OUTPUT_NAME="$1"
TMP_dir="$2"
OUTPUT_DIR="$3"
SCRIPTS_DIR="$4"

ls $TMP_dir/*.sam > "${OUTPUT_DIR}/sam_list.txt"

SAM_FILE="${OUTPUT_DIR}/${OUTPUT_NAME}.sam"
BAM_FILE="${OUTPUT_DIR}/${OUTPUT_NAME}.bam"

module load SAMtools

echo "Concatenating SAM Chunks"

samtools merge -@ 4 -o "${SAM_FILE}" -b "${OUTPUT_DIR}/sam_list.txt"

echo "Converting SAM to sorted BAM and indexing..."
samtools view -@ 4 -bS "${SAM_FILE}" | samtools sort -@ 4 -o "${BAM_FILE}" && samtools index -@ 4 "${BAM_FILE}"

echo "Creating knee plot and detecting real cells..."

samtools view -@ 4 "${BAM_FILE}" | python "${SCRIPTS_DIR}/bin/readcounts.py" -o "${OUTPUT_DIR}/readcounts.csv"

cat "${OUTPUT_DIR}/readcounts.csv" | python "${SCRIPTS_DIR}/bin/detect_cells.py" --plot "${OUTPUT_DIR}/kneeplot.png" > "${OUTPUT_DIR}/real_cells.txt"
