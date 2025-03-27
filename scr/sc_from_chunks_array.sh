#!/bin/bash
#SBATCH -J sc_array
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out

##########################################################################################################################
# This job extracts reads from single cells from each chunk, then compiles into single cell bams                         #
##########################################################################################################################

chunk_indices="$1"
TMP_dir="$2"
barcode_file="$4"
scripts_DIR="$4"

chunk=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$chunk_indices")

sam_file="${TMP_DIR}/${chunk}.sam"

cell_count=$(wc -l "${barcode_file}")

array_ID=$(sbatch --parsable --array:1-$cell_count "${scripts_DIR}/scr/extract_sc_array.sh" "${sam_file}" "${barcode_file}")

# Wait for the subjobs to finish
srun --dependency=afterok:$array_ID sleep 1