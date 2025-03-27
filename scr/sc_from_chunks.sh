#!/bin/bash
#SBATCH -J sc_chunks
#SBATCH -o SLURM_outs/%x_%A.out

##########################################################################################################################
# This job extracts reads from single cells from each chunk, then compiles into single cell bams                         #
##########################################################################################################################

chunk_indices="$1"
TMP_dir="$2"
barcode_file="$3"
scripts_DIR="$4"

# Create an empty string to store job IDs
job_ids=""


# Count barcodes
cell_count=$(wc -l < "$barcode_file")

# Iterate over each chunk
while IFS= read -r chunk; do
    sam_file="${TMP_dir}/${chunk}.sam"
    

    echo "Submitting array job for chunk: $chunk with $cell_count barcodes"
    
    # Submit an array job for processing all barcodes in this chunk
    array_ID=$(sbatch --parsable --array=1-"$cell_count" "${scripts_DIR}/scr/extract_sc_array.sh" "$sam_file" "$barcode_file")
    
    # Add the array job ID to the list of job IDs
    job_ids="$job_ids$array_ID,"
    
    echo "Submitted job $array_ID for chunk $chunk"
done < "$chunk_indices"

# Remove the trailing comma
job_ids=${job_ids%,}

echo "Waiting for all array jobs to finish: $job_ids"

# Wait for all the submitted array jobs to finish
srun --dependency=afterok:$job_ids sleep 1

echo "All chunk jobs completed."