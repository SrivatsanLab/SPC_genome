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

max_jobs=50000 #single user limit on rhino

# Function to check current job count
check_job_limit() {
    # Count the number of jobs the user has in the queue (running, pending, etc.)
    job_count=$(squeue --array -u $USER | grep -v "squeue" | wc -l)
    echo "cumulative total jobs: $job_count"
    job_count=$((job_count + cell_count))

    # If the number of jobs exceeds the max limit, wait until there is space
    if [ "$job_count" -ge "$max_jobs" ]; then
        echo "Max job limit reached ($job_count jobs). Waiting for jobs to finish..."
        while [ "$job_count" -ge "$max_jobs" ]; do
            # Check job count again and wait if needed
            job_count=$(squeue -u $USER | grep -v "squeue" | wc -l)
			job_count=$((job_count + cell_count))
            sleep 60  # Wait for 1 minute before checking again
        done
        echo "Space available. Continuing job submissions."
    fi
}

# Iterate over each chunk
while IFS= read -r chunk; do
    bam_file="${TMP_dir}/${chunk}.bam"

    echo "Submitting array job for chunk: $chunk with $cell_count barcodes"

    check_job_limit

    # Submit an array job for processing all barcodes in this chunk
    array_ID=$(sbatch --parsable --array=1-"$cell_count" "${scripts_DIR}/scr/extract_sc_array.sh" "$bam_file" "$barcode_file")

    # Add the array job ID to the list of job IDs
    job_ids="$job_ids$array_ID,"

    echo "Submitted job $array_ID for chunk $chunk, Job ID: ${array_ID}"
done < "$chunk_indices"

# Remove the trailing comma
job_ids=${job_ids%,}

echo "Waiting for all array jobs to finish: $job_ids"

# Wait for all the submitted array jobs to finish
srun --dependency=afterok:$job_ids sleep 1

echo "All chunk jobs completed."