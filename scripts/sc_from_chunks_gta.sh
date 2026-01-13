#!/bin/bash
#SBATCH -J sc_chunks_gta
#SBATCH -o SLURM_outs/%x_%A.out

##########################################################################################################################
# This job extracts DNA and RNA reads from single cells from each chunk, then compiles into single cell BAMs            #
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
    echo "Submitting array job for chunk: $chunk with $cell_count barcodes"

    check_job_limit

    # Submit an array job for processing all barcodes in this chunk
    # This will extract both DNA and RNA reads for each barcode
    array_ID=$(sbatch --parsable --array=1-"$cell_count" "${scripts_DIR}/scripts/extract_sc_array_gta.sh" "$chunk" "$TMP_dir" "$barcode_file")

    # Add the array job ID to the list of job IDs
    job_ids="$job_ids$array_ID,"

    echo "Submitted job $array_ID for chunk $chunk, Job ID: ${array_ID}"
done < "$chunk_indices"

# Remove the trailing comma
job_ids=${job_ids%,}

echo "Waiting for all array jobs to finish: $job_ids"

# Wait for all the submitted array jobs to finish
srun --dependency=afterok:$job_ids sleep 1

echo "All chunk extraction jobs completed."
echo ""
echo "Merging per-chunk BAMs into final per-cell BAMs..."

# Merge the chunked BAMs into final single-cell BAMs
FINAL_SC_OUTPUT="${TMP_dir}/sc_outputs_final"
mkdir -p "${FINAL_SC_OUTPUT}"

merge_job_id=$(sbatch --parsable --dependency=afterok:$job_ids "${scripts_DIR}/scripts/merge_sc_chunks_gta.sh" \
    "$barcode_file" \
    "${TMP_dir}/sc_outputs" \
    "${FINAL_SC_OUTPUT}")

echo "Merge job submitted (ID: ${merge_job_id}). Final BAMs will be in: ${FINAL_SC_OUTPUT}"
echo ""

# Copy BAMs to final data directory
DATA_OUTPUT="${scripts_DIR}/data/PolE_worm_pilot/sc_outputs"
mkdir -p "${DATA_OUTPUT}"

copy_job_id=$(sbatch --parsable --dependency=afterok:${merge_job_id} --wrap="cp -v ${FINAL_SC_OUTPUT}/*.bam ${DATA_OUTPUT}/" -J copy_sc_bams -o SLURM_outs/copy_sc_bams_%j.out)

echo "Copy job submitted (ID: ${copy_job_id}). BAMs will be copied to: ${DATA_OUTPUT}"
echo ""

# Submit count matrix job
GTF="/shared/biodata/ngs/Reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"
COUNT_OUTPUT="${scripts_DIR}/results/PolE_worm_pilot"

count_job_id=$(sbatch --parsable --dependency=afterok:${copy_job_id} "${scripts_DIR}/scripts/create_rna_count_matrix.sh" \
    "${DATA_OUTPUT}" \
    "${GTF}" \
    "${COUNT_OUTPUT}/rna_counts")

echo "Count matrix job submitted (ID: ${count_job_id}). Results will be in: ${COUNT_OUTPUT}/rna_counts"
