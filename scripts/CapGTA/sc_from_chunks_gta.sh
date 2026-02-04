#!/bin/bash
#SBATCH -J sc_chunks_gta
#SBATCH -o SLURM_outs/%x_%A.out

##########################################################################################################################
# This job extracts DNA and RNA reads from single cells from each chunk, then compiles into single cell BAMs            #
# Then performs variant calling and RNA counting                                                                        #
##########################################################################################################################

chunk_indices="$1"
TMP_dir="$2"
barcode_file="$3"
scripts_DIR="$4"
SC_OUTPUTS_DIR="$5"
RESULTS_DIR="$6"
OUTPUT_NAME="$7"

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
    array_ID=$(sbatch --parsable --array=1-"$cell_count" "${scripts_DIR}/scripts/CapGTA/extract_sc_array_gta.sh" "$chunk" "$TMP_dir" "$barcode_file")

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

# Merge the chunked BAMs into final single-cell BAMs and copy to sc_outputs directory
FINAL_SC_OUTPUT="${TMP_dir}/sc_outputs_final"
mkdir -p "${FINAL_SC_OUTPUT}"
mkdir -p "${SC_OUTPUTS_DIR}"

merge_job_id=$(sbatch --parsable --dependency=afterok:$job_ids "${scripts_DIR}/scripts/CapGTA/merge_sc_chunks_gta.sh" \
    "$barcode_file" \
    "${TMP_dir}/sc_outputs" \
    "${FINAL_SC_OUTPUT}")

echo "Merge job submitted (ID: ${merge_job_id}). Final BAMs will be in: ${FINAL_SC_OUTPUT}"
echo ""

# Copy BAMs to final data directory
copy_job_id=$(sbatch --parsable --dependency=afterok:${merge_job_id} --wrap="cp -v ${FINAL_SC_OUTPUT}/*.bam ${SC_OUTPUTS_DIR}/" -J copy_sc_bams -o SLURM_outs/copy_sc_bams_%j.out)

echo "Copy job submitted (ID: ${copy_job_id}). BAMs will be copied to: ${SC_OUTPUTS_DIR}"
echo ""

# Submit count matrix job
GTF="/shared/biodata/ngs/Reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"

count_job_id=$(sbatch --parsable --dependency=afterok:${copy_job_id} "${scripts_DIR}/scripts/CapGTA/create_rna_count_matrix.sh" \
    "${SC_OUTPUTS_DIR}" \
    "${GTF}" \
    "${RESULTS_DIR}/rna_counts")

echo "Count matrix job submitted (ID: ${count_job_id}). Results will be in: ${RESULTS_DIR}/rna_counts"
echo ""

# Submit single-cell variant calling (BCFtools)
REFERENCE_FA="/shared/biodata/ngs/Reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa"

vcall_job_id=$(sbatch --parsable --dependency=afterok:${copy_job_id} --array=1-${cell_count} "${scripts_DIR}/scripts/CapGTA/sc_variant_calling_bcftools_array.sh" \
    "${barcode_file}" \
    "${SC_OUTPUTS_DIR}" \
    "${REFERENCE_FA}" \
    "${SC_OUTPUTS_DIR}")

echo "Single-cell variant calling job submitted (ID: ${vcall_job_id}). VCFs will be in: ${SC_OUTPUTS_DIR}"
echo ""

# Merge single-cell VCFs
MERGED_VCF="${RESULTS_DIR}/sc_variants_merged.vcf.gz"

merge_vcf_job_id=$(sbatch --parsable --dependency=afterok:${vcall_job_id} "${scripts_DIR}/scripts/CapGTA/merge_sc_vcfs.sh" \
    "${SC_OUTPUTS_DIR}" \
    "${MERGED_VCF}")

echo "VCF merge job submitted (ID: ${merge_vcf_job_id}). Merged VCF will be: ${MERGED_VCF}"
echo ""
echo "================================"
echo "Pipeline submission complete!"
echo "================================"
echo "Final outputs:"
echo "  - RNA count matrix: ${RESULTS_DIR}/rna_counts_matrix.csv"
echo "  - Merged VCF: ${MERGED_VCF}"
