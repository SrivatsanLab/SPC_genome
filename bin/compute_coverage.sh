#!/bin/bash
#SBATCH --job-name=compute_cov_master

# TODO: 
# How to handle alt chromosomes

# Check command-line arguments
if [ "$#" -lt 2 ] || [ "$#" -gt 5 ]; then
    echo "Usage: $0 <list_of_cells.txt> <input_bam.bam> [OUTPUT_DIR (default: .)] [path/to/sci_Genome_PP (default: ./sci_Genome_PP)] [path/to/ref_genome_chrom_sizes (default: ./hg38/hg38.chrom.sizes)]"
    exit 1
fi

CELLS="$1"
INPUT_BAM="$2"
OUTPUT_DIR="${3:-.}"
sci_Genome_PP="${4:-./sci_Genome_PP}"
chrom_sizes="${5:-./hg38/hg38.chrom.sizes}"

echo "setting job array to compute coverage- cells : ${CELLS}"
echo "input bam : ${INPUT_BAM}"
echo "output dir : ${OUTPUT_DIR}"
echo "path/to/sci_Genom_PP : ${sci_Genome_PP}"
echo "path/to/ref_genome_chrom_sizes : ${chrom_sizes}"

# Create a unique temporary directory
# temp_dir=$(mktemp -d -t jobarray-XXXXXX)
# chmod a+rwx $temp_dir
# echo "directory for cell specific outputs created at $temp_dir"

mkdir -p "${OUTPUT_DIR}/sc_output"
temp_dir="${OUTPUT_DIR}/sc_output"
echo "directory for cell specific outputs created at $temp_dir"

# Get the number of lines (barcodes) in the barcode file
num_barcodes=$(wc -l < "$CELLS")
echo "Submitting job array with $num_barcodes tasks, each processing one barcode."

# Submit job array, passing the temp directory and other variables
array_job_id=$(sbatch --array=1-"$num_barcodes" "${sci_Genome_PP}/scripts/sc_coverage_array.sh" "${CELLS}" "${INPUT_BAM}" "${temp_dir}" "${chrom_sizes}" "${sci_Genome_PP}" | awk '{print $4}')
echo "Job array submitted with job ID $array_job_id"

# wait for array to finish, concatenate results from temp dir
echo "barcode,coverage" > "${OUTPUT_DIR}/coverages.csv"
cat_job_id=$(sbatch --dependency=afterok:${array_job_id} --wrap="cat \"${temp_dir}/*.csv\" >> \"${OUTPUT_DIR}/coverages.csv\"")
echo "results will be concatenated into ${OUTPUT_DIR}/coverages.csv"

# Submit python job that depends on the completion of the array to plot the distribution of cell coverages and write to the ouput_dir
# sbatch --dependency=afterok:${cat_job_id} "${sci_Genome_PP}/bin/python_batch.sh ${sci_Genome_PP}/scripts/coverage_dist.py ${CELLS}"
# echo "Python plotting job submitted to run after arry completes"


# Submit a cleanup job that depends on the completion of the python job
#sbatch --dependency=afterok:$array_job_id --wrap="rm -rf '$temp_dir'"
#echo "Cleanup job submitted to run after array job $python_job_id completes"

