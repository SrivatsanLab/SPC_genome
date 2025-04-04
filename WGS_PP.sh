#!/bin/bash
#sbatch --job-name=Preprocessing
#sbatch --output=SLURM_outs/%x_%j.out
#sbatch -c 2

mkdir -p SLURM_outs/

# Set Defaults
SCRIPTS_DIR="./SPC_genome/"
OUTPUT_DIR="."
N_CHUNKS=500
TMP_DIR=$(mktemp -d)


mkdir -p "${OUTPUT_DIR}"

# Function to display help message
show_help() {
    cat << EOF
Usage: $0 -o <OUTPUT_NAME> -1 <read1.fastq.gz> -2 <read2.fastq.gz> -g <reference_genome.fa> -r <read_count>

Required arguments:
  -o    <output_name>           Desired sample name (prefix for outputs)
  -1    <read1.fastq.gz>        Read 1 FASTQ file
  -2    <read2.fastq.gz>        Read 2 FASTQ file
  -g    <reference_genome.fa>   Path to a **BWA INDEXED** reference genome fasta (for alignment)
  -r    <READ_COUNT>			Number of reads (from sequencing run info)
  
Optional arguments:
  -s    <scripts_DIR> 			path to the SPC_genome directory (default: ./SPC_genome)
  -O    <output_dir>            desired output directory (default: .)
  -n	<N_CHUNKS>				Number of subjobs for SLURM arrays. Default: 500
  -t    <TMP_DIR>				Temp directory for fastq chunks. Provide in order to override mktemp (e.g. when scratch space is limited)
  -h                            Show this help message and exit
EOF
}

# Parse options using getopts
while getopts ":o:1:2:g:s:O:n:t:h" option; do
  case $option in
    o) OUTPUT_NAME=$OPTARG ;;
    1) READ1=$OPTARG ;;
    2) READ2=$OPTARG ;;
    g) REFERENCE_GENOME=$OPTARG ;;
	r) READ_COUNT=$OPTARG ;;
    s) scripts_DIR=$OPTARG ;;  # Optional argument with default "./SPC_scWGS"
	O) OUTPUT_DIR=$OPTARG ;; # Optional argument with default "."
	n) N_CHUNKS=$OPTARG ;; # Optional argument with default 5000
	t) TMP_DIR=$OPTARG ;; # Optional argument, provide to override defualt `mktemp`
    h) show_help; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; show_help; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; show_help; exit 1 ;;
  esac
done

# Check that required arguments are provided
if [ -z "$OUTPUT_NAME" ] || [ -z "$READ1" ] || [ -z "$READ2" ] || [ -z "$REFERENCE_GENOME" ] || [ -z "$READ_COUNT" ]; then
    echo "Missing required arguments." >&2
    show_help
    exit 1
fi

# Print inputs
echo "processing $READ_COUNT reads with input parameters:"
echo "Sample Name: ${OUTPUT_NAME}"
echo "Read 1: ${READ1}"
echo "Read 2: ${READ2}"
echo "Reference Genome: ${REFERENCE_GENOME}"

echo "Scripts Directory: ${SCRIPTS_DIR}"
echo "Barcode Directory: ${BARCODE_DIR}"
echo "Output Directory: ${OUTPUT_DIR}"

######################################################################################################
#### chunk the fastqs

echo "Splitting fastq's into chunks for parallel processing, this might take a while..."

# Compute total lines (each read has 4 lines)
total_lines=$(( READ_COUNT * 4 ))

# Compute lines per chunk
CHUNK_LINES=$(( total_lines / N_CHUNKS ))

# Ensure CHUNK_LINES is a multiple of 4
CHUNK_LINES=$(( CHUNK_LINES / 4 * 4 ))

# Step 3: Split each FASTQ file into chunks based on respective chunk size
zcat "$READ1" | split -l $CHUNK_LINES - "$TMP_DIR/read1_chunk_" &
zcat "$READ2" | split -l $CHUNK_LINES - "$TMP_DIR/read2_chunk_" &
wait

ls "$TMP_DIR"/read1_chunk_* | sed 's/.*chunk_//' > "${OUTPUT_DIR}/chunk_indices.txt"

chunk_count=$(wc -l "${OUTPUT_DIR}/chunk_indices.txt")

######################################################################################################
#### Submit First Job array

PP_array_ID=$(sbatch --parsable --array=1-$chunk_count "${scripts_DIR}/PP_array.sh" "${OUTPUT_DIR}/chunk_indices.txt" "${REFERENCE_GENOME}" "${scripts_DIR}" "${TMP_DIR}")

echo "Preprocessing array job ID: ${PP_array_ID}"

######################################################################################################
#### Make Knee Plot and detect real cells

knee_plot_array_ID=$()

######################################################################################################
#### Submit single cell extraction arrays

sc_from_chunks_array_ID=$(sbatch --parsable --dependency=afterok:$knee_plot_array_ID "${scripts_DIR}/sc_from_chunks.sh" "${OUTPUT_DIR}/chunk_indices.txt" "${TMP_DIR}" "${OUTPUT_DIR}/real_cells.txt" "${SCRIPTS_DIR}")

echo "Compiling single cell bam files with array job ID: ${sc_array_ID}"

######################################################################################################
#### Submit single cell variant calling array
cell_count=$(wc -l "${OUTPUT_DIR}/real_cells.txt")
mkdir -p "${OUTPUT_DIR}/sc_ouputs"

sc_var_array_ID=$(sbatch --parsable --array=1-$cell_count --dependency=afterok:$sc_comp_array_ID "${scripts_DIR}/sc_var_array.sh" "${TMP_DIR}" "${OUTPUT_DIR}/real_cells.txt" "${REFERENCE_GENOME}" "${OUTPUT_DIR}/sc_ouputs")

echo "Compiling single cell bam files with array job ID: ${sc_array_ID}"

######################################################################################################
#### Submit matched bulk variant calling job

# bulk_v_job_ID=$(sbatch --parsable --dependency=afterok:$comp_job_ID )

# echo "Bulk calling job array ID: ${bulk_v_job_ID}"

######################################################################################################
#### Submit joint calling job

# ls "${OUTPUT_DIR}/sc_output/*.bam" > bam_list.txt

# joint_calling_job_id=$(sbatch --parsable --dependency=afterok:$--dependency=afterok:joint_calling_job_id)

# echo "Joint Calling Job ID: ${joint_calling_job_id}"

######################################################################################################
#### Submit anndata job