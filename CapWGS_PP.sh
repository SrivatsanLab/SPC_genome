#!/bin/bash
#sbatch --job-name=Preprocessing
#sbatch --output=SLURM_outs/%x_%j.out
#sbatch -c 2

mkdir -p SLURM_outs/

# Function to parse YAML config file
parse_yaml() {
    local yaml_file=$1
    local prefix=$2
    local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
    sed -ne "s|^\($s\):|\1|" \
         -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
         -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p" "$yaml_file" |
    awk -F$fs '{
        indent = length($1)/2;
        vname[indent] = $2;
        for (i in vname) {if (i > indent) {delete vname[i]}}
        if (length($3) > 0) {
            vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
            printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
        }
    }'
}

# Load defaults from config.yaml if it exists
CONFIG_FILE="./config.yaml"
if [ -f "$CONFIG_FILE" ]; then
    eval $(parse_yaml "$CONFIG_FILE" "CONFIG_")

    # Set defaults from config
    SCRIPTS_DIR="."  # Project root directory (contains bin/, scripts/, barcodes/)
    OUTPUT_DIR="${CONFIG_output_base_dir:-.}"
    N_CHUNKS="${CONFIG_processing_n_chunks:-500}"
    TMP_DIR="${CONFIG_processing_tmp_dir:-/hpc/temp/srivatsan_s/SPC_genome_preprocessing}"
    REFERENCE_GENOME="${CONFIG_reference_genome_dir:-/shared/biodata/reference/GATK/hg38}"
    READ1="${CONFIG_data_read1}"
    READ2="${CONFIG_data_read2}"
    READ_COUNT="${CONFIG_data_read_count}"
    OUTPUT_NAME="${CONFIG_output_sample_name}"
else
    # Hardcoded defaults if config.yaml doesn't exist
    SCRIPTS_DIR="."  # Project root directory
    OUTPUT_DIR="./data/K562_tree"
    N_CHUNKS=500
    TMP_DIR="/hpc/temp/srivatsan_s/SPC_genome_preprocessing"
    REFERENCE_GENOME="/shared/biodata/reference/GATK/hg38"
fi

mkdir -p "${OUTPUT_DIR}"

# Function to display help message
show_help() {
    cat << EOF
Usage: $0 -o <OUTPUT_NAME> -1 <read1.fastq.gz> -2 <read2.fastq.gz> -g <reference_genome.fa> -r <read_count>

This script uses defaults from config.yaml if present. Command-line arguments override config.yaml values.

Required arguments:
  -o    <output_name>           Desired sample name (prefix for outputs)
  -1    <read1.fastq.gz>        Read 1 FASTQ file
  -2    <read2.fastq.gz>        Read 2 FASTQ file
  -g    <reference_genome>      Path to directory containing genome fasta, fasta index, and BWA index folder
  -r    <READ_COUNT>            Number of reads (from sequencing run info)

Optional arguments:
  -s    <scripts_DIR>           Path to the SPC_genome directory (default: from config.yaml or ./SPC_genome)
  -O    <output_dir>            Desired output directory (default: from config.yaml or .)
  -n    <N_CHUNKS>              Number of subjobs for SLURM arrays (default: from config.yaml or 500)
  -t    <TMP_DIR>               Temp directory for fastq chunks (default: from config.yaml or /hpc/temp/srivatsan_s/SPC_genome_preprocessing)
  -h                            Show this help message and exit

Note: Values are applied in this order: hardcoded defaults < config.yaml < command-line arguments
EOF
}

# Parse command-line options (these override config.yaml values)
while getopts ":o:1:2:g:r:s:O:n:t:h" option; do
  case $option in
    o) OUTPUT_NAME=$OPTARG ;;
    1) READ1=$OPTARG ;;
    2) READ2=$OPTARG ;;
    g) REFERENCE_GENOME=$OPTARG ;;
	r) READ_COUNT=$OPTARG ;;
    s) SCRIPTS_DIR=$OPTARG ;;
	O) OUTPUT_DIR=$OPTARG ;;
	n) N_CHUNKS=$OPTARG ;;
	t) TMP_DIR=$OPTARG ;;
    h) show_help; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; show_help; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; show_help; exit 1 ;;
  esac
done

# Check that required arguments are provided
if [ -z "$OUTPUT_NAME" ] || [ -z "$READ1" ] || [ -z "$READ2" ] || [ -z "$READ_COUNT" ]; then
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
# Limit to specified read count, then split into chunks
zcat "$READ1" | head -n $total_lines | split -l $CHUNK_LINES - "$TMP_DIR/read1_chunk_" &
zcat "$READ2" | head -n $total_lines | split -l $CHUNK_LINES - "$TMP_DIR/read2_chunk_" &
wait

ls "$TMP_DIR"/read1_chunk_* | sed 's/.*chunk_//' > "${OUTPUT_DIR}/chunk_indices.txt"

chunk_count=$(wc -l < "${OUTPUT_DIR}/chunk_indices.txt")

######################################################################################################
#### Submit First Job array

PP_array_ID=$(sbatch --parsable --array=1-$chunk_count "${SCRIPTS_DIR}/scripts/PP_array.sh" "${OUTPUT_DIR}/chunk_indices.txt" "${REFERENCE_GENOME}" "${SCRIPTS_DIR}" "${TMP_DIR}")

echo "Preprocessing array job ID: ${PP_array_ID}"

######################################################################################################
#### Concatenate SAM files, create BAM, and detect real cells

concat_job_ID=$(sbatch --parsable --dependency=afterok:$PP_array_ID "${SCRIPTS_DIR}/scripts/concatenate.sh" "${OUTPUT_NAME}" "${TMP_DIR}" "${OUTPUT_DIR}" "${SCRIPTS_DIR}")

echo "Concatenation and cell detection job ID: ${concat_job_ID}"

######################################################################################################
#### Submit single cell extraction arrays

sc_from_chunks_array_ID=$(sbatch --parsable --dependency=afterok:$concat_job_ID "${SCRIPTS_DIR}/scripts/sc_from_chunks.sh" "${OUTPUT_DIR}/chunk_indices.txt" "${TMP_DIR}" "${OUTPUT_DIR}/real_cells.txt" "${SCRIPTS_DIR}")

echo "Compiling single cell bam files with array job ID: ${sc_array_ID}"

######################################################################################################
#### Submit single cell variant calling array
cell_count=$(wc -l < "${OUTPUT_DIR}/real_cells.txt")
mkdir -p "${OUTPUT_DIR}/sc_ouputs"

sc_var_array_ID=$(sbatch --parsable --array=1-$cell_count --dependency=afterok:$sc_comp_array_ID "${SCRIPTS_DIR}/sc_var_array.sh" "${TMP_DIR}" "${OUTPUT_DIR}/real_cells.txt" "${REFERENCE_GENOME}" "${OUTPUT_DIR}/sc_ouputs")

echo "Compiling single cell bam files with array job ID: ${sc_array_ID}"

######################################################################################################
#### Submit matched bulk variant calling job

# bulk_v_job_ID=$(sbatch --parsable --dependency=afterok:$comp_job_ID )

# echo "Bulk calling job array ID: ${bulk_v_job_ID}"

######################################################################################################
#### joint calling

#### Prep for joint calling:

for f in sc_outputs/*.g.vcf.gz; do
  [[ -f "${f}.tbi" ]] || continue  # skip if no index file
  bc=$(basename "$f" .g.vcf.gz)
  echo -e "${bc}\t${f}"
done > "${OUTPUT_DIR}/barcodes.map"



# joint_calling_array_id=$(sbatch --parsable --dependency=afterok:$sc_var_array_ID "${SCRIPTS_DIR}/joint_calling_array.sh")

# echo "Joint Calling Job ID: ${joint_calling_job_id}"

######################################################################################################
#### Submit anndata job