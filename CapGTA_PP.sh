#!/bin/bash
#SBATCH --job-name=GTA_Preprocessing
#SBATCH --output=SLURM_outs/%x_%j.out
#SBATCH -c 2

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
    REFERENCE_GENOME="${CONFIG_reference_genome_dir:-/shared/biodata/reference/iGenomes/Caenorhabditis_elegans/UCSC/ce10/Sequence}"
    READ1="${CONFIG_data_read1}"
    READ2="${CONFIG_data_read2}"
    READ_COUNT="${CONFIG_data_read_count}"
    OUTPUT_NAME="${CONFIG_output_sample_name}"
else
    # Hardcoded defaults if config.yaml doesn't exist
    SCRIPTS_DIR="."  # Project root directory
    OUTPUT_DIR="./data/PolE_worm_pilot"
    N_CHUNKS=500
    TMP_DIR="/hpc/temp/srivatsan_s/SPC_genome_preprocessing"
    REFERENCE_GENOME="/shared/biodata/reference/iGenomes/Caenorhabditis_elegans/UCSC/ce10/Sequence"
fi

mkdir -p "${OUTPUT_DIR}"

# Function to display help message
show_help() {
    cat << EOF
Usage: $0 -o <OUTPUT_NAME> -1 <read1.fastq.gz> -2 <read2.fastq.gz> -g <reference_genome> -r <read_count>

This script processes genome-transcriptome coassay (CapGTA) data, separating DNA and RNA reads.
Uses defaults from config.yaml if present. Command-line arguments override config.yaml values.

Required arguments:
  -o    <output_name>           Desired sample name (prefix for outputs)
  -1    <read1.fastq.gz>        Read 1 FASTQ file(s) - can be single file or quoted pattern (e.g., "lane*_R1*.fastq.gz")
  -2    <read2.fastq.gz>        Read 2 FASTQ file(s) - can be single file or quoted pattern (e.g., "lane*_R2*.fastq.gz")
  -g    <reference_genome>      Path to directory containing genome sequence and indices (BWAIndex, STARIndex)
  -r    <READ_COUNT>            Number of reads (from sequencing run info)

Optional arguments:
  -s    <scripts_DIR>           Path to the SPC_genome directory (default: from config.yaml or .)
  -O    <output_dir>            Desired output directory (default: from config.yaml or ./data/PolE_worm_pilot)
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

# Create sample-specific TMP_DIR to avoid collisions when running multiple samples in parallel
# Save the base TMP_DIR first
TMP_DIR_BASE="${TMP_DIR}"
TMP_DIR="${TMP_DIR_BASE}/${OUTPUT_NAME}"

# Print inputs
echo "Processing $READ_COUNT reads for genome-transcriptome coassay with input parameters:"
echo "Sample Name: ${OUTPUT_NAME}"
echo "Read 1: ${READ1}"
echo "Read 2: ${READ2}"
echo "Reference Genome: ${REFERENCE_GENOME}"

echo "Scripts Directory: ${SCRIPTS_DIR}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "Temp Directory: ${TMP_DIR}"

# Create necessary directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TMP_DIR}"

# Set up directory structure
BIN_DIR="${SCRIPTS_DIR}/bin/${OUTPUT_NAME}"
ALIGNED_DIR="${SCRIPTS_DIR}/data/${OUTPUT_NAME}/aligned"
RESULTS_DIR="${SCRIPTS_DIR}/results/${OUTPUT_NAME}"

mkdir -p "${BIN_DIR}"
mkdir -p "${ALIGNED_DIR}"
mkdir -p "${RESULTS_DIR}"
echo "Bin Directory: ${BIN_DIR}"
echo "Aligned Data Directory: ${ALIGNED_DIR}"
echo "Results Directory: ${RESULTS_DIR}"

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
# Use -d for numeric suffixes (00, 01, 02...) to avoid file extension collisions (e.g., 'gz')
# Note: READ1 and READ2 are intentionally unquoted to allow shell expansion for multi-lane patterns
zcat $READ1 | head -n $total_lines | split -d -l $CHUNK_LINES - "$TMP_DIR/read1_chunk_" &
zcat $READ2 | head -n $total_lines | split -d -l $CHUNK_LINES - "$TMP_DIR/read2_chunk_" &
wait

ls "$TMP_DIR"/read1_chunk_* | sed 's/.*chunk_//' > "${BIN_DIR}/chunk_indices.txt"

chunk_count=$(wc -l < "${BIN_DIR}/chunk_indices.txt")

######################################################################################################
#### Submit STAR-only alignment array job
#### Aligns all reads with STAR, then splits by splice junctions

PP_array_ID=$(sbatch --parsable --array=1-$chunk_count "${SCRIPTS_DIR}/scripts/PP_array_gta_star_only.sh" "${BIN_DIR}/chunk_indices.txt" "${REFERENCE_GENOME}" "${SCRIPTS_DIR}" "${TMP_DIR}")

echo "STAR-only alignment preprocessing array job ID: ${PP_array_ID}"

######################################################################################################
#### Concatenate DNA and RNA SAM files, create BAMs, and detect real cells

concat_job_ID=$(sbatch --parsable --dependency=afterok:$PP_array_ID "${SCRIPTS_DIR}/scripts/concatenate_gta.sh" "${OUTPUT_NAME}" "${TMP_DIR}" "${ALIGNED_DIR}" "${BIN_DIR}" "${SCRIPTS_DIR}")

echo "Concatenation and cell detection job ID: ${concat_job_ID}"

######################################################################################################
#### Extract single cells (DNA and RNA BAMs separately)

# This job extracts reads for each detected cell from the chunked SAM files
sc_from_chunks_job_ID=$(sbatch --parsable --dependency=afterok:$concat_job_ID "${SCRIPTS_DIR}/scripts/sc_from_chunks_gta.sh" "${BIN_DIR}/chunk_indices.txt" "${TMP_DIR}" "${BIN_DIR}/real_cells.txt" "${SCRIPTS_DIR}")

echo "Single cell extraction job ID: ${sc_from_chunks_job_ID}"

######################################################################################################
#### Pipeline completion

echo ""
echo "=========================================="
echo "CapGTA pipeline submitted successfully!"
echo "=========================================="
echo "Pipeline will:"
echo "  1. Demultiplex and trim reads"
echo "  2. Align DNA reads with BWA-MEM"
echo "  3. Align RNA reads with STAR"
echo "  4. Detect real cells using combined DNA+RNA counts"
echo "  5. Extract per-cell DNA and RNA BAMs"
echo ""
echo "Output locations:"
echo "  - Bulk DNA BAM: ${ALIGNED_DIR}/${OUTPUT_NAME}_dna.bam"
echo "  - Bulk RNA BAM: ${ALIGNED_DIR}/${OUTPUT_NAME}_rna.bam"
echo "  - Knee plot: ${BIN_DIR}/kneeplot.png"
echo "  - Real cells list: ${BIN_DIR}/real_cells.txt"
echo "  - Read counts: ${BIN_DIR}/readcounts.csv"
echo "  - SC outputs: ${ALIGNED_DIR}/sc_outputs/"
echo ""
