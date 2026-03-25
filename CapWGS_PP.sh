#!/bin/bash
#SBATCH --job-name=Preprocessing
#SBATCH --output=SLURM_outs/%x_%j.out
#SBATCH -c 2

mkdir -p SLURM_outs/

# Function to wait for a SLURM job to complete
wait_for_job() {
    local job_id=$1
    while squeue -j "$job_id" -h -t PENDING,RUNNING 2>/dev/null | grep -q .; do
        sleep 60
    done
    # Check if any tasks failed
    if sacct -j "$job_id" --format=State -n | grep -q FAILED; then
        echo "Job $job_id had failures" >&2
        return 1
    fi
}

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
    SCRIPTS_DIR="."  # Project root directory
    N_CHUNKS="${CONFIG_processing_n_chunks:-500}"
    TMP_DIR_BASE="${CONFIG_processing_tmp_dir:-/hpc/temp/srivatsan_s/SPC_genome_preprocessing}"
    REFERENCE_GENOME="${CONFIG_reference_genome_dir:-/shared/biodata/reference/GATK/hg38}"
    VARIANT_CALLER="${CONFIG_variant_caller:-bcftools}"
    READ1="${CONFIG_data_read1}"
    READ2="${CONFIG_data_read2}"
    READ_COUNT="${CONFIG_data_read_count}"
    OUTPUT_NAME="${CONFIG_output_sample_name}"
else
    # Hardcoded defaults if config.yaml doesn't exist
    SCRIPTS_DIR="."  # Project root directory
    N_CHUNKS=500
    TMP_DIR_BASE="/hpc/temp/srivatsan_s/SPC_genome_preprocessing"
    REFERENCE_GENOME="/shared/biodata/reference/GATK/hg38"
    VARIANT_CALLER="bcftools"
fi

# Function to display help message
show_help() {
    cat << EOF
Usage: $0 -o <OUTPUT_NAME> -1 <read1.fastq.gz> -2 <read2.fastq.gz> -g <reference_genome.fa> -r <read_count>

This script processes CapWGS data: alignment, cell detection, and variant calling.
Uses defaults from config.yaml if present. Command-line arguments override config.yaml values.

Required arguments:
  -o    <output_name>           Sample name (creates data/{sample}/ and results/{sample}/)
  -1    <read1.fastq.gz>        Read 1 FASTQ file(s) - can be single file or quoted pattern (e.g., "lane*_R1*.fastq.gz")
  -2    <read2.fastq.gz>        Read 2 FASTQ file(s) - can be single file or quoted pattern (e.g., "lane*_R2*.fastq.gz")
  -g    <reference_genome>      Path to directory containing genome fasta, fasta index, and BWA index folder
  -r    <READ_COUNT>            Number of reads (from sequencing run info)

Optional arguments:
  -s    <scripts_DIR>           Path to the SPC_genome directory (default: from config.yaml or .)
  -n    <N_CHUNKS>              Number of subjobs for SLURM arrays (default: from config.yaml or 500)
  -t    <TMP_DIR>               Temp directory for fastq chunks (default: from config.yaml or /hpc/temp/srivatsan_s/SPC_genome_preprocessing/{sample}/)
  -c    <cell_count>            Override automatic cell detection - select top N cells by read count (optional)
  -v    <variant_caller>        Variant caller: bcftools, gatk, or none (default: from config.yaml or bcftools)
                                  bcftools: Faster, suitable for shallow pilot runs
                                  gatk: More rigorous, follows GATK best practices (includes mark duplicates and BQSR)
                                  none: Skip variant calling, run QC only
  --use-existing-chunks         Skip fastq splitting, use existing chunks in temp directory
                                  Requires chunk_indices.txt to exist in results directory
                                  Useful for development/testing and retry scenarios
  -h                            Show this help message and exit

Directory structure created:
  - data/{sample}/                 Bulk alignments
  - data/{sample}/sc_outputs/      Single-cell BAMs and VCFs
  - results/{sample}/              Final results and QC metrics
  - {TMP_DIR}/{sample}/            Temporary chunks (deleted after completion)

Note: Values are applied in this order: hardcoded defaults < config.yaml < command-line arguments
EOF
}

# Initialize flag for using existing chunks
USE_EXISTING_CHUNKS=false

# Parse long options manually before getopts
for arg in "$@"; do
    if [ "$arg" = "--use-existing-chunks" ]; then
        USE_EXISTING_CHUNKS=true
    fi
done

# Parse command-line options (these override config.yaml values)
while getopts ":o:1:2:g:r:s:n:t:c:v:h" option; do
  case $option in
    o) OUTPUT_NAME=$OPTARG ;;
    1) READ1=$OPTARG ;;
    2) READ2=$OPTARG ;;
    g) REFERENCE_GENOME=$OPTARG ;;
	r) READ_COUNT=$OPTARG ;;
    s) SCRIPTS_DIR=$OPTARG ;;
	n) N_CHUNKS=$OPTARG ;;
	t) TMP_DIR_BASE=$OPTARG ;;
	c) CELL_COUNT=$OPTARG ;;
	v) VARIANT_CALLER=$OPTARG ;;
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

# Set default for CELL_COUNT if not provided
CELL_COUNT="${CELL_COUNT:-}"

# Validate variant caller option (convert to lowercase for comparison)
VARIANT_CALLER_LOWER=$(echo "$VARIANT_CALLER" | tr '[:upper:]' '[:lower:]')
if [ "$VARIANT_CALLER_LOWER" != "bcftools" ] && [ "$VARIANT_CALLER_LOWER" != "gatk" ] && [ "$VARIANT_CALLER_LOWER" != "none" ]; then
    echo "Error: Invalid variant caller '${VARIANT_CALLER}'. Must be 'bcftools', 'gatk', or 'none'." >&2
    show_help
    exit 1
fi

# Normalize variant caller to lowercase
VARIANT_CALLER="$VARIANT_CALLER_LOWER"

# Extract sample name from output path (handles paths like "benchmarking_coverage/HSC2_enzyme")
SAMPLE_NAME=$(basename "${OUTPUT_NAME}")

# Set up sample-specific directory structure
DATA_DIR="${SCRIPTS_DIR}/data/${OUTPUT_NAME}"
SC_OUTPUTS_DIR="${DATA_DIR}/sc_outputs"
RESULTS_DIR="${SCRIPTS_DIR}/results/${OUTPUT_NAME}"
TMP_DIR="${TMP_DIR_BASE}/${SAMPLE_NAME}"

# Create necessary directories
mkdir -p "${DATA_DIR}"
mkdir -p "${SC_OUTPUTS_DIR}"
mkdir -p "${RESULTS_DIR}"
mkdir -p "${TMP_DIR}"

# Print inputs
echo "=========================================="
echo "CapWGS Pipeline Configuration"
echo "=========================================="
echo "Sample Name: ${OUTPUT_NAME}"
echo "Read 1: ${READ1}"
echo "Read 2: ${READ2}"
echo "Read Count: ${READ_COUNT}"
echo "Reference Genome: ${REFERENCE_GENOME}"
echo "Variant Caller: ${VARIANT_CALLER}"
echo ""
echo "Scripts Directory: ${SCRIPTS_DIR}"
echo "Data Directory: ${DATA_DIR}"
echo "Single Cell Outputs: ${SC_OUTPUTS_DIR}"
echo "Results Directory: ${RESULTS_DIR}"
echo "Temp Directory: ${TMP_DIR}"
echo "Number of Chunks: ${N_CHUNKS}"
echo "=========================================="

######################################################################################################
#### chunk the fastqs (or use existing chunks)

if [ "$USE_EXISTING_CHUNKS" = true ]; then
    echo "=========================================="
    echo "Using existing fastq chunks"
    echo "=========================================="
    echo "Temp directory: ${TMP_DIR}"
    echo "Results directory: ${RESULTS_DIR}"
    echo ""

    # Verify chunk_indices.txt exists
    if [ ! -f "${RESULTS_DIR}/chunk_indices.txt" ]; then
        echo "ERROR: chunk_indices.txt not found at ${RESULTS_DIR}/chunk_indices.txt"
        echo "This file is required when using --use-existing-chunks"
        echo ""
        echo "chunk_indices.txt should contain the list of chunk suffixes (one per line)"
        echo "Example: 00, 01, 02, etc."
        exit 1
    fi

    chunk_count=$(wc -l < "${RESULTS_DIR}/chunk_indices.txt")
    echo "Found ${chunk_count} chunks in chunk_indices.txt"
    echo "Skipping fastq splitting step"
    echo "=========================================="
    echo ""
else
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

    ls "$TMP_DIR"/read1_chunk_* | sed 's/.*chunk_//' > "${RESULTS_DIR}/chunk_indices.txt"

    chunk_count=$(wc -l < "${RESULTS_DIR}/chunk_indices.txt")
fi

######################################################################################################
#### Submit preprocessing array and wait for completion

PP_array_ID=$(sbatch --parsable --array=1-$chunk_count "${SCRIPTS_DIR}/scripts/CapWGS/PP_array.sh" "${RESULTS_DIR}/chunk_indices.txt" "${REFERENCE_GENOME}" "${SCRIPTS_DIR}" "${TMP_DIR}" "${SAMPLE_NAME}")

echo "Preprocessing array job ID: ${PP_array_ID}"

# Wait for preprocessing array to complete
echo "Waiting for preprocessing array to complete..."
wait_for_job "${PP_array_ID}"
echo "Preprocessing array complete!"

######################################################################################################
#### Concatenate SAM files, create BAM, and detect real cells (runs inline)

echo "=========================================="
echo "Concatenating SAM files and detecting cells..."
echo "=========================================="

# Activate conda environment for python scripts
eval "$(micromamba shell hook --shell bash)"
micromamba activate spc_genome

# Give filesystem a moment to sync (NFS delay)
sleep 5

# Debug: show what we're looking for
echo "Looking for SAM files in: ${TMP_DIR}"
sam_count=$(ls "${TMP_DIR}"/*.sam 2>/dev/null | wc -l)
echo "Found ${sam_count} SAM files"

ls "${TMP_DIR}"/*.sam > "${RESULTS_DIR}/sam_list.txt"

BAM_FILE="${DATA_DIR}/${SAMPLE_NAME}.bam"

module load SAMtools

echo "Merging SAM chunks directly to sorted BAM..."
# Merge and sort in one streaming operation
# -c: Combine @RG headers with colliding IDs (don't alter them to be distinct)
# -f 0x2: Filter for properly paired reads to remove duplicate records from chunking
samtools merge -@ 4 -c -b "${RESULTS_DIR}/sam_list.txt" -O SAM - | \
    samtools view -@ 4 -f 0x2 -b - | \
    samtools sort -@ 4 -o "${BAM_FILE}"

echo "Indexing BAM file..."
samtools index -@ 4 "${BAM_FILE}"

echo "Computing read counts per barcode..."
samtools view -@ 4 "${BAM_FILE}" | python "${SCRIPTS_DIR}/scripts/utils/readcounts.py" -o "${RESULTS_DIR}/readcounts.csv"

# Calculate barcode assignment statistics
echo ""
echo "Calculating barcode assignment statistics..."
total_reads_in_bam=$(tail -n +2 "${RESULTS_DIR}/readcounts.csv" | cut -d',' -f2 | awk '{s+=$1} END {print s}')

ASSIGNMENT_STATS="${RESULTS_DIR}/barcode_assignment_stats.txt"
{
    echo "=========================================="
    echo "Barcode Assignment Statistics"
    echo "=========================================="
    echo ""

    if [ -n "${READ_COUNT}" ] && [ "${READ_COUNT}" -gt 0 ]; then
        # Multiply READ_COUNT by 2 for paired-end reads
        total_input_reads=$((READ_COUNT * 2))
        assignment_rate=$(awk "BEGIN {printf \"%.2f\", ($total_reads_in_bam / $total_input_reads) * 100}")
        echo "Total input reads (paired-end): ${total_input_reads}"
        echo "Reads assigned to valid barcodes: ${total_reads_in_bam}"
        echo "Assignment rate: ${assignment_rate}%"
    else
        echo "Reads assigned to valid barcodes: ${total_reads_in_bam}"
        echo "(Input read count not provided, cannot calculate assignment rate)"
    fi

    echo ""
    echo "=========================================="
} > "${ASSIGNMENT_STATS}"

cat "${ASSIGNMENT_STATS}"
echo "Barcode assignment statistics written to: ${ASSIGNMENT_STATS}"
echo ""

# Cell detection
if [ -n "${CELL_COUNT}" ]; then
    echo "Using user-provided cell count: ${CELL_COUNT}"
    echo "Selecting top ${CELL_COUNT} cells by read count..."
    # Generate knee plot for QC
    cat "${RESULTS_DIR}/readcounts.csv" | python "${SCRIPTS_DIR}/scripts/utils/detect_cells.py" --plot "${RESULTS_DIR}/kneeplot.png" > /dev/null || true
    # Select top N cells
    tail -n +2 "${RESULTS_DIR}/readcounts.csv" | \
        sort -t',' -k2 -nr | \
        head -${CELL_COUNT} | \
        cut -d',' -f1 > "${RESULTS_DIR}/real_cells.txt"
else
    echo "Using automatic cell detection (knee plot)..."
    cat "${RESULTS_DIR}/readcounts.csv" | \
        python "${SCRIPTS_DIR}/scripts/utils/detect_cells.py" \
        --plot "${RESULTS_DIR}/kneeplot.png" > "${RESULTS_DIR}/real_cells.txt"
fi

DETECTED_CELL_COUNT=$(wc -l < "${RESULTS_DIR}/real_cells.txt")
echo "Selected ${DETECTED_CELL_COUNT} cells"
echo "Concatenation and cell detection complete!"
echo "=========================================="
echo ""

######################################################################################################
#### Submit unified single cell processing array (extraction + preprocessing + QC)

# Single array job that does EVERYTHING for each cell:
# 1. Extract reads from bulk BAM
# 2. MarkDuplicates + BQSR
# 3. Generate bigwig
# 4. Generate Lorenz curve
# 5. Collect QC metrics

BULK_BAM="${DATA_DIR}/${SAMPLE_NAME}.bam"
QC_METRICS_DIR="${RESULTS_DIR}/qc_metrics"
mkdir -p "${QC_METRICS_DIR}"

echo "Submitting unified single-cell processing array for ${DETECTED_CELL_COUNT} cells..."

sc_unified_job_ID=$(sbatch --parsable \
    --array=1-${DETECTED_CELL_COUNT} \
    "${SCRIPTS_DIR}/scripts/CapWGS/sc_extract_preprocess_qc_array.sh" \
    "${BULK_BAM}" \
    "${RESULTS_DIR}/real_cells.txt" \
    "${SC_OUTPUTS_DIR}" \
    "${QC_METRICS_DIR}" \
    "${REFERENCE_GENOME}" \
    "${SCRIPTS_DIR}" \
    1000)

echo "Single-cell unified processing job ID: ${sc_unified_job_ID}"
echo "Array size: 1-${DETECTED_CELL_COUNT}"

# Wait for unified processing array to complete
echo "Waiting for single-cell processing array to complete..."
wait_for_job "${sc_unified_job_ID}"
echo "Single-cell processing array complete!"

# All outputs (BAMs, bigwigs, Lorenz curves, QC metrics) created by this single array
FINAL_SC_BAM_DIR="${SC_OUTPUTS_DIR}"
preprocessing_dependency=$sc_unified_job_ID

######################################################################################################
#### Submit single cell variant calling array (skip if variant_caller=none)

if [ "$VARIANT_CALLER" != "none" ]; then
    # This wrapper script will:
    # 1. Read real_cells.txt to determine how many cells were detected
    # 2. Submit the variant calling array job with the correct array size
    # Both GATK and bcftools modes use preprocessed BAMs (after MarkDuplicates + BQSR)

    submit_sc_var_job_ID=$(sbatch --parsable "${SCRIPTS_DIR}/scripts/CapWGS/submit_sc_variant_calling.sh" "${FINAL_SC_BAM_DIR}" "${RESULTS_DIR}/real_cells.txt" "${REFERENCE_GENOME}" "${FINAL_SC_BAM_DIR}" "${SCRIPTS_DIR}" "${VARIANT_CALLER}")

    echo "Single cell variant calling submission job ID: ${submit_sc_var_job_ID}"

    # Wait for variant calling to complete
    echo "Waiting for single-cell variant calling to complete..."
    wait_for_job "${submit_sc_var_job_ID}"
    echo "Single-cell variant calling complete!"

    ######################################################################################################
    #### Joint calling

    # This wrapper script will:
    # 1. Generate genomic intervals for parallelization (GATK mode)
    # 2. Create barcodes.map from VCF/GVCF files
    # 3. Submit joint calling array job
    #    - GATK mode: GenomicsDBImport + GenotypeGVCFs per interval
    #    - BCFtools mode: bcftools merge + normalize

    submit_jc_job_ID=$(sbatch --parsable "${SCRIPTS_DIR}/scripts/CapWGS/submit_joint_calling.sh" "${REFERENCE_GENOME}" "${RESULTS_DIR}" "${SC_OUTPUTS_DIR}" "${RESULTS_DIR}" "${SCRIPTS_DIR}" "${SAMPLE_NAME}" "${VARIANT_CALLER}" "${TMP_DIR}")

    echo "Joint calling submission job ID: ${submit_jc_job_ID}"
    echo "(Joint calling is the final output - no wait needed)"
else
    echo "Variant calling skipped (VARIANT_CALLER=none)"
    echo "Running QC-only mode"
fi

######################################################################################################
#### Compile QC results (runs inline)

echo "=========================================="
echo "Compiling QC Results"
echo "=========================================="
echo "Single-cell BAM directory: ${FINAL_SC_BAM_DIR}"
echo "QC metrics directory: ${QC_METRICS_DIR}"
echo "Results directory: ${RESULTS_DIR}"
echo "=========================================="
echo ""

# Activate python environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

# Compile Lorenz curves
echo "Compiling Lorenz curves..."
LORENZ_OUTPUT="${RESULTS_DIR}/compiled_lorenz_curves.csv"

if ls "${SC_OUTPUTS_DIR}"/*_lorenz.csv 1> /dev/null 2>&1; then
    python scripts/CapWGS_QC/compile_lorenz.py \
        "${SC_OUTPUTS_DIR}" \
        "${LORENZ_OUTPUT}"
    echo "✓ Lorenz curves compiled: ${LORENZ_OUTPUT}"
else
    echo "⚠ No Lorenz curve files found in ${SC_OUTPUTS_DIR}"
fi

echo ""

# Compile benchmarking QC metrics
echo "Compiling benchmarking QC metrics..."
QC_OUTPUT="${RESULTS_DIR}/compiled_qc_metrics.csv"

if ls "${QC_METRICS_DIR}"/*_alignment_metrics.txt 1> /dev/null 2>&1; then
    python scripts/CapWGS_QC/compile_qc_metrics.py \
        "${QC_METRICS_DIR}" \
        "${QC_OUTPUT}"
    echo "✓ QC metrics compiled: ${QC_OUTPUT}"
else
    echo "⚠ No QC metrics files found in ${QC_METRICS_DIR}"
fi

echo ""

# Generate Lander-Waterman plot
echo "Generating Lander-Waterman coverage plot..."
LW_PLOT="${RESULTS_DIR}/lander_waterman_coverage.png"

if [ -f "${QC_OUTPUT}" ]; then
    # Get genome length from reference
    GENOME_LENGTH=$(samtools faidx "${REFERENCE_GENOME}" 2>/dev/null && awk '{sum+=$2} END {print sum}' "${REFERENCE_GENOME}.fai")

    if [ -n "${GENOME_LENGTH}" ] && [ "${GENOME_LENGTH}" -gt 0 ]; then
        python scripts/utils/plot_lander_waterman.py \
            "${QC_OUTPUT}" \
            "${GENOME_LENGTH}" \
            "${LW_PLOT}"
        echo "✓ Lander-Waterman plot generated: ${LW_PLOT}"
    else
        echo "⚠ Could not determine genome length - skipping Lander-Waterman plot"
    fi
else
    echo "⚠ No compiled QC metrics - skipping Lander-Waterman plot"
fi

echo ""
echo "=========================================="
echo "QC Compilation Complete!"
echo "=========================================="
echo "Output files:"
if [ -f "${LORENZ_OUTPUT}" ]; then
    echo "  ✓ ${LORENZ_OUTPUT}"
fi
if [ -f "${QC_OUTPUT}" ]; then
    echo "  ✓ ${QC_OUTPUT}"
fi
if [ -f "${LW_PLOT}" ]; then
    echo "  ✓ ${LW_PLOT}"
fi
echo "=========================================="