#!/bin/bash
#SBATCH --job-name=GTA_Preprocessing
#SBATCH --output=SLURM_outs/%x_%j.out
#SBATCH -c 2
#SBATCH -t 12-00:00:00

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
    REFERENCE_GENOME="${CONFIG_reference_genome_dir:-data/reference/worm_GCA_028201515.1_combined}"
    GTF_FILE="${CONFIG_reference_gtf}"
    READ1="${CONFIG_data_read1}"
    READ2="${CONFIG_data_read2}"
    READ_COUNT="${CONFIG_data_read_count}"
    OUTPUT_NAME="${CONFIG_output_sample_name}"
else
    # Hardcoded defaults if config.yaml doesn't exist
    SCRIPTS_DIR="."  # Project root directory
    N_CHUNKS=500
    TMP_DIR_BASE="/hpc/temp/srivatsan_s/SPC_genome_preprocessing"
    REFERENCE_GENOME="data/reference/worm_GCA_028201515.1_combined"
    GTF_FILE=""
fi

# Function to display help message
show_help() {
    cat << EOF
Usage: $0 -o <OUTPUT_NAME> -1 <read1.fastq.gz> -2 <read2.fastq.gz> -g <reference_genome> -a <gtf_file> -r <read_count>

This script processes genome-transcriptome coassay (CapGTA) data, separating DNA and RNA reads.
Includes alignment, cell detection, single-cell extraction, QC metrics, variant calling, and RNA count matrix generation.
Uses defaults from config.yaml if present. Command-line arguments override config.yaml values.

Required arguments:
  -o    <output_name>           Sample name (creates data/{sample}/ and results/{sample}/)
  -1    <read1.fastq.gz>        Read 1 FASTQ file(s) - can be single file or quoted pattern (e.g., "lane*_R1*.fastq.gz")
  -2    <read2.fastq.gz>        Read 2 FASTQ file(s) - can be single file or quoted pattern (e.g., "lane*_R2*.fastq.gz")
  -g    <reference_genome>      Path to directory containing genome sequence and indices (BWAIndex, STARIndex)
  -a    <gtf_file>              GTF annotation file for exonic enrichment calculation
  -r    <READ_COUNT>            Number of reads (from sequencing run info)

Optional arguments:
  -s    <scripts_DIR>           Path to the SPC_genome directory (default: from config.yaml or .)
  -n    <N_CHUNKS>              Number of subjobs for SLURM arrays (default: from config.yaml or 500)
  -t    <TMP_DIR>               Temp directory for fastq chunks (default: from config.yaml or /hpc/temp/srivatsan_s/SPC_genome_preprocessing/{sample}/)
  -c    <cell_count>            Override automatic cell detection - select top N cells by read count (optional)
  --use-existing-chunks         Skip fastq splitting, use existing chunks in temp directory
                                Requires chunk_indices.txt to exist in results directory
                                Useful for development/testing and retry scenarios
  -h                            Show this help message and exit

Directory structure created:
  - data/{sample}/                 Bulk DNA and RNA alignments
  - data/{sample}/sc_outputs/      Single-cell DNA/RNA BAMs and VCFs
  - data/{sample}/qc_metrics/      Per-cell QC metrics
  - results/{sample}/              Joint calling VCF, QC metrics, RNA count matrices
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
while getopts ":o:1:2:g:a:r:s:n:t:c:h" option; do
  case $option in
    o) OUTPUT_NAME=$OPTARG ;;
    1) READ1=$OPTARG ;;
    2) READ2=$OPTARG ;;
    g) REFERENCE_GENOME=$OPTARG ;;
    a) GTF_FILE=$OPTARG ;;
	r) READ_COUNT=$OPTARG ;;
    s) SCRIPTS_DIR=$OPTARG ;;
	n) N_CHUNKS=$OPTARG ;;
	t) TMP_DIR_BASE=$OPTARG ;;
	c) CELL_COUNT=$OPTARG ;;
    h) show_help; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; show_help; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; show_help; exit 1 ;;
  esac
done

# Check that required arguments are provided
if [ -z "$OUTPUT_NAME" ] || [ -z "$READ1" ] || [ -z "$READ2" ] || [ -z "$READ_COUNT" ] || [ -z "$GTF_FILE" ]; then
    echo "Missing required arguments." >&2
    show_help
    exit 1
fi

# Set default for CELL_COUNT if not provided
CELL_COUNT="${CELL_COUNT:-}"

# Extract sample name from output path (handles paths like "worm_CapGTA_L4_pilot/L4_3_UDI7")
SAMPLE_NAME=$(basename "${OUTPUT_NAME}")

# Set up sample-specific directory structure
DATA_DIR="${SCRIPTS_DIR}/data/${OUTPUT_NAME}"
SC_OUTPUTS_DIR="${DATA_DIR}/sc_outputs"
QC_METRICS_DIR="${DATA_DIR}/qc_metrics"
RESULTS_DIR="${SCRIPTS_DIR}/results/${OUTPUT_NAME}"
TMP_DIR="${TMP_DIR_BASE}/${SAMPLE_NAME}"

# Create necessary directories
mkdir -p "${DATA_DIR}"
mkdir -p "${SC_OUTPUTS_DIR}"
mkdir -p "${QC_METRICS_DIR}"
mkdir -p "${RESULTS_DIR}"
mkdir -p "${TMP_DIR}"

# Print inputs
echo "=========================================="
echo "CapGTA Pipeline Configuration"
echo "=========================================="
echo "Sample Name: ${OUTPUT_NAME}"
echo "Read 1: ${READ1}"
echo "Read 2: ${READ2}"
echo "Read Count: ${READ_COUNT}"
echo "Reference Genome: ${REFERENCE_GENOME}"
echo "GTF File: ${GTF_FILE}"
echo ""
echo "Scripts Directory: ${SCRIPTS_DIR}"
echo "Data Directory: ${DATA_DIR}"
echo "Single Cell Outputs: ${SC_OUTPUTS_DIR}"
echo "QC Metrics Directory: ${QC_METRICS_DIR}"
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
#### Submit STAR-only alignment array job and wait for completion

PP_array_ID=$(sbatch --parsable --array=1-$chunk_count "${SCRIPTS_DIR}/scripts/CapGTA/PP_array_gta_star_only.sh" "${RESULTS_DIR}/chunk_indices.txt" "${REFERENCE_GENOME}" "${SCRIPTS_DIR}" "${TMP_DIR}")

echo "STAR-only alignment preprocessing array job ID: ${PP_array_ID}"

# Wait for preprocessing array to complete
echo "Waiting for alignment array to complete..."
wait_for_job "${PP_array_ID}"
echo "Alignment array complete!"

######################################################################################################
#### Two-stage parallel SAM merge (DNA and RNA separately, runs inline)

echo "=========================================="
echo "Concatenating SAM files and creating BAMs..."
echo "=========================================="

# Activate conda environment for python scripts
eval "$(micromamba shell hook --shell bash)"
micromamba activate spc_genome

# Give filesystem a moment to sync (NFS delay)
sleep 5

module load SAMtools

echo "=========================================="
echo "Two-Stage Parallel SAM Merge (DNA)"
echo "=========================================="

# Stage 1: Create groups of DNA SAM files for parallel merging
MERGE_GROUP_SIZE=20  # Number of SAMs per group (configurable)
GROUP_DIR_DNA="${TMP_DIR}/merge_groups_dna"
INTERMEDIATE_DIR_DNA="${TMP_DIR}/intermediate_bams_dna"

echo "Stage 1: Preparing parallel merge groups (DNA)..."
mkdir -p "${GROUP_DIR_DNA}"
mkdir -p "${INTERMEDIATE_DIR_DNA}"

ls "${TMP_DIR}"/*_dna.sam > "${RESULTS_DIR}/sam_list_dna.txt"

# Split sam_list.txt into groups
split -l ${MERGE_GROUP_SIZE} -d -a 2 "${RESULTS_DIR}/sam_list_dna.txt" "${GROUP_DIR_DNA}/group_"

# Count groups created
N_GROUPS_DNA=$(ls ${GROUP_DIR_DNA}/group_* 2>/dev/null | wc -l)
echo "  Created ${N_GROUPS_DNA} merge groups of ~${MERGE_GROUP_SIZE} SAMs each"
echo ""

# Stage 2: Submit parallel merge array for DNA
echo "Stage 2: Submitting parallel merge array for DNA (${N_GROUPS_DNA} jobs)..."
parallel_merge_job_dna=$(sbatch --parsable --array=1-${N_GROUPS_DNA} \
    "${SCRIPTS_DIR}/scripts/utils/parallel_merge_array.sh" \
    "${GROUP_DIR_DNA}" \
    "${INTERMEDIATE_DIR_DNA}" \
    "${SAMPLE_NAME}_dna")

echo "  Parallel merge job ID (DNA): ${parallel_merge_job_dna}"

echo "=========================================="
echo "Two-Stage Parallel SAM Merge (RNA)"
echo "=========================================="

# Stage 1: Create groups of RNA SAM files for parallel merging
GROUP_DIR_RNA="${TMP_DIR}/merge_groups_rna"
INTERMEDIATE_DIR_RNA="${TMP_DIR}/intermediate_bams_rna"

echo "Stage 1: Preparing parallel merge groups (RNA)..."
mkdir -p "${GROUP_DIR_RNA}"
mkdir -p "${INTERMEDIATE_DIR_RNA}"

ls "${TMP_DIR}"/*_rna.sam > "${RESULTS_DIR}/sam_list_rna.txt"

# Split sam_list.txt into groups
split -l ${MERGE_GROUP_SIZE} -d -a 2 "${RESULTS_DIR}/sam_list_rna.txt" "${GROUP_DIR_RNA}/group_"

# Count groups created
N_GROUPS_RNA=$(ls ${GROUP_DIR_RNA}/group_* 2>/dev/null | wc -l)
echo "  Created ${N_GROUPS_RNA} merge groups of ~${MERGE_GROUP_SIZE} SAMs each"
echo ""

# Stage 2: Submit parallel merge array for RNA
echo "Stage 2: Submitting parallel merge array for RNA (${N_GROUPS_RNA} jobs)..."
parallel_merge_job_rna=$(sbatch --parsable --array=1-${N_GROUPS_RNA} \
    "${SCRIPTS_DIR}/scripts/utils/parallel_merge_array.sh" \
    "${GROUP_DIR_RNA}" \
    "${INTERMEDIATE_DIR_RNA}" \
    "${SAMPLE_NAME}_rna")

echo "  Parallel merge job ID (RNA): ${parallel_merge_job_rna}"
echo ""

# Wait for both DNA and RNA parallel merge arrays to complete
echo "Waiting for DNA and RNA parallel merge arrays to complete..."
wait_for_job "${parallel_merge_job_dna}"
echo "  ✓ DNA parallel merge complete!"
wait_for_job "${parallel_merge_job_rna}"
echo "  ✓ RNA parallel merge complete!"
echo ""

# Stage 3: Final merge and sort for DNA
echo "Stage 3: Final merge and sort (DNA)..."
ls "${INTERMEDIATE_DIR_DNA}"/*.bam > "${RESULTS_DIR}/intermediate_bam_list_dna.txt"

N_INTERMEDIATE_DNA=$(wc -l < "${RESULTS_DIR}/intermediate_bam_list_dna.txt")
echo "  Merging ${N_INTERMEDIATE_DNA} intermediate BAMs into final sorted DNA BAM..."

# Final merge and sort (unfiltered first for accurate flagstat)
# -c: Combine @RG headers with colliding IDs (preserve read groups)
DNA_BAM="${DATA_DIR}/${SAMPLE_NAME}_dna.bam"
UNFILTERED_DNA_BAM="${DNA_BAM%.bam}.unfiltered.bam"

samtools merge -@ 4 -c -b "${RESULTS_DIR}/intermediate_bam_list_dna.txt" -O BAM - | \
    samtools sort -@ 4 -o "${UNFILTERED_DNA_BAM}"

echo "  ✓ Unfiltered DNA BAM created: ${UNFILTERED_DNA_BAM}"

# Run flagstat on unfiltered DNA BAM to get true mapping rate
echo "Running samtools flagstat (DNA, unfiltered)..."
FLAGSTAT_DNA="${RESULTS_DIR}/flagstat_dna.txt"
samtools flagstat -@ 4 "${UNFILTERED_DNA_BAM}" > "${FLAGSTAT_DNA}"
echo "DNA flagstat written to: ${FLAGSTAT_DNA}"
echo ""

# Filter for properly paired reads and create final DNA BAM
echo "Filtering for properly paired reads (DNA)..."
samtools view -@ 4 -f 0x2 -b "${UNFILTERED_DNA_BAM}" | \
    samtools sort -@ 4 -o "${DNA_BAM}"
rm "${UNFILTERED_DNA_BAM}"

echo "  ✓ Final DNA BAM created: ${DNA_BAM}"
echo ""

# Stage 3: Final merge and sort for RNA
echo "Stage 3: Final merge and sort (RNA)..."
ls "${INTERMEDIATE_DIR_RNA}"/*.bam > "${RESULTS_DIR}/intermediate_bam_list_rna.txt"

N_INTERMEDIATE_RNA=$(wc -l < "${RESULTS_DIR}/intermediate_bam_list_rna.txt")
echo "  Merging ${N_INTERMEDIATE_RNA} intermediate BAMs into final sorted RNA BAM..."

# Final merge and sort (unfiltered first for accurate flagstat)
RNA_BAM="${DATA_DIR}/${SAMPLE_NAME}_rna.bam"
UNFILTERED_RNA_BAM="${RNA_BAM%.bam}.unfiltered.bam"

samtools merge -@ 4 -c -b "${RESULTS_DIR}/intermediate_bam_list_rna.txt" -O BAM - | \
    samtools sort -@ 4 -o "${UNFILTERED_RNA_BAM}"

echo "  ✓ Unfiltered RNA BAM created: ${UNFILTERED_RNA_BAM}"

# Run flagstat on unfiltered RNA BAM to get true mapping rate
echo "Running samtools flagstat (RNA, unfiltered)..."
FLAGSTAT_RNA="${RESULTS_DIR}/flagstat_rna.txt"
samtools flagstat -@ 4 "${UNFILTERED_RNA_BAM}" > "${FLAGSTAT_RNA}"
echo "RNA flagstat written to: ${FLAGSTAT_RNA}"
echo ""

# Filter for properly paired reads and create final RNA BAM
echo "Filtering for properly paired reads (RNA)..."
samtools view -@ 4 -f 0x2 -b "${UNFILTERED_RNA_BAM}" | \
    samtools sort -@ 4 -o "${RNA_BAM}"
rm "${UNFILTERED_RNA_BAM}"

echo "  ✓ Final RNA BAM created: ${RNA_BAM}"
echo ""

echo "=========================================="
echo "Two-Stage Merge Complete!"
echo "=========================================="
echo ""

# Index both BAMs
echo "Indexing BAM files..."
samtools index -@ 4 "${DNA_BAM}"
samtools index -@ 4 "${RNA_BAM}"
echo ""

######################################################################################################
#### Compute read counts and detect cells (runs inline)

echo "Computing read counts per barcode (DNA and RNA)..."
samtools view -@ 4 "${DNA_BAM}" | python "${SCRIPTS_DIR}/scripts/utils/readcounts.py" -o "${RESULTS_DIR}/readcounts_dna.csv"
samtools view -@ 4 "${RNA_BAM}" | python "${SCRIPTS_DIR}/scripts/utils/readcounts.py" -o "${RESULTS_DIR}/readcounts_rna.csv"

# Combine DNA and RNA read counts
echo "Combining DNA and RNA read counts..."
python "${SCRIPTS_DIR}/scripts/utils/combine_readcounts_gta.py" \
    "${RESULTS_DIR}/readcounts_dna.csv" \
    "${RESULTS_DIR}/readcounts_rna.csv" \
    "${RESULTS_DIR}/readcounts.csv"

# Calculate barcode assignment statistics
echo ""
echo "Calculating barcode assignment statistics..."
total_dna_reads=$(tail -n +2 "${RESULTS_DIR}/readcounts_dna.csv" | cut -d',' -f2 | awk '{s+=$1} END {print s}')
total_rna_reads=$(tail -n +2 "${RESULTS_DIR}/readcounts_rna.csv" | cut -d',' -f2 | awk '{s+=$1} END {print s}')
total_reads_in_bams=$((total_dna_reads + total_rna_reads))

ASSIGNMENT_STATS="${RESULTS_DIR}/barcode_assignment_stats.txt"
{
    echo "=========================================="
    echo "Barcode Assignment Statistics"
    echo "=========================================="
    echo ""

    if [ -n "${READ_COUNT}" ] && [ "${READ_COUNT}" -gt 0 ]; then
        # Multiply READ_COUNT by 2 for paired-end reads
        total_input_reads=$((READ_COUNT * 2))
        assignment_rate=$(awk "BEGIN {printf \"%.2f\", ($total_reads_in_bams / $total_input_reads) * 100}")
        echo "Total input reads (paired-end): ${total_input_reads}"
        echo "Reads assigned to valid barcodes (DNA + RNA): ${total_reads_in_bams}"
        echo "  - DNA reads: ${total_dna_reads}"
        echo "  - RNA reads: ${total_rna_reads}"
        echo "Assignment rate: ${assignment_rate}%"
    else
        echo "Reads assigned to valid barcodes (DNA + RNA): ${total_reads_in_bams}"
        echo "  - DNA reads: ${total_dna_reads}"
        echo "  - RNA reads: ${total_rna_reads}"
        echo "(Input read count not provided, cannot calculate assignment rate)"
    fi

    echo ""
    echo "=========================================="
} > "${ASSIGNMENT_STATS}"

cat "${ASSIGNMENT_STATS}"
echo "Barcode assignment statistics written to: ${ASSIGNMENT_STATS}"
echo ""

# Cell detection (using combined DNA+RNA read counts)
if [ -n "${CELL_COUNT}" ]; then
    echo "Using user-provided cell count: ${CELL_COUNT}"
    echo "Selecting top ${CELL_COUNT} cells by total read count..."
    # Generate knee plot for QC
    cat "${RESULTS_DIR}/readcounts.csv" | python "${SCRIPTS_DIR}/scripts/utils/detect_cells.py" --plot "${RESULTS_DIR}/kneeplot.png" > /dev/null || true
    # Select top N cells (use read_count column, which equals total_reads)
    tail -n +2 "${RESULTS_DIR}/readcounts.csv" | \
        sort -t',' -k5 -nr | \
        head -${CELL_COUNT} | \
        cut -d',' -f1 > "${RESULTS_DIR}/real_cells.txt"
else
    echo "Using automatic cell detection (knee plot)..."
    # Use read_count column (same as total_reads) for cell detection
    cat "${RESULTS_DIR}/readcounts.csv" | \
        python "${SCRIPTS_DIR}/scripts/utils/detect_cells.py" \
        --plot "${RESULTS_DIR}/kneeplot.png" > "${RESULTS_DIR}/real_cells.txt"
fi

DETECTED_CELL_COUNT=$(wc -l < "${RESULTS_DIR}/real_cells.txt")
echo "Selected ${DETECTED_CELL_COUNT} cells"
echo "Cell detection complete!"
echo "=========================================="
echo ""

######################################################################################################
#### Submit unified single cell processing array (extraction + preprocessing + QC)

echo "Submitting unified single-cell processing array for ${DETECTED_CELL_COUNT} cells..."

sc_unified_job_ID=$(sbatch --parsable \
    --array=1-${DETECTED_CELL_COUNT} \
    "${SCRIPTS_DIR}/scripts/CapGTA/sc_extract_preprocess_qc_array.sh" \
    "${DNA_BAM}" \
    "${RNA_BAM}" \
    "${RESULTS_DIR}/real_cells.txt" \
    "${SC_OUTPUTS_DIR}" \
    "${QC_METRICS_DIR}" \
    "${REFERENCE_GENOME}" \
    "${GTF_FILE}" \
    "${SCRIPTS_DIR}" \
    1000)

echo "Single-cell unified processing job ID: ${sc_unified_job_ID}"
echo "Array size: 1-${DETECTED_CELL_COUNT}"

# Wait for unified processing array to complete
echo "Waiting for single-cell processing array to complete..."
wait_for_job "${sc_unified_job_ID}"
echo "Single-cell processing array complete!"

######################################################################################################
#### Compile QC results (runs inline)

echo "=========================================="
echo "Compiling QC Results"
echo "=========================================="
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

if ls "${QC_METRICS_DIR}"/*_lorenz.csv 1> /dev/null 2>&1; then
    if python scripts/utils/compile_lorenz.py \
        "${QC_METRICS_DIR}" \
        "${LORENZ_OUTPUT}"; then
        if [ -f "${LORENZ_OUTPUT}" ]; then
            echo "✓ Lorenz curves compiled: ${LORENZ_OUTPUT}"
        else
            echo "⚠ Lorenz curve compilation completed but output file not found"
        fi
    else
        echo "⚠ Lorenz curve compilation failed (check log above for errors)"
    fi
else
    echo "⚠ No Lorenz curve files found in ${QC_METRICS_DIR}"
fi

echo ""

# Compile QC metrics
echo "Compiling QC metrics..."
QC_OUTPUT="${RESULTS_DIR}/compiled_qc_metrics.csv"

if ls "${QC_METRICS_DIR}"/*_alignment_metrics.txt 1> /dev/null 2>&1; then
    if python scripts/utils/compile_qc_metrics.py \
        "${QC_METRICS_DIR}" \
        "${QC_OUTPUT}"; then
        if [ -f "${QC_OUTPUT}" ]; then
            echo "✓ QC metrics compiled: ${QC_OUTPUT}"
        else
            echo "⚠ QC metrics compilation completed but output file not found"
        fi
    else
        echo "⚠ QC metrics compilation failed (check log above for errors)"
    fi
else
    echo "⚠ No QC metrics files found in ${QC_METRICS_DIR}"
fi

echo ""

# Compile exonic enrichment metrics
echo "Compiling exonic enrichment metrics..."
EXONIC_ENRICHMENT_OUTPUT="${RESULTS_DIR}/compiled_exonic_enrichment.csv"

if ls "${QC_METRICS_DIR}"/*_exonic_enrichment_dna.txt 1> /dev/null 2>&1; then
    if python scripts/utils/compile_exonic_enrichment.py \
        "${QC_METRICS_DIR}" \
        "${EXONIC_ENRICHMENT_OUTPUT}"; then
        if [ -f "${EXONIC_ENRICHMENT_OUTPUT}" ]; then
            echo "✓ Exonic enrichment metrics compiled: ${EXONIC_ENRICHMENT_OUTPUT}"
        else
            echo "⚠ Exonic enrichment compilation completed but output file not found"
        fi
    else
        echo "⚠ Exonic enrichment compilation failed (check log above for errors)"
    fi
else
    echo "⚠ No exonic enrichment files found in ${QC_METRICS_DIR}"
fi

echo ""

# Generate Lander-Waterman plot for DNA coverage
echo "Generating Lander-Waterman coverage plot (DNA)..."
LW_PLOT="${RESULTS_DIR}/lander_waterman_coverage_dna.png"

if [ -f "${QC_OUTPUT}" ]; then
    if python scripts/utils/plot_lander_waterman.py \
        "${QC_OUTPUT}" \
        "${REFERENCE_GENOME}" \
        "${LW_PLOT}"; then
        echo "✓ Lander-Waterman plot generated: ${LW_PLOT}"
    else
        echo "⚠ Lander-Waterman plot generation failed (check log above for errors)"
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
if [ -f "${EXONIC_ENRICHMENT_OUTPUT}" ]; then
    echo "  ✓ ${EXONIC_ENRICHMENT_OUTPUT}"
fi
if [ -f "${LW_PLOT}" ]; then
    echo "  ✓ ${LW_PLOT}"
fi
echo "=========================================="
echo ""

######################################################################################################
#### Submit single cell variant calling and RNA count matrix generation

# This wrapper script will submit variant calling and RNA count matrix generation
submit_sc_final_job_ID=$(sbatch --parsable "${SCRIPTS_DIR}/scripts/CapGTA/submit_final_processing.sh" \
    "${SC_OUTPUTS_DIR}" \
    "${RESULTS_DIR}/real_cells.txt" \
    "${REFERENCE_GENOME}" \
    "${GTF_FILE}" \
    "${RESULTS_DIR}" \
    "${SCRIPTS_DIR}" \
    "${OUTPUT_NAME}")

echo "Single cell variant calling and RNA count matrix job ID: ${submit_sc_final_job_ID}"
echo "(Final processing - no wait needed)"

echo ""
echo "=========================================="
echo "CapGTA Pipeline Submitted Successfully!"
echo "=========================================="
echo "Pipeline will:"
echo "  1. ✓ Demultiplex and trim reads"
echo "  2. ✓ Align all reads with STAR (separates DNA and RNA by splice junctions)"
echo "  3. ✓ Parallel merge SAM files (DNA and RNA)"
echo "  4. ✓ Detect real cells using combined DNA+RNA counts"
echo "  5. ✓ Extract per-cell DNA and RNA BAMs"
echo "  6. ✓ Generate comprehensive QC metrics (Picard, Lorenz, exonic enrichment)"
echo "  7. → Call variants on DNA BAMs with BCFtools"
echo "  8. → Generate RNA count matrix from RNA BAMs"
echo ""
echo "Output locations:"
echo "  - Bulk DNA BAM: ${DNA_BAM}"
echo "  - Bulk RNA BAM: ${RNA_BAM}"
echo "  - Knee plot: ${RESULTS_DIR}/kneeplot.png"
echo "  - Real cells list: ${RESULTS_DIR}/real_cells.txt"
echo "  - Read counts: ${RESULTS_DIR}/readcounts.csv"
echo "  - QC metrics: ${RESULTS_DIR}/compiled_qc_metrics.csv"
echo "  - Exonic enrichment: ${RESULTS_DIR}/compiled_exonic_enrichment.csv"
echo "  - SC BAMs and VCFs: ${SC_OUTPUTS_DIR}/"
echo "  - Merged VCF: ${RESULTS_DIR}/sc_variants_merged.vcf.gz"
echo "  - RNA count matrix: ${RESULTS_DIR}/rna_counts_matrix.csv"
echo ""
