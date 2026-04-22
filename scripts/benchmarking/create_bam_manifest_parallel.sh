#!/bin/bash
##########################################################################################################################
# Parallelized BAM Manifest Builder
#
# Purpose: Create BAM manifest by parallelizing the slow read counting step
# Step 1: Creates basic metadata file
# Step 2: Submits SLURM array to count reads in parallel
# Step 3: Compiles results into final manifest
##########################################################################################################################

set -euo pipefail

# Check arguments
if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    echo "Usage: $0 <input_bam_list.txt> <output_manifest.tsv> [--wait]"
    echo ""
    echo "Options:"
    echo "  --wait    Wait for SLURM jobs to complete before returning"
    echo ""
    echo "Example:"
    echo "  $0 data/benchmarking/manifests/input_bams.txt \\"
    echo "     data/benchmarking/manifests/bam_manifest.tsv --wait"
    exit 1
fi

INPUT_LIST="$1"
OUTPUT_TSV="$2"
WAIT_FLAG="${3:-}"

# Validate input file exists
if [ ! -f "$INPUT_LIST" ]; then
    echo "Error: Input file not found: $INPUT_LIST"
    exit 1
fi

# Temporary directory for intermediate files
TEMP_DIR="$(dirname "$OUTPUT_TSV")/tmp_manifest_$(date +%s)"
mkdir -p "$TEMP_DIR"

echo "==========================================="
echo "Parallelized BAM Manifest Builder"
echo "==========================================="
echo "Input list: $INPUT_LIST"
echo "Output TSV: $OUTPUT_TSV"
echo "Temp dir: $TEMP_DIR"
echo ""

# Count BAMs
TOTAL_BAMS=$(grep -v "^#" "$INPUT_LIST" | grep -v "^[[:space:]]*$" | wc -l)
echo "Total BAMs to process: $TOTAL_BAMS"
echo ""

##########################################################################################################################
# Step 1: Create basic metadata file (fast - no read counting)
##########################################################################################################################

echo "Step 1: Creating basic metadata file..."
echo ""

METADATA_FILE="${TEMP_DIR}/metadata.tsv"

# Write header
echo -e "task_id\tbam_path\tdataset\tcell_id" > "$METADATA_FILE"

TASK_ID=0
while IFS= read -r BAM_PATH || [ -n "$BAM_PATH" ]; do
    # Skip empty lines and comments
    [[ -z "$BAM_PATH" || "$BAM_PATH" =~ ^[[:space:]]*# ]] && continue

    TASK_ID=$((TASK_ID + 1))

    # Validate BAM exists
    if [ ! -f "$BAM_PATH" ]; then
        echo "  WARNING: BAM file not found, skipping: $BAM_PATH"
        continue
    fi

    # Extract dataset (parent directory name)
    PARENT_DIR=$(dirname "$BAM_PATH")
    DATASET=$(basename "$PARENT_DIR")

    # If parent is 'sc_outputs', go up one more level
    if [ "$DATASET" = "sc_outputs" ]; then
        DATASET=$(basename "$(dirname "$PARENT_DIR")")
    fi

    # Extract cell_id (basename without .bam)
    CELL_ID=$(basename "$BAM_PATH" .bam)

    # Write to metadata file
    echo -e "${TASK_ID}\t${BAM_PATH}\t${DATASET}\t${CELL_ID}" >> "$METADATA_FILE"

done < "$INPUT_LIST"

VALID_TASKS=$(tail -n +2 "$METADATA_FILE" | wc -l)
echo "✓ Metadata file created: $VALID_TASKS valid BAMs"
echo ""

##########################################################################################################################
# Step 2: Submit SLURM array job to count reads in parallel
##########################################################################################################################

echo "Step 2: Submitting SLURM array job to count reads..."
echo ""

# Create the array job script
ARRAY_SCRIPT="${TEMP_DIR}/count_reads_array.sh"

cat > "$ARRAY_SCRIPT" <<'ARRAY_EOF'
#!/bin/bash
#SBATCH -J count_reads
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 30:00

set -euo pipefail

METADATA_FILE="$1"
TEMP_DIR="$2"

# Get task info
TASK_LINE=$(tail -n +2 "$METADATA_FILE" | sed -n "${SLURM_ARRAY_TASK_ID}p")
IFS=$'\t' read -r TASK_ID BAM_PATH DATASET CELL_ID <<< "$TASK_LINE"

echo "Task ${TASK_ID}: Processing ${CELL_ID}"

# Load samtools
module load SAMtools

# Extract reference from BAM header
HEADER=$(samtools view -H "$BAM_PATH" 2>/dev/null || echo "")

# Try to find reference in @PG lines
REFERENCE=$(echo "$HEADER" | grep "^@PG" | grep -oP '(?<=-x\s|--reference\s|-R\s)\S+' | head -1 || echo "")

# If that didn't work, make educated guess from @SQ lines
if [ -z "$REFERENCE" ]; then
    FIRST_SEQ=$(echo "$HEADER" | grep "^@SQ" | head -1 | grep -oP 'SN:\K\S+' || echo "")

    if [[ "$FIRST_SEQ" =~ ^chr ]]; then
        if [ -f "/shared/biodata/reference/GATK/hg38/Homo_sapiens_assembly38.fasta" ]; then
            REFERENCE="/shared/biodata/reference/GATK/hg38/Homo_sapiens_assembly38.fasta"
        elif [ -f "/shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa" ]; then
            REFERENCE="/shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa"
        else
            REFERENCE="UNKNOWN_hg38"
        fi
    else
        REFERENCE="UNKNOWN_${FIRST_SEQ}"
    fi
fi

# Count primary mapped reads
MAPPED_READS=$(samtools view -c -F 260 "$BAM_PATH" 2>/dev/null || echo "0")

# Write result
OUTPUT_FILE="${TEMP_DIR}/task_${TASK_ID}.tsv"
echo -e "${DATASET}\t${CELL_ID}\t${BAM_PATH}\t${REFERENCE}\t${MAPPED_READS}" > "$OUTPUT_FILE"

echo "✓ Task ${TASK_ID} complete: ${MAPPED_READS} reads"
ARRAY_EOF

chmod +x "$ARRAY_SCRIPT"

# Submit array job
cd /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome

JOB=$(sbatch --array=1-${VALID_TASKS} "$ARRAY_SCRIPT" "$METADATA_FILE" "$TEMP_DIR")
JOB_ID=$(echo "$JOB" | awk '{print $NF}')

echo "$JOB"
echo ""

##########################################################################################################################
# Step 3: Create compilation script (runs after array completes)
##########################################################################################################################

echo "Step 3: Creating compilation job..."
echo ""

COMPILE_SCRIPT="${TEMP_DIR}/compile_manifest.sh"

cat > "$COMPILE_SCRIPT" <<COMPILE_EOF
#!/bin/bash
#SBATCH -J compile_manifest
#SBATCH -o SLURM_outs/compile_manifest_%j.out
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 10:00

set -euo pipefail

echo "Compiling manifest from parallel results..."
echo ""

# Create output with header
echo -e "dataset\tcell_id\tbam_path\treference\ttotal_mapped_reads" > "$OUTPUT_TSV"

# Concatenate all task results
cat ${TEMP_DIR}/task_*.tsv >> "$OUTPUT_TSV"

# Count results
TOTAL_ROWS=\$(tail -n +2 "$OUTPUT_TSV" | wc -l)

echo "✓ Manifest compiled: \$TOTAL_ROWS cells"
echo "Output: $OUTPUT_TSV"
echo ""

# Show preview
echo "Preview:"
head -5 "$OUTPUT_TSV" | column -t -s \$'\t'
echo ""

# Clean up temp directory
echo "Cleaning up temporary files..."
rm -rf "$TEMP_DIR"
echo "✓ Done!"
COMPILE_EOF

chmod +x "$COMPILE_SCRIPT"

# Submit compilation job (depends on array completion)
COMPILE_JOB=$(sbatch --dependency=afterok:${JOB_ID} "$COMPILE_SCRIPT")
COMPILE_JOB_ID=$(echo "$COMPILE_JOB" | awk '{print $NF}')

echo "$COMPILE_JOB"
echo ""

##########################################################################################################################
# Summary
##########################################################################################################################

echo "==========================================="
echo "✓ Jobs submitted successfully!"
echo "==========================================="
echo ""
echo "Job chain:"
echo "  1. Count reads:  $JOB_ID ($VALID_TASKS array tasks)"
echo "  2. Compile:      $COMPILE_JOB_ID (waits for job 1)"
echo ""
echo "Monitor progress:"
echo "  squeue -j $JOB_ID"
echo "  squeue -j $COMPILE_JOB_ID"
echo ""
echo "Final output will be: $OUTPUT_TSV"
echo ""

# Wait for completion if requested
if [ "$WAIT_FLAG" = "--wait" ]; then
    echo "Waiting for jobs to complete..."
    echo ""

    # Wait for compilation job to finish
    while squeue -j $COMPILE_JOB_ID 2>/dev/null | grep -q $COMPILE_JOB_ID; do
        sleep 10
    done

    echo "✓ Jobs complete!"
    echo ""

    if [ -f "$OUTPUT_TSV" ]; then
        echo "Manifest created successfully:"
        head -5 "$OUTPUT_TSV" | column -t -s $'\t'
        echo ""
    else
        echo "Error: Manifest file not found. Check SLURM logs."
        exit 1
    fi
fi

echo "==========================================="
