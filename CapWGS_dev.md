# CapWGS Pipeline Development Notes

## Implementation Plan

### Overview
This document outlines the planned improvements to the CapWGS and bulk GATK pipelines to improve QC coverage, performance, and efficiency.

### Key Improvements
1. **QC Integration**: Add comprehensive QC to `CapWGS_PP.sh` for both GATK and bcftools modes
2. **Architecture Change**: Move MarkDuplicates/BQSR to single-cell level (GATK mode only)
3. **Parallelization**: Replace chromosome-level joint calling with ~100MB interval-based parallelization

---

## TODO #1: Consolidate QC into CapWGS_PP.sh

### Goal
Add comprehensive QC checks to `CapWGS_PP.sh` so it always runs QC alongside variant calling for both GATK and bcftools modes.

### Components

#### 1.1: Add Barcode Assignment Metric
**File**: `scripts/CapWGS/concatenate.sh`

**Changes**:
- After `readcounts.py` runs, sum reads from `readcounts.csv` (much faster than re-counting BAM)
- Calculate assignment rate: `(reads_in_bam / READ_COUNT) * 100`
- Print statistics to log with clear formatting
- Pass `READ_COUNT` parameter from `CapWGS_PP.sh` to `concatenate.sh`

**Implementation**:
```bash
# After line 38 (readcounts.py call)
total_reads_in_bam=$(tail -n +2 "${BIN_DIR}/readcounts.csv" | cut -d',' -f2 | awk '{s+=$1} END {print s}')
assignment_rate=$(awk "BEGIN {printf \"%.2f\", ($total_reads_in_bam / ${READ_COUNT}) * 100}")

echo "=========================================="
echo "Barcode Assignment Statistics"
echo "=========================================="
echo "Total reads (input): ${READ_COUNT}"
echo "Reads assigned to valid barcodes: ${total_reads_in_bam}"
echo "Assignment rate: ${assignment_rate}%"
echo "=========================================="
```

**Files to modify**:
- `scripts/CapWGS/concatenate.sh`: Add READ_COUNT parameter and statistics output
- `CapWGS_PP.sh` line 192: Pass READ_COUNT to concatenate.sh

#### 1.2: Add QC Job Submissions to CapWGS_PP.sh
**File**: `CapWGS_PP.sh`

**Changes**:
Add QC array job submissions after single-cell extraction (after line 233):

1. **Bigwig and Lorenz curves**: Submit `submit_bigwig_lorenz.sh`
2. **Benchmarking QC metrics**: Submit `submit_benchmarking_qc.sh`
3. **QC compilation**: Submit compilation job after both arrays complete

**Implementation**:
```bash
######################################################################################################
#### Submit QC analysis arrays

# Bigwig and Lorenz curve generation (for coverage uniformity assessment)
bigwig_lorenz_job_ID=$(sbatch --parsable \
    --dependency=afterok:$preprocessing_dependency \
    "${SCRIPTS_DIR}/scripts/CapWGS_QC/submit_bigwig_lorenz.sh" \
    "${RESULTS_DIR}/real_cells.txt" \
    "${SC_OUTPUTS_DIR}" \
    "${SC_OUTPUTS_DIR}" \
    1000)

echo "Bigwig and Lorenz curve generation job ID: ${bigwig_lorenz_job_ID}"

# Benchmarking QC metrics (alignment, GC bias, duplicates, coverage)
benchmarking_qc_job_ID=$(sbatch --parsable \
    --dependency=afterok:$preprocessing_dependency \
    "${SCRIPTS_DIR}/scripts/CapWGS_QC/submit_benchmarking_qc.sh" \
    "${RESULTS_DIR}/real_cells.txt" \
    "${SC_OUTPUTS_DIR}/qc_metrics" \
    "${REFERENCE_GENOME}/Homo_sapiens_assembly38.fasta" \
    "${SC_OUTPUTS_DIR}")

echo "Benchmarking QC metrics collection job ID: ${benchmarking_qc_job_ID}"
```

**Note**: QC runs in parallel with variant calling (independent operations)

#### 1.3: Add QC Compilation Step
**File**: `CapWGS_PP.sh`

**Changes**:
Add compilation job after joint calling submission (after line 259):

```bash
######################################################################################################
#### Compile QC results

compile_qc_job_ID=$(sbatch --parsable \
    --dependency=afterok:${bigwig_lorenz_job_ID}:${benchmarking_qc_job_ID} \
    "${SCRIPTS_DIR}/scripts/CapWGS_QC/compile_qc_results.sh" \
    "${SC_OUTPUTS_DIR}" \
    "${SC_OUTPUTS_DIR}/qc_metrics" \
    "${RESULTS_DIR}")

echo "QC compilation job ID: ${compile_qc_job_ID}"
```

#### 1.4: Create QC Compilation Script
**New file**: `scripts/CapWGS_QC/compile_qc_results.sh`

**Purpose**: Aggregate Lorenz curves and QC metrics into summary CSV files

**Implementation**:
```bash
#!/bin/bash
#SBATCH -J compile_qc
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 30

set -euo pipefail

SC_OUTPUTS_DIR="$1"
QC_METRICS_DIR="$2"
RESULTS_DIR="$3"

eval "$(micromamba shell hook --shell bash)"
micromamba activate default_jupyter

echo "Compiling QC results..."

# Compile Lorenz curves
python scripts/CapWGS_QC/compile_lorenz.py \
    "${SC_OUTPUTS_DIR}" \
    "${RESULTS_DIR}/compiled_lorenz_curves.csv"

# Compile benchmarking QC metrics
python scripts/CapWGS_QC/compile_qc_metrics.py \
    "${QC_METRICS_DIR}" \
    "${RESULTS_DIR}/compiled_qc_metrics.csv"

echo "QC compilation complete!"
echo "Outputs:"
echo "  - ${RESULTS_DIR}/compiled_lorenz_curves.csv"
echo "  - ${RESULTS_DIR}/compiled_qc_metrics.csv"
```

#### 1.5: Update submit_benchmarking_qc.sh
**File**: `scripts/CapWGS_QC/submit_benchmarking_qc.sh`

**Changes**:
Update to accept 4 positional parameters matching new calling convention:
- $1: SAMPLE_LIST (real_cells.txt)
- $2: OUTPUT_DIR (qc_metrics directory)
- $3: REFERENCE (full path to fasta)
- $4: BAM_DIR (single-cell BAM directory)

---

## TODO #2: Move MarkDuplicates/BQSR to Single-Cell Level (GATK mode only)

### Rationale
Currently, MarkDuplicates and BQSR run on the bulk concatenated BAM (GATK mode only). This is extremely slow for large datasets. Moving these operations to single-cell level provides:

1. **Massive parallelization**: N parallel jobs instead of 1 slow bulk job
2. **Better duplicate detection**: Duplicates within single cells are more biologically meaningful
3. **Cell-specific BQSR**: Base quality recalibration tailored to each cell's coverage
4. **Faster pipeline**: Array jobs complete much faster than bulk processing
5. **Better resource utilization**: Many small jobs vs one huge memory-intensive job

### Architecture Changes

#### Current GATK Pipeline Flow:
```
1. PP_array (alignment)
2. concatenate (merge chunks, detect cells)
3. sc_from_bam (extract single cells)
4. markdup_bqsr (bulk BAM preprocessing) ← SLOW
5. sc_from_bam again (extract from preprocessed BAM)
6. sc_variant_calling (HaplotypeCaller)
7. joint_calling (GenomicsDB + GenotypeGVCFs)
```

#### New GATK Pipeline Flow (UNIFIED ARCHITECTURE):
```
Within CapWGS_PP.sh:
1. Parse args (no change)
2. Split fastqs (no change)
3. Submit PP_array → WAIT
4. Concatenate SAMs, detect cells (moved from concatenate.sh, runs in CapWGS_PP.sh)
5. Submit unified SC array (extraction + preprocessing + QC) → WAIT
6. Submit SC variant calling array → WAIT
7. Submit joint calling (submit_joint_calling.sh keeps internal arrays, no wait needed)
8. Compile QC metrics (moved from compile_qc_results.sh, runs in CapWGS_PP.sh)

Key Changes:
- CapWGS_PP.sh submits arrays directly (no wrapper scripts)
- Uses `srun --dependency=afterok:$array_ID --wait=0 true` to wait for arrays
- Single-threaded steps run directly in CapWGS_PP.sh
- Eliminates race conditions from wrapper completion
```

#### BCFtools Pipeline Flow (unchanged):
```
1. PP_array (alignment)
2. concatenate (merge chunks, detect cells)
3. sc_from_bam (extract single cells)
4. QC arrays (bigwig, lorenz, benchmarking)
5. sc_variant_calling (bcftools mpileup/call)
6. joint_calling (bcftools merge + normalize)
```

### Implementation: Refactor Main Pipeline to Eliminate Race Conditions

#### 2.1: Core Principle
**CapWGS_PP.sh serves as the main scheduler:**
- Runs single-threaded housekeeping tasks directly (fastq splitting, concatenation, QC compilation)
- Submits array jobs directly (no wrapper scripts)
- Waits for arrays to complete using `srun --dependency=afterok:$array_ID --wait=0 true`
- Eliminates all wrapper script race conditions

#### 2.2: Changes to CapWGS_PP.sh

**Move concatenation inline** (from `concatenate.sh`):
- After PP_array wait completes, run SAM merging directly
- Run cell detection directly
- Read cell count: `CELL_COUNT=$(wc -l < "${RESULTS_DIR}/real_cells.txt")`

**Submit unified SC array directly**:
```bash
sc_unified_job_ID=$(sbatch --parsable \
    --array=1-${CELL_COUNT} \
    "${SCRIPTS_DIR}/scripts/CapWGS/sc_extract_preprocess_qc_array.sh" \
    "${BULK_BAM}" \
    "${RESULTS_DIR}/real_cells.txt" \
    "${SC_OUTPUTS_DIR}" \
    "${QC_METRICS_DIR}" \
    "${REFERENCE_GENOME}" \
    "${SCRIPTS_DIR}" \
    1000)

# Wait for array to complete
srun --dependency=afterok:${sc_unified_job_ID} --wait=0 true
```

**Move QC compilation inline** (from `compile_qc_results.sh`):
- After SC unified array wait completes, run compilation directly

#### 2.3: Create Single-Cell Preprocessing Array Script
**New file**: `scripts/CapWGS/sc_preprocessing_array.sh`

**Purpose**: Run MarkDuplicates and BQSR on individual single-cell BAMs

**Implementation**:
```bash
#!/bin/bash
#SBATCH -J sc_preprocess
#SBATCH --output=SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 4:00:00

set -euo pipefail

SAMPLE_LIST="$1"
BAM_DIR="$2"
OUTPUT_DIR="$3"
REFERENCE_DIR="$4"
SCRIPTS_DIR="$5"

# Get cell barcode from array task ID
BARCODE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

# Input/output files
INPUT_BAM="${BAM_DIR}/${BARCODE}.bam"
MARKDUP_BAM="${OUTPUT_DIR}/${BARCODE}.markdup.bam"
MARKDUP_METRICS="${OUTPUT_DIR}/${BARCODE}.markdup_metrics.txt"
RECAL_TABLE="${OUTPUT_DIR}/${BARCODE}.recal_data.table"
FINAL_BAM="${OUTPUT_DIR}/${BARCODE}.preprocessed.bam"

REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"

echo "Preprocessing cell: ${BARCODE}"

# Detect known sites for BQSR
source "${SCRIPTS_DIR}/scripts/utils/detect_known_sites.sh"
detect_known_sites "${REFERENCE_DIR}"

# Mark duplicates
module load picard
java -Xmx6g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I="${INPUT_BAM}" \
    O="${MARKDUP_BAM}" \
    M="${MARKDUP_METRICS}" \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=true

# BQSR (if known sites available)
if [ -n "${KNOWN_SITES_ARGS}" ]; then
    module load GATK

    gatk BaseRecalibrator \
        -R "${REFERENCE}" \
        -I "${MARKDUP_BAM}" \
        ${KNOWN_SITES_ARGS} \
        -O "${RECAL_TABLE}"

    gatk ApplyBQSR \
        -R "${REFERENCE}" \
        -I "${MARKDUP_BAM}" \
        --bqsr-recal-file "${RECAL_TABLE}" \
        -O "${FINAL_BAM}"

    # Remove intermediate markdup BAM
    rm "${MARKDUP_BAM}" "${MARKDUP_BAM}.bai"
else
    # No BQSR, just rename markdup BAM
    mv "${MARKDUP_BAM}" "${FINAL_BAM}"
    mv "${MARKDUP_BAM}.bai" "${FINAL_BAM}.bai"
fi

echo "Preprocessing complete for ${BARCODE}"
```

#### 2.2: Create Preprocessing Submission Script
**New file**: `scripts/CapWGS/submit_sc_preprocessing.sh`

**Purpose**: Submit preprocessing array after single-cell extraction

**Implementation**:
```bash
#!/bin/bash
#SBATCH -J submit_preprocess
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH -t 10

set -euo pipefail

REAL_CELLS_FILE="$1"
RAW_BAM_DIR="$2"
PREPROCESSED_BAM_DIR="$3"
REFERENCE_DIR="$4"
SCRIPTS_DIR="$5"

# Count cells
cell_count=$(wc -l < "$REAL_CELLS_FILE")

if [ "$cell_count" -eq 0 ]; then
    echo "ERROR: No cells detected"
    exit 1
fi

mkdir -p "${PREPROCESSED_BAM_DIR}"

echo "Submitting preprocessing array for ${cell_count} cells..."

preprocess_array_ID=$(sbatch --parsable \
    --array=1-${cell_count} \
    "${SCRIPTS_DIR}/scripts/CapWGS/sc_preprocessing_array.sh" \
    "${REAL_CELLS_FILE}" \
    "${RAW_BAM_DIR}" \
    "${PREPROCESSED_BAM_DIR}" \
    "${REFERENCE_DIR}" \
    "${SCRIPTS_DIR}")

echo "Preprocessing array submitted: ${preprocess_array_ID}"
echo "Array size: ${cell_count}"
```

#### 2.3: Update CapWGS_PP.sh Pipeline Logic
**File**: `CapWGS_PP.sh`

**Changes**:

1. **Remove bulk preprocessing section** (lines 199-211):
   - Delete `markdup_bqsr_job_ID` submission
   - Delete conditional `preprocessing_dependency`

2. **Add single-cell preprocessing** (after line 233, before QC):
```bash
######################################################################################################
#### Preprocessing for GATK mode (MarkDuplicates + BQSR on single cells)

if [ "$VARIANT_CALLER" = "gatk" ]; then
    # Create directory for preprocessed single-cell BAMs
    PREPROCESSED_SC_DIR="${SC_OUTPUTS_DIR}/preprocessed"

    # Run MarkDuplicates and BQSR on each single cell in parallel
    sc_preprocessing_job_ID=$(sbatch --parsable \
        --dependency=afterok:$sc_extraction_job_ID \
        "${SCRIPTS_DIR}/scripts/CapWGS/submit_sc_preprocessing.sh" \
        "${RESULTS_DIR}/real_cells.txt" \
        "${SC_OUTPUTS_DIR}" \
        "${PREPROCESSED_SC_DIR}" \
        "${REFERENCE_GENOME}" \
        "${SCRIPTS_DIR}")

    echo "Single-cell preprocessing array job ID: ${sc_preprocessing_job_ID}"

    # Set BAM directory for downstream steps
    FINAL_SC_BAM_DIR="${PREPROCESSED_SC_DIR}"
    preprocessing_dependency=$sc_preprocessing_job_ID
else
    # bcftools mode: no preprocessing needed
    FINAL_SC_BAM_DIR="${SC_OUTPUTS_DIR}"
    preprocessing_dependency=$sc_extraction_job_ID
fi
```

3. **Update QC and variant calling dependencies**:
   - QC jobs depend on `$preprocessing_dependency`
   - Variant calling uses `$FINAL_SC_BAM_DIR`

4. **Remove second sc_extraction** (lines 216-233):
   - No longer need to extract from preprocessed bulk BAM
   - Single extraction happens before preprocessing

#### 2.4: Remove Obsolete Script
**File to remove**: `scripts/CapWGS/markdup_bqsr.sh`

This script is replaced by the array-based approach.

#### 2.5: Update Variant Calling Scripts
**Files**: `scripts/CapWGS/submit_sc_variant_calling.sh`

**Changes**:
Ensure it uses `$FINAL_SC_BAM_DIR` for BAM location (passed from CapWGS_PP.sh)

---

## TODO #3: Parallelize Joint Calling with ~100MB Chunks

### Goal
Replace chromosome-level parallelization with interval-based parallelization using ~100 genomic chunks for faster joint calling. This avoids timeouts on large chromosomes (especially chr1).

### Approach
Use GATK's `SplitIntervals` tool with `--scatter-count 100` to generate ~100 balanced genomic intervals of roughly equal size (~30MB each for human genome).

### Components

#### 3.1: Create Interval Generation Script
**New file**: `scripts/utils/generate_genomic_intervals.sh`

**Purpose**: Generate balanced genomic intervals using GATK SplitIntervals

**Implementation**:
```bash
#!/bin/bash
#SBATCH -J split_intervals
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 1:00:00

set -euo pipefail

REFERENCE_DIR="$1"
INTERVAL_LIST="$2"     # Output: list of interval file paths
TMP_DIR="$3"           # Temp location for interval files
SCATTER_COUNT="${4:-100}"  # Number of chunks (default: 100)

REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
CALLING_REGIONS="${REFERENCE_DIR}/wgs_calling_regions.hg38.interval_list"
INTERVALS_DIR="${TMP_DIR}/intervals"

mkdir -p "${INTERVALS_DIR}"

module load GATK

echo "Generating ${SCATTER_COUNT} balanced genomic intervals..."
echo "Reference: ${REFERENCE}"
echo "Calling regions: ${CALLING_REGIONS}"
echo "Output directory: ${INTERVALS_DIR}"

gatk SplitIntervals \
    -R "${REFERENCE}" \
    -L "${CALLING_REGIONS}" \
    --scatter-count "${SCATTER_COUNT}" \
    --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
    -O "${INTERVALS_DIR}"

# Create list of interval files (sorted by genomic position)
ls "${INTERVALS_DIR}"/*.interval_list | sort -V > "${INTERVAL_LIST}"

N_INTERVALS=$(wc -l < "${INTERVAL_LIST}")

echo ""
echo "Generated ${N_INTERVALS} interval files"
echo "Interval list: ${INTERVAL_LIST}"
echo ""
echo "First few intervals:"
head -5 "${INTERVAL_LIST}"
```

**Note**: `BALANCING_WITHOUT_INTERVAL_SUBDIVISION` creates chunks with similar base pair counts without splitting individual callable regions.

#### 3.2: Create Interval-Based GenomicsDB Import Array
**New file**: `scripts/CapWGS/gatk_genomicsdb_import_interval_array.sh`

**Purpose**: Import GVCFs to GenomicsDB for one genomic interval

**Implementation**:
```bash
#!/bin/bash
#SBATCH -J genomicsdb_interval
#SBATCH --output=SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 12:00:00

set -euo pipefail

GVCF_DIR="$1"
GENOMICSDB_BASE_DIR="$2"
SAMPLE_MAP="$3"
REFERENCE_DIR="$4"
INTERVAL_LIST="$5"

# Get interval file for this task
INTERVAL_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${INTERVAL_LIST}")
INTERVAL_NAME=$(basename "${INTERVAL_FILE}" .interval_list)

REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
GENOMICSDB_DIR="${GENOMICSDB_BASE_DIR}/${INTERVAL_NAME}"
TMP_DIR="${GENOMICSDB_BASE_DIR}/tmp_${INTERVAL_NAME}"

mkdir -p "${TMP_DIR}"

echo "GenomicsDB import for interval: ${INTERVAL_NAME}"
echo "Interval file: ${INTERVAL_FILE}"
echo "Output: ${GENOMICSDB_DIR}"

module load GATK

gatk --java-options "-Xmx28g" GenomicsDBImport \
    --sample-name-map "${SAMPLE_MAP}" \
    --genomicsdb-workspace-path "${GENOMICSDB_DIR}" \
    -L "${INTERVAL_FILE}" \
    --tmp-dir "${TMP_DIR}" \
    --reader-threads 4

rm -rf "${TMP_DIR}"

echo "GenomicsDB import complete for ${INTERVAL_NAME}"
```

#### 3.3: Create Interval-Based Joint Calling Array
**New file**: `scripts/CapWGS/gatk_joint_calling_interval_array.sh`

**Purpose**: Run GenotypeGVCFs on one genomic interval

**Implementation**:
```bash
#!/bin/bash
#SBATCH -J joint_call_interval
#SBATCH --output=SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 2-00:00:00

set -euo pipefail

GENOMICSDB_BASE_DIR="$1"
OUTPUT_DIR="$2"
REFERENCE_DIR="$3"
INTERVAL_LIST="$4"

# Get interval file for this task
INTERVAL_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${INTERVAL_LIST}")
INTERVAL_NAME=$(basename "${INTERVAL_FILE}" .interval_list)

REFERENCE="${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta"
GENOMICSDB_DIR="${GENOMICSDB_BASE_DIR}/${INTERVAL_NAME}"
OUTPUT_VCF="${OUTPUT_DIR}/${INTERVAL_NAME}.vcf.gz"
TMP_DIR="${OUTPUT_DIR}/tmp_${INTERVAL_NAME}"

mkdir -p "${TMP_DIR}"

echo "Joint calling interval: ${INTERVAL_NAME}"
echo "GenomicsDB: ${GENOMICSDB_DIR}"
echo "Output VCF: ${OUTPUT_VCF}"

# Verify GenomicsDB exists
if [ ! -f "${GENOMICSDB_DIR}/callset.json" ]; then
    echo "ERROR: GenomicsDB not found for ${INTERVAL_NAME}"
    exit 1
fi

module load GATK

gatk --java-options "-Xmx28g" GenotypeGVCFs \
    -R "${REFERENCE}" \
    -V gendb://"${GENOMICSDB_DIR}" \
    -L "${INTERVAL_FILE}" \
    -O "${OUTPUT_VCF}" \
    --tmp-dir "${TMP_DIR}"

rm -rf "${TMP_DIR}"

echo "Joint calling complete for ${INTERVAL_NAME}"
```

#### 3.4: Update Merge Script for Intervals
**File**: `scripts/CapWGS/gatk_joint_calling_parallel_merge.sh`

**Changes**:
Update to merge interval-based VCFs (instead of chromosome-based):

```bash
#!/bin/bash
#SBATCH -J merge_joint_vcfs
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 12:00:00

set -euo pipefail

PER_INTERVAL_DIR="$1"
OUTPUT_VCF="$2"
INTERVAL_LIST="$3"

echo "Merging interval VCFs into final joint VCF..."
echo "Input directory: ${PER_INTERVAL_DIR}"
echo "Output VCF: ${OUTPUT_VCF}"

# Create input list for GatherVcfs
VCF_LIST="${PER_INTERVAL_DIR}/vcf_list.txt"
> "${VCF_LIST}"

while read interval_file; do
    interval_name=$(basename "${interval_file}" .interval_list)
    vcf_file="${PER_INTERVAL_DIR}/${interval_name}.vcf.gz"

    if [ ! -f "${vcf_file}" ]; then
        echo "ERROR: Missing VCF for interval ${interval_name}"
        exit 1
    fi

    echo "${vcf_file}" >> "${VCF_LIST}"
done < "${INTERVAL_LIST}"

n_vcfs=$(wc -l < "${VCF_LIST}")
echo "Found ${n_vcfs} interval VCFs to merge"

module load picard

# Use Picard GatherVcfs (preserves genomic order)
java -Xmx28g -jar $EBROOTPICARD/picard.jar GatherVcfs \
    I="${VCF_LIST}" \
    O="${OUTPUT_VCF}"

echo "Merge complete: ${OUTPUT_VCF}"

# Generate statistics
module load BCFtools/1.18-GCC-12.2.0

bcftools stats "${OUTPUT_VCF}" > "${OUTPUT_VCF%.vcf.gz}_stats.txt"

echo "Variant statistics:"
bcftools view -H "${OUTPUT_VCF}" | wc -l | awk '{print "Total variant sites: " $1}'
bcftools view -H -v snps "${OUTPUT_VCF}" | wc -l | awk '{print "SNP sites: " $1}'
bcftools view -H -v indels "${OUTPUT_VCF}" | wc -l | awk '{print "Indel sites: " $1}'

echo "Merge complete!"
```

#### 3.5: Update CapWGS submit_joint_calling.sh
**File**: `scripts/CapWGS/submit_joint_calling.sh`

**Changes** (GATK mode section, lines 62-192):

Replace chromosome-based logic with interval-based:

```bash
elif [ "$VARIANT_CALLER" = "gatk" ]; then
    ##################################################################################################################
    # GATK mode: Interval-based GenomicsDBImport + GenotypeGVCFs + merge
    ##################################################################################################################

    echo "GATK mode: Running parallelized joint calling with genomic intervals..."
    echo ""

    # Define paths
    GENOMICSDB_DIR="${RESULTS_DIR}/genomicsdb"
    SAMPLE_MAP="${RESULTS_DIR}/sample_map.txt"
    INTERVAL_LIST="${TMP_DIR}/interval_list.txt"
    PER_INTERVAL_DIR="${TMP_DIR}/per_interval_vcfs"
    FINAL_VCF="${RESULTS_DIR}/${OUTPUT_NAME}_joint.vcf.gz"

    mkdir -p "${GENOMICSDB_DIR}"
    mkdir -p "${PER_INTERVAL_DIR}"
    mkdir -p "${TMP_DIR}"

    ##############################################################################################################
    # Step 1: Generate genomic intervals
    ##############################################################################################################

    echo "Step 1: Generating ~100 balanced genomic intervals..."

    interval_job_ID=$(sbatch --parsable \
        "${SCRIPTS_DIR}/scripts/utils/generate_genomic_intervals.sh" \
        "${REFERENCE_GENOME}" \
        "${INTERVAL_LIST}" \
        "${TMP_DIR}" \
        100)

    echo "Interval generation job submitted: ${interval_job_ID}"
    echo ""

    ##############################################################################################################
    # Step 2: Create sample map from GVCFs
    ##############################################################################################################

    echo "Step 2: Creating sample map from GVCFs..."

    > "${SAMPLE_MAP}"
    for gvcf in "${ALIGNED_DIR}"/*.g.vcf.gz; do
        if [ -f "$gvcf" ]; then
            sample_name=$(basename "$gvcf" .g.vcf.gz)
            echo -e "${sample_name}\t${gvcf}" >> "${SAMPLE_MAP}"
        fi
    done

    n_cells=$(wc -l < "${SAMPLE_MAP}")
    if [ "$n_cells" -eq 0 ]; then
        echo "ERROR: No GVCF files found in ${ALIGNED_DIR}"
        exit 1
    fi

    echo "Found ${n_cells} cells"
    echo ""

    ##############################################################################################################
    # Step 3: Submit interval-based GenomicsDB import array (waits for intervals)
    ##############################################################################################################

    echo "Step 3: Submitting interval-based GenomicsDB import array..."

    # This wrapper script reads interval_list.txt after it's created
    genomicsdb_array_ID=$(sbatch --parsable \
        --dependency=afterok:${interval_job_ID} \
        "${SCRIPTS_DIR}/scripts/CapWGS/submit_genomicsdb_import_array.sh" \
        "${ALIGNED_DIR}" \
        "${GENOMICSDB_DIR}" \
        "${SAMPLE_MAP}" \
        "${REFERENCE_GENOME}" \
        "${INTERVAL_LIST}")

    echo "GenomicsDB import array submitted: ${genomicsdb_array_ID}"
    echo ""

    ##############################################################################################################
    # Step 4: Submit interval-based joint calling array
    ##############################################################################################################

    echo "Step 4: Submitting interval-based joint calling array..."

    jc_array_ID=$(sbatch --parsable \
        --dependency=afterok:${genomicsdb_array_ID} \
        "${SCRIPTS_DIR}/scripts/CapWGS/submit_joint_calling_array.sh" \
        "${GENOMICSDB_DIR}" \
        "${PER_INTERVAL_DIR}" \
        "${REFERENCE_GENOME}" \
        "${INTERVAL_LIST}")

    echo "Joint calling array submitted: ${jc_array_ID}"
    echo ""

    ##############################################################################################################
    # Step 5: Submit merge job
    ##############################################################################################################

    echo "Step 5: Submitting merge job..."

    merge_job_ID=$(sbatch --parsable \
        --dependency=afterok:${jc_array_ID} \
        "${SCRIPTS_DIR}/scripts/CapWGS/gatk_joint_calling_parallel_merge.sh" \
        "${PER_INTERVAL_DIR}" \
        "${FINAL_VCF}" \
        "${INTERVAL_LIST}")

    echo "Merge job submitted: ${merge_job_ID}"
    echo ""

    ##############################################################################################################
    # Summary
    ##############################################################################################################

    echo "GATK joint calling pipeline submitted!"
    echo ""
    echo "Job chain:"
    echo "  1. Generate intervals: ${interval_job_ID}"
    echo "  2. GenomicsDB import: ${genomicsdb_array_ID} (~100 intervals)"
    echo "  3. Joint calling: ${jc_array_ID} (~100 intervals)"
    echo "  4. Merge VCFs: ${merge_job_ID}"
    echo ""
    echo "Output files:"
    echo "  GenomicsDB (per-interval): ${GENOMICSDB_DIR}/"
    echo "  Per-interval VCFs (temp): ${PER_INTERVAL_DIR}/"
    echo "  Final merged VCF: ${FINAL_VCF}"
    echo ""
fi
```

#### 3.6: Create Wrapper Scripts for Array Submission
**New file**: `scripts/CapWGS/submit_genomicsdb_import_array.sh`

**Purpose**: Submit GenomicsDB import array after determining interval count

```bash
#!/bin/bash
#SBATCH -J submit_genomicsdb
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH -t 10

set -euo pipefail

GVCF_DIR="$1"
GENOMICSDB_DIR="$2"
SAMPLE_MAP="$3"
REFERENCE_DIR="$4"
INTERVAL_LIST="$5"

# Count intervals
n_intervals=$(wc -l < "${INTERVAL_LIST}")

echo "Submitting GenomicsDB import array for ${n_intervals} intervals..."

genomicsdb_array_ID=$(sbatch --parsable \
    --array=1-${n_intervals} \
    ./scripts/CapWGS/gatk_genomicsdb_import_interval_array.sh \
    "${GVCF_DIR}" \
    "${GENOMICSDB_DIR}" \
    "${SAMPLE_MAP}" \
    "${REFERENCE_DIR}" \
    "${INTERVAL_LIST}")

echo "GenomicsDB import array submitted: ${genomicsdb_array_ID}"
```

**New file**: `scripts/CapWGS/submit_joint_calling_array.sh`

**Purpose**: Submit joint calling array after determining interval count

```bash
#!/bin/bash
#SBATCH -J submit_joint_call
#SBATCH -o SLURM_outs/%x_%j.out
#SBATCH -c 1
#SBATCH -t 10

set -euo pipefail

GENOMICSDB_DIR="$1"
OUTPUT_DIR="$2"
REFERENCE_DIR="$3"
INTERVAL_LIST="$4"

# Count intervals
n_intervals=$(wc -l < "${INTERVAL_LIST}")

echo "Submitting joint calling array for ${n_intervals} intervals..."

jc_array_ID=$(sbatch --parsable \
    --array=1-${n_intervals} \
    ./scripts/CapWGS/gatk_joint_calling_interval_array.sh \
    "${GENOMICSDB_DIR}" \
    "${OUTPUT_DIR}" \
    "${REFERENCE_DIR}" \
    "${INTERVAL_LIST}")

echo "Joint calling array submitted: ${jc_array_ID}"
```

#### 3.7: Update Bulk GATK Pipeline
**File**: `scripts/bulk/gatk_pipeline.sh`

**Changes**:
Apply the same interval-based approach:

1. Replace chromosome list generation (lines 167-193) with interval generation
2. Update GenomicsDB import to use interval arrays
3. Update joint calling to use interval arrays
4. Update merge to handle interval-based VCFs
5. Copy/adapt the interval-based array scripts to `scripts/bulk/`

**Note**: `gatk_single_sample.sh` doesn't need changes (single-sample calling, no joint calling)

---

## Implementation Order

### Phase 1: QC Integration (TODO #1)
1. Update `concatenate.sh` with barcode metrics
2. Add QC job submissions to `CapWGS_PP.sh`
3. Create `compile_qc_results.sh` script
4. Update `submit_benchmarking_qc.sh` parameters
5. Test on small dataset

### Phase 2: Architecture Change (TODO #2)
1. Create `sc_preprocessing_array.sh`
2. Create `submit_sc_preprocessing.sh`
3. Update `CapWGS_PP.sh` pipeline logic
4. Remove `markdup_bqsr.sh` (obsolete)
5. Update variant calling scripts to use preprocessed BAMs
6. Test on small dataset

### Phase 3: Interval Parallelization (TODO #3)
1. Create `generate_genomic_intervals.sh`
2. Create interval-based GenomicsDB import array script
3. Create interval-based joint calling array script
4. Create wrapper submission scripts
5. Update merge script for intervals
6. Update `submit_joint_calling.sh`
7. Apply same changes to `scripts/bulk/gatk_pipeline.sh`
8. Test on small dataset

---

## Testing Strategy

### For TODO #1 (QC Integration):
1. Run `CapWGS_PP.sh` on small test dataset
2. Verify barcode assignment statistics in concatenation log
3. Verify bigwig/Lorenz jobs submit and complete
4. Verify benchmarking QC jobs submit and complete
5. Verify compiled QC CSV files appear in results directory
6. Test with both GATK and bcftools modes

### For TODO #2 (Architecture Change):
1. Run GATK mode on small dataset
2. Verify preprocessing array submits after cell extraction
3. Verify preprocessed BAMs created in correct directory
4. Verify QC runs on preprocessed BAMs
5. Verify variant calling uses preprocessed BAMs
6. Compare runtime: new approach should be much faster
7. Verify bcftools mode unchanged (no preprocessing)

### For TODO #3 (Interval Parallelization):
1. Test interval generation on reference genome
2. Verify ~100 intervals created
3. Run small pilot with few cells
4. Monitor GenomicsDB import array (should complete much faster)
5. Monitor joint calling array (should complete much faster)
6. Verify merged VCF correctness and completeness
7. Compare total runtime: should be significantly faster than chromosome-level
8. Test on bulk pipeline as well

---

## Expected Performance Improvements

### TODO #1 (QC Integration):
- **No performance impact**: QC runs in parallel with variant calling
- **Benefit**: Comprehensive QC metrics without separate pipeline runs

### TODO #2 (Architecture Change):
- **Current**: 1 slow bulk markdup/BQSR job (hours to days for large datasets)
- **New**: N parallel single-cell jobs (minutes to hours)
- **Expected speedup**: 5-10x faster preprocessing step

### TODO #3 (Interval Parallelization):
- **Current**: Chromosome-level (25 jobs, chr1 can take days)
- **New**: Interval-level (~100 jobs, each much smaller)
- **Expected speedup**: 3-5x faster joint calling
- **Reason**: Better load balancing, smaller chunks, avoids chr1 bottleneck

---

## Files to Create

### TODO #1:
- `scripts/CapWGS_QC/compile_qc_results.sh`

### TODO #2:
- `scripts/CapWGS/sc_preprocessing_array.sh`
- `scripts/CapWGS/submit_sc_preprocessing.sh`

### TODO #3:
- `scripts/utils/generate_genomic_intervals.sh`
- `scripts/CapWGS/gatk_genomicsdb_import_interval_array.sh`
- `scripts/CapWGS/gatk_joint_calling_interval_array.sh`
- `scripts/CapWGS/submit_genomicsdb_import_array.sh`
- `scripts/CapWGS/submit_joint_calling_array.sh`
- `scripts/bulk/gatk_genomicsdb_import_interval_array.sh` (copy/adapt)
- `scripts/bulk/gatk_joint_calling_interval_array.sh` (copy/adapt)
- `scripts/bulk/submit_genomicsdb_import_array.sh` (copy/adapt)
- `scripts/bulk/submit_joint_calling_array.sh` (copy/adapt)

## Files to Modify

### TODO #1:
- `scripts/CapWGS/concatenate.sh`
- `CapWGS_PP.sh`
- `scripts/CapWGS_QC/submit_benchmarking_qc.sh`

### TODO #2:
- `CapWGS_PP.sh` (major refactor of GATK pipeline logic)
- `scripts/CapWGS/submit_sc_variant_calling.sh`

### TODO #3:
- `scripts/CapWGS/submit_joint_calling.sh` (major refactor)
- `scripts/CapWGS/gatk_joint_calling_parallel_merge.sh`
- `scripts/bulk/gatk_pipeline.sh` (major refactor)

## Files to Remove

### TODO #2:
- `scripts/CapWGS/markdup_bqsr.sh` (replaced by array-based preprocessing)

---

## Notes

- All new scripts should follow existing conventions (SLURM headers, error handling, logging)
- Intermediate files (per-interval VCFs, GenomicsDB) stored in user-provided temp directory
- Compiled QC results saved to results directory for easy access
- Pipeline maintains backward compatibility for bcftools mode (no preprocessing)
- Interval generation is deterministic (same intervals for same reference) 