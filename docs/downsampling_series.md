# Downsampling Series Coverage Analysis

## Overview

Compare CapWGS single-cell data to publicly available single-cell WGA methods by analyzing coverage uniformity across different sequencing depths. The pipeline downsamples BAM files to standardized read depths (1M, 5M, 10M, 50M, 100M), collects WGS metrics, and generates Lander-Waterman saturation curves.

## Design Decisions

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| **Random seed** | 42 | Fixed for reproducibility; ensures nested subsampling (1M ⊂ 5M ⊂ 10M...) |
| **Target depths** | 1M, 5M, 10M, 50M, 100M | Standard benchmarking range for single-cell WGA |
| **Input BAMs** | `*.bam` (NOT `*.preprocessed.bam`) | Non-deduplicated, consistent with public data |
| **Samtools flags** | `-F 260` | Primary mapped reads only (excludes unmapped + secondary) |
| **Reference detection** | From BAM `@SQ` header | Automatic, flexible to any genome |
| **Cell selection** | User-provided BAM list | Explicit control over which cells to analyze |

## Directory Structure

```
data/benchmarking/
├── manifests/
│   ├── input_bams.txt              # User-provided: one BAM path per line
│   ├── bam_manifest.tsv            # Processed: dataset, cell_id, bam_path, reference, total_mapped_reads
│   └── downsample_tasks.tsv        # Expanded: all (cell × depth) combinations with fractions
├── downsampled/
│   └── {dataset}/
│       └── {cell_id}_{depth}.bam   # e.g., AAAA-BBBB-CCCC-DDDD_5M.bam (+ .bai)
└── qc_metrics/
    └── {dataset}/
        └── {cell_id}_{depth}_wgs_metrics.txt

results/benchmarking/
├── compiled_metrics.csv            # Long-format: dataset, cell, depth, all Picard metrics
└── downsampling_analysis.ipynb     # Lander-Waterman fitting & plots
```

## Scripts

### 1. Manifest Builder
**Path:** `scripts/benchmarking/create_bam_manifest.sh`

**Purpose:** Process user-provided BAM list and extract metadata.

**Input:**
- `data/benchmarking/manifests/input_bams.txt` (one BAM path per line)

**Output:**
- `data/benchmarking/manifests/bam_manifest.tsv`

**Columns:**
- `dataset`: Inferred from parent directory (e.g., `HSC4`, `public_WGA`)
- `cell_id`: Basename without `.bam` extension
- `bam_path`: Full absolute path
- `reference`: Extracted from BAM `@SQ SN:` header lines
- `total_mapped_reads`: From `samtools view -c -F 260` (primary mapped only)

**Method:**
- For each BAM in input list:
  - Extract reference genome path from BAM header using `samtools view -H`
  - Count primary mapped reads with `samtools view -c -F 260`
  - Parse dataset/cell_id from file path

### 2. Task Expander
**Path:** `scripts/benchmarking/create_downsample_tasks.py`

**Purpose:** Generate downsampling task list with calculated fractions.

**Input:**
- `data/benchmarking/manifests/bam_manifest.tsv`
- Target depths: 1M, 5M, 10M, 50M, 100M

**Output:**
- `data/benchmarking/manifests/downsample_tasks.tsv`

**Columns:**
- `task_id`: Sequential integer for SLURM array indexing
- `dataset`: From manifest
- `cell_id`: From manifest
- `input_bam`: Full path to source BAM
- `reference`: Reference genome path
- `target_depth`: Target read count (e.g., 5000000)
- `depth_label`: Human-readable (e.g., "5M")
- `fraction`: `target_depth / total_mapped_reads` (capped at 1.0)
- `output_bam`: Path to downsampled BAM

**Logic:**
- Skip tasks where `target_depth > total_mapped_reads` (insufficient reads)
- For remaining tasks, calculate exact fraction for `samtools view -s`

### 3. Downsampling Array
**Path:** `scripts/benchmarking/downsample_array.sh`

**Purpose:** SLURM array job to downsample BAMs.

**Input:**
- Task file: `data/benchmarking/manifests/downsample_tasks.tsv`
- Array task ID: `$SLURM_ARRAY_TASK_ID`

**Command:**
```bash
samtools view -s 42.${FRACTION} -b -o ${OUTPUT_BAM} ${INPUT_BAM}
samtools index ${OUTPUT_BAM}
```

**Resources:**
- CPUs: 1
- Memory: 8GB
- Time: 2 hours

### 4. CollectWgsMetrics Array
**Path:** `scripts/benchmarking/collect_wgs_metrics_array.sh`

**Purpose:** SLURM array job to collect Picard WGS metrics on downsampled BAMs.

**Input:**
- Task file: `data/benchmarking/manifests/downsample_tasks.tsv`
- Array task ID: `$SLURM_ARRAY_TASK_ID`

**Command:**
```bash
module load picard
java -jar $EBROOTPICARD/picard.jar CollectWgsMetrics \
    I=${DOWNSAMPLED_BAM} \
    O=${WGS_METRICS} \
    R=${REFERENCE}
```

**Output:**
- `data/benchmarking/qc_metrics/{dataset}/{cell_id}_{depth}_wgs_metrics.txt`

**Resources:**
- CPUs: 1
- Memory: 16GB
- Time: 4 hours

**Borrowed from:** `scripts/CapWGS/sc_extract_preprocess_qc_array.sh:246-250`

### 5. Metrics Compiler
**Path:** `scripts/benchmarking/compile_metrics.py`

**Purpose:** Parse all Picard WGS metrics into single CSV.

**Input:**
- Directory: `data/benchmarking/qc_metrics/`
- Pattern: `**/*_wgs_metrics.txt`

**Output:**
- `results/benchmarking/compiled_metrics.csv`

**Key Columns:**
- Metadata: `dataset`, `cell_id`, `target_depth`, `depth_label`
- Coverage: `MEAN_COVERAGE`, `MEDIAN_COVERAGE`, `MAD_COVERAGE`
- Breadth: `PCT_1X`, `PCT_5X`, `PCT_10X`, `PCT_20X`, `PCT_30X`
- Exclusions: `PCT_EXC_DUPE`, `PCT_EXC_MAPQ`, `PCT_EXC_BASEQ`, `PCT_EXC_TOTAL`
- Other: `GENOME_TERRITORY`, `SD_COVERAGE`

**Method:**
- For each metrics file:
  - Skip comment lines (starting with `#`)
  - Parse header and data rows
  - Extract dataset/cell/depth from filename
  - Append to master DataFrame

### 6. Analysis Notebook
**Path:** `results/benchmarking/downsampling_analysis.ipynb`

**Purpose:** Fit Lander-Waterman curves and generate comparison plots.

**Analyses:**
1. Per-cell saturation curves (coverage vs. reads)
2. Dataset-level coverage uniformity (Gini, Lorenz)
3. Method comparison (CapWGS vs. public methods)
4. Depth recommendations (reads needed for target coverage)

**Dependencies:**
- `cellspec` (installed in micromamba environment)
- pandas, matplotlib, seaborn, scipy

## Workflow Execution

### Step 1: Create Input BAM List
```bash
# User manually creates this file with desired BAMs
# Example for HSC4:
ls /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/data/HSC4/sc_outputs/*.bam | \
  grep -v "preprocessed" > data/benchmarking/manifests/input_bams.txt

# Add public data:
ls /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome/data/public_WGA/*.bam >> \
  data/benchmarking/manifests/input_bams.txt
```

### Step 2: Build Manifest
```bash
bash scripts/benchmarking/create_bam_manifest.sh \
    data/benchmarking/manifests/input_bams.txt \
    data/benchmarking/manifests/bam_manifest.tsv
```

### Step 3: Generate Downsampling Tasks
```bash
python scripts/benchmarking/create_downsample_tasks.py \
    --input data/benchmarking/manifests/bam_manifest.tsv \
    --depths 1000000,5000000,10000000,50000000,100000000 \
    --output data/benchmarking/manifests/downsample_tasks.tsv
```

### Step 4: Submit Downsampling Array
```bash
NTASKS=$(tail -n +2 data/benchmarking/manifests/downsample_tasks.tsv | wc -l)
JOB1=$(sbatch --array=1-${NTASKS} \
    scripts/benchmarking/downsample_array.sh \
    data/benchmarking/manifests/downsample_tasks.tsv)
echo "Downsampling job ID: ${JOB1##* }"
```

### Step 5: Submit CollectWgsMetrics Array
```bash
JOB2=$(sbatch --dependency=afterok:${JOB1##* } \
    --array=1-${NTASKS} \
    scripts/benchmarking/collect_wgs_metrics_array.sh \
    data/benchmarking/manifests/downsample_tasks.tsv)
echo "WGS metrics job ID: ${JOB2##* }"
```

### Step 6: Compile Metrics
```bash
# Run after all array tasks complete
python scripts/benchmarking/compile_metrics.py \
    --metrics-dir data/benchmarking/qc_metrics \
    --output results/benchmarking/compiled_metrics.csv
```

### Step 7: Analysis
```bash
# Launch Jupyter and run downsampling_analysis.ipynb
jupyter notebook results/benchmarking/downsampling_analysis.ipynb
```

## Technical Notes

### Reference Genome Detection
The manifest builder extracts reference paths from BAM headers using:
```bash
samtools view -H ${BAM} | grep "^@SQ" | head -1 | awk '{print $2}' | cut -d: -f2
```

This assumes standard BAM headers with `@SQ SN:chr1 LN:248956422` format. The reference FASTA path is reconstructed from the sequence dictionary.

### Nested Downsampling
Using a fixed seed (42) ensures that reads selected at 1M depth are also present in the 5M sample. This creates nested subsets, making saturation curves monotonic and improving model fitting stability.

### Handling Low-Depth Cells
If a cell has fewer reads than a target depth (e.g., 500K reads but target is 1M), the task expander skips that (cell, depth) combination. This prevents `samtools view -s` errors and ensures `fraction ≤ 1.0`.

### Memory Requirements
- Downsampling: 8GB sufficient for most single-cell BAMs (<200M reads)
- CollectWgsMetrics: 16GB recommended for human genome references
- Scale memory if working with larger genomes or ultra-deep cells

## Expected Outputs

### Compiled Metrics CSV
~500-5000 rows (depending on # cells × depths), with columns:
- `dataset`, `cell_id`, `target_depth`, `depth_label`
- `MEAN_COVERAGE`, `MEDIAN_COVERAGE`, `SD_COVERAGE`, `MAD_COVERAGE`
- `PCT_1X`, `PCT_5X`, `PCT_10X`, `PCT_20X`, `PCT_30X`
- `PCT_EXC_DUPE`, `PCT_EXC_MAPQ`, `PCT_EXC_BASEQ`, `PCT_EXC_TOTAL`
- `GENOME_TERRITORY`

### Lander-Waterman Curves
Per-cell plots showing:
- X-axis: Total reads (1M to 100M)
- Y-axis: Mean coverage
- Fitted curve: `C = L * (1 - e^(-reads/genome_size))`
- Key parameter: Effective library complexity `L`

### Method Comparison
Box plots or violin plots showing:
- Coverage uniformity (Gini coefficient) per method
- Breadth of coverage (PCT_10X) per method
- Duplication rate vs. depth per method
