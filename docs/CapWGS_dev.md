# CapWGS Pipeline Development Notes

## Recent Improvements (April 2026)

### 5. HTML Dashboard Generation ✅ COMPLETED
**Goal**: Create self-contained interactive HTML dashboard for pipeline QC results

**Problem**: QC metrics spread across multiple CSV files and static PNG plots, difficult to share and explore

**Solution**: Single-file HTML dashboard with embedded Plotly.js visualizations
- Self-contained: all data, plots, and code in one HTML file
- Interactive charts: zoom, pan, hover tooltips on all plots
- Tabbed interface: Summary, Per-Cell QC, Coverage tabs
- Styled after 10x Genomics Cell Ranger web_summary.html

**Features**:
- **Summary Tab**:
  - Metric cards (cells, reads, mapping rates, duplication, etc.)
  - Read retention Sankey diagram (demux → trimming → alignment)
  - Preprocessing statistics table (with flagstat integration)
  - Adapter position histograms (R1 + R2 overlaid)
  - Interactive knee plot (colored by condition if provided)
- **Per-Cell QC Tab**:
  - Alignment rate, duplication rate, coverage distributions
  - Gini coefficient histogram
  - Coverage vs genome fraction scatter
  - All plots support condition-based coloring
- **Coverage Tab**:
  - Lorenz curves (up to 50 cells sampled)
  - Lander-Waterman plot (embedded PNG)

**Optional Features**:
- Condition assignment via JSON config + barcode plate map
- Color-coded plots by experimental condition
- Barcode-to-condition mapping for interactive exploration

**Files Created**:
- `scripts/CapWGS/generate_dashboard.py`: Dashboard generator (927 lines)
- `scripts/CapWGS/preprocessing_parsers.py`: Parse SLURM logs for stats
- `scripts/CapWGS/parse_slurm_logs.py`: Extract preprocessing stats from old runs

**Integration**:
- Reads existing pipeline outputs (compiled_qc_metrics.csv, lorenz curves, etc.)
- Uses static Lander-Waterman PNG from plot_lander_waterman.py
- Embeds all plots and data as base64/JSON in HTML
- Can be generated post-hoc for old runs by parsing SLURM logs

**Usage**:
```bash
# For new runs (with preprocessing_summary.json)
python3 scripts/CapWGS/generate_dashboard.py results/sample/

# For old runs (parse SLURM logs first)
python3 scripts/CapWGS/parse_slurm_logs.py SLURM_outs/array_outs results/sample/ --job-id 12345
python3 scripts/CapWGS/generate_dashboard.py results/sample/

# With experimental conditions
python3 scripts/CapWGS/generate_dashboard.py results/sample/ \
    --conditions conditions.json --barcode-map barcodes/map.csv
```

**Commit**: `adcc504` (Srivatsan, PI)

---

### 6. Lander-Waterman Plot Formula Fix ✅ COMPLETED
**Goal**: Fix theoretical curve falling 2 orders of magnitude below data points

**Problem**: Original implementation missing read length factor in coverage calculation
- Formula used: `coverage = reads_per_mb / 1e6` (assumed 1bp reads!)
- Curve appeared ~100x too low on log-scale plots

**Solution**: Include mean read length in Lander-Waterman formula
- Correct formula: `coverage = (reads_per_mb * read_length) / 1e6`
- Extracts mean read length from QC metrics (typically ~142bp for Illumina)
- Curve now properly overlays with empirical data

**Additional Improvements**:
- **Flexible input format**: Accept either genome length (integer) OR FASTA path
  - `plot_lander_waterman.py 12000000 ...` (simple, original behavior)
  - `plot_lander_waterman.py /path/to/genome.fa ...` (pipeline compatible)
  - Auto-detects format and parses accordingly
- **Simplified code**: Removed complex directory search logic (~25 lines)
  - Now expects explicit FASTA file path, not directory
  - Cleaner, more maintainable implementation
- **Better output**: Shows mean read length used in summary statistics

**Formula Details**:
```python
# Old (incorrect):
y_theory = 1 - np.exp(-reads_per_mb / 1e6)

# New (correct):
mean_read_length = df['mean_read_length'].mean()  # ~142 bp
y_theory = 1 - np.exp(-reads_per_mb * mean_read_length / 1e6)
```

**Files Modified**:
- `scripts/utils/plot_lander_waterman.py`: Fixed formula, simplified input handling

**Testing**:
- Regenerated plot for HSC4_pilot (500 cells, hg38 reference)
- Theoretical curve now aligns with data across 4 orders of magnitude
- Mean read length: 141.8 bp (from compiled QC metrics)

**Commit**: `79f79f3`

---

### 7. Barcode Annotation Refactoring ✅ COMPLETED
**Goal**: Eliminate code duplication between notebook and pipeline barcode annotation tools

**Problem**: Two scripts with ~150 lines of duplicated logic
- `scripts/utils/barcode_annotation.py`: DataFrame API for notebooks
- `scripts/CapWGS/assign_conditions.py`: CLI tool for dashboards
- Both implemented same barcode slicing, well parsing, lookup building

**Solution**: Layered architecture with shared core utilities

**Architecture**:
```
┌─────────────────────────────────────────────────────────┐
│  HIGH-LEVEL APIs (Use-case specific)                    │
├─────────────────────────────────────────────────────────┤
│  barcode_annotation.py   │  assign_conditions.py        │
│  (for notebooks)         │  (for pipeline/CLI)          │
├─────────────────────────────────────────────────────────┤
│  CORE LOGIC (Shared utilities)                          │
├─────────────────────────────────────────────────────────┤
│  barcode_core.py - Constants, plate loading, lookup     │
│                    building, condition assignment        │
└─────────────────────────────────────────────────────────┘
```

**Files Created**:
- `scripts/utils/barcode_core.py` (228 lines): Shared utilities
  - Constants: `BC_SLICES`, `POSITION_COLUMNS` for all 4 positions (A/B/C/D)
  - `load_barcode_map_dict/dataframe()`: Support both dict and DataFrame
  - `parse_well_spec()`: Flexible well parsing (rows, columns, specific wells)
  - `build_barcode_lookup()`: Condition → barcode subsequence mapping
  - `assign_conditions()`: Cell barcode → condition assignment
- `scripts/utils/BARCODE_ANNOTATION.md` (214 lines): Comprehensive docs
  - Architecture diagrams, usage examples for both APIs
  - Barcode structure reference, well specification formats
  - Design rationale, testing instructions, migration guide

**Files Refactored**:
- `scripts/utils/barcode_annotation.py`: Now imports from barcode_core
  - Added detailed examples and docstrings
  - Maintains all functionality (DataFrame in/out, flexible wells)
  - Cross-references CLI tool
- `scripts/CapWGS/assign_conditions.py`: Reduced from 169 → 109 lines (40% reduction!)
  - Now imports from barcode_core
  - Added `--position` flag for flexibility (not just position A)
  - Maintains CLI interface and dashboard compatibility

**Enhanced Features**:
- **Flexible well specification** (both APIs):
  - Row letters: `'A'`, `'B'` → all valid columns in that row
  - Column numbers: `1`, `2`, `3` → all rows in that column
  - Specific wells: `'A1'`, `'B3'` → single well
- **Multi-position support**: Can annotate by any barcode position (A/B/C/D)
- **Better error handling**: Descriptive exceptions (notebooks) or warnings (pipeline)
- **Backward compatible**: Existing notebooks and configs still work

**Benefits**:
- ✅ Eliminated ~150 lines of duplicated code (DRY principle)
- ✅ Single source of truth for barcode logic
- ✅ Both use cases well-served (interactive + automation)
- ✅ Easier testing and maintenance
- ✅ Better documented with comprehensive guide

**Net Change**: -40 lines of duplication, +442 lines of core utilities and docs

**Commit**: `fb3a85c`

---

## Improvements from March 2026

### 1. QC Pipeline Integration ✅ COMPLETED
**Goal**: Consolidate all QC metrics into single compiled outputs

**Changes Made**:
- Integrated Picard QC metrics, Lorenz curves, and gini coefficients into unified pipeline
- Added Lander-Waterman coverage plot (`scripts/utils/plot_lander_waterman.py`)
- All QC outputs now compile to `results/{sample}/compiled_qc_metrics.csv`
- Column naming standardized to use 'barcode' instead of 'sample'
- Added file existence checks to prevent misleading success messages
- QC metrics directory: `data/{sample}/qc_metrics/`

**Files Modified**:
- `CapWGS_PP.sh`: Added QC compilation step with proper output verification
- `scripts/CapWGS/sc_extract_preprocess_qc_array.sh`: Unified single-cell processing
- `scripts/CapWGS_QC/compile_qc_metrics.py`: Added gini coefficient support, changed column names
- `scripts/CapWGS_QC/lorenz.py`: Already had --gini-output flag, now used

**Files Created**:
- `scripts/utils/plot_lander_waterman.py`: Generates coverage vs sequencing depth plots

---

### 2. Reference Genome Path Flexibility ✅ COMPLETED
**Goal**: Handle non-standard reference genome filenames and directory structures

**Problem**: Pipeline failed on C. elegans data using `GCA_028201515.1_genomic.fna` instead of standard `genome.fa`

**Solution**: Implemented flexible FASTA file detection
- Accepts both file and directory paths as input
- Searches for any `.fa` or `.fna` file with BWA index (`.amb` file)
- Checks both root directory and `BWAIndex/` subdirectory
- Auto-detects FASTA file and constructs proper paths

**Files Modified**:
- `scripts/CapWGS/PP_array.sh`: Added flexible BWA index search
- `scripts/CapWGS/sc_extract_preprocess_qc_array.sh`: Added flexible reference path handling
- `scripts/utils/plot_lander_waterman.py`: Handles directory input for genome length calculation

---

### 3. SAM Merge Filtering ✅ COMPLETED
**Goal**: Remove improperly paired reads caused by FASTQ chunking artifacts

**Problem**: Duplicate SAM records from chunking caused issues with downstream tools

**Solution**: Added `-f 0x2` filtering during SAM merge
- Filters for properly paired reads only
- Removes ~4% of reads (chunking artifacts)
- Applied during bulk BAM creation step

**Files Modified**:
- `CapWGS_PP.sh`: Added `-f 0x2` flag to samtools merge pipeline

---

### 4. Parallel Two-Stage SAM/BAM Merge ✅ COMPLETED
**Goal**: Reduce 3-5 day merge bottleneck to ~3-6 hours (20-40x speedup)

**Problem**: Serial merge of 500 SAM chunks takes 3-5 days for large datasets (10TB, 25B reads)

**Solution**: Two-stage parallel merge using SLURM arrays
- **Stage 1**: Split SAMs into groups, merge in parallel (default: 20 SAMs per group)
- **Stage 2**: Merge intermediate BAMs into final sorted BAM
- Configurable group size via `MERGE_GROUP_SIZE` variable
- Filtering (`-f 0x2`) and RG preservation (`-c`) in both stages

**Implementation**:
- Stage 1 submits array job to short partition
- Each task merges one group of SAMs → intermediate BAM
- Stage 2 merges all intermediate BAMs → final BAM
- Pipeline waits for arrays using job polling

**Files Created**:
- `scripts/CapWGS/parallel_merge_array.sh`: Array script for group merging

**Files Modified**:
- `CapWGS_PP.sh`: Replaced serial merge with two-stage parallel approach

**Testing**:
- Tested with 501 yeast SAM files
- Created 26 parallel tasks (25 × 20 SAMs, 1 × 1 SAM remainder)
- 25/26 tasks completed successfully
- Final BAM: 1.2GB, 20.7M reads (validated with samtools quickcheck)
- 4% read reduction from proper-pair filtering as expected

---

## Known Issues

### None Currently

All previously identified issues have been resolved:
- ✅ Array indexing edge case (fixed in `parallel_merge_array.sh`)
- ✅ Lander-Waterman curve offset (fixed with read length factor)
- ✅ Barcode annotation code duplication (refactored into layered architecture)

---

## Pending Work

### TODO: Dashboard Enhancements
- Add Lander-Waterman plot as interactive Plotly chart (currently static PNG)
- Consider adding condition coloring to Lander-Waterman points
- Add filtering/selection controls for interactive plots
- Export dashboard data to JSON for custom analysis

### TODO: Test Variant Calling Features
The following features are implemented but not yet tested:
- Issue #8: Add `-v none` option for alignment-only mode
- Issue #9: Add `--use-existing-chunks` flag to skip FASTQ splitting

### TODO: Consider Future Optimizations
- Convert SAM to BAM before Stage 1 merge (minor speedup)
- Adjust `MERGE_GROUP_SIZE` based on dataset characteristics
- Profile parallel merge performance on NovaSeq scale datasets
- Investigate whether paired-end reads need 2× read length factor in Lander-Waterman formula

---

## Pipeline Architecture (Current)

### CapWGS_PP.sh Flow:
```
1. Parse arguments
2. Split FASTQs into chunks (N_CHUNKS, default 500)
3. Submit PP_array (alignment) → WAIT
4. Two-stage parallel SAM merge:
   - Create merge groups (default: 20 SAMs per group)
   - Submit parallel merge array → WAIT
   - Final merge of intermediate BAMs
5. Index BAM, compute readcounts
6. Detect cells (knee plot)
7. Submit single-cell extraction array → WAIT
8. Submit single-cell processing array (MarkDup, BQSR, bigwig, Lorenz, QC) → WAIT
9. Compile QC metrics (Lorenz curves, Picard metrics, Lander-Waterman plot)
10. Generate HTML dashboard (optional, via generate_dashboard.py)
11. Submit variant calling (if requested)
```

### Key Directories:
- Bulk BAM: `data/{sample}/{sample}.bam`
- Single-cell BAMs: `data/{sample}/sc_outputs/{barcode}.bam`
- QC metrics: `data/{sample}/qc_metrics/`
- Results: `results/{sample}/`
  - `compiled_qc_metrics.csv`: Per-cell QC metrics
  - `compiled_lorenz_curves.csv`: Coverage uniformity curves
  - `lander_waterman_coverage.png`: Coverage vs depth plot
  - `dashboard.html`: Interactive HTML QC dashboard
  - `real_cells.txt`: List of detected cell barcodes
  - `readcounts.csv`: Read counts per barcode
  - `kneeplot.png`: Barcode rank plot
- Temp files: `/hpc/temp/{user}/SPC_genome_preprocessing/{sample}/`

---

## Testing Summary

### Multispecies Test (Yeast + Worm):
- **Yeast (S. cerevisiae)**: ✅ All 120 cells processed, QC compiled, Lander-Waterman plot generated
- **Worm (C. elegans)**: ✅ All 120 cells processed, full QC suite generated including Lander-Waterman plot

### HSC Deep Sequencing Test:
- **HSC4_pilot**: ✅ 500 cells processed, complete QC metrics
- Lander-Waterman plot regenerated with corrected formula
- Theoretical curve properly aligns with data (4 orders of magnitude range)
- Dashboard integration tested and working

### Parallel Merge Test:
- **Dataset**: 501 yeast SAM files (~1.2GB final BAM)
- **Performance**: 26 parallel tasks completed in ~5 minutes
- **Validation**: BAM validated, read counts match expected filtering results

---

## Git Branch Status

**Current Branch**: `main`

**Recent Commits (April 2026)**:
- `fb3a85c`: Refactor barcode annotation into layered architecture
- `79f79f3`: Fix Lander-Waterman plot curve calculation and simplify input handling
- `adcc504`: Add HTML dashboard and preprocessing stats capture for CapWGS pipeline (PI)

**Previous Commits (March 2026, merged from `pipeline-improvements`)**:
- `000f936`: Implement parallel two-stage SAM/BAM merge for 20-40x speedup
- `c03aaa6`: Fix misleading success message when Lander-Waterman plot fails
- `8c33f57`: Add file existence checks to all QC compilation success messages

**Status**: All improvements tested and pushed to main. Dashboard system fully integrated.
