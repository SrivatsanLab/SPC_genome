# CapGTA Pipeline Enhancement Plan
## Recent CapWGS Enhancements - Compatibility Analysis

**Date**: April 9, 2026
**Purpose**: Assess recent CapWGS improvements for integration into CapGTA pipeline overhaul

---

## Summary

Recent CapWGS improvements (March-April 2026) include 10 major enhancements. Of these, **7 are directly applicable** to CapGTA and should be integrated into the pipeline overhaul. This document provides a comprehensive compatibility analysis and updated implementation roadmap.

---

## Enhancement Compatibility Matrix

| # | Enhancement | CapGTA Compatible? | Priority | Notes |
|---|-------------|-------------------|----------|-------|
| 1 | YAML config file support | ✅ Already Implemented | N/A | CapGTA_PP.sh lines 8-47 |
| 2 | `wait_for_job()` error handling | ✅ **HIGH PRIORITY** | P1 | Essential for robust pipeline |
| 3 | `--use-existing-chunks` flag | ✅ **HIGH PRIORITY** | P1 | Critical for development/testing |
| 4 | Sample name extraction | ✅ Already Implemented | N/A | CapGTA_PP.sh line 187 |
| 5 | Preprocessing stats compilation | ✅ **MEDIUM PRIORITY** | P2 | Valuable QC, but not critical |
| 6 | Unfiltered BAM + flagstat | ✅ **HIGH PRIORITY** | P1 | DNA & RNA separately |
| 7 | Barcode assignment statistics | ✅ **HIGH PRIORITY** | P1 | DNA + RNA combined stats |
| 8 | HTML Dashboard generation | ✅ **VERY HIGH PRIORITY** | P1 | Must-have for modern QC |
| 9 | Variant caller option (-v) | ❌ Not Applicable | N/A | CapGTA uses bcftools only |
| 10 | Known sites directory (-k) | ❌ Not Applicable | N/A | CapGTA has no GATK/BQSR |

---

## Detailed Enhancement Analysis

### 1. ✅ YAML Config File Support
**Status**: Already implemented in CapGTA
**Action**: None required

CapGTA_PP.sh already has identical YAML parsing (lines 8-47).

---

### 2. ✅ `wait_for_job()` Function with Error Handling
**Status**: **NOT in CapGTA - MUST ADD**
**Priority**: P1 (HIGH)
**Compatibility**: 100% compatible

**What it does**:
```bash
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
```

**Why it's critical**:
- Current CapGTA uses `--dependency=afterok:$JOB_ID` for job dependencies
- This can cause silent failures to propagate through pipeline
- `wait_for_job()` allows main script to detect failures and exit early

**CapGTA-specific considerations**:
- Must wait for: alignment array, concatenation, single-cell arrays, variant calling
- DNA/RNA parallel merge arrays can be launched simultaneously with `wait` for both

**Implementation**:
- Add function to top of CapGTA_PP.sh (after helper functions)
- Replace all `sbatch --dependency=afterok:$JOB_ID` with inline job submission + `wait_for_job`

---

### 3. ✅ `--use-existing-chunks` Flag
**Status**: **NOT in CapGTA - MUST ADD**
**Priority**: P1 (HIGH)
**Compatibility**: 100% compatible

**What it does**:
- Skips FASTQ splitting step (expensive for large files)
- Reuses existing chunks in temp directory
- Requires `chunk_indices.txt` to exist in results directory

**Why it's critical**:
- Essential for development and testing (avoid re-chunking)
- Useful for pipeline restarts after failures
- Already mentioned in CAPGTA_IMPROVEMENTS_ROADMAP.md line 19

**Implementation**:
```bash
# Add flag parsing
USE_EXISTING_CHUNKS=false
for arg in "$@"; do
    if [ "$arg" = "--use-existing-chunks" ]; then
        USE_EXISTING_CHUNKS=true
    fi
done

# Add conditional chunking logic (copy from CapWGS_PP.sh lines 197-243)
if [ "$USE_EXISTING_CHUNKS" = true ]; then
    # Verify chunk_indices.txt exists
    # Load chunk count
else
    # Execute normal FASTQ chunking
fi
```

---

### 4. ✅ Sample Name Extraction
**Status**: Already implemented in CapGTA
**Action**: None required

CapGTA_PP.sh line 187 already uses `basename "${OUTPUT_NAME}"` to handle paths.

---

### 5. ✅ Preprocessing Statistics Compilation
**Status**: **NOT in CapGTA - RECOMMENDED**
**Priority**: P2 (MEDIUM)
**Compatibility**: 100% compatible

**What it does**:
- Aggregates per-chunk demux, trimming, and alignment stats
- Produces:
  - `preprocessing_summary.json` (overall stats)
  - `adapter_histogram_r1.csv` (adapter position distribution)
  - `adapter_histogram_r2.csv`

**Why it's valuable**:
- Quantifies barcode mapping rate across chunks
- Tracks adapter contamination patterns
- Essential input for HTML dashboard
- Debugging tool for preprocessing issues

**CapGTA-specific considerations**:
- STAR alignment array (`PP_array_gta_star_only.sh`) would need to output per-chunk JSON
- Requires modifications to array script to capture stats
- Dashboard integration requires this feature

**Implementation**:
1. Modify `scripts/CapGTA/PP_array_gta_star_only.sh`:
   - Add JSON stats output (copy from CapWGS version)
   - Create `${STATS_DIR}/chunk_${chunk}.json` per task
2. Add compilation step to CapGTA_PP.sh (after alignment array completes)
3. Use `scripts/CapWGS/compile_preprocessing_stats.py` (works for both pipelines)

**Decision**: Include in Phase 1 if dashboard is prioritized, otherwise defer to Phase 2

---

### 6. ✅ Unfiltered BAM + Flagstat
**Status**: **NOT in CapGTA - MUST ADD**
**Priority**: P1 (HIGH)
**Compatibility**: 100% compatible (DNA & RNA separately)

**What it does**:
- Creates unfiltered BAM before applying `-f 0x2` filtering
- Runs `samtools flagstat` on unfiltered BAM
- Provides true mapping rate before filtering artifacts
- Deletes unfiltered BAM after flagstat

**Why it's critical**:
- `-f 0x2` filtering removes ~4% of reads (chunking artifacts)
- Flagstat on filtered BAM shows misleading mapping rates
- Unfiltered flagstat = true alignment quality

**CapGTA-specific considerations**:
- Must create TWO unfiltered BAMs (DNA and RNA)
- Must run flagstat on both
- Output: `${RESULTS_DIR}/flagstat_dna.txt` and `${RESULTS_DIR}/flagstat_rna.txt`

**Implementation** (copy from CapWGS_PP.sh lines 336-349):
```bash
# For DNA BAM:
UNFILTERED_DNA_BAM="${DATA_DIR}/${SAMPLE_NAME}_dna.unfiltered.bam"
samtools merge -@ 4 -c -b "${RESULTS_DIR}/intermediate_bam_list_dna.txt" -O BAM - | \
    samtools sort -@ 4 -o "${UNFILTERED_DNA_BAM}"
samtools flagstat -@ 4 "${UNFILTERED_DNA_BAM}" > "${RESULTS_DIR}/flagstat_dna.txt"
samtools view -@ 4 -f 0x2 -b "${UNFILTERED_DNA_BAM}" | \
    samtools sort -@ 4 -o "${DATA_DIR}/${SAMPLE_NAME}_dna.bam"
rm "${UNFILTERED_DNA_BAM}"

# Repeat for RNA BAM
```

---

### 7. ✅ Barcode Assignment Statistics
**Status**: **NOT in CapGTA - MUST ADD**
**Priority**: P1 (HIGH)
**Compatibility**: 100% compatible (with DNA+RNA modification)

**What it does**:
- Calculates percentage of input reads assigned to valid barcodes
- Reports: total input reads, total assigned reads, assignment rate
- Output: `${RESULTS_DIR}/barcode_assignment_stats.txt`

**Why it's critical**:
- Key QC metric for library preparation quality
- Identifies barcode design or sequencing issues
- Dashboard summary metric

**CapGTA-specific considerations**:
- Must calculate for **combined DNA + RNA reads**
- Input reads = `READ_COUNT * 2` (paired-end)
- Assigned reads = DNA reads + RNA reads

**Implementation** (adapted from CapWGS_PP.sh lines 370-400):
```bash
# After DNA and RNA read counting
total_dna_reads=$(tail -n +2 "${RESULTS_DIR}/readcounts_dna.csv" | cut -d',' -f2 | awk '{s+=$1} END {print s}')
total_rna_reads=$(tail -n +2 "${RESULTS_DIR}/readcounts_rna.csv" | cut -d',' -f2 | awk '{s+=$1} END {print s}')
total_reads_in_bams=$((total_dna_reads + total_rna_reads))

ASSIGNMENT_STATS="${RESULTS_DIR}/barcode_assignment_stats.txt"
{
    echo "=========================================="
    echo "Barcode Assignment Statistics"
    echo "=========================================="
    echo ""

    total_input_reads=$((READ_COUNT * 2))
    assignment_rate=$(awk "BEGIN {printf \"%.2f\", ($total_reads_in_bams / $total_input_reads) * 100}")
    echo "Total input reads (paired-end): ${total_input_reads}"
    echo "Reads assigned to valid barcodes (DNA + RNA): ${total_reads_in_bams}"
    echo "  - DNA reads: ${total_dna_reads}"
    echo "  - RNA reads: ${total_rna_reads}"
    echo "Assignment rate: ${assignment_rate}%"
    echo ""
    echo "=========================================="
} > "${ASSIGNMENT_STATS}"
```

---

### 8. ✅ HTML Dashboard Generation
**Status**: **NOT in CapGTA - VERY HIGH PRIORITY**
**Priority**: P1 (VERY HIGH)
**Compatibility**: 95% compatible (needs CapGTA-specific adaptations)

**What it does**:
- Single-file self-contained HTML dashboard (like 10x Cell Ranger)
- Interactive Plotly.js visualizations (zoom, pan, hover)
- Tabbed interface: Summary, Per-Cell QC, Coverage
- Embeds all data, plots, and code in one shareable file

**Dashboard Features**:
- **Summary Tab**:
  - Metric cards (cells, reads, mapping rates, duplication)
  - Read retention Sankey diagram (demux → trimming → alignment)
  - Preprocessing stats table (with flagstat integration)
  - Adapter position histograms (R1 + R2)
  - Interactive knee plot
- **Per-Cell QC Tab**:
  - Alignment rate, duplication rate, coverage distributions
  - Gini coefficient histogram
  - Coverage vs genome fraction scatter
- **Coverage Tab**:
  - Lorenz curves (sampled cells)
  - Lander-Waterman plot

**Why it's critical**:
- Modern standard for pipeline QC (10x, Parse, etc. all use HTML dashboards)
- Enables easy sharing with collaborators
- Interactive exploration of QC metrics
- Professional presentation for papers/presentations

**CapGTA-specific adaptations required**:
1. **DNA/RNA separation metrics**:
   - Show DNA and RNA read counts separately in summary
   - Separate alignment rate metrics for DNA and RNA
   - Display exonic enrichment metrics (pre/post rescue)
2. **RNA-specific plots**:
   - RNA count histogram (genes detected per cell)
   - Optionally: RNA PCA plot (if count matrix available)
3. **CapGTA-specific summary cards**:
   - DNA reads per cell
   - RNA reads per cell (pre/post rescue if applicable)
   - RNA counts per cell (median genes detected)
   - Exonic enrichment (DNA should be ~0, RNA >0.2)
4. **Dual Lander-Waterman plots**:
   - One for DNA coverage
   - One for RNA coverage (may not be applicable)

**Implementation path**:
1. **Phase 1**: Use existing `generate_dashboard.py` with DNA BAM only
   - Dashboard shows DNA QC metrics only (like CapWGS)
   - RNA metrics displayed as summary cards but not plotted
2. **Phase 2**: Create `scripts/CapGTA/generate_dashboard_gta.py`
   - Extend with RNA-specific visualizations
   - Add DNA/RNA comparison plots
   - Integrate exonic enrichment metrics

**Files to create/modify**:
- Option A: Extend `scripts/CapWGS/generate_dashboard.py` with `--mode gta` flag
- Option B: Create separate `scripts/CapGTA/generate_dashboard_gta.py` (cleaner)

**Recommendation**: Start with Option A (extend existing script) for faster deployment, refactor to Option B if divergence becomes significant.

---

### 9. ❌ Variant Caller Option (-v flag)
**Status**: Not applicable to CapGTA
**Reason**: CapGTA uses bcftools only, no GATK variant calling option

CapWGS supports `bcftools`, `gatk`, or `none`. CapGTA pipeline only uses bcftools for variant calling and does not require this flexibility.

**Decision**: Do not implement

---

### 10. ❌ Known Sites Directory (-k flag)
**Status**: Not applicable to CapGTA
**Reason**: CapGTA has no BQSR step (bcftools variant calling only)

BQSR (Base Quality Score Recalibration) is a GATK best practice requiring known variant sites. CapGTA uses bcftools which does not require BQSR.

**Decision**: Do not implement

---

## Additional CapWGS Improvements to Consider

### A. Reference Genome Path Flexibility
**Status**: Should verify CapGTA implementation
**Priority**: P2 (MEDIUM)

CapWGS can auto-detect FASTA files with non-standard names (e.g., `GCA_028201515.1_genomic.fna`).

**Action**: Verify CapGTA STAR alignment array handles this correctly.

---

### B. Lander-Waterman Plot Formula Fix
**Status**: Should be in shared utilities
**Priority**: P1 (HIGH)

CapWGS fixed the theoretical curve formula to include read length factor (commit `79f79f3`).

**Current formula** (CORRECT):
```python
mean_read_length = df['mean_read_length'].mean()  # ~142 bp
y_theory = 1 - np.exp(-reads_per_mb * mean_read_length / 1e6)
```

**Action**: Verify `scripts/utils/plot_lander_waterman.py` is up-to-date and will be used by CapGTA.

---

### C. Barcode Annotation Refactoring
**Status**: Shared utilities
**Priority**: P3 (LOW)

CapWGS refactored barcode annotation into layered architecture (commit `fb3a85c`):
- `scripts/utils/barcode_core.py`: Shared utilities
- `scripts/utils/barcode_annotation.py`: DataFrame API (notebooks)
- `scripts/CapWGS/assign_conditions.py`: CLI tool (dashboards)

**Action**: If CapGTA dashboard supports experimental conditions, use these shared utilities.

---

## Updated CapGTA Roadmap

### Phase 1: Core Pipeline Refactor + Essential Enhancements

**Goal**: Bring CapGTA to parity with CapWGS structure and essential improvements

**Tasks**:
1. **Add `wait_for_job()` function** (P1 - CRITICAL)
   - Replace all `--dependency=afterok:` with inline waiting
   - Enable early failure detection
2. **Add `--use-existing-chunks` flag** (P1 - CRITICAL)
   - Essential for development/testing
   - Copy logic from CapWGS_PP.sh lines 197-243
3. **Implement two-stage parallel merge** (P1 - already in roadmap)
   - DNA and RNA BAMs separately (can run in parallel)
   - Copy from CapWGS implementation
4. **Add unfiltered BAM + flagstat** (P1 - NEW)
   - Create `flagstat_dna.txt` and `flagstat_rna.txt`
   - Accurate mapping rate metrics
5. **Calculate barcode assignment statistics** (P1 - NEW)
   - DNA + RNA combined statistics
   - Output to `barcode_assignment_stats.txt`
6. **Move shared scripts to `scripts/utils/`** (P1 - already in roadmap)
   - `parallel_merge_array.sh`
   - `compile_lorenz.py`
   - `compile_qc_metrics.py`
   - `lorenz.py`
7. **Create unified single-cell processing array** (P1 - already in roadmap)
   - Extract DNA + RNA BAMs
   - MarkDuplicates (QC only)
   - Generate bigwig, Lorenz curve, gini coefficient
   - Collect Picard metrics
   - Calculate exonic enrichment (pre-rescue)
8. **Rewrite CapGTA_PP.sh main script** (P1 - already in roadmap)
   - Match CapWGS structure
   - Inline cell detection
   - Inline QC compilation
   - Add all new enhancements above

**Estimated Time**: 2-3 days

---

### Phase 2: Advanced QC and Dashboard

**Goal**: Add comprehensive QC outputs and modern dashboard

**Tasks**:
1. **Preprocessing statistics compilation** (P2)
   - Modify alignment array to output per-chunk JSON
   - Add compilation step to main script
   - Output `preprocessing_summary.json`, adapter histograms
2. **Generate HTML Dashboard** (P1 - VERY HIGH PRIORITY)
   - Start with existing `generate_dashboard.py` for DNA QC
   - Add CapGTA-specific summary metrics (DNA/RNA reads, exonic enrichment)
   - Display RNA count summary statistics
   - Add dashboard generation to main pipeline
3. **Create exonic enrichment compilation** (P2 - already in roadmap)
   - `scripts/utils/compile_exonic_enrichment.py`
   - Input to dashboard
4. **Add RNA rescue integration stubs** (P2 - already in roadmap)
   - `--enable-rna-rescue` flag
   - RNA rescue array stub
   - Post-rescue QC array
   - Conditional execution logic

**Estimated Time**: 2-3 days

---

### Phase 3: RNA Rescue Algorithm (Future Work)

**Goal**: Implement statistical RNA reclassification (unchanged from original roadmap)

**Tasks**:
1. Implement algorithm in `rna_rescue_array.sh`
2. Validate on test cells
3. Full production run
4. Enhance dashboard with post-rescue metrics

**Estimated Time**: 3-5 days

---

## Summary of New Additions to Roadmap

**High Priority (Phase 1)**:
- ✅ `wait_for_job()` function with error handling
- ✅ `--use-existing-chunks` flag
- ✅ Unfiltered BAM + flagstat (DNA and RNA)
- ✅ Barcode assignment statistics (DNA + RNA combined)

**Medium/High Priority (Phase 2)**:
- ✅ Preprocessing statistics compilation
- ✅ HTML Dashboard generation (adapted for CapGTA)

**Total New Items**: 6 enhancements (4 critical, 2 important)

---

## Implementation Priority Order

### Absolutely Critical (Do First):
1. `wait_for_job()` function - prevents silent failure propagation
2. `--use-existing-chunks` flag - essential for development workflow
3. Two-stage parallel merge - already in roadmap, core performance improvement

### Very Important (Do Second):
4. Barcode assignment statistics - key QC metric
5. Unfiltered BAM + flagstat - accurate mapping rate reporting
6. HTML Dashboard - modern QC standard, high impact

### Nice to Have (Do Third):
7. Preprocessing statistics compilation - valuable for dashboard, but not critical

---

## Testing Plan

**Phase 1 Testing**:
- Run on existing worm CapGTA data (e.g., `worm_CapGTA_UDI_5`)
- Verify:
  - `wait_for_job()` detects failures correctly
  - `--use-existing-chunks` works for retry scenarios
  - Parallel merge achieves expected speedup (20-40x)
  - Flagstat shows accurate mapping rates
  - Barcode assignment stats are reasonable (typically 60-80%)

**Phase 2 Testing**:
- Verify dashboard displays correctly
- Validate all interactive features work
- Confirm CapGTA-specific metrics display properly
- Test with/without experimental conditions

---

## Files to Create

**New Scripts (Phase 1)**:
1. *(None - all changes are modifications to existing scripts)*

**New Scripts (Phase 2)**:
1. `scripts/utils/compile_exonic_enrichment.py` (already in roadmap)
2. `scripts/CapGTA/generate_dashboard_gta.py` (or extend existing CapWGS version)

**Modified Scripts (Phase 1)**:
1. `CapGTA_PP.sh` - major rewrite incorporating all enhancements
2. `scripts/CapGTA/PP_array_gta_star_only.sh` - add stats JSON output (Phase 2)
3. `scripts/CapGTA/sc_extract_preprocess_qc_array.sh` - create new (already in roadmap)

**Scripts to Move (Phase 1)**:
1. `scripts/CapWGS/parallel_merge_array.sh` → `scripts/utils/`
2. `scripts/CapWGS/compile_lorenz.py` → `scripts/utils/`
3. `scripts/CapWGS/compile_qc_metrics.py` → `scripts/utils/`
4. `scripts/CapWGS/lorenz.py` → `scripts/utils/`

---

## Conclusion

Recent CapWGS improvements (March-April 2026) provide **7 valuable enhancements** applicable to CapGTA:

**Critical (P1)**:
- `wait_for_job()` function (reliability)
- `--use-existing-chunks` flag (development efficiency)
- Unfiltered BAM + flagstat (QC accuracy)
- Barcode assignment statistics (library QC)
- HTML Dashboard (modern QC standard)

**Important (P2)**:
- Preprocessing statistics compilation (detailed QC)

**Already Planned**:
- Parallel two-stage merge (performance)

These enhancements should be integrated into the CapGTA pipeline overhaul alongside the originally planned improvements (parallel merge, comprehensive QC, RNA rescue integration). The combined effort will bring CapGTA to full parity with CapWGS while adding CapGTA-specific features (exonic enrichment, RNA rescue).

**Recommendation**: Proceed with Phase 1 implementation incorporating all P1 enhancements. The incremental effort is minimal (1-2 additional days) compared to the substantial QC and reliability benefits.
