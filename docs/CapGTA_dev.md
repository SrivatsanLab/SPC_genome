# CapGTA Pipeline Development

---

## Known Bugs and Fixes

### Bug: Path duplication in single-cell extraction (March 31, 2026)

**Issue:** When `OUTPUT_NAME` contains path separators (e.g., `worm_CapGTA_L4_pilot/L4_3_UDI7`), the BAM path construction in `CapGTA_PP.sh` incorrectly duplicates the path component.

**Root cause:** Lines 186-187 in `CapGTA_PP.sh`:
```bash
DNA_BAM="${DATA_DIR}/${OUTPUT_NAME}_dna.bam"
RNA_BAM="${DATA_DIR}/${OUTPUT_NAME}_rna.bam"
```

When:
- `OUTPUT_NAME="worm_CapGTA_L4_pilot/L4_3_UDI7"`
- `DATA_DIR="./data/worm_CapGTA_L4_pilot/L4_3_UDI7"`

This expands to:
- `DNA_BAM="./data/worm_CapGTA_L4_pilot/L4_3_UDI7/worm_CapGTA_L4_pilot/L4_3_UDI7_dna.bam"` ❌

But actual file is at:
- `./data/worm_CapGTA_L4_pilot/L4_3_UDI7/L4_3_UDI7_dna.bam` ✓

**Workaround:** Use `basename` to extract just the sample name:
```bash
SAMPLE_NAME=$(basename "${OUTPUT_NAME}")
DNA_BAM="${DATA_DIR}/${SAMPLE_NAME}_dna.bam"
RNA_BAM="${DATA_DIR}/${SAMPLE_NAME}_rna.bam"
```

**Status:** ✅ **FIXED** (commit dd539af on fix-overhang-validation branch, March 31, 2026)

**Affected jobs:** L4 worm pilot experiments failed at single-cell extraction step. Resubmission script created at `bin/resubmit_sc_extraction_L4_pilot.sh` with corrected paths.

---

## Phase 2: Statistical Algorithm for RNA Classification (Feb 25, 2026)

### Algorithm Design

**Goal:** Estimate per-gene RNA fraction and probabilistically reclassify reads from DNA BAMs to RNA BAMs

**Statistical Approach:**

1. **Estimate DNA baseline coverage (λ_DNA)**
   - Compute per-base coverage in large intergenic regions (>10kb from genes)
   - Model as Poisson(λ_DNA) or Negative Binomial if overdispersed
   - Provides expected DNA coverage rate

2. **Test each gene for exonic enrichment**
   - For each gene, aggregate coverage across all exons
   - Expected DNA reads: E[R_DNA] = λ_DNA × total_exonic_length
   - Observed reads: R_obs
   - **Poisson rate test:** P(R ≥ R_obs | λ_DNA × L)
   - Apply Benjamini-Hochberg FDR correction across genes

3. **Estimate per-gene RNA proportion (f_RNA)**
   - Excess reads: R_excess = max(0, R_obs - E[R_DNA])
   - RNA fraction: f_RNA = R_excess / R_obs
   - Apply **Bayesian shrinkage** (Beta prior) for low-coverage genes:
     - Prior: Beta(α=3, β=2) ≈ 60% RNA in expressed genes
     - Posterior: f_RNA = (R_excess + α) / (R_obs + α + β)

4. **Probabilistically reassign reads**
   - For each read overlapping exons:
     - Get gene's f_RNA
     - Assign to RNA with probability f_RNA, else keep in DNA
   - Preserves expected proportions, avoids hard thresholds

5. **Handle edge cases**
   - Reads with splice junctions: keep in RNA (already classified correctly)
   - Intronic reads: keep in DNA (likely true DNA)
   - Intergenic reads: keep in DNA
   - Multi-gene overlaps: use max(f_RNA) (conservative)

### Variant Calling Strategy

**Issue:** Reclassification may move true DNA reads into RNA BAM, reducing variant calling sensitivity

**Solution:** BCFtools supports multiple input BAMs per sample:
```bash
# Call variants using BOTH DNA and RNA BAMs for same cell
bcftools mpileup -f ref.fa cell_dna.bam cell_rna.bam | \
    bcftools call --gvcf 1 -mv -Oz -o cell.g.vcf.gz
```

**Advantages:**
- Recovers coverage from reclassified reads
- No loss of variant calling sensitivity
- RNA BAM provides additional support for true variants
- BCFtools handles duplicate reads automatically

**Implementation:**
- Modify `scripts/CapGTA/sc_variant_calling_bcftools_array.sh` to accept both BAMs
- Update after RNA reclassification is complete

### Implementation Plan

**Scripts to create:**
1. `bin/classify_rna_from_dna.py` - Core algorithm
   - Parse GTF for genes, exons, intergenic regions
   - Estimate λ_DNA from intergenic coverage
   - Compute per-gene f_RNA with statistical tests
   - Probabilistically reassign reads
   - Output new DNA and RNA BAMs

2. `bin/reclassify_rna_array.sh` - SLURM array wrapper
   - Process cells in parallel
   - 1 cell per array task

3. Test scripts for validation on small barcode subset

**Testing strategy:**
- Start with 5-10 test cells
- Validate at each step:
  - λ_DNA estimates reasonable (~0.01-0.1 reads/base)
  - f_RNA distributions make sense (0-0.9 range)
  - RNA counts improve significantly
  - DNA intergenic coverage unchanged

**Parameters to tune:**
```python
min_intergenic_size = 10000  # bp
fdr_threshold = 0.05
alpha_prior = 3  # RNA reads prior
beta_prior = 2   # DNA reads prior
min_gene_coverage = 10  # minimum reads for reliable estimate
```

### Validation Metrics

After reclassification, check:
1. **RNA counts per cell:** Should increase from ~2,141 to ~5,000-10,000
2. **DNA intergenic coverage:** Should remain unchanged (sanity check)
3. **Gene expression correlation:** With WormSeq reference atlas
4. **Variant calling:** No quality degradation with dual-BAM approach

### Note on Cell Count

When regenerating count matrices, use **top 1200 barcodes** instead of 1000 (previous selection was slightly conservative based on knee plot analysis). 
---

## Phase 1: Pipeline Improvements - Parallel Merge + Comprehensive QC (March 2026)

### Overview
Bring CapGTA pipeline up to parity with recent CapWGS improvements:
- **Parallel two-stage SAM/BAM merge** (20-40x speedup for bulk BAM creation)
- **Comprehensive QC outputs** (Picard metrics, Lorenz curves, gini coefficients, Lander-Waterman plots)
- **Exonic enrichment metrics** (DNA vs RNA, pre/post RNA rescue)
- **Unified single-cell processing array** (avoid race conditions)
- **RNA rescue integration** (optional, algorithm in Phase 2)

### Motivation
Current CapGTA pipeline lacks:
1. QC metrics for assessing DNA coverage uniformity and quality
2. Quantitative metrics for RNA/DNA separation quality (exonic enrichment)
3. Fast parallel merge (currently serial, slow for large datasets)
4. Structured integration point for RNA rescue algorithm

### Implementation Roadmap

**Full detailed roadmap**: See `CAPGTA_IMPROVEMENTS_ROADMAP.md`

**Summary:**

**Phase 1a: Refactor Core Pipeline** (Priority 1)
1. Move shared scripts to `scripts/utils/` (parallel_merge_array.sh, compile_lorenz.py, compile_qc_metrics.py, lorenz.py)
2. Rewrite `CapGTA_PP.sh` to match CapWGS structure:
   - Inline two-stage parallel merge (DNA and RNA separately)
   - Inline cell detection and barcode statistics
   - Inline QC compilation
3. Create unified single-cell processing array (`scripts/CapGTA/sc_extract_preprocess_qc_array.sh`):
   - Extract DNA and RNA BAMs
   - MarkDuplicates on DNA (QC only)
   - Generate DNA bigwig, Lorenz curve, gini coefficient
   - Collect Picard alignment metrics (DNA)
   - Calculate exonic enrichment (DNA and RNA, pre-rescue)
4. Add QC compilation section (Picard metrics, Lorenz curves, exonic enrichment, Lander-Waterman plot)
5. Delete obsolete scripts (concatenate_gta.sh, sc_from_bam_gta.sh)

**Phase 1b: RNA Rescue Integration** (Priority 2)
1. Add `--enable-rna-rescue` flag to CapGTA_PP.sh
2. Create RNA rescue array stub (`scripts/CapGTA/rna_rescue_array.sh`)
3. Create post-rescue QC array (`scripts/CapGTA/post_rescue_qc_array.sh`)
4. Add conditional execution logic (skip rescue steps if flag not set)

**Phase 2: RNA Rescue Algorithm Implementation** (Future Work)
- Implement statistical algorithm in rna_rescue_array.sh
- Validate on test cells
- Full production runs

### Exonic Enrichment Metric

**Definition:**
```
Exonic Enrichment = (Exonic Reads / Total Reads) - (Exonic Bases / Total Reference Bases)
```

**Expected Values:**
- **DNA (pre-rescue)**: ~0 (no enrichment, random genomic sampling)
- **RNA (pre-rescue)**: >0.2 (20%+ enrichment for expressed genes)
- **DNA (post-rescue)**: ~0 (should remain unchanged)
- **RNA (post-rescue)**: >0.4 (increased enrichment after rescue)

**Purpose:**
- Quantify RNA/DNA separation quality
- Validate RNA rescue algorithm effectiveness
- QC metric for library preparation quality

### Directory Structure Changes

**Matches CapWGS structure:**
```
data/{sample}/
├── {sample}_dna.bam              # Bulk DNA BAM
├── {sample}_rna.bam              # Bulk RNA BAM
├── sc_outputs/                   # Single-cell BAMs and VCFs
│   ├── {barcode}_dna.bam         # Overwritten by rescue if enabled
│   ├── {barcode}_rna.bam         # Overwritten by rescue if enabled
│   ├── {barcode}_dna.bw
│   └── {barcode}.vcf.gz
└── qc_metrics/                   # Per-cell QC outputs
    ├── {barcode}_alignment_metrics.txt
    ├── {barcode}_duplication_metrics.txt
    ├── {barcode}_lorenz.csv
    ├── {barcode}_gini.txt
    ├── {barcode}_exonic_enrichment_dna.txt
    ├── {barcode}_exonic_enrichment_rna.txt
    ├── {barcode}_exonic_enrichment_dna_postrescue.txt
    └── {barcode}_exonic_enrichment_rna_postrescue.txt

results/{sample}/
├── compiled_qc_metrics.csv           # Picard + gini
├── compiled_lorenz_curves.csv
├── compiled_exonic_enrichment.csv    # Pre/post rescue
├── lander_waterman_coverage.png
├── kneeplot.png
├── readcounts.csv
├── real_cells.txt
├── barcode_assignment_stats.txt
├── sc_variants_merged.vcf.gz
└── rna_counts_matrix.csv
```

### New Scripts Required

**scripts/CapGTA/:**
1. `sc_extract_preprocess_qc_array.sh` - Unified SC processing
2. `rna_rescue_array.sh` - RNA reclassification (stub initially)
3. `post_rescue_qc_array.sh` - Re-run QC after rescue
4. `calculate_exonic_enrichment.py` - Compute exonic enrichment

**scripts/utils/:**
1. `combine_readcounts_gta.py` - Combine DNA + RNA counts
2. `compile_exonic_enrichment.py` - Compile enrichment metrics

### Expected Improvements

**Performance:**
- Bulk BAM creation: 3-5 days → 3-6 hours (20-40x speedup)
- No race conditions from sequential job dependencies

**QC Coverage:**
- Per-cell Picard alignment metrics (mapping quality, insert size, etc.)
- Coverage uniformity (Lorenz curves, gini coefficients)
- Lander-Waterman coverage plots
- Exonic enrichment quantification

**Pipeline Structure:**
- Matches CapWGS for maintainability
- Clear integration point for RNA rescue algorithm
- Consolidated QC outputs

### Testing Plan

**Phase 1a Testing:**
- Run on existing worm CapGTA data (worm_CapGTA_UDI_5)
- Verify parallel merge speedup
- Validate all QC outputs generated
- Confirm exonic enrichment: DNA ~0, RNA >0.2

**Phase 1b Testing:**
- Run with `--enable-rna-rescue` flag (stub)
- Verify QC outputs correctly overwritten
- Confirm pipeline completes successfully

**Phase 2 Testing:**
- Validate algorithm on 5-10 test cells
- Check RNA count improvement (2,141 → 5,000-10,000)
- Verify DNA intergenic coverage unchanged
- Confirm variant calling quality maintained

### Related Issues

See GitHub issue #11 for detailed implementation tracking.
