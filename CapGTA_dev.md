# CapGTA Pipeline Development

Raw Data: `./SRA_upload/20260206_LH00740_0181_B232F3VLT4/Unaligned/L4/`

# Initial Issues (Feb 11, 2026)

I am processing new data using the @CapGTA_PP.sh pipeline.
I submitted this new data with @bin/submit_worm_big_run.sh
Read that script to see where the data is, how many reads, etc. There are 4 separate groups of fastqs that are being processed together (UDI_5,UDI_6,UDI_7,UDI_8)
So, there were 4 pipeline jobs, and you can read there output logs:
UDI_5 : @SLURM_outs/GTA_Preprocessing_46347497.out
UDI_6 : @SLURM_outs/GTA_Preprocessing_46347498.out
UDI_7 : @SLURM_outs/GTA_Preprocessing_46347499.out
UDI_8 : @SLURM_outs/GTA_Preprocessing_46347500.out
2 of the jobs, corresponding to UDI_5 and UDI_8, successfully made it through the concatenation step, but were failing on single cell extraction, so I cancelled their downstream jobs.
The jobs corresponding to UDI_6 and UDI_7 failed before concatenation, leaving the concatenation jobs (46366324, 46358195) hung on `DependencyNeverSatisfied`
You may need to read all of the CapGTA scripts in `scripts/CapGTA/` to understand the pipeline.
I need to:
1. Understand why the UDI_6 and UDI_7 jobs failed before concatenation.
    - Given that 2 jobs made it through concatenation, and I have used versions of these scripts in the past, I think it was probably not related to my code, and instead due to something simple like using too much memory, an incidental node failure, or a timeout
    - If it was actually a code issue, propose a fix.
2. Understand why the jobs were failing at the single cell extraction step.
    - Check the output logs of the array jobs
3. Fix a known bug in single cell extraction:
    - Lines 33 and 40 do a `grep` operation on the bam files to extract single cells
    `samtools view -h "${dna_sam}" | grep -E "^@|CB:Z:${barcode}" | samtools view -b -o "${dna_bam}"`
    - This should be changed to
    `samtools view -h ${input_file} | grep -E "^@|CB:Z:${barcode}" | samtools sort -o "${output_bam}"`
    - piping to sort ensures the output bams are sorted, preventing downstream issues.
    - Also add a line to index the bams:
    `samtools index "${output_bam}"`
4. Once we have resolved the above issues, I would like to cleanup the existing temp and data directories, and resubmit the jobs.
    - edit the submissions script to use -c 1000 to override the cell prediction and use the top 1000 barcodes. based on the data I saw from UDI_5, the predictions seem off when I run them on this coassay data
    - also remove the argument for number of chunks and just use the default 500.

---

## Session Summary (Feb 12, 2026)

### Issues Diagnosed

**1. UDI_6 & UDI_7 Preprocessing Failures**
- Investigation revealed only 1 task failed per job (46358194_1024, 46366323_1629) out of 5000 total tasks
- Exit code 127 indicates sporadic node failures, not code issues
- 4998/5000 tasks completed successfully for each run
- Conclusion: Transient infrastructure issues, not bugs in code

**2. Single Cell Extraction Bug**
- Root cause: `grep` filtering disrupts coordinate sort order in BAM files
- Affects both DNA and RNA single-cell BAMs in `scripts/CapGTA/extract_sc_array_gta.sh`
- Would cause downstream failures in indexing and variant calling

### Fixes Implemented

**1. Fixed Single Cell Extraction (`scripts/CapGTA/extract_sc_array_gta.sh`)**
- Lines 33 & 40: Changed `samtools view -b -o` to `samtools sort -o`
- Added lines 34 & 41: Added `samtools index` for both DNA and RNA BAMs
- Ensures all extracted single-cell BAMs are coordinate-sorted and indexed

**2. Implemented gVCF Mode for Variant Calling**
- **Issue discovered**: BCFtools variant calling only reported variant sites, not reference calls
- Created [GitHub Issue #3](https://github.com/SrivatsanLab/SPC_genome/issues/3) to track this bug across all pipelines

**Modified `scripts/CapGTA/sc_variant_calling_bcftools_array.sh`:**
- Added `--gvcf 1` flag to `bcftools call` to report all sites with DP≥1
- Changed output filename from `${barcode}_variants.vcf.gz` to `${barcode}.g.vcf.gz`
- Updated intermediate filename pattern from `_raw.vcf.gz` to `_raw.g.vcf.gz`

**Modified `scripts/CapGTA/merge_sc_vcfs.sh`:**
- Updated glob pattern from `*_variants.vcf.gz` to `*.g.vcf.gz`
- Updated barcode extraction from `${filename%_variants.vcf.gz}` to `${filename%.g.vcf.gz}`

**3. Updated Submission Script (`bin/submit_worm_big_run.sh`)**
- Removed `-n 5000` argument (now uses default 500 chunks)
- Added `-c 1000` argument to override cell detection and select top 1000 barcodes by read count
- Rationale: Automatic cell detection appears biased on co-assay data based on UDI_5 results

### Jobs Resubmitted

Cleaned up ~9TB of temp data and all partial outputs, then resubmitted all 4 jobs:
- UDI_5: Job 46407150
- UDI_6: Job 46407151
- UDI_7: Job 46407152
- UDI_8: Job 46407153

All jobs submitted with:
- Corrected extraction script (sorted + indexed BAMs)
- gVCF mode enabled for variant calling
- Top 1000 cells by read count
- Default 500 chunks for parallelization

### Files Modified
- `scripts/CapGTA/extract_sc_array_gta.sh`
- `scripts/CapGTA/sc_variant_calling_bcftools_array.sh`
- `scripts/CapGTA/merge_sc_vcfs.sh`
- `bin/submit_worm_big_run.sh` (not tracked in git)

### Related Issues
- Issue #3: BCFtools variant calling should use gVCF mode across all pipelines (CapGTA, CapWGS, Bulk)

---

## Session Summary (Feb 17, 2026)

### Problem Diagnosed

**Job Explosion Causing SLURM Failures**
- Resubmitted jobs (46407150-153) completed with exit code 0, but `sc_outputs/` directories were empty
- Investigation revealed single cells were extracted to temp but never merged or copied to final location
- Root cause: `sc_from_chunks_gta.sh` hit SLURM's QOSMaxSubmitJobPerUserLimit
  - 500 chunks × 1000 cells = 500,000 array tasks submitted
  - SLURM limit: ~50,000 jobs per user
  - Result: Some job submissions failed, creating malformed dependency strings
  - All downstream jobs (merge, copy, variant calling, RNA counting) failed with "Job dependency problem"

### Solution Implemented

**Consolidated Single-Cell Extraction (Addresses Issue #4)**

Created unified extraction utilities in `scripts/utils/`:
- `extract_sc_from_bam_array.sh`: Core array job for extracting cells from bulk BAMs
  - Accepts optional suffix parameter for dual-output pipelines (e.g., "_dna", "_rna")
  - Uses `samtools sort` to ensure coordinate-sorted output
  - 24-hour timeout for large NovaSeq runs
- `sc_from_bam.sh`: Generic wrapper for submitting extraction jobs

**CapGTA-Specific Changes:**
- Created `scripts/CapGTA/sc_from_bam_gta.sh`: Orchestrates dual extraction (DNA + RNA) from bulk BAMs
  - Replaces chunk-based approach with bulk BAM extraction
  - Submits downstream variant calling and RNA counting jobs
- Updated `CapGTA_PP.sh` to use new bulk BAM extraction approach
  - Extracts from `data/{sample}/{sample}_dna.bam` and `data/{sample}/{sample}_rna.bam`
  - No longer processes per-chunk SAM files

**Impact:**
- Task count reduced: 500,000 → 2,000 array tasks (1000 DNA + 1000 RNA)
- Eliminates SLURM job submission failures
- Also updated CapWGS_PP.sh and CapWGS_PP_QC_only.sh to use unified utilities

### Jobs Resubmitted

Cleaned up temp `sc_outputs/` directories and manually resubmitted extraction for all 4 samples:
- UDI_5: Job 47693348 → spawned jobs 47693376-47693380
- UDI_6: Job 47694133
- UDI_7: Job 47694134
- UDI_8: Job 47694140

Each wrapper job submits:
1. DNA extraction array (1000 tasks)
2. RNA extraction array (1000 tasks)
3. Variant calling array (1000 tasks, depends on DNA extraction)
4. VCF merge job (depends on variant calling)
5. RNA count matrix job (depends on RNA extraction)

### Files Modified
- Created: `scripts/utils/extract_sc_from_bam_array.sh`
- Created: `scripts/utils/sc_from_bam.sh`
- Created: `scripts/CapGTA/sc_from_bam_gta.sh`
- Modified: `CapGTA_PP.sh`
- Modified: `CapWGS_PP.sh`
- Modified: `CapWGS_PP_QC_only.sh`

### Commits
- Merged `2-vcaller_selection` branch into main
- Merged `CapGTA_dev` branch into main
- Commit fd4da94: "Consolidate single-cell extraction into unified utilities (fixes #4)"

---

## Session Summary (Feb 18, 2026)

### Problem Diagnosed

**Variant Calling Complete Failure (100% Failure Rate)**
- DNA/RNA extraction completed successfully for all 4 samples (UDI_5-8), 1000 cells each
- All 4000 variant calling jobs failed (exit codes 11 and 127)
- Merge jobs stuck in PENDING with `DependencyNeverSatisfied`
- All gVCF output files were 0 bytes

**Root Causes Identified:**

1. **Reference Genome Mismatch**
   - BAMs aligned to N2 strain-specific reference with renamed chromosomes (I, II, III, IV, V, X)
   - Variant calling script used wrong Ensembl WBcel235 reference
   - **Correct reference**: `data/reference/worm_GCA_028201515.1_STAR/GCA_028201515.1_genomic_renamed.fna`

2. **bcftools Bug: `-A` Flag Incompatible with gVCF Mode**
   - Combining `bcftools call -A --gvcf 1` triggers segmentation faults
   - Two types of failures observed:
     - Exit code 11: Segmentation fault (484/1000 jobs)
     - Exit code 127: "Unexpected number of PL fields" (516/1000 jobs)
   - Removing `-A` flag completely resolves crashes

### Testing & Validation

**Manual Testing:**
- Tested bcftools with correct reference + `-A` flag: ❌ Segfault
- Tested bcftools with correct reference, no `-A`: ✅ Success (~57s, 95MB RAM, 1.3MB gVCF)
- Validated gVCF format: Proper MIN_DP reference blocks and variant calls present
- Tested merge on 3-cell subset: ✅ Success (29,243 records, multi-sample format)

### Fixes Implemented

**1. `scripts/CapGTA/sc_variant_calling_bcftools_array.sh`**
   - Removed `-A` flag from `bcftools call` (line 65)
   - Added `-p short` partition for faster scheduling
   - Runtime: ~30 seconds per cell (well under 2-hour limit)

**2. `scripts/CapGTA/sc_from_bam_gta.sh`**
   - Updated reference path to correct N2 genome (line 100)

**3. `CapGTA_PP.sh`**
   - Updated default reference in config parsing (line 36)
   - Updated hardcoded default reference (line 46)

**4. `scripts/utils/extract_sc_from_bam_array.sh`**
   - Added `-p short` partition
   - Reduced timeout from 24h to 4h (observed max: 1h 47m)

**5. Created `bin/resubmit_variant_calling.sh`**
   - Automated resubmission script for all 4 samples
   - Generates barcode files from existing DNA BAMs
   - Submits variant calling arrays and merge jobs with correct dependencies

### Jobs Resubmitted

Cleaned up 1000 empty gVCF files per sample, then resubmitted:
- UDI_5: Variant calling 47730222, Merge 47730223
- UDI_6: Variant calling 47730224, Merge 47730225
- UDI_7: Variant calling 47730226, Merge 47730227
- UDI_8: Variant calling 47730428, Merge 47730429

**Status:** Jobs running successfully on `short` partition with corrected parameters.

### Files Modified
- `scripts/CapGTA/sc_variant_calling_bcftools_array.sh`
- `scripts/CapGTA/sc_from_bam_gta.sh`
- `CapGTA_PP.sh`
- `scripts/utils/extract_sc_from_bam_array.sh`
- Created: `bin/resubmit_variant_calling.sh`

### Related Issues
- Issue #3: Fixed bcftools gVCF mode failures (removed incompatible `-A` flag) 

# New tasks (2/25/26)

I debugged and ran the CapGTA pipeline on 4 samples, as described above. I have been analyzing the data, and the RNA counts are extremely low. I suspect that only taking reads spanning splice junctions is leaving quite a bit of trnascript reads on the table. I need to test that hypothesis by examing coverage- is there enrichment over exons in the "DNA" bam files? Then, if so, I need to devise an algorithm to estimate the proportion of reads that are likely rna by comparing exonic read depth to background non-coding DNA read depth. There are scripts to create bigwigs and plot coverage tracks in scripts/CapWGS_QC/ that may be useful for examining enrichment.

---

## Plan for RNA Read Recovery (Feb 25, 2026)

### Background

**Current Pipeline Behavior:**
- STAR aligns all reads
- Reads split by splice junctions:
  - **RNA BAM**: Only reads with `N` in CIGAR (splice junctions) → high confidence transcripts
  - **DNA BAM**: All other reads (includes unspliced RNA reads)
- **Problem**: RNA counts extremely low (mean ~2,141 transcripts/cell vs reference ~43 genes detected/cell)
- **Hypothesis**: Many exonic RNA reads without splice junctions are being classified as DNA

### Phase 1: Investigate Exonic Enrichment in DNA BAMs (Diagnostic)

**1.1 Generate BigWig files for visual inspection**
- Use existing `scripts/CapWGS_QC/generate_bigwig_array.sh` to create coverage tracks for subset of DNA BAMs
- Visually inspect coverage over gene models using `scripts/CapWGS_QC/coverage_tracks.sh`

**1.2 Quantitative enrichment analysis**
- Create script: `scripts/CapGTA/analyze_exonic_enrichment.py`
- For each cell's DNA BAM, compute:
  - Exonic coverage (reads overlapping exons)
  - Intronic coverage
  - Intergenic (non-coding) coverage
  - Enrichment ratios (exonic/intergenic, exonic/intronic)
- Run on representative cells from UDI_5
- Output: Summary statistics showing whether DNA BAMs contain RNA signal

### Phase 2: Develop RNA Classification Algorithm (if enrichment confirmed)

**2.1 Statistical approach (preferred):**
- For each cell:
  1. Estimate background DNA coverage from intergenic regions
  2. Model expected DNA coverage in exons based on background
  3. Identify reads exceeding expected DNA coverage as likely RNA
  4. Use statistical test (e.g., Poisson test) to classify regions as DNA vs RNA

**2.2 Implementation:**
- Create: `scripts/CapGTA/classify_rna_reads.py`
- For each DNA BAM:
  - Calculate background (intergenic) coverage rate
  - For each exon, test if coverage significantly exceeds background
  - Extract reads from high-coverage exons as additional RNA reads
  - Merge with existing spliced RNA reads

**2.3 Alternative approach (if statistical fails):**
- Machine learning classifier using features:
  - Coverage depth relative to background
  - Proximity to splice junctions
  - Gene expression patterns from reference data

### Phase 3: Integrate into Pipeline (if algorithm works)

**3.1 Update pipeline scripts:**
- Create: `scripts/CapGTA/reclassify_rna_from_dna.sh`
- Update `CapGTA_PP.sh` to optionally run reclassification
- Add pipeline flag: `-e/--extract-exonic-rna` (optional, off by default)

**3.2 Recount RNA matrices:**
- Merge reclassified RNA reads with existing spliced RNA reads
- Regenerate count matrices using `featureCounts` or similar
- **Note**: When regenerating, use top 1200 barcodes instead of 1000 (was slightly conservative in cell selection)

**3.3 Validation:**
- Compare new counts to current counts and reference datasets
- Check if cell type classification improves
- Verify DNA variant calling not affected by read reclassification

### Immediate Next Steps

1. **Create quantitative enrichment analysis script** (Phase 1.2)
   - Run on UDI_5 representative cells
   - Determine if exonic enrichment exists in DNA BAMs

2. **If enrichment confirmed** → Proceed to Phase 2.1

3. **If no enrichment** → Investigate other causes:
   - Alignment parameters
   - GTF annotation mismatch
   - featureCounts parameters
   - Library quality issues

---

## Phase 1 Results: Exonic Enrichment Analysis (Feb 25, 2026)

### Analysis Completed

**Created scripts:**
- `bin/analyze_exonic_enrichment.py` - Python script to compute exonic/intronic/intergenic coverage
- `bin/enrichment_analysis_array.sh` - SLURM array job for parallel processing
- `bin/compile_enrichment_results.py` - Compile individual results into summary

**Analyzed:** 46 cells from worm_CapGTA_UDI_5 using parallel array job (Job ID: 48531452)

### Key Findings: **STRONG EXONIC ENRICHMENT CONFIRMED**

**Mean Read Distribution in DNA BAMs:**
- **Exonic: 45.8%** (±2.4%)
- Intronic: 22.0% (±1.7%)
- Intergenic: 32.2% (±2.8%)

**Mean Enrichment Ratios:**
- **Exonic/Intergenic: 1.44x** (±0.19x)
- **Exonic/Intronic: 2.10x** (±0.22x)
- **Genic/Intergenic: 2.13x** (±0.28x)

### Interpretation

**The DNA BAMs contain a massive RNA signal:**
1. **~46% of reads map to exons** - far higher than expected for random genomic DNA
2. Expected DNA coverage should be roughly proportional to genome composition (~2-3% exonic in C. elegans)
3. **15-20x enrichment over expected** suggests most exonic reads are likely RNA
4. This explains the extremely low RNA counts (only counting spliced reads = ~2,141 transcripts/cell)
5. By only taking reads with splice junctions, we're missing the vast majority of transcript reads

### Conclusion

**Phase 1 hypothesis CONFIRMED.** The current pipeline is severely underestimating RNA content by:
- Classifying unspliced RNA reads as "DNA"
- Only counting ~5-10% of actual RNA reads (those spanning introns)
- Missing reads from:
  - Short exons that don't span introns
  - Intronless genes
  - Reads that happen to fall within single exons

**Next step:** Proceed to Phase 2 - develop algorithm to reclassify RNA reads from DNA BAMs.

### Technical Notes

- 4 jobs timed out (30 min limit) - likely cells with very high coverage
- Should increase timeout to 1 hour for full dataset
- Results saved to: `results/worm_CapGTA_UDI_5/enrichment_summary.csv`
- Individual cell results in: `results/worm_CapGTA_UDI_5/enrichment_analysis/`

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