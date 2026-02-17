# CapGTA Pipeline Development

## Initial Issues (Feb 11, 2026)

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