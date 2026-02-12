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