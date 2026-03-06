# CapWGS Pipeline Development Notes

## Per-Chromosome GenomicsDB Implementation (LATEST)

### Date: 2026-03-06

**Problem:** GenomicsDB import was a severe bottleneck. Single job importing all chromosomes for all samples took 26+ hours and timed out.

**Solution:** Parallelized GenomicsDB import by chromosome, creating 24 independent databases instead of one monolithic database.

### Architecture Change

**Before:**
```
Single GenomicsDB import (all chromosomes) → Per-chromosome GenotypeGVCFs
     [BOTTLENECK: 26+ hours, TIMEOUT]        [Fast: 24 parallel jobs]
```

**After:**
```
Per-chromosome GenomicsDB import → Per-chromosome GenotypeGVCFs → Merge VCFs
     [24 parallel jobs, 2-4 hours]      [24 parallel jobs]           [1 job]
```

### Implementation

**Created:**
- `scripts/bulk/gatk_genomicsdb_import_array.sh` - Per-chromosome GenomicsDB for bulk pipeline
- `scripts/CapWGS/gatk_genomicsdb_import_array.sh` - Per-chromosome GenomicsDB for CapWGS pipeline
- `bin/submit_K562_joint_calling.sh` - Bulk joint calling submission script (moved from scripts/)

**Modified:**
- `scripts/bulk/gatk_joint_calling_parallel_array.sh` - Updated to read per-chromosome GenomicsDB
- `scripts/CapWGS/gatk_joint_calling_parallel_array.sh` - Updated to read per-chromosome GenomicsDB
- `scripts/CapWGS/submit_joint_calling.sh` - Now submits per-chromosome GenomicsDB array
- `scripts/bulk/submit_parallel_joint_calling.sh` - Moved to `bin/submit_K562_joint_calling.sh`

**Key Changes:**

1. **GenomicsDB structure changed from monolithic to per-chromosome:**
   - Before: `genomicsdb_${COHORT_NAME}/` (single database with all chromosomes)
   - After: `genomicsdb/${CHROMOSOME}/` (24 separate databases, e.g., `genomicsdb/chr1/`)

2. **Array script arguments updated:**
   - Before: `GENOMICS_DB` (single path)
   - After: `GENOMICSDB_DIR` (base directory containing per-chromosome databases)
   - Script constructs path: `GENOMICS_DB="${GENOMICSDB_DIR}/${CHROMOSOME}"`

3. **Job submission now has 3-step chain:**
   - Step 1: Per-chromosome GenomicsDB import array (24 jobs, 4 hours each)
   - Step 2: Per-chromosome joint calling array (24 jobs, depends on step 1)
   - Step 3: Merge VCFs (1 job, depends on step 2)

4. **Each chromosome processes independently:**
   - GenomicsDB import only reads that chromosome's data from GVCFs
   - Dramatically reduces I/O and memory requirements per job
   - Enables true parallelization of the bottleneck step

### Performance Impact

- **GenomicsDB import:** 26+ hours (timeout) → 2-4 hours expected (24x parallelization)
- **Total pipeline time:** Significantly reduced, no more timeout issues
- **Resource efficiency:** Better cluster utilization with smaller parallel jobs
- **Scalability:** Handles large cohorts (17+ samples) that previously timed out

### Testing Status

- **K562 bulk samples (17 samples):**
  - Job 48940067: GenomicsDB array (24 chromosomes) - RUNNING
  - Job 48940068: Joint calling array (24 chromosomes) - PENDING
  - Job 48940069: Merge VCFs - PENDING

- **HSC2_enzyme CapWGS (80 cells expected):**
  - Job 48940070: Preprocessing - RUNNING
  - Joint calling will use new per-chromosome GenomicsDB when preprocessing completes

### Deprecated

The old single-database approach is no longer used:
- `scripts/bulk/gatk_genomicsdb_import.sh` - Replaced by array version
- `scripts/CapWGS/gatk_genomicsdb_import.sh` - Replaced by array version

---

## Parallel Joint Calling Implementation

### Background

**Current CapWGS GATK joint calling architecture:**
1. Generates 1Mb intervals across entire genome (thousands of intervals)
2. For each interval: runs GenomicsDBImport + GenotypeGVCFs in parallel
3. Outputs per-interval BCF files to temp directory
4. **NO merge step** - ends up with thousands of separate BCF files

**Issues with current approach:**
- GenomicsDBImport runs thousands of times (once per 1Mb interval) - massive I/O overhead
- Creates thousands of small GenomicsDB workspaces instead of one unified database
- No final merge step to create a single joint VCF
- Each 1Mb interval imports all GVCFs independently (very inefficient)

### Proposed Architecture

**Goal**: Match the efficient bulk pipeline approach (chromosome-based parallelization)

**New Architecture** (similar to bulk pipeline):

1. **Single GenomicsDB import** (one-time, run once)
   - Import all single-cell GVCFs into one GenomicsDB
   - Use GATK `wgs_calling_regions.hg38.interval_list` to restrict regions
   - Job: `scripts/CapWGS/gatk_genomicsdb_import.sh`

2. **Parallel per-chromosome GenotypeGVCFs** (24 array jobs)
   - Extract chromosome list from interval list (chr1-22, X, Y)
   - Call variants on each chromosome in parallel
   - Job: `scripts/CapWGS/gatk_joint_calling_parallel_array.sh`

3. **Merge per-chromosome VCFs** (single job after array completes)
   - Concatenate 24 chromosome VCFs into final joint VCF
   - Normalize and generate statistics
   - Job: `scripts/CapWGS/gatk_joint_calling_parallel_merge.sh`

### Implementation Plan

#### New Scripts to Create:

**1. `scripts/CapWGS/gatk_genomicsdb_import.sh`**
- Create GenomicsDB from barcodes.map
- Output to `${RESULTS_DIR}/genomicsdb_${SAMPLE_NAME}`
- SLURM config: 12 hours, 64GB RAM, 8 cores
- Auto-generate barcodes.map from GVCFs if not provided
- Robust cleanup of incomplete GenomicsDB

**2. `scripts/CapWGS/gatk_joint_calling_parallel_array.sh`**
- Per-chromosome GenotypeGVCFs
- Read chromosome from chromosome_list.txt using SLURM_ARRAY_TASK_ID
- Output to `${RESULTS_DIR}/vcfs/per_chromosome/${CHROMOSOME}.vcf.gz`
- SLURM config: 2 days, 32GB RAM, 4 cores
- Array job (1-24 for each chromosome)

**3. `scripts/CapWGS/gatk_joint_calling_parallel_merge.sh`**
- Concatenate per-chromosome VCFs with bcftools concat
- Normalize and generate statistics
- Output final `${RESULTS_DIR}/${SAMPLE_NAME}_joint.vcf.gz`
- SLURM config: 12 hours, 64GB RAM, 16 cores
- Generate variant summary statistics

#### Modified Scripts:

**1. `scripts/CapWGS/submit_joint_calling.sh`** (lines 61-163)
Replace current GATK mode logic with parallelized approach:
- Extract chromosome list from `wgs_calling_regions.hg38.interval_list` (NOT from make_intervals.sh)
  ```bash
  grep -v "^@" "${INTERVALS}" | cut -f1 | sort -u > "${CHROMOSOME_FILE}"
  ```
- Submit GenomicsDB import job
- Submit chromosome array job with dependency on GenomicsDB: `--dependency=afterok:${genomicsdb_job_id}`
- Submit merge job with dependency on array job: `--dependency=afterok:${array_job_id}`

**2. Retain for BCFtools mode:**
- Keep existing bcftools logic (lines 39-60) unchanged
- `bcftools_joint_calling.sh` works fine as-is

#### Scripts to Deprecate:
- `scripts/CapWGS/joint_calling_array.sh` - replaced by chromosome-based approach
- `scripts/CapWGS/make_intervals.sh` - no longer needed (using chromosomes instead of 1Mb windows)

### Benefits

1. **90% reduction** in GenomicsDB operations (1 vs thousands)
2. **Much faster**: 24 parallel chromosome jobs instead of thousands of tiny intervals
3. **Single merged VCF** output instead of thousands of BCF fragments
4. **Consistent** with bulk pipeline architecture
5. **Easier to monitor** and debug (24 jobs vs thousands)
6. **Follows GATK best practices**: Uses curated calling regions interval list

### GATK Calling Regions and Alt Contigs

**Reference Information:**
- GRCh38/hg38 full reference: 3,366 contigs (includes alt contigs, unplaced scaffolds)
- GATK `wgs_calling_regions.hg38.interval_list`: 24 chromosomes only (chr1-22, X, Y)
- **Note**: chrM is excluded from GATK calling regions

**Alt Contigs Explained:**
- Alt contigs are alternative haplotypes with very long flanking sequences nearly identical to primary assembly
- BWA-MEM runs in alt-aware mode when `.alt` file exists in index
- Alt-aware mode prevents alt contigs from degrading mapping quality for primary assembly reads
- Reads preferentially map to primary assembly; alt contig hits become supplementary alignments
- GATK **intentionally excludes** alt contigs from variant calling (best practices)

**Why exclude alt contigs from variant calling:**
- Standard variant callers cannot take advantage of alt-aware mapping
- Including alt contigs causes reduced variant calling sensitivity due to zero mapping quality assignments
- Alt contigs help BWA-MEM assign better mapping quality scores, but variants are only called on primary assembly

**Source**: [Heng Li's guide on human reference genomes](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)

### Implementation Status

- [x] Create `scripts/CapWGS/gatk_genomicsdb_import.sh`
- [x] Create `scripts/CapWGS/gatk_joint_calling_parallel_array.sh`
- [x] Create `scripts/CapWGS/gatk_joint_calling_parallel_merge.sh`
- [x] Modify `scripts/CapWGS/submit_joint_calling.sh` GATK mode
- [x] Modify `CapWGS_PP.sh` to pass TMP_DIR to submit_joint_calling.sh
- [x] Configure per-chromosome VCFs to write to temp directory
- [ ] Test full pipeline with small dataset
- [ ] Deprecate old scripts after validation

### File Locations Summary

**Permanent files (in results directory):**
- `${RESULTS_DIR}/genomicsdb_${OUTPUT_NAME}/` - GenomicsDB workspace
- `${RESULTS_DIR}/sample_map.txt` - Sample map for GenomicsDB
- `${RESULTS_DIR}/chromosome_list.txt` - List of chromosomes to process
- `${RESULTS_DIR}/${OUTPUT_NAME}_joint.vcf.gz` - Final merged joint VCF
- `${RESULTS_DIR}/${OUTPUT_NAME}_joint_stats.txt` - Variant statistics

**Temporary files (in temp directory):**
- `${TMP_DIR}/per_chromosome_vcfs/*.vcf.gz` - Per-chromosome VCFs (intermediate)
- `${TMP_DIR}/per_chromosome_vcfs/*_stats.txt` - Per-chromosome statistics
- Automatically cleaned up when temp directory is removed

### Created Files

**scripts/CapWGS/gatk_genomicsdb_import.sh**
- Arguments: GVCF_DIR, REFERENCE_DIR, GENOMICS_DB, SAMPLE_MAP
- Auto-generates sample map from GVCFs in sc_outputs/
- Removes incomplete GenomicsDB and recreates
- 12 hours, 64GB RAM, 8 cores

**scripts/CapWGS/gatk_joint_calling_parallel_array.sh**
- Arguments: GENOMICS_DB, OUTPUT_DIR, REFERENCE_DIR, CHROMOSOME_FILE
- Per-chromosome GenotypeGVCFs using SLURM_ARRAY_TASK_ID
- Generates per-chromosome statistics
- 2 days, 32GB RAM, 4 cores

**scripts/CapWGS/gatk_joint_calling_parallel_merge.sh**
- Arguments: PER_CHR_DIR, OUTPUT_VCF, CHROMOSOME_FILE
- Validates all per-chromosome VCFs exist before merging
- Uses bcftools concat for efficient chromosome-ordered concatenation
- Generates final statistics
- 12 hours, 64GB RAM, 16 cores

**Modified: scripts/CapWGS/submit_joint_calling.sh**
- Added TMP_DIR as 8th argument (defaults to /hpc/temp/srivatsan_s/joint_calling_${OUTPUT_NAME})
- GATK mode now uses 4-step process:
  1. Extract chromosomes from wgs_calling_regions.hg38.interval_list
  2. Submit GenomicsDB import job
  3. Submit per-chromosome array job (depends on GenomicsDB)
  4. Submit merge job (depends on array completion)
- Per-chromosome VCFs written to temp directory: `${TMP_DIR}/per_chromosome_vcfs/`
- Final merged VCF written to results directory: `${RESULTS_DIR}/${OUTPUT_NAME}_joint.vcf.gz`
- BCFtools mode unchanged

**Modified: CapWGS_PP.sh**
- Line 257: Added TMP_DIR as 8th argument to submit_joint_calling.sh call
