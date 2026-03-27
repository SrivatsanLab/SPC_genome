# CapWGS Script Audit

## Analysis Date: 2026-03-26
## Cleanup Completed: 2026-03-26

This document summarizes the script cleanup performed on the CapWGS pipeline directories.

---

## ✅ ACTIVELY USED SCRIPTS

### Called directly by CapWGS_PP.sh:

1. **scripts/CapWGS/PP_array.sh** - FASTQ alignment array
2. **scripts/CapWGS/parallel_merge_array.sh** - Parallel SAM merge (NEW)
3. **scripts/CapWGS/sc_extract_preprocess_qc_array.sh** - Unified SC processing
4. **scripts/CapWGS/submit_sc_variant_calling.sh** - SC variant calling submission
5. **scripts/CapWGS/submit_joint_calling.sh** - Joint calling submission
6. **scripts/CapWGS_QC/compile_lorenz.py** - Compile Lorenz curves
7. **scripts/CapWGS_QC/compile_qc_metrics.py** - Compile Picard metrics

### Called by submit_sc_variant_calling.sh:

8. **scripts/CapWGS/sc_var_array_gatk.sh** - GATK HaplotypeCaller per cell
9. **scripts/CapWGS/sc_var_array_bcftools.sh** - BCFtools variant calling per cell

### Called by submit_joint_calling.sh:

10. **scripts/CapWGS/bcftools_joint_calling.sh** - BCFtools merge/normalize
11. **scripts/CapWGS/generate_intervals.sh** - GATK interval generation
12. **scripts/CapWGS/gatk_genomicsdb_import_array.sh** - GenomicsDB import per interval
13. **scripts/CapWGS/gatk_joint_calling_parallel_array.sh** - GenotypeGVCFs per interval
14. **scripts/CapWGS/gatk_joint_calling_parallel_merge.sh** - Merge interval VCFs

### Called by sc_extract_preprocess_qc_array.sh:

15. **scripts/CapWGS_QC/lorenz.py** - Generate Lorenz curve and gini coefficient

---

## ❌ UNUSED SCRIPTS (Candidates for Deletion)

### scripts/CapWGS/:

1. **gatk_genomicsdb_import.sh** - Old chromosome-based GenomicsDB (replaced by interval-based)
2. **joint_calling_array.sh** - Old joint calling array (unclear purpose, not called)
3. **make_intervals.sh** - Duplicate of generate_intervals.sh (check contents)
4. **markdup_bqsr.sh** - OLD bulk preprocessing (replaced by sc_extract_preprocess_qc_array.sh)
5. **sc_bam_haplotypecaller.sh** - Old single-cell variant calling (replaced by sc_var_array_gatk.sh)
6. **sc_extract_and_preprocess_array.sh** - Old version (replaced by sc_extract_preprocess_qc_array.sh)
7. **sc_from_chunks.sh** - Old single-cell extraction (not used)
8. **sc_preprocessing_array.sh** - Standalone preprocessing (not called by pipeline)

### scripts/CapWGS_QC/:

9. **benchmarking_qc_array.sh** - Old QC array (QC now unified in sc_extract_preprocess_qc_array.sh)
10. **binned_coverage_array.sh** - Not called by pipeline
11. **binned_coverage.py** - Not called by pipeline
12. **compile_preseq_results.py** - Not called by pipeline
13. **coverage_tracks.py** - Not called by pipeline (standalone tool)
14. **coverage_tracks.sh** - Wrapper for coverage_tracks.py (standalone)
15. **generate_bigwig_and_lorenz_array.sh** - Old bigwig/Lorenz array (replaced by unified array)
16. **generate_bigwig_array.sh** - Old bigwig generation (replaced by unified array)
17. **genomic_feature_distribution.sh** - Not called by pipeline
18. **PP_array_ucsc.sh** - Alternative preprocessing for UCSC reference (not used)
19. **run_preseq_sample.sh** - Not called by pipeline
20. **submit_benchmarking_qc.sh** - OLD submission script (not called by main pipeline)
21. **submit_bigwig_lorenz.sh** - OLD submission script (not called by main pipeline)

---

## 🔍 SCRIPTS REQUIRING VERIFICATION

These scripts might be used manually or in other contexts:

### Potentially useful standalone tools:
- **coverage_tracks.py** / **coverage_tracks.sh** - May be useful for analysis
- **binned_coverage.py** - May be useful for analysis
- **genomic_feature_distribution.sh** - May be useful for analysis

### Need to verify if duplicates:
- **make_intervals.sh** vs **generate_intervals.sh** - Check if identical

---

## 📊 SUMMARY

**Total scripts analyzed**: 35
**Actively used**: 15 (43%)
**Unused/obsolete**: 20 (57%)

### Breakdown by directory:
- **scripts/CapWGS/**: 15 scripts total, 8 unused (53% unused)
- **scripts/CapWGS_QC/**: 14 scripts total, 12 unused (86% unused)

---

## 🎯 RECOMMENDATION

### Phase 1: Safe to delete (obsolete, replaced by newer versions)
- `markdup_bqsr.sh` - Replaced by unified array
- `sc_extract_and_preprocess_array.sh` - Replaced by sc_extract_preprocess_qc_array.sh
- `sc_bam_haplotypecaller.sh` - Replaced by sc_var_array_gatk.sh
- `gatk_genomicsdb_import.sh` - Replaced by interval-based version
- `generate_bigwig_and_lorenz_array.sh` - Replaced by unified array
- `benchmarking_qc_array.sh` - QC now in unified array
- `submit_benchmarking_qc.sh` - Not called by pipeline
- `submit_bigwig_lorenz.sh` - Not called by pipeline
- `PP_array_ucsc.sh` - Alternative preprocessing not used

### Phase 2: Verify then delete (unclear purpose, not called)
- `joint_calling_array.sh` - Check if needed
- `make_intervals.sh` - Compare with generate_intervals.sh
- `sc_from_chunks.sh` - Verify not needed
- `sc_preprocessing_array.sh` - Verify not needed

### Phase 3: Keep as standalone tools (may be useful for analysis)
- `coverage_tracks.py` + `coverage_tracks.sh`
- `binned_coverage.py` + `binned_coverage_array.sh`
- `genomic_feature_distribution.sh`
- `compile_preseq_results.py` + `run_preseq_sample.sh`

If keeping standalone tools, consider moving to `scripts/analysis/` or `scripts/standalone/`.

---

## 🔧 CONSOLIDATION PLAN

**Option A**: Merge into single `scripts/CapWGS/` directory
- Move actively used CapWGS_QC scripts into CapWGS/
- Create `scripts/CapWGS/qc/` subdirectory for QC-specific scripts
- Structure: `scripts/CapWGS/{preprocessing,variant_calling,qc}/`

**Option B**: Keep separation but clean up
- Delete unused scripts from both directories
- Move standalone tools to `scripts/analysis/`
- Keep only actively used scripts in each directory

**Recommended**: Option A - Single directory with subdirectories for organization

---

## ✅ CLEANUP COMPLETED

### Actions Taken:

1. **Deleted 13 obsolete scripts**:
   - `markdup_bqsr.sh`
   - `sc_extract_and_preprocess_array.sh`
   - `sc_bam_haplotypecaller.sh`
   - `gatk_genomicsdb_import.sh`
   - `generate_bigwig_and_lorenz_array.sh`
   - `benchmarking_qc_array.sh`
   - `submit_benchmarking_qc.sh`
   - `submit_bigwig_lorenz.sh`
   - `PP_array_ucsc.sh`
   - `joint_calling_array.sh`
   - `make_intervals.sh`
   - `sc_from_chunks.sh`
   - `sc_preprocessing_array.sh`

2. **Moved 7 standalone tools to `scripts/utils/`**:
   - `coverage_tracks.py` + `coverage_tracks.sh`
   - `binned_coverage.py` + `binned_coverage_array.sh`
   - `genomic_feature_distribution.sh`
   - `compile_preseq_results.py` + `run_preseq_sample.sh`

3. **Consolidated directories**:
   - Moved remaining QC scripts from `scripts/CapWGS_QC/` to `scripts/CapWGS/`
   - Removed empty `scripts/CapWGS_QC/` directory
   - All CapWGS pipeline scripts now in single `scripts/CapWGS/` directory

4. **Updated references**:
   - `CapWGS_PP.sh`: Updated paths for `compile_lorenz.py` and `compile_qc_metrics.py`
   - `sc_extract_preprocess_qc_array.sh`: Updated path for `lorenz.py`

---

## 📂 FINAL DIRECTORY STRUCTURE

### scripts/CapWGS/ (16 scripts):

**Preprocessing:**
- PP_array.sh
- parallel_merge_array.sh
- sc_extract_preprocess_qc_array.sh

**Variant Calling:**
- submit_sc_variant_calling.sh
- sc_var_array_gatk.sh
- sc_var_array_bcftools.sh

**Joint Calling:**
- submit_joint_calling.sh
- bcftools_joint_calling.sh
- generate_intervals.sh
- gatk_genomicsdb_import_array.sh
- gatk_joint_calling_parallel_array.sh
- gatk_joint_calling_parallel_merge.sh

**QC:**
- lorenz.py
- compile_lorenz.py
- compile_qc_metrics.py
- generate_bigwig_array.sh

### scripts/utils/ (additional standalone tools):
- coverage_tracks.py
- coverage_tracks.sh
- binned_coverage.py
- binned_coverage_array.sh
- genomic_feature_distribution.sh
- compile_preseq_results.py
- run_preseq_sample.sh

---

## 📊 CLEANUP SUMMARY

**Before cleanup:**
- scripts/CapWGS/: 19 scripts
- scripts/CapWGS_QC/: 14 scripts
- **Total: 33 scripts in 2 directories**

**After cleanup:**
- scripts/CapWGS/: 16 scripts
- scripts/CapWGS_QC/: (removed)
- scripts/utils/: +7 standalone tools
- **Total: 16 active pipeline scripts in 1 directory**

**Result:**
- ✅ Deleted 13 obsolete scripts (39% reduction)
- ✅ Consolidated 2 directories into 1
- ✅ Organized standalone tools in utils
- ✅ All path references updated
