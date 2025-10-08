# Preprocessing Pipeline Test Results

## Test Run Summary

**Date**: October 8, 2025
**Test Dataset**: 10 million read subset
**Status**: Partial success - identified bugs

## Bugs Found and Fixed

### 1. Missing `-r` option in getopts (Line 39)
**Issue**: The `-r` (read count) option was documented but not included in the getopts string.
**Fix**: Changed `getopts ":o:1:2:g:s:O:n:t:h"` to `getopts ":o:1:2:g:r:s:O:n:t:h"`

### 2. `wc -l` command includes filename (Lines 95, 118)
**Issue**: `wc -l file` returns "count filename" instead of just count, causing sbatch array syntax error.
**Fix**: Changed to `wc -l < file` to get only the count.

### 3. Variable name inconsistency (Multiple lines)
**Issue**: Mixed use of `SCRIPTS_DIR` and `scripts_DIR` throughout WGS_PP.sh
**Fix**: Standardized all references to use `SCRIPTS_DIR`

### 4. Incorrect scripts directory path (scr/PP_array.sh)
**Issue**: PP_array.sh expects project root directory but receives `./scr` path
- Script looks for `${scripts_DIR}/bin/atrandi_demux.py`
- But bin is at `./bin/`, not `./scr/bin/`

**Status**: NOT YET FIXED - needs design decision:
- Option A: Change WGS_PP.sh to pass project root instead of `./scr`
- Option B: Change PP_array.sh to expect `./scr` and adjust bin path

## Test Progress

### Completed Steps:
1. ✅ Created 10M read subset (609MB R1, 579MB R2)
2. ✅ FASTQ chunking (50 chunks created successfully)
3. ✅ Array job submission (job 36228894, 50 tasks)
4. ✅ Array jobs started

### Failed Steps:
1. ❌ Demultiplexing - script not found at expected path
2. ❌ Subsequent processing steps not reached

## Files Generated

```
test_output/
├── chunk_indices.txt           # 50 chunk IDs
├── barcodes.map               # Empty (expected at this stage)
└── sc_ouputs/                 # Empty directory

SLURM_outs/array_outs/
└── PP_array_36228894_*.out    # 50 log files
```

## Next Steps

1. Fix the scripts directory path issue
2. Rerun test to complete preprocessing
3. Verify outputs match expected structure
4. Compare with original run results

## Recommendations

### For Full Pipeline Run:
- Default temp directory is now `/hpc/temp/srivatsan_s/SPC_genome_preprocessing`
- This directory has been pre-created and is ready to use
- ~4TB of space will be used for FASTQ chunks during processing
- Monitor disk space: `du -sh /hpc/temp/srivatsan_s/SPC_genome_preprocessing`
- Clean up after completion: `rm -rf /hpc/temp/srivatsan_s/SPC_genome_preprocessing/*`

### Pipeline Improvements Needed:
1. Better error handling for missing scripts/directories
2. Input validation for directory paths
3. Check for required dependencies before starting
4. Progress reporting during long-running steps

## Bug Fix Commits

All fixes committed:
- ✅ WGS_PP.sh: Add `-r` to getopts
- ✅ WGS_PP.sh: Fix `wc -l` commands
- ✅ WGS_PP.sh: Standardize SCRIPTS_DIR variable
- ✅ WGS_PP.sh: Fix bin path issue (pass project root instead of ./scr)
- ✅ WGS_PP.sh: Set default temp directory to /hpc/temp/srivatsan_s/SPC_genome_preprocessing
- ✅ config.yaml: Updated with new temp directory default
- ✅ Documentation: Updated all references to temp directory
