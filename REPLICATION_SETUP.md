# Data Replication Setup - PREPROCESSING MODE

This document describes the setup for **preprocessing pipeline testing**. All processed data symlinks have been removed, and only raw input data is available. The goal is to regenerate all processed outputs using the preprocessing pipeline.

## Setup Summary

All key data files from `/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq/` have been symlinked into the `data/` directory structure.

## Current Data Structure

```
data/
└── raw/                                   # Raw input FASTQs (symlinked)
    ├── k562_tree_1_S1_R1_001.fastq.gz    # Read 1 (1.8 TB)
    └── k562_tree_1_S1_R2_001.fastq.gz    # Read 2 (1.8 TB)
```

**Note**: All processed data has been removed. The preprocessing pipeline will generate:
- `sc_outputs/` - Per-cell BAM and VCF files
- `vcf_chunks/` - VCF chunks
- `*.vcf.gz` - Joint-called variants
- `*.h5ad` - AnnData objects
- `real_cells.txt` - Detected cell barcodes

## Configuration

The `config.yaml` file has been updated with raw input paths:

```yaml
data:
  base_dir: "./data"
  read1: "./data/raw/k562_tree_1_S1_R1_001.fastq.gz"
  read2: "./data/raw/k562_tree_1_S1_R2_001.fastq.gz"
  read_count: 9000000000  # Estimate
```

## Running the Analysis

### Option 1: Use existing environment (recommended)

```bash
# Activate your existing jupyter environment
micromamba activate /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/envs/default_jupyter

# Start JupyterLab
jupyter lab

# Open and run sc_analysis.ipynb
```

### Option 2: Use the new local environment

The repository includes `environment.yml` and `setup.sh` for creating a new environment. Note: Setup may be in progress or you may need to complete it.

```bash
# If not already done, complete the setup
./setup.sh

# Activate the local environment
micromamba activate ./env

# Start JupyterLab
jupyter lab
```

## Analysis Notebooks

The repository includes the following analysis notebooks:

1. **`sc_analysis.ipynb`** - Main single-cell analysis notebook
   - Uses: `data/processed/anndata/sc_PolE_filtered.h5ad`
   - Contains: QC, filtering, mutational spectra, phylogenetic analysis

2. **`bulk_spectra_analysis.ipynb`** - Bulk WGS analysis
   - Uses: Bulk sequencing data
   - Contains: Mutational signature analysis

## Quick Start for Replication

```python
import scanpy as sc
import numpy as np
import pandas as pd

# Load the filtered AnnData object
adata = sc.read_h5ad('data/processed/anndata/sc_PolE_filtered.h5ad')

# View basic info
print(f"Cells: {adata.n_obs}")
print(f"Variants: {adata.n_vars}")
print(f"Observations: {adata.obs.columns.tolist()}")
print(f"Variables: {adata.var.columns.tolist()}")

# Continue with your analysis...
```

## Data Files

All symlinked files are read-only links to your original analysis directory. No data has been copied, so:

- Changes to files in either location will affect both
- Deleting a symlink won't delete the original file
- No additional disk space is used

## Verification

To verify the symlinks are working:

```bash
# Check symlinks
ls -lh data/processed/anndata/

# Check file accessibility
file data/processed/anndata/sc_PolE_filtered.h5ad

# Count VCF chunks (should be ~800 directories)
ls data/vcf_chunks/ | wc -l
```

## Notes

- The main filtered AnnData file is 55 GB, so loading may take 1-2 minutes
- Use `sc_PolE_filtered_no_singletons.h5ad` (2.7 GB) for faster testing
- All paths in notebooks should work with relative paths from repository root
- Environment contains all required packages: scanpy, pandas, matplotlib, etc.

## Next Steps

1. Open `sc_analysis.ipynb` in JupyterLab
2. Update any hardcoded paths to use the new `data/` structure
3. Run cells to replicate your analysis
4. Compare outputs with original results in `../sc_PolE_novaseq/figures/`

## Troubleshooting

**If notebooks can't find data:**
- Check you're in the repository root: `pwd` should show `.../SPC_genome`
- Verify symlinks: `ls -lh data/processed/anndata/`

**If Python packages are missing:**
- Activate the environment: `micromamba activate /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/envs/default_jupyter`
- Or use the local environment after completing setup

**If loading is too slow:**
- Use the smaller filtered file without singletons for testing
- Consider working on a compute node instead of rhino
