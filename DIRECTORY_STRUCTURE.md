# Directory Structure

This document clarifies the directory structure used by the SPC_genome pipeline.

## Overview

```
SPC_genome/
├── bin/                           # Metadata files (indices, cell lists)
├── data/                          # All sequencing data (raw + processed)
├── results/                       # Analysis outputs only (.h5ad, CSVs, plots)
├── scripts/                       # Pipeline scripts
└── notebooks/                     # Analysis notebooks
```

## Detailed Structure

### bin/ - Metadata and Indices
**Purpose:** Generated metadata files that describe the data

```
bin/
├── test/                          # Test run metadata
│   ├── chunk_indices.txt         # List of chunk identifiers
│   ├── real_cells.txt            # List of detected cell barcodes
│   ├── sam_list.txt              # List of SAM files to merge
│   ├── readcounts.csv            # Reads per barcode
│   ├── kneeplot.png              # Cell detection plot
│   └── barcodes.map              # Barcode to VCF mapping
└── [sample_name]/                # Production run metadata
```

### data/ - All Sequencing Data
**Purpose:** Raw FASTQ files AND processed BAM/VCF files

```
data/
├── K562_tree/
│   ├── raw/                      # Raw sequencing data
│   │   ├── sample_R1.fastq.gz
│   │   └── sample_R2.fastq.gz
│   ├── test/                     # Test data
│   │   ├── test_R1.fastq.gz     # Test FASTQ files
│   │   ├── test_R2.fastq.gz
│   │   └── aligned/              # Processed test data
│   │       ├── sample.bam        # Bulk merged BAM
│   │       ├── sample.sam        # Bulk merged SAM
│   │       ├── BARCODE1.bam      # Single-cell BAMs
│   │       ├── BARCODE1.g.vcf.gz # Single-cell VCFs
│   │       └── ...
│   └── aligned/                  # Production aligned data
│       ├── sample.bam
│       └── BARCODE*.vcf.gz
└── benchmarking/
    └── aligned/                  # Public dataset aligned files
```

### results/ - Analysis Outputs ONLY
**Purpose:** Final outputs from analysis notebooks (.h5ad, CSVs, plots, trees)

```
results/
├── K562_tree/
│   ├── sample.h5ad               # AnnData object
│   ├── sc_spectra.csv            # Mutational spectra
│   ├── sc_test.newick            # Phylogenetic tree
│   ├── snp_counts.csv            # Variant counts
│   └── plots/                    # Generated plots
└── benchmarking/
    └── benchmarking_results.csv
```

### scripts/ - Pipeline Scripts
**Purpose:** Executable scripts for the pipeline

```
scripts/
├── PP_array.sh                   # Parallel preprocessing
├── concatenate.sh                # Merge and cell detection
├── sc_from_chunks.sh             # Single-cell extraction
├── sc_var_array.sh               # Variant calling
├── submit_sc_variant_calling.sh  # Wrapper for variant calling
└── *.py                          # Python helper scripts
```

## Key Principles

1. **bin/** = Metadata about the data (lists, indices, mappings)
2. **data/** = All sequencing data (raw FASTQ + aligned BAM/VCF)
3. **results/** = Analysis outputs from notebooks (.h5ad, CSVs, plots)
4. **Intermediate files** (temp SAMs, unsorted BAMs) go to `/hpc/temp/` and are deleted

## Pipeline I/O Flow

```
Input (data/):
  └─ FASTQ files

Processing (/hpc/temp/):
  └─ Chunked FASTQs → Aligned SAMs → [deleted after merge]

Metadata (bin/):
  └─ chunk_indices.txt, real_cells.txt, etc.

Aligned Data (data/aligned/):
  └─ Bulk BAM, Single-cell BAMs, VCFs

Analysis (notebooks → results/):
  └─ h5ad files, CSVs, plots, trees
```

## Command Line Usage

When running the pipeline, the `-O` flag is currently used to determine test vs production mode:

```bash
# Test mode (uses bin/test/ and data/K562_tree/test/aligned/)
./CapWGS_PP.sh -O "results/test" ...

# Production mode (uses bin/[sample]/ and data/[sample]/aligned/)
./CapWGS_PP.sh -O "results/[sample]" ...
```

Note: Despite the flag name, outputs go to `data/`, not `results/`. The flag is used for mode detection only.
