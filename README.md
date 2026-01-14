# SPC Genome Pipelines

Pipelines for processing single-cell whole genome sequencing (CapWGS) and genome-transcriptome co-assay (CapGTA) data. Includes tools for alignment, cell detection, variant calling, and AnnData generation.

## Installation

```bash
git clone https://github.com/SrivatsanLab/SPC_genome
cd SPC_genome

# Create the environment
micromamba create -n spc_genome -f environment.yml
micromamba activate spc_genome
```

## Pipelines

### CapWGS_PP.sh - Whole Genome Sequencing with Variant Calling

Complete pipeline for single-cell whole genome sequencing data:
- BWA-MEM alignment
- Cell detection from barcode counts
- Single-cell BAM extraction
- GATK variant calling (HaplotypeCaller + GenotypeGVCFs)
- AnnData generation with variant calls

**Output:** `data/{sample}/`, `data/{sample}/sc_outputs/`, `results/{sample}/variants.h5ad`

### CapWGS_PP_QC_only.sh - Coverage QC Pipeline

QC-only pipeline for coverage benchmarking (no variant calling):
- BWA-MEM alignment
- Cell detection
- Single-cell BAM extraction
- BigWig generation
- Lorenz curve analysis for coverage uniformity

**Output:** `data/{sample}/sc_outputs/`, `results/{sample}/` (QC metrics, Lorenz curves)

### CapGTA_PP.sh - Genome-Transcriptome Co-Assay

Pipeline for simultaneous genome and transcriptome profiling:
- STAR alignment (separates DNA/RNA by splice junctions)
- Dual-modality cell detection
- Single-cell DNA and RNA BAM extraction
- BCFtools variant calling on DNA
- RNA count matrix generation
- Dual AnnData output (variants + gene expression)

**Output:** `data/{sample}/sc_outputs/`, `results/{sample}/variants.h5ad`, `results/{sample}/rna_counts.h5ad`

## Usage Example

```bash
./CapWGS_PP.sh \
  -o sample_name \
  -1 /path/to/read1.fastq.gz \
  -2 /path/to/read2.fastq.gz \
  -r 1000000000 \
  -g /shared/biodata/reference/GATK/hg38/
```

**Required arguments:**
- `-o` Sample name (creates `data/{sample}/` and `results/{sample}/`)
- `-1` Read 1 FASTQ file(s)
- `-2` Read 2 FASTQ file(s)
- `-r` Total read count
- `-g` Reference genome directory

**Optional arguments:**
- `-s` Scripts directory (default: `.`)
- `-n` Number of chunks for parallelization (default: 500)
- `-t` Temporary directory (default: `/hpc/temp/srivatsan_s/SPC_genome_preprocessing/{sample}/`)
- `-h` Show help message

All pipelines support the same argument structure and can optionally use `config.yaml` for default values.

## Configuration (Optional)

Create `config.yaml` to set defaults:

```yaml
processing:
  n_chunks: 500
  tmp_dir: /hpc/temp/srivatsan_s/SPC_genome_preprocessing

reference:
  genome_dir: /shared/biodata/reference/GATK/hg38/

data:
  read1: /path/to/read1.fastq.gz
  read2: /path/to/read2.fastq.gz
  read_count: 1000000000

output:
  sample_name: my_sample
```

Command-line arguments override config values.

## Directory Structure

```
SPC_genome/
├── CapWGS_PP.sh              # Main WGS pipeline
├── CapWGS_PP_QC_only.sh      # QC-only pipeline
├── CapGTA_PP.sh              # Genome-transcriptome pipeline
├── scripts/                  # Pipeline scripts
│   ├── CapWGS/              # WGS + GATK scripts
│   ├── CapWGS_QC/           # Coverage QC scripts
│   ├── CapGTA/              # Genome-transcriptome scripts
│   ├── bulk/                # Bulk processing utilities
│   └── utils/               # Shared utilities
├── bin/                      # Experiment metadata and one-off scripts
├── data/                     # Alignments and single-cell outputs (gitignored)
│   └── {sample}/
│       ├── {sample}.bam     # Bulk alignment
│       └── sc_outputs/      # Single-cell BAMs, VCFs, bigwigs
├── results/                  # Final outputs (gitignored)
│   └── {sample}/
│       ├── variants.h5ad    # Variant AnnData
│       ├── rna_counts.h5ad  # RNA AnnData (CapGTA only)
│       └── *_qc_summary.csv # QC metrics
├── notebooks/                # Analysis notebooks
│   ├── K562_tree/           # K562 tree experiment analysis
│   ├── benchmarking/        # CapWGS benchmarking analysis
│   └── PolE_worm_pilot/     # C. elegans CapGTA analysis
└── environment.yml           # Micromamba environment

```

## Analysis Notebooks

- **K562_tree/** - K562 lineage tree experiment
  - `sc_analysis.ipynb` - Single-cell variant analysis
  - `bulk_spectra_analysis.ipynb` - Mutational spectra
- **benchmarking/** - CapWGS validation
  - `benchmarking.ipynb` - Comparison with public scWGS datasets
- **PolE_worm_pilot/** - C. elegans CapGTA analysis
  - `PolE_worm_pilot_analysis.ipynb` - Multi-modal analysis

## Citation

Capsule-based Whole Genome Sequencing and lineage tracing - [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.03.14.643253v1)
