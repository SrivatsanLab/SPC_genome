# CapWGS
This repository provides a pipeline for preprocesssing CapWGS data, as well as the code required to reproduce all of the analysis in our paper: Capsule-based Whole Genome Sequencing and lineage tracing ([bioRxiv](https://www.biorxiv.org/content/10.1101/2025.03.14.643253v1))

### 1. Installation

```bash
git clone https://github.com/SrivatsanLab/SPC_genome
cd SPC_genome

# Create the environment
micromamba create -n spc_genome -f environment.yml

# Activate the environment
micromamba activate spc_genome
```

### 2. Configuration

Edit `config.yaml` with your data paths:
```yaml
data:
  base_dir: "../sc_PolE_novaseq"
  read1: "/path/to/sample_R1.fastq.gz"
  read2: "/path/to/sample_R2.fastq.gz"

reference:
  genome_dir: "/shared/biodata/reference/GATK/hg38/"

output:
  base_dir: "./output"
  sample_name: "my_sample"
```

See [`DATA_README.md`](DATA_README.md) for detailed data requirements.

### 3. Run Preprocessing

Process CapWGS data using the main pipeline script:

```bash
./CapWGS_PP.sh \
  -o sample_name \
  -1 /path/to/read1.fastq.gz \
  -2 /path/to/read2.fastq.gz \
  -r 1000000000 \
  -g /shared/biodata/reference/GATK/hg38/
```

**Required arguments:**
- `-o` Sample name (output prefix)
- `-1` Read 1 FASTQ file
- `-2` Read 2 FASTQ file
- `-r` Total read count
- `-g` Reference genome directory

**Optional arguments:**
- `-s` Scripts directory (default: `./SPC_genome`)
- `-O` Output directory (default: `.`)
- `-n` Number of chunks for parallelization (default: 500)
- `-t` Temporary directory (overrides mktemp)

### Analysis notebooks:
- `notebooks/K562_tree/` - K562 tree experiment analysis notebooks
  - `sc_analysis.ipynb` - Main single-cell CapWGS analysis
  - `bulk_spectra_analysis.ipynb` - Bulk WGS mutational spectra
- `notebooks/benchmarking/` - CapWGS benchmarking analysis
  - `benchmarking.ipynb` - Comparison with public scWGS datasets
- `species_mixing/` - Species mixing experiment analysis

### Repository Structure

```
SPC_genome/
├── CapWGS_PP.sh                   # Main preprocessing script
├── scripts/                       # Helper scripts for pipeline
│   ├── PP_array.sh               # Parallel preprocessing array job
│   ├── sc_var_array.sh           # Single-cell variant calling
│   ├── make_variant_anndata.py   # Create AnnData from VCF
│   └── ...
├── test/                          # Test suite
│   ├── test_runner.sh            # Main test runner
│   ├── unit/                     # Unit tests
│   ├── integration/              # Integration tests
│   └── fixtures/                 # Test data
├── amplicon_analysis/             # Amplicon sequencing scripts
├── bin/                           # Generated files (indices, cell lists, metadata)
│   └── benchmarking/             # Benchmarking metadata files
├── data/                          # Project data (gitignored)
│   ├── K562_tree/                # K562 tree experiment data
│   │   ├── raw/                  # Raw FASTQ files
│   │   └── test/                 # Test data
│   └── benchmarking/             # Benchmarking datasets
│       └── aligned/              # Aligned public scWGS data
├── docs/summaries/                # Documentation summaries (gitignored)
├── notebooks/                     # Analysis notebooks
│   ├── K562_tree/                # K562 tree experiment analysis
│   │   ├── sc_analysis.ipynb    # Main single-cell analysis
│   │   └── bulk_spectra_analysis.ipynb
│   └── benchmarking/             # CapWGS benchmarking analysis
│       └── benchmarking.ipynb   # Public data comparison
├── results/                       # Generated outputs (gitignored)
│   ├── K562_tree/                # K562 tree results
│   └── benchmarking/             # Benchmarking results
├── environment.yml                # Micromamba environment
├── config.yaml                    # Configuration file
├── setup.sh                       # Setup verification script
└── DATA_README.md                 # Data organization guide
```
