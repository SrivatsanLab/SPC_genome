# SPC Genome - Single Cell Genomics Pipeline

Computational pipeline for processing and analyzing Capsule-based Whole Genome Sequencing (CapWGS) data, with a focus on somatic mutation detection in PolE-P286R K562 cells.

## Overview

This repository contains the computational pipeline for the CapWGS method, enabling high-quality whole genome sequencing of single cells through microfluidic capsule-based cell isolation and linear amplification.

**Key components:**
- **Preprocessing pipeline:** Parallel processing of CapWGS data using SLURM arrays
- **Variant calling:** Single-cell and bulk variant detection using GATK
- **Benchmarking:** Comparison with existing scWGS methods (LIANTI, PTA)
- **Analysis notebooks:** Downstream analysis including mutational spectra and phylogenetic reconstruction

**Reference:** Capsule-based Whole Genome Sequencing enables high-quality single cell genomics ([bioRxiv](https://www.biorxiv.org/content/10.1101/2025.03.14.643253v1))

## Quick Start

### 1. Installation

**Prerequisites:**
- Access to Fred Hutch HPC cluster (Rhino/Gizmo)
- Micromamba package manager
- SLURM job scheduler

**Install the environment:**
```bash
git clone https://github.com/SrivatsanLab/SPC_genome
cd SPC_genome

# Create the environment
micromamba create -n spc_genome -f environment.yml

# Activate the environment
micromamba activate spc_genome
```

You can also run `./setup.sh` for a guided setup that checks dependencies.

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

### 4. Testing (Optional)

Run the test suite to verify your installation:

```bash
# Run all tests
./test/test_runner.sh

# Run only unit tests (fast, no SLURM jobs)
./test/test_runner.sh unit

# Run only integration tests (slower, submits SLURM jobs)
./test/test_runner.sh integration
```

See [`test/README.md`](test/README.md) for details on writing and running tests.

### 5. Run Analysis

Start JupyterLab to run analysis notebooks:
```bash
# On rhino node
micromamba activate spc_genome
jupyter lab
```

**Analysis notebooks:**
- `notebooks/K562_tree/` - K562 tree experiment analysis notebooks
  - `sc_analysis.ipynb` - Main single-cell CapWGS analysis
  - `bulk_spectra_analysis.ipynb` - Bulk WGS mutational spectra
- `notebooks/benchmarking/` - CapWGS benchmarking analysis
  - `benchmarking.ipynb` - Comparison with public scWGS datasets
- `species_mixing/` - Species mixing experiment analysis

## Repository Structure

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

## Pipeline Workflow

1. **FASTQ Splitting** - Split input FASTQs into chunks for parallel processing
2. **Alignment** - BWA-MEM alignment to reference genome
3. **Cell Detection** - Knee plot analysis to identify real cells
4. **Single Cell Extraction** - Extract reads per cell barcode
5. **Variant Calling** - GATK HaplotypeCaller per cell
6. **Joint Calling** - Combine variants across cells
7. **AnnData Generation** - Create AnnData object for analysis

## Requirements

### Computational Resources

**Recommended SLURM resources:**
- CPUs: 2-4 per job
- Memory: 16-32 GB per job
- Time: 2-4 days for full pipeline
- Partitions: `campus` or `short`

**Storage:**
- ~60-125 GB per sample
- Use `/fh/fast/` for active analysis
- Archive to S3 for long-term storage

### Software Dependencies

Core tools (installed via micromamba):
- Python 3.11+
- BWA, Samtools, BCFtools
- GATK4, Picard
- Scanpy, AnnData
- NumPy, Pandas, Matplotlib

See [`environment.yml`](environment.yml) for complete list.

## Data

**Not included in repository:**
- Raw sequencing data (`.fastq.gz`)
- Aligned BAM files
- VCF files
- Processed AnnData objects

This repository supports multiple datasets:
- **K562_tree**: PolE-P286R K562 phylogenetic tree experiment
- **benchmarking**: CapWGS comparison with public scWGS methods (LIANTI, PTA)
- **species_mixing**: Species mixing validation experiment

Data is organized in `data/` with subdirectories for each experiment (see [`DATA_README.md`](DATA_README.md))

## Cluster-Specific Information

This pipeline is optimized for the **Fred Hutch HPC cluster** (Rhino/Gizmo):

- Login nodes: `rhino01-03`
- Job scheduler: SLURM
- Reference genomes: `/shared/biodata/reference/GATK/`
- Module system: Lmod

For other HPC systems, you may need to:
- Adjust SLURM directives in scripts
- Update reference genome paths
- Modify temporary storage locations

## Troubleshooting

### Common Issues

**Job fails with "out of memory":**
```bash
# Increase memory allocation
sbatch --mem=32G script.sh
```

**Permission denied:**
```bash
chmod +x *.sh scripts/*.sh
```

**Module not found:**
```bash
# Ensure environment is activated
micromamba activate spc_genome
```

**Can't find data files:**
- Check paths in `config.yaml`
- Use absolute paths if needed
- See [`DATA_README.md`](DATA_README.md)

### Getting Help

- **Pipeline issues:** Open an issue on GitHub
- **Fred Hutch cluster:** Email scicomp@fredhutch.org
- **SLURM jobs:** Use `squeue -u $USER` to check status

## Citation

If you use this pipeline, please cite:
[Citation information to be added]

## License

See [`LICENSE`](LICENSE) for details.

## Contact

**Srivatsan Lab**
Fred Hutchinson Cancer Center

For questions about this repository, please open an issue on GitHub.
