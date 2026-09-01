# SPC Genome

This repository contains the code used for preprocessing and analysis of genomic data from our paper: Capsule-based Whole Genome Sequencing and lineage tracing - [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.03.14.643253v1)

Here you will find pipelines and tools for processing and analyzing capsule-based single-cell whole genome sequencing (CapWGS) and genome-transcriptome co-assay (CapGTA) data, as well as the code used to produce all data and figures for our paper. This includes jupyter notebooks for figures and statistical analyses, and scripts to process public data for benchmarking, and bulk WGS data from the K562 mutation accumulation experiment and HSC donor samples. 

## Setup

```bash
git clone https://github.com/SrivatsanLab/SPC_genome
cd SPC_genome

# Create the environment
micromamba create -n spc_genome -f environment.yml
micromamba activate spc_genome
```

#  Capsule-based single-cell whole genome sequencing (CapWGS)

### preprocessing: `CapWGS_PP.sh`

Complete pipeline for single-cell whole genome sequencing data:
- BWA-MEM alignment
- Cell detection from per-barcode read counts
- Single-cell BAM extraction
- GATK variant calling (HaplotypeCaller + GenotypeGVCFs)

**Outputs:**
- single cell bam files
- single cell .g.vcf files
- joint variant call set in .vcf format
- thorough per cell qc metrics, including lorenz curves

### Usage Example

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

### Configuration (Optional)

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

# capsule-based single-cell genome-transcriptome co-assay (CapGTA)

### Upstream preprocessing

1. [spc-demultiplex](https://github.com/SrivatsanLab/SPC-demultiplex)- for demultiplexing combinatorial index barcodes. Outputs per-barcode fastq files. 
2. [spc-align](https://github.com/SrivatsanLab/SPC-align)- using `STAR` splice aware alignment to WBcel235. Outputs per-barcode bam files.

Downstream pipelines (below) require a `real_cells.csv`, which has columns barcode,path. User should examine qc metrics from `spc-align` to determine which barcodes are associated with occupied capsules (real cells) to create this csv for the downstream gene expression and variant calling pipelines. 

### Gene expression preprocessing

`CapGTA_gex_soupx_pipeline.sh` runs the full spliced-RNA GEX pipeline end-to-end (SLURM-chained via `afterok`):

1. Spliced-read filter (CIGAR `N`) → featureCounts on real cells → per-cell integer count matrix.
2. Scanpy preprocessing: QC, WBID→gene-symbol rename, per-gene scaling by `1 / (n_junctions * exonic_length_kb)` to correct for junction count and cDNA-length amplification bias, HVG, PCA, UMAP, Leiden clustering.
3. Same count pipeline on a random subset of empty droplets to build an ambient-RNA soup profile — with a WGA-aware QC filter that drops empties dominated by a single amplified transcript.
4. `SoupX` ambient decontamination against real-cell clusters.
5. Re-preprocess on the decontaminated counts (with cuticle-gene drop) to get final clusters.

```bash
bash CapGTA_gex_soupx_pipeline.sh \
    results/<sample>/real_cells.csv \
    /shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf \
    results/<sample> \
    [--n-empties 5000] [--resolution 2]
```

**Required arguments:**
- `<real_cells.csv>` cells to process (barcode,bam_path from spc-align)
- `<annotation.gtf>` reference GTF (WBcel235 iGenomes)
- `<output_dir>` all outputs land here

**Key outputs** (under `<output_dir>/`):
- `adata_gex_spliced_soupx_processed.h5ad` — **final** h5ad: SoupX-decontaminated counts, preprocessed with clusters, ready for annotation
- `adata_gex_spliced.h5ad` — preprocessed pre-SoupX (useful for comparison)
- `rna_counts/` and `soupx_empties/rna_counts/` — raw featureCounts output for real cells and empties
- `soupx_out/soupx_summary.csv` — per-cell contamination fractions (`rho`)
- `gene_junctions.csv` / `gene_lengths.csv` — per-gene splice-junction counts and union-exonic length used for the combined splice+length normalization

### Variant Calling

`CapGTA_joint_variant_calling.sh` performs true joint variant calling with `bcftools mpileup -b bam_list | call -mv`, parallelized by fixed-size regions and concatenated into a single VCF.

```bash
bash CapGTA_joint_variant_calling.sh \
    results/<sample>/real_cells.csv \
    /shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa \
    results/<sample>/joint_variants \
    1000000
```

**Required arguments:**
- `<real_cells.csv>` cells to include (barcode,bam_path)
- `<reference.fa>` indexed reference FASTA (`.fai` required alongside)
- `<output_dir>` all outputs land here
- `<region_size>` bp per region for parallelization (1000000 = 1 Mb → ~104 tasks on WBcel235)

**Key outputs** (under `<output_dir>/`):
- `joint_variants.vcf.gz` (+ `.csi`) — single joint VCF across all cells; `FORMAT/AD, FORMAT/DP` retained for downstream cellspec analysis
- `per_region/region_*.vcf.gz` — per-region intermediates
- `bam_list.txt`, `regions.txt` — inputs echoed for reproducibility


# Phylogenetic and mutation spectrum analysis


# Bulk WGS preprocessing


# Public data preprocessing


# Analysis Notebooks

