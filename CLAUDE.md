# SPC genome preprocessing

### Project Overview

This repository contains the code for processing sequencing data associated with the preprint [Capsule Based Genome Sequencing and Lineage Tracing](https://www.biorxiv.org/content/10.1101/2025.03.14.643253v1.full.pdf). This paper has been updated, so some of what I outline below may not pertain to the contents of the paper, but new data that we have added. In this paper we describe the use of semi-permeable capsules to perform single cell whole genome amplification at scale, and prepare illumina sequencing libraries. We use single cell combinatorial indexing, wherein encapsulated cells are sequentially split and short oligos are added to their gDNA such that each cell recieves a random combination of 4 oligos to create a unique barcode. We term this new assay CapWGS. In the paper, we benchmark the performance of CapWGS, and then apply it to perform lineage tracing on cells harboring hypermutator alleles of DNA polymerase epsilon. The referenced preprint has been reviewed, and we are in the process of adding a substantial amount of new data to address reviewer comments. We have performed more rigorous benchmarking using deep sequencing of HSCs, and we've also added a genome-transcriptome coassay - CapGTA - which we have applied to perform lineage tracing on cells isolated from c elegans.  

There are 3 main pipeline scripts:

* `CapWGS_PP.sh`: the main script for Capsule based whole genome sequencing- its primary purpose is to perform joint variant calling with GATK. I typically use the grch38 for precise variant calling with many alt contigs. 
* `CapWGS_PP_QC_only.sh`: This pipeline is for when we are benchmarking coverage uniformity and basic sequencing run stats, and not necessarily performing variant calling. It outputs single cell bams, bigwigs, and some per-cell qc metrics. I typically use the standard grch37 (fewer haplotypes/alt contigs) to perform coverage uniformity benchmarking. 
* `CapGTA_PP.sh`: For joint genome and transcriptome co-assay data. This uses star for alignment, and separates reads containing splice junctions into a separate RNA bam, from which it generates count matrices. 

Additionally, we often need to process bulk WGS data for benchmarking or other experiments, so there are bulk processing scripts in `scripts/bulk/`.

Each of these pipelines are built for the fred hutch cluster, which uses `SLURM`. You can find context about this cluster and `SLURM` in `../scicomp_context/`.

Each of these scripts take raw fastq's as input, with arguments for things like the reference genome, and information about how to chunk the files for parallelization. 

These main scripts mostly serve as job schedulers, and call scripts from `scripts/{pipeline}/`. They all chunk their input files and use SLURM job arrays. 

`scripts/utils/` contains scripts that are used by multiple pipelines, like `atrandi_demux.py`, which parses the barcodes from the fastqs. 

As a convention output alignments to a sample specific folder in `data/`, and output results (csv's, figures, etc) to corresponding sample specific folder in `results/`.

In data, the bulked alignments should be written to `data/{sample}/`, and all files corresponding to single cells should be stored in a subdirectory: `data/{sample}/sc_outputs/`.

You are working in a micromamba environment containing any python libraries needed for any of the scripts used by the pipeline, or that you might need for basic scientific computing, so you do not need to load the fred hutch python module. 

For basic bioinformatics tools like `samtools`, `bcftools`, `GATK`, `picard`, `deeptools`, you can load them from modules. Job error and output logs are always written to `SLURM_outs` and `SLURM_outs/array_outs`. 

`bin` is a place to store miscellaeneous scripts that are not required for the main pipelines, such as job submission scripts, and one off analysis scripts. 

I perform analysis of this data using a python package I am developing called [`cellspec`](https://github.com/harrispopgen/cellspec), which is installed in this environment. I do all of this analysis in jupyter notebooks found in `notebooks`, and save my results to `results`. 

**Commonly used reference genomes**:

When processing human data for variant calling, I typically use grch38:
`/shared/biodata/reference/GATK/hg38/Homo_sapiens_assembly38.fasta`
BWA index located here:
`/shared/biodata/reference/GATK/hg38/BWAindex`

When processing human data for evaluating coverage uniformity, I typically use grch37:
`/shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa`

For c elegans data, I use a bristol N2 strain specific reference I have downloaded and built indices for here:
`data/reference/worm_GCA_028201515.1_combined` (`STARIndex/` for CapGTA, `BWAIndex/` for bulk WGS)