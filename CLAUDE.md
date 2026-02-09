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

### Current tasks

**quick fixes:**

**Bigger jobs**

1. In `scripts/bulk/` directory, you'll notice there are scripts to perform joint calling with bcftools (less rigorous) or following the GATK best practices. I would like to implement the ability to perform variant calling with either of these tools in `CapWGS_PP.sh` and `CapGTA_PP.sh`. I often use bcftools for shallower pilot runs where I can tolerate imprecise variant calling, but prefer GATK for deeper runs. It is most important that `CapWGS_PP.sh` implements the GATK best practices, including marking duplicates. You can read about GATK best practices for data preprocessing [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery) and for germline SNP calling [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels) (I follow the germline best practices becuase even though these are technically somatic variants we are calling, becuase it is single cell data it resembles germline discovery more)


**Data to process**
1. I have bulk WGS fastqs that need to be processed using `scripts/bulk/gatk_pipeline.sh`. This a 4 billion read version of K562_mut_accumulation_pilot. Fastq's are located in `/fh/fast/srivatsan_s/SR/ngs/illumina/sanjay/20260130_LH00740_0176_A23F3HLLT4/Unaligned/MutationAccumulationHuman`. Samples are named `AAVS_Clone_[4-6]_P[1-3]` or `PolE_Clone_[4-6]_P[1-3]`. PolE_Clone_5 does not have a P3, so there are 34 fastq's corresponding to 17 total samples. Run `scripts/bulk/gatk_pipeline.sh` and put individual bam files in `data/K562_mut_accumulation/bams/` and individual .g.vcf files in `data/K562_mut_accumulation/gvcfs/`. Output the final joint calling vcf to `data/K562_mut_accumulation/`. Use grch38 for alignment and variant calling `/shared/biodata/reference/GATK/hg38/Homo_sapiens_assembly38.fasta`

2. There are several sequencing runs of HSC's, both CapWGS and bulk WGS, that need to be processed. This is for the benchmarking portion of th paper. We are using the same datasets to benchmark coverage and variant calling. The data is processed with `CapWGS_PP_QC_only.sh` and the UCSC grch37 reference (see above) for coverage benchmarking, and `CapWGS_PP.sh ` using our new GATK best practices calling mode and the GATK grch38 reference (see above). Bulk samples are included for comparison, and should be processed with the bulk pipeline (for benchmarking coverage, use `scripts/bulk/bcftools_align_single.sh`, for variant benchmakring use the GATK pipeline we adjusted in `scripts/bulk/`). Below is a list of the samples, and checkmarks indicate whether or not they've been processed for both benchmarking_coverage and benchmarking_variant. For each, I have provided a directory where you can find the raw fastq's. All of these have a corresponding folder in `data/` and `results/`

    * benchmarking_coverage
        - [ ] HSC2_bulk : `/fh/fast/srivatsan_s/SR/ngs/illumina/sanjay/20260130_LH00740_0176_A23F3HLLT4/Unaligned/HSC2_bulk_single_cell/HSC2_bulk_rep*` (2 replicates to process separately, multiple lanes each)
        - [ ] HSC2_enzyme : `/fh/fast/srivatsan_s/SR/ngs/illumina/sanjay/20260130_LH00740_0176_A23F3HLLT4/Unaligned/HSC2_bulk_single_cell/HSC2_bigSPC_mMDA_PTA_*` (multiple lanes)
        - [x] HSC_bulk
        - [x] HSC_CellGrowth
        - [x] HSC_enzyme
        - [x] public
    * benchmarking_variant
        - [ ] HSC2_bulk : `/fh/fast/srivatsan_s/SR/ngs/illumina/sanjay/20260130_LH00740_0176_A23F3HLLT4/Unaligned/HSC2_bulk_single_cell/HSC2_bulk_rep*` (2 replicates to process separately, multiple lanes each)
        - [ ] HSC2_enzyme : `/fh/fast/srivatsan_s/SR/ngs/illumina/sanjay/20260130_LH00740_0176_A23F3HLLT4/Unaligned/HSC2_bulk_single_cell/HSC2_bigSPC_mMDA_PTA_*` (multiple lanes)
        - [ ] HSC_bulk : `/fh/working/srivatsan_s/CapWGS/HSC_bulk_single_cell/HSC_bulk*` (multiple lanes)
        - [ ] HSC_CellGrowth : `/fh/working/srivatsan_s/CapWGS/HSC_bulk_single_cell/spcHSC_CellGrowth*` (multiple lanes)
        - [ ] HSC_enzyme : `/fh/working/srivatsan_s/CapWGS/HSC_bulk_single_cell/spcHSC_enzyme*` (multiple lanes)

Read counts for samples in `/fh/fast/srivatsan_s/SR/ngs/illumina/sanjay/20260130_LH00740_0176_A23F3HLLT4/Unaligned/HSC2_bulk_single_cell/`:
    Sample	Reads
    HSC2_bigSPC_mMDA_PTA	17202550439
    HSC2_bulk_rep1	66530138
    HSC2_bulk_rep2	10083006

Read counts for samples in `/fh/working/srivatsan_s/CapWGS/HSC_bulk_single_cell/`:
    Sample	Reads
    spcHSC_CellGrowth	10654389129
    spcHSC_enzyme	21433544624
    HSC_bulk	1283653389
---

## Recent Session Summary (Feb 6-9, 2026)

### Completed Fixes and Implementations

**1. Variant Caller Selection (Issue #2)**
- Implemented `-v/--variant-caller` parameter in `CapWGS_PP.sh` to select between bcftools (default) and gatk
- GATK mode follows best practices: mark duplicates (Picard) + BQSR + HaplotypeCaller in GVCF mode
- Created supporting scripts:
  - `scripts/CapWGS/markdup_bqsr.sh`: Preprocessing bulk BAM with mark duplicates and BQSR
  - `scripts/CapWGS/sc_from_bam.sh` + `extract_sc_from_bam_array.sh`: Extract single cells from preprocessed BAM
  - `scripts/CapWGS/sc_var_array_gatk.sh`: GATK HaplotypeCaller for single cells
  - `scripts/CapWGS/sc_var_array_bcftools.sh`: BCFtools variant calling for single cells
  - `scripts/CapWGS/bcftools_joint_calling.sh`: Merge BCFtools VCFs
  - `scripts/utils/detect_known_sites.sh`: Auto-detect dbSNP and known indels in bundle/ subdirectory

**2. Read Group Addition**
- Fixed GATK BQSR error: "Number of read groups must be >= 1, but is 0"
- Added read groups to BWA alignments with `-R` flag in:
  - `scripts/CapWGS/PP_array.sh` (passes SAMPLE_NAME parameter)
  - `scripts/bulk/gatk_single_sample.sh` (new single-sample GATK pipeline)
  - `bin/HSC2_bulk_variant.sh` (standalone HSC2 bulk variant script)
- Each read group includes: ID, SM (sample), PL (ILLUMINA), LB (library)

**3. Script Path Fixes for Directory Reorganization**
- Fixed 7 script path references to match CapWGS/ and CapWGS_QC/ directory structure:
  - `sc_from_chunks.sh`: extract_sc_array.sh → utils/extract_sc_array.sh
  - `submit_bigwig_lorenz.sh`: generate_bigwig_and_lorenz_array.sh → CapWGS_QC/
  - `submit_joint_calling.sh`: joint_calling_array.sh → CapWGS/ (2 occurrences)
  - `submit_benchmarking_qc.sh`: benchmarking_qc_array.sh → CapWGS_QC/
  - `generate_bigwig_and_lorenz_array.sh`: lorenz.py → CapWGS_QC/
  - `coverage_tracks.sh`: coverage_tracks.py → CapWGS_QC/
  - `binned_coverage_array.sh`: binned_coverage.py → CapWGS_QC/

**4. Pipeline Path Handling**
- Fixed `CapWGS_PP.sh` and `CapWGS_PP_QC_only.sh` to handle nested output paths
- Added `SAMPLE_NAME=$(basename "${OUTPUT_NAME}")` to support directory structures like `benchmarking_coverage/HSC2_enzyme`

**5. Bulk GATK Pipeline Improvements**
- Updated `scripts/bulk/gatk_align_call_array.sh` to use detect_known_sites.sh utility
- Fixed sample list generation to preserve individual passages (not merge them)
- Increased walltime from 24 hours to 3 days for large datasets
- Added SCRIPTS_ROOT parameter for proper path resolution in SLURM jobs
- Created `scripts/bulk/gatk_single_sample.sh` for single-sample GATK processing (produces final VCF, not GVCF)

**6. PP_array.sh BWA Index Fix**
- Modified `scripts/CapWGS/PP_array.sh` to handle both directory and FASTA file paths for genome parameter
- Automatically detects if genome is a directory or file and constructs BWA index path accordingly

### Current Job Status (as of Feb 9, 2026)

**K562 Mutation Accumulation (Bulk GATK):**
- ✅ Jobs completed successfully (46110750, 46110773)
- 17 samples processed with GATK best practices
- Location: `data/K562_mut_accumulation/`

**HSC2 Benchmarking Data:**

1. **HSC2_bulk coverage** (bcftools)
   - ✅ COMPLETED (job 46114040, 2h 4m)
   - Output: `data/benchmarking_coverage/HSC2_bulk/HSC2_bulk.bam`

2. **HSC2_bulk variant** (GATK)
   - ✅ COMPLETED (job 46157680, 14h 15m)
   - Output: `data/benchmarking_variant/HSC2_bulk/gvcfs/HSC2_bulk.g.vcf.gz` (4.4GB)
   - Note: GVCF created; need to run GenotypeGVCFs for final VCF

3. **HSC2_enzyme coverage** (CapWGS QC-only, 17B reads)
   - 🔄 IN PROGRESS - Resubmitted with fixes
   - Bulk BAM completed: `data/benchmarking_coverage/HSC2_enzyme/HSC2_enzyme.bam` (408GB)
   - Currently extracting top 80 single cells (job 46221307)
   - Pending: Bigwig and Lorenz curve generation

4. **HSC2_enzyme variant** (CapWGS + GATK, 17B reads)
   - ❌ FAILED - All 501 preprocessing tasks failed
   - Issue: BWA index path problem when passing FASTA file instead of directory
   - PP_array.sh fix has been implemented
   - **TODO**: Clean up and resubmit with corrected paths

**Remaining Benchmarking Data to Process:**
- HSC_bulk (variant only)
- HSC_CellGrowth (variant only)
- HSC_enzyme (variant only)

### Known Issues and Next Steps

1. **HSC2_enzyme variant pipeline** needs resubmission after PP_array.sh BWA index fix
2. **HSC2_bulk GVCF** needs GenotypeGVCFs step to produce final VCF (can use new `gatk_single_sample.sh` approach in future)
3. **HSC2_enzyme coverage** awaiting single cell extraction completion, then submit bigwig/Lorenz generation
4. Process remaining HSC samples for variant benchmarking using `CapWGS_PP.sh -v gatk`

### Git Commits on Branch `2-vcaller_selection`

- `6e9da20` Fix script paths for reorganized directory structure
- `2a5f62d` Implement variant caller selection for CapWGS pipeline (#2)
- `9783ccb` Add read groups to BWA alignments for GATK compatibility
- Earlier commits: Remove AnnData generation, fix BAM paths, etc.
