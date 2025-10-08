# Preprocessing Pipeline Test Setup

This document describes the setup for testing the preprocessing pipeline with the original NovaSeq data.

## Data Setup

### Raw Input Data (Symlinked)

Located in `data/raw/`:
- **Read 1**: `k562_tree_1_S1_R1_001.fastq.gz` (1.8 TB)
- **Read 2**: `k562_tree_1_S1_R2_001.fastq.gz` (1.8 TB)

**Source**: `/fh/fast/srivatsan_s/SR/ngs/illumina/sanjay/20250228_LH00740_0051_B22T53HLT4/Unaligned/SS_TruSeq_FS_Ultra/`

These are symlinks to the original NovaSeq sequencing run data.

### Processed Data (To Be Generated)

The pipeline will generate:
- `sc_outputs/` - Per-cell BAM and VCF files
- `vcf_chunks/` - VCF chunks for processing
- `*.vcf.gz` - Joint-called VCF
- `*.h5ad` - AnnData objects with variant calls
- `real_cells.txt` - List of detected real cell barcodes

## Pipeline Overview

The preprocessing pipeline (`WGS_PP.sh`) performs the following steps:

1. **FASTQ Chunking**: Split large FASTQ files into chunks for parallel processing
2. **Parallel Preprocessing** (`PP_array.sh`): For each chunk:
   - Barcode error correction
   - Adapter trimming (Trim Galore)
   - BWA-MEM alignment
   - SAM to BAM conversion
3. **Cell Detection**: Knee plot analysis to identify real cells vs. empty droplets
4. **Single Cell Extraction** (`sc_from_chunks.sh`): Extract reads per cell barcode
5. **Variant Calling** (`sc_var_array.sh`): GATK HaplotypeCaller per cell
6. **Joint Calling** (`joint_calling_array.sh`): Combine variants across cells
7. **AnnData Generation**: Create scanpy-compatible AnnData object

## Running the Pipeline

### Method 1: Full Pipeline (Recommended for First Test)

This will process all ~3.6 TB of FASTQ data:

```bash
# Estimate read count (or use value from sequencing report)
# For this dataset, approximate read count can be estimated from file size
READ_COUNT=9000000000  # Approximate, adjust as needed

# Run the full pipeline
# Note: Temp directory defaults to /hpc/temp/srivatsan_s/SPC_genome_preprocessing
sbatch WGS_PP.sh \
  -o sc_PolE_test \
  -1 data/raw/k562_tree_1_S1_R1_001.fastq.gz \
  -2 data/raw/k562_tree_1_S1_R2_001.fastq.gz \
  -r $READ_COUNT \
  -g /shared/biodata/reference/GATK/hg38/ \
  -s ./scr \
  -O . \
  -n 500
```

**Parameters:**
- `-o sc_PolE_test` - Output name prefix
- `-1` / `-2` - Read 1 and Read 2 FASTQ files
- `-r` - Total read count (estimate from file size or sequencing report)
- `-g` - Reference genome directory (hg38)
- `-s ./scr` - Scripts directory
- `-O .` - Output directory (current directory)
- `-n 500` - Number of chunks for parallelization

### Method 2: Subset Test (Faster)

For a quick test, you can create a small subset of reads:

```bash
# Create a test subset (e.g., first 10 million reads)
mkdir -p data/test
zcat data/raw/k562_tree_1_S1_R1_001.fastq.gz | head -n 40000000 | gzip > data/test/test_R1.fastq.gz
zcat data/raw/k562_tree_1_S1_R2_001.fastq.gz | head -n 40000000 | gzip > data/test/test_R2.fastq.gz

# Run pipeline on subset
sbatch WGS_PP.sh \
  -o sc_PolE_subset_test \
  -1 data/test/test_R1.fastq.gz \
  -2 data/test/test_R2.fastq.gz \
  -r 10000000 \
  -g /shared/biodata/reference/GATK/hg38/ \
  -s ./scr \
  -O ./test_output \
  -n 50
```

## Resource Requirements

### Computational Resources

**For full pipeline:**
- **Time**: 2-4 days
- **CPUs**: 2-36 per job (varies by step)
- **Memory**: 16-64 GB per job
- **Storage**: ~125 GB for outputs
- **Temp Storage**: Significant space needed for FASTQ chunks

**SLURM partitions:**
- Main job: `campus` partition
- Array jobs: Configure in individual scripts

### Storage Breakdown

Expected disk usage:
- Raw FASTQs: 3.6 TB (symlinked, no additional space)
- FASTQ chunks (temp): ~4 TB (in `/hpc/temp/srivatsan_s/SPC_genome_preprocessing`)
- Aligned BAMs: ~30-60 GB
- VCF files: ~5-10 GB
- Single-cell outputs: ~50 GB
- AnnData: ~1-60 GB (depending on filtering)

**Total temporary**: ~4 TB (stored in `/hpc/temp/srivatsan_s/`)
**Total permanent**: ~125 GB

**Note**: The temp directory `/hpc/temp/srivatsan_s/SPC_genome_preprocessing` is now the default and has been pre-created for this project.

## Pipeline Scripts

Located in `scr/`:

1. **`PP_array.sh`** - Preprocessing array job
   - Barcode correction, trimming, alignment
   - Runs as SLURM array with `--array=1-N_CHUNKS`

2. **`sc_from_chunks.sh`** - Single cell extraction
   - Extracts reads by barcode from chunks
   - Creates per-cell BAM files

3. **`sc_var_array.sh`** - Variant calling array
   - GATK HaplotypeCaller per cell
   - Generates GVCF files

4. **`joint_calling_array.sh`** - Joint variant calling
   - Combines GVCFs across cells
   - Creates final multi-sample VCF

5. **`make_variant_anndata.py`** - AnnData creation
   - Converts VCF to AnnData format
   - Adds cell and variant annotations

## Barcodes

Located in `barcodes/`:
- `bcA_24_exD.txt`
- `bcB_24_exC.txt`
- `bcC_24_exB.txt`
- `bcD_24_exA.txt`

These contain the barcode sequences for the 4-level combinatorial indexing used in the experiment.

## Expected Outputs

### Intermediate Files

- `chunk_indices.txt` - List of chunk identifiers
- `SLURM_outs/` - SLURM job output logs
- `readcounts.csv` - Read counts per barcode
- `real_cells.txt` - Barcodes passing cell detection threshold
- `knee_plot.png` - Knee plot for cell detection

### Final Outputs

- `sc_PolE_test.vcf.gz` - Joint-called VCF
- `sc_PolE_test.vcf.gz.csi` - VCF index
- `sc_PolE_all_SNPs.h5ad` - All detected variants
- `sc_PolE_filtered.h5ad` - Quality-filtered variants
- `sc_outputs/` - Per-cell BAM and GVCF files
- `vcf_chunks/` - Intermediate VCF chunks

## Monitoring Progress

### Check SLURM Jobs

```bash
# View all your jobs
squeue -u $USER

# Check specific job
squeue -j <job_id>

# View job details
scontrol show job <job_id>
```

### Check Output Logs

```bash
# View preprocessing log
tail -f SLURM_outs/Preprocessing_<job_id>.out

# View array job logs
ls SLURM_outs/array_outs/
tail SLURM_outs/array_outs/PP_array_<job_id>_<task_id>.out
```

### Track File Generation

```bash
# Count chunks processed
ls sc_outputs/*.bam 2>/dev/null | wc -l

# Check VCF generation
ls sc_outputs/*.g.vcf.gz 2>/dev/null | wc -l

# Monitor disk usage
du -sh sc_outputs/ vcf_chunks/
```

## Troubleshooting

### Common Issues

**1. Out of disk space in temp directory**
```bash
# The default temp directory is now /hpc/temp/srivatsan_s/SPC_genome_preprocessing
# If you need a different location, use the -t flag:
sbatch WGS_PP.sh ... -t /path/to/alternative/temp/dir

# To clean the default temp directory:
rm -rf /hpc/temp/srivatsan_s/SPC_genome_preprocessing/*
```

**2. Job killed due to memory**
```bash
# Edit the array script to request more memory
# In scr/PP_array.sh, increase:
#SBATCH --mem=32G  # or higher
```

**3. Reference genome not found**
```bash
# Verify reference genome exists and is indexed
ls -lh /shared/biodata/reference/GATK/hg38/

# Should contain:
# - genome.fa (or similar)
# - genome.fa.fai
# - BWA index files (*.amb, *.ann, *.bwt, *.pac, *.sa)
```

**4. Module loading errors**
```bash
# If modules fail to load, check available versions
module spider BWA
module spider SAMtools
module spider GATK

# Update module load commands in scripts if needed
```

## Comparison with Original Run

To verify your preprocessing matches the original:

```bash
# Compare cell counts
wc -l real_cells.txt
wc -l /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq/real_cells.txt

# Compare VCF sizes (should be similar)
ls -lh *.vcf.gz
ls -lh /fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq/*.vcf.gz

# Compare AnnData shapes
python -c "import scanpy as sc; \
  adata = sc.read_h5ad('sc_PolE_filtered.h5ad'); \
  print('New:', adata.shape)"

python -c "import scanpy as sc; \
  adata = sc.read_h5ad('/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/sc_PolE_novaseq/sc_PolE_filtered.h5ad'); \
  print('Original:', adata.shape)"
```

## Next Steps

After successful preprocessing:

1. **Quality Control**: Run `sc_analysis.ipynb` to perform QC
2. **Variant Filtering**: Apply quality filters to variants
3. **Mutational Spectra**: Analyze mutation patterns
4. **Phylogenetic Analysis**: Build cell lineage trees
5. **Comparison**: Compare with original analysis results

## Notes

- The original preprocessing used read count from the sequencing report
- Chunk size is automatically calculated based on total reads and number of chunks
- The pipeline uses SLURM job dependencies to ensure correct execution order
- Temporary files in `$TMP_DIR` are cleaned up automatically (unless job fails)
- Consider archiving intermediate files to S3 for long-term storage

## Getting Read Count

If you need the exact read count from the FASTQ files:

```bash
# WARNING: This will take a long time for 1.8TB files!
READ_COUNT=$(zcat data/raw/k562_tree_1_S1_R1_001.fastq.gz | wc -l | awk '{print $1/4}')
echo "Total reads: $READ_COUNT"
```

Or estimate from file size:
```bash
# Typical FASTQ compression: ~4 bytes per base
# Estimate: file_size_bytes / (read_length * 2)  # *2 for forward and reverse
# For 150bp reads: 1.8TB / (150 * 2) ≈ 6 billion reads
```
