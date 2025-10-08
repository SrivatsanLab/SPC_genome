# Data Organization and Requirements

This document describes the expected data structure and requirements for the SPC Genome project.

## Data Location

The main data directory should be located at:
```
../sc_PolE_novaseq/
```

This is relative to the project root directory.

## Expected Data Structure

```
sc_PolE_novaseq/
├── raw_data/
│   ├── sample_R1.fastq.gz    # Read 1 sequencing files
│   ├── sample_R2.fastq.gz    # Read 2 sequencing files
│   └── ...
├── processed/
│   ├── bam_files/
│   ├── vcf_files/
│   └── anndata/
└── metadata/
    └── sample_info.csv        # Sample metadata (if applicable)
```

## Input Requirements

### Raw Sequencing Data

**Format:** Paired-end FASTQ files (gzipped)

**Files needed:**
- `*_R1.fastq.gz` or `*_read1.fastq.gz` - Read 1
- `*_R2.fastq.gz` or `*_read2.fastq.gz` - Read 2

**Quality metrics:**
- Minimum read length: 50 bp
- Recommended depth: 0.1-1x coverage per cell for CapWGS

### Reference Genome

**Location:** `/shared/biodata/reference/GATK/hg38/` (default on Fred Hutch cluster)

**Required files:**
- `*.fa` or `*.fasta` - Reference genome FASTA
- `*.fa.fai` - FASTA index (created by samtools faidx)
- BWA index files (`*.amb`, `*.ann`, `*.bwt`, `*.pac`, `*.sa`)

**Note:** If using a different reference, update `config.yaml`

### Metadata (Optional)

Sample information file with:
- Sample names
- Barcodes
- Cell counts
- Sequencing run information
- Read counts per sample

## Storage Considerations

### Fred Hutch Storage Options

1. **`/fh/fast/`** - High-performance storage (recommended for active analysis)
   - 5TB free, then charged
   - Best for raw data and intermediate files

2. **`/home/`** - User home directory
   - 100GB limit
   - For code and small files only

3. **`/hpc/temp`** - Temporary storage
   - For intermediate processing files
   - Files may be deleted after 30 days

4. **`s3://fh-pi-*/`** - Economy cloud storage
   - 100TB free
   - For long-term archival

### Space Requirements

Estimated space needed per sample:
- Raw FASTQ (compressed): ~20-50 GB
- Aligned BAM files: ~30-60 GB
- VCF files: ~5-10 GB
- AnnData objects: ~1-5 GB

**Total:** ~60-125 GB per sample

## Data Preparation Checklist

Before running the pipeline:

- [ ] Raw FASTQ files are located in accessible storage
- [ ] Files are gzipped (`.fastq.gz`)
- [ ] Read 1 and Read 2 files are paired correctly
- [ ] Know the total read count (from sequencing report)
- [ ] Reference genome is indexed and accessible
- [ ] Sufficient storage space is available
- [ ] `config.yaml` is updated with correct paths

## Getting Read Count

Read count can be obtained from:

1. **Sequencing run report** (recommended)
2. **Count manually** (slow):
   ```bash
   zcat sample_R1.fastq.gz | wc -l | awk '{print $1/4}'
   ```

## Troubleshooting

### Common Issues

**Issue:** "No such file or directory"
- Check that paths in `config.yaml` are correct
- Use absolute paths if relative paths fail

**Issue:** "Permission denied"
- Ensure you have read access to data directory
- Check file permissions: `ls -l data_file`

**Issue:** "Disk quota exceeded"
- Check storage usage: `du -sh /fh/fast/pi_name/`
- Move old files to S3 or delete temporary files

**Issue:** "Reference genome not found"
- Verify reference path exists
- Check that BWA index files are present
- Consider using default: `/shared/biodata/reference/GATK/hg38/`
