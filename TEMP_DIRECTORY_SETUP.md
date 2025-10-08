# Temp Directory Configuration

## Location

**Default temp directory**: `/hpc/temp/srivatsan_s/SPC_genome_preprocessing`

This directory has been created and is ready for use by the preprocessing pipeline.

## Purpose

During preprocessing, the pipeline will:
1. Split 3.6 TB of FASTQ files into 500 chunks
2. Store chunks temporarily (~4 TB total)
3. Process each chunk in parallel
4. Clean up chunks as processing completes

## Permissions

```bash
$ ls -la /hpc/temp/srivatsan_s/SPC_genome_preprocessing
drwxrws--- 2 dmullane srivatsan_s_grp 4096 Oct  8 15:31 .
```

- Owner: dmullane
- Group: srivatsan_s_grp
- Permissions: rwxrws--- (group writable with setgid)

## Configuration

The temp directory is now the default in:

1. **WGS_PP.sh** (Line 12):
   ```bash
   TMP_DIR="/hpc/temp/srivatsan_s/SPC_genome_preprocessing"
   ```

2. **config.yaml** (Line 29):
   ```yaml
   tmp_dir: "/hpc/temp/srivatsan_s/SPC_genome_preprocessing"
   ```

3. **run_preprocessing_test.sh** (Line 18):
   ```bash
   TMP_DIR="/hpc/temp/srivatsan_s/SPC_genome_preprocessing"
   ```

## Monitoring Disk Usage

```bash
# Check current usage
du -sh /hpc/temp/srivatsan_s/SPC_genome_preprocessing

# Monitor in real-time during run
watch -n 60 'du -sh /hpc/temp/srivatsan_s/SPC_genome_preprocessing'

# Check available space on /hpc/temp
df -h /hpc/temp
```

## Cleaning Up

After successful pipeline completion:

```bash
# Remove all temporary files
rm -rf /hpc/temp/srivatsan_s/SPC_genome_preprocessing/*

# Verify cleanup
ls -la /hpc/temp/srivatsan_s/SPC_genome_preprocessing
```

## Override Default

If you need to use a different temp directory:

**Command line:**
```bash
sbatch WGS_PP.sh ... -t /path/to/alternative/temp/dir
```

**Script:**
```bash
./run_preprocessing_test.sh -t /path/to/alternative/temp/dir
```

## Why /hpc/temp?

- **Large capacity**: Suitable for multi-TB temporary storage
- **Shared access**: Group members can access if needed
- **Cluster-optimized**: Better performance than `/loc/scratch` for large files
- **Persistent**: Files remain until explicitly deleted (unlike `/loc/scratch` which cleans up after job completion)

## Troubleshooting

### "No space left on device"

Check available space:
```bash
df -h /hpc/temp
```

If full, either:
1. Clean up old files in other directories
2. Use alternative temp directory with `-t` flag

### Permission denied

Ensure you're in the srivatsan_s_grp group:
```bash
groups | grep srivatsan_s_grp
```

### Stale files from previous runs

Clean up manually:
```bash
rm -rf /hpc/temp/srivatsan_s/SPC_genome_preprocessing/*
```

## Notes

- The directory persists between runs
- Clean up manually after each run to free space
- Consider setting up a cron job for automatic cleanup of old files
- Files are NOT automatically deleted after 30 days like in some scratch spaces
