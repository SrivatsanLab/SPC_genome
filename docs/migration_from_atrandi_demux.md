# Migration: `atrandi_demux.py` → `spc-demux` (Rust)

Guide for switching `SPC_genome/CapWGS_PP.sh` and `SPC_genome/CapGTA_PP.sh`
from the legacy Python demultiplexer to the new Rust binary.

---

## What changes

| Aspect | `atrandi_demux.py` | `spc-demux` |
|--------|--------------------|-------------|
| Input | One R1/R2 chunk per array task | Full R1/R2 files, one demux pass per sample |
| Output (FASTQ mode) | Paired corrected R1 + R2 FASTQ | 576 bucketed FASTQ files on local scratch, each interleaved R1 (R2 only when `>15bp` remain after barcode strip) |
| Output (SAM mode) | n/a | Single unaligned SAM/BAM with `CB:Z`/`BC:Z`/`RG:Z` tags |
| `CB:Z:` length | 45 char (4 × 8 bp barcode + 3 × 4 bp linker + 1 bp tail) | **32 char** (4 × 8 bp barcode only, linkers dropped) |
| Trimming | External (`Trim Galore`) | Inline R1 adapter + R2-tail-RC trimming |
| Parallelism | Per-chunk SLURM array (the orchestrator splits FASTQs) | Per-bucket (576 array tasks) or single-task STAR |
| Throughput | ~50–100 k reads/s | ~1–2 M reads/s |

The 45 → 32 char `CB:Z:` change is the breakage hotspot. Anything downstream
that regex-matches a 45-base CB, or that slices the CB to extract individual
round barcodes, needs updating. The new layout is:

```
CB[ 0: 8] = bcD
CB[ 8:16] = bcC
CB[16:24] = bcB
CB[24:32] = bcA
```

Well-map lookups using `resources/barcodes/bc*_well_map.tsv` still work — they
map individual 8 bp barcodes to wells, not the concatenated CB.

---

## CapWGS_PP.sh (BWA-MEM)

### Target architecture

```
spc-demux ──► 576 bucketed FASTQs on local scratch
                │
                └──► SLURM array (1 task per bucket): BWA-MEM ──► per-bucket BAM
                                                                   │
                                                                   └──► samtools merge ──► final BAM
```

The fastq chunker, the per-chunk Python demux, and Trim Galore all collapse
into one `spc-demux` call. The downstream BWA alignment becomes a 576-way
array (one task per bucket) instead of the current `N_CHUNKS=500` array.

### Orchestrator changes (`CapWGS_PP.sh`)

**Replace** the chunking block at lines 240–263:

```bash
# OLD
zcat $READ1 | head -n $total_lines | split -d -l $CHUNK_LINES - "$TMP_DIR/read1_chunk_" &
zcat $READ2 | head -n $total_lines | split -d -l $CHUNK_LINES - "$TMP_DIR/read2_chunk_" &
wait
ls "$TMP_DIR"/read1_chunk_* | sed 's/.*chunk_//' > "${RESULTS_DIR}/chunk_indices.txt"
chunk_count=$(wc -l < "${RESULTS_DIR}/chunk_indices.txt")
```

with:

```bash
# NEW — single demux pass writes 576 buckets to scratch_dir from config.yaml
module load SAMtools  # only needed if you set --sam-output later; harmless to keep
SPC_DEMUX="${SCRIPTS_DIR}/scripts/utils/spc-demux"   # or wherever you install it
SPC_CONFIG="${SCRIPTS_DIR}/config/spc_demux.yaml"

${SPC_DEMUX} \
    --config "${SPC_CONFIG}" \
    --r1 ${READ1} --r2 ${READ2} \
    --sample "${SAMPLE_NAME}" \
    --stats "${RESULTS_DIR}/demux_stats.json" \
    --cell-counts "${RESULTS_DIR}/cell_counts.tsv"

# scratch_dir is set in the demux config.yaml. Enumerate non-empty buckets:
SCRATCH_DIR=$(grep -E '^\s*scratch_dir:' "${SPC_CONFIG}" | awk '{print $2}')
find "${SCRATCH_DIR}" -name "${SAMPLE_NAME}_*.fastq.gz" -size +28c \
    > "${RESULTS_DIR}/bucket_list.txt"   # 28 bytes = empty gzip member
bucket_count=$(wc -l < "${RESULTS_DIR}/bucket_list.txt")
```

**Replace** the SLURM array submission at line 268:

```bash
# OLD: PP_array.sh, called per-chunk
PP_array_ID=$(sbatch --parsable --array=1-$chunk_count "${SCRIPTS_DIR}/scripts/CapWGS/PP_array.sh" \
    "${RESULTS_DIR}/chunk_indices.txt" "${REFERENCE_GENOME}" "${SCRIPTS_DIR}" \
    "${TMP_DIR}" "${SAMPLE_NAME}" "${RESULTS_DIR}" "${DEBUG}")

# NEW: align_bucket_array.sh, called per-bucket (see below)
PP_array_ID=$(sbatch --parsable --array=1-$bucket_count "${SCRIPTS_DIR}/scripts/CapWGS/align_bucket_array.sh" \
    "${RESULTS_DIR}/bucket_list.txt" "${REFERENCE_GENOME}" "${SCRIPTS_DIR}" \
    "${TMP_DIR}" "${SAMPLE_NAME}" "${RESULTS_DIR}" "${DEBUG}")
```

**Drop the `-f 0x2` filter** at line 393. Bucketed FASTQ is R1-only for the
standard SPC read structure (R2 = 45 bp barcode region, nothing left after
strip), so no read will be flagged as properly paired:

```bash
# OLD
samtools view -@ 36 -f 0x2 -b "${UNFILTERED_BAM}" | samtools sort -@ 36 -o "${BAM_FILE}"

# NEW — keep mapped reads only
samtools view -@ 36 -F 4 -b "${UNFILTERED_BAM}" | samtools sort -@ 36 -o "${BAM_FILE}"
```

Update the read-count math at line 425 (now SE):

```bash
# OLD
total_input_reads=$((READ_COUNT * 2))   # paired-end

# NEW
total_input_reads=${READ_COUNT}         # SE — one BAM record per input read pair
```

You can also delete:
- `N_CHUNKS` (always 576 buckets now — comes from `bucket_rounds: [round_A, round_B]` in the demux config)
- `--use-existing-chunks` mode (demux is one pass over the whole input; not worth caching)
- `compile_preprocessing_stats.py` invocation (the demux emits one stats JSON; no per-chunk aggregation needed)

### Per-task alignment script (`align_bucket_array.sh`)

Replaces `PP_array.sh`. Single bucket → single SAM. No demux, no trimming.

```bash
#!/bin/bash
#SBATCH -J align_bucket
#SBATCH -o SLURM_outs/array_outs/%x_%A_%a.out
#SBATCH -c 4
#SBATCH -p short
#SBATCH -t 4:00:00

set -euo pipefail
mkdir -p SLURM_outs/array_outs

bucket_list="$1"; genome="$2"; scripts_DIR="$3"
TMP_DIR="$4"; SAMPLE_NAME="$5"; RESULTS_DIR="$6"; DEBUG="${7:-false}"

BUCKET_FQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${bucket_list}")
NAME=$(basename "${BUCKET_FQ}" .fastq.gz)

module load BWA SAMtools

# (Re-use the BWA_INDEX discovery block from the existing PP_array.sh — omitted here for brevity.)
BWA_INDEX=...

RG="@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:ILLUMINA\tLB:${SAMPLE_NAME}"

# Key flags:
#   -C   propagate FASTQ header comment (CB:Z:) into SAM tag
#   no -p — bucket FASTQ is single-end R1 (R2 dropped by demux for standard SPC structure)
bwa mem -t 4 -c 1 -C -R "${RG}" "${BWA_INDEX}" "${BUCKET_FQ}" \
    > "${TMP_DIR}/${NAME}.sam"
```

`bwa mem -C` is already used in the legacy `PP_array.sh` (line 190) for the
exact same reason — it just continues to work, since `spc-demux` annotates the
FASTQ header as `@... CB:Z:<32bp>` (same form as before, just shorter CB).

> **If you have true PE WGS data** (R1 and R2 both genomic, barcode somewhere
> else): `spc-demux` does not currently write paired interleaved bucketed
> FASTQ. You either need to add that feature, or keep using `atrandi_demux.py`
> for that read structure. The 1×250 SPC structure (R2 is barcode-only) is the
> supported case.

### Bucketed output is a partition of the cell space — exploit it

Bucketing isn't just a parallelism trick. A cell's 32 bp CB decomposes as
`bcD(8) | bcC(8) | bcB(8) | bcA(8)`, and its bucket file is
`{sample}_{bcA}_{bcB}.bam`. **Every read for a given cell lives in exactly
one bucket.** This changes two stages of the current CapWGS pipeline.

**1. Cell detection — skip the BAM readcount round-trip**

Today, `CapWGS_PP.sh` scans the merged BAM to compute per-barcode read counts
(lines 408–409):

```bash
# OLD
samtools view -@ 36 "${BAM_FILE}" | \
    python "${SCRIPTS_DIR}/scripts/utils/readcounts.py" -o "${RESULTS_DIR}/readcounts.csv"
```

`spc-demux --cell-counts` already emits this exact data, computed inline
during demux (pre-alignment) and pre-sorted by read count. Skip the BAM
round-trip:

```bash
# NEW
python "${SCRIPTS_DIR}/scripts/utils/detect_cells.py" \
    --plot "${RESULTS_DIR}/kneeplot.png" \
    < "${RESULTS_DIR}/cell_counts.tsv" \
    > "${RESULTS_DIR}/real_cells.txt"
```

> Caveat: the spc-demux counts are pre-alignment, so they include reads that
> won't ultimately align. For most samples the kneeplot shape is identical
> (alignment rate is roughly uniform across cells), but for samples with
> per-cell alignment-rate skew you may want to keep the BAM-derived readcounts
> as the cell-detection input. A/B-check on the first run.

**2. Per-cell BAM extraction — go to the right bucket, don't scan the bulk**

`sc_extract_preprocess_qc_array.sh` line 113 currently extracts each cell by
grep-streaming the entire bulk BAM:

```bash
# OLD — one full BAM scan per cell
samtools view -h -F 0x900 "${BULK_BAM}" | \
    grep -E "^@|CB:Z:${BARCODE}" | \
    samtools sort -@ 4 -o "${RAW_BAM}"
```

For *N* cells this is *N* full BAM scans. With bucketing, the cell's reads
are confined to one bucket BAM (up to 576 cells per bucket, vs. all cells in
the bulk). Derive the bucket from the CB:

```bash
# NEW — extract from the one bucket that contains this cell
BARCODE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${BARCODE_FILE}")
BCA="${BARCODE:24:8}"   # CB[24:32]
BCB="${BARCODE:16:8}"   # CB[16:24]
BUCKET_BAM="${BUCKET_BAM_DIR}/${SAMPLE_NAME}_${BCA}_${BCB}.bam"

samtools view -h -F 0x900 "${BUCKET_BAM}" | \
    grep -E "^@|CB:Z:${BARCODE}" | \
    samtools sort -@ 4 -o "${RAW_BAM}"
```

The bucket BAM is ~1/576 the size of the bulk (for an evenly-loaded sample),
so per-cell extraction speeds up by roughly that factor. Pass
`BUCKET_BAM_DIR` to the array instead of (or in addition to) `BULK_BAM`.

**Trade-off: keep buckets or merge to a single bulk BAM**

| Strategy | Per-cell extraction | Bulk QC (flagstat, IGV) | Inode count |
|----------|--------------------|-----|--------|
| Keep 576 bucket BAMs only | fast (one small BAM) | needs on-demand merge | high (576×2 per sample) |
| Merge all to one bulk BAM | slow (full scan per cell) | trivial | low |
| **Hybrid (recommended)** | fast | trivial | high |

The hybrid: keep the indexed per-bucket BAMs for the single-cell array, AND
produce a merged bulk BAM (using the existing two-stage parallel merge at
lines 312–401) for flagstat and any IGV/debug work. Disk cost ≈ 2×; the one
extra merge pass is cheap relative to alignment.

**Empty buckets and load imbalance**

- For small samples most of the 576 buckets are empty (28-byte gzip stubs).
  The bucket-list enumeration in step 1 already filters with
  `find -size +28c`; the BWA array tolerates empty input naturally
  (`bwa mem` on an empty FASTQ emits a header-only SAM, which merges cleanly).
- Bucket sizes are very heterogeneous — wells with many cells produce fat
  buckets, sparse wells produce thin ones. Expect 10–100× variation in array
  task duration. Set the BWA array `--time` to the worst-case bucket (not the
  average), or cap concurrency with `--array=1-N%K` if you're sharing a
  partition.

---

## CapGTA_PP.sh (STAR)

### Target architecture

```
spc-demux --sam-output sample.bam ──► one unaligned BAM with CB:Z tags
                                       │
                                       └──► STAR --readFilesType SAM SE ──► aligned BAM
                                                                             │
                                                                             └──► split DNA / RNA by N in CIGAR
```

The current SE script (`PP_array_gta_star_only_SE.sh`) already does
demux → trim → `fastq_to_unaligned_sam.py` → STAR. The first three steps
collapse into a single `spc-demux --sam-output` call.

### Orchestrator changes (`CapGTA_PP.sh`)

**Replace** the chunking block (lines 233–256) and the array submission (lines
261–270) with a single demux call followed by a single STAR job:

```bash
# NEW — one demux pass, one STAR call
SPC_DEMUX="${SCRIPTS_DIR}/scripts/utils/spc-demux"
SPC_CONFIG="${SCRIPTS_DIR}/config/spc_demux.yaml"

UNALIGNED_BAM="${TMP_DIR}/${SAMPLE_NAME}_unaligned.bam"

module load SAMtools
${SPC_DEMUX} \
    --config "${SPC_CONFIG}" \
    --r1 ${READ1} --r2 ${READ2} \
    --sample "${SAMPLE_NAME}" \
    --sam-output "${UNALIGNED_BAM}" \
    --stats "${RESULTS_DIR}/demux_stats.json" \
    --cell-counts "${RESULTS_DIR}/cell_counts.tsv"
# add --paired-output if R1 *and* R2 both carry genomic content
# (1x250 SPC: omit; both reads insert: add)
```

STAR on the single unaligned BAM:

```bash
module load STAR SAMtools

if [ "$SINGLE_END" = "true" ]; then
    READ_FILES_TYPE="SAM SE"
else
    READ_FILES_TYPE="SAM PE"
fi

STAR --runThreadN 16 \
     --genomeDir "${REFERENCE_GENOME}/STARIndex" \
     --readFilesIn "${UNALIGNED_BAM}" \
     --readFilesType ${READ_FILES_TYPE} \
     --readFilesCommand "samtools view" \
     --outFileNamePrefix "${TMP_DIR}/${SAMPLE_NAME}_" \
     --outSAMtype BAM Unsorted \
     --outSAMattributes NH HI AS nM NM MD CB BC RG \
     --outFilterMultimapNmax 1 \
     --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 \
     --alignIntronMax 15000
```

The two changes that matter:

1. `--readFilesCommand "samtools view"` so STAR pipes the BAM through samtools
   into its SAM parser. With `.sam` input you can drop this flag.
2. **Add `CB BC RG` to `--outSAMattributes`** so STAR propagates the tags into
   the aligned output. Without this the tags are stripped — the whole reason
   we're feeding STAR a SAM in the first place is lost.

### What can be deleted

- `PP_array_gta_star_only.sh` and `PP_array_gta_star_only_SE.sh` (per-chunk
  STAR — replaced by one STAR call)
- `fastq_to_unaligned_sam.py` (replaced by `spc-demux --sam-output`)
- The two-stage parallel SAM merge (lines 280–400) — no merge needed when STAR
  emits a single output BAM. Keep the DNA/RNA splitting awk block; it operates
  on the STAR BAM regardless of how that BAM was produced.
- `N_CHUNKS`, `--use-existing-chunks`, `chunk_indices.txt`

### If STAR throughput is the bottleneck

A single STAR call on a large sample can take a long time. If you need to
parallelize, the cleanest split is on the bucketed FASTQ side and then convert
per-bucket back to BAM — but at that point you're rebuilding the chunk-and-
align flow. Easier: bump `--runThreadN` to the node's full CPU count. STAR
parallelizes well within a single job.

---

## Common gotchas

1. **`scratch_dir` must exist and have headroom.** `spc-demux` writes
   uncompressed FASTQ buckets, then compresses them in parallel with gzip.
   Peak disk = ~2× the size of the compressed input. Set
   `output.scratch_dir` in the demux config to a local NVMe path.

2. **CB tag length changed from 45 → 32.** Audit any code that grep-matches
   `CB:Z:[ACGT]{45}`. Look in `assign_conditions.py` (uses `BC_SLICES` from
   `barcode_core` — already 32-char aware in the new repo, but the
   `SPC_genome` copy may be stale) and any per-cell extraction scripts that
   parse the CB into round barcodes.

3. **`-f 0x2` filter drops everything when R2 has no genomic content.**
   Replace with `-F 4` (mapped only) for the SE-flavored SPC structure.

4. **`samtools` is required at runtime for `.bam` output.** `spc-demux` shells
   out to `samtools view -b`. Either `module load SAMtools` before running, or
   write `.sam` and pipe yourself.

5. **`barcode_matching.ignore_linkers` matters for the assignment-rate stat.**
   Default in the demux repo is `true` (faster — only validates barcodes, not
   linkers). The old `atrandi_demux.py` `--ignore-overhangs` flag is the same
   knob. Set to `false` only if linker QC is part of your acceptance criteria.

6. **Trim Galore is gone.** `spc-demux` does R1 adapter trimming + R2-tail-RC
   trimming inline. The trimming config (`min_length`, `max_errors`,
   `min_overlap`, `quality_cutoff`) lives in the demux config.yaml. Note that
   `quality_cutoff` is currently **not implemented** in the Rust binary — only
   adapter trimming runs. If quality trimming matters for your samples, either
   keep a separate Trim Galore step downstream or wait for the trimming gap to
   be closed.

7. **Read-group propagation.** `spc-demux --sam-output` writes one `@RG` line
   (`ID:{sample} SM:{sample} PL:ILLUMINA`). The bucketed FASTQ output does
   **not** carry RG info — you need `-R` on the BWA-MEM command line as
   before. Both pipelines already do this.

8. **Bucket layout matters for per-cell work.** See the "Bucketed output is a
   partition of the cell space" subsection — the single-cell extraction array
   should look up the bucket from the CB instead of scanning the merged BAM.

---

## Validation steps before flipping the switch

1. **CB length:** `samtools view sample.bam | awk '{for(i=12;i<=NF;i++) if($i~/^CB:Z:/) print length($i)-5}' | sort -u` should print only `32`.
2. **Tag round-trip through aligner:** after BWA / STAR, `samtools view aligned.bam | head -10000 | grep -c 'CB:Z:'` should equal record count. If it's 0, you forgot `-C` (BWA) or `CB` in `--outSAMattributes` (STAR).
3. **Assignment rate sanity:** `demux_stats.json` reports `mapped_reads / total_reads`. Compare against the legacy `barcode_assignment_stats.txt` rate on the same input — should be within a couple percent (the new demux is stricter about ambiguous variants, so a small drop is expected).
4. **Cell count:** `wc -l demux_stats cell_counts.tsv` vs the legacy per-chunk demux-stats sum. Order-of-magnitude check only — exact match isn't expected because of the ambiguous-variant handling difference.
