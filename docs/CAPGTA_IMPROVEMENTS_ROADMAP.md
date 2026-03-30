# CapGTA Pipeline Improvements - Implementation Roadmap

## Goal
Bring CapGTA pipeline up to parity with CapWGS improvements:
- Parallel two-stage SAM/BAM merge (20-40x speedup)
- Comprehensive QC outputs (Picard metrics, Lorenz curves, Lander-Waterman plots)
- Exonic enrichment metrics (pre/post RNA rescue)
- Unified single-cell processing array (avoid race conditions)
- Integrate RNA rescue algorithm (optional, in development)

---

## Updated Pipeline Flow

### 1. Parse Arguments and Setup
**Script**: `CapGTA_PP.sh` (inline)
**Status**: **MODIFY EXISTING**

- Add `--use-existing-chunks` flag (like CapWGS)
- Add `--enable-rna-rescue` flag (optional, default: disabled)
- Match CapWGS argument parsing structure
- Create directory structure matching CapWGS (data/, sc_outputs/, qc_metrics/, results/)

---

### 2. Chunk FASTQs
**Script**: `CapGTA_PP.sh` (inline)
**Status**: **MODIFY EXISTING** (copy logic from CapWGS)

- Reuse FASTQ chunking logic from CapWGS_PP.sh (lines 191-237)
- Support `--use-existing-chunks` flag
- Create `chunk_indices.txt` for array job submission

---

### 3. STAR Alignment Array
**Script**: `scripts/CapGTA/PP_array_gta_star_only.sh`
**Status**: **EXISTS - NO CHANGES NEEDED**

- Submit array job: `sbatch --array=1-$chunk_count`
- Wait for completion using `wait_for_job()` function
- Outputs: DNA and RNA SAM chunks in temp directory

---

### 4. Two-Stage Parallel SAM Merge (DNA)
**Script**: `CapGTA_PP.sh` (inline) + `scripts/utils/parallel_merge_array.sh`
**Status**: **MOVE SHARED SCRIPT + ADD TO CapGTA_PP.sh**

**Stage 1: Create merge groups**
- List DNA SAM files: `ls ${TMP_DIR}/*_dna.sam > ${RESULTS_DIR}/sam_list_dna.txt`
- Split into groups of 20: `split -l 20 -d -a 2 sam_list_dna.txt ${TMP_DIR}/merge_groups_dna/group_`

**Stage 2: Submit parallel merge array**
- Submit: `sbatch --array=1-${N_GROUPS} scripts/utils/parallel_merge_array.sh ${TMP_DIR}/merge_groups_dna ${TMP_DIR}/intermediate_bams_dna ${SAMPLE_NAME}_dna`
- Wait for completion

**Stage 3: Final merge**
- List intermediate BAMs: `ls ${TMP_DIR}/intermediate_bams_dna/*.bam > ${RESULTS_DIR}/intermediate_bam_list_dna.txt`
- Final merge: `samtools merge -@ 4 -c -b intermediate_bam_list_dna.txt -O BAM - | samtools view -@ 4 -f 0x2 -b - | samtools sort -@ 4 -o ${DATA_DIR}/${SAMPLE_NAME}_dna.bam`
- Index: `samtools index ${DATA_DIR}/${SAMPLE_NAME}_dna.bam`

---

### 5. Two-Stage Parallel SAM Merge (RNA)
**Script**: `CapGTA_PP.sh` (inline) + `scripts/utils/parallel_merge_array.sh`
**Status**: **MOVE SHARED SCRIPT + ADD TO CapGTA_PP.sh**

**Same as Step 4, but for RNA SAMs:**
- List RNA SAM files: `ls ${TMP_DIR}/*_rna.sam > ${RESULTS_DIR}/sam_list_rna.txt`
- Split, merge, and sort → `${DATA_DIR}/${SAMPLE_NAME}_rna.bam`

**Note**: DNA and RNA parallel merge arrays can be submitted simultaneously and run in parallel.

---

### 6. Compute Read Counts and Detect Cells
**Script**: `CapGTA_PP.sh` (inline)
**Status**: **MODIFY EXISTING** (currently in concatenate_gta.sh)

- Compute DNA read counts: `samtools view ${DNA_BAM} | python scripts/utils/readcounts.py -o ${RESULTS_DIR}/readcounts_dna.csv`
- Compute RNA read counts: `samtools view ${RNA_BAM} | python scripts/utils/readcounts.py -o ${RESULTS_DIR}/readcounts_rna.csv`
- Combine DNA + RNA counts: `python scripts/utils/combine_readcounts_gta.py`
- Detect cells: `python scripts/utils/detect_cells.py --plot ${RESULTS_DIR}/kneeplot.png > ${RESULTS_DIR}/real_cells.txt`
- Calculate barcode assignment statistics (like CapWGS)

---

### 7. Unified Single-Cell Processing Array (Pre-Rescue)
**Script**: `scripts/CapGTA/sc_extract_preprocess_qc_array.sh`
**Status**: **CREATE NEW** (based on CapWGS version)

**Submitted as**: `sbatch --array=1-${DETECTED_CELL_COUNT}`

**Per-cell tasks (in order):**
1. Extract DNA BAM from bulk: `samtools view -b ${BULK_DNA_BAM} -d CB:${BARCODE} > ${SC_OUTPUTS_DIR}/${BARCODE}_dna.bam`
2. Extract RNA BAM from bulk: `samtools view -b ${BULK_RNA_BAM} -d CB:${BARCODE} > ${SC_OUTPUTS_DIR}/${BARCODE}_rna.bam`
3. Run MarkDuplicates on DNA BAM (QC only, output to qc_metrics):
   - Input: `${BARCODE}_dna.bam`
   - Output: `${QC_METRICS_DIR}/${BARCODE}_duplication_metrics.txt`
   - Keep original BAM (no REMOVE_DUPLICATES)
4. Generate DNA bigwig: `bamCoverage -b ${BARCODE}_dna.bam -o ${SC_OUTPUTS_DIR}/${BARCODE}_dna.bw`
5. Generate DNA Lorenz curve + gini: `python scripts/utils/lorenz.py ${BARCODE}_dna.bw -o ${QC_METRICS_DIR}/${BARCODE}_lorenz.csv --gini-output ${QC_METRICS_DIR}/${BARCODE}_gini.txt`
6. Collect Picard alignment metrics (DNA only): `picard CollectAlignmentSummaryMetrics`
7. Calculate exonic enrichment (pre-rescue):
   - DNA: `python scripts/CapGTA/calculate_exonic_enrichment.py ${BARCODE}_dna.bam ${GTF} ${GENOME_SIZE} > ${QC_METRICS_DIR}/${BARCODE}_exonic_enrichment_dna.txt`
   - RNA: `python scripts/CapGTA/calculate_exonic_enrichment.py ${BARCODE}_rna.bam ${GTF} ${GENOME_SIZE} > ${QC_METRICS_DIR}/${BARCODE}_exonic_enrichment_rna.txt`

**Outputs per cell:**
- `sc_outputs/${BARCODE}_dna.bam` (original, will be overwritten if rescue runs)
- `sc_outputs/${BARCODE}_rna.bam` (original, will be overwritten if rescue runs)
- `sc_outputs/${BARCODE}_dna.bw`
- `qc_metrics/${BARCODE}_duplication_metrics.txt`
- `qc_metrics/${BARCODE}_alignment_metrics.txt`
- `qc_metrics/${BARCODE}_lorenz.csv`
- `qc_metrics/${BARCODE}_gini.txt`
- `qc_metrics/${BARCODE}_exonic_enrichment_dna.txt`
- `qc_metrics/${BARCODE}_exonic_enrichment_rna.txt`

Wait for completion.

---

### 8. RNA Rescue Array [OPTIONAL]
**Script**: `scripts/CapGTA/rna_rescue_array.sh`
**Status**: **CREATE NEW STUB** (implement algorithm later)

**Condition**: Only run if `--enable-rna-rescue` flag is set

**Submitted as**: `sbatch --array=1-${DETECTED_CELL_COUNT}`

**Per-cell tasks:**
1. Run statistical reclassification algorithm:
   - Input: `${BARCODE}_dna.bam`, `${BARCODE}_rna.bam`, GTF annotation
   - Estimate λ_DNA from intergenic coverage
   - Compute per-gene RNA fraction (f_RNA)
   - Probabilistically reassign unspliced exonic reads
   - Output: **Overwrite** `${BARCODE}_dna.bam` and `${BARCODE}_rna.bam`

**For Phase 1 (stub)**: Just copy files to simulate completion
```bash
# Stub version - just verify files exist
if [ -f "${SC_OUTPUTS_DIR}/${BARCODE}_dna.bam" ] && [ -f "${SC_OUTPUTS_DIR}/${BARCODE}_rna.bam" ]; then
    echo "RNA rescue would run here (stub)"
fi
```

Wait for completion.

---

### 9. Post-Rescue QC Array [OPTIONAL]
**Script**: `scripts/CapGTA/post_rescue_qc_array.sh`
**Status**: **CREATE NEW**

**Condition**: Only run if RNA rescue was enabled and completed

**Submitted as**: `sbatch --array=1-${DETECTED_CELL_COUNT}`

**Per-cell tasks (re-run QC on rescued BAMs):**
1. Re-run MarkDuplicates on rescued DNA BAM → **Overwrite** `${QC_METRICS_DIR}/${BARCODE}_duplication_metrics.txt`
2. Re-generate DNA bigwig → **Overwrite** `${SC_OUTPUTS_DIR}/${BARCODE}_dna.bw`
3. Re-generate DNA Lorenz curve + gini → **Overwrite** `${QC_METRICS_DIR}/${BARCODE}_lorenz.csv` and `${QC_METRICS_DIR}/${BARCODE}_gini.txt`
4. Re-collect Picard alignment metrics → **Overwrite** `${QC_METRICS_DIR}/${BARCODE}_alignment_metrics.txt`
5. Calculate exonic enrichment (post-rescue):
   - DNA: `python scripts/CapGTA/calculate_exonic_enrichment.py ${BARCODE}_dna.bam ${GTF} ${GENOME_SIZE} > ${QC_METRICS_DIR}/${BARCODE}_exonic_enrichment_dna_postrescue.txt`
   - RNA: `python scripts/CapGTA/calculate_exonic_enrichment.py ${BARCODE}_rna.bam ${GTF} ${GENOME_SIZE} > ${QC_METRICS_DIR}/${BARCODE}_exonic_enrichment_rna_postrescue.txt`

Wait for completion.

---

### 10. Compile QC Metrics
**Script**: `CapGTA_PP.sh` (inline)
**Status**: **ADD NEW SECTION** (copy from CapWGS)

**Uses shared scripts from scripts/utils/:**
- `compile_lorenz.py` - Compile Lorenz curves
- `compile_qc_metrics.py` - Compile Picard metrics + gini coefficients
- `plot_lander_waterman.py` - Generate Lander-Waterman plot
- `compile_exonic_enrichment.py` - **CREATE NEW** - Compile exonic enrichment metrics

**Outputs:**
- `results/${SAMPLE}/compiled_qc_metrics.csv` (includes gini coefficient)
- `results/${SAMPLE}/compiled_lorenz_curves.csv`
- `results/${SAMPLE}/compiled_exonic_enrichment.csv` (columns: barcode, dna_pre, rna_pre, dna_post, rna_post)
- `results/${SAMPLE}/lander_waterman_coverage.png`

---

### 11. Variant Calling Array
**Script**: `scripts/CapGTA/sc_variant_calling_bcftools_array.sh`
**Status**: **EXISTS - NO CHANGES NEEDED**

- Calls variants on DNA BAMs (post-rescue if available, else original)
- BCFtools only (no GATK option)
- Submit: `sbatch --array=1-${DETECTED_CELL_COUNT}`
- Outputs: `sc_outputs/${BARCODE}.vcf.gz`

Wait for completion.

---

### 12. Merge VCFs and Generate RNA Count Matrix
**Script**: `scripts/CapGTA/merge_sc_vcfs.sh` + `scripts/CapGTA/create_rna_count_matrix.sh`
**Status**: **EXISTS - NO CHANGES NEEDED**

**Submitted as two separate jobs** (or combined into one submission script):
1. Merge per-cell VCFs → `results/${SAMPLE}/sc_variants_merged.vcf.gz`
2. Generate RNA count matrix from RNA BAMs → `results/${SAMPLE}/rna_counts_matrix.csv`

These are final outputs, no wait needed.

---

## Script Organization Changes

### Scripts to Move to scripts/utils/:

**From scripts/CapWGS/ → scripts/utils/:**
1. `parallel_merge_array.sh` - Already generic, shared by both pipelines
2. `compile_lorenz.py` - Generic QC compilation
3. `compile_qc_metrics.py` - Generic QC compilation
4. `lorenz.py` - Generic Lorenz curve generation

**Update references in:**
- `CapWGS_PP.sh` - Update paths to scripts/utils/
- `scripts/CapWGS/sc_extract_preprocess_qc_array.sh` - Update lorenz.py path

---

### New Scripts to Create:

**scripts/CapGTA/:**
1. `sc_extract_preprocess_qc_array.sh` - Unified single-cell processing (based on CapWGS version)
2. `rna_rescue_array.sh` - RNA reclassification (stub for Phase 1)
3. `post_rescue_qc_array.sh` - Re-run QC after rescue
4. `calculate_exonic_enrichment.py` - Calculate exonic enrichment metric

**scripts/utils/:**
1. `combine_readcounts_gta.py` - Combine DNA + RNA read counts for cell detection
2. `compile_exonic_enrichment.py` - Compile exonic enrichment metrics into CSV

---

### Scripts to Modify:

**Major modifications:**
1. `CapGTA_PP.sh` - Complete rewrite to match CapWGS structure

**Minor modifications:**
2. `CapWGS_PP.sh` - Update script paths to scripts/utils/
3. `scripts/CapWGS/sc_extract_preprocess_qc_array.sh` - Update lorenz.py path

---

### Scripts to Delete:

1. `scripts/CapGTA/concatenate_gta.sh` - Logic moved inline to main script
2. `scripts/CapGTA/sc_from_bam_gta.sh` - Replaced by unified array

---

## Exonic Enrichment Metric Definition

**Formula:**
```
Exonic Enrichment = (Exonic Reads / Total Reads) - (Exonic Bases / Total Reference Bases)
```

**Interpretation:**
- **Positive value**: Enrichment for exonic regions (expected for RNA)
- **~Zero**: No enrichment (expected for DNA pre-rescue)
- **Negative value**: Depletion of exonic regions (unexpected)

**Example:**
- Cell has 100,000 total reads, 30,000 map to exons
- Reference genome: 100 Mb total, 2 Mb exonic (2%)
- Exonic reads proportion: 30,000 / 100,000 = 0.30 (30%)
- Expected exonic proportion: 2 Mb / 100 Mb = 0.02 (2%)
- **Exonic enrichment = 0.30 - 0.02 = 0.28 (28% enrichment)**

**Implementation requirements:**
- GTF file for exon coordinates
- Reference genome size (can be calculated from FASTA index)

---

## Directory Structure (Final)

```
data/{sample}/
├── {sample}_dna.bam              # Bulk DNA BAM
├── {sample}_rna.bam              # Bulk RNA BAM
├── sc_outputs/                   # Single-cell outputs
│   ├── {barcode}_dna.bam         # DNA BAM (post-rescue if rescue run)
│   ├── {barcode}_rna.bam         # RNA BAM (post-rescue if rescue run)
│   ├── {barcode}_dna.bw          # DNA bigwig (post-rescue if rescue run)
│   └── {barcode}.vcf.gz          # Variant calls
└── qc_metrics/                   # QC outputs (overwritten if rescue run)
    ├── {barcode}_alignment_metrics.txt           # Picard metrics
    ├── {barcode}_duplication_metrics.txt
    ├── {barcode}_lorenz.csv
    ├── {barcode}_gini.txt
    ├── {barcode}_exonic_enrichment_dna.txt       # Pre-rescue DNA
    ├── {barcode}_exonic_enrichment_rna.txt       # Pre-rescue RNA
    ├── {barcode}_exonic_enrichment_dna_postrescue.txt  # Post-rescue DNA
    └── {barcode}_exonic_enrichment_rna_postrescue.txt  # Post-rescue RNA

results/{sample}/
├── compiled_qc_metrics.csv           # Picard metrics + gini coefficients
├── compiled_lorenz_curves.csv
├── compiled_exonic_enrichment.csv    # Pre/post rescue for DNA and RNA
├── lander_waterman_coverage.png      # DNA coverage plot
├── kneeplot.png
├── readcounts.csv                    # Combined DNA + RNA counts
├── real_cells.txt
├── barcode_assignment_stats.txt
├── sc_variants_merged.vcf.gz         # Joint VCF
└── rna_counts_matrix.csv             # RNA count matrix

/hpc/temp/.../SPC_genome_preprocessing/{sample}/
├── read1_chunk_*, read2_chunk_*      # FASTQ chunks
├── *_dna.sam, *_rna.sam              # Alignment chunks
├── merge_groups_dna/                 # DNA merge groups
├── merge_groups_rna/                 # RNA merge groups
├── intermediate_bams_dna/            # DNA intermediate BAMs
└── intermediate_bams_rna/            # RNA intermediate BAMs
```

---

## Implementation Phases

### Phase 1: Parallel Merge + Basic QC (Priority 1)
**Goal**: Bring CapGTA to parity with CapWGS improvements

**Tasks:**
1. Move shared scripts to `scripts/utils/` and update references
2. Rewrite `CapGTA_PP.sh` main script
3. Create `scripts/CapGTA/sc_extract_preprocess_qc_array.sh`
4. Create `scripts/CapGTA/calculate_exonic_enrichment.py`
5. Create `scripts/utils/combine_readcounts_gta.py`
6. Create `scripts/utils/compile_exonic_enrichment.py`
7. Add QC compilation section to CapGTA_PP.sh
8. Delete obsolete scripts (concatenate_gta.sh, sc_from_bam_gta.sh)

**Testing:**
- Run on existing worm CapGTA data
- Verify parallel merge speedup
- Validate QC outputs
- Compare exonic enrichment pre-rescue (DNA should be ~0, RNA should be positive)

**Estimated time**: 1-2 days of development + testing

---

### Phase 2: RNA Rescue Integration (Priority 2)
**Goal**: Add RNA rescue pipeline flow (stub implementation)

**Tasks:**
1. Create `scripts/CapGTA/rna_rescue_array.sh` (stub - just copies files)
2. Add `--enable-rna-rescue` flag to CapGTA_PP.sh
3. Create `scripts/CapGTA/post_rescue_qc_array.sh`
4. Add conditional execution logic to main pipeline
5. Test pipeline flow with stubbed rescue step

**Testing:**
- Run with `--enable-rna-rescue` flag
- Verify QC outputs are overwritten correctly
- Confirm post-rescue metrics are generated
- Pipeline should complete successfully with stub

**Estimated time**: 4-6 hours

---

### Phase 3: RNA Rescue Algorithm (Priority 3 - Future Work)
**Goal**: Implement statistical reclassification algorithm

**Tasks:**
1. Implement algorithm in `rna_rescue_array.sh`:
   - Parse GTF for genes, exons, intergenic regions
   - Estimate λ_DNA from intergenic coverage
   - Compute per-gene RNA fraction (f_RNA) with statistical tests
   - Probabilistically reassign unspliced exonic reads
2. Validate on test cells (5-10 cells first)
3. Full production run with all cells

**Testing:**
- λ_DNA estimates reasonable (0.01-0.1 reads/base)
- f_RNA distributions sensible (0-0.9 range)
- RNA counts improve significantly (2,141 → 5,000-10,000)
- DNA intergenic coverage unchanged (sanity check)
- Exonic enrichment post-rescue shows improvement
- Variant calling quality maintained

**Estimated time**: 3-5 days (algorithm development + validation)

---

## Success Metrics

### Phase 1 Completion:
- ✅ Parallel merge reduces bulk BAM creation time by 20-40x
- ✅ All Picard QC metrics generated for DNA BAMs
- ✅ Lorenz curves + gini coefficients compiled
- ✅ Lander-Waterman plot generated
- ✅ Exonic enrichment metrics show DNA ~0, RNA >0.2

### Phase 2 Completion:
- ✅ Pipeline runs successfully with `--enable-rna-rescue` flag
- ✅ QC outputs correctly overwritten after rescue
- ✅ Post-rescue exonic enrichment metrics generated

### Phase 3 Completion:
- ✅ RNA counts per cell increase from ~2,000 to ~5,000-10,000
- ✅ DNA exonic enrichment remains near zero post-rescue
- ✅ RNA exonic enrichment increases post-rescue
- ✅ Variant calling sensitivity maintained
- ✅ Gene expression correlates with reference atlas
