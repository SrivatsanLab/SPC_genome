# CapGTA Pipeline Development

Status of `CapGTA_PP.sh` and supporting scripts. Updated 2026-05-26.

---

## Implemented

### Pipeline (`CapGTA_PP.sh`)
- `wait_for_job()` SLURM error handling (replaces `--dependency=afterok:`)
- `--use-existing-chunks` flag (skip FASTQ split, reuse `chunk_indices.txt`)
- `--debug` mode: capture unmapped reads with CB tags into a separate BAM, write per-barcode unmapped readcounts + knee plot, keep all intermediates
- Two-stage parallel SAM merge, DNA and RNA streams submitted concurrently (`submit_parallel_merge` / `finalize_stream_bam`)
- Unfiltered BAM + `flagstat_{dna,rna}.txt` before `-f 0x2` filtering
- DNA+RNA combined barcode assignment statistics → `barcode_assignment_stats.txt`
- Inline cell detection (knee plot on combined DNA+RNA counts)
- Inline QC compilation (Lorenz, Picard alignment metrics, exonic enrichment, Lander-Waterman)
- `SAMPLE_NAME=$(basename "${OUTPUT_NAME}")` to prevent BAM-path duplication when `OUTPUT_NAME` contains `/`

### Scripts
**`scripts/CapGTA/`**
- `sc_extract_preprocess_qc_array.sh` — unified per-cell extraction + MarkDuplicates + bigwig + Lorenz/gini + Picard + exonic enrichment (DNA & RNA)
- `calculate_exonic_enrichment.py` — `(exonic reads / total) − (exonic bases / genome bases)`
- `submit_final_processing.sh` — submits variant calling, VCF merge, RNA count matrix
- `PP_array_gta_star_only.sh` — STAR alignment array (DNA/RNA separation by splice junctions)
- `sc_variant_calling_bcftools_array.sh`, `merge_sc_vcfs.sh`, `create_rna_count_matrix.sh`

**`scripts/utils/`** (shared with CapWGS)
- `parallel_merge_array.sh`, `compile_lorenz.py`, `compile_qc_metrics.py`, `lorenz.py`
- `combine_readcounts_gta.py`, `compile_exonic_enrichment.py`

### Output structure
```
data/{sample}/
├── {sample}_dna.bam, {sample}_rna.bam
├── sc_outputs/{barcode}_{dna,rna}.bam, {barcode}_dna.bw, {barcode}.vcf.gz
└── qc_metrics/{barcode}_{alignment,duplication,gc,wgs}_metrics.txt
                {barcode}_lorenz.csv, {barcode}_gini.txt
                {barcode}_exonic_enrichment_{dna,rna}.txt

results/{sample}/
├── compiled_qc_metrics.csv, compiled_lorenz_curves.csv, compiled_exonic_enrichment.csv
├── lander_waterman_coverage_dna.png, kneeplot.png
├── readcounts{,_dna,_rna}.csv, real_cells.txt
├── barcode_assignment_stats.txt, flagstat_{dna,rna}.txt
└── sc_variants_merged.vcf.gz, rna_counts_matrix.csv
```

---

## TODO

### Phase 1 validation (still pending)
- Run on `worm_CapGTA_UDI_5` end-to-end
- Benchmark parallel-merge speedup vs old serial merge
- Verify exonic enrichment values (DNA ≈ 0, RNA > 0.2)
- Confirm barcode assignment rate (expect 60–80%)

### Phase 1b — RNA rescue plumbing (stubs)
- `--enable-rna-rescue` flag in `CapGTA_PP.sh`
- `scripts/CapGTA/rna_rescue_array.sh` (stub initially)
- `scripts/CapGTA/post_rescue_qc_array.sh` — re-run QC on rescued BAMs, write `*_exonic_enrichment_{dna,rna}_postrescue.txt`
- Conditional execution + post-rescue columns in `compiled_exonic_enrichment.csv`

### Phase 2 — Advanced QC
- Per-chunk preprocessing stats JSON from `PP_array_gta_star_only.sh` (demux/trim/alignment stats); inline compilation → `preprocessing_summary.json`, `adapter_histogram_{r1,r2}.csv`
- HTML dashboard (`generate_dashboard_gta.py`, or extend `scripts/CapWGS/generate_dashboard.py` with `--mode gta`): DNA/RNA read counts, exonic enrichment, RNA counts/cell, knee plot, Lorenz, Lander-Waterman

### Phase 3 — RNA rescue algorithm
Statistical reclassification of reads from DNA BAMs into RNA BAMs:

1. Estimate DNA baseline coverage λ_DNA from intergenic regions (>10 kb from genes), modeled as Poisson (or NB if overdispersed).
2. Per-gene exonic-enrichment test: Poisson rate test on `R_obs` vs `λ_DNA × exonic_length`, BH FDR across genes.
3. Per-gene RNA fraction `f_RNA = (R_excess + α) / (R_obs + α + β)` with Beta(α=3, β=2) shrinkage.
4. Probabilistic reassignment of exonic reads using gene's `f_RNA`. Keep splice-junction reads in RNA; intronic/intergenic in DNA; multi-gene overlaps use `max(f_RNA)`.

Implementation:
- `bin/classify_rna_from_dna.py` (core algorithm)
- `bin/reclassify_rna_array.sh` (per-cell SLURM array)
- Tunables: `min_intergenic_size=10000`, `fdr_threshold=0.05`, `alpha_prior=3`, `beta_prior=2`, `min_gene_coverage=10`

Variant calling after rescue: modify `sc_variant_calling_bcftools_array.sh` to accept both DNA and RNA BAMs per cell (`bcftools mpileup -f ref.fa cell_dna.bam cell_rna.bam | bcftools call …`) to recover coverage from reassigned reads without losing sensitivity.

Validation: RNA counts/cell should rise from ~2,141 → ~5,000–10,000; DNA intergenic coverage unchanged; gene expression correlates with WormSeq atlas; variant-calling quality maintained. When regenerating count matrices, use **top 1200 barcodes** (previous 1000 was conservative on the knee plot).
