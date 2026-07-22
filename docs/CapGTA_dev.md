# CapGTA Pipeline Development

Status of the CapGTA downstream processing. Updated 2026-07-22.

The upstream demultiplexing + STAR alignment steps have moved to the standalone
[`spc-demultiplex`](../../spc-demultiplex) and [`spc-align`](../../spc-align)
pipelines (Modules 2 and 3). This repo now handles only the *downstream* steps:
building an RNA count matrix and doing joint variant calling on the per-cell
BAMs produced by `spc-align`.

---

## Downstream pipeline (from spc-align outputs)

### Inputs
- A `real_cells.csv` with columns `barcode,bam_path` (one row per called cell).
  `bam_path` points at a coordinate-sorted, indexed per-cell BAM from
  `spc-align` (STAR-aligned, MarkDuplicates flagged, RG SM = barcode).
- The reference FASTA and GTF that `spc-align` used (iGenomes layout).

### Objective 1 — RNA count matrix

Two-stage: filter each per-cell BAM to spliced reads (CIGAR contains `N`), then
run `featureCounts` across all filtered BAMs.

Scripts (in `scripts/CapGTA/`):
- `filter_spliced_reads_array.sh` — SLURM array, one task per cell.
- `create_rna_count_matrix.sh` — single multi-threaded featureCounts job.
- `submit_rna_counts.sh` — orchestrator (submits the array + count job with `afterok`).

Submit:
```bash
sbatch --wrap='bash scripts/CapGTA/submit_rna_counts.sh \
    results/<sample>/real_cells.csv \
    /shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf \
    results/<sample>/rna_counts'
```

Output layout:
```
results/<sample>/rna_counts/
├── bam_list_source.txt
├── spliced_bams/{barcode}.bam[.bai]
├── bam_list_spliced.txt
├── rna_counts_raw.txt              # featureCounts raw
├── rna_counts_matrix.csv           # gene_id x barcode
└── rna_counts_summary.csv          # per-cell totals + genes detected
```

Notes:
- The BAMs are DNA+RNA coassay; spliced-only is a conservative first-pass RNA
  estimator. Expect low counts per cell (~0.2% of aligned reads are spliced).
- A more sophisticated estimator (Poisson exonic-enrichment model) is planned.

### Objective 2 — Joint variant calling

True joint calling with `bcftools mpileup -b bam_list -r region | bcftools call -mv`
parallelized by 1 Mb regions, then concatenated.

Scripts (in `scripts/CapGTA/`):
- `make_regions.sh` — split a `.fai` into fixed-size regions.
- `joint_variant_calling_array.sh` — SLURM array, one region per task, all cells jointly.
- `concat_region_vcfs.sh` — concat per-region VCFs into one joint VCF.
- `submit_joint_variant_calling.sh` — orchestrator (submits the array + concat with `afterok`).

Submit:
```bash
sbatch --wrap='bash scripts/CapGTA/submit_joint_variant_calling.sh \
    results/<sample>/real_cells.csv \
    /shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa \
    results/<sample>/joint_variants \
    1000000'
```

Output layout:
```
results/<sample>/joint_variants/
├── bam_list.txt
├── regions.txt
├── per_region/region_000001.vcf.gz ...
└── joint_variants.vcf.gz[.csi]     # final
```

Notes:
- Per-sample `FORMAT/AD,FORMAT/DP` are kept in the VCF for downstream cellspec analysis.
- Duplicate reads are skipped by default (`bcftools mpileup` default skip-flags).
- Mapping quality filter `-q 20` passes STAR uniquely-mapped reads (MAPQ 255)
  while excluding anything sneakier that gets in.
- `ulimit -n 8192` is raised inside the array script — 1660+ BAMs open at once
  per task.
- WBcel235 with 1 Mb regions → 104 array tasks.

---

## Deprecated (superseded by spc-align)

The following scripts implement the old in-repo demux/align/QC pipeline and are
obsolete now that Modules 2 and 3 handle that work. Kept for historical
reference; safe to delete once no notebook references them:

- `PP_array_gta.sh`, `PP_array_gta_star_only.sh` — STAR alignment array
- `sc_extract_preprocess_qc_array.sh` — per-cell extraction + MarkDuplicates + Picard + bigwig
- `calculate_exonic_enrichment.py` — exonic enrichment metric
- `build_wbcel235_star_index.sh` — STAR index build
- `submit_final_processing.sh` — old orchestrator (called per-cell variant calling + merge)
- `sc_variant_calling_bcftools_array.sh` — per-cell (not joint) variant calling
- `merge_sc_vcfs.sh` — `bcftools merge` union of per-cell VCFs (not true joint)

---

## Future work

### RNA rescue algorithm (Phase 3, previously planned)
Statistical reclassification of exonic reads from DNA→RNA using per-gene Poisson
enrichment tests. Was scoped against the old two-BAM (DNA/RNA) layout; needs to
be re-derived against the unified single-BAM layout from spc-align.

### HTML dashboard
Analogous to spc-align's `align_report.html`, but summarizing downstream
outputs: cells x reads, exonic enrichment, RNA counts per cell, variants per
cell, per-cell VAF distributions.
