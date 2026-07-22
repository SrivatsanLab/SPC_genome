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
- `counts_matrix_to_h5ad.py` — convert the matrix CSV into a sparse `.h5ad` for scanpy.
- `submit_rna_counts.sh` — orchestrator (array → count job → h5ad, chained with `afterok`).

Submit:
```bash
# Spliced-only (default; conservative RNA proxy)
bash scripts/CapGTA/submit_rna_counts.sh \
    results/<sample>/real_cells.csv \
    /shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf \
    results/<sample>/rna_counts

# All reads (skip spliced filter; higher counts, DNA contamination)
bash scripts/CapGTA/submit_rna_counts.sh --all-reads \
    results/<sample>/real_cells.csv \
    /shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf \
    results/<sample>/rna_counts_all_reads
```

Output layout:
```
results/<sample>/rna_counts/
├── bam_list_source.txt
├── spliced_bams/{barcode}.bam[.bai]
├── bam_list_spliced.txt
├── rna_counts_raw.txt              # featureCounts raw
├── rna_counts_matrix.csv           # gene_id x barcode
├── rna_counts_summary.csv          # per-cell totals + genes detected
└── rna_counts.h5ad                 # scanpy-ready AnnData (cells x genes, sparse)
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

### Objective 3 — RNA expression estimate (excess over NoFeatures baseline)

Neither spliced-only (too sparse) nor all-reads (dominated by WGA DNA background) is a
great per-gene expression proxy on its own. This step models the DNA baseline explicitly
and returns the excess:

```
non_exonic_bp = sum(chromosome lengths) - sum(featureCounts Length column)
lambda[c]     = Unassigned_NoFeatures[c] / non_exonic_bp
R_exp[c, g]   = lambda[c] * exonic_length[g]
R_rna[c, g]   = max(0, R_obs[c, g] - R_exp[c, g])
```

Where `R_obs`, `exonic_length`, and per-cell `Unassigned_NoFeatures` all come from the
same all-reads featureCounts run. Using NoFeatures as the baseline (rather than a
separately-counted intergenic bed) matches R_obs's fragment definition exactly and
includes intronic sequence — important for WGA data where amplification is
biased toward gene bodies vs. intergenic.

Scripts (in `scripts/CapGTA/`):
- `assemble_expression_matrix.py` — reads featureCounts raw + .summary + .fai; writes CSV + h5ad.
- `submit_rna_expression.sh` — one-line sbatch wrapper around the assembler.

Submit (requires all-reads featureCounts to have already run):
```bash
bash scripts/CapGTA/submit_rna_expression.sh \
    /shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa \
    results/<sample>/rna_counts_all_reads/rna_counts_raw.txt \
    results/<sample>/rna_expression
```

Output layout:
```
results/<sample>/rna_expression/
├── rna_expression_matrix.csv       # gene_id x barcode, float R_rna
└── rna_expression.h5ad             # X = R_rna, layers = raw_exon / expected
```

Notes:
- The h5ad `.uns` records `genome_bp`, `exonic_bp`, `non_exonic_bp`, and `model = noFeatures_excess`.
- Room to grow (v2): per-gene overdispersion (NB rather than Poisson-implied), Beta shrinkage on `f_RNA = (R_obs - R_exp) / R_obs`, per-gene FDR flag.

---

## Future work

### Expression model v2
Extend the simple-excess estimator: Beta shrinkage on the per-gene RNA fraction,
per-gene NB fit for overdispersion, per-gene FDR flag for "expressed above baseline".
See `assemble_expression_matrix.py` for the extension point.

### HTML dashboard
Analogous to spc-align's `align_report.html`, but summarizing downstream
outputs: cells x reads, exonic enrichment, RNA counts per cell, variants per
cell, per-cell VAF distributions.
