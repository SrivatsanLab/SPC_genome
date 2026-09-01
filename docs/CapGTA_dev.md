# CapGTA Pipeline Development

Status of the CapGTA downstream processing. Updated 2026-09-01.

The upstream demultiplexing + STAR alignment steps have moved to the standalone
[`spc-demultiplex`](../../spc-demultiplex) and [`spc-align`](../../spc-align)
pipelines (Modules 2 and 3). This repo now handles only the *downstream* steps:
building an RNA count matrix (with SoupX decontamination) and doing joint variant
calling on the per-cell BAMs produced by `spc-align`.

Two orchestrators live at the repo root, one per objective:

- `CapGTA_gex_soupx_pipeline.sh`    — Objective 1 (GEX + SoupX + preprocess)
- `CapGTA_joint_variant_calling.sh` — Objective 2 (single-cell bcftools joint calling)

Objective 3 (RNA expression estimator) is a small post-hoc job on the outputs
of Objective 1's all-reads mode.

Per-worm pseudobulk variant calling (post-demultiplex, GATK) is a separate
concern; see `scripts/CapGTA/pseudobulk_gatk/` and `docs/GATK_PSEUDOBULK_PLAN.md`.

---

## Downstream pipeline (from spc-align outputs)

### Inputs
- A `real_cells.csv` with columns `barcode,bam_path` (one row per called cell).
  `bam_path` points at a coordinate-sorted, indexed per-cell BAM from
  `spc-align` (STAR-aligned, MarkDuplicates flagged, RG SM = barcode).
- The reference FASTA and GTF that `spc-align` used (iGenomes layout).

### Objective 1 — GEX count matrix with SoupX decontamination

Eight-stage chain that goes from real-cell BAMs to a SoupX-decontaminated,
Leiden-clustered scanpy h5ad. Uses spliced-only reads (CIGAR contains `N`) as
the RNA proxy. A dev variant that also accepts polyA/T-containing reads lives
under `CapGTA_dev/` — see § Future work.

Stages (chained via SLURM `afterok`):

0. **Junction + length tables** — `gtf_gene_junctions.py` + `gtf_gene_lengths.py` (inline).
1. **Real-cell spliced counts** — `submit_rna_counts.sh` (array: `filter_spliced_reads_array.sh` → `create_rna_count_matrix.sh` → `counts_matrix_to_h5ad.py`).
2. **Preprocess (pre-SoupX)** — `preprocess_scanpy.py`: QC, WBID→symbol rename (WBID kept in `.var['wormbase_id']`), junction/length correction on `.X`, log1p, HVG, PCA, neighbors, UMAP, Leiden. The Leiden partition is what SoupX consumes.
3. **Enumerate empties** — `enumerate_empties.py` samples N empty-droplet BAMs from the same run directories as the real cells.
4. **Empty-droplet spliced counts** — `submit_rna_counts.sh` again, on the empties CSV.
5. **Prep SoupX inputs** — `prep_soupx_inputs.py`: WGA-aware QC on empties + gene-order alignment between the real and empty matrices.
6. **SoupX** — `soupx_decontaminate.R` runs `autoEstCont` + `adjustCounts` (R via `ml fhR`).
7. **Reimport** — `reimport_soupx.py` writes SoupX-adjusted counts back into an h5ad.
8. **Preprocess (post-SoupX)** — `preprocess_scanpy.py` again with `--drop-cuticle` (removes `col-*`, `dpy-7`, `sqt-1/3`, `rol-6`, `bli-1/2` before renormalization); writes the terminal h5ad.

Submit:
```bash
bash CapGTA_gex_soupx_pipeline.sh \
    results/<sample>/real_cells.csv \
    /shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf \
    results/<sample>/gex \
    --n-empties 5000 --resolution 2
```

Terminal output: `results/<sample>/gex/adata_gex_spliced_soupx_processed.h5ad`

Full layout under `<output_dir>`:
```
gex/
├── gene_junctions.csv                                # exon–exon junctions per gene
├── gene_lengths.csv                                  # exonic length per gene
├── rna_counts/                                       # real cells (spliced)
│   └── rna_counts.h5ad
├── adata_gex_spliced.h5ad                            # stage 2 (preprocess pre-SoupX)
├── soupx_empties/
│   ├── empty_cells.csv
│   └── rna_counts/                                   # empty-droplet spliced counts
│       └── rna_counts.h5ad
├── soupx_in/, soupx_out/                             # SoupX I/O
├── adata_gex_spliced_soupx.h5ad                      # reimported SoupX counts
└── adata_gex_spliced_soupx_processed.h5ad           # terminal (cuticle dropped)
```

Notes:
- `.layers['counts']` on every preprocess output holds raw integer counts as-loaded — this is what SoupX and downstream count-model tools expect. The junction/length correction is applied to `.X` only.
- BAMs are DNA+RNA coassay; spliced-only is a conservative RNA proxy (~0.2% of aligned reads are spliced). See Future work for the polyA/T-inclusive variant.
- `submit_rna_counts.sh` is still callable directly for one-off count matrices, including its `--all-reads` mode (skips the spliced filter, feeds Objective 3).

### Objective 2 — Joint variant calling (single-cell, bcftools)

True joint calling with `bcftools mpileup -b bam_list -r region | bcftools call -mv`,
parallelized by 1 Mb regions and concatenated. Coverage per cell is ~1×, which
is why GATK is **not** used here — HaplotypeCaller's assembly-based calling
needs materially more depth than a single cell provides. GATK re-enters the
picture only after cells are demultiplexed into per-worm pseudobulks (see
`scripts/CapGTA/pseudobulk_gatk/`).

Stages (chained via SLURM `afterok`):
1. `make_regions.sh` — split the reference `.fai` into fixed-size regions.
2. `joint_variant_calling_array.sh` — SLURM array, one region per task, all cells jointly.
3. `concat_region_vcfs.sh` — concat per-region VCFs into one joint VCF.

Submit:
```bash
bash CapGTA_joint_variant_calling.sh \
    results/<sample>/real_cells.csv \
    /shared/biodata/reference/iGenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa \
    results/<sample>/joint_variants \
    1000000
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
- `ulimit -n 8192` is raised inside the array script — 1660+ BAMs open at once per task.
- WBcel235 with 1 Mb regions → 104 array tasks.

### Objective 3 — RNA expression estimate (excess over NoFeatures baseline)

Neither spliced-only (too sparse) nor all-reads (dominated by WGA DNA background)
is a great per-gene expression proxy on its own. This step models the DNA baseline
explicitly and returns the excess:

```
non_exonic_bp = sum(chromosome lengths) - sum(featureCounts Length column)
lambda[c]     = Unassigned_NoFeatures[c] / non_exonic_bp
R_exp[c, g]   = lambda[c] * exonic_length[g]
R_rna[c, g]   = max(0, R_obs[c, g] - R_exp[c, g])
```

Where `R_obs`, `exonic_length`, and per-cell `Unassigned_NoFeatures` all come from the
same all-reads featureCounts run. Using NoFeatures as the baseline (rather than a
separately-counted intergenic bed) matches R_obs's fragment definition exactly and
includes intronic sequence — important for WGA data where amplification is biased
toward gene bodies vs. intergenic.

Scripts (in `scripts/CapGTA/`):
- `assemble_expression_matrix.py` — reads featureCounts raw + .summary + .fai; writes CSV + h5ad.
- `submit_rna_expression.sh` — one-line sbatch wrapper around the assembler.

Submit (requires an all-reads featureCounts run — see § Objective 1's `submit_rna_counts.sh --all-reads`):
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
- Room to grow (v2): per-gene overdispersion (NB rather than Poisson-implied), Beta shrinkage on `f_RNA = (R_obs - R_exp) / R_obs`, per-gene FDR flag. The `CapGTA_dev/` sweep (§ Future work) is where v2 is being calibrated.

---

## Analysis utilities

Standalone scripts under `scripts/CapGTA/` that operate on a terminal GEX h5ad
(post-SoupX preprocess). Not part of the SLURM chain — invoked per analysis.

### polyA/T per-cell QC
- `count_polyat_per_cell.py` (+ `.sh` sbatch wrapper) — count reads with ≥N-bp A/T run per cell (samtools, process pool). Same rule as the `CapGTA_dev/` spliced+polyA/T RNA filter, so this is the natural QC. Two-step so the 6 GB h5ad is never fully rewritten: (1) counting → resumable CSV, (2) `--inject` writes the counts as a new `.obs` dataset via h5py (~20 KB touched).

### WormBase annotation
- `build_wormbase_anatomy_dict.py` — download the current WormBase anatomy GAF + OBO and pivot to a WBID × WBbt-name binary matrix. Drop-in replacement for the `tissue_enrichment_analysis` package's broken `fetch_dictionary()` (the Caltech URL 404s).
- `build_expanded_markers.py` — build tissue-specific marker sets from Cao 2017 L2 sci-RNA-seq (117 fine cell types → 14 broad tissues; anything unmapped → "neuron"). Wilcoxon vs rest, cross-tissue specificity filter, top-N WBGene IDs per tissue → JSON.
- `deg_wormbase_annotate.py` — `sc.tl.rank_genes_groups` per cluster, pull top-N genes, query WormBase REST for gene name / class / concise description → CSV.
- `tea_cluster_enrichment.py` — WormBase tissue / phenotype / GO enrichment (TEA) on top-N DEGs per cluster via the `tissue_enrichment_analysis` package. Pair with `build_wormbase_anatomy_dict.py`'s output via `--dict-csv` because the package's built-in fetcher is broken.

---

## Future work

### Expression model v2 — baseline-calibration sweep

Under `CapGTA_dev/scripts/`, an experimental framework compares baseline choices
and models against a per-gene spliced-read ground truth. The point is to pick
the successor to Objective 3's simple NoFeatures-excess estimator.

- `filter_polyat_reads.py`, `filter_sl_reads.py`, `filter_spliced_polyat_reads.py` (each with a `_array.sh` SLURM array wrapper) — per-cell BAM filters for polyA/T-containing reads (≥N-bp A/T run in SEQ, including soft-clip), SL1/SL2 seed-matched reads (trans-spliced leader), and spliced ∪ polyA/T reads.
- `build_local_baselines.py` — per-gene local gDNA baselines (merged intron span; ±window flanking region minus any exon/intron overlap) → featureCounts-compatible SAFs → counts.
- `build_calibration_table.py` — assemble a cells × genes AnnData with `.X = exon_count` (all-reads featureCounts) and `.layers[{spliced, polyat, sl, intron, flank}]`, plus `.var[{exon_bp, intron_bp, flanking_bp, n_junc}]` and `.obs[{spliced_libsize, polyat_libsize, sl_libsize}]`.
- `models.py` — expression models sharing a `fit(adata) / predict(adata) → (n_cells, n_genes)` interface.
- `run_sweep.py` — sweep baseline × floor_alpha × model × features. **Primary metric is per-gene Pearson/Spearman** between predicted counts summed over cells and observed spliced-read counts, broken down by coverage tier (`pct_1x` from the QC metrics table). Cell-level metrics are diagnostic-only — Poisson thinning of spliced counts makes cell-level correlations dominated by noise, not model quality (see `sweep_v1_notes.md`).

### Spliced ∪ polyA/T GEX pipeline (dev)

`CapGTA_dev/CapGTA_gex_soupx_pipeline_spliced_polyat.sh` is a parallel version
of the main Objective 1 pipeline that broadens the RNA-read definition to
spliced ∪ polyA/T (via `CapGTA_dev/scripts/submit_rna_counts_spliced_polyat.sh`
+ `filter_spliced_polyat_reads*`). All other stages — preprocess, empties,
SoupX, reimport, post-SoupX preprocess — reuse the main-pipeline scripts
unchanged. Outputs are named `adata_gex_splicedpolyat_soupx_processed.h5ad` etc.
so they don't collide with the spliced-only outputs.

### HTML dashboard

Analogous to spc-align's `align_report.html`, but summarizing downstream
outputs: cells x reads, exonic enrichment, RNA counts per cell, variants per
cell, per-cell VAF distributions.
