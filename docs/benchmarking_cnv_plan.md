# CapWGS benchmarking: CNV / aneuploidy calls vs. public WGA methods

## Goal

Per-cell CNV / aneuploidy / large-SV calling on CapWGS HSCs and the public
single-cell WGA datasets, benchmarked against paired bulk WGS. Recreates the
spirit of Gonzalez-Pena et al. (PNAS 2021,
[10.1073/pnas.2024176118](https://www.pnas.org/doi/10.1073/pnas.2024176118))
Fig. 2: per-cell coverage uniformity and spurious CNV calls as a quality
measure of single-cell whole-genome amplification.

Companion to the variant-calling benchmark
(`docs/benchmarking_variant_calling_plan.md`); shares manifest, donor groups,
and depth-matched (25M-read) BAMs.

## Why AneuFinder

* Same tool the user has previously used on CapWGS data — directly comparable.
* Designed for shallow single-cell sequencing (works at <1× coverage); the
  25M-read depth-matched BAMs sit around 0.7–0.9× post-dedup.
* HMM-based segmentation that accommodates the bimodal "amplified vs
  drop-out" bin-count pattern of WGA data, which off-the-shelf bulk CNV
  callers (CNVkit, GATK gCNV) handle poorly.
* Returns per-cell, per-bin copy-state calls, segment tables, and built-in
  quality metrics (`spikiness`, `bhattacharyya`, `entropy`, etc.) for QC.

## Donor groups (reuses variant-calling manifest)

| donor_id | method | cells | bulk(s) | notes |
|---|---|---:|---|---|
| `HSC4` | CapWGS | 50 | `bulkHSC4.bqsr.bam` | this study; HSCs, expected near-diploid |
| `CD34_cord_blood` | PTA | 29 | SRR8438271 | normal HSCs, expected near-diploid |
| `B_lymphocyte` | PTA + SCMDA | 20 | (missing) | no bulk truth — analyzed but excluded from bulk-comparison metrics |
| `BJ_fibroblast` | LIANTI/MDA/MALBAC/DOP-PCR | 44 | SRR5365377, SRR5365378 | primary fibroblasts; bulk may show mild clonal events |

No K562/HepG2-type known-aneuploid samples available, so we **cannot** measure
sensitivity to true CNVs. The comparison is asymmetric:

* **False-positive CNVs** (spurious calls in cells that bulk does not show):
  directly measurable, primary endpoint.
* **Coverage uniformity / bin-variance**: directly measurable, independent of
  truth, key WGA quality metric (analogous to MAPD).
* **Aneuploidy index** (total fraction of genome called non-diploid per cell):
  directly measurable.
* **True-positive sensitivity**: not addressable with current donors; flag as
  a limitation and skip.

## Inputs

Per-cell BAMs already exist:

* Public cells: `data/benchmarking/bams_25M/<donor>/<sample>_25M.bam` (post
  MarkDuplicates, 25M reads).
* HSC4 cells: native depth in `data/HSC4/sc_outputs/`. **Need 25M-downsampled
  variant** for depth-matched comparison — see Open Prereqs.
* Bulks: `data/benchmarking/bulks/<donor>/*.bqsr.bam` and
  `data/HSC4_bulk/bams/bulkHSC4.bqsr.bam` (bulks stay at native depth; bin
  variance is not the metric for bulks).

AneuFinder accepts MarkDup'd BAMs directly — no BQSR or GVCF step needed.

## Pipeline outline

```
                  +-------------------------+
                  |  Per-cell 25M BAMs      |
                  |  + per-donor bulk BAM(s)|
                  +-----------+-------------+
                              |
                              v
              +-------------------------------+
              | scripts/benchmarking/cnv/     |
              | aneufinder_array.sh           |
              |   per cell (and per bulk)     |
              |   - binReads (variable+fixed) |
              |   - GC + mappability correct  |
              |   - findCNVs (HMM, edivisive) |
              |   - karyotypeMeasures         |
              +---------------+---------------+
                              |
                              v
              +-------------------------------+
              | per-sample AneuFinder RDS:    |
              | data/benchmarking/cnv/        |
              |   <donor>/<sample>.rds        |
              | + segment TSV + QC TSV        |
              +---------------+---------------+
                              |
                              v
              +-------------------------------+
              | scripts/benchmarking/cnv/     |
              | compile_cnv_metrics.py        |
              |   join with donor manifest    |
              |   per-cell FP CNV count vs    |
              |   matched bulk                |
              +---------------+---------------+
                              |
                              v
              +-------------------------------+
              | notebooks/benchmarking_cnv.   |
              | ipynb                         |
              |  - bin-variance / spikiness   |
              |  - genome-wide CN heatmaps    |
              |  - per-method boxplots        |
              +-------------------------------+
```

## Scripts to write

1. **`scripts/benchmarking/cnv/build_cnv_task_list.py`** — emits a TSV of
   `(donor, sample, bam_path, is_bulk)` rows over all cells + bulks, similar
   to `build_donor_manifest.py`. Reuses
   `data/benchmarking/manifests/donor_manifest.tsv` plus the 25M-BAM paths.
2. **`scripts/benchmarking/cnv/aneufinder_array.sh`** — SLURM array, one task
   per (donor, sample). Runs an R driver `aneufinder_one_sample.R` that:
   * `binReads()` at variable-width bins (default ~1Mb effective; let
     AneuFinder choose given read count) and **also** at 1Mb fixed bins for
     visualization consistency across cells.
   * `correctGC()` + `correctMappability()` against hg38 blacklists.
   * `findCNVs(method=c('edivisive','HMM'))` — both, write both calls.
   * `karyotypeMeasures()` for per-cell QC.
   * Save `.rds` + flattened segment TSV + a one-line QC TSV.
3. **`scripts/benchmarking/cnv/aneufinder_one_sample.R`** — the actual R
   driver invoked by the array. Single sample in, all outputs out.
4. **`scripts/benchmarking/cnv/compile_cnv_metrics.py`** — assembles per-cell
   QC + segment summaries into one long-form TSV, joined with the donor
   manifest. Computes per-cell FP CNV count = (cell segments off-diploid) ∖
   (bulk segments off-diploid), with a tolerance on breakpoint matching
   (e.g. ±2 bins).
5. **`notebooks/benchmarking_cnv.ipynb`** — figures:
   * Per-cell bin spikiness / `bhattacharyya` boxplots by method (analog of
     PTA Fig 2: coverage uniformity).
   * Per-cell number of non-diploid segments by method.
   * Per-cell aneuploidy index (`fraction-of-genome-not-2N`) by method.
   * Cell-vs-bulk side-by-side CN heatmaps for representative cells per
     method.
   * Optional: clustering by AneuFinder QC features (`clusterByQuality`).

Convention: outputs go under `data/benchmarking/cnv/` and follow the same
per-donor subdir layout used by the variant pipeline.

## R environment

Use the Hutch `fhR` module — AneuFinder 1.32.0 ships pre-installed there.

```
ml fhR        # R 4.4.1, AneuFinder 1.32.0
Rscript scripts/benchmarking/cnv/aneufinder_one_sample.R ...
```

The SLURM wrapper just does `module load fhR` at the top; no conda env
needed.

## Open prerequisite: HSC4 25M BAMs

HSC4 cells are native depth (170M – 1.7B reads). The 25M downsampled BAMs
exist for the public donors but **not** for HSC4. Same situation as the
variant-calling plan. Options:

1. Run `scripts/benchmarking/downsample_array.sh` for HSC4 at 25M, MarkDup,
   then feed those BAMs to the AneuFinder array. (Preferred.)
2. Run AneuFinder on HSC4 at native depth and note that lower-depth public
   cells will look noisier purely from binomial sampling. Bin-variance
   comparison gets confounded in that case, so this option weakens the
   primary metric.

Pick option 1 before launching the array.

## Open prerequisite: bin size + mappability tracks

AneuFinder needs an hg38 mappability BigWig for `correctMappability()`. Lab
may already have one; if not, generate or download (GEM / ENCODE 100mer
mappability for hg38). Decide bin size at run time: aim for ~200 reads/bin
median, which at 25M reads is roughly 1Mb variable bins.

## Metrics (primary endpoints)

| metric | unit | source | what it tests |
|---|---|---|---|
| `spikiness` | per cell | AneuFinder QC | bin-to-bin read-count instability — analog of MAPD |
| `bhattacharyya` | per cell | AneuFinder QC | separation of euploid bin-count peaks — segmentation tractability |
| `n_segments_nondiploid` | per cell | segment TSV | spurious large CNV burden |
| `aneuploidy_index` | per cell | derived from segments | fraction of genome called non-2N (excluding bulk-supported regions) |
| `cell_vs_bulk_FP_count` | per cell | segment TSV ∖ bulk segments | direct FP rate vs paired bulk |
| `cell_vs_bulk_overlap` | per cell, % | segment intersection | analog of CN concordance |

Aggregation: median ± IQR per WGA method (CapWGS, PTA, LIANTI, MDA-pooled
via the `wga_method_grouping` memory rule, MALBAC, DOP-PCR).

## Limitations to document up front

* No K562/HepG2/etc. known-aneuploid samples → cannot estimate sensitivity
  to true CNVs; only FP and uniformity.
* B_lymphocyte has no bulk → excluded from cell-vs-bulk FP analysis;
  retained for inter-method uniformity comparison.
* BJ_fib bulk is from a different ATCC vial run than the LIANTI bulks of
  origin — small clonal differences possible; treat marginal bulk events
  with caution.
