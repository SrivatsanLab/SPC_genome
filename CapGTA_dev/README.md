# CapGTA_dev — local-baseline expression model

Sandbox for experimenting with a per-locus, per-cell expression model
calibrated against spliced-read ground truth. One-off for the current
`worm6_final` batch; if it works, migrate the winning approach back into
the main pipeline.

## Problem

The current pipeline has two count paths, both flawed:

| path | script | signal | flaw |
|---|---|---|---|
| spliced-only | `filter_spliced_reads_array.sh` → `create_rna_count_matrix.sh` | reads with `N` in CIGAR | zero for 24,498 junction-less genes; Poisson-thinned even where non-zero |
| noFeatures excess | `assemble_expression_matrix.py` | `R_obs - λ_cell · L_gene` | `λ` is global per cell, ignores MDA's per-locus bias; no calibration |

MDA amplification is highly locus-variable (orders of magnitude across
neighboring windows), so a global per-cell `λ` mis-scales the expectation
at every gene. The excess model still produces a count for junction-less
genes, but we've never checked whether that number is meaningful.

## Reframe

Spliced counts, upsampled by junction-spanning probability
`p_junc(g) ≈ n_junc·(R-1)/L`, are an approximate RNA ground truth for
junction-containing genes. Use them to calibrate a locally-baselined
excess model, then apply the calibrated model to junction-less genes.

## Model

Hierarchical / generative form (all cells × all genes):

```
R_obs[c,g]   ~ Poisson( λ_gDNA[c,g] · L_g  +  μ_RNA[c,g] )
S_obs[c,g]   ~ Binomial( μ_RNA[c,g], p_junc(g) )    # only genes w/ junctions
λ_gDNA[c,g]  = local baseline (per-cell, per-locus)
μ_RNA[c,g]   = size_factor[c] · rate_g · covariate correction
```

Spliced observations pull `μ_RNA` toward the truth for junction genes;
shared per-cell / per-gene-covariate parameters carry that constraint to
the junction-less genes.

For a "one-off, single batch" scope we don't need full Bayes. Two
tractable variants of the same idea to sweep:

- **NB-GLM** with `log(λ_local·L_g)` as offset, size factor per cell, gene
  random effect, covariates. Gives calibrated point estimates + rough SEs.
- **GBDT residual correction** — predict `upsampled_S - excess_estimate`
  from gene/cell covariates. Fast, flexible; no principled uncertainty.

## Local baseline options (sweep dimension)

- **intron**: coverage over per-gene merged intron intervals. Best local
  match for multi-exon genes; unavailable for junction-less genes → fall
  back to flanking.
- **flanking**: ±5–10 kb windows around gene body, subtracting any
  overlap with other annotated features. Available for all genes;
  chromatin around expressed genes may differ from within-gene gDNA.
- **hybrid**: intron where available, flanking elsewhere.
- **global** (current): keep as the baseline-of-baselines to beat.

## Sweep dimensions (as of sweep v2)

| dim | values |
|---|---|
| baseline | **global** (winner) / flanking with `floor_alpha ∈ {0, 0.25, 0.5, 1}` |
| model | excess (implemented) / NB-GLM (stub) / GBDT (stub) |
| features | minimal (length, `n_junc`) / full (+ `spliced_libsize`, `expected_gdna`) |
| RNA signals | spliced (kept) + polyA/T (added, ~15% more evidence) + SL (dropped, too sparse) |

**Dropped from earlier plans:** intron/hybrid baselines (pre-mRNA
contamination), read-position-adjusted `p_junc` (would require bulk
statistics — see sweep_v1_notes.md), flanking-window sweep (fixed at 5 kb).

## Validation

No external bulk RNA-seq — validate internally:

1. **Held-out junction genes** (20%, stratified by length / expression /
   biotype): compare `μ_RNA` to upsampled spliced count. Metrics: Pearson
   on `log(1+·)`, Poisson deviance.
2. **Propensity-weighted extrapolation error**: reweight the junction
   held-out set toward the junction-less covariate distribution. This is
   the number that actually matters.
3. **Silent negative controls**: junction-less genes expected to be off
   in this tissue/sex context — model must give near-zero.
4. **Marker recovery** on the final expression matrix: does adding
   junction-less genes surface new markers (tRNAs, ncRNAs, single-exon
   TFs) that the spliced-only matrix couldn't see?
5. **Distribution QC**: per-cell RNA-count distributions from each
   config — plot against the current excess and spliced-only outputs.

## Layout

```
CapGTA_dev/
├── README.md                              # this file
├── scripts/
│   ├── build_local_baselines.py           # GTF → intron/flanking SAF → featureCounts
│   ├── build_calibration_table.py         # merge counts + features → calibration.h5ad
│   ├── models.py                          # Excess, NBGLM, GBDT — fit(adata)/predict(adata) → (C,G) float32
│   └── run_sweep.py                       # grid over baseline × model × features
├── notebooks/
│   └── 01_explore_and_sweep.ipynb         # baseline distributions, sweep results, count-distribution plots
└── results/
    └── worm6_final/
        ├── calibration.h5ad                # merged input (exon/spliced/intron/flank as layers)
        └── sweep/<tag>/                    # per-config predictions.h5ad + metrics.json + config.json
```

## Scope

- Runs only against `results/worm6_final`. No portability requirements.
- Uncertainty on the RNA estimate is nice-to-have, not required.
- No independent bulk RNA-seq available.
- If a winning combination emerges, port `build_local_baselines.py` + the
  chosen model into `scripts/CapGTA/` and wire it into the main pipeline.
