# Genetic demultiplexing of pooled single-cell WGS — analysis summary

Prepared for review by another agent. This documents the problem, the methods
tried, what was empirically established, what was **retracted**, and what
remains open. Numeric claims are labelled by their source: **[data]** = the
user's real data, **[sim]** = simulation run during this analysis.

---

## 1. The problem

Single-cell whole-genome amplification and sequencing of ~2,400 capsules from a
pool of ~16–20 *C. elegans* individuals. All individuals are selfing descendants
of one genetically engineered founder carrying a **pole-1 P270R hypermutator**
allele. Goal: assign each capsule to its worm of origin, so that somatic
mutations can subsequently be called per worm.

**Data.** AnnData object, `.layers["AD"]` (alt depth) and `.layers["DP"]` (total
depth), 963,336 called sites (bcftools joint calling), ~2,404 capsules.
`.var["bulk_vaf"]` is a precomputed pseudobulk VAF using hypergeometric
downsampling (`cellspec.tl.compute_bulk_vaf`, `target_dp=1000`).

**Properties that drive every design decision:**

| property | value | consequence |
|---|---|---|
| depth at covered sites | **~1×** (86% of covered sites at DP=1) **[data]** | per-cell VAF is 0 or 1; marginal per-variant statistics carry no structural information |
| coverage breadth | spans ~5%–80% of the genome per cell **[data]** | any method requiring a minimum per-cell site count discards a large fraction of capsules |
| zygosity | selfing ⇒ worms are effectively **homozygous** | no heterozygotes; allelic dropout is irrelevant; a somatic mutation (het, VAF≈0.5) is distinguishable from germline (hom, VAF≈1) *in aggregate* |
| relatedness | all worms share one founder | large common variant background; closely related worms may be unresolvable |
| contamination | **pervasive and graded** | every capsule contains ambient ("soup") DNA; no discrete empty-droplet population |
| somatic load | ~495k sites below the germline band vs ~225k in it **[data]** | somatic variants numerically dominate the panel |

---

## 2. Key quantitative facts established

### 2.1 Read-level bulk VAF is contamination-invariant

For a homozygous variant private to worm *k*, pooled bulk VAF equals **worm
*k*'s share of total reads**, independent of contamination α. Contamination
redistributes alt reads among capsules without changing the pool-wide alt
fraction.

Verified **[sim]**: correlation between worm read share and module VAF =
**0.9997** at median α = 0.31, ratios 1.00–1.08 across a 25× range of worm
sizes.

**Consequences:**
- Worm-private variants sit at bulk VAF ≈ (that worm's cell share), **not** at a
  common 1/K unless worms contributed equally. With cell counts spanning
  539→7 **[data]**, expected read shares span ~0.22→0.003.
- The germline VAF spectrum should be **multimodal**, one peak per worm.
- Carrier fraction (cells with ≥1 alt read ÷ cells covered) is **NOT**
  equivalent and should not be substituted: it is inflated by depth ×
  contamination (1.00× at DP=1, 1.30× at 1.5×, 1.99× at 3× **[sim]**).

### 2.2 Marginal feature selection cannot work at 1× coverage

A variant carried by 1/K of cells at VAF≈1 is **marginally indistinguishable**
from a variant occurring at rate 1/K in every cell — both give Bernoulli(1/K)
per-cell observations. Verified **[sim]**: median overdispersion (chi-square vs
binomial null) = **1.078 for germline vs 1.080 for artifact**, null = 1.0.

This is why scRNA-seq highly-variable-gene selection does not transfer. All
structural information lives in **co-occurrence across variants**.

### 2.3 Per-variant centering is what makes the method work

`colmean` = unweighted mean of per-cell VAF over covered cells. Numerically
≈ bulk VAF at this depth (ratio 1.000 **[sim]**), but **not** the carrier
fraction (0.083 vs 0.094 **[sim]**). It is the correct center because it makes
each column sum to exactly zero over observed entries.

Centering removes the shared founder background (a variant carried by all worms
has colmean ≈ 1, so its centered value ≈ 0 for everyone). This is why
per-variant-centered Pearson beat the alternatives on real data:

| metric | within-worm | cross-worm | gap **[sim]** |
|---|---|---|---|
| Hamming | 0.98 | **+0.59** | 0.35 |
| phi (Pearson on binary) | 0.94 | −0.16 | 1.10 |
| **Pearson on centered AD/DP** | 0.99 | −0.54 | **1.54** |

Hamming leaves founder/clade variants counted as agreements → compressed range →
renders as a smooth gradient, which is exactly what the user observed **[data]**.
phi centers each cell by its own marginal, never by site frequency, so shared
ancestry persists; its marginals are also unstable in low-coverage cells.

Cross-worm correlation is **negative** by construction: per-variant centering
forces `r_cross ≈ −r_within/(K−1)`, confirmed at K=4 (−0.20) and K=8 (−0.09)
**[sim]**.

### 2.4 Why observed within-worm correlation is ~0.4 rather than ~0.95

Attributable to somatic load, false-positive calls, and unfixed heterozygous
germline variants **[sim]**:

| scenario | within | cross | gap |
|---|---|---|---|
| germline only, no error | 0.965 | −0.112 | 1.08 |
| + somatic 6× germline | 0.623 | −0.078 | 0.70 |
| + somatic 6× germline, germline band only | 0.952 | −0.113 | 1.07 |
| + 2% false-positive calls | 0.523 | −0.065 | 0.59 |
| + 30% unfixed (het) germline | 0.738 | −0.089 | 0.83 |

The **gap** is what clustering uses; the absolute level is not the figure of
merit.

---

## 3. Final method

### Pipeline

1. **Site panel** — filter on `bulk_vaf` window (three windows swept:
   germline_peak 0.018–0.06, high_vaf 0.018–1.0, all_variants 0–1.0), plus
   `bulk_dp > 400`. `MIN_DP_VAR = 1` (every read counts).
2. **Centered matrix** — `X` = cells × variants, per-cell VAF minus per-variant
   mean over covered cells; uncovered entries imputed at that mean, hence
   exactly 0 after centering, hence dropping out of the factorization.
3. **Provisional decomposition** — `TruncatedSVD(n_components=ncomp)` on `X`;
   variant loadings `components_.T` are L2-normalized; K-means into K modules;
   `d_all` = distance from each variant to its assigned centroid.
4. **Core selection** — keep the `N_CORE` (12,000) variants with lowest `d_all`;
   refit SVD + K-means on that subset.
5. **Cell scoring** — for each module, `S[i,j]` = (alt reads cell *i* has at
   module *j*'s variants) ÷ (total reads there) = depth-weighted pooled VAF.
   `W = S / S.sum(1)`. No cell filtering anywhere; a low-coverage cell simply
   contributes fewer reads.
6. **Calling** — `donor = argmax(W)`; `purity = max(W)`; `alpha_hap = 1−purity`;
   `ratio = second/first`; `eff_modules = 1/Σw²`.

### Why the score vector is interpretable

For a cell from worm *k* with contamination α: `S[i,k] ≈ (1−α) + α/K`,
`S[i,j≠k] ≈ α/K`, rows sum to ≈1 (verified **[sim]**). Hence `1 − purity`
estimates α directly (α=0.4 → purity 0.620).

### Two SVD directions (a persistent source of confusion)

| | object | shape | used for |
|---|---|---|---|
| `svd.components_.T` | **variants** | n_var × ncomp | K-means → modules; `d_own` |
| `svd.fit_transform(X)` | **cells** | n_cells × ncomp | UMAP embedding |

Cells are **never clustered**. `d` is always a **K-means** centroid distance in
the L2-normalized variant-loading space.

---

## 4. Metric pitfalls discovered (important for review)

Several metrics were found to be **mechanically confounded**. Each was a real
bug that changed conclusions.

| metric | confound | status |
|---|---|---|
| `median d` (distance to centroid) | rises with `ncomp` (0.177 at 11 comp → 0.909 at 60, with clustering accuracy flat at ARI 1.000 **[sim]**); falls with K | compare **down columns** at fixed ncomp only |
| `excess_purity` = (p−1/K)/(1−1/K) | rises with K for free as the 1/K floor shrinks (~2% over K=15→26) | never use alone to choose K |
| `size_cv`, `smallest` | computed as `sz = np.bincount(donor); sz = sz[sz>0]` — **empty modules are invisible**, so over-specifying K costs nothing **[sim]** | **fixed**: added `occupancy()` reporting modules with ≥1 cell / ≥5 cells / ≥1% of cells, plus `size_cv_full` including empties |
| `occupied_ge1pct` (core arm) | denominator was assigned-cell count, which shrinks as K rises (more background) → the 1% bar drops → more modules clear it | **fixed**: denominator is now total capsules in both arms |
| core selection at every K in the sweep | each K evaluated on a panel optimized for that K → over-splitting masked | **fixed**: added a **no-selection arm** (`noselect_*` columns) which should choose K |
| held-out silhouette | circular if computed on the same variants used to assign | use A/B variant split |

### The scree cliff is at **K−1**, not K

Per-variant centering forces the K group means to sum to zero, dropping the rank
to K−1. Verified at K = 6, 8, 10, 12, and robust to 25% coverage, 1× depth,
somatic structure, and unequal cell counts **[sim]**.

**But the cliff is undetectable with unequal cell counts**: gap-at-K ÷
flat-region gap falls from ~150 (ideal) to **~1.1** when cell counts are
proportional to the real data **[sim]**. The user's real scree indeed shows no
cliff, with the largest gap at component 2 **[data]** — itself evidence of the
size imbalance.

---

## 5. Core-selection debate (unresolved but characterized)

The user challenged core selection as circular. Testing produced a mixed
picture:

**Against the concern:** at the correct K, core selection *purged* artifacts
entirely (0.0% vs 9.6% baseline **[sim]**), because artifact variants fit no
worm module and get high `d`.

**Supporting the concern:** at K=14 and 16 (true K=12), an entire module became
**100% artifact variants** **[sim]** — artifacts do claim their own module when
K is over-specified. This is quarantine rather than contamination (cell ARI
stayed 0.81), but cells assigned to that module are misassigned.

**Against the method:** at the correct K, core selection was *harmful* —
cell ARI 0.856 with selection vs **0.952 without** **[sim]**.

**For the method:** at over-specified K it was strongly protective — 0.42
(no selection) → **0.82** (core) at K=20 **[sim]**.

**Circularity appears not to be doing practical harm:** selecting the core from
a random half of the cells vs all cells gave 0.879/0.856, 0.899/0.810,
0.801/0.824, with ~70% panel overlap **[sim]**.

**Alternatives tested and rejected [sim]:**

| selection | % artifact retained | cell ARI |
|---|---|---|
| none | 9.6% | **0.952** |
| core (lowest d) | 0.0% | 0.856 |
| most-covered | 8.9% | 0.942 |
| communality (‖loading‖ at ncomp=50) | **20.0%** | **0.384** |
| marginal overdispersion | 9.8% | 0.941 |

Communality — proposed as a K-free structural criterion — *enriched* artifacts
2× because artifacts are highly structured. It should not be used.

**Current resolution:** run both arms. Choose K from the no-selection arm
(where over-specification is penalized); use core selection for the final fit
(where it confers robustness). Written `.obs`/`.var` come from the core arm.

---

## 6. Capsule classification

### The soup/doublet confusion (diagnosed)

Original rule: `soup = purity < 2/K`; `doublet = ~soup & (ratio > 0.6)`.

The user observed few soup calls, many doublets, and two distinct UMAP clusters
of "doublets" with low purity **[data]**.

**The two rules flag identical cell sets.** Geometrically, weights sum to 1 and
the other K−2 modules cannot exceed the second, so `ratio ≥ (1−p)/(p(K−1))`,
meaning `ratio ≤ 0.6` is only possible when `purity ≥ 1/(1+0.6(K−1))` = 0.0806
at K=20 — just below `2/K` = 0.10. Verified: 0.00% of capsules fall in the
region where the rules could differ **[sim]**. The rules differ only in
**labeling**, splitting the flagged set by purity.

**That split may be a coverage artifact.** With the *same* soup capsules, only
varying weight-estimate noise **[sim]**:

| soup weight noise | called "soup" | called "doublet" |
|---|---|---|
| near-uniform | 100% | 0% |
| noisy | 68% | 30% |
| very noisy | 7% | **82%** |

Noisier weights inflate the max above 2/K → fails soup test → falls through to
doublet. Since weight noise is set by covered sites per cell, **the soup/doublet
split may track coverage rather than capsule content**.

**Soup purity does not equal 1/K in practice** — it is inflated by the noise
floor: median 0.286 at 6 covered sites/module, 0.157 at 25, 0.106 at 100, 0.092
at 200 **[sim]**, against a theoretical 0.0625 at K=16.

### Recommended statistic: effective modules

`eff = 1/Σw²` (inverse Simpson) reads out how many genomes contributed **[sim]**:

| capsule | purity | ratio | eff |
|---|---|---|---|
| clean singlet | 1.000 | 0.000 | 1.00 |
| singlet, α=0.85 | 0.193 | 0.221 | 14.0 |
| clean doublet | 0.500 | 1.000 | **2.00** |
| doublet, α=0.6 | 0.230 | 1.000 | 8.20 |
| pure soup | 0.050 | 1.000 | **20.0** |

Unlike purity, its scale is fixed by the model rather than by the noise floor.
Note `ratio` cannot distinguish soup from doublets (both ≈1) and *misses*
heavily contaminated singlets (α=0.85 → ratio 0.221).

**Current user setting [data]:** `background = (ratio > 0.6) | (purity < 0.2)`.
At K=20, purity 0.2 corresponds to α ≈ 0.84. Note 0.2 is 4/K at K=20 but 3.2/K
at K=16 — scale it if comparing across K.

### A retracted claim, and a real bug found

The user reported cells with purity=1.0 being called background under a
ratio-only rule. This is **mathematically impossible** (`second ≤ 1−purity`, so
`ratio > 0.6` requires `purity < 0.625`); verified across 200,000 simulated
weight vectors, max purity among `ratio>0.6` cells = 0.567 **[sim]**. The user
confirmed it was a bug in their code, since fixed.

**Genuine edge cases in `_score` that remain unguarded:**
- A capsule covering **zero** core sites → `S` row all zero → `W` all zero →
  `purity = 0`, `argmax` returns module 0 arbitrarily.
- A capsule covering sites in **only one** module → `W = (1,0,…)` → **purity =
  1.0** on almost no evidence. Purity measures the *shape* of evidence, never
  its *quantity*; `n_core_sites` is the necessary companion.

---

## 7. Open questions

1. **True K.** The corrected `median d` sweep gave a **plateau across K=14–18**
   rather than a sharp optimum **[data]**. Scree suggested K≈15–19. User
   selected K=20 for a later run. The occupancy fix has not yet been rerun.

2. **Rising occupancy at high K [data].** `occupied` and `ge5` hug the y=K line
   to K=22, and `ge1pct` for germline_peak rises through K=26. Two candidate
   explanations, not yet distinguished:
   - more worms than assumed;
   - **sub-worm somatic clonal structure** — plausible under a hypermutator, and
     coherent enough to form modules holding tens of cells.

   **Proposed discriminator (not yet run):** within-module VAF. Germline
   (homozygous) modules sit at ≈(1−α); somatic (het) sub-clones at ≈half that.
   The ratio stays ~0.5 across α **[sim]**, so normalizing by the maximum removes
   α. Complementary test: VAF-vs-cell-share — germline modules on the identity
   line, somatic sub-clones well below it.

3. **Relatedness ceiling.** Sibling pairs sharing >90% of mutations fuse
   (ARI 0.842 at 90%, 0.602 at 98% **[sim]**), while clade-level structure
   remains perfectly recoverable. If K resolves below the worm count, reporting
   clades is the honest output.

4. **`N_CORE` = 12,000 was fixed by hand** and never swept. It interacts with
   both K and the VAF window (the germline_peak panel is far smaller than
   all_variants, so the core is a much harsher distillation in one than the
   other).

5. **`ncomp` vs K.** Best combination reported as `ncomp=16, K=20` **[data]**,
   but K worms span K−1 = 19 dimensions after centering, so 16 components cannot
   fully separate 20 worms. The extended grid (`ncomp` to 30) now covers
   `ncomp ≥ K−1` across the K range; needs rechecking after the rerun.

---

## 8. Code inventory

| file | purpose |
|---|---|
| `sweep_one.py` | one (window, ncomp, K) combination; two arms; writes `metrics.csv`, `obs.csv`, `var.csv`, `weights.npz`, `scree.png/csv`, `distances.png`, `params.json` |
| `make_params.sh` | writes `params.tsv` (3 windows × 11 ncomp × 9 K = 297 rows) |
| `sweep_array.sh` | SLURM array task; pulls its row by `SLURM_ARRAY_TASK_ID`; skips completed runs; pins BLAS threads |
| `submit_sweep.sh` | generates the grid and submits; `--per-window` sizes memory per window |
| `collect_sweep.py` | concatenates all `metrics.csv`; ranks by `noselect_combined`; prints occupancy saturation |
| `final_cell.py` / `cells678.py` | final assignment and diagnostics |
| `module_sizes.py` | sites-per-module, sites-vs-cells, VAF-vs-cell-share validation |
| `splitcheck2.py` | permutation test for fused (two-worm) modules; validated z=38.5 fused vs z=2.1 single **[sim]** |
| `viz.py` | sweep result visualization |

### Known caveats in the simulator used throughout

- Somatic mutations were modelled as **homozygous** in early simulations; this
  was wrong (they arise on one chromosome, VAF≈0.5) and made separation look
  artificially clean. Corrected in later tests.
- Donors drawn as independent HWE genotypes in early sims (maximally distinct);
  relatedness added later.
- Dropout is uniform-random; real WGA coverage is spatially correlated.
- Ambient contamination injected at uniform pool frequency — the cleanest
  possible version.

**All simulation-derived thresholds should be treated as order-of-magnitude
guidance, not calibrated values.**

---

## 9. Summary of retractions

Claims made during the analysis and later withdrawn on evidence:

1. "Carrier fraction ≈ colmean at 1× coverage" — **wrong**, 0.083 vs 0.094, and
   they diverge sharply with depth.
2. "Worm-private variants sit at 1/K" — **only if worms contribute equally**;
   correct statement is (that worm's read share).
3. "L2 normalization makes clustering key on pattern rather than coverage" —
   **overstated**; normalization changed nothing across every regime tested
   (ARI identical to 3 decimals). Its real value is making `d` interpretable on
   a fixed scale.
4. "The scree cliff lands at K" — **it is at K−1**.
5. "Core selection will preferentially retain artifacts" — **the opposite** at
   the correct K.
6. "Communality is a good K-free selection criterion" — **enriched artifacts 2×**.
7. "Reselecting the core at each K fixes circularity" — it does not; it flatters
   every K with its own optimized panel.
8. Per-module bulk VAF as a module-quality diagnostic — **confounded** when the
   VAF filter window is narrow, since every module then reports the window
   center.
