# CellPhy tree inference pipeline

Genotype-likelihood inference on somatic panels, per
`docs/PANEL_DIAGNOSTICS_AND_ENCODING_PLAN.md` §3.

Reads per-cell `PL` straight from the joint VCF — covered-ref becomes a soft
`0` (~0.3 log10 nudge toward ref at DP=1) instead of a hard `?` or hard `0`.

## Layout

Panel candidates live at:
```
results/worm6_final/DNA_analysis/phylogeny/panels/<panel_tag>/worm_<W>.h5ad
```

Outputs are grouped by panel tag:
```
results/worm6_final/DNA_analysis/phylogeny/inference/<panel_tag>/
├── vcf/worm_<W>/
│   ├── subset.bcf                    positional × sample subset of joint VCF
│   ├── real.vcf.gz                   exact-allele filter of subset
│   ├── covmask.vcf.gz                soft-0 at every DP>0 (allele signal wiped)
│   └── shuffled.vcf.gz               per-variant permutation of (GT,AD,PL) among DP>0 cells
├── cellphy/worm_<W>/{real,covmask,shuffled}/
│   ├── run.raxml.bestTree            ML tree
│   ├── run.raxml.bestModel
│   ├── bs.raxml.bootstraps           100 bootstrap replicates
│   ├── sup.raxml.supportFBP          FBP support on best tree
│   └── sup.raxml.supportTBE          TBE support on best tree
└── summary/worm_<W>/
    ├── bipartitions.tsv
    ├── overlap.tsv                   max-Jaccard vs cellphy:real
    └── permutation_null.tsv          leaf-shuffle null for max-Jaccard
```

## Running

### 0. Build a panel

Write a YAML in `scripts/phylogeny/configs/` (see `tct_kept.yaml` for the
retain-TCT-only-C>A example), then:

```
sbatch --job-name=panel_<tag> \
    --export=ALL,CONFIG=scripts/phylogeny/configs/<tag>.yaml,PANEL_TAG=<tag> \
    scripts/phylogeny/inference/sbatch_build_panel.sh
```

### 1. Prepare VCFs

```
sbatch --job-name=prep_<tag> \
    --export=ALL,PANEL_TAG=<tag> \
    scripts/phylogeny/inference/sbatch_prepare.sh
```

16-task array (one per worm). ~5 min per worm.

### 2. CellPhy

```
sbatch --job-name=cp_<tag> \
    --export=ALL,PANEL_TAG=<tag> \
    scripts/phylogeny/inference/sbatch_cellphy.sh
```

48-task array (16 worms × 3 matrices). `cellphy.sh SEARCH` in the default GL
mode (`-l` opts into ML — we do not want that), then bootstrap and FBP+TBE
support mapping via the bundled raxml-ng. Each task runs ~30–60 min for SEARCH
plus another 20–40 min for bootstrap on 100 replicates.

To restrict to a subset:
```
sbatch --export=ALL,PANEL_TAG=<tag>,WORMS="worm07 worm10",MATRICES="real covmask" \
    --array=0-3 \
    scripts/phylogeny/inference/sbatch_cellphy.sh
```

### 3. Diagnostics

```
micromamba activate cellphy
for W in worm07 worm10 worm20; do
    python scripts/phylogeny/inference/diagnostics.py \
        --worm $W \
        --inference-dir results/worm6_final/DNA_analysis/phylogeny/inference/<tag> \
        --panel-h5ad    results/worm6_final/DNA_analysis/phylogeny/panels/<tag>/worm_${W}.h5ad \
        --out-dir       results/worm6_final/DNA_analysis/phylogeny/inference/<tag>/summary/worm_${W}
done
```

## Design notes

- **GL mode is the default.** `cellphy.sh` reads `PL` from the VCF and builds
  tip likelihood vectors directly. `-l` opts into the ML/GT mode we do *not*
  want (plan §3.2).
- **Model.** `cellphy.sh` defaults to `GT16+FO` (16-state genotype with ML
  frequencies). Do not add `-a` (approximate 10-state) — the plan pins full GT16.
- **covmask semantics.** For every sample × variant: DP>0 → `GT=0/0`,
  `AD=(DP,0)`, `PL=(0,3,20)`. DP=0 → `GT=./.`, `AD=(0,0)`, `PL=(0,0,0)`. This
  preserves the coverage pattern exactly while removing all allele signal —
  any tree recovered from it is coverage confound.
- **shuffled semantics.** Per variant, permute the `(GT,AD,PL)` triple among
  DP>0 cells in this worm. Preserves marginal PL distribution while breaking
  the cell-lineage association. Seed via `SEED` env var (default 0).
- **Cell scope.** Per-worm cells only (plan Q1 answer A). Permutation null in
  `diagnostics.py` shuffles leaf labels within the worm's cell set.
- **SVD-on-covmask control was dropped.** The production
  `SVDBipartitionBuilder` centers per variant in `build_matrix`
  (`svd_kmeans.py:62`), so a matrix with `AD=DP` at all covered cells has
  X ≡ 0 after centering — the builder is architecturally immune to
  coverage-only signal. Recorded as a known property of the tool; no test
  needed.
- **Panel filters.** `phylo.pp.build_panels` supports
  `exclude_mutation_types` (list of substitution types, e.g. `[C>A]`) with
  an override `retain_trinuc_contexts` (list of `anc>der` 3-mer pairs, e.g.
  `[TCT>TAT]`) that keeps matching variants even when their sub-type is in
  the exclusion list. Used by `configs/tct_kept.yaml` to salvage the pole-1
  TCT>TAT hotspot from the C>A exclusion.
- **Missing pieces.** Winner criterion pre-registration (plan §5) is TBD —
  the user will provide before the first real evaluation.
