# wormdemux

Simulation and analysis code for **demultiplexing pooled single-cell WGS into
donor individuals**, specialised to the case of selfing *C. elegans* lines
descended from a single (hypermutator) founder.

The package exists to answer one question before you spend effort on method
tuning: **is there enough signal in the data to separate the individuals at
all?** Everything else follows from that.

---

## Install

```bash
pip install -r requirements.txt      # numpy, scipy, matplotlib, scikit-learn, seaborn
```

No install step is needed for the package itself; scripts add the repository
root to `sys.path`. To use it from elsewhere:

```bash
pip install -e .        # if you add a pyproject.toml, or just:
export PYTHONPATH=/path/to/wormdemux:$PYTHONPATH
```

## Layout

```
wormdemux/
  simulate.py      simulate_pool(): worms, germline + somatic mutations, coverage, error
  cluster.py       site selection, distances (hamming/jaccard/phi/corr), UPGMA, EM
  diagnostics.py   feasibility_report(), carrier-fraction histogram - for REAL data
scripts/
  01_feasibility_mutation_load.py   how much mutation load is needed
  02_somatic_load.py                do somatic mutations break it
  03_correlation_attenuation.py     why within-group correlation looks low
  04_relatedness.py                 how much shared ancestry can be tolerated
  05_metric_comparison.py           does the distance metric matter
  06_real_data_diagnostics.py       run diagnostics on your own .h5ad
```

## Quickstart

```python
from wormdemux import simulate_pool, demultiplex, within_cross
from sklearn.metrics import adjusted_rand_score

d = simulate_pool(n_worms=10, cells_per_worm=60,
                  germline_per_worm=200, somatic_per_worm=800)
r = demultiplex(d["AD"], d["DP"], K=10, metric="phi")
print(adjusted_rand_score(d["truth"], r["labels"]))
print(within_cross(r["similarity"], d["truth"]))
```

On your own data:

```bash
python scripts/06_real_data_diagnostics.py mydata.h5ad --K 20 --plot carrier.png
```

```python
from wormdemux import feasibility_report, demultiplex
feasibility_report(adata.layers["AD"], adata.layers["DP"], K=20)
res = demultiplex(adata.layers["AD"], adata.layers["DP"], K=20,
                  metric="phi", max_sites=20000)
adata.obs["worm"] = res["labels"]
```

`max_sites` bounds memory: the pairwise step densifies `cells x sites`, so
3000 x 20000 float32 is ~240 MB per matrix. Raise it if you have room.

---

## The model

Selfing drives heterozygosity to ~0, so worms are treated as **homozygous**
(genotype 0/1 per site). Consequences that shape the whole approach:

- **Allelic dropout is irrelevant** - there are no heterozygotes to drop. A
  single read at a site is near-deterministic evidence.
- **Binary presence/absence loses nothing.** The usual objection to binarising
  genotypes (it collapses het and hom-alt) does not apply here.
- Sites separate cleanly by **carrier fraction** = (cells with an alt read) /
  (cells covered):

  | class | carrier fraction | meaning |
  |---|---|---|
  | somatic | `< 1/K` | clonal subset of one worm's cells |
  | germline | `~ j/K` | all cells of *j* worms |
  | founder | `~ 1.0` | strain background vs the reference |

  A homozygous variant private to one worm sits at exactly that worm's share of
  the sample. This is why the carrier-fraction histogram is diagnostic: discrete
  peaks at multiples of `1/K` confirm the sites track worm identity and
  independently estimate `K`.

All pairwise quantities are **pairwise-complete** - computed only over sites
covered in both cells. Treating uncovered sites as reference makes distance
track sequencing depth instead of genotype; this is the standard failure mode
and the reason PCA on a zero-filled matrix gives ambiguous structure.

---

## What the experiments show

Numbers below are from the defaults (10 worms x 60 cells, breadth 0.25,
~1.8x depth at covered sites). Rerun with `--worms/--cells/--breadth` to
match your own design; thresholds shift with pool size.

### 01 - Mutation load sets feasibility

| mutations/lineage | regime | shared sites/pair | ARI |
|---|---|---|---|
| 5 | wild-type, ~15 generations | 3 | 0.38 |
| 15 | wild-type, ~50 generations | 9 | 0.58 |
| 50 | mild elevation | 31 | 0.85 |
| 200 | hypermutator | 125 | 1.00 |
| 1000 | strong hypermutator | 624 | 1.00 |

At wild-type *C. elegans* rates (~0.3 mutations/genome/generation) a lineage
accumulates too few private variants to demultiplex over realistic generation
counts. A hypermutator allele moves the problem into the trivially-solvable
regime. **The median shared-informative-sites-per-pair is the number to
check** - roughly >=150 is comfortable, ~40 is degraded, below that the
information is absent.

### 02 - Somatic mutations do not break it (and are partly informative)

| germline/worm | somatic/worm | no filter | germline band only |
|---|---|---|---|
| 300 | 300 | 1.000 | 1.000 |
| 100 | 800 | 1.000 | 1.000 |
| 30 | 1500 | 0.996 | 0.682 |
| 10 | 2000 | 1.000 | 0.280 |
| 3 | 2500 | 0.331 | 0.074 |

Somatic mutations create nested clonal structure *inside* each worm, but they
arise in one worm only, so they carry partial worm-identity information rather
than being pure noise. **Filtering them out can therefore cost accuracy when
germline signal is scarce.** Their real cost is to the *appearance* of the
similarity matrix (see 03), not to clustering accuracy.

### 03 - Why within-group correlation looks low

| scenario | panel | within | cross | gap |
|---|---|---|---|---|
| germline only, no error | all sites | 0.965 | -0.112 | 1.08 |
| + somatic 2x germline | all sites | 0.815 | -0.097 | 0.91 |
| + somatic 6x germline | all sites | 0.623 | -0.078 | 0.70 |
| + somatic 6x germline | germline band | 0.952 | -0.113 | 1.07 |
| + 2% false-positive calls | all sites | 0.523 | -0.065 | 0.59 |
| + 30% unfixed (het) germline | all sites | 0.738 | -0.089 | 0.83 |

Clean homozygous worms should give within-group correlation ~0.95. Observed
values well below that are explained by somatic load, false-positive calls, or
germline variants not yet fixed by selfing (heterozygous, so a coin flip per
read at 1x depth). Restricting to the germline band restores the correlation
almost fully when somatic load is the cause - which makes this a **diagnostic**
for distinguishing the three explanations.

Note the `cross` column: it stays near `-1/(K-1)` in every scenario. That is a
consequence of per-variant centring, not a bug. **The within/cross gap is what
clustering uses, and it is preserved even when the absolute level drops.**

### 04 - Relatedness is the hard limit

| sibling sharing | ARI at K=n_worms | ARI at K=n_worms/2 (clades) |
|---|---|---|
| 0% | 1.000 | 0.179 |
| 50% | 1.000 | 1.000 |
| 90% | 0.842 | 1.000 |
| 98% | 0.602 | 1.000 |

Worms separated by few selfing generations share most of their mutations and
fuse. No method recovers a distinction the data does not contain. The right
response is to report **clades** rather than forcing `K` = number of worms:
the last column shows clade structure remains perfectly recoverable when
individuals do not.

### 05 - Metric choice matters only when coverage is sparse

| scenario | hamming | jaccard | phi | corr |
|---|---|---|---|---|
| clean | 1.000 | 1.000 | 1.000 | 1.000 |
| somatic 4x | 1.000 | 1.000 | 1.000 | 1.000 |
| somatic 4x + 2% FP | 1.000 | 1.000 | 1.000 | 1.000 |
| sparse (breadth 0.05) | 0.108 | 0.088 | **0.882** | 0.299 |

With adequate coverage all metrics tie. At low breadth **phi (binary Pearson)
is dramatically better** - it discounts co-absence against marginal frequency,
so it does not drown in the co-absence that dominates a sparse matrix. Use
`metric="phi"` as the default; it costs nothing when coverage is good.

---

## Choosing K

`K` is a claim about the pool, not a tuning knob. If you know how many
individuals went in, pass it. Otherwise:

- Sweep `K` and check `np.corrcoef(Q)` between cluster consensus genotypes
  (`consensus(...)`); off-diagonal above ~0.95 means two clusters are one
  individual over-split.
- Read `K` off the carrier-fraction histogram: peaks at multiples of `1/K`.
- Expect fewer resolvable clusters than worms if any are closely related (04).

## Caveats - what is NOT modelled

- **Mapping and alignment artefacts.** Recurrent false variants in repeats or
  paralogous regions are not simulated; on real data these are a leading
  alternative explanation for an unexpectedly large informative site set.
- **Amplification chimeras** and reference-mapping bias.
- **Locus-level coverage structure.** Dropout is uniform-random per site here;
  real WGA coverage is spatially correlated along the genome, which reduces the
  effective number of independent observations.
- **Selection** and any non-neutral dynamics in the lineages.
- **Doublets and ambient capsules** are supported in `assign_em` (ambient
  component) but are not part of the worm simulator; they matter more for
  outbred multi-donor pools than for this design.
- The mutation-load thresholds in 01 are for the *default* pool geometry.
  Larger pools need more mutations per lineage to keep pairs distinguishable.

Results are from simulation. They establish the *shape* of the problem and the
relative ordering of methods; absolute accuracy on real data will be lower.
