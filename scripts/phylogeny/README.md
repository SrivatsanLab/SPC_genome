# phylo — phylogenetic inference for single-cell WGS

Scanpy-style API for building somatic-mutation trees on the demultiplexed
CapWGS/CapGTA AnnData object. Structured to fold into `cellspec.tl.phylo`.

Companion doc: `results/worm6_final/DNA_analysis/phylogeny/PLAN_AMENDMENTS.md`.

## Layout

```
phylo/
  pp/       QC predicates, per-worm somatic panel construction
  tl/       G0 diagnostics, tree builders, DNA-internal validation
  pl/       tree / VAF / covariate plots
  utils/    config + run management
cli.py      python -m phylo <stage> --config <config.yaml>
configs/    YAML configs; one per experiment
```

## Usage

```bash
# from scripts/phylogeny/
pip install -e .

# run one stage
python -m phylo g0      --config configs/default.yaml
python -m phylo panels  --config configs/default.yaml
python -m phylo trees   --config configs/default.yaml

# rerun (skips completed stages unless --force)
python -m phylo all --config configs/default.yaml --force
```

Each run writes to `<output_root>/runs/<name>/`. The resolved config is
persisted as `config.resolved.yaml` in the run dir; every stage marks
completion with `.done_<stage>`.

## Design

- **AnnData in, AnnData out** — matches cellspec conventions
- **Per-worm h5ad panels** — canonical intermediate; small (few MB each),
  self-contained, notebook-inspectable
- **Composable QC** — each filter is a named predicate; the config lists
  which to apply
- **Swappable tree builders** — `phylo.tl.trees.TreeBuilder` ABC; drop-in
  implementations for SVD-bipartition (default), NJ, ...
- **Cheap re-runs** — stages skip if their marker exists; drop the marker
  or pass `--force` to redo
