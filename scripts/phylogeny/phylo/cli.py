"""CLI: ``python -m phylo <stage> --config <path> [--force] [--name NAME]``.

Stages: ``g0``, ``panels``, ``trees``, ``validate``, ``covariates``, ``bwm``,
``all``. Each stage checks its ``.done_<stage>`` marker in the run dir; if
present and ``--force`` is not set, the stage is skipped. Trees / validate /
covariates / bwm are placeholders until implemented.
"""

from __future__ import annotations

import argparse
import sys

import os

import anndata as ad
import pandas as pd

from .pp.panels import build_panels
from .pp.qc import apply_qc
from .tl.diagnostics import run_g0
from .tl.trees import SVDBipartitionBuilder, build_worm_tree
from .utils.config import load_config
from .utils.run import STAGES, RunDir, clear_stage_marker, stage_done, write_stage_marker

ALL_STAGES = STAGES


def _load_adata(cfg: dict) -> ad.AnnData:
    path = cfg["inputs"]["adata"]
    print(f"[cli] loading {path}")
    return ad.read_h5ad(path)


def _stage_g0(cfg: dict, run: RunDir) -> None:
    adata = _load_adata(cfg)
    apply_qc(adata, cfg["qc"]["predicates"], inplace=True)
    run_g0(adata, cfg, run)


def _stage_panels(cfg: dict, run: RunDir) -> None:
    adata = _load_adata(cfg)
    apply_qc(adata, cfg["qc"]["predicates"], inplace=True)
    build_panels(adata, cfg, run)


def _stage_trees(cfg: dict, run: RunDir) -> None:
    panels_dir = run.path / "panels"
    if not panels_dir.exists():
        raise FileNotFoundError(f"No panels/ dir under {run.path}; run 'panels' stage first.")
    panel_files = sorted(panels_dir.glob("worm_*.h5ad"))
    if not panel_files:
        raise FileNotFoundError(f"No worm_*.h5ad panels in {panels_dir}")

    trees_dir = run.subdir("trees")
    tcfg = cfg["trees"]
    method = tcfg.get("method", "svd_bipart")
    if method != "svd_bipart":
        raise NotImplementedError(f"tree method '{method}' not implemented; only 'svd_bipart' for now.")

    # nproc: config > SLURM_CPUS_PER_TASK > 1
    nproc = tcfg.get("mquad_nproc")
    if nproc in (None, "auto"):
        nproc = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    builder = SVDBipartitionBuilder(
        ncomp=int(tcfg.get("ncomp", 2)),
        n_core=(int(tcfg["n_core"]) if tcfg.get("n_core") else None),
        min_dp_var=int(tcfg.get("min_dp_var", 1)),
        min_variants_to_split=int(tcfg.get("min_variants_to_split", 20)),
        seed=int(tcfg.get("seed", 0)),
        core_selection=str(tcfg.get("core_selection", "none")),
        mquad_delta_bic_threshold=float(tcfg.get("mquad_delta_bic_threshold", 10.0)),
        min_carriers_at_node=int(tcfg.get("min_carriers_at_node", 2)),
        mquad_nproc=int(nproc),
    )
    print(f"[trees] builder: {tcfg.get('method')} | core_selection={builder.core_selection} | "
          f"ΔBIC≥{builder.mquad_delta_bic_threshold} | mquad_nproc={builder.mquad_nproc}")

    summaries = []
    all_assignments = []
    for p in panel_files:
        worm = p.stem.replace("worm_", "")
        adata = ad.read_h5ad(p)
        print(f"[trees] {worm}: {adata.n_obs} cells × {adata.n_vars} panel variants")
        tree = build_worm_tree(
            adata,
            builder,
            min_clade_size=int(tcfg.get("min_clade_size", 10)),
            max_depth=int(tcfg.get("max_depth", 6)),
            worm=worm,
            params={"builder": builder.__dict__, **{k: v for k, v in tcfg.items()}},
        )
        wdir = trees_dir / f"worm_{worm}"
        wdir.mkdir(parents=True, exist_ok=True)
        tree.to_json(wdir / "tree.json")
        pd.DataFrame(tree.cell_assignments(), columns=["cell_id", "leaf_path"]).to_csv(
            wdir / "cell_assignments.csv", index=False
        )
        for cell_id, leaf in tree.cell_assignments():
            all_assignments.append({"worm": worm, "cell_id": cell_id, "leaf_path": leaf})
        stats = tree.stats()
        summaries.append(stats)
        print(f"[trees] {worm}: {stats['n_splits']} splits, {stats['n_leaves']} leaves, max_depth={stats['max_depth']}")

    pd.DataFrame(summaries).to_csv(trees_dir / "summary.csv", index=False)
    pd.DataFrame(all_assignments).to_csv(trees_dir / "cell_assignments.csv", index=False)
    print(f"[trees] {len(summaries)} trees written → {trees_dir}")

    from .pl.trees import plot_all_trees
    from .tl.bipartitions import extract_bipartitions
    extract_bipartitions(run.path, write=True)
    plot_all_trees(run.path)


def _stage_clade_matching(cfg: dict, run: RunDir) -> None:
    from .tl.clade_matching import run as run_matching
    ccfg = cfg.get("clade_matching", {})
    if not ccfg.get("rna_path"):
        raise ValueError("clade_matching.rna_path is required")
    run_matching(
        run.path, ccfg["rna_path"],
        similarity_methods=tuple(ccfg.get("similarity_methods", ("pearson", "spearman", "cosine"))),
        n_null=int(ccfg.get("n_null", 500)),
        min_total_counts=int(ccfg.get("min_total_counts", 20)),
        counts_layer=str(ccfg.get("counts_layer", "counts")),
        min_worms_per_depth=int(ccfg.get("min_worms_per_depth", 2)),
        residualize_by_worm=bool(ccfg.get("residualize_by_worm", True)),
        trustworthy_only=bool(ccfg.get("trustworthy_only", False)),
        seed=int(ccfg.get("seed", 0)),
    )


def _stage_todo(name: str) -> None:
    print(f"[cli] stage '{name}' not implemented yet; skipping.")


DISPATCH = {
    "g0": _stage_g0,
    "panels": _stage_panels,
    "trees": _stage_trees,
    "validate": lambda cfg, run: _stage_todo("validate"),
    "covariates": lambda cfg, run: _stage_todo("covariates"),
    "bwm": lambda cfg, run: _stage_todo("bwm"),
    "clade_matching": _stage_clade_matching,
}


def _run_stage(stage: str, cfg: dict, run: RunDir, force: bool) -> None:
    if force:
        clear_stage_marker(run, stage)
    if stage_done(run, stage):
        print(f"[cli] {stage}: already done ({run.path}/.done_{stage}) — skipping. --force to redo.")
        return
    print(f"[cli] === stage: {stage} ===")
    DISPATCH[stage](cfg, run)
    write_stage_marker(run, stage)


def main(argv: list[str] | None = None) -> int:
    """Entrypoint for ``python -m phylo``."""
    parser = argparse.ArgumentParser(prog="phylo")
    parser.add_argument("stage", choices=[*ALL_STAGES, "all"])
    parser.add_argument("--config", required=True)
    parser.add_argument("--name", default=None, help="Run name (default: config file stem)")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    if args.name:
        cfg["name"] = args.name
    run = RunDir.create(cfg["output_root"], cfg["name"], cfg)

    stages = list(ALL_STAGES) if args.stage == "all" else [args.stage]
    for s in stages:
        _run_stage(s, cfg, run, force=args.force)
    return 0


if __name__ == "__main__":
    sys.exit(main())
