#!/usr/bin/env python3
"""Build a per-worm panel set and write it to ``panels/<tag>/``.

Wraps ``phylo.pp.build_panels`` so its ``run.subdir('panels')`` writes straight
into the tag directory instead of the ``runs/<name>/panels/`` layout used by the
main CLI. Cheaper: no stray run directory carrying every other stage's config.

Usage
-----
  build_panel_at_tag.py \
      --config  scripts/phylogeny/configs/tct_kept.yaml \
      --panel-tag tct_kept
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import yaml

REPO = Path("/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome")
sys.path.insert(0, str(REPO / "scripts" / "phylogeny"))

from phylo.pp.panels import build_panels  # noqa: E402
from phylo.pp.qc import apply_qc  # noqa: E402
from phylo.utils.config import load_config  # noqa: E402


class RunDirTagged:
    """RunDir stand-in that routes ``subdir('panels')`` to a specific tag dir."""

    def __init__(self, root: Path, panels_target: Path):
        self.path = root
        self._panels_target = panels_target

    def subdir(self, stage: str) -> Path:
        if stage == "panels":
            self._panels_target.mkdir(parents=True, exist_ok=True)
            return self._panels_target
        d = self.path / stage
        d.mkdir(parents=True, exist_ok=True)
        return d


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--panel-tag", required=True)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    adata_path = Path(cfg["inputs"]["adata"])
    if not adata_path.is_absolute():
        adata_path = (args.config.parent / adata_path).resolve()
    print(f"[build_panel] loading adata {adata_path}")
    adata = ad.read_h5ad(adata_path)

    print(f"[build_panel] applying qc predicates")
    apply_qc(adata, cfg["qc"]["predicates"], inplace=True)

    phylo_root = Path(cfg["output_root"])
    if not phylo_root.is_absolute():
        phylo_root = (args.config.parent / phylo_root).resolve()
    target = phylo_root / "panels" / args.panel_tag
    run = RunDirTagged(phylo_root, target)

    print(f"[build_panel] writing panels to {target}")
    build_panels(adata, cfg, run)

    with open(target / "panel_config.resolved.yaml", "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    print(f"[build_panel] wrote panel_config.resolved.yaml")


if __name__ == "__main__":
    main()
