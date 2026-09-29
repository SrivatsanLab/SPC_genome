#!/usr/bin/env python3
"""MP encoding-stringency + parsimony-model sweep on one worm.

Sweeps configurations of (dp_threshold K, Camin-Sokal reversal cost k), running
phangorn MP + bootstrap for each and collecting stats side-by-side.

See docs/mp_encoding_and_model_sweep.md for design.

Usage
-----
  micromamba run -n cellphy python \
      scripts/phylogeny/inference/mp_encoding_model_sweep.py \
      --panel-tag tct_kept_relax03__svdImpSym_mquad10_cs3 \
      --worm worm10 \
      --out-root scratch/mp_sweep \
      [--bootstrap 100] [--threads 4]

Environment
-----------
Requires the `cellphy` micromamba env for the encoder, and the `fhR` module
for the R runner. This launcher runs those in-process (blocking) — suitable
for single-worm smell tests, not a scale-up sweep.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path("/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome")
PANELS = REPO / "results/worm6_final/DNA_analysis/phylogeny/panels"

# Config sweep. Cost pair is (gain_cost = 0->1 transition, cs_cost = 1->0 transition).
# gain_cost = 1, cs_cost = ∞  -> strict Camin-Sokal (irreversible gains)
# gain_cost = 1, cs_cost = 1  -> Fitch (symmetric)
# gain_cost = 100, cs_cost = 1 -> Dollo approximation (gains heavily penalized,
#                                 losses cheap; the optimizer strongly prefers
#                                 "one gain, many losses" per character).
CONFIGS = [
    dict(label="baseline_K1_kInf",   dp_threshold=1, gain_cost=1, cs_cost=1e9),
    dict(label="relax_K1_k20",       dp_threshold=1, gain_cost=1, cs_cost=20),
    dict(label="relax_K1_k10",       dp_threshold=1, gain_cost=1, cs_cost=10),
    dict(label="relax_K1_k5",        dp_threshold=1, gain_cost=1, cs_cost=5),
    dict(label="relax_K1_k3",        dp_threshold=1, gain_cost=1, cs_cost=3),
    dict(label="fitch_K1_k1",        dp_threshold=1, gain_cost=1, cs_cost=1),
    dict(label="stringent_K3_kInf", dp_threshold=3, gain_cost=1, cs_cost=1e9),
    dict(label="combined_K3_k5",    dp_threshold=3, gain_cost=1, cs_cost=5),
    # Dollo-approximation (asymmetric with high gain cost + cheap loss)
    dict(label="dollo_g10_l1",       dp_threshold=1, gain_cost=10,   cs_cost=1),
    dict(label="dollo_g100_l1",      dp_threshold=1, gain_cost=100,  cs_cost=1),
    dict(label="dollo_g1000_l1",     dp_threshold=1, gain_cost=1000, cs_cost=1),
]


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-tag", required=True)
    ap.add_argument("--worm", required=True, help="e.g. worm10")
    ap.add_argument("--out-root", type=Path, required=True,
                    help="Scratch directory root; per-config subdirs will be created.")
    ap.add_argument("--bootstrap", type=int, default=100)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--configs", nargs="+", default=None,
                    help="Optional subset of config labels to run. Default runs all.")
    return ap.parse_args()


def run(cmd, log_prefix=""):
    print(f"{log_prefix}$ {' '.join(str(x) for x in cmd)}", flush=True)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise SystemExit(f"command failed with code {proc.returncode}: {' '.join(str(x) for x in cmd)}")
    return proc.stdout, proc.stderr


def main():
    args = parse_args()
    panel_h5 = PANELS / args.panel_tag / f"worm_{args.worm}.h5ad"
    if not panel_h5.exists():
        raise SystemExit(f"panel h5ad not found: {panel_h5}")

    configs = CONFIGS
    if args.configs is not None:
        wanted = set(args.configs)
        configs = [c for c in CONFIGS if c["label"] in wanted]
        if not configs:
            raise SystemExit(f"no configs matched {args.configs}")

    print(f"[sweep] panel={args.panel_tag}  worm={args.worm}  n_configs={len(configs)}")
    print(f"[sweep] panel_h5={panel_h5}")

    args.out_root.mkdir(parents=True, exist_ok=True)
    summary_rows = []

    for cfg in configs:
        label = cfg["label"]
        K = cfg["dp_threshold"]
        gain = cfg.get("gain_cost", 1)
        loss = cfg["cs_cost"]
        cfg_dir = args.out_root / f"{args.panel_tag}__{args.worm}__{label}"
        cfg_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n=== {label} (K={K}, gain={gain}, loss={loss}) ===")

        # 1. Encode with this K threshold.
        run([
            "python", str(REPO / "scripts/phylogeny/inference/encode_mp_matrix.py"),
            "--panel-h5ad", str(panel_h5),
            "--out-dir",    str(cfg_dir),
            "--dp-threshold", str(K),
            "--tag", "real",
        ], log_prefix=f"[{label}] encode ")

        # 2. Run phangorn MP + bootstrap with the (gain, loss) cost pair.
        rcmd = (
            f"module load fhR && "
            f"Rscript {REPO}/scripts/phylogeny/inference/mp_run.R "
            f"--tsv {cfg_dir}/real.tsv "
            f"--out-dir {cfg_dir} "
            f"--tag real "
            f"--bootstrap {args.bootstrap} "
            f"--k 10 "
            f"--seed 0 "
            f"--threads {args.threads} "
            f"--gain-cost {gain} "
            f"--cs-cost {loss}"
        )
        run(["bash", "-lc", rcmd], log_prefix=f"[{label}] mp ")

        # 3. Read stats
        stats_p = cfg_dir / "real.stats.json"
        if not stats_p.exists():
            print(f"[{label}] WARNING: no stats.json produced")
            continue
        stats = json.loads(stats_p.read_text())

        # 4. Bootstrap-support distribution from mp_support.newick
        sup_stats = _support_stats(cfg_dir / "real.mp_support.newick")

        row = dict(label=label, K=K, gain=gain, loss=loss,
                   **{f"stats_{k2}": v for k2, v in stats.items()}, **sup_stats)
        summary_rows.append(row)

    # Print comparison table
    print("\n" + "=" * 130)
    print(f"SWEEP SUMMARY  panel={args.panel_tag}  worm={args.worm}")
    print("=" * 130)
    header = ["label", "K", "gain", "loss", "mp_score", "n_mpts", "ci", "ri", "hi",
              "res_strict", "res_maj", "mean_sup", "median_sup", "pct_ge50", "pct_ge70", "pct_ge90",
              "elapsed_s"]
    print("\t".join(header))
    for r in summary_rows:
        print("\t".join(str(_fmt(r.get(k, r.get(f"stats_{k}", "")))) for k in [
            "label", "K", "gain", "loss",
            "stats_mp_score", "stats_n_mpts", "stats_ci", "stats_ri", "stats_hi",
            "stats_resolution_strict", "stats_resolution_majority",
            "mean_sup", "median_sup", "pct_ge50", "pct_ge70", "pct_ge90",
            "stats_elapsed_search_s",
        ]))

    # Also write a machine-readable summary.
    out_summary = args.out_root / f"summary__{args.panel_tag}__{args.worm}.json"
    out_summary.write_text(json.dumps(summary_rows, indent=2, default=str))
    print(f"\n[sweep] wrote {out_summary}")


def _fmt(v):
    if isinstance(v, float):
        return f"{v:.3f}"
    return str(v)


def _support_stats(newick_path: Path) -> dict:
    """Return mean/median/percentile bootstrap support from an mp_support.newick."""
    if not newick_path.exists():
        return dict(mean_sup=None, median_sup=None, pct_ge50=None, pct_ge70=None, pct_ge90=None)
    import ete3
    import numpy as np
    nw = newick_path.read_text()
    nw = re.sub(r"\)NA(?=[,\)])", ")0", nw)
    t = ete3.Tree(nw, format=0)
    for lf in t.get_leaves(): lf.name = lf.name.strip("'")
    vals = []
    for n in t.traverse():
        if n.is_leaf() or n.is_root(): continue
        s = getattr(n, "support", None)
        try:
            f = float(s)
            if 0 <= f <= 100 and not np.isnan(f): vals.append(f)
        except (TypeError, ValueError): pass
    if not vals:
        return dict(mean_sup=None, median_sup=None, pct_ge50=None, pct_ge70=None, pct_ge90=None)
    v = np.array(vals)
    return dict(
        mean_sup=float(v.mean()),
        median_sup=float(np.median(v)),
        pct_ge50=float((v >= 50).mean() * 100),
        pct_ge70=float((v >= 70).mean() * 100),
        pct_ge90=float((v >= 90).mean() * 100),
    )


if __name__ == "__main__":
    main()
