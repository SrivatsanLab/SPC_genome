"""Build per-worm cell/BAM lists and summary from the assignment h5ad.

Reads adata.obs["donor"]/["bam_path"]/["alpha_hap"]/["purity"]/["n_core_sites"];
writes one <worm>.cells.txt and one <worm>.bampaths per worm, plus
worm_summary.csv with n_cells / med_alpha / med_purity / est_depth. Excludes
`background` capsules and cells with alpha_hap >= 0.9.

Usage:
    python build_worm_lists.py <h5ad> <output_dir> [--min-cells 50]

Writes:
    <output_dir>/worm_summary.csv
    <output_dir>/worm_lists/<worm>.cells.txt
    <output_dir>/worm_lists/<worm>.bampaths
    <output_dir>/worm_lists/passing_worms.txt   # worms with n_cells >= min_cells
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("h5ad")
    p.add_argument("output_dir")
    p.add_argument("--min-cells", type=int, default=50)
    p.add_argument("--alpha-cutoff", type=float, default=0.9)
    args = p.parse_args()

    out = Path(args.output_dir)
    lists_dir = out / "worm_lists"
    lists_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.h5ad} (backed='r')…")
    adata = ad.read_h5ad(args.h5ad, backed="r")
    obs = adata.obs

    for col in ("donor", "bam_path", "alpha_hap", "purity", "n_core_sites"):
        if col not in obs.columns:
            print(f"Error: adata.obs missing column '{col}'", file=sys.stderr)
            return 1

    keep = obs["donor"].astype(str).str.startswith("worm")
    keep &= obs["alpha_hap"].astype(float) < args.alpha_cutoff
    kept = obs[keep].copy()
    print(f"  {int(keep.sum())} cells pass (background excluded; alpha_hap < {args.alpha_cutoff})")

    # Spot-check that bam_path entries exist on disk.
    n_spot = min(10, len(kept))
    spot = kept["bam_path"].sample(n_spot, random_state=0)
    missing = [p for p in spot if not Path(str(p)).is_file()]
    if missing:
        print(f"Error: {len(missing)}/{n_spot} spot-checked BAM paths missing, e.g. {missing[:3]}",
              file=sys.stderr)
        return 1
    print(f"  spot-checked {n_spot} BAM paths on disk — all present")

    summary = (
        kept.groupby("donor", observed=True)
        .agg(
            n_cells=("purity", "size"),
            med_alpha=("alpha_hap", "median"),
            med_purity=("purity", "median"),
            med_sites=("n_core_sites", "median"),
        )
        .sort_values("n_cells", ascending=False)
    )
    # Rough per-worm depth estimate: (mean reads per cell) * (mean pct aligned) * (read_length / genome_size).
    # Real depth is measured post-merge (see merge_pseudobulk_array.sh), so this is order-of-magnitude only.
    summary["est_depth"] = summary["n_cells"] * 1.15 * 0.25
    summary["passes"] = summary["n_cells"] >= args.min_cells
    summary.to_csv(out / "worm_summary.csv")
    print(f"\n{summary.to_string()}\n")
    print(f"  wrote {out / 'worm_summary.csv'}")

    passing = summary.index[summary["passes"]].tolist()
    dropped = summary.index[~summary["passes"]].tolist()

    for worm, grp in kept.groupby("donor", observed=True):
        (lists_dir / f"{worm}.cells.txt").write_text("\n".join(grp.index.astype(str)) + "\n")
        (lists_dir / f"{worm}.bampaths").write_text("\n".join(grp["bam_path"].astype(str)) + "\n")

    (lists_dir / "passing_worms.txt").write_text("\n".join(passing) + "\n")

    print(f"  wrote per-worm cell/bampath files to {lists_dir}")
    print(f"  passing worms (n_cells >= {args.min_cells}): {len(passing)} — {passing}")
    if dropped:
        print(f"  dropped worms (< {args.min_cells} cells): {dropped}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
