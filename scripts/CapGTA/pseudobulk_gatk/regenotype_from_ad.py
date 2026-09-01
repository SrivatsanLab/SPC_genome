"""Re-genotype the joint VCF from AD/DP using per-worm contamination.

GATK's diploid GT is contamination-blind: a homozygous germline site at pooled VAF
(1-alpha)*1 + alpha*(1/K) reads at 0.62 when alpha=0.4 and would be miscalled 0/1.
This script does the plan §5 correction:

    hom_lo = 0.75 * (1 - alpha_worm)   # >= hom_lo    -> 1/1
    het_lo = 0.25 * (1 - alpha_worm)   # >= het_lo    -> 0/1
    else                                -> 0/0
    dp < MIN_DP                         -> ./.

Also writes a per-worm VAF histogram — the sanity check that must pass before any
downstream interpretation (bimodal at ~(1-alpha) and ~(1-alpha)/2).

Requires GATK on PATH (calls `gatk VariantsToTable`) and a Python env with pandas
+ matplotlib.

Usage:
    python regenotype_from_ad.py <output_dir>
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MIN_DP = 10


def _ensure_tbi(vcf: Path) -> None:
    """GATK requires .tbi; bcftools defaults to .csi. Build .tbi if missing."""
    tbi = vcf.with_suffix(vcf.suffix + ".tbi")
    if tbi.exists():
        return
    cmd = ["gatk", "IndexFeatureFile", "-I", str(vcf)]
    print("Building .tbi index:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_variants_to_table(vcf: Path, out_tsv: Path) -> None:
    """Extract CHROM/POS/REF/ALT/TYPE/QUAL/FILTER + per-worm AD, DP, GQ."""
    _ensure_tbi(vcf)
    cmd = [
        "gatk",
        "VariantsToTable",
        "-V", str(vcf),
        "-O", str(out_tsv),
        "-F", "CHROM", "-F", "POS", "-F", "REF", "-F", "ALT",
        "-F", "TYPE", "-F", "QUAL", "-F", "FILTER",
        "-GF", "AD", "-GF", "DP", "-GF", "GQ",
        "--show-filtered",
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_ad(col: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Parse a comma-joined AD column into (ref_dp, alt_dp) numeric Series."""
    parts = col.astype(str).str.split(",", n=1, expand=True)
    ref = pd.to_numeric(parts[0], errors="coerce").fillna(0).astype(int)
    alt = pd.to_numeric(parts[1], errors="coerce").fillna(0).astype(int)
    return ref, alt


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("output_dir")
    p.add_argument("--min-dp", type=int, default=MIN_DP,
                   help=f"Min per-worm depth to genotype at all (default {MIN_DP})")
    args = p.parse_args()

    out = Path(args.output_dir)
    vcf = out / "joint.vcf.gz"
    summary_csv = out / "worm_summary.csv"
    for f in (vcf, summary_csv):
        if not f.exists():
            print(f"Error: required input not found: {f}", file=sys.stderr)
            return 1

    analyses = out / "analyses"
    analyses.mkdir(parents=True, exist_ok=True)

    tsv = out / "joint_table.tsv"
    if tsv.exists() and tsv.stat().st_mtime >= vcf.stat().st_mtime:
        print(f"Reusing existing {tsv}")
    else:
        if tsv.exists():
            print(f"{tsv} is older than {vcf} — regenerating")
        run_variants_to_table(vcf, tsv)

    print(f"Loading {tsv}…")
    tab = pd.read_csv(tsv, sep="\t")

    alpha = pd.read_csv(summary_csv, index_col=0)["med_alpha"]
    worms = sorted({c[:-3] for c in tab.columns if c.endswith(".AD")})
    print(f"  {len(tab)} sites, {len(worms)} worms")

    keep_cols = ["CHROM", "POS", "REF", "ALT", "TYPE", "QUAL", "FILTER"]
    out_frame = tab[keep_cols].copy()

    fig, axes = plt.subplots(
        (len(worms) + 3) // 4, 4,
        figsize=(16, 3 * ((len(worms) + 3) // 4)),
        squeeze=False,
    )

    for i, w in enumerate(worms):
        ref, alt = parse_ad(tab[f"{w}.AD"])
        dp = ref + alt
        vaf = alt / dp.replace(0, np.nan)
        a = float(alpha.get(w, 0.4))
        hom_lo = 0.75 * (1 - a)
        het_lo = 0.25 * (1 - a)

        gt = np.where(dp < args.min_dp, "./.",
             np.where(vaf >= hom_lo, "1/1",
             np.where(vaf >= het_lo, "0/1", "0/0")))

        out_frame[f"{w}_dp"] = dp
        out_frame[f"{w}_ad"] = alt
        out_frame[f"{w}_vaf"] = vaf
        out_frame[f"{w}_gt"] = gt

        # VAF histogram: non-reference calls only (0/1 or 1/1), depth-passing.
        mask = (gt != "0/0") & (gt != "./.")
        ax = axes[i // 4, i % 4]
        v = vaf[mask].dropna()
        ax.hist(v, bins=50, range=(0, 1), color="steelblue", edgecolor="none")
        ax.axvline(1 - a, color="crimson", ls="--", lw=1, label=f"1-α={1-a:.2f} (hom)")
        ax.axvline((1 - a) / 2, color="orange", ls="--", lw=1, label=f"(1-α)/2={((1-a)/2):.2f} (het)")
        ax.set_title(f"{w}  n={mask.sum()}  α={a:.2f}", fontsize=10)
        ax.set_xlabel("VAF")
        ax.set_ylabel("sites")
        ax.legend(fontsize=7, loc="upper center")

    # Blank any unused subplots.
    for j in range(len(worms), axes.size):
        axes.flat[j].axis("off")

    fig.tight_layout()
    fig_path = analyses / "vaf_histograms.png"
    fig.savefig(fig_path, dpi=140)
    plt.close(fig)
    print(f"  wrote {fig_path}")

    regeno = out / "regenotyped.tsv"
    out_frame.to_csv(regeno, sep="\t", index=False)
    print(f"  wrote {regeno}   ({len(out_frame)} rows, {out_frame.shape[1]} cols)")

    # Quick per-worm hom / het counts.
    print("\nPer-worm genotype counts (variant sites only):")
    counts_rows = []
    for w in worms:
        gt = out_frame[f"{w}_gt"]
        n_hom = int((gt == "1/1").sum())
        n_het = int((gt == "0/1").sum())
        n_ref = int((gt == "0/0").sum())
        n_miss = int((gt == "./.").sum())
        counts_rows.append((w, n_hom, n_het, n_ref, n_miss,
                            (n_hom / n_het) if n_het else np.nan))
    counts = pd.DataFrame(counts_rows,
                          columns=["worm", "n_hom", "n_het", "n_ref", "n_miss", "hom_over_het"])
    counts.to_csv(analyses / "worm_genotype_counts.csv", index=False)
    print(counts.to_string(index=False))
    print(f"  wrote {analyses / 'worm_genotype_counts.csv'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
