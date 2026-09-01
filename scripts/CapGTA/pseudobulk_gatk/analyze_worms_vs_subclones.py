"""Downstream analyses on the regenotyped pseudobulk callset.

Runs plan §6:
    6.1 per-worm hom vs het counts (already in worm_genotype_counts.csv from regenotype step)
    6.2 private-variant zygosity: fraction of each worm's private variants that are hom
        (worm vs somatic-subclone discriminator)
    6.3 new sites vs existing bcftools callset, by TYPE
    6.4 germline blacklist BED (union of non-reference positions across worms)
    +   pairwise homozygous-site sharing matrix (fused worms detection)
    +   REPORT.md summarizing the findings

Usage:
    python analyze_worms_vs_subclones.py <output_dir> <existing_vcf>
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


def load_regenotyped(path: Path) -> tuple[pd.DataFrame, list[str]]:
    tab = pd.read_csv(path, sep="\t")
    worms = sorted({c[:-3] for c in tab.columns if c.endswith("_gt")})
    return tab, worms


def gt_matrix(tab: pd.DataFrame, worms: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Return (is_hom, is_alt) boolean matrices, shape (n_sites, n_worms)."""
    is_hom = np.column_stack([tab[f"{w}_gt"].values == "1/1" for w in worms])
    is_alt = np.column_stack([tab[f"{w}_gt"].isin(["0/1", "1/1"]).values for w in worms])
    return is_hom, is_alt


def variant_ids(tab: pd.DataFrame) -> pd.Series:
    """Match cellspec's chr-pos-ref>alt variant_id convention (see load_vcf.py)."""
    return (tab["CHROM"].astype(str) + "-"
            + tab["POS"].astype(str) + "-"
            + tab["REF"].astype(str) + ">"
            + tab["ALT"].astype(str))


def existing_variant_ids(vcf: Path) -> set[str]:
    """Extract chrom-pos-ref>alt IDs from an existing bcftools VCF (worm6 joint call)."""
    print(f"Extracting variant IDs from {vcf}…")
    proc = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(vcf)],
        check=True, capture_output=True, text=True,
    )
    ids: set[str] = set()
    for line in proc.stdout.splitlines():
        chrom, pos, ref, alt = line.rstrip("\n").split("\t")
        # bcftools norm -m -both was already applied to worm6 joint VCF → ALT is single-allele.
        ids.add(f"{chrom}-{pos}-{ref}>{alt}")
    return ids


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("output_dir")
    p.add_argument("existing_vcf")
    args = p.parse_args()

    out = Path(args.output_dir)
    existing_vcf = Path(args.existing_vcf)
    for f in (out / "regenotyped.tsv", existing_vcf):
        if not f.exists():
            print(f"Error: required input not found: {f}", file=sys.stderr)
            return 1

    analyses = out / "analyses"
    analyses.mkdir(parents=True, exist_ok=True)

    tab, worms = load_regenotyped(out / "regenotyped.tsv")
    print(f"Loaded {len(tab)} sites × {len(worms)} worms")

    is_hom, is_alt = gt_matrix(tab, worms)
    vids = variant_ids(tab)

    # --- 6.2 Private-variant zygosity ---------------------------------------
    print("\n[6.2] Private-variant zygosity per worm…")
    n_alt_across = is_alt.sum(axis=1)
    priv = n_alt_across == 1                      # exactly one worm carries it
    rows = []
    for j, w in enumerate(worms):
        priv_j = priv & is_alt[:, j]
        n_priv = int(priv_j.sum())
        n_priv_hom = int((priv_j & is_hom[:, j]).sum())
        frac_hom = (n_priv_hom / n_priv) if n_priv else np.nan
        rows.append((w, n_priv, n_priv_hom, frac_hom))
    priv_df = pd.DataFrame(rows,
                           columns=["worm", "n_private", "n_private_hom", "frac_hom"])
    priv_df.to_csv(analyses / "private_variant_zygosity.csv", index=False)
    print(priv_df.to_string(index=False))
    #  Interpretation: high frac_hom → distinct germline lineage (genuine worm).
    #                  low  frac_hom → somatic subclone wrongly split from a parent worm.

    # --- 6.4 + fused-worm detection: pairwise hom sharing --------------------
    print("\n[6.4+] Pairwise hom-site sharing matrix…")
    #  For each worm, size = n hom sites. Sharing[i,j] = |hom_i ∩ hom_j| / min(|hom_i|,|hom_j|).
    hom_counts = is_hom.sum(axis=0)
    n = len(worms)
    share = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            inter = int((is_hom[:, i] & is_hom[:, j]).sum())
            denom = min(hom_counts[i], hom_counts[j])
            share[i, j] = (inter / denom) if denom else 0.0
    share_df = pd.DataFrame(share, index=worms, columns=worms)
    share_df.to_csv(analyses / "hom_sharing_matrix.csv")

    fig, ax = plt.subplots(figsize=(0.5 * n + 2, 0.5 * n + 2))
    im = ax.imshow(share, cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(n)); ax.set_xticklabels(worms, rotation=90, fontsize=8)
    ax.set_yticks(range(n)); ax.set_yticklabels(worms, fontsize=8)
    ax.set_title("Pairwise hom-site sharing (|A∩B|/min(|A|,|B|))")
    fig.colorbar(im, ax=ax, shrink=0.7)
    fig.tight_layout()
    fig.savefig(analyses / "hom_sharing_matrix.png", dpi=140)
    plt.close(fig)

    # Flag fused pairs (>90% overlap) — these are one worm split into two, or
    # unresolvably related siblings.
    fused = []
    for i in range(n):
        for j in range(i + 1, n):
            if share[i, j] > 0.9:
                fused.append((worms[i], worms[j], float(share[i, j])))
    if fused:
        print("  *** fused worm pairs (>90% hom sharing) ***")
        for a, b, s in fused:
            print(f"    {a} <-> {b}   share={s:.3f}")
    pd.DataFrame(fused, columns=["worm_a", "worm_b", "share"]).to_csv(
        analyses / "fused_worm_pairs.csv", index=False)

    # --- 6.3 New sites vs existing bcftools callset --------------------------
    print("\n[6.3] Comparing to existing bcftools callset…")
    existing = existing_variant_ids(existing_vcf)
    print(f"  existing callset: {len(existing)} variant IDs")

    new_mask = ~vids.isin(existing)
    #  "new" restricted to sites where at least one worm has an alt call — GATK
    #  emits reference-confidence sites during genotyping, but those are not
    #  discoveries.
    interesting = new_mask & (n_alt_across > 0)
    summary_by_type = (
        tab.loc[interesting, "TYPE"].value_counts().rename("new_sites")
    )
    existing_by_type = (
        tab.loc[~new_mask & (n_alt_across > 0), "TYPE"].value_counts().rename("overlap_sites")
    )
    type_summary = pd.concat([summary_by_type, existing_by_type], axis=1).fillna(0).astype(int)
    type_summary.to_csv(analyses / "new_sites_summary.csv")
    (analyses / "new_sites_summary.txt").write_text(
        f"Existing callset: {len(existing)} variant IDs\n"
        f"HC callset:       {int((n_alt_across > 0).sum())} variant sites (>=1 alt worm)\n"
        f"New in HC:        {int(interesting.sum())}\n"
        f"Overlap:          {int((~new_mask & (n_alt_across > 0)).sum())}\n"
        f"\nBy TYPE:\n{type_summary.to_string()}\n"
    )
    print(type_summary.to_string())

    # --- 6.4 germline blacklist BED (union of non-reference positions) ------
    print("\n[6.4] Germline blacklist BED…")
    #  A site is on the blacklist if any worm was called non-reference at it.
    #  BED is 0-based half-open — VCF POS is 1-based, so use (POS-1, POS).
    bl = tab.loc[n_alt_across > 0, ["CHROM", "POS"]].copy()
    bl["start"] = bl["POS"].astype(int) - 1
    bl["end"] = bl["POS"].astype(int)
    bl[["CHROM", "start", "end"]].to_csv(
        analyses / "germline_blacklist.bed", sep="\t", header=False, index=False)
    print(f"  wrote {analyses / 'germline_blacklist.bed'}   ({len(bl)} positions)")

    # --- REPORT.md ----------------------------------------------------------
    report = analyses / "REPORT.md"
    report.write_text(
        f"# Worm6 pseudobulk GATK callset — analysis report\n\n"
        f"Generated by `analyze_worms_vs_subclones.py`.\n\n"
        f"## Inputs\n"
        f"- Regenotyped table: `regenotyped.tsv` ({len(tab)} sites × {len(worms)} worms)\n"
        f"- Existing callset:  `{existing_vcf}` ({len(existing)} variant IDs)\n\n"
        f"## 6.1 Per-worm hom vs het counts\n"
        f"See `worm_genotype_counts.csv` (written by `regenotype_from_ad.py`).\n\n"
        f"## 6.2 Private-variant zygosity — worm vs subclone discriminator\n"
        f"See `private_variant_zygosity.csv`.\n\n"
        f"Interpretation: high `frac_hom` → distinct germline lineage (a real worm). "
        f"Low `frac_hom` → somatic subclone that was wrongly split from a parent worm; "
        f"its private variants are somatic (het), and it shares most of its hom sites "
        f"with the parent worm.\n\n"
        f"## 6.3 New sites vs existing bcftools callset (by TYPE)\n\n"
        f"```\n{type_summary.to_string()}\n```\n\n"
        f"See `new_sites_summary.txt` for the counts and `new_sites_summary.csv` for the table.\n\n"
        f"## 6.4 Germline blacklist\n"
        f"`germline_blacklist.bed` — union of {len(bl)} non-reference positions across "
        f"all worms. This is **necessary but not sufficient** for downstream somatic "
        f"calling — the ambient pool also carries other cells' somatic mutations, so "
        f"read-count thresholds are still required at the somatic-calling stage.\n\n"
        f"## Fused-worm detection\n"
        f"Pairs with hom-site overlap > 0.9 in `fused_worm_pairs.csv`:\n\n"
        f"```\n{pd.DataFrame(fused, columns=['worm_a', 'worm_b', 'share']).to_string(index=False) if fused else 'none'}\n```\n\n"
        f"See heatmap: `hom_sharing_matrix.png`.\n\n"
        f"## Caveats (from GATK_PSEUDOBULK_PLAN.md §10)\n"
        f"- Cell assignments are inferred; every pseudobulk inherits assignment error.\n"
        f"- This analysis cannot validate the demultiplexing that produced its input.\n"
        f"- Contamination is graded, not discrete. Per-worm α is a summary of a distribution; "
        f"sites near the hom/het boundary are genuinely ambiguous.\n"
        f"- Worms differ in cell count by ~30×; per-worm depth and sensitivity differ likewise. "
        f"Never compare raw variant counts between worms without depth normalization.\n"
    )
    print(f"\nWrote {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
