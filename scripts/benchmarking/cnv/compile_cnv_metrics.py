#!/usr/bin/env python
"""Compile per-sample AneuFinder outputs into long-form benchmark tables.

Reads:
    data/benchmarking/cnv/<donor>/<sample>.qc.tsv
    data/benchmarking/cnv/<donor>/<sample>.segments.tsv

Writes:
    <out>/cnv_qc.tsv          one row per sample, joined to donor manifest +
                              is_bulk flag from the task list.
    <out>/cnv_segments.tsv    long-form segment table, all samples concatenated.
    <out>/cnv_fp_segments.tsv per-cell vs paired-bulk FP count (cells only).

Run:
    python scripts/benchmarking/cnv/compile_cnv_metrics.py \
        --cnv-root data/benchmarking/cnv \
        --task-list data/benchmarking/manifests/cnv_tasks_25M.tsv \
        --manifest data/benchmarking/manifests/donor_manifest.tsv \
        --out data/benchmarking/cnv/_compiled
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path

import pandas as pd


def read_qcs(cnv_root: Path) -> pd.DataFrame:
    rows = []
    for f in cnv_root.glob("*/*.qc.tsv"):
        df = pd.read_csv(f, sep="\t")
        df["_donor_dir"] = f.parent.name
        rows.append(df)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def read_segs(cnv_root: Path) -> pd.DataFrame:
    rows = []
    for f in cnv_root.glob("*/*.segments.tsv"):
        df = pd.read_csv(f, sep="\t")
        if df.empty:
            continue
        df["_donor_dir"] = f.parent.name
        rows.append(df)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def overlap_bp(a: pd.DataFrame, b: pd.DataFrame) -> int:
    """Total base-pairs of overlap between two same-chromosome segment lists."""
    if a.empty or b.empty:
        return 0
    total = 0
    for chrom in set(a["seqnames"]).intersection(b["seqnames"]):
        aa = a[a["seqnames"] == chrom][["start", "end"]].values
        bb = b[b["seqnames"] == chrom][["start", "end"]].values
        for sa, ea in aa:
            for sb, eb in bb:
                lo = max(sa, sb)
                hi = min(ea, eb)
                if hi > lo:
                    total += hi - lo
    return total


def compute_fp_segments(segs: pd.DataFrame, task: pd.DataFrame,
                        method: str = "HMM") -> pd.DataFrame:
    """Cell-vs-bulk FP non-diploid segment count, per cell, per donor.

    FP segment = a non-diploid cell segment that does NOT overlap any
    non-diploid bulk segment (in the same donor). We measure both the count
    and the bp coverage, since a single long FP is qualitatively different
    from many short ones.
    """
    if segs.empty or task.empty:
        return pd.DataFrame()
    segs = segs[segs["method"] == method].copy()
    segs = segs[segs["copy.number"] != 2]
    if segs.empty:
        return pd.DataFrame(columns=["donor_id", "sample", "n_fp_segs",
                                     "fp_bp", "n_bulk_nondip_segs"])

    bulks_by_donor = {}
    for donor in task["donor_id"].unique():
        bulk_samples = task.loc[(task["donor_id"] == donor) &
                                (task["is_bulk"] == 1), "sample_id"].tolist()
        if not bulk_samples:
            continue
        bulk_segs = segs[segs["sample"].isin(bulk_samples)]
        bulks_by_donor[donor] = bulk_segs

    out_rows = []
    cell_task = task[task["is_bulk"] == 0]
    for _, t in cell_task.iterrows():
        donor, sample = t["donor_id"], t["sample_id"]
        cell_segs = segs[segs["sample"] == sample]
        if cell_segs.empty:
            continue
        bulk_segs = bulks_by_donor.get(donor)
        if bulk_segs is None or bulk_segs.empty:
            # No bulk truth for this donor: report cell counts, NA for FP
            out_rows.append({
                "donor_id": donor, "sample": sample,
                "n_cell_nondip_segs": len(cell_segs),
                "cell_nondip_bp": int(cell_segs["width"].sum()),
                "n_bulk_nondip_segs": pd.NA,
                "n_fp_segs": pd.NA,
                "fp_bp": pd.NA,
            })
            continue

        n_fp = 0
        fp_bp = 0
        for _, c in cell_segs.iterrows():
            same_chrom_bulk = bulk_segs[bulk_segs["seqnames"] == c["seqnames"]]
            overlap = any(
                (c["start"] <= eb and c["end"] >= sb)
                for sb, eb in zip(same_chrom_bulk["start"], same_chrom_bulk["end"])
            )
            if not overlap:
                n_fp += 1
                fp_bp += int(c["width"])
        out_rows.append({
            "donor_id": donor, "sample": sample,
            "n_cell_nondip_segs": len(cell_segs),
            "cell_nondip_bp": int(cell_segs["width"].sum()),
            "n_bulk_nondip_segs": len(bulk_segs),
            "n_fp_segs": n_fp,
            "fp_bp": fp_bp,
        })
    return pd.DataFrame(out_rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cnv-root", required=True,
                    help="Root dir containing per-donor AneuFinder outputs")
    ap.add_argument("--task-list", required=True,
                    help="data/benchmarking/manifests/cnv_tasks_25M.tsv")
    ap.add_argument("--manifest", required=True,
                    help="data/benchmarking/manifests/donor_manifest.tsv")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--method", default="HMM", choices=["HMM", "edivisive"],
                    help="Which segmentation method to use for FP analysis")
    args = ap.parse_args()

    cnv_root = Path(args.cnv_root)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    qc = read_qcs(cnv_root)
    segs = read_segs(cnv_root)
    task = pd.read_csv(args.task_list, sep="\t")
    manifest = pd.read_csv(args.manifest, sep="\t")[
        ["donor_id", "sample_id", "sample_type", "method", "tissue", "sex"]
    ]

    if qc.empty:
        print(f"ERROR: no QC files under {cnv_root}", file=sys.stderr)
        return 1

    qc_full = (qc.rename(columns={"sample": "sample_id"})
                 .merge(task[["donor_id", "sample_id", "is_bulk"]],
                        on="sample_id", how="left")
                 .merge(manifest, on=["donor_id", "sample_id"], how="left"))
    qc_full.to_csv(out_dir / "cnv_qc.tsv", sep="\t", index=False)
    print(f"Wrote {len(qc_full)} QC rows -> {out_dir/'cnv_qc.tsv'}")

    if not segs.empty:
        seg_full = (segs.rename(columns={"sample": "sample_id"})
                       .merge(task[["donor_id", "sample_id", "is_bulk"]],
                              on="sample_id", how="left")
                       .merge(manifest, on=["donor_id", "sample_id"], how="left"))
        seg_full.to_csv(out_dir / "cnv_segments.tsv", sep="\t", index=False)
        print(f"Wrote {len(seg_full)} segment rows -> {out_dir/'cnv_segments.tsv'}")

        fp = compute_fp_segments(
            segs.rename(columns={"sample": "sample"}), task, method=args.method)
        if not fp.empty:
            fp = fp.merge(manifest, left_on=["donor_id", "sample"],
                          right_on=["donor_id", "sample_id"], how="left")
            fp.to_csv(out_dir / "cnv_fp_segments.tsv", sep="\t", index=False)
            print(f"Wrote {len(fp)} FP-segment rows -> {out_dir/'cnv_fp_segments.tsv'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
