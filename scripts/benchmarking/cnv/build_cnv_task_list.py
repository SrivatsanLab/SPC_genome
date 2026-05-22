#!/usr/bin/env python
"""Build the AneuFinder array task list from the donor manifest.

Maps each manifest row to the BAM that AneuFinder should use:

  * public cells -> data/benchmarking/downsampled/public_WGA/<sample>_25M.bam
  * HSC4 cells   -> manifest bam_path (native depth; 25M variant TBD)
  * public bulks -> data/benchmarking/bulks/<donor>/<sample>.bqsr.bam
  * HSC4 bulk    -> manifest bam_path (already MarkDup'd + BQSR'd)

Rows whose target BAM doesn't exist on disk are dropped (with a warning).

Output:
    <out_tsv>: donor_id  sample_id  bam_path  is_bulk

Run:
    python scripts/benchmarking/cnv/build_cnv_task_list.py \
        --manifest data/benchmarking/manifests/donor_manifest.tsv \
        --out data/benchmarking/manifests/cnv_tasks_25M.tsv
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[3]


def target_bam(row) -> Path:
    """Map a manifest row to the BAM we want AneuFinder to read."""
    donor = row["donor_id"]
    sample = row["sample_id"]
    is_bulk = row["sample_type"] == "bulk"

    if is_bulk:
        if donor == "HSC4":
            return REPO / "data/HSC4_bulk/bams/bulkHSC4.bqsr.bam"
        return REPO / f"data/benchmarking/bulks/{donor}/{sample}.bqsr.bam"

    # Cells
    if donor == "HSC4":
        # Native depth (25M HSC4 BAMs don't exist yet — see plan prereq).
        return Path(row["bam_path"])
    return REPO / f"data/benchmarking/downsampled/public_WGA/{sample}_25M.bam"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True,
                    help="data/benchmarking/manifests/donor_manifest.tsv")
    ap.add_argument("--out", required=True, help="Output task TSV")
    ap.add_argument("--include-hsc4-native", action="store_true",
                    help="Include HSC4 cells at native depth (default: skip "
                         "until 25M HSC4 BAMs are built — keeps comparison "
                         "depth-matched).")
    args = ap.parse_args()

    df = pd.read_csv(args.manifest, sep="\t")
    df["bam_target"] = df.apply(target_bam, axis=1).astype(str)
    df["is_bulk"] = (df["sample_type"] == "bulk").astype(int)

    if not args.include_hsc4_native:
        before = len(df)
        df = df[~((df["donor_id"] == "HSC4") & (df["sample_type"] == "cell"))]
        dropped = before - len(df)
        if dropped:
            print(f"Skipped {dropped} HSC4 cells at native depth (use "
                  f"--include-hsc4-native to keep them).", file=sys.stderr)

    missing = ~df["bam_target"].map(lambda p: Path(p).is_file())
    if missing.any():
        for _, r in df[missing].iterrows():
            print(f"WARN missing BAM: {r['donor_id']}/{r['sample_id']} "
                  f"-> {r['bam_target']}", file=sys.stderr)
        df = df[~missing]

    out_df = df[["donor_id", "sample_id", "bam_target", "is_bulk"]].rename(
        columns={"bam_target": "bam_path"})

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)

    print(f"Wrote {len(out_df)} tasks to {args.out}")
    print(out_df.groupby(["donor_id", "is_bulk"]).size()
          .reset_index(name="n").to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
