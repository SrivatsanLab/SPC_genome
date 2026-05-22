#!/usr/bin/env python
"""
Build a donor manifest mapping every benchmarking sample (HSC4 cells/bulk +
public_WGA SRR runs) to a donor_id and a sample_type ("cell" or "bulk").

The donor groups are:
  HSC4              - CapWGS HSCs (this study)                       cells + bulk
  BJ_fibroblast     - LIANTI paper, BJ foreskin fibroblast donor     cells + 2 bulks
  CD34_cord_blood   - PTA paper, CD34+ cord blood (matched to HSC4)  cells + 1 bulk
  B_lymphocyte      - PTA paper, B-lymphocyte donor                  cells + 1 bulk

Cells from public_WGA point at the 100M downsampled BAM
(data/benchmarking/downsampled/public_WGA/{SRR}_100M.bam) so that all
single-cell joint calling is at a matched depth.
"""

from pathlib import Path
import argparse
import csv
import sys

REPO = Path("/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/SPC_genome")
PTA_META = REPO / "results/benchmarking/PTA_meta.tsv"
LIANTI_META = REPO / "results/benchmarking/LIANTI_meta.tsv"
LIANTI_META_ADDL = REPO / "results/benchmarking/LIANTI_meta_additional.tsv"
HSC4_SC_DIR = REPO / "data/HSC4/sc_outputs"
HSC4_BULK_BAM = REPO / "data/HSC4_bulk/bams/bulkHSC4.bqsr.bam"
DOWNSAMPLED_DIR = REPO / "data/benchmarking/downsampled/public_WGA"

REF_GATK = "/shared/biodata/reference/GATK/hg38"


def read_tsv(path):
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def classify_pta_method(row):
    """Map PTA_meta amplification_method to canonical method string, or None to drop."""
    am = (row.get("amplification_method") or "").strip()
    if am == "PTA":
        return "PTA"
    if am == "SCMDA":
        return "MDA"
    if am == "":
        # Cord blood cells in PTA paper: blank == PTA
        return "PTA"
    if am == "None":
        # Bulk control samples
        return "bulk"
    return None


def lianti_method(row):
    # In LIANTI_meta.tsv the human-readable label is in `replicate`, not `isolate`.
    iso = row["replicate"]
    if "Bulk" in iso:
        return "bulk"
    if iso.startswith("LIANTI"):
        return "LIANTI"
    if "MDA (Qiagen)" in iso:
        return "MDA_Qiagen"
    if "MDA (GE)" in iso:
        return "MDA_GE"
    if "MALBAC (Yikon)" in iso:
        return "MALBAC_Yikon"
    if "MALBAC-like (Rubicon)" in iso:
        return "MALBAC_Rubicon"
    if "DOP-PCR" in iso:
        return "DOP_PCR"
    # UV / S-phase / UDG cells - keep but mark by isolate
    if "UV " in iso or "UDG" in iso or "S-Phase" in iso:
        return "LIANTI"  # All UV/S-phase/UDG used LIANTI per the paper
    return "unknown"


def main(out_path):
    rows = []

    # ---- HSC4 (CapWGS) ----
    for bam in sorted(HSC4_SC_DIR.glob("*.preprocessed.bam")):
        cell_id = bam.stem.replace(".preprocessed", "")
        gvcf = HSC4_SC_DIR / f"{cell_id}.g.vcf.gz"
        rows.append({
            "donor_id": "HSC4",
            "dataset": "HSC4",
            "sample_id": cell_id,
            "sample_type": "cell",
            "method": "CapWGS",
            "bam_path": str(bam),
            "gvcf_path": str(gvcf) if gvcf.exists() else "",
            "tissue": "HSC",
            "sex": "",
            "reference_dir": REF_GATK,
        })

    rows.append({
        "donor_id": "HSC4",
        "dataset": "HSC4",
        "sample_id": "bulkHSC4",
        "sample_type": "bulk",
        "method": "bulk",
        "bam_path": str(HSC4_BULK_BAM),
        "gvcf_path": "",  # exists as non-GVCF; will recall in GVCF mode
        "tissue": "HSC",
        "sex": "",
        "reference_dir": REF_GATK,
    })

    # ---- PTA paper ----
    for r in read_tsv(PTA_META):
        srr = r["SRA_Run"].strip()
        method = classify_pta_method(r)
        if method is None:
            continue
        tissue = r["tissue"].strip()
        # Determine donor by tissue/cell_type
        if "Cord Blood" in tissue or r.get("cell_subtype", "").strip() == "CD34+":
            donor = "CD34_cord_blood"
        elif "B-Lymphocyte" in (r.get("cell_type") or ""):
            donor = "B_lymphocyte"
        else:
            donor = "PTA_other"

        if method == "bulk":
            sample_type = "bulk"
            # Bulks for PTA paper donors: use the full (non-downsampled) BAM
            bam = REPO / f"data/public_WGA/{srr}.bam"
        else:
            sample_type = "cell"
            bam = DOWNSAMPLED_DIR / f"{srr}_100M.bam"

        rows.append({
            "donor_id": donor,
            "dataset": "public_WGA",
            "sample_id": srr,
            "sample_type": sample_type,
            "method": method,
            "bam_path": str(bam),
            "gvcf_path": "",
            "tissue": tissue,
            "sex": r["sex"].strip(),
            "reference_dir": REF_GATK,
        })

    # ---- LIANTI paper (BJ fibroblast donor) ----
    seen_runs = set(r["sample_id"] for r in rows)
    lianti_rows = read_tsv(LIANTI_META) + read_tsv(LIANTI_META_ADDL)
    for r in lianti_rows:
        srr = r["SRA_Run"].strip()
        if srr in seen_runs:
            continue
        seen_runs.add(srr)
        method = lianti_method(r)
        if method == "bulk":
            sample_type = "bulk"
            bam = REPO / f"data/public_WGA/{srr}.bam"
        else:
            sample_type = "cell"
            bam = DOWNSAMPLED_DIR / f"{srr}_100M.bam"

        rows.append({
            "donor_id": "BJ_fibroblast",
            "dataset": "public_WGA",
            "sample_id": srr,
            "sample_type": sample_type,
            "method": method,
            "bam_path": str(bam),
            "gvcf_path": "",
            "tissue": r["tissue"].strip(),
            "sex": r["sex"].strip(),
            "reference_dir": REF_GATK,
        })

    # ---- write ----
    cols = [
        "donor_id", "dataset", "sample_id", "sample_type", "method",
        "bam_path", "gvcf_path", "tissue", "sex", "reference_dir",
    ]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # ---- summary ----
    print(f"Wrote {len(rows)} samples to {out_path}", file=sys.stderr)
    from collections import Counter
    by_donor = Counter((r["donor_id"], r["sample_type"], r["method"]) for r in rows)
    print("\nDonor breakdown (donor / sample_type / method : count):", file=sys.stderr)
    for k, v in sorted(by_donor.items()):
        print(f"  {k[0]:<18} {k[1]:<5} {k[2]:<15} {v}", file=sys.stderr)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out", default=str(REPO / "data/benchmarking/manifests/donor_manifest.tsv"))
    args = p.parse_args()
    main(args.out)
