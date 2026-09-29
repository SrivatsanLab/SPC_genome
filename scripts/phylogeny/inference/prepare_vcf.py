#!/usr/bin/env python3
"""Build per-worm CellPhy input VCFs.

For one panel h5ad, emits three VCFs, all with PL kept intact:

  real.vcf.gz       joint VCF subset to (panel variants) x (worm cells)
  covmask.vcf.gz    same coverage pattern; every covered site is a soft 0
                    (GT=0/0, AD=(DP,0), PL=(0,3,20) at DP>0; PL=(0,0,0) at DP=0)
  shuffled.vcf.gz   per-variant permutation of (GT,AD,PL) among DP>0 cells

Usage
-----
  prepare_vcf.py \
      --panel-h5ad phylogeny/panels/<tag>/worm_<W>.h5ad \
      --joint-vcf  DNA_analysis/joint_variants/joint_variants.vcf.gz \
      --out-dir    phylogeny/inference/<tag>/vcf/worm_<W> \
      --seed 0
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pysam


# PL used for a soft-0 covered site — matches the DP=1 flat prior from GATK.
# Weight = -log10(P(alt|read))*10 ~ 3 (het) vs 20 (hom-alt). Ref preferred by
# ~0.3 bits, which is the informational content of a single covered-ref read.
COVMASK_PL_COV = (0, 3, 20)
COVMASK_PL_UNCOV = (0, 0, 0)


# --- diploid biallelic PL model (err = 1e-3) --------------------------------
# Standard model per bcftools: PL_g = -10 * log10(P(reads | g)).
#   P(reads|0/0) = (1-ε)^(DP-AD) * ε^AD
#   P(reads|0/1) = 0.5^DP
#   P(reads|1/1) = ε^(DP-AD) * (1-ε)^AD
_LOG10_1MEPS = -0.000434  # log10(1 - 1e-3)
_LOG10_EPS = -3.0         # log10(1e-3)
_LOG10_HALF = -0.30103    # log10(0.5)

def _pl_from_ad_dp(ad_alt: int, dp: int) -> tuple[int, int, int]:
    """Return (PL_00, PL_01, PL_11) normalized so min == 0. Integer-rounded."""
    if dp <= 0:
        return (0, 0, 0)  # uninformative
    ad_ref = max(0, dp - ad_alt)
    ad_alt = max(0, min(ad_alt, dp))
    l00 = ad_ref * _LOG10_1MEPS + ad_alt * _LOG10_EPS
    l01 = dp * _LOG10_HALF
    l11 = ad_ref * _LOG10_EPS + ad_alt * _LOG10_1MEPS
    pl00 = int(round(-10 * l00))
    pl01 = int(round(-10 * l01))
    pl11 = int(round(-10 * l11))
    m = min(pl00, pl01, pl11)
    return (pl00 - m, pl01 - m, pl11 - m)


def _gt_from_ad_dp(ad_alt: int, dp: int) -> tuple[int | None, int | None]:
    """Naive GT call from a PL tuple. Uses argmin-PL heuristic."""
    if dp <= 0:
        return (None, None)
    p00, p01, p11 = _pl_from_ad_dp(ad_alt, dp)
    m = min(p00, p01, p11)
    if m == p01:
        return (0, 1)
    if m == p11:
        return (1, 1)
    return (0, 0)


def _write_real_imputed(subset_bcf: Path, panel_variants, out_vcf: Path,
                        imputed_panel_h5ad: Path) -> tuple[int, int]:
    """Copy subset records but override per-sample GT/AD/DP/PL from imputed panel."""
    import anndata as ad_mod
    from scipy.sparse import issparse
    a = ad_mod.read_h5ad(imputed_panel_h5ad)
    AD = a.layers["AD"].toarray() if issparse(a.layers["AD"]) else a.layers["AD"]
    DP = a.layers["DP"].toarray() if issparse(a.layers["DP"]) else a.layers["DP"]
    cell_idx = {c: i for i, c in enumerate(a.obs_names.astype(str))}
    # Panel var IDs like "I-140649-C>T" — build a lookup keyed by (chrom, pos, ref, alt)
    var_by_key: dict[tuple[str, int, str, str], int] = {}
    for j, row in a.var.iterrows():
        anc = str(row["anc"]); der = str(row["der"])
        ref = anc[1] if len(anc) == 3 else anc
        alt = der[1] if len(der) == 3 else der
        var_by_key[(str(row["chrom"]), int(row["pos"]), ref, alt)] = int(a.var_names.get_loc(j))

    vin = pysam.VariantFile(str(subset_bcf))
    vout = pysam.VariantFile(str(out_vcf), "wz", header=vin.header)
    n_seen = 0
    n_written = 0
    header_samples = list(vin.header.samples)
    for rec in vin:
        n_seen += 1
        matched_alt_idx = -1
        for i, alt in enumerate(rec.alts or ()):
            key = (rec.chrom, rec.pos, rec.ref, alt)
            if key in panel_variants and key in var_by_key:
                matched_alt_idx = i
                break
        if matched_alt_idx < 0:
            continue
        # Collapse multiallelic if needed
        if len(rec.alts) > 1:
            _keep_single_alt(rec, matched_alt_idx)
        key = (rec.chrom, rec.pos, rec.ref, rec.alts[0])
        v_i = var_by_key[key]
        # Override every sample
        for s_name in header_samples:
            c_i = cell_idx.get(s_name)
            if c_i is None:
                continue
            ad_alt = int(AD[c_i, v_i])
            dp = int(DP[c_i, v_i])
            ad_ref = max(0, dp - ad_alt)
            sample = rec.samples[s_name]
            sample["DP"] = dp
            sample["AD"] = (ad_ref, min(ad_alt, dp))
            sample["GT"] = _gt_from_ad_dp(ad_alt, dp)
            sample["PL"] = _pl_from_ad_dp(ad_alt, dp)
        vout.write(rec)
        n_written += 1
    vin.close()
    vout.close()
    return n_seen, n_written



def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-h5ad", type=Path, required=True)
    ap.add_argument("--joint-vcf", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--imputed-panel-h5ad", type=Path, default=None,
                    help="Panel h5ad with imputed AD/DP layers. If set, overrides GT/AD/PL "
                         "from this panel using standard diploid PL formula; joint VCF is "
                         "still used for site records/topology.")
    ap.add_argument("--skip-existing", action="store_true",
                    help="Skip a matrix if its .vcf.gz already exists.")
    return ap.parse_args()


def write_tabix(vcf_path: Path) -> None:
    """Tabix-index a bgzipped VCF, overwriting an existing index."""
    subprocess.run(["tabix", "-f", "-p", "vcf", str(vcf_path)], check=True)


def subset_joint_vcf(joint_vcf: Path, out_bcf: Path, panel_variants, cell_ids) -> None:
    """Subset joint VCF to (panel positions) x (this worm's cells).

    Positional subset with bcftools (fast, uses tabix index); exact-allele
    filtering happens in the pysam pass.
    """
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        regions_bed = td / "regions.bed"
        with open(regions_bed, "w") as fh:
            for chrom, pos, _ref, _alt in sorted(panel_variants):
                fh.write(f"{chrom}\t{pos-1}\t{pos}\n")
        samples_txt = td / "samples.txt"
        with open(samples_txt, "w") as fh:
            fh.write("\n".join(cell_ids) + "\n")
        cmd = [
            "bcftools", "view",
            "-R", str(regions_bed),
            "-S", str(samples_txt),
            "--force-samples",
            "-Ob", "-o", str(out_bcf),
            str(joint_vcf),
        ]
        subprocess.run(cmd, check=True)
    subprocess.run(["bcftools", "index", "-f", str(out_bcf)], check=True)


def load_panel(panel_h5ad: Path):
    """Return (cells, variants, worm).

    Panel var stores trinucleotide contexts (``anc``, ``der``); the actual
    ref/alt bases at ``pos`` are the middle characters. The VCF carries
    single-base ref/alt, so we match on those.
    """
    a = ad.read_h5ad(panel_h5ad)
    cells = list(a.obs_names.astype(str))
    variants = set()
    for _, row in a.var.iterrows():
        anc = str(row["anc"])
        der = str(row["der"])
        ref = anc[1] if len(anc) == 3 else anc
        alt = der[1] if len(der) == 3 else der
        variants.add((str(row["chrom"]), int(row["pos"]), ref, alt))
    worm = a.uns.get("panel_worm", panel_h5ad.stem.replace("worm_", ""))
    return cells, variants, str(worm)


def _write_real(subset_bcf: Path, panel_variants, out_vcf: Path) -> tuple[int, int]:
    """Copy subset records whose (chrom,pos,ref,alt) match a panel variant."""
    vin = pysam.VariantFile(str(subset_bcf))
    vout = pysam.VariantFile(str(out_vcf), "wz", header=vin.header)
    n_seen = 0
    n_written = 0
    for rec in vin:
        n_seen += 1
        for i, alt in enumerate(rec.alts or ()):
            key = (rec.chrom, rec.pos, rec.ref, alt)
            if key in panel_variants:
                # If multi-allelic, keep only the matching alt. This matters
                # for GT16 state-space assumptions downstream.
                if len(rec.alts) > 1:
                    _keep_single_alt(rec, i)
                vout.write(rec)
                n_written += 1
                break
    vin.close()
    vout.close()
    return n_seen, n_written


def _keep_single_alt(rec: pysam.VariantRecord, keep_idx: int) -> None:
    """Collapse a multiallelic record to the single alt at keep_idx.

    Only touches fields we consume downstream (GT, AD, PL, plus alts). Other
    format fields are left untouched; CellPhy uses PL alone in GL mode.
    """
    rec.alts = (rec.alts[keep_idx],)
    n_alt = 1
    n_g = 3  # diploid, biallelic
    for sample in rec.samples.values():
        # GT: re-map indices; keep original if 0 or matches keep_idx+1, else missing
        gt = sample.get("GT")
        if gt is not None:
            new_gt = tuple(
                None if a is None else (0 if a == 0 else (1 if a == keep_idx + 1 else None))
                for a in gt
            )
            sample["GT"] = new_gt
        # AD: keep (ref, alt_keep)
        ad_vals = sample.get("AD")
        if ad_vals is not None and len(ad_vals) > keep_idx + 1:
            sample["AD"] = (ad_vals[0], ad_vals[keep_idx + 1])
        # PL: keep the three entries for genotypes {0/0, 0/keep, keep/keep}
        # Diploid PL ordering is F(j,k) = k*(k+1)/2 + j for genotype j/k.
        pl_vals = sample.get("PL")
        if pl_vals is not None:
            k = keep_idx + 1
            i00 = 0
            i0k = k * (k + 1) // 2 + 0
            ikk = k * (k + 1) // 2 + k
            if max(i00, i0k, ikk) < len(pl_vals):
                sample["PL"] = (pl_vals[i00], pl_vals[i0k], pl_vals[ikk])


def _write_covmask(real_vcf: Path, out_vcf: Path) -> int:
    """Rewrite PL/GT/AD to a soft-0 wherever DP>0."""
    vin = pysam.VariantFile(str(real_vcf))
    vout = pysam.VariantFile(str(out_vcf), "wz", header=vin.header)
    n = 0
    for rec in vin:
        for sample in rec.samples.values():
            dp = sample.get("DP") or 0
            if dp > 0:
                sample["GT"] = (0, 0)
                sample["AD"] = (dp, 0)
                sample["PL"] = COVMASK_PL_COV
            else:
                sample["GT"] = (None, None)
                sample["AD"] = (0, 0)
                sample["PL"] = COVMASK_PL_UNCOV
        vout.write(rec)
        n += 1
    vin.close()
    vout.close()
    return n


def _write_shuffled(real_vcf: Path, out_vcf: Path, seed: int) -> int:
    """Per-variant, permute (GT,AD,PL) among DP>0 samples (G0.2 shuffle)."""
    rng = np.random.default_rng(seed)
    vin = pysam.VariantFile(str(real_vcf))
    vout = pysam.VariantFile(str(out_vcf), "wz", header=vin.header)
    sample_names = list(vin.header.samples)
    n = 0
    for rec in vin:
        cov_idx = [i for i, s in enumerate(sample_names)
                   if (rec.samples[s].get("DP") or 0) > 0]
        if len(cov_idx) < 2:
            vout.write(rec)
            n += 1
            continue
        payloads = [(
            tuple(rec.samples[sample_names[i]].get("GT") or (None, None)),
            tuple(rec.samples[sample_names[i]].get("AD") or (0, 0)),
            tuple(rec.samples[sample_names[i]].get("PL") or COVMASK_PL_UNCOV),
        ) for i in cov_idx]
        perm = rng.permutation(len(cov_idx))
        for dst_pos, src_pos in enumerate(perm):
            s = rec.samples[sample_names[cov_idx[dst_pos]]]
            gt, ad_v, pl = payloads[src_pos]
            s["GT"] = gt
            s["AD"] = ad_v
            s["PL"] = pl
        vout.write(rec)
        n += 1
    vin.close()
    vout.close()
    return n


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    cells, panel_variants, worm = load_panel(args.panel_h5ad)
    print(f"[prepare_vcf] worm={worm}  n_cells={len(cells)}  n_panel_vars={len(panel_variants)}")

    subset_bcf = args.out_dir / "subset.bcf"
    if args.skip_existing and subset_bcf.exists():
        print(f"[prepare_vcf] subset exists, skipping bcftools view")
    else:
        subset_joint_vcf(args.joint_vcf, subset_bcf, panel_variants, cells)
        print(f"[prepare_vcf] subset.bcf written")

    real = args.out_dir / "real.vcf.gz"
    covmask = args.out_dir / "covmask.vcf.gz"
    shuffled = args.out_dir / "shuffled.vcf.gz"

    if args.skip_existing and real.exists():
        print(f"[prepare_vcf] real exists, skipping")
    else:
        if args.imputed_panel_h5ad is not None:
            n_seen, n_written = _write_real_imputed(
                subset_bcf, panel_variants, real, args.imputed_panel_h5ad)
            print(f"[prepare_vcf] real (IMPUTED from {args.imputed_panel_h5ad.name}): "
                  f"{n_seen} subset records → {n_written} panel matches")
        else:
            n_seen, n_written = _write_real(subset_bcf, panel_variants, real)
            print(f"[prepare_vcf] real: {n_seen} subset records → {n_written} panel matches")
        write_tabix(real)

    if args.skip_existing and covmask.exists():
        print(f"[prepare_vcf] covmask exists, skipping")
    else:
        n = _write_covmask(real, covmask)
        write_tabix(covmask)
        print(f"[prepare_vcf] covmask: {n} records")

    if args.skip_existing and shuffled.exists():
        print(f"[prepare_vcf] shuffled exists, skipping")
    else:
        n = _write_shuffled(real, shuffled, args.seed)
        write_tabix(shuffled)
        print(f"[prepare_vcf] shuffled: {n} records (seed={args.seed})")


if __name__ == "__main__":
    main()
