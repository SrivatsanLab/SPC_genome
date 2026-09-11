#!/usr/bin/env python
"""Export the worm6 haplotype UMAP embedding for plotting.

Reproduces cells 56-58 of notebooks/worm6_final_haplotype_assignment.ipynb and
writes the coordinates to CSV so worm6_haplotype_umap.R can render them. The
notebook itself only saves the rendered PNG, so without this the embedding
exists nowhere on disk.

The NMF/haplotype-deconvolution step upstream does not need re-running: the
notebook's last cell writes adata back to joint_variants_process.h5ad, so
obs['purity'], obs['donor'] and var['is_core'] are already persisted there.

Usage:
    python export_worm6_umap_coords.py <joint_variants_process.h5ad> <out.csv>
"""

import sys

import numpy as np
import pandas as pd
import anndata
from sklearn.decomposition import TruncatedSVD
import umap

# notebook constants
N_EMB = 18
SEED = 1
MIN_DP_VAR = 1
BG_PURITY = 0.4
DROP = ["worm02", "worm04", "worm14", "worm19"]


def main(h5ad_path: str, out_path: str) -> None:
    print(f"reading {h5ad_path}", flush=True)
    adata = anndata.read_h5ad(h5ad_path)
    print(f"  {adata.n_obs} obs x {adata.n_vars} vars", flush=True)

    for key, where in (("purity", adata.obs), ("donor", adata.obs)):
        if key not in where.columns:
            raise SystemExit(f"missing obs['{key}'] - was the notebook's final "
                             "write_h5ad cell run?")
    if "is_core" not in adata.var.columns:
        raise SystemExit("missing var['is_core'] - was the notebook's final "
                         "write_h5ad cell run?")

    purity_all = adata.obs["purity"].values
    donor_all = adata.obs["donor"].values
    background = purity_all < BG_PURITY

    keep = ~background & ~adata.obs["donor"].isin(DROP).values
    assert len(keep) == adata.n_obs, f"keep {len(keep)} vs {adata.n_obs}"
    print(f"keeping {keep.sum()} of {adata.n_obs} "
          f"({(~keep).sum()} dropped: background + {DROP})", flush=True)

    core_idx = np.flatnonzero(adata.var["is_core"].values)
    print(f"{len(core_idx)} core variants", flush=True)

    a = adata.layers["AD"][:, core_idx].toarray()[keep].astype(np.float32)
    d = adata.layers["DP"][:, core_idx].toarray()[keep].astype(np.float32)
    m = d >= MIN_DP_VAR
    with np.errstate(invalid="ignore", divide="ignore"):
        f = np.where(m, a / np.maximum(d, 1), np.nan)
    Xk = np.where(m, f - np.nanmean(f, 0), 0.0).astype(np.float32)
    del a, d, f

    print("TruncatedSVD ...", flush=True)
    emb = TruncatedSVD(n_components=N_EMB, random_state=SEED).fit_transform(Xk)
    print("UMAP ...", flush=True)
    U = umap.UMAP(n_neighbors=30, min_dist=0.3,
                  random_state=SEED).fit_transform(emb)

    # renumber the surviving worms to consecutive labels, as the notebook does
    labels = donor_all[keep]
    present = sorted(set(labels))
    worm_lbls = [w for w in present if w.startswith("worm")]
    renum = {w: f"worm{i + 1:02d}" for i, w in enumerate(worm_lbls)}
    renum.update({w: w for w in present if not w.startswith("worm")})
    call_k = np.array([renum[l] for l in labels])

    print(pd.Series(renum).to_string(), flush=True)

    out = pd.DataFrame({
        "UMAP1": U[:, 0],
        "UMAP2": U[:, 1],
        "donor": call_k,
        "purity": purity_all[keep],
    })
    out.to_csv(out_path, index=False)
    print(f"wrote {len(out)} rows to {out_path}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2])
