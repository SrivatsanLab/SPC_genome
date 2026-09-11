#!/usr/bin/env python
"""Export the worm6 co-assay per-cell counts for plotting.

Reproduces the data assembly in cells 13 and 16 of
notebooks/worm6_final_GEX.ipynb and writes the three columns the scatter
needs to CSV, so worm6_coassay_scatter.R can render it.

Unlike the UMAP, nothing here is stochastic: these are per-cell values already
stored in the two inputs, so the rendered figure reproduces the published one
exactly.

Usage:
    python export_worm6_coassay_counts.py <gex.h5ad> <dna_obs.csv> <out.csv>
"""

import sys

import anndata
import pandas as pd


def main(gex_path: str, dna_obs_path: str, out_path: str) -> None:
    print(f"reading {gex_path}", flush=True)
    adata = anndata.read_h5ad(gex_path)
    print(f"  {adata.n_obs} cells", flush=True)

    dna = pd.read_csv(dna_obs_path, index_col=0)
    missing = adata.obs_names.difference(dna.index)
    if len(missing):
        raise SystemExit(f"{len(missing)} GEX cells absent from {dna_obs_path}")

    obs = adata.obs.copy()
    # cell 13: the DNA metrics are joined onto the GEX obs by barcode
    obs["mean_coverage"] = dna.loc[obs.index, "mean_coverage"]

    for col in ("total_counts", "genes_detected"):
        if col not in obs.columns:
            raise SystemExit(f"missing obs['{col}'] in {gex_path}")

    # cell 16: sorted by total_counts descending with the top cell dropped
    out = (
        obs[["total_counts", "genes_detected", "mean_coverage"]]
        .sort_values("total_counts", ascending=False)
        .iloc[1:]
    )
    out.to_csv(out_path)
    print(f"wrote {len(out)} rows to {out_path} "
          f"(dropped the single highest total_counts cell, as the notebook does)",
          flush=True)
    print(out.describe().to_string(), flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
