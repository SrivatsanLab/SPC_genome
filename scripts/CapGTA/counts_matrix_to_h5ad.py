#!/usr/bin/env python3
"""Convert a gene_id x barcode counts CSV (from create_rna_count_matrix.sh)
into a sparse .h5ad AnnData file for scanpy.

The AnnData follows the standard scRNA-seq layout:
    obs (rows) = cells (barcodes)
    var (cols) = genes

Usage:
    counts_matrix_to_h5ad.py <matrix.csv> <output.h5ad> [--summary <summary.csv>]

If --summary is passed and the file exists, the per-cell `total_counts` and
`genes_detected` columns are attached to obs.
"""

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


def build_anndata(matrix_csv: Path, summary_csv: Path | None) -> ad.AnnData:
    df = pd.read_csv(matrix_csv)
    if 'gene_id' not in df.columns:
        raise ValueError(f"{matrix_csv} is missing a 'gene_id' column")

    gene_ids = df['gene_id'].astype(str).values
    counts = df.drop(columns=['gene_id'])
    barcodes = counts.columns.astype(str).values

    # counts is genes x cells; AnnData wants cells x genes.
    X = csr_matrix(counts.to_numpy(dtype=np.int32).T)

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=pd.Index(barcodes, name='barcode')),
        var=pd.DataFrame(index=pd.Index(gene_ids, name='gene_id')),
    )

    if summary_csv is not None and summary_csv.exists():
        summary = pd.read_csv(summary_csv).set_index('barcode')
        # Align to obs_names, tolerate missing rows (fill NaN).
        aligned = summary.reindex(adata.obs_names)
        for col in ('total_counts', 'genes_detected'):
            if col in aligned.columns:
                adata.obs[col] = aligned[col].values

    return adata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('matrix_csv', type=Path, help='gene_id x barcode counts CSV')
    parser.add_argument('output_h5ad', type=Path, help='output .h5ad path')
    parser.add_argument('--summary', type=Path, default=None, help='optional per-cell summary CSV (barcode,total_counts,genes_detected)')
    args = parser.parse_args()

    if not args.matrix_csv.is_file():
        print(f"Error: matrix CSV not found: {args.matrix_csv}", file=sys.stderr)
        return 1

    adata = build_anndata(args.matrix_csv, args.summary)

    args.output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output_h5ad)

    print(f"Wrote {args.output_h5ad}")
    print(f"  shape:   {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"  X dtype: {adata.X.dtype} ({type(adata.X).__name__})")
    print(f"  obs:     {list(adata.obs.columns)}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
