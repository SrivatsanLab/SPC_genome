#!/usr/bin/env python3
"""Assemble per-cell per-gene RNA expression estimates for a coassay run.

Model (simple excess):
    lambda[c] = intergenic_fragments[c] / intergenic_bp     # per-cell DNA baseline
    R_exp[c, g] = lambda[c] * exonic_length[g]
    R_rna[c, g] = max(0, R_obs[c, g] - R_exp[c, g])

Where R_obs is the featureCounts fragment count per (cell, gene), and
exonic_length is featureCounts's Length column (union of exons per gene).

Inputs:
    - featureCounts raw output (from create_rna_count_matrix.sh --all-reads)
    - intergenic BED (from make_intergenic_bed.py)
    - directory of per-cell intergenic-count TSVs (from count_intergenic_reads_array.sh)

Outputs at <output_prefix>:
    _matrix.csv     gene_id x barcode   float R_rna
    .h5ad           AnnData: X = R_rna (sparse float32), layers = raw_exon / expected,
                    obs = intergenic_fragments / lambda_intergenic,
                    var = exonic_length
"""

import argparse
import glob
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


def total_bed_bp(bed_path: Path) -> int:
    total = 0
    with bed_path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            total += int(fields[2]) - int(fields[1])
    return total


def load_per_cell_counts(per_cell_dir: Path) -> pd.DataFrame:
    rows = []
    for tsv in sorted(glob.glob(str(per_cell_dir / '*.tsv'))):
        with open(tsv) as fh:
            for line in fh:
                barcode, count = line.rstrip('\n').split('\t')
                rows.append((barcode, int(count)))
    if not rows:
        raise RuntimeError(f'No per-cell TSVs found in {per_cell_dir}')
    return pd.DataFrame(rows, columns=['barcode', 'intergenic_fragments']).set_index('barcode')


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument('featurecounts_raw', type=Path, help='featureCounts raw output (rna_counts_raw.txt)')
    p.add_argument('intergenic_bed', type=Path, help='intergenic BED (for total bp)')
    p.add_argument('per_cell_dir', type=Path, help='directory of {barcode}.tsv intergenic counts')
    p.add_argument('output_prefix', type=Path, help='output prefix (writes _matrix.csv and .h5ad)')
    args = p.parse_args()

    for f in (args.featurecounts_raw, args.intergenic_bed):
        if not f.is_file():
            print(f'Error: not found: {f}', file=sys.stderr)
            return 1
    if not args.per_cell_dir.is_dir():
        print(f'Error: per-cell dir not found: {args.per_cell_dir}', file=sys.stderr)
        return 1

    intergenic_bp = total_bed_bp(args.intergenic_bed)
    print(f'Intergenic bases: {intergenic_bp:,}')

    fc = pd.read_csv(args.featurecounts_raw, sep='\t', comment='#')
    gene_ids = fc['Geneid'].astype(str).to_numpy()
    exonic_length = fc['Length'].astype(np.int32).to_numpy()
    bam_cols = list(fc.columns[6:])
    barcodes = np.array([Path(c).stem for c in bam_cols])
    r_obs = fc[bam_cols].to_numpy(dtype=np.int32).T  # cells x genes
    print(f'featureCounts matrix: {r_obs.shape[0]} cells x {r_obs.shape[1]} genes')

    per_cell = load_per_cell_counts(args.per_cell_dir)
    per_cell['lambda_intergenic'] = per_cell['intergenic_fragments'] / intergenic_bp
    print(f'Per-cell intergenic counts: {len(per_cell)} cells')

    aligned = per_cell.reindex(barcodes)
    missing = aligned.index[aligned['intergenic_fragments'].isna()]
    if len(missing) > 0:
        print(f'Error: missing intergenic counts for {len(missing)} cells (first 5: {list(missing[:5])})', file=sys.stderr)
        return 1

    lam = aligned['lambda_intergenic'].to_numpy()

    r_exp = np.outer(lam, exonic_length).astype(np.float32)
    r_rna = np.maximum(0, r_obs.astype(np.float32) - r_exp)

    n_zero_lambda = int((lam == 0).sum())
    if n_zero_lambda:
        print(f'Warning: {n_zero_lambda} cells with zero intergenic fragments — R_rna == R_obs for those.')

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)

    csv_path = args.output_prefix.with_name(args.output_prefix.name + '_matrix.csv')
    pd.DataFrame(r_rna.T, columns=barcodes, index=pd.Index(gene_ids, name='gene_id')).to_csv(csv_path)
    print(f'Wrote {csv_path}')

    adata = ad.AnnData(
        X=csr_matrix(r_rna),
        obs=pd.DataFrame({
            'intergenic_fragments': aligned['intergenic_fragments'].to_numpy(),
            'lambda_intergenic':    lam,
        }, index=pd.Index(barcodes, name='barcode')),
        var=pd.DataFrame({
            'exonic_length': exonic_length,
        }, index=pd.Index(gene_ids, name='gene_id')),
    )
    adata.layers['raw_exon'] = csr_matrix(r_obs)
    adata.layers['expected'] = csr_matrix(r_exp)
    adata.uns['intergenic_bp'] = intergenic_bp
    adata.uns['model'] = 'simple_excess'  # R_rna = max(0, R_obs - lambda * exonic_length)

    h5ad_path = args.output_prefix.with_suffix('.h5ad')
    adata.write_h5ad(h5ad_path)
    print(f'Wrote {h5ad_path}')

    print()
    print('Summary:')
    print(f'  shape:                   {adata.n_obs} cells x {adata.n_vars} genes')
    print(f'  median lambda:           {np.median(lam):.3e} reads/bp')
    print(f'  mean R_obs sum / cell:   {r_obs.sum(axis=1).mean():.0f}')
    print(f'  mean R_exp sum / cell:   {r_exp.sum(axis=1).mean():.0f}')
    print(f'  mean R_rna sum / cell:   {r_rna.sum(axis=1).mean():.0f}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
