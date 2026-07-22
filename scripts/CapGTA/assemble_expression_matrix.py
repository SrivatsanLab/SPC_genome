#!/usr/bin/env python3
"""Assemble per-cell per-gene RNA expression estimates for a coassay run.

Model (simple-excess, NoFeatures baseline):

    non_exonic_bp    = sum(chromosome lengths) - sum(featureCounts Length column)
    lambda[c]        = Unassigned_NoFeatures[c] / non_exonic_bp
    R_exp[c, g]      = lambda[c] * exonic_length[g]
    R_rna[c, g]      = max(0, R_obs[c, g] - R_exp[c, g])

The per-cell DNA baseline `lambda` comes from featureCounts' own
`Unassigned_NoFeatures` column (fragments that passed all filters but did
not overlap any exon feature by >= --fracOverlap). This matches R_obs's
fragment definition exactly and, importantly, includes intronic sequence
in the baseline — which is critical for WGA data where amplification
is biased toward gene bodies vs. intergenic.

Inputs:
    - featureCounts raw output (from create_rna_count_matrix.sh --all-reads)
    - the sibling .summary file
    - reference FASTA index (.fai) for the total genome size

Outputs at <output_prefix>:
    _matrix.csv     gene_id x barcode   float R_rna
    .h5ad           AnnData: X = R_rna (sparse float32), layers = raw_exon / expected,
                    obs = background_fragments / lambda_background,
                    var = exonic_length
"""

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


def total_genome_bp(fai_path: Path) -> int:
    total = 0
    with fai_path.open() as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            total += int(fields[1])
    return total


def load_no_features(summary_path: Path) -> pd.Series:
    """Return a Series indexed by BAM path, values = Unassigned_NoFeatures count."""
    df = pd.read_csv(summary_path, sep='\t')
    if 'Status' not in df.columns:
        raise ValueError(f'{summary_path}: first column is not "Status"')
    row = df[df['Status'] == 'Unassigned_NoFeatures']
    if row.empty:
        raise ValueError(f'{summary_path}: no Unassigned_NoFeatures row')
    return row.iloc[0].drop('Status').astype(int)


def main() -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument('featurecounts_raw', type=Path, help='featureCounts raw output (rna_counts_raw.txt)')
    p.add_argument('featurecounts_summary', type=Path, help='featureCounts .summary file')
    p.add_argument('fai', type=Path, help='reference FASTA index (.fai)')
    p.add_argument('output_prefix', type=Path, help='output prefix (writes _matrix.csv and .h5ad)')
    args = p.parse_args()

    for f in (args.featurecounts_raw, args.featurecounts_summary, args.fai):
        if not f.is_file():
            print(f'Error: not found: {f}', file=sys.stderr)
            return 1

    genome_bp = total_genome_bp(args.fai)

    fc = pd.read_csv(args.featurecounts_raw, sep='\t', comment='#')
    gene_ids = fc['Geneid'].astype(str).to_numpy()
    exonic_length = fc['Length'].astype(np.int32).to_numpy()
    bam_cols = list(fc.columns[6:])
    barcodes = np.array([Path(c).stem for c in bam_cols])
    r_obs = fc[bam_cols].to_numpy(dtype=np.int32).T  # cells x genes

    exonic_bp = int(exonic_length.sum())
    non_exonic_bp = genome_bp - exonic_bp
    if non_exonic_bp <= 0:
        print(f'Error: non-exonic bp <= 0 (genome={genome_bp:,}, exonic_sum={exonic_bp:,})', file=sys.stderr)
        return 1

    print(f'Genome:          {genome_bp:,} bp')
    print(f'Exonic (sum L):  {exonic_bp:,} bp')
    print(f'Non-exonic:      {non_exonic_bp:,} bp')
    print(f'R_obs matrix:    {r_obs.shape[0]} cells x {r_obs.shape[1]} genes')

    no_features = load_no_features(args.featurecounts_summary)
    # Reindex by BAM path to match fc columns exactly, then align to our barcode order.
    no_features = no_features.reindex(bam_cols)
    if no_features.isna().any():
        missing = [c for c, v in zip(bam_cols, no_features) if pd.isna(v)]
        print(f'Error: summary missing NoFeatures for {len(missing)} BAMs (first 5: {missing[:5]})', file=sys.stderr)
        return 1
    bg_fragments = no_features.to_numpy(dtype=np.int64)

    lam = bg_fragments / non_exonic_bp

    r_exp = np.outer(lam, exonic_length).astype(np.float32)
    r_rna = np.maximum(0, r_obs.astype(np.float32) - r_exp)

    n_zero = int((bg_fragments == 0).sum())
    if n_zero:
        print(f'Warning: {n_zero} cells with zero background fragments — R_rna == R_obs for those.')

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)

    csv_path = args.output_prefix.with_name(args.output_prefix.name + '_matrix.csv')
    pd.DataFrame(r_rna.T, columns=barcodes, index=pd.Index(gene_ids, name='gene_id')).to_csv(csv_path)
    print(f'Wrote {csv_path}')

    adata = ad.AnnData(
        X=csr_matrix(r_rna),
        obs=pd.DataFrame({
            'background_fragments': bg_fragments,
            'lambda_background':    lam,
        }, index=pd.Index(barcodes, name='barcode')),
        var=pd.DataFrame({
            'exonic_length': exonic_length,
        }, index=pd.Index(gene_ids, name='gene_id')),
    )
    adata.layers['raw_exon'] = csr_matrix(r_obs)
    adata.layers['expected'] = csr_matrix(r_exp)
    adata.uns['genome_bp'] = genome_bp
    adata.uns['exonic_bp'] = exonic_bp
    adata.uns['non_exonic_bp'] = non_exonic_bp
    adata.uns['model'] = 'noFeatures_excess'

    h5ad_path = args.output_prefix.with_suffix('.h5ad')
    adata.write_h5ad(h5ad_path)
    print(f'Wrote {h5ad_path}')

    print()
    print('Summary:')
    print(f'  shape:                   {adata.n_obs} cells x {adata.n_vars} genes')
    print(f'  median bg fragments:     {int(np.median(bg_fragments)):,}')
    print(f'  median lambda:           {np.median(lam):.3e} fragments/bp')
    print(f'  mean R_obs sum / cell:   {r_obs.sum(axis=1).mean():.0f}')
    print(f'  mean R_exp sum / cell:   {r_exp.sum(axis=1).mean():.0f}')
    print(f'  mean R_rna sum / cell:   {r_rna.sum(axis=1).mean():.0f}')
    print(f'  median R_rna / R_obs:    {np.median(r_rna.sum(axis=1) / np.maximum(r_obs.sum(axis=1), 1)):.3f}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
