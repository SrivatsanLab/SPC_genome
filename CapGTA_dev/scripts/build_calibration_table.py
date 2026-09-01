#!/usr/bin/env python3
"""Assemble the (cell, gene) modeling table for the baseline-calibration sweep.

Output: an AnnData written to <output>.h5ad, cells x genes, with

    X                  = exon_count       (sparse int32)  (all-reads featureCounts)
    layers['spliced']  = spliced_count    (sparse int32)  (spliced-only featureCounts, no -B)
    layers['polyat']   = polyat_count     (sparse int32)  (>=10 A/T run reads)
    layers['sl']       = sl_count         (sparse int32)  (SL1/SL2 seed-matched reads)
    layers['intron']   = intron_count     (sparse int32)  (from build_local_baselines.py)
    layers['flank']    = flanking_count   (sparse int32)  (from build_local_baselines.py)
    var:  exon_bp, intron_bp, flanking_bp, n_junc, chrom
    obs:  spliced_libsize, polyat_libsize, sl_libsize

Every "count" input is a featureCounts raw output file (--*-fc). Any of
the RNA-signal inputs (spliced/polyat/sl) may be omitted; missing layers
just don't appear in the output h5ad.

Usage:
    build_calibration_table.py \\
        --exonic-fc <all_reads_rna_counts.txt> \\
        --intron-fc <intron_counts.tsv> \\
        --flanking-fc <flanking_counts.tsv> \\
        --intron-lengths <intron_lengths.tsv> \\
        --flanking-lengths <flanking_lengths.tsv> \\
        --junctions-csv <gene_junctions.csv> \\
        --output <calibration.h5ad> \\
        [--spliced-fc <spliced_no_B_raw.txt>] \\
        [--polyat-fc <polyat_counts_raw.txt>] \\
        [--sl-fc <sl_counts_raw.txt>]
"""

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


def load_featurecounts(path: Path) -> tuple[pd.DataFrame, pd.Series]:
    """Return (counts_df: gene x barcode, length_series indexed by gene)."""
    fc = pd.read_csv(path, sep='\t', comment='#')
    gene_ids = fc['Geneid'].astype(str)
    length = fc.set_index('Geneid')['Length']
    length.name = 'Length'
    bam_cols = list(fc.columns[6:])
    barcodes = [Path(c).stem for c in bam_cols]
    counts = fc[bam_cols].astype(np.int32)
    counts.columns = barcodes
    counts.index = gene_ids
    counts.index.name = 'gene_id'
    return counts, length


def align_wide(dfs: list[pd.DataFrame], genes, barcodes) -> list[pd.DataFrame]:
    """Reindex every wide gene x barcode DataFrame to (genes, barcodes)."""
    return [d.reindex(index=genes, columns=barcodes, fill_value=0) for d in dfs]


def wide_to_cells_by_genes(df: pd.DataFrame) -> csr_matrix:
    """gene x barcode -> csr(cells x genes)."""
    return csr_matrix(df.to_numpy(dtype=np.int32).T)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--exonic-fc',        type=Path, required=True)
    ap.add_argument('--intron-fc',        type=Path, required=True)
    ap.add_argument('--flanking-fc',      type=Path, required=True)
    ap.add_argument('--intron-lengths',   type=Path, required=True)
    ap.add_argument('--flanking-lengths', type=Path, required=True)
    ap.add_argument('--junctions-csv',    type=Path, required=True)
    ap.add_argument('--output',           type=Path, required=True)
    ap.add_argument('--spliced-fc',       type=Path, default=None)
    ap.add_argument('--polyat-fc',        type=Path, default=None)
    ap.add_argument('--sl-fc',            type=Path, default=None)
    args = ap.parse_args()

    print('Loading exonic (all-reads) featureCounts…')
    exon_counts, exon_length = load_featurecounts(args.exonic_fc)
    print(f'  {exon_counts.shape[0]:,} genes x {exon_counts.shape[1]:,} cells')

    print('Loading intron baseline…')
    intron_counts, _ = load_featurecounts(args.intron_fc)
    intron_len = pd.read_csv(args.intron_lengths, sep='\t').set_index('gene_id')['total_bp']

    print('Loading flanking baseline…')
    flank_counts, _ = load_featurecounts(args.flanking_fc)
    flank_len = pd.read_csv(args.flanking_lengths, sep='\t').set_index('gene_id')['total_bp']

    print('Loading gene junctions…')
    junc = pd.read_csv(args.junctions_csv).set_index('gene_id')

    # Optional RNA-signal layers
    optional_layers = {}
    for name, path in [('spliced', args.spliced_fc),
                       ('polyat',  args.polyat_fc),
                       ('sl',      args.sl_fc)]:
        if path is not None:
            print(f'Loading {name} featureCounts…')
            counts, _ = load_featurecounts(path)
            optional_layers[name] = counts

    # Core layers (exon, intron, flank) must all agree on barcodes/genes.
    # RNA signal layers (spliced, polyat, sl) get reindexed to the core set,
    # filling missing cells (e.g. empty filtered BAMs) with zero.
    core_frames = [exon_counts, intron_counts, flank_counts]
    genes    = sorted(set.intersection(*(set(f.index) for f in core_frames)))
    barcodes = sorted(set.intersection(*(set(f.columns) for f in core_frames)))
    print(f'Core (exon∩intron∩flank): {len(genes):,} genes x {len(barcodes):,} cells')

    exon_counts, intron_counts, flank_counts = align_wide(
        [exon_counts, intron_counts, flank_counts], genes, barcodes
    )
    for name, m in optional_layers.items():
        aligned = align_wide([m], genes, barcodes)[0]
        n_missing_cells = m.shape[1] - len(set(m.columns) & set(barcodes))
        core_missing = len(barcodes) - len(set(m.columns) & set(barcodes))
        if core_missing:
            print(f'  {name}: {core_missing} core cells missing from input (filled 0)')
        optional_layers[name] = aligned

    print('Assembling AnnData…')
    X      = wide_to_cells_by_genes(exon_counts)
    intron = wide_to_cells_by_genes(intron_counts)
    flank  = wide_to_cells_by_genes(flank_counts)
    layers = {'intron': intron, 'flank': flank}
    for name, m in optional_layers.items():
        layers[name] = wide_to_cells_by_genes(m)

    var = pd.DataFrame({
        'exon_bp':     pd.Series(exon_length).reindex(genes).fillna(0).astype('int32').values,
        'intron_bp':   pd.Series(intron_len).reindex(genes).fillna(0).astype('int32').values,
        'flanking_bp': pd.Series(flank_len).reindex(genes).fillna(0).astype('int32').values,
        'n_junc':      pd.Series(junc['n_junctions']).reindex(genes).fillna(0).astype('int32').values,
        'chrom':       pd.Series(junc['chrom']).reindex(genes).fillna('').astype(str).values,
    }, index=pd.Index(genes, name='gene_id'))

    obs_data = {}
    for name in ('spliced', 'polyat', 'sl'):
        if name in layers:
            obs_data[f'{name}_libsize'] = np.asarray(layers[name].sum(axis=1)).ravel().astype(np.int32)
    obs = pd.DataFrame(obs_data, index=pd.Index(barcodes, name='barcode'))

    adata = ad.AnnData(X=X, obs=obs, var=var, layers=layers)
    adata.uns['baseline_sources'] = ['global', 'flanking']
    adata.uns['rna_signals'] = list(optional_layers.keys())

    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)
    print(f'Wrote {args.output}')
    print(f'  X (exon):  nnz={adata.X.nnz:,}')
    for name, m in layers.items():
        print(f'  {name:8s}   nnz={m.nnz:,}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
