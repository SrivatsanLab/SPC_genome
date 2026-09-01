#!/usr/bin/env python3
"""Sweep baseline × floor_alpha × model × features on calibration.h5ad.

Primary metric: **per-gene** Pearson/Spearman between predictions summed
over cells and the raw spliced-read count (or upsampled spliced count if
requested). Also breaks down per coverage tier using pct_1x from the QC
metrics table.

Cell-level metrics are still computed but reported as diagnostic-only —
Poisson thinning of spliced counts makes cell-level Pearson dominated by
noise, not model quality (see sweep_v1_notes.md).

Usage:
    run_sweep.py --calibration <calibration.h5ad> --output-dir <out> \\
        --qc-csv <pooled_sc_qc_metrics.csv> \\
        [--read-len 100] [--seed 0] [--test-frac 0.2]
"""

import argparse
import itertools
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.stats as ss
from scipy.sparse import csr_matrix

from models import MODELS, p_junc


# Baseline / floor_alpha grid.
# floor_alpha only matters for the flanking baseline; global is a scalar per cell.
CONFIGS = []
CONFIGS += [{'baseline': 'global',   'floor_alpha': 0.0}]
for a in [0.0, 0.25, 0.5, 1.0]:
    CONFIGS.append({'baseline': 'flanking', 'floor_alpha': a})

MODEL_CHOICES    = ['excess', 'nb_glm', 'gbdt']
FEATURE_CHOICES  = ['minimal', 'full']

FEATURE_SETS = {
    'minimal': ('exon_bp', 'n_junc'),
    'full':    ('exon_bp', 'n_junc', 'spliced_libsize', 'expected_gdna'),
}

COVERAGE_TIERS = [
    ('<5',   0.00, 0.05),
    ('5-20', 0.05, 0.20),
    ('20-50', 0.20, 0.50),
    ('50-80', 0.50, 0.80),
    ('>=80', 0.80, 1.01),
]


def split_test_genes(adata, seed: int, test_frac: float) -> np.ndarray:
    n_junc = adata.var['n_junc'].to_numpy()
    junc_idx = np.where(n_junc > 0)[0]
    rng = np.random.default_rng(seed)
    n_test = int(len(junc_idx) * test_frac)
    test_idx = rng.choice(junc_idx, size=n_test, replace=False)
    mask = np.zeros(adata.n_vars, dtype=bool)
    mask[test_idx] = True
    return mask


def _dense(m):
    return m.toarray() if hasattr(m, 'toarray') else np.asarray(m)


def per_gene_metrics(pred, cal, held_mask, cell_mask=None):
    """Gene-level metrics restricted to held-out junction genes.

    pred, cal.layers['spliced']: cells x genes matrices.
    Returns dict with per-gene Pearson (log1p) and Spearman on
    prediction-sum vs raw spliced-sum, plus a scale-invariant sanity.
    """
    S = _dense(cal.layers['spliced']).astype(np.float32)
    n_junc = cal.var['n_junc'].to_numpy()
    if cell_mask is None:
        cell_mask = np.ones(cal.n_obs, dtype=bool)
    if cell_mask.sum() == 0:
        return {'n_cells': 0}
    per_gene_pred = pred[cell_mask, :].sum(axis=0)
    per_gene_S    = S[cell_mask, :].sum(axis=0)
    gene_mask = held_mask & (n_junc > 0)
    if gene_mask.sum() < 3:
        return {'n_cells': int(cell_mask.sum()), 'n_genes': int(gene_mask.sum())}
    x = per_gene_pred[gene_mask]
    y = per_gene_S[gene_mask]
    try:
        pearson = float(np.corrcoef(np.log1p(x), np.log1p(y))[0, 1])
    except Exception:
        pearson = float('nan')
    try:
        spearman = float(ss.spearmanr(x, y).statistic)
    except Exception:
        spearman = float('nan')
    return {
        'n_cells':            int(cell_mask.sum()),
        'n_genes':            int(gene_mask.sum()),
        'pearson_log1p':      pearson,
        'spearman':           spearman,
        'mean_pred_per_gene': float(x.mean()),
        'mean_raw_spliced':   float(y.mean()),
    }


def per_cell_diagnostic(pred, cal, held_mask):
    """Diagnostic-only: cell-level Pearson (noisy due to Poisson thinning)."""
    S = _dense(cal.layers['spliced']).astype(np.float32)
    n_junc = cal.var['n_junc'].to_numpy()
    exon_bp = cal.var['exon_bp'].to_numpy().astype(np.float32)
    p = p_junc(n_junc, exon_bp, 100)
    with np.errstate(divide='ignore', invalid='ignore'):
        target = np.where(p[None, :] > 0, S / np.maximum(p[None, :], 1e-9), np.nan)
    m = held_mask
    ok = np.isfinite(target[:, m]) & np.isfinite(pred[:, m])
    if ok.sum() == 0:
        return {'cell_gene_pearson_log1p': float('nan')}
    p_flat = pred[:, m][ok]
    t_flat = target[:, m][ok]
    return {
        'cell_gene_pearson_log1p':  float(np.corrcoef(np.log1p(p_flat), np.log1p(t_flat))[0, 1]),
        'cell_gene_n':              int(ok.sum()),
    }


def run_one(cal, config, out_dir, seed, test_frac, coverage_tiers):
    out_dir.mkdir(parents=True, exist_ok=True)
    ModelClass = MODELS[config['model']]
    model = ModelClass(
        baseline=config['baseline'],
        floor_alpha=config['floor_alpha'],
        read_len=config['read_len'],
        features=FEATURE_SETS[config['features']],
    )

    test_mask = split_test_genes(cal, seed=seed, test_frac=test_frac)
    train_view = cal[:, ~test_mask]
    model.fit(train_view)
    preds = model.predict(cal).astype(np.float32)

    metrics = {'overall': per_gene_metrics(preds, cal, test_mask)}
    for tier_name, lo, hi in coverage_tiers:
        cm = ((cal.obs['pct_1x'] >= lo) & (cal.obs['pct_1x'] < hi)).to_numpy()
        metrics[f'tier_{tier_name}'] = per_gene_metrics(preds, cal, test_mask, cm)
    metrics['diag_cell_gene'] = per_cell_diagnostic(preds, cal, test_mask)

    # Save predictions as an h5ad with held_out mask.
    out = ad.AnnData(X=csr_matrix(preds), obs=cal.obs.copy(), var=cal.var.copy())
    out.var['held_out'] = test_mask
    out.uns['config'] = config
    out.write_h5ad(out_dir / 'predictions.h5ad')

    (out_dir / 'metrics.json').write_text(json.dumps(metrics, indent=2))
    (out_dir / 'config.json').write_text(json.dumps(config, indent=2))

    ovl = metrics['overall']
    print(f'  overall: pearson_log1p={ovl.get("pearson_log1p"):.3f}  '
          f'spearman={ovl.get("spearman"):.3f}  '
          f'mean_pred={ovl.get("mean_pred_per_gene"):.0f}  '
          f'mean_raw_S={ovl.get("mean_raw_spliced"):.0f}')


def attach_qc(cal, qc_csv):
    qc = pd.read_csv(qc_csv, usecols=['barcode', 'pct_1x', 'gini', 'mean_coverage'])
    qc = qc.set_index('barcode').reindex(cal.obs_names)
    for col in ['pct_1x', 'gini', 'mean_coverage']:
        cal.obs[col] = qc[col].to_numpy()
    n_missing = cal.obs['pct_1x'].isna().sum()
    if n_missing:
        print(f'WARNING: {n_missing} cells missing QC — filling pct_1x=0 for those.')
        cal.obs['pct_1x'] = cal.obs['pct_1x'].fillna(0.0)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--calibration', type=Path, required=True)
    ap.add_argument('--output-dir',  type=Path, required=True)
    ap.add_argument('--qc-csv',      type=Path, required=True)
    ap.add_argument('--read-len',    type=int, default=100)
    ap.add_argument('--seed',        type=int, default=0)
    ap.add_argument('--test-frac',   type=float, default=0.2)
    ap.add_argument('--only',        type=str, default=None,
                    help='JSON dict {baseline, floor_alpha, model, features} to run one config')
    args = ap.parse_args()

    print(f'Loading calibration: {args.calibration}')
    cal = ad.read_h5ad(args.calibration)
    print(f'  {cal.n_obs:,} cells x {cal.n_vars:,} genes')
    attach_qc(cal, args.qc_csv)

    if args.only:
        configs = [json.loads(args.only)]
    else:
        configs = [
            {**c, 'model': m, 'features': f, 'read_len': args.read_len}
            for c in CONFIGS
            for m in MODEL_CHOICES
            for f in FEATURE_CHOICES
        ]

    for c in configs:
        c.setdefault('read_len', args.read_len)
        tag = f'{c["baseline"]}_a{c["floor_alpha"]:.2f}__{c["model"]}__{c["features"]}'
        print(f'\n=== {tag} ===')
        run_one(cal, c, args.output_dir / tag, args.seed, args.test_frac, COVERAGE_TIERS)

    return 0


if __name__ == '__main__':
    sys.exit(main())
