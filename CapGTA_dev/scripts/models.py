"""Expression models for the baseline-calibration sweep.

All models share:

    fit(adata)      → self
    predict(adata)  → np.ndarray (n_cells, n_genes) of RNA count estimates

`adata` is the AnnData produced by build_calibration_table.py:
    X                 = exon_count
    layers['spliced'] = spliced_count
    layers['intron']  = intron_count
    layers['flank']   = flanking_count
    var:  exon_bp, intron_bp, flanking_bp, n_junc, chrom
    obs:  spliced_libsize

Read length is passed at construction (needed for p_junc).
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _dense(mat) -> np.ndarray:
    """Return a dense float32 numpy array from a numpy/scipy-sparse matrix."""
    if hasattr(mat, 'toarray'):
        return mat.toarray().astype(np.float32)
    return np.asarray(mat, dtype=np.float32)


def p_junc(n_junc: np.ndarray, exon_bp: np.ndarray, read_len: int) -> np.ndarray:
    """P(a random RNA read spans a splice junction), uniform-read model.

    ≈ n_junc * (read_len - 1) / exon_bp    (clipped to [0, 0.999])
    Returns a 1D array of length n_genes.
    """
    n_junc  = np.asarray(n_junc,  dtype=np.float32)
    exon_bp = np.asarray(exon_bp, dtype=np.float32)
    p = n_junc * (read_len - 1) / np.maximum(exon_bp, 1.0)
    return np.clip(p, 0.0, 0.999)


def _global_per_cell_rate(adata) -> np.ndarray:
    """Per-cell scalar gDNA rate: total (intron+flank) reads / total (intron+flank) bp."""
    intron_bp = adata.var['intron_bp'].to_numpy(dtype=np.float32)
    flank_bp  = adata.var['flanking_bp'].to_numpy(dtype=np.float32)
    intron_ct = _dense(adata.layers['intron'])
    flank_ct  = _dense(adata.layers['flank'])
    bg_reads = intron_ct.sum(axis=1) + flank_ct.sum(axis=1)  # (n_cells,)
    bg_bp    = intron_bp.sum() + flank_bp.sum()              # scalar
    return (bg_reads / max(bg_bp, 1.0)).astype(np.float32)   # (n_cells,)


def local_baseline_lambda(adata, source: str, floor_alpha: float = 0.0) -> np.ndarray:
    """Per-(cell, gene) gDNA rate (reads / bp) from a chosen local source.

    Returns (n_cells, n_genes) float32.

    source ∈ {'global', 'flanking'}
      global   — per-cell scalar rate, broadcast to all genes.
      flanking — flank_count / flanking_bp per (cell, gene), with an
                 optional global-floor rule:
                     λ[c,g] = max(flank_rate[c,g], floor_alpha · global_rate[c])
                 floor_alpha=0 → pure per-locus (zeros stay zero);
                 floor_alpha=1 → local rate never drops below per-cell rate.

    Intron / hybrid sources are removed — pre-mRNA contamination in
    introns breaks them (see sweep_v1_notes.md).
    """
    flank_bp = adata.var['flanking_bp'].to_numpy(dtype=np.float32)
    flank_ct = _dense(adata.layers['flank'])
    global_cell = _global_per_cell_rate(adata)  # (n_cells,)

    if source == 'global':
        return np.broadcast_to(global_cell[:, None], flank_ct.shape).astype(np.float32)

    if source == 'flanking':
        flank_rate = np.where(flank_bp[None, :] > 0,
                              flank_ct / np.maximum(flank_bp[None, :], 1.0),
                              0.0).astype(np.float32)
        if floor_alpha > 0:
            floor = (floor_alpha * global_cell)[:, None].astype(np.float32)
            return np.maximum(flank_rate, floor)
        return flank_rate

    raise ValueError(f'unknown baseline source: {source} (supported: global, flanking)')


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------

class ExcessModel:
    """Direct excess: max(0, exon_count - λ_local * exon_bp)."""
    def __init__(self, baseline: str = 'global', read_len: int = 100,
                 floor_alpha: float = 0.0, **_):
        self.baseline = baseline
        self.read_len = read_len
        self.floor_alpha = floor_alpha

    def fit(self, adata):
        return self

    def predict(self, adata) -> np.ndarray:
        lam = local_baseline_lambda(adata, self.baseline, self.floor_alpha)  # (C, G)
        expected = lam * adata.var['exon_bp'].to_numpy(dtype=np.float32)[None, :]
        return np.maximum(0.0, _dense(adata.X) - expected)


class NBGLMModel:
    """Placeholder — falls back to excess. Fill in composite likelihood later.

    Target model:
        exon_count[c,g]     ~ Poisson(λ_local[c,g] * L_g + μ_RNA[c,g])
        spliced_count[c,g]  ~ Binomial(μ_RNA[c,g], p_junc[g])           (n_junc[g] > 0 only)
        log μ_RNA[c,g]      = size_factor[c] + β · x_g
    """
    def __init__(self, baseline: str = 'global', read_len: int = 100,
                 floor_alpha: float = 0.0,
                 features: tuple[str, ...] = ('exon_bp', 'n_junc'), **_):
        self.baseline = baseline
        self.read_len = read_len
        self.floor_alpha = floor_alpha
        self.features = features
        self._params: dict | None = None

    def fit(self, adata):
        # TODO(dustin): fit composite likelihood on junction genes only.
        self._params = {'baseline': self.baseline}
        return self

    def predict(self, adata) -> np.ndarray:
        return ExcessModel(self.baseline, self.read_len, self.floor_alpha).predict(adata)


class GBDTModel:
    """Placeholder GBDT residual regressor. Falls back to excess for now.

    Target: upsampled spliced count = spliced_count / p_junc(gene).
    """
    def __init__(self, baseline: str = 'global', read_len: int = 100,
                 floor_alpha: float = 0.0,
                 features: tuple[str, ...] = (
                     'exon_count', 'expected_gdna', 'exon_bp', 'n_junc', 'spliced_libsize',
                 ), **_):
        self.baseline = baseline
        self.read_len = read_len
        self.floor_alpha = floor_alpha
        self.features = features
        self._model = None

    def fit(self, adata):
        # TODO(dustin): melt to long, subset n_junc > 0, fit XGBoost/LightGBM regressor.
        return self

    def predict(self, adata) -> np.ndarray:
        return ExcessModel(self.baseline, self.read_len, self.floor_alpha).predict(adata)


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

MODELS = {
    'excess':  ExcessModel,
    'nb_glm':  NBGLMModel,
    'gbdt':    GBDTModel,
}
