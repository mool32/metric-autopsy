"""Reference metrics — worked examples of things people compute and then over-trust.

Every metric here has signature ``metric(data, *, ...params) -> float`` so that binding
the params with ``functools.partial`` yields the ``metric(data) -> float`` contract the
gates consume:

    >>> from functools import partial
    >>> m = partial(mi_3bin, gene_a="Smad3", gene_b="Col1a1")
    >>> m(data)   # -> float

`mi_3bin` is the antihero: mutual information with a zero/low/high discretization. It is
exactly the metric that failed 0/7 gates, and it is here so the harness can demonstrate a
metric *failing* — a validator that can only bless things is useless.
"""
from __future__ import annotations

import numpy as np

from .core import gene_column


def _to_3bin(x: np.ndarray) -> np.ndarray:
    """Zero / low / high discretization. Bin 0 is 'not detected' — the confound."""
    z = x <= 0
    out = np.zeros_like(x, dtype=int)
    nz = x[~z]
    if nz.size == 0:
        return out
    med = np.median(nz)
    out[~z] = np.where(nz > med, 2, 1)
    return out


def mi_3bin(data, *, gene_a: str, gene_b: str) -> float:
    """Mutual information (nats) between two genes under 0/low/high discretization.

    CONFOUNDED BY SPARSITY: bin 0 == 'not detected', so anything that changes the
    detection rate (library size, dropout, sequencing depth) moves this metric with no
    change in biology. This is the metric that motivated the whole project.
    """
    a = _to_3bin(gene_column(data, gene_a))
    b = _to_3bin(gene_column(data, gene_b))
    n = len(a)
    if n == 0:
        return 0.0
    mi = 0.0
    for ia in (0, 1, 2):
        pa = np.mean(a == ia)
        if pa == 0:
            continue
        for ib in (0, 1, 2):
            pb = np.mean(b == ib)
            if pb == 0:
                continue
            pab = np.mean((a == ia) & (b == ib))
            if pab > 0:
                mi += pab * np.log(pab / (pa * pb))
    return float(mi)


def pearson(data, *, gene_a: str, gene_b: str) -> float:
    """Plain Pearson correlation of two genes across cells (includes zeros)."""
    a = gene_column(data, gene_a).astype(float)
    b = gene_column(data, gene_b).astype(float)
    if a.std() == 0 or b.std() == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def codetected_spearman(data, *, gene_a: str, gene_b: str) -> float:
    """Spearman correlation using only cells where BOTH genes are detected.

    A deliberately more robust reference: by conditioning on co-detection it removes the
    zero-bin dependence that sinks `mi_3bin`, so it is the metric that *passes* GATE 0 in
    the tests. Not a universal fix — just an illustration that the confound is avoidable.
    """
    a = gene_column(data, gene_a).astype(float)
    b = gene_column(data, gene_b).astype(float)
    both = (a > 0) & (b > 0)
    if both.sum() < 3:
        return 0.0
    ra = _rankdata(a[both])
    rb = _rankdata(b[both])
    if ra.std() == 0 or rb.std() == 0:
        return 0.0
    return float(np.corrcoef(ra, rb)[0, 1])


def norm_pearson(data, *, gene_a: str, gene_b: str) -> float:
    """Library-normalized (CP10k + log1p) Pearson on co-detected cells.

    The robust reference: per-cell normalization divides out library size, so unlike
    `pearson`/`codetected_spearman` this is invariant to the per-cell depth nuisance, and
    conditioning on co-detection removes the zero-bin dependence. It is the metric that
    *passes* GATE 0 in the test suite — proof that the confounds are avoidable, not that
    all metrics are doomed.
    """
    from .core import as_dense, unique_col_index

    X = as_dense(data.X)
    var_names = list(data.var_names)
    ia, ib = unique_col_index(var_names, gene_a), unique_col_index(var_names, gene_b)
    raw_a, raw_b = X[:, ia], X[:, ib]
    both = (raw_a > 0) & (raw_b > 0)
    if both.sum() < 3:
        return 0.0
    tot = X.sum(axis=1)
    tot = np.where(tot == 0, 1.0, tot)
    a = np.log1p(raw_a / tot * 1e4)[both]
    b = np.log1p(raw_b / tot * 1e4)[both]
    if a.std() == 0 or b.std() == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def spectral_entropy(data) -> float:
    """Shannon entropy of the normalized eigenvalue spectrum of the gene-gene correlation.

    A dataset-level scalar (no genes to bind). CONFOUNDED BY VARIANCE STRUCTURE: inflate
    the variance of a few genes and the spectrum — hence this 'metric' — moves, with no
    change in biology. Included to exercise GATE 0's variance perturbation.
    """
    from .core import as_dense

    X = as_dense(data.X)
    # keep genes with nonzero variance
    v = X.var(axis=0)
    X = X[:, v > 0]
    if X.shape[1] < 2:
        return 0.0
    C = np.corrcoef(X, rowvar=False)
    C = np.nan_to_num(C)
    w = np.linalg.eigvalsh(C)
    w = np.clip(w, 0, None)
    s = w.sum()
    if s == 0:
        return 0.0
    p = w / s
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())


def _rankdata(x: np.ndarray) -> np.ndarray:
    """Average-rank of x (ties averaged), scipy-free."""
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    ranks[order] = np.arange(1, len(x) + 1)
    # average ties
    _, inv, counts = np.unique(x, return_inverse=True, return_counts=True)
    sums = np.zeros(len(counts))
    np.add.at(sums, inv, ranks)
    return (sums / counts)[inv]
