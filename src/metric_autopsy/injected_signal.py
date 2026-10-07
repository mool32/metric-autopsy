"""Injected signal: plant a known change in the construct, by binomial thinning of counts.

A metric earns a ``metric_validity`` PASS only by *responding* to signal (a positive
control, or an injected signal) as well as resisting nuisance. These constructors return
an ``inject(data, rng) -> data`` callable for ``run_autopsy(signal_test=...)``. Thinning
only removes molecules, so the injected data stay valid counts with real technical noise
(the approach of seqgendiff; Gerard 2020).

Each injector carries ``inject.sham``: the same thinning, of the same intensity, without the
signal. GATE 4 contrasts the metric on injected data with the metric on sham data. Contrasting
with the untouched data instead confuses the signal with the thinning noise: injecting coupling
into a pair that is already strongly coupled *lowers* its correlation against the original
(0.541 -> 0.524) while raising it against the sham (0.370 -> 0.524).

Not to be confused with ERCC *spike-ins*, the external RNA standard used for content
estimands (see ``equalize``).
"""
from __future__ import annotations

from typing import Sequence

import numpy as np

from .core import SimpleData, as_dense, unique_col_index


def _materialize(data) -> SimpleData:
    if isinstance(data, SimpleData):
        return SimpleData(data.X.copy(), data.obs.copy(), data.var_names)
    return SimpleData(as_dense(data.X).copy(), data.obs.copy(), list(data.var_names))


def _require_counts(X: np.ndarray, what: str):
    if not (np.all(X >= 0) and np.allclose(X, np.round(X))):
        raise ValueError(f"{what} needs raw counts (non-negative integers)")


def module(genes: Sequence[str], fold: float = 2.0, frac: float = 0.3):
    """Raise the *relative* expression of a gene module ``fold``-times in a random ``frac``
    of cells, by thinning every other gene to 1/fold in those cells (a composition change).
    """
    genes = list(genes)
    if fold <= 1:
        raise ValueError("fold must be > 1")

    def inject(data, rng):
        sd = _materialize(data)
        _require_counts(sd.X, "module injection")
        idx = [unique_col_index(sd.var_names, g) for g in genes]
        other = np.setdiff1d(np.arange(sd.n_vars), idx)
        cells = np.where(rng.random(sd.n_obs) < frac)[0]
        X = sd.X
        before = X[np.ix_(cells, other)]
        after = rng.binomial(before.astype(np.int64), 1.0 / fold).astype(float)
        X[np.ix_(cells, other)] = after
        obs = sd.obs
        if "total_counts" in obs:
            # gene-subset objects carry library size in obs: remove the expected loss
            tot = np.asarray(obs["total_counts"], dtype=float)
            mod = X[cells][:, idx].sum(axis=1)
            tot[cells] = mod + (tot[cells] - mod) / fold
            obs["total_counts"] = tot
        return SimpleData(X, obs, sd.var_names)

    def sham(data, rng):
        """Thin *every* gene to 1/fold in the same share of cells: the same loss of depth,
        no change in composition."""
        sd = _materialize(data)
        _require_counts(sd.X, "module injection")
        cells = np.where(rng.random(sd.n_obs) < frac)[0]
        X = sd.X
        X[cells] = rng.binomial(X[cells].astype(np.int64), 1.0 / fold).astype(float)
        obs = sd.obs
        if "total_counts" in obs:
            tot = np.asarray(obs["total_counts"], dtype=float)
            tot[cells] = tot[cells] / fold
            obs["total_counts"] = tot
        return SimpleData(X, obs, sd.var_names)

    inject.description = f"module({len(genes)} genes, fold={fold}, frac={frac})"
    inject.sham = sham
    return inject


def coupling(gene_a: str, gene_b: str, strength: float = 1.0):
    """Induce co-variation between two genes: both are thinned with the same per-cell
    keep-probability ``sigmoid(strength * z)``, ``z ~ N(0, 1)`` (a shared factor). The sham
    thins each gene with its own, independent factor of the same distribution."""

    def _thin(data, rng, shared: bool):
        sd = _materialize(data)
        _require_counts(sd.X, "coupling injection")
        cols = [unique_col_index(sd.var_names, gene_a), unique_col_index(sd.var_names, gene_b)]
        z = rng.normal(size=(sd.n_obs, 1)) if shared else rng.normal(size=(sd.n_obs, 2))
        p = 1.0 / (1.0 + np.exp(-strength * z))
        before = sd.X[:, cols]
        after = rng.binomial(before.astype(np.int64), np.broadcast_to(p, before.shape)).astype(float)
        sd.X[:, cols] = after
        obs = sd.obs
        if "total_counts" in obs:
            obs["total_counts"] = (np.asarray(obs["total_counts"], dtype=float)
                                   - (before - after).sum(axis=1))
        return SimpleData(sd.X, obs, sd.var_names)

    def inject(data, rng):
        return _thin(data, rng, shared=True)

    def sham(data, rng):
        return _thin(data, rng, shared=False)

    inject.description = f"coupling({gene_a}, {gene_b}, strength={strength})"
    inject.sham = sham
    return inject


def describe(inject) -> str:
    return getattr(inject, "description", getattr(inject, "__name__", "custom injection"))
