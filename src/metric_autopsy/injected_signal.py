"""Injected signal: plant a known change in the construct, by binomial thinning of counts.

A metric earns a ``metric_validity`` PASS only by *responding* to signal (a positive
control, or an injected signal) as well as resisting nuisance. These constructors return
an ``inject(data, rng) -> data`` callable for ``run_autopsy(signal_test=...)``. Thinning
only removes molecules, so the injected data stay valid counts with real technical noise
(the approach of seqgendiff; Gerard 2020).

The coupling injector carries ``inject.sham``: the same thinning, of the same intensity, without
the signal. GATE 4 contrasts the metric on injected data with the metric on sham data. Contrasting
with the untouched data instead confuses the signal with the thinning noise: injecting coupling
into a pair that is already strongly coupled *lowers* its correlation against the original
(0.541 -> 0.524) while raising it against the sham (0.370 -> 0.524).

The module injector has no sham (v0.4 change 2): it lowers the module by thinning the module's own
genes, so the loss of counts is the signal itself, as it is for a real decrease of expression, and
the reference is the untouched data. Each injector names the direction in which it must move a
metric of its construct (``inject.direction``) and the genes it touches (``inject.genes``), and
rewrites only those genes' columns (``inject.columns(before, rng)``), which lets GATE 4 inject into
one private copy of each replicate's cells.

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
    """Non-negative and integer up to float error (the thinning rounds before it draws)."""
    X = np.asarray(X)
    if X.size == 0:
        return
    if X.dtype.kind in "iu":
        ok = X.min() >= 0
    else:
        r = np.round(X)
        ok = X.min() >= 0 and np.abs(X - r).max() <= 1e-6 * max(1.0, float(r.max()))
    if not ok:
        raise ValueError(f"{what} needs raw counts (non-negative integers)")


def module(genes: Sequence[str], fold: float = 2.0, frac: float = 1.0):
    """Lower the expression of a gene module ``fold``-times, by thinning the module's own genes to
    1/fold (binomial thinning: the counts stay real) in a random ``frac`` of the cells (all by
    default). No sham: the loss of the module's counts is the signal itself, as it is for a real
    decrease of expression, and the cells' depth falls only by the module's share of their counts
    (``inject.depth_loss(data)`` gives that share x (1 - 1/fold) x frac). A metric of the module's
    expression must decrease (``inject.direction``).

    v0.3 raised the module by thinning every other gene against a sham that thinned every gene; on
    shallow data the sham's loss of single counts then made half of the response (v0.4 change 2).
    """
    genes = list(genes)
    if fold <= 1:
        raise ValueError("fold must be > 1")
    if not 0 < frac <= 1:
        raise ValueError("frac must be in (0, 1]")

    def _columns(before, rng):
        before = np.round(np.asarray(before)).astype(np.int64)
        after = before.astype(float)
        rows = np.where(rng.random(before.shape[0]) < frac)[0] if frac < 1 else np.arange(before.shape[0])
        after[rows] = rng.binomial(before[rows], 1.0 / fold).astype(float)
        return after

    def inject(data, rng):
        sd = _materialize(data)
        idx = [unique_col_index(sd.var_names, g) for g in genes]
        before = sd.X[:, idx]
        _require_counts(before, "module injection")
        after = _columns(before, rng)
        sd.X[:, idx] = after
        obs = sd.obs
        if "total_counts" in obs:
            # gene-subset objects carry library size in obs: remove the counts thinned (a copy:
            # pandas' copy-on-write hands out read-only arrays)
            tot = np.array(obs["total_counts"], dtype=float)
            obs["total_counts"] = tot - (before - after).sum(axis=1)
        return SimpleData(sd.X, obs, sd.var_names)

    def depth_loss(data) -> float:
        """The expected share of the cells' counts the injection removes."""
        X = as_dense(data.X)
        idx = [unique_col_index(list(data.var_names), g) for g in genes]
        total = float(X.sum())
        return float(frac * (1.0 - 1.0 / fold) * X[:, idx].sum() / total) if total > 0 else 0.0

    inject.description = f"module({len(genes)} genes, 1/{fold:g} by thinning, frac={frac:g})"
    inject.genes = tuple(genes)
    inject.columns = _columns
    inject.direction = "decrease"
    inject.depth_loss = depth_loss
    inject.sham = None
    return inject


def coupling(gene_a: str, gene_b: str, strength: float = 1.0):
    """Induce co-variation between two genes: both are thinned with the same per-cell
    keep-probability ``sigmoid(strength * z)``, ``z ~ N(0, 1)`` (a shared factor). The sham
    thins each gene with its own, independent factor of the same distribution.

    The injection touches only the two genes: ``inject.genes`` names them and
    ``inject.columns(before, rng)`` (and ``inject.sham.columns``) returns their thinned counts
    from their counts, with the same random draws as the injection itself. GATE 4 and GATE 5 use
    it to rewrite the two columns of one private copy of the data instead of copying the whole
    matrix for every injection (``gates.injected_deltas``)."""
    genes = (gene_a, gene_b)

    def _columns(before, rng, shared: bool):
        n = before.shape[0]
        z = rng.normal(size=(n, 1)) if shared else rng.normal(size=(n, 2))
        p = 1.0 / (1.0 + np.exp(-strength * z))
        return rng.binomial(np.round(np.asarray(before)).astype(np.int64), np.broadcast_to(p, before.shape)).astype(float)

    def _thin(data, rng, shared: bool):
        sd = _materialize(data)
        cols = [unique_col_index(sd.var_names, g) for g in genes]
        before = sd.X[:, cols]
        _require_counts(before, "coupling injection")
        after = _columns(before, rng, shared)
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
    inject.strength = strength
    inject.direction = "increase"
    inject.genes = sham.genes = genes
    inject.columns = lambda before, rng: _columns(before, rng, True)
    sham.columns = lambda before, rng: _columns(before, rng, False)
    inject.sham = sham
    return inject


def describe(inject) -> str:
    return getattr(inject, "description", getattr(inject, "__name__", "custom injection"))
