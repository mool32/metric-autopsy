"""Injected signal: plant a known change in the construct, by binomial thinning of counts.

A metric earns a ``metric_validity`` PASS only by *responding* to signal (a positive
control, or an injected signal) as well as resisting nuisance. These constructors return
an ``inject(data, rng) -> data`` callable for ``run_autopsy(signal_test=...)``. Thinning
only removes molecules, so the injected data stay valid counts with real technical noise
(the approach of seqgendiff; Gerard 2020).

The coupling injectors (``coupling`` for a pair, ``coupling_set`` for a set of genes, v0.4 change 9)
carry ``inject.sham``: the same thinning, of the same intensity, without the signal. GATE 4
contrasts the metric on injected data with the metric on sham data. Contrasting
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

from .core import SimpleData, as_dense, canonical_csr, column_block, issparse, unique_col_index, with_columns


def _materialize(data) -> SimpleData:
    """A private SimpleData to inject into. A sparse X is not copied: it is never written in place, the
    injected columns go into a new matrix (``_replaced``; v0.4 change 10)."""
    if isinstance(data, SimpleData):
        return SimpleData(data.X if issparse(data.X) else data.X.copy(), data.obs.copy(), data.var_names)
    if issparse(data.X):
        return SimpleData(canonical_csr(data.X), data.obs.copy(), list(data.var_names))
    return SimpleData(as_dense(data.X).copy(), data.obs.copy(), list(data.var_names))


def _replaced(sd: SimpleData, cols, after):
    """`sd`'s matrix with the columns `cols` set to `after`: in place for a dense X (a private copy), a new
    matrix for a sparse one."""
    if issparse(sd.X):
        return with_columns(sd.X, cols, after)
    sd.X[:, cols] = after
    return sd.X


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
        before = column_block(sd.X, idx)
        _require_counts(before, "module injection")
        after = _columns(before, rng)
        X = _replaced(sd, idx, after)
        obs = sd.obs
        if "total_counts" in obs:
            # gene-subset objects carry library size in obs: remove the counts thinned (a copy:
            # pandas' copy-on-write hands out read-only arrays)
            tot = np.array(obs["total_counts"], dtype=float)
            obs["total_counts"] = tot - (before - after).sum(axis=1)
        return SimpleData(X, obs, sd.var_names)

    def depth_loss(data) -> float:
        """The expected share of the cells' counts the injection removes."""
        idx = [unique_col_index(list(data.var_names), g) for g in genes]
        if issparse(data.X):
            total, part = float(data.X.sum()), float(column_block(data.X, idx).sum())
        else:
            X = as_dense(data.X)
            total, part = float(X.sum()), X[:, idx].sum()
        return float(frac * (1.0 - 1.0 / fold) * part / total) if total > 0 else 0.0

    inject.description = f"module({len(genes)} genes, 1/{fold:g} by thinning, frac={frac:g})"
    inject.kind = "module"
    inject.fold, inject.frac = float(fold), float(frac)
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
        before = column_block(sd.X, cols)
        _require_counts(before, "coupling injection")
        after = _columns(before, rng, shared)
        X = _replaced(sd, cols, after)
        obs = sd.obs
        if "total_counts" in obs:
            obs["total_counts"] = (np.asarray(obs["total_counts"], dtype=float)
                                   - (before - after).sum(axis=1))
        return SimpleData(X, obs, sd.var_names)

    def inject(data, rng):
        return _thin(data, rng, shared=True)

    def sham(data, rng):
        return _thin(data, rng, shared=False)

    inject.description = f"coupling({gene_a}, {gene_b}, strength={strength:g})"
    inject.kind = "coupling"
    inject.strength = strength
    inject.direction = "increase"
    inject.genes = sham.genes = genes
    inject.columns = lambda before, rng: _columns(before, rng, True)
    sham.columns = lambda before, rng: _columns(before, rng, False)
    inject.sham = sham
    return inject


COUPLE = ("all", "half")


def coupling_set(genes: Sequence[str], strength: float = 1.0, couple: str = "all"):
    """Induce co-variation in a set of genes (v0.4 change 9): the genes the probe couples are thinned
    with the same per-cell keep-probability ``sigmoid(strength * z)``, ``z ~ N(0, 1)`` (a shared
    factor); the sham thins each gene of the set with its own, independent factor of the same
    distribution, so the injection and the sham remove the same expected counts.

    `couple` names the genes the factor couples, fixed by what the metric reads: ``"all"`` (the
    default) the whole set; ``"half"`` a random half of it (``len(genes) // 2`` genes, drawn anew at
    every injection; the other half is thinned as the sham's). A factor shared by every gene of the set
    cancels in what is relative within the set (shares, proportions), which ``"half"`` moves. The set
    of two with ``couple="all"`` is ``coupling``'s probe, draw for draw.

    Column-local, as ``coupling``: ``inject.genes`` names the set, ``inject.columns(before, rng)`` and
    ``inject.sham.columns`` return its thinned counts."""
    genes = tuple(genes)
    k = len(genes)
    if k < 2:
        raise ValueError("a set needs at least two genes")
    if len(set(genes)) != k:
        raise ValueError("the set's genes must be distinct")
    if couple not in COUPLE:
        raise ValueError(f"couple must be one of {COUPLE}, not {couple!r}")
    if couple == "half" and k < 4:
        raise ValueError("couple='half' needs at least 4 genes (a half of at least two)")

    def _columns(before, rng, shared: bool):
        n = before.shape[0]
        if shared and couple == "all":
            z = rng.normal(size=(n, 1))
        else:
            z = rng.normal(size=(n, k))
            if shared:  # one random half shares a factor; the other keeps its own (the sham's)
                z[:, rng.permutation(k)[: k // 2]] = rng.normal(size=(n, 1))
        p = 1.0 / (1.0 + np.exp(-strength * z))
        return rng.binomial(np.round(np.asarray(before)).astype(np.int64), np.broadcast_to(p, before.shape)).astype(float)

    def _thin(data, rng, shared: bool):
        sd = _materialize(data)
        cols = [unique_col_index(sd.var_names, g) for g in genes]
        before = column_block(sd.X, cols)
        _require_counts(before, "coupling injection")
        after = _columns(before, rng, shared)
        X = _replaced(sd, cols, after)
        obs = sd.obs
        if "total_counts" in obs:
            obs["total_counts"] = (np.asarray(obs["total_counts"], dtype=float)
                                   - (before - after).sum(axis=1))
        return SimpleData(X, obs, sd.var_names)

    def inject(data, rng):
        return _thin(data, rng, shared=True)

    def sham(data, rng):
        return _thin(data, rng, shared=False)

    inject.description = f"coupling_set({k} genes, couple={couple}, strength={strength:g})"
    inject.kind = "coupling_set"
    inject.strength = strength
    inject.couple = couple
    inject.direction = "increase"
    inject.genes = sham.genes = genes
    inject.columns = lambda before, rng: _columns(before, rng, True)
    sham.columns = lambda before, rng: _columns(before, rng, False)
    inject.sham = sham
    return inject


FRONT_END_PROBES = ("coupling", "module")


def named(kind: str | None, *, gene_pair=None, module_genes=None, module_fold: float = 2.0):
    """The probe a front end names (the CLI's ``--inject-signal``, the MCP tool's ``inject_signal``; v0.4 change
    11): ``"coupling"`` of the gene pair (none without a pair, as before), ``"module"`` of `module_genes` lowered
    `module_fold`-times by thinning (with a SESOI on the construct's scale GATE 4 doses it to exp(SESOI)-fold);
    None for no probe."""
    if not kind:
        return None
    if kind == "coupling":
        return coupling(*gene_pair) if gene_pair else None
    if kind == "module":
        if not module_genes:
            raise ValueError("the module probe needs its genes (--module-genes, module_genes)")
        return module(list(module_genes), fold=float(module_fold))
    raise ValueError(f"unknown probe {kind!r}: one of {FRONT_END_PROBES}")


def describe(inject) -> str:
    return getattr(inject, "description", getattr(inject, "__name__", "custom injection"))


# --------------------------------------------------------------------------- #
# the probe's dose from a SESOI on the construct's scale (v0.4 change 3)
# --------------------------------------------------------------------------- #
BIO_VAR_FLOOR = 0.01  # the smallest variance of log expected expression the coupling dose assumes (SD 0.1)


def log_sigmoid_variance(strength: float, n_nodes: int = 96) -> float:
    """Var(log sigmoid(strength x z)), z ~ N(0, 1): the variance that the shared factor of
    ``coupling`` adds to each gene's log expected expression (Gauss-Hermite quadrature)."""
    z, w = np.polynomial.hermite_e.hermegauss(n_nodes)
    w = w / w.sum()
    f = -np.logaddexp(0.0, -float(strength) * z)
    m = float((w * f).sum())
    return float(max((w * (f - m) ** 2).sum(), 0.0))


def biological_log_variance(counts, totals) -> float:
    """The variance of a gene's log expected expression across cells, the counting noise removed: a
    Poisson-lognormal method of moments on the counts scaled to the cells' mean depth,
    exp(var) - 1 = (Var(x) - mean(x) mean(1/s)) / mean(x)^2 with x = count / s and s = depth / mean depth.
    NaN without counts."""
    y = np.asarray(counts, float)
    t = np.asarray(totals, float)
    ok = t > 0
    y, t = y[ok], t[ok]
    if len(y) < 3:
        return float("nan")
    s = t / t.mean()
    x = y / s
    mu = float(x.mean())
    if mu <= 0:
        return float("nan")
    excess = (float(x.var(ddof=1)) - mu * float(np.mean(1.0 / s))) / mu ** 2
    return float(np.log1p(max(excess, 0.0)))


def _strength_for_variance(v: float) -> float:
    """The strength whose shared factor adds `v` to a gene's log expected expression
    (``log_sigmoid_variance`` inverted by bisection)."""
    a, b = 0.0, 1.0
    while log_sigmoid_variance(b) < v:
        b *= 2.0
    for _ in range(200):
        mid = 0.5 * (a + b)
        a, b = (mid, b) if log_sigmoid_variance(mid) < v else (a, mid)
    return float(0.5 * (a + b))


def coupling_strength_for(construct_change: float, var_a: float, var_b: float) -> dict:
    """The strength at which ``coupling``'s shared factor raises the correlation of the two genes' log
    expected expression by `construct_change` against its sham. The factor adds v = Var(log sigmoid(s z))
    to both genes' variances and to their covariance; the sham adds v to the variances only, so the
    correlation rises by v / sqrt((var_a + v)(var_b + v)) whatever the genes' own correlation. The
    variances are floored at ``BIO_VAR_FLOOR``."""
    S = float(construct_change)
    if not 0 < S < 1:
        raise ValueError(f"a change of correlation must lie in (0, 1), not {S!r}")
    va, vb = (max(float(v), BIO_VAR_FLOOR) if np.isfinite(v) else BIO_VAR_FLOOR for v in (var_a, var_b))

    def rise(v):
        return v / np.sqrt((va + v) * (vb + v))
    lo, hi = 0.0, 1.0
    while rise(hi) < S:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if rise(mid) < S else (lo, mid)
    v = 0.5 * (lo + hi)
    return dict(strength=_strength_for_variance(v), added_variance=float(v), var_a=float(va), var_b=float(vb),
                floored=[g for g, x in (("a", var_a), ("b", var_b)) if not np.isfinite(x) or x < BIO_VAR_FLOOR])


def coupling_set_strength_for(construct_change: float, variances: Sequence[float]) -> dict:
    """The strength at which ``coupling_set``'s shared factor raises the mean correlation of the coupled
    genes' log expected expression, over their pairs, by `construct_change` against its sham (v0.4 change
    9; for a set of two, ``coupling_strength_for``). The factor adds v = Var(log sigmoid(s z)) to every
    coupled gene's variance and to every pair's covariance, so pair (i, j) rises by w_i w_j with
    w_i = sqrt(v / (var_i + v)); the mean over the pairs is ((sum w)^2 - sum w^2) / (k (k - 1)). With
    ``couple="half"`` the coupled half is random, so the mean over all the set's pairs is the expected
    mean over the coupled ones. The variances are floored at ``BIO_VAR_FLOOR``."""
    S = float(construct_change)
    if not 0 < S < 1:
        raise ValueError(f"a change of correlation must lie in (0, 1), not {S!r}")
    raw = [float(x) if x is not None else float("nan") for x in variances]
    if len(raw) < 2:
        raise ValueError("a set needs at least two genes' variances")
    var = np.array([max(x, BIO_VAR_FLOOR) if np.isfinite(x) else BIO_VAR_FLOOR for x in raw])
    k = len(var)

    def rise(v):
        w = np.sqrt(v / (var + v))
        return float((w.sum() ** 2 - (w ** 2).sum()) / (k * (k - 1)))
    lo, hi = 0.0, 1.0
    while rise(hi) < S:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if rise(mid) < S else (lo, mid)
    v = 0.5 * (lo + hi)
    return dict(strength=_strength_for_variance(v), added_variance=float(v), variances=[float(x) for x in var],
                n_pairs=k * (k - 1) // 2,
                floored=[i for i, x in enumerate(raw) if not np.isfinite(x) or x < BIO_VAR_FLOOR])


def at_construct_change(inject, construct_change: float, *, var_a: float | None = None, var_b: float | None = None,
                        variances: Sequence[float] | None = None):
    """The probe of `inject`'s kind and genes whose construct change is `construct_change`, and a
    record of the dose: a module lowered exp(change)-fold (a change of its log expression), a coupling
    at the strength that raises the correlation of the pair's log expected expression by the change
    (given the genes' biological variances, ``biological_log_variance``), a set's coupling at the
    strength that raises the mean correlation over the pairs it couples by the change (`variances`, one
    per gene of the set; change 9). Returns (probe, info); the probe is None for an injector of
    another kind."""
    kind = getattr(inject, "kind", None)
    S = float(construct_change)
    if kind == "module":
        fold = float(np.exp(S))
        return module(list(inject.genes), fold=fold, frac=getattr(inject, "frac", 1.0)), dict(
            kind=kind, construct_change=S, scale="log expression (natural log)", fold=fold)
    if kind == "coupling":
        info = coupling_strength_for(S, var_a if var_a is not None else float("nan"),
                                     var_b if var_b is not None else float("nan"))
        return coupling(*inject.genes, strength=info["strength"]), dict(
            kind=kind, construct_change=S, scale="correlation of the pair's log expected expression", **info)
    if kind == "coupling_set":
        genes = list(inject.genes)
        info = coupling_set_strength_for(S, variances if variances is not None else [float("nan")] * len(genes))
        info["floored"] = [genes[i] for i in info["floored"]]
        couple = getattr(inject, "couple", "all")
        return coupling_set(genes, strength=info["strength"], couple=couple), dict(
            kind=kind, construct_change=S, couple=couple,
            scale=("mean correlation of the set's log expected expression over its pairs" if couple == "all" else
                   "mean correlation of log expected expression over the pairs of the half it couples"),
            **info)
    return None, dict(kind=kind, construct_change=S, note="no rule translates this injector's dose")
