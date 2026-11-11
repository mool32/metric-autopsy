"""The claim's estimate within the cell types of its domain, the types weighted equally, with its
uncertainty over the replicates (v0.4, decided 2026-10-18).

A *unit* is the cells of one replicate in one cell type (``composition.MIN_CELLS`` of them at least).

* **Nested design.** Within each stratum of ``within`` and each type, the difference of the two
  groups' mean unit values, groups[0] minus groups[1] (the effect test's convention); the stratum's
  difference is the mean of its types' differences with equal weights; the strata are combined with
  the effect test's precision weights 1 / (1/n_A + 1/n_B) over their replicates.
* **Paired design** (every replicate holds both groups). Each unit is a replicate's difference
  within a type, and the estimate is the mean over the types of their mean differences.

The standard error is cluster-robust over the replicates: a replicate that contributes to several
types is one unit of inference. With one type and one stratum it is Welch's (nested) or the paired
t's standard error; the degrees of freedom are Welch-Satterthwaite's over the (stratum, group)
components. The test relabels the replicates within the strata (sign flips for a paired design), as
the effect test does; a type left with an empty group under a relabeling drops out of that
relabeling's mean.

Every function works on a matrix of values, one row per statistic (the claim, or each of K
controls), one column per unit, so that the controls of the check of specificity share one pass.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np

from .stats import PermResult, TEstimate, nested_relabelings


@dataclass
class Units:
    """The units of a comparison and the bookkeeping the estimator needs.

    ``reps`` the replicates (in a fixed order), ``rep_first`` 1 for a replicate of groups[0] (nested),
    ``rep_stratum`` each replicate's stratum index; ``unit_rep`` and ``unit_type`` each unit's
    replicate and type indices; ``cells`` each unit's cell indices (nested) or the pair of index
    arrays of its two groups' cells (paired)."""

    kind: str
    reps: list
    rep_first: np.ndarray
    rep_stratum: np.ndarray
    strata: list
    types: list
    unit_rep: np.ndarray
    unit_type: np.ndarray
    cells: list = field(default_factory=list)

    @property
    def n_units(self) -> int:
        return int(len(self.unit_rep))

    def subset(self, keep: np.ndarray) -> "Units":
        """The units where `keep` is True (the replicates and strata are kept as they are)."""
        keep = np.asarray(keep, bool)
        return Units(self.kind, self.reps, self.rep_first, self.rep_stratum, self.strata, self.types,
                     self.unit_rep[keep], self.unit_type[keep], [c for c, k in zip(self.cells, keep) if k])


def _onehot(index: np.ndarray, n: int) -> np.ndarray:
    E = np.zeros((len(index), n))
    if len(index):
        E[np.arange(len(index)), index] = 1.0
    return E


class Estimator:
    """The combined estimate for one set of units; `estimate` takes a (C, U) matrix of unit values
    (for a paired design, the units' within-replicate differences) and the replicates' labels."""

    def __init__(self, units: Units):
        self.u = units
        R = len(units.reps)
        self.R = R
        self.unit_stratum = units.rep_stratum[units.unit_rep] if units.kind == "nested" else np.zeros(units.n_units, int)
        S = len(units.strata) if units.kind == "nested" else 1
        self.S = S
        K = len(units.types)
        # the cells (stratum, type) that hold units
        cell_of = self.unit_stratum * K + units.unit_type
        cells, inverse = np.unique(cell_of, return_inverse=True)
        self.cell_index = inverse
        self.J = len(cells)
        self.cell_stratum = (cells // K).astype(int)
        self.E = _onehot(inverse, self.J)               # unit -> cell
        self.F = _onehot(self.cell_stratum, S)          # cell -> stratum
        self.Rm = _onehot(units.unit_rep, R)            # unit -> replicate
        # the replicates that contribute (hold at least one unit), and the effect test's weights
        self.has_unit = self.Rm.sum(axis=0) > 0
        if units.kind == "nested":
            self.weights = np.zeros(S)
            for s in range(S):
                in_s = self.has_unit & (units.rep_stratum == s)
                n_a, n_b = int((in_s & (units.rep_first == 1)).sum()), int((in_s & (units.rep_first == 0)).sum())
                self.weights[s] = 1.0 / (1.0 / n_a + 1.0 / n_b) if n_a and n_b else 0.0

    # ------------------------------------------------------------------ nested
    def _nested(self, V: np.ndarray, g: np.ndarray, with_se: bool):
        M = np.isfinite(V)
        V0 = np.where(M, V, 0.0)
        gu = g[self.u.unit_rep].astype(float)
        S1, N1 = (V0 * gu) @ self.E, (M * gu) @ self.E
        S0, N0 = (V0 * (1.0 - gu)) @ self.E, (M * (1.0 - gu)) @ self.E
        with np.errstate(invalid="ignore", divide="ignore"):
            m1, m0 = S1 / N1, S0 / N0
        defined = (N1 > 0) & (N0 > 0)
        d = np.where(defined, m1 - m0, 0.0)
        den = defined.astype(float) @ self.F
        with np.errstate(invalid="ignore", divide="ignore"):
            Ds = (d @ self.F) / den
        ok = den > 0
        w = np.where(ok, self.weights[None, :], 0.0)
        W = w.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            D = np.where(W > 0, (np.where(ok, Ds, 0.0) * w).sum(axis=1) / W, np.nan)
        if not with_se:
            return D, None, None
        mean_u = np.where(gu[None, :] == 1, m1[:, self.cell_index], m0[:, self.cell_index])
        n_u = np.where(gu[None, :] == 1, N1[:, self.cell_index], N0[:, self.cell_index])
        k_u = den[:, self.cell_stratum[self.cell_index]]
        valid = M & defined[:, self.cell_index]
        with np.errstate(invalid="ignore", divide="ignore"):
            contrib = np.where(valid, (V - mean_u) / (k_u * n_u), 0.0)
        a = contrib @ self.Rm                                   # (C, R)
        var = np.zeros(V.shape[0])
        comps, dfs = [], []
        for s in range(self.S):
            for h in (1, 0):
                members = self.has_unit & (self.u.rep_stratum == s) & (g == h)
                m = int(members.sum())
                if m == 0:
                    continue
                scale = (w[:, s] / np.where(W > 0, W, 1.0)) ** 2
                if m < 2:
                    comp = np.full(V.shape[0], np.nan)
                else:
                    comp = scale * (m / (m - 1.0)) * (a[:, members] ** 2).sum(axis=1)
                comps.append(comp)
                dfs.append(m - 1)
        if not comps:
            return D, np.full(V.shape[0], np.nan), np.zeros(V.shape[0])
        C = np.array(comps)                                     # (n_comp, C)
        dfs = np.array(dfs, float)[:, None]
        var = C.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            df = var ** 2 / np.where(dfs > 0, C ** 2 / dfs, np.inf).sum(axis=0)
        return D, np.sqrt(var), df

    # ------------------------------------------------------------------ paired
    def _paired(self, V: np.ndarray, sign: np.ndarray, with_se: bool):
        M = np.isfinite(V)
        Vs = np.where(M, V * sign[self.u.unit_rep][None, :], 0.0)
        S, N = Vs @ self.E, M.astype(float) @ self.E          # cells are the types here
        with np.errstate(invalid="ignore", divide="ignore"):
            mk = S / N
        defined = N > 0
        K = defined.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            D = np.where(K > 0, np.where(defined, mk, 0.0).sum(axis=1) / K, np.nan)
        if not with_se:
            return D, None, None
        mean_u = mk[:, self.cell_index]
        n_u = N[:, self.cell_index]
        with np.errstate(invalid="ignore", divide="ignore"):
            contrib = np.where(M, (Vs - mean_u) / (K[:, None] * n_u), 0.0)
        a = contrib @ self.Rm
        m = int(self.has_unit.sum())
        if m < 2:
            return D, np.full(V.shape[0], np.nan), np.zeros(V.shape[0])
        var = (m / (m - 1.0)) * (a[:, self.has_unit] ** 2).sum(axis=1)
        return D, np.sqrt(var), np.full(V.shape[0], float(m - 1))

    def estimate(self, V: np.ndarray, labels: np.ndarray | None = None, with_se: bool = True):
        """(estimate, standard error, degrees of freedom), one per row of `V`. `labels` are the
        replicates' group labels (nested: 1 for groups[0]) or signs (paired); by default the observed."""
        V = np.atleast_2d(np.asarray(V, float))
        if self.u.kind == "nested":
            g = self.u.rep_first if labels is None else np.asarray(labels)
            return self._nested(V, g, with_se)
        sign = np.ones(self.R) if labels is None else np.asarray(labels, float)
        return self._paired(V, sign, with_se)

    # ------------------------------------------------------------------ relabelings
    def relabelings(self, *, n_perm: int, rng, max_exact: int) -> tuple[np.ndarray, str, float]:
        """The relabelings of the contributing replicates (the observed one among the exact ones), the
        method, and the smallest attainable two-sided p."""
        idx = np.where(self.has_unit)[0]
        if self.u.kind == "nested":
            g = self.u.rep_first[idx]
            s = self.u.rep_stratum[idx]
            sub, method = nested_relabelings(g, s, n_perm=n_perm, rng=rng, max_exact=max_exact)
            full = np.tile(self.u.rep_first, (len(sub), 1))
            full[:, idx] = sub
            strata = [st for st in np.unique(s) if (g[s == st] == 1).any() and (g[s == st] == 0).any()]
            balanced = all(2 * int((g[s == st] == 1).sum()) == int((s == st).sum()) for st in strata)
            total = len(sub)
            min_p = ((2.0 if balanced else 1.0) / total) if method == "exact" else 1.0 / (1 + n_perm)
            return full, method, min_p
        n = len(idx)
        if 2 ** n <= max_exact:
            flips = np.array(list(itertools.product((1.0, -1.0), repeat=n)))
            method = "exact"
            min_p = 2.0 / len(flips)
        else:
            flips = rng.choice((1.0, -1.0), size=(n_perm, n))
            method = "monte-carlo"
            min_p = 1.0 / (1 + n_perm)
        full = np.ones((len(flips), self.R))
        full[:, idx] = flips
        return full, method, min_p

    def permutation(self, V: np.ndarray, *, n_perm: int, rng, max_exact: int) -> list:
        """The two-sided relabeling test of each row of `V` (the PermResult of ``stats``)."""
        V = np.atleast_2d(np.asarray(V, float))
        obs, _, _ = self.estimate(V, with_se=False)
        labs, method, min_p = self.relabelings(n_perm=n_perm, rng=rng, max_exact=max_exact)
        null = np.array([self.estimate(V, lab, with_se=False)[0] for lab in labs])  # (L, C)
        out = []
        for c in range(V.shape[0]):
            o = obs[c]
            if not np.isfinite(o):
                out.append(PermResult(float("nan"), "none", 0, 0.0, float("nan"), float(o)))
                continue
            col = null[:, c]
            col = col[np.isfinite(col)]
            eps = 1e-12 * max(1.0, abs(o))
            count = int(np.sum(np.abs(col) >= abs(o) - eps))
            if method == "exact":
                p = count / len(col) if len(col) else float("nan")
                out.append(PermResult(float(p), "exact", int(len(col)), 0.0, float(min_p), float(o)))
            else:
                p = (1.0 + count) / (1.0 + len(col))
                out.append(PermResult(float(p), "monte-carlo", int(len(col)), float(np.sqrt(p * (1 - p) / max(len(col), 1))),
                                      float(min_p), float(o)))
        return out

    def counts(self) -> dict:
        """The contributing replicates per group (nested) or the pairs (paired)."""
        if self.u.kind == "paired":
            return dict(pairs=int(self.has_unit.sum()))
        return dict(first=int((self.has_unit & (self.u.rep_first == 1)).sum()),
                    second=int((self.has_unit & (self.u.rep_first == 0)).sum()))


def t_estimate(est: float, se: float, df: float, kind: str) -> TEstimate:
    method = "domain: cluster-robust over replicates" + (" (paired)" if kind == "paired" else "")
    return TEstimate(float(est), float(se) if np.isfinite(se) else float("nan"),
                     float(df) if np.isfinite(df) else 0.0, method)


def one_group_mean(R: np.ndarray, units: Units, group: int, weights: np.ndarray):
    """GATE 4's response on the domain: the mean over the injected group's units of their responses
    (`R`, (C, U) for the units of that group only; NaN where a unit gave none), combined as the effect
    is (types equally within a stratum, the strata by the effect test's `weights`), with its
    cluster-robust standard error over the group's replicates and Satterthwaite's degrees of freedom.
    Returns (mean, se, df, n_replicates)."""
    R = np.atleast_2d(np.asarray(R, float))
    M = np.isfinite(R)
    K = len(units.types)
    stratum = units.rep_stratum[units.unit_rep] if units.kind == "nested" else np.zeros(units.n_units, int)
    S = int(max(len(units.strata), 1)) if units.kind == "nested" else 1
    cell_of = stratum * K + units.unit_type
    cells, inv = np.unique(cell_of, return_inverse=True)
    E = _onehot(inv, len(cells))
    F = _onehot((cells // K).astype(int), S)
    Rm = _onehot(units.unit_rep, len(units.reps))
    Sum, N = np.where(M, R, 0.0) @ E, M.astype(float) @ E
    with np.errstate(invalid="ignore", divide="ignore"):
        mk = Sum / N
    defined = N > 0
    den = defined.astype(float) @ F
    with np.errstate(invalid="ignore", divide="ignore"):
        Ms = (np.where(defined, mk, 0.0) @ F) / den
    w = np.where(den > 0, np.asarray(weights, float)[None, :S], 0.0)
    W = w.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(W > 0, (np.where(den > 0, Ms, 0.0) * w).sum(axis=1) / W, np.nan)
        k_u = den[:, (cells // K).astype(int)[inv]]
        contrib = np.where(M, (R - mk[:, inv]) / (k_u * N[:, inv]), 0.0)
    a = contrib @ Rm
    rep_stratum = units.rep_stratum if units.kind == "nested" else np.zeros(len(units.reps), int)
    has = Rm.sum(axis=0) > 0
    comps, dfs = [], []
    for s in range(S):
        members = has & (rep_stratum == s)
        m = int(members.sum())
        if m == 0:
            continue
        scale = (w[:, s] / np.where(W > 0, W, 1.0)) ** 2
        comps.append(np.full(R.shape[0], np.nan) if m < 2 else scale * (m / (m - 1.0)) * (a[:, members] ** 2).sum(axis=1))
        dfs.append(m - 1)
    n_reps = int(has.sum())
    if not comps:
        return mean, np.full(R.shape[0], np.nan), np.zeros(R.shape[0]), n_reps
    Cm = np.array(comps)
    dfs = np.array(dfs, float)[:, None]
    var = Cm.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        df = var ** 2 / np.where(dfs > 0, Cm ** 2 / dfs, np.inf).sum(axis=0)
    return mean, np.sqrt(var), df, n_reps
