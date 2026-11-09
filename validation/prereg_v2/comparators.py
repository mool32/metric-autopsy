"""Validation v2's comparators, without the engine (v2.md sections 3.2 and 3.7; requirement 18):

(a) the test over cells: the claim's metric on all cells of each group, a two-sided permutation test that
    relabels cells (the group sizes kept); donors ignored. Reported, not compared.
(b) the test over donors without gates: the metric on each donor's cells (each donor's difference between
    the groups in a paired design), compared by the engine's graded rule coded apart — an exact
    permutation over donors within strata (Monte Carlo above `max_exact` relabelings) at four or more a
    group, the t interval at three, no verdict below — with a TOST for absence and EFFECT BELOW SESOI
    (v0.4 change 14), so that it makes the same kinds of claims as the engine (v2_requirements.md
    section 1). No correction, no gate.

The metrics are recomputed from per-cell vectors read once (`PerCell`): the reference metrics the engine
offers by name (`norm_pearson`, `mi_3bin`, `set_norm_pearson`, the arithmetic of `metrics.py` of v0.4),
the panel's marker metrics (the mean of log1p CP10k of a gene, and of its raw counts; v2.md section 3.2)
and the useless ones (a random number, a constant). `test_v2_w1.py` checks every metric and every test
against the engine, number for number; this file imports numpy and this directory's `num` only.

The difference is groups[0] - groups[1], the engine's convention: a positive difference is a decrease
from groups[0] to groups[1].
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field
from math import comb

import numpy as np

from num import t_ppf, t_two_sided_p

MIN_CODETECTED = 20          # metrics.MIN_CODETECTED of v0.4
FLOOR_POWER = 0.95           # effect.FLOOR_POWER of v0.4 (change 4)
PAIR_METRICS = ("norm_pearson", "mi_3bin")
SET_METRICS = ("set_norm_pearson",)
MARKER_METRICS = ("marker_log1p_cp10k", "marker_mean_counts")
USELESS_METRICS = ("random", "constant")
DIRECTIONS = ("increase", "decrease", "two-sided")


# --------------------------------------------------------------------------- #
# the metrics, from per-cell vectors
# --------------------------------------------------------------------------- #
def _column(X, j: int) -> np.ndarray:
    if hasattr(X, "tocsc") and not isinstance(X, np.ndarray):
        return np.asarray(X[:, j].toarray(), dtype=float).ravel()
    return np.asarray(X[:, j], dtype=float).ravel()


def row_sums(X) -> np.ndarray:
    """Every cell's total over all genes: the engine's `core.row_sums` (a sparse matrix of integer
    counts summed on its entries, other sparse values row by row as numpy sums the dense rows)."""
    if hasattr(X, "tocsr") and not isinstance(X, np.ndarray):
        d = X.data
        if np.all(d == np.round(d)):
            return np.asarray(X.sum(axis=1), dtype=float).ravel()
        out = np.empty(X.shape[0], dtype=float)
        for r0 in range(0, X.shape[0], 256):
            out[r0:r0 + 256] = np.asarray(X[r0:r0 + 256].toarray(), dtype=float).sum(axis=1)
        return out
    return np.asarray(X, dtype=float).sum(axis=1)


@dataclass
class PerCell:
    """What a claim's metric reads from each cell, read once: the raw counts of the claim's genes
    (cells x genes) and the cells' totals over all genes."""

    counts: np.ndarray
    totals: np.ndarray

    @classmethod
    def from_matrix(cls, X, var_names, genes) -> "PerCell":
        names = list(var_names)
        idx = []
        for g in genes:
            hits = [i for i, v in enumerate(names) if v == g]
            if len(hits) != 1:
                raise ValueError(f"gene {g!r} is in the data {len(hits)} times")
            idx.append(hits[0])
        cols = np.column_stack([_column(X, j) for j in idx]) if idx else np.zeros((X.shape[0], 0))
        return cls(cols, row_sums(X))

    def subset(self, mask) -> "PerCell":
        return PerCell(self.counts[mask], self.totals[mask])


def _cp10k_log1p(raw: np.ndarray, tot: np.ndarray) -> np.ndarray:
    tot = np.where(tot == 0, 1.0, tot)
    return np.log1p(raw / tot * 1e4)


def _norm_pearson(raw_a, raw_b, tot, min_codetected: int = MIN_CODETECTED) -> float:
    both = (raw_a > 0) & (raw_b > 0)
    if both.sum() < max(3, int(min_codetected)):
        return float("nan")
    a = _cp10k_log1p(raw_a, tot)[both]
    b = _cp10k_log1p(raw_b, tot)[both]
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _to_3bin(x: np.ndarray) -> np.ndarray:
    z = x <= 0
    out = np.zeros_like(x, dtype=int)
    nz = x[~z]
    if nz.size == 0:
        return out
    med = np.median(nz)
    out[~z] = np.where(nz > med, 2, 1)
    return out


def _mi_3bin(raw_a, raw_b) -> float:
    a, b = _to_3bin(np.asarray(raw_a, float)), _to_3bin(np.asarray(raw_b, float))
    if len(a) == 0:
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


@dataclass
class Metric:
    """A claim's metric by name, computed on a `PerCell` (a pair's two genes, a set's genes, a marker's
    one gene). `random` draws from its own stream (keyed by the caller), `constant` is 0."""

    name: str
    rng: np.random.Generator | None = None
    min_codetected: int = MIN_CODETECTED

    def __post_init__(self):
        known = PAIR_METRICS + SET_METRICS + MARKER_METRICS + USELESS_METRICS
        if self.name not in known:
            raise ValueError(f"unknown metric {self.name!r}")
        if self.name == "random" and self.rng is None:
            raise ValueError("the random metric needs its stream (rng)")

    def __call__(self, cells: PerCell) -> float:
        c, tot = cells.counts, cells.totals
        if self.name == "norm_pearson":
            return _norm_pearson(c[:, 0], c[:, 1], tot, self.min_codetected)
        if self.name == "mi_3bin":
            return _mi_3bin(c[:, 0], c[:, 1])
        if self.name == "set_norm_pearson":
            if c.shape[1] < 2:
                raise ValueError("a set needs at least two genes")
            values = []
            for i, j in itertools.combinations(range(c.shape[1]), 2):
                v = _norm_pearson(c[:, i], c[:, j], tot, self.min_codetected)
                if not np.isfinite(v):
                    return float("nan")
                values.append(v)
            return float(np.mean(values))
        if self.name == "marker_log1p_cp10k":
            return float(np.mean(_cp10k_log1p(c[:, 0], tot))) if len(tot) else float("nan")
        if self.name == "marker_mean_counts":
            return float(np.mean(c[:, 0])) if len(tot) else float("nan")
        if self.name == "random":
            return float(self.rng.standard_normal())
        return 0.0


# --------------------------------------------------------------------------- #
# (a) the test over cells
# --------------------------------------------------------------------------- #
def _direction_call(est: float, direction: str) -> str:
    """'effect' or 'against' for a detected difference (groups[0] - groups[1]) and the claim's direction."""
    if direction == "two-sided":
        return "effect"
    observed = "decrease" if est > 0 else "increase"
    return "effect" if observed == direction else "against"


def cell_test(metric: Metric, cells: PerCell, in_first, *, direction: str, alpha: float = 0.05,
              n_perm: int = 999, rng: np.random.Generator) -> dict:
    """(a): the metric on the cells of groups[0] (`in_first`) and of groups[1], a two-sided permutation
    test that relabels the cells, (1 + count) / (1 + n_perm); an effect or one against the claim where
    p < alpha, else no claim (it makes no claim of absence)."""
    if direction not in DIRECTIONS:
        raise ValueError(f"direction must be one of {DIRECTIONS}")
    in_first = np.asarray(in_first, bool)
    obs = metric(cells.subset(in_first)) - metric(cells.subset(~in_first))
    if not np.isfinite(obs):
        return dict(call="none", est=float("nan"), p=float("nan"), n_perm=0)
    count = 0
    lab = in_first.copy()
    eps = 1e-12 * max(1.0, abs(obs))
    for _ in range(n_perm):
        rng.shuffle(lab)
        d = metric(cells.subset(lab)) - metric(cells.subset(~lab))
        if np.isfinite(d) and abs(d) >= abs(obs) - eps:
            count += 1
    p = (1.0 + count) / (1.0 + n_perm)
    call = _direction_call(obs, direction) if p < alpha else "none"
    return dict(call=call, est=float(obs), p=float(p), n_perm=n_perm)


# --------------------------------------------------------------------------- #
# (b) the test over donors: the engine's graded rule, coded apart
# --------------------------------------------------------------------------- #
@dataclass
class TEst:
    est: float
    se: float
    df: float
    method: str

    def ci(self, level: float) -> tuple[float, float]:
        if not (np.isfinite(self.se) and self.df > 0):
            return (float("nan"), float("nan"))
        q = t_ppf(0.5 + level / 2.0, self.df)
        return (self.est - q * self.se, self.est + q * self.se)

    @property
    def p(self) -> float:
        return t_two_sided_p(self.est, self.se, self.df)


@dataclass
class Perm:
    p: float
    method: str
    n: int
    min_attainable_p: float
    observed: float


def _weights(g: np.ndarray, s: np.ndarray) -> dict:
    out = {}
    for st in np.unique(s):
        n_a = int(np.sum((s == st) & (g == 1)))
        n_b = int(np.sum((s == st) & (g == 0)))
        if n_a and n_b:
            out[st] = 1.0 / (1.0 / n_a + 1.0 / n_b)
    return out


def _sdiff(y: np.ndarray, g: np.ndarray, s: np.ndarray, weights: dict) -> float:
    if not weights:
        return float("nan")
    num = 0.0
    for st, w in weights.items():
        m = s == st
        num += w * (y[m & (g == 1)].mean() - y[m & (g == 0)].mean())
    return float(num / sum(weights.values()))


def _nested_t(y, g, s) -> TEst:
    weights = _weights(g, s)
    est = _sdiff(y, g, s, weights)
    if len(weights) == 1:
        st = next(iter(weights))
        a, b = y[(s == st) & (g == 1)], y[(s == st) & (g == 0)]
        if len(a) < 2 or len(b) < 2:
            return TEst(est, float("nan"), 0.0, "welch")
        va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
        se = math.sqrt(va + vb)
        if se == 0:
            return TEst(est, 0.0, float(len(a) + len(b) - 2), "welch")
        df = (va + vb) ** 2 / ((va ** 2) / (len(a) - 1) + (vb ** 2) / (len(b) - 1)) if (va or vb) else 1.0
        return TEst(est, se, float(df), "welch")
    ss, df = 0.0, 0
    for st in weights:
        for grp in (0, 1):
            v = y[(s == st) & (g == grp)]
            if len(v) >= 2:
                ss += float(((v - v.mean()) ** 2).sum())
                df += len(v) - 1
    if df == 0:
        return TEst(est, float("nan"), 0.0, "pooled-stratified")
    se = math.sqrt(ss / df / sum(weights.values()))
    return TEst(est, se, float(df), "pooled-stratified")


def _paired_t(d) -> TEst:
    d = np.asarray(d, float)
    if len(d) < 2:
        return TEst(float(d.mean()) if len(d) else float("nan"), float("nan"), 0.0, "paired")
    return TEst(float(d.mean()), float(d.std(ddof=1) / math.sqrt(len(d))), float(len(d) - 1), "paired")


def _outer_sums(per_stratum: list) -> np.ndarray:
    acc = np.zeros(1)
    for arr in per_stratum:
        acc = (acc[:, None] + arr[None, :]).ravel()
    return acc


def _perm_nested(y, g, s, *, n_perm, rng, max_exact) -> Perm:
    weights = _weights(g, s)
    obs = _sdiff(y, g, s, weights)
    if not weights:
        return Perm(float("nan"), "none", 0, float("nan"), obs)
    strata = list(weights)
    total = 1
    for st in strata:
        m = s == st
        total *= comb(int(m.sum()), int((g[m] == 1).sum()))
    balanced = all(2 * int((g[s == st] == 1).sum()) == int((s == st).sum()) for st in strata)
    wsum = sum(weights.values())
    eps = 1e-12 * max(1.0, abs(obs))
    if total <= max_exact:
        per = []
        for st in strata:
            idx = np.where(s == st)[0]
            k = int((g[idx] == 1).sum())
            vals = y[idx]
            tot_sum, n_s = vals.sum(), len(idx)
            diffs = []
            for c in itertools.combinations(range(n_s), k):
                a_sum = vals[list(c)].sum()
                diffs.append(a_sum / k - (tot_sum - a_sum) / (n_s - k))
            per.append(weights[st] / wsum * np.asarray(diffs))
        stats = _outer_sums(per)
        p = float(np.mean(np.abs(stats) >= abs(obs) - eps))
        return Perm(p, "exact", int(total), (2.0 if balanced else 1.0) / total, obs)
    count = 0
    g_perm = g.copy()
    for _ in range(n_perm):
        for st in strata:
            idx = np.where(s == st)[0]
            g_perm[idx] = rng.permutation(g[idx])
        if abs(_sdiff(y, g_perm, s, weights)) >= abs(obs) - eps:
            count += 1
    p = (1.0 + count) / (1.0 + n_perm)
    return Perm(p, "monte-carlo", int(n_perm), 1.0 / (1 + n_perm), obs)


def _perm_paired(d, *, n_perm, rng, max_exact) -> Perm:
    d = np.asarray(d, float)
    obs = float(d.mean()) if len(d) else float("nan")
    n = len(d)
    if n == 0:
        return Perm(float("nan"), "none", 0, float("nan"), obs)
    eps = 1e-12 * max(1.0, abs(obs))
    total = 2 ** n
    if total <= max_exact:
        signs = np.array(list(itertools.product((1.0, -1.0), repeat=n)))
        stats = (signs * d[None, :]).mean(axis=1)
        p = float(np.mean(np.abs(stats) >= abs(obs) - eps))
        return Perm(p, "exact", int(total), 2.0 / total, obs)
    flips = rng.choice((1.0, -1.0), size=(n_perm, n))
    stats = (flips * d[None, :]).mean(axis=1)
    count = int(np.sum(np.abs(stats) >= abs(obs) - eps))
    p = (1.0 + count) / (1.0 + n_perm)
    return Perm(p, "monte-carlo", int(n_perm), 1.0 / (1 + n_perm), obs)


def design_floor(se: float, df: float, alpha: float, power: float = FLOOR_POWER) -> float:
    """The design's floor (v0.4 change 4): se x (t(1 - alpha, df) + t((1 + power) / 2, df))."""
    if not (np.isfinite(se) and se > 0 and df > 0):
        return float("nan")
    return float(se * (t_ppf(1.0 - alpha, df) + t_ppf((1.0 + power) / 2.0, df)))


@dataclass
class Design:
    """The replicates of a card: ``kind`` "nested" (each replicate in one group; ``first`` names those of
    groups[0]) or "paired" (each replicate's value is its difference groups[0] - groups[1]); ``reps`` in
    the order the engine forms them; ``stratum`` each replicate's stratum (one stratum if omitted)."""

    kind: str
    reps: list
    first: set = field(default_factory=set)
    stratum: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.kind not in ("nested", "paired"):
            raise ValueError("a design is 'nested' or 'paired'")


def _graded(values: dict, design: Design, *, alpha, n_perm, max_exact, min_replicates, rng) -> dict:
    """The engine's `_infer` and `_tiered` on replicate values (NaN or None left out)."""
    reps = [r for r in design.reps if values.get(r) is not None and np.isfinite(values[r])]
    if design.kind == "paired":
        d = np.array([values[r] for r in reps], float)
        te = _paired_t(d)
        perm = _perm_paired(d, n_perm=n_perm, rng=rng, max_exact=max_exact)
        n_eff = len(d)
        counts = dict(pairs=n_eff)
    else:
        g = np.array([1 if r in design.first else 0 for r in reps])
        labels = {s: i for i, s in enumerate(dict.fromkeys(design.stratum.get(r, "") for r in reps))}
        s = np.array([labels[design.stratum.get(r, "")] for r in reps])
        y = np.array([values[r] for r in reps], float)
        te = _nested_t(y, g, s)
        perm = _perm_nested(y, g, s, n_perm=n_perm, rng=rng, max_exact=max_exact)
        counts = dict(first=int(g.sum()), second=int((1 - g).sum()))
        n_eff = int(min(g.sum(), (1 - g).sum()))
    threshold = max(3, int(min_replicates or 3))
    tier = "insufficient" if (n_eff <= 2 or n_eff < threshold) else ("parametric" if n_eff == 3 else "permutation")
    if tier == "permutation" and not (perm.min_attainable_p < alpha):
        tier = "parametric"
    if tier == "permutation":
        detected = bool(perm.p < alpha)
    elif tier == "parametric":
        detected = bool(np.isfinite(te.p) and te.p < alpha)
    else:
        detected = False
    return dict(te=te, perm=perm, tier=tier, detected=detected, n_eff=n_eff, counts=counts,
                measured=len(reps), unmeasured=len(design.reps) - len(reps))


def donor_test(values: dict, design: Design, *, direction: str, margin: float | None = None,
               margin_source: str | None = None, alpha: float = 0.05, n_perm: int = 1000,
               max_exact: int = 20000, min_replicates: int = 3, rng: np.random.Generator | None = None) -> dict:
    """(b): the graded test on the replicate values, a TOST at `margin` (the SESOI in the metric's units,
    ``margin_source`` "declared"; or, with no margin, the design's floor from this test's own standard
    error, "floor") and the engine's verdict order on the effect field alone: a detected effect against
    the claimed direction is "against"; one in it whose (1 - 2 alpha) interval lies within ±margin, with a
    declared margin, is "below"; any other detected effect is "effect"; an undetected one whose interval
    lies within ±margin is "absent"; the rest, and too few replicates, "none"."""
    if direction not in DIRECTIONS:
        raise ValueError(f"direction must be one of {DIRECTIONS}")
    rng = rng if rng is not None else np.random.default_rng(0)
    inf = _graded(values, design, alpha=alpha, n_perm=n_perm, max_exact=max_exact,
                  min_replicates=min_replicates, rng=rng)
    te, perm = inf["te"], inf["perm"]
    if margin is None:
        margin, source = design_floor(te.se, te.df, alpha), "floor"
    else:
        source = margin_source or "declared"
    lo2, hi2 = te.ci(1 - 2 * alpha)
    margin_ok = margin is not None and np.isfinite(margin) and margin > 0
    equivalent = bool(margin_ok and np.isfinite(lo2) and np.isfinite(hi2) and lo2 > -margin and hi2 < margin)
    if inf["tier"] == "insufficient":
        call = "none"
    elif inf["detected"]:
        call = _direction_call(te.est, direction)
        if call == "effect" and equivalent and source == "declared":
            call = "below"
    else:
        call = "absent" if equivalent else "none"
    return dict(call=call, est=te.est, se=te.se, df=te.df, t_method=te.method, p_t=te.p,
                ci=te.ci(1 - alpha), ci_2alpha=(lo2, hi2), p_perm=perm.p, perm_method=perm.method,
                min_attainable_p=perm.min_attainable_p, tier=inf["tier"], detected=inf["detected"],
                equivalent=equivalent, margin=margin if margin_ok else None, margin_source=source,
                n_eff=inf["n_eff"], counts=inf["counts"], measured=inf["measured"], unmeasured=inf["unmeasured"])


def donor_values(metric: Metric, cells: PerCell, rep_of_cell, in_first, design: Design) -> dict:
    """Each replicate's value of the metric: on its cells (nested), or its difference groups[0] -
    groups[1] (paired), as the engine's `_rep_values` forms them; NaN where the metric gives none."""
    rep_of_cell = np.asarray(rep_of_cell).astype(str)
    in_first = np.asarray(in_first, bool)
    out = {}
    for r in design.reps:
        m = rep_of_cell == str(r)
        if design.kind == "nested":
            out[r] = metric(cells.subset(m))
        else:
            va, vb = metric(cells.subset(m & in_first)), metric(cells.subset(m & ~in_first))
            out[r] = va - vb if np.isfinite(va) and np.isfinite(vb) else float("nan")
    return out
