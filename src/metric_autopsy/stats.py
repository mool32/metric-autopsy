"""Small statistics used by the effect, QC and control gates.

Works with numpy alone; scipy is used for t quantiles when installed and a pure-python
fallback takes over otherwise, so the core install stays numpy + pandas.

Conventions: two-sided tests; intervals are equal-tailed; ``alpha`` is the family level the
caller has already adjusted (e.g. Bonferroni across strata).
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from math import comb

import numpy as np


# --------------------------------------------------------------------------- #
# distributions
# --------------------------------------------------------------------------- #
def _betacf(a: float, b: float, x: float, max_iter: int = 400, eps: float = 3e-15) -> float:
    """Continued fraction for the regularized incomplete beta (Numerical Recipes)."""
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def _betainc(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    front = math.exp(lbeta + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def t_cdf(t: float, df: float) -> float:
    """Student t CDF."""
    if not np.isfinite(t):
        return 1.0 if t > 0 else 0.0
    try:
        from scipy.stats import t as _t  # type: ignore
        return float(_t.cdf(t, df))
    except Exception:
        tail = 0.5 * _betainc(df / 2.0, 0.5, df / (df + t * t))
        return 1.0 - tail if t > 0 else tail


def t_ppf(q: float, df: float) -> float:
    """Student t quantile (inverse CDF)."""
    try:
        from scipy.stats import t as _t  # type: ignore
        return float(_t.ppf(q, df))
    except Exception:
        lo, hi = -1e3, 1e3
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if t_cdf(mid, df) < q:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)


def norm_sf(z: float) -> float:
    """Standard normal survival function."""
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def t_two_sided_p(est: float, se: float, df: float) -> float:
    if not (np.isfinite(se) and se > 0 and df > 0):
        return float("nan")
    return float(2.0 * (1.0 - t_cdf(abs(est) / se, df)))


# --------------------------------------------------------------------------- #
# estimates and intervals for replicate-level values
# --------------------------------------------------------------------------- #
@dataclass
class TEstimate:
    """A difference (or mean difference) with its t-based uncertainty."""

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


def stratum_weights(g: np.ndarray, s: np.ndarray) -> dict:
    """Precision weights 1 / (1/n_A + 1/n_B) for strata that contain both groups."""
    out = {}
    for st in np.unique(s):
        n_a = int(np.sum((s == st) & (g == 1)))
        n_b = int(np.sum((s == st) & (g == 0)))
        if n_a and n_b:
            out[st] = 1.0 / (1.0 / n_a + 1.0 / n_b)
    return out


def stratified_diff(y: np.ndarray, g: np.ndarray, s: np.ndarray, weights: dict | None = None) -> float:
    """Precision-weighted mean over strata of (mean of group 1 - mean of group 0)."""
    weights = weights if weights is not None else stratum_weights(g, s)
    if not weights:
        return float("nan")
    num = 0.0
    for st, w in weights.items():
        m = s == st
        num += w * (y[m & (g == 1)].mean() - y[m & (g == 0)].mean())
    return float(num / sum(weights.values()))


def nested_t(y: np.ndarray, g: np.ndarray, s: np.ndarray) -> TEstimate:
    """t estimate of the group difference for replicates nested in groups (and strata).

    One stratum: Welch. Several strata: precision-weighted stratified difference with the
    pooled within-(stratum, group) variance (assumes a common replicate-level variance).
    """
    weights = stratum_weights(g, s)
    est = stratified_diff(y, g, s, weights)
    if len(weights) == 1:
        st = next(iter(weights))
        a, b = y[(s == st) & (g == 1)], y[(s == st) & (g == 0)]
        if len(a) < 2 or len(b) < 2:
            return TEstimate(est, float("nan"), 0.0, "welch")
        va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
        se = math.sqrt(va + vb)
        if se == 0:
            return TEstimate(est, 0.0, float(len(a) + len(b) - 2), "welch")
        df = (va + vb) ** 2 / ((va ** 2) / (len(a) - 1) + (vb ** 2) / (len(b) - 1)) if (va or vb) else 1.0
        return TEstimate(est, se, float(df), "welch")
    ss, df = 0.0, 0
    for st in weights:
        for grp in (0, 1):
            v = y[(s == st) & (g == grp)]
            if len(v) >= 2:
                ss += float(((v - v.mean()) ** 2).sum())
                df += len(v) - 1
    if df == 0:
        return TEstimate(est, float("nan"), 0.0, "pooled-stratified")
    se = math.sqrt(ss / df / sum(weights.values()))
    return TEstimate(est, se, float(df), "pooled-stratified")


def paired_t(d: np.ndarray) -> TEstimate:
    """One-sample t on within-replicate differences."""
    d = np.asarray(d, float)
    if len(d) < 2:
        return TEstimate(float(d.mean()) if len(d) else float("nan"), float("nan"), 0.0, "paired")
    return TEstimate(float(d.mean()), float(d.std(ddof=1) / math.sqrt(len(d))), float(len(d) - 1), "paired")


def mde(se: float, df: float, alpha: float, power: float) -> float:
    """Minimum detectable effect of a two-sided t test at the given power."""
    if not (np.isfinite(se) and df > 0):
        return float("nan")
    return float((t_ppf(1 - alpha / 2, df) + t_ppf(power, df)) * se)


# --------------------------------------------------------------------------- #
# permutation tests over replicates
# --------------------------------------------------------------------------- #
@dataclass
class PermResult:
    """Two-sided permutation test; ``mc_se`` is 0 for exact enumeration."""

    p: float
    method: str            # "exact" | "monte-carlo"
    n: int                 # assignments enumerated, or permutations drawn
    mc_se: float
    min_attainable_p: float
    observed: float


def _sum_distribution(per_stratum: list[np.ndarray]) -> np.ndarray:
    """All sums choosing one value from each array (outer-sum, flattened)."""
    acc = np.zeros(1)
    for arr in per_stratum:
        acc = (acc[:, None] + arr[None, :]).ravel()
    return acc


def permutation_test_nested(y, g, s, *, n_perm: int = 1000, rng=None, max_exact: int = 20000) -> PermResult:
    """Permute group labels of replicates within strata; statistic = stratified difference."""
    y, g, s = np.asarray(y, float), np.asarray(g, int), np.asarray(s)
    rng = rng if rng is not None else np.random.default_rng(0)
    weights = stratum_weights(g, s)
    obs = stratified_diff(y, g, s, weights)
    if not weights:
        return PermResult(float("nan"), "none", 0, 0.0, float("nan"), obs)
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
        stats = _sum_distribution(per)
        p = float(np.mean(np.abs(stats) >= abs(obs) - eps))
        return PermResult(p, "exact", int(total), 0.0, (2.0 if balanced else 1.0) / total, obs)
    count = 0
    g_perm = g.copy()
    for _ in range(n_perm):
        for st in strata:
            idx = np.where(s == st)[0]
            g_perm[idx] = rng.permutation(g[idx])
        if abs(stratified_diff(y, g_perm, s, weights)) >= abs(obs) - eps:
            count += 1
    p = (1.0 + count) / (1.0 + n_perm)
    return PermResult(p, "monte-carlo", int(n_perm), math.sqrt(p * (1 - p) / n_perm), 1.0 / (1 + n_perm), obs)


def permutation_test_paired(d, *, n_perm: int = 1000, rng=None, max_exact: int = 20000) -> PermResult:
    """Sign-flip test on within-replicate differences; statistic = mean difference."""
    d = np.asarray(d, float)
    rng = rng if rng is not None else np.random.default_rng(0)
    obs = float(d.mean()) if len(d) else float("nan")
    n = len(d)
    if n == 0:
        return PermResult(float("nan"), "none", 0, 0.0, float("nan"), obs)
    eps = 1e-12 * max(1.0, abs(obs))
    total = 2 ** n
    if total <= max_exact:
        signs = np.array(list(itertools.product((1.0, -1.0), repeat=n)))
        stats = (signs * d[None, :]).mean(axis=1)
        p = float(np.mean(np.abs(stats) >= abs(obs) - eps))
        return PermResult(p, "exact", int(total), 0.0, 2.0 / total, obs)
    flips = rng.choice((1.0, -1.0), size=(n_perm, n))
    stats = (flips * d[None, :]).mean(axis=1)
    count = int(np.sum(np.abs(stats) >= abs(obs) - eps))
    p = (1.0 + count) / (1.0 + n_perm)
    return PermResult(p, "monte-carlo", int(n_perm), math.sqrt(p * (1 - p) / n_perm), 1.0 / (1 + n_perm), obs)


# --------------------------------------------------------------------------- #
# empirical nulls
# --------------------------------------------------------------------------- #
def empirical_two_sided_p(x: float, null: np.ndarray, alpha: float) -> tuple[float, str]:
    """Two-sided p of ``x`` against an empirical null, centred on the null median.

    Rank-based ``(1 + #more extreme) / (1 + n)`` when ``alpha`` is attainable with ``n``
    draws; otherwise a normal approximation (mean, sd of the null), reported as such.
    """
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    n = len(null)
    if n < 5 or not np.isfinite(x):
        return float("nan"), "none"
    if 1.0 / (n + 1) < alpha:
        c = float(np.median(null))
        more = int(np.sum(np.abs(null - c) >= abs(x - c)))
        return (1.0 + more) / (1.0 + n), "empirical"
    sd = float(null.std(ddof=1))
    if sd == 0:
        return (0.0 if x != float(null.mean()) else 1.0), "normal"
    return float(2.0 * norm_sf(abs(x - float(null.mean())) / sd)), "normal"
