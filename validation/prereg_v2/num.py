"""The numbers validation v2's own code needs, without the engine (v2.md section 2: the panel's code
never imports the engine): the regularized incomplete beta, Student's t, and the Clopper-Pearson
interval, also at an effective size (Korn and Graubard 1998).

Student's t takes scipy's functions where scipy is installed and otherwise the same continued fraction
and bisection as the engine's ``stats.py``, so that a test over donors computed here and the engine's own
test give the same p value and interval in either environment (`test_v2_w1.py` checks both against the
engine).
"""
from __future__ import annotations

import math


def _betacf(a: float, b: float, x: float, max_iter: int = 400, eps: float = 3e-15) -> float:
    """The continued fraction of the regularized incomplete beta (Numerical Recipes, betacf)."""
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


def betainc(a: float, b: float, x: float) -> float:
    """The regularized incomplete beta I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    front = math.exp(lbeta + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def beta_ppf(q: float, a: float, b: float) -> float:
    """The q quantile of Beta(a, b), by bisection on the regularized incomplete beta."""
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if betainc(a, b, mid) < q:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-12:
            break
    return 0.5 * (lo + hi)


def t_cdf(t: float, df: float) -> float:
    """Student's t distribution function (scipy's where installed, as the engine's)."""
    if not math.isfinite(t):
        return 1.0 if t > 0 else 0.0
    try:
        from scipy.stats import t as _t  # type: ignore
        return float(_t.cdf(t, df))
    except Exception:
        tail = 0.5 * betainc(df / 2.0, 0.5, df / (df + t * t))
        return 1.0 - tail if t > 0 else tail


def t_ppf(q: float, df: float) -> float:
    """Student's t quantile (scipy's where installed, as the engine's)."""
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


def t_two_sided_p(est: float, se: float, df: float) -> float:
    """The two-sided p of an estimate with its standard error on `df` degrees of freedom."""
    if not (math.isfinite(se) and se > 0 and df > 0):
        return float("nan")
    return float(2.0 * (1.0 - t_cdf(abs(est) / se, df)))


def clopper_pearson(k: float, n: float, alpha: float = 0.05) -> tuple[float, float]:
    """The two-sided 1 - alpha Clopper-Pearson interval of k / n; with a non-integer (effective) k and n
    it is the Korn-Graubard interval. scipy's beta quantiles where installed, else the bisection above."""
    if n <= 0:
        return (float("nan"), float("nan"))
    try:
        from scipy.stats import beta as _beta  # type: ignore

        def ppf(q, a, b):
            return float(_beta.ppf(q, a, b))
    except Exception:
        ppf = beta_ppf
    lo = 0.0 if k <= 0 else ppf(alpha / 2.0, k, n - k + 1)
    hi = 1.0 if k >= n else ppf(1.0 - alpha / 2.0, k + 1, n - k)
    return (float(lo), float(hi))
