"""Operating characteristics of the success criteria in validation/prereg/v1.md.

One principle for every criterion (decided by the project owner, 2026-10-08): a sound validator
passes with probability >= 0.90, and a validator whose error is double the nominal passes with
probability <= 0.05. A criterion is a bound on a binomial rate from n datasets: "the two-sided
95% Clopper-Pearson upper bound of the error rate is <= T" (equivalently: at most k* errors of
n), or for decisiveness "the lower bound is >= T". For each criterion this script finds the
smallest n for which some threshold meets the principle, and the threshold that the planned n
gives (the most lenient threshold at which the doubled-error validator still passes with
probability <= 0.05).

    python validation/prereg/oc.py > validation/prereg/oc.log
"""
from __future__ import annotations

import math


def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Two-sided exact binomial interval (scipy; the scoring never imports the engine)."""
    from scipy import stats
    lo = 0.0 if k == 0 else float(stats.beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(stats.beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi

P_PASS_SOUND = 0.90
P_PASS_DOUBLED = 0.05


def binom_cdf(k: int, n: int, p: float) -> float:
    """P(X <= k) for X ~ Binomial(n, p)."""
    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    total = 0.0
    for i in range(k + 1):
        total += math.exp(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
                          + (i * math.log(p) if i else 0.0) + (n - i) * math.log1p(-p))
    return min(total, 1.0)


def error_rule(n: int, p0: float):
    """For an error-rate criterion: the largest k* with P(X <= k* | 2 p0) <= 0.05, and
    P(pass | p0), P(pass | 2 p0) for the rule 'at most k* errors'."""
    p1 = 2 * p0
    k = -1
    while binom_cdf(k + 1, n, p1) <= P_PASS_DOUBLED:
        k += 1
    return k, binom_cdf(k, n, p0), binom_cdf(k, n, p1)


def decisiveness_rule(n: int, d0: float):
    """For decisiveness (a share that should be high): nominal non-decisiveness 1 - d0, doubled
    1 - d1 with d1 = 1 - 2 (1 - d0). The smallest k* with P(X >= k* | d1) <= 0.05, and
    P(pass | d0), P(pass | d1) for the rule 'at least k* definite verdicts'."""
    d1 = 1 - 2 * (1 - d0)
    k = n + 1
    while k - 1 >= 0 and 1 - binom_cdf(k - 2, n, d1) <= P_PASS_DOUBLED:
        k -= 1
    p_sound = 1 - binom_cdf(k - 1, n, d0)
    p_doubled = 1 - binom_cdf(k - 1, n, d1)
    return k, p_sound, p_doubled


def smallest_n(fn, start: int, nominal: float, step: int = 10) -> int:
    n = start
    while fn(n, nominal)[1] < P_PASS_SOUND:
        n += step
    lo = max(1, n - step)
    for m in range(lo, n + 1):
        if fn(m, nominal)[1] >= P_PASS_SOUND:
            return m
    return n


def show_error(name: str, n: int, p0: float):
    k, ps, pd = error_rule(n, p0)
    hi = clopper_pearson(k, n)[1]
    hi_next = clopper_pearson(k + 1, n)[1]
    ok = "meets" if ps >= P_PASS_SOUND and pd <= P_PASS_DOUBLED else "FAILS"
    print(f"  {name}: n = {n}, pass with at most {k}/{n} errors (upper 95% CP bound <= {hi:.4f}; any "
          f"T in [{hi:.4f}, {hi_next:.4f}) is the same rule); P(pass | {p0:.3f}) = {ps:.3f}, "
          f"P(pass | {2 * p0:.3f}) = {pd:.3f} — {ok} the principle")


def show_decisive(name: str, n: int, d0: float):
    k, ps, pd = decisiveness_rule(n, d0)
    lo = clopper_pearson(k, n)[0]
    lo_prev = clopper_pearson(k - 1, n)[0]
    ok = "meets" if ps >= P_PASS_SOUND and pd <= P_PASS_DOUBLED else "FAILS"
    d1 = 1 - 2 * (1 - d0)
    print(f"  {name}: n = {n}, pass with at least {k}/{n} definite verdicts (lower 95% CP bound >= "
          f"{lo:.4f}; any T in ({lo_prev:.4f}, {lo:.4f}] is the same rule); P(pass | {d0:.2f}) = {ps:.3f}, "
          f"P(pass | {d1:.2f}) = {pd:.3f} — {ok} the principle")


def main():
    print("# operating characteristics of the v1 success criteria (two-sided 95% Clopper-Pearson)")
    print(f"# principle: P(pass | sound) >= {P_PASS_SOUND}, P(pass | doubled nominal error) <= {P_PASS_DOUBLED}")
    print()
    print("## Nominal rates")
    print("  false SUPPORTED, directional claims: a two-sided test at alpha = 0.05 whose SUPPORTED needs")
    print("  the pre-registered direction gives alpha/2 = 0.025 on a null with a valid metric (less where")
    print("  the metric is invalid); doubled: 0.05. The same nominal 0.025 is used for the share of verdicts")
    print("  outside the allowed set (dominated by the same false SUPPORTED). Decisiveness: nominal 0.85")
    print("  (establishable cases have oracle power >= 0.9; a validator must also demonstrate the metric's")
    print("  response, which the oracle is told); doubled non-decisiveness 0.30, i.e. decisiveness 0.70.")
    print()
    print("## Smallest number of datasets that can meet the principle")
    n_err = smallest_n(error_rule, 200, 0.025)
    n_dec = smallest_n(decisiveness_rule, 20, 0.85, step=1)
    print(f"  error rate, nominal 0.025: n >= {n_err}")
    print(f"  error rate, nominal 0.050 (non-directional): n >= {smallest_n(error_rule, 100, 0.05)}")
    print(f"  decisiveness, nominal 0.85: n >= {n_dec}")
    print()
    print("## The criteria at the planned numbers of datasets")
    show_error("S1 false SUPPORTED, each key null condition", 600, 0.025)
    show_error("S1 at the earlier plan of 400 per condition", 400, 0.025)
    show_error("S2 false SUPPORTED, pooled over N1-N7", 3550, 0.025)
    show_decisive("S3 decisiveness, pooled establishable cases", 200, 0.85)
    show_decisive("S3 at 1000 establishable datasets", 1000, 0.85)
    show_error("S4 outside the allowed set, pooled N and E", 4450, 0.025)
    print()
    print("## Precision of a rate from n datasets (half-width of the 95% CP interval)")
    for n in (100, 200, 400, 600, 1000):
        row = ", ".join(f"p={p:.3f}: ±{(clopper_pearson(round(p * n), n)[1] - clopper_pearson(round(p * n), n)[0]) / 2:.3f}"
                        for p in (0.025, 0.05, 0.5))
        print(f"  n={n:4d}: {row}")


if __name__ == "__main__":
    main()
