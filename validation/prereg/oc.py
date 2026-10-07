"""Operating characteristics of the success criteria in validation/prereg/v1.md.

For a criterion on a binomial rate with n independent datasets, this prints the probability that
the criterion is met as a function of the true rate. Intervals are two-sided 95% Clopper-Pearson,
as in every report of the project.

    python validation/prereg/oc.py > validation/prereg/oc.log
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from metric_autopsy.stats import clopper_pearson  # noqa: E402


def binom_pmf(k: int, n: int, p: float) -> float:
    if p <= 0:
        return 1.0 if k == 0 else 0.0
    if p >= 1:
        return 1.0 if k == n else 0.0
    return math.exp(math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
                    + k * math.log(p) + (n - k) * math.log1p(-p))


def p_upper_below(n: int, p: float, threshold: float) -> float:
    """P(upper 95% CP bound <= threshold) when the true rate is p."""
    return sum(binom_pmf(k, n, p) for k in range(n + 1) if clopper_pearson(k, n)[1] <= threshold)


def p_lower_above(n: int, p: float, threshold: float) -> float:
    """P(lower 95% CP bound >= threshold) when the true rate is p."""
    return sum(binom_pmf(k, n, p) for k in range(n + 1) if clopper_pearson(k, n)[0] >= threshold)


def largest_k_passing_upper(n: int, threshold: float) -> int:
    return max((k for k in range(n + 1) if clopper_pearson(k, n)[1] <= threshold), default=-1)


def smallest_k_passing_lower(n: int, threshold: float) -> int:
    return min((k for k in range(n + 1) if clopper_pearson(k, n)[0] >= threshold), default=n + 1)


def main():
    print("# operating characteristics of the v1 success criteria (two-sided 95% Clopper-Pearson)")
    print()
    print("## False SUPPORTED: criterion 'upper bound <= T'")
    print("A calibrated two-sided test at alpha = 0.05 gives ~5% SUPPORTED on null datasets when the")
    print("metric is valid and judgment is resolved; requiring the effect in the pre-registered")
    print("direction halves that to ~2.5%.")
    for n in (200, 1000, 2000, 2750):
        for T in (0.07, 0.10):
            k = largest_k_passing_upper(n, T)
            row = ", ".join(f"p={p:.3f}: {p_upper_below(n, p, T):.2f}" for p in (0.01, 0.025, 0.05, 0.075, 0.10))
            print(f"  n={n:4d}, T={T:.2f} (passes with at most {k}/{n}): P(pass) — {row}")
    print()
    print("## Criterion 'upper bound <= 2 x nominal' (discriminates nominal from doubled error)")
    for nominal in (0.025, 0.05):
        T = 2 * nominal
        for n in (200, 300, 400, 600, 1000):
            k = largest_k_passing_upper(n, T)
            print(f"  nominal {nominal:.3f}, T={T:.3f}, n={n:4d} (passes with at most {k}/{n}): "
                  f"P(pass | nominal) = {p_upper_below(n, nominal, T):.2f}, "
                  f"P(pass | 1.5 x nominal) = {p_upper_below(n, 1.5 * nominal, T):.2f}, "
                  f"P(pass | 2 x nominal) = {p_upper_below(n, T, T):.2f}")
    print()
    print("## Outside the allowed set, pooled over N and E: criterion 'upper bound <= 0.10'")
    n = 3650
    k = largest_k_passing_upper(n, 0.10)
    row = ", ".join(f"p={p:.3f}: {p_upper_below(n, p, 0.10):.2f}" for p in (0.025, 0.05, 0.075, 0.09, 0.10))
    print(f"  n={n:4d} (passes with at most {k}/{n}): P(pass) — {row}")
    print()
    print("## Decisiveness: criterion 'lower bound >= 0.70'")
    for n in (100, 200, 400):
        k = smallest_k_passing_lower(n, 0.70)
        row = ", ".join(f"d={d:.2f}: {p_lower_above(n, d, 0.70):.2f}" for d in (0.70, 0.75, 0.80, 0.85, 0.90))
        print(f"  n={n:4d} (passes with at least {k}/{n}): P(pass) — {row}")
    print()
    print("## Precision of a rate from n datasets (half-width of the 95% CP interval)")
    for n in (30, 100, 200, 400, 1000):
        row = ", ".join(f"p={p:.2f}: ±{(clopper_pearson(round(p * n), n)[1] - clopper_pearson(round(p * n), n)[0]) / 2:.3f}"
                        for p in (0.05, 0.5))
        print(f"  n={n:4d}: {row}")


if __name__ == "__main__":
    main()
