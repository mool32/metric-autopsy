"""Power of the real claims' decision rule (validation/plans/real_claims.md), before any selection.

Exact binomial arithmetic, no data.

1. The second question, rightness: a one-sided sign test at 0.05 on the deciding votes. For each
   number of votes n, the critical value (the fewest votes for the engine that reject a coin at
   0.05), the test's exact size, and its power when the engine is right in 80%, 85% or 90% of the
   disagreements; then the smallest n at which an engine right in 85% passes with probability at
   least 0.8, and the numbers above it where it does not.
2. The first question, frequency: the share of papers where the engine changes a donor-level
   test's conclusion, against a threshold set in advance. For the proposed threshold, the
   probability that the rule passes at a true frequency of 0.05 to 0.40, for 10 to 40 papers.

    python validation/exploratory/real_claims_power.py > validation/exploratory/real_claims_power.log
"""
from __future__ import annotations

from math import comb

ALPHA = 0.05
THRESHOLD = 0.20  # the proposed frequency threshold: one paper in five


def upper(n: int, k: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p)."""
    return sum(comb(n, j) * p ** j * (1 - p) ** (n - j) for j in range(max(k, 0), n + 1))


def critical(n: int) -> int | None:
    """The fewest votes for the engine, out of n, that reject a fair coin one-sided at ALPHA."""
    return next((k for k in range(n + 1) if upper(n, k, 0.5) <= ALPHA), None)


def main() -> int:
    out = ["# Power of the real claims' decision rule (exact binomial; no data)", "",
           "## 1. Rightness: a one-sided sign test at 0.05 on the deciding votes",
           "  n votes: critical value c (votes for the engine needed), exact size P(X >= c | 0.5);",
           "  power P(X >= c | p) for an engine right in a share p of the disagreements",
           "     n   c   size    p=0.80  p=0.85  p=0.90"]
    power85 = {}
    for n in range(5, 41):
        c = critical(n)
        pw = [upper(n, c, p) for p in (0.80, 0.85, 0.90)]
        power85[n] = pw[1]
        if n <= 30 or n % 5 == 0:
            out.append(f"  {n:>4} {c:>3}  {upper(n, c, 0.5):.4f}   {pw[0]:.3f}   {pw[1]:.3f}   {pw[2]:.3f}")
    first = next(n for n in sorted(power85) if power85[n] >= 0.8)
    fails = [n for n in sorted(power85) if n > first and power85[n] < 0.8]
    low = min((power85[n], n) for n in power85 if n >= first)
    out += ["", f"  the smallest n at which an engine right in 85% passes with probability >= 0.8: {first} "
            f"(power {power85[first]:.3f})",
            f"  numbers of votes above it, up to 40, where it does not: {fails if fails else 'none'}; the lowest "
            f"power from {first} up: {low[0]:.3f} at n = {low[1]}",
            f"  at 15 votes: power {power85[15]:.3f}", ""]

    out += ["## 2. Frequency: the share of papers where the engine changes the donor-level test's conclusion",
            f"  rule: the share observed is at least {THRESHOLD:.2f} (at least ceil({THRESHOLD:.2f} N) of N papers);",
            "  the probability that the rule passes at a true frequency f",
            "     N  needed   f=0.05  f=0.10  f=0.20  f=0.30  f=0.40"]
    for n in (10, 15, 20, 25, 30, 40):
        k = -(-int(round(THRESHOLD * n * 1000)) // 1000)  # ceil without float error
        out.append(f"  {n:>4} {k:>6}    " + "  ".join(f"{upper(n, k, f):.3f}" for f in (0.05, 0.10, 0.20, 0.30, 0.40)))
    out.append("")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
