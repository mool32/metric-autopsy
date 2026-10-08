"""Operating characteristics of the success criteria in validation/prereg/v1.md.

One principle for every criterion (decided by the project owner, 2026-10-08): a sound validator
passes the criterion with probability >= 0.90, and a validator whose error is double the
nominal passes with probability <= 0.05. A criterion is a bound on binomial rates from n
datasets: "at most k* errors of n" (equivalently: the two-sided 95% Clopper-Pearson upper bound
of the rate is <= T), or for decisiveness "at least k* correct definite verdicts". k* is the
most lenient threshold at which the doubled-error validator still passes with probability
<= 0.05. S1 is a conjunction over the key null conditions, so the principle applies to S1 as a
whole: P(every condition passes | sound) >= 0.90, while a doubled error in any one condition
passes that condition with probability <= 0.05.

    python validation/prereg/oc.py > validation/prereg/oc.log
"""
from __future__ import annotations

import math

import numpy as np

import panel as P

P_PASS_SOUND = 0.90
P_PASS_DOUBLED = 0.05
NOMINAL_ERROR = 0.025         # false SUPPORTED with directional claims: alpha / 2
NOMINAL_DECISIVENESS = 0.85   # correct definite verdicts on establishable cards


def clopper_pearson(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Two-sided exact binomial interval (scipy; the scoring never imports the engine)."""
    from scipy import stats
    lo = 0.0 if k == 0 else float(stats.beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(stats.beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


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


def joint_error_rule(n: int, p0: float, m: int):
    """S1 as a whole over m conditions of n datasets each: the per-condition k* (as in
    `error_rule`), P(all m pass | p0 everywhere), P(one condition passes | 2 p0 there)."""
    k, ps, pd = error_rule(n, p0)
    return k, ps ** m, pd, ps


def decisiveness_rule(n: int, d0: float):
    """For decisiveness (a share that should be high): nominal non-decisiveness 1 - d0, doubled
    1 - d1 with d1 = 1 - 2 (1 - d0). The smallest k* with P(X >= k* | d1) <= 0.05, and
    P(pass | d0), P(pass | d1) for the rule 'at least k* correct definite verdicts'."""
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


def design_counts(dropped=()) -> dict:
    """Cards per class in the panel's design (panel.CONDITIONS)."""
    key = [c for c in P.CONDITIONS if c.key]
    n_key = {c.name: dict(c.variants)[c.key_variant] for c in key}
    return dict(key=n_key, null=P.n_cards(dropped, null=True), effect=P.n_cards(dropped, null=False),
                all=P.n_cards(dropped), datasets=P.n_datasets(dropped),
                oracle=sum(n * c.cards for c in P.CONDITIONS for v, n in c.variants
                           if c.oracle and not P._dropped(c.name, v, dropped)))


def joint_pass_probability(counts: dict, n_est: int, sims: int = 200_000, seed: int = 0) -> float:
    """P(S1-S4 all pass) for a sound validator, by simulation: every null card is a false
    SUPPORTED with probability 0.025 (and then outside the allowed set), every real-effect card
    outside with 0.025, and an establishable card correct and definite with 0.85 (a share of
    n_est spread over null and effect cards in proportion). The criteria share cards, so they
    are simulated together."""
    rng = np.random.default_rng(seed)
    key = list(counts["key"].values())
    k1 = [error_rule(n, NOMINAL_ERROR)[0] for n in key]
    n_null, n_eff, n_all = counts["null"], counts["effect"], counts["all"]
    k2 = error_rule(n_null, NOMINAL_ERROR)[0]
    k4 = error_rule(n_all, NOMINAL_ERROR)[0]
    k3 = decisiveness_rule(n_est, NOMINAL_DECISIVENESS)[0]
    key_sup = np.stack([rng.binomial(n, NOMINAL_ERROR, sims) for n in key])
    s1 = np.all(key_sup <= np.asarray(k1)[:, None], axis=0)
    rest_null = rng.binomial(n_null - sum(key), NOMINAL_ERROR, sims)
    null_sup = key_sup.sum(axis=0) + rest_null
    s2 = null_sup <= k2
    eff_out = rng.binomial(n_eff, NOMINAL_ERROR, sims)
    s4 = (null_sup + eff_out) <= k4
    # among establishable cards, outside ones are not correct; conditional on not outside, a
    # card is correct and definite with 0.85 / 0.975
    est_out = rng.binomial(n_est, NOMINAL_ERROR, sims)
    correct = rng.binomial(n_est - est_out, NOMINAL_DECISIVENESS / (1 - NOMINAL_ERROR))
    s3 = correct >= k3
    return float(np.mean(s1 & s2 & s3 & s4))


def _ranges(xs: list) -> str:
    out, start, prev = [], None, None
    for x in xs + [None]:
        if start is None:
            start = prev = x
        elif x is not None and x == prev + 1:
            prev = x
        else:
            out.append(f"{start}" if start == prev else f"{start}-{prev}")
            start = prev = x
    return ", ".join(out)


def show_error(name: str, n: int, p0: float):
    k, ps, pd = error_rule(n, p0)
    hi = clopper_pearson(k, n)[1]
    hi_next = clopper_pearson(k + 1, n)[1]
    ok = "meets" if ps >= P_PASS_SOUND and pd <= P_PASS_DOUBLED else "FAILS"
    print(f"  {name}: n = {n}, pass with at most {k}/{n} errors (upper 95% CP bound <= {hi:.4f}; any "
          f"T in [{hi:.4f}, {hi_next:.4f}) is the same rule); P(pass | {p0:.3f}) = {ps:.3f}, "
          f"P(pass | {2 * p0:.3f}) = {pd:.3f} — {ok} the principle")


def show_joint(name: str, n: int, p0: float, m: int):
    k, pall, pd, pone = joint_error_rule(n, p0, m)
    hi = clopper_pearson(k, n)[1]
    ok = "meets" if pall >= P_PASS_SOUND and pd <= P_PASS_DOUBLED else "FAILS"
    print(f"  {name}: {m} conditions x n = {n}, each passes with at most {k}/{n} false SUPPORTED (upper "
          f"95% CP bound <= {hi:.4f}); P(one passes | {p0:.3f}) = {pone:.4f}, P(all {m} pass | "
          f"{p0:.3f}) = {pall:.4f}, P(a condition passes | {2 * p0:.3f} there) = {pd:.4f} — {ok} the principle")


def show_decisive(name: str, n: int, d0: float):
    k, ps, pd = decisiveness_rule(n, d0)
    lo = clopper_pearson(k, n)[0]
    lo_prev = clopper_pearson(k - 1, n)[0]
    ok = "meets" if ps >= P_PASS_SOUND and pd <= P_PASS_DOUBLED else "FAILS"
    d1 = 1 - 2 * (1 - d0)
    print(f"  {name}: n = {n}, pass with at least {k}/{n} correct definite verdicts (lower 95% CP "
          f"bound >= {lo:.4f}; any T in ({lo_prev:.4f}, {lo:.4f}] is the same rule); P(pass | {d0:.2f}) = "
          f"{ps:.3f}, P(pass | {d1:.2f}) = {pd:.3f} — {ok} the principle")


def main():
    counts = design_counts()
    m = len(counts["key"])
    print("# operating characteristics of the v1 success criteria (two-sided 95% Clopper-Pearson)")
    print(f"# principle: P(pass | sound) >= {P_PASS_SOUND}, P(pass | doubled nominal error) <= {P_PASS_DOUBLED};")
    print("# S1, a conjunction over the key null conditions, meets it as a whole")
    print()
    print("## Nominal rates")
    print("  false SUPPORTED, directional claims: the engine tests at alpha = 0.05 two-sided and checks")
    print("  the direction afterwards (test_a_directional_claim_is_tested_two_sided_at_alpha), so a null")
    print("  with a valid metric gives alpha/2 = 0.025 (less where the metric is invalid); doubled: 0.05.")
    print("  The same nominal 0.025 is used for the share of verdicts outside the allowed set. Decisiveness")
    print("  counts correct definite verdicts (definite and in the allowed set): nominal 0.85 (establishable")
    print("  cases have oracle power >= 0.9; a validator must also demonstrate the metric's response, which")
    print("  the oracle is told); doubled shortfall 0.30, i.e. 0.70.")
    print()
    print("## The design (panel.py)")
    print(f"  key null conditions (S1): {', '.join(f'{k} {v}' for k, v in counts['key'].items())}")
    print(f"  null cards {counts['null']}, real-effect cards {counts['effect']}, all cards {counts['all']} "
          f"({counts['datasets']} datasets); cards the oracle can establish (not N4): {counts['oracle']}")
    print()
    print(f"## S1 as a whole: the number of datasets per key condition ({m} conditions)")
    first = next(n for n in range(400, 2001) if joint_error_rule(n, NOMINAL_ERROR, m)[1] >= P_PASS_SOUND)
    fails = [n for n in range(first, 1001) if joint_error_rule(n, NOMINAL_ERROR, m)[1] < P_PASS_SOUND]
    print(f"  the rule first holds at n = {first}; k* moves in steps, so it fails again at n = "
          f"{_ranges(fails)} (up to 1000); the plan's {P.KEY_N} meets it")
    for n in (600, 700, 763, 770, 785, 790, 800, 850):
        k, pall, pd, pone = joint_error_rule(n, NOMINAL_ERROR, m)
        mark = "meets" if pall >= P_PASS_SOUND and pd <= P_PASS_DOUBLED else "fails"
        print(f"    n = {n}: k* = {k}, P(one | sound) = {pone:.4f}, P(all {m} | sound) = {pall:.4f}, "
              f"P(one | doubled) = {pd:.4f} — {mark}")
    print()
    print("## The criteria at the planned numbers of datasets")
    show_joint("S1 false SUPPORTED, key null conditions", P.KEY_N, NOMINAL_ERROR, m)
    show_error("S2 false SUPPORTED, all null cards", counts["null"], NOMINAL_ERROR)
    for n in (200, 1000, 3000, counts["oracle"]):
        show_decisive("S3 correct definite verdicts, establishable cards", n, NOMINAL_DECISIVENESS)
    show_error("S4 outside the allowed set, all cards", counts["all"], NOMINAL_ERROR)
    print()
    print("## S1-S4 together, for a sound validator (simulated; the criteria share cards)")
    for n_est in (200, 1000, 3000, counts["oracle"]):
        print(f"  n_est = {n_est}: P(S1, S2, S3 and S4 all pass) = {joint_pass_probability(counts, n_est):.3f}")
    print()
    print("## Precision of a rate from n datasets (half-width of the 95% CP interval)")
    for n in (100, 200, 400, 790, 1000):
        row = ", ".join(f"p={p:.3f}: ±{(clopper_pearson(round(p * n), n)[1] - clopper_pearson(round(p * n), n)[0]) / 2:.3f}"
                        for p in (0.025, 0.05, 0.5))
        print(f"  n={n:4d}: {row}")


if __name__ == "__main__":
    main()
