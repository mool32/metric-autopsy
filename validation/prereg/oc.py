"""Operating characteristics of the success criteria in validation/prereg/v1.md.

One principle for every criterion (decided by the project owner, 2026-10-08): a sound validator
passes the criterion with probability >= 0.90, and a validator whose error is double the
nominal passes with probability <= 0.05. A criterion is a bound on binomial rates: "at most k*
errors of n" (equivalently: the two-sided 95% Clopper-Pearson upper bound of the rate is <= T),
or for decisiveness "at least k* correct definite outcomes". k* is the most lenient threshold at
which the doubled-error validator still passes with probability <= 0.05. S1 is a conjunction over
the key null conditions, so the principle applies to S1 as a whole.

Errors are judged by the cause of a verdict, not its label (decided 2026-10-08): every card has
its allowed outcomes (``panel.allowed``) and its nominal error, the sum of the designed sizes of
the routes to a wrong outcome that are open on it (``nominal_error``). The n of S2-S5 depend on
the truth about the metric on each pair (pilot.json) and on the key's assignment, so ``score.py``
computes every threshold with the rules here from the realized cards; this file shows them on
the design and on scenarios for the truth, and, given pilot.json, on its expected cards.

    python validation/prereg/oc.py > validation/prereg/oc.log
    python validation/prereg/oc.py --pilot pilot.json > oc_pilot.log
"""
from __future__ import annotations

import argparse
import json

import numpy as np

import panel as P

P_PASS_SOUND = 0.90
P_PASS_DOUBLED = 0.05
ALPHA = 0.05
E_SUPPORTED = ALPHA / 2       # false SUPPORTED: the two-sided effect test at alpha, one direction
E_INVALID = 2 * ALPHA         # false "metric invalid": GATE 4 (alpha / 2) + GATE 5's negative control
                              # (alpha) + a silent positive control shown blind (alpha / 2)
E_NDE = ALPHA                 # false NO DETECTABLE EFFECT: the TOST's size at the SESOI
D_NOMINAL = 0.85              # correct definite outcomes on establishable cards


def nominal_error(ok: frozenset) -> float:
    """A card's nominal error: the designed sizes of the routes to an outcome outside its allowed
    set `ok`: SUPPORTED alpha/2, NOT SUPPORTED (metric invalid) 2 alpha, NO DETECTABLE EFFECT alpha."""
    return ((E_SUPPORTED if P.SUPPORTED not in ok else 0.0) + (E_INVALID if P.NS_INVALID not in ok else 0.0)
            + (E_NDE if P.NDE not in ok else 0.0))


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
    if p <= 0:
        return 1.0
    if p >= 1:
        return 0.0
    from scipy import stats
    return float(stats.binom.cdf(k, n, p))


def error_rule(n: int, p0: float):
    """For an error-rate criterion: the largest k* with P(X <= k* | 2 p0) <= 0.05, and
    P(pass | p0), P(pass | 2 p0) for the rule 'at most k* errors'. For a mixture of cards with
    different nominal rates, p0 is their mean; the count is then Poisson-binomial, more
    concentrated than the binomial (Hoeffding 1956), so the binomial rule is conservative on
    both sides."""
    p1 = min(1.0, 2 * p0)
    if n <= 0:
        return -1, 1.0, 1.0
    from scipy import stats
    k = int(stats.binom.ppf(P_PASS_DOUBLED, n, p1))
    while k >= 0 and binom_cdf(k, n, p1) > P_PASS_DOUBLED:
        k -= 1
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
    P(pass | d0), P(pass | d1) for the rule 'at least k* correct definite outcomes'."""
    d1 = 1 - 2 * (1 - d0)
    if n <= 0:
        return 1, 0.0, 0.0
    from scipy import stats
    k = int(stats.binom.isf(P_PASS_DOUBLED, n, d1))
    while k - 1 >= 0 and 1 - binom_cdf(k - 2, n, d1) <= P_PASS_DOUBLED:
        k -= 1
    while 1 - binom_cdf(k - 1, n, d1) > P_PASS_DOUBLED:
        k += 1
    return k, 1 - binom_cdf(k - 1, n, d0), 1 - binom_cdf(k - 1, n, d1)


def joint_pass_probability(rows: list[dict], sims: int = 20_000, seed: int = 0) -> float:
    """P(S1-S5 all pass) for a sound validator, by simulation on the cards `rows` (as
    ``score.card_rows`` gives them: key, sup_error_possible, valid, establishable, condition and
    the allowed routes): on every card a false SUPPORTED with E_SUPPORTED where SUPPORTED is not
    allowed, a false "metric invalid" with E_INVALID where it is not, a false NO DETECTABLE
    EFFECT with E_NDE where it is not (mutually exclusive), and an establishable card without an
    error correct and definite with D_NOMINAL / (1 - its nominal error). The criteria share cards,
    so they are simulated together."""
    if not rows:
        return float("nan")
    rng = np.random.default_rng(seed)
    n = len(rows)
    p_sup = np.array([E_SUPPORTED if r["sup_error_possible"] else 0.0 for r in rows])
    p_inv = np.array([E_INVALID if r["invalid_error_possible"] else 0.0 for r in rows])
    p_nde = np.array([E_NDE if r["nde_error_possible"] else 0.0 for r in rows])
    nominal = p_sup + p_inv + p_nde
    est = np.array([r["establishable"] for r in rows])
    valid = np.array([r["valid"] for r in rows])
    key_groups = {}
    for i, r in enumerate(rows):
        if r["key"]:
            key_groups.setdefault(r["condition"], []).append(i)
    k1 = {c: error_rule(len(ix), E_SUPPORTED)[0] for c, ix in key_groups.items()}
    k2 = error_rule(int((p_sup > 0).sum()), E_SUPPORTED)[0]
    k3 = decisiveness_rule(int(est.sum()), D_NOMINAL)[0]
    k4 = error_rule(n, float(nominal.mean()))[0]
    k5 = error_rule(int(valid.sum()), E_INVALID)[0]
    passed = 0
    chunk = 500
    p_ok = np.clip(D_NOMINAL / np.maximum(1e-12, 1 - nominal), 0, 1)
    for start in range(0, sims, chunk):
        m = min(chunk, sims - start)
        u = rng.random((m, n))
        sup = u < p_sup
        inv = (u >= p_sup) & (u < p_sup + p_inv)
        nde = (u >= p_sup + p_inv) & (u < nominal)
        err = sup | inv | nde
        good = (~err) & est & (rng.random((m, n)) < p_ok)
        s1 = np.ones(m, bool)
        for c, ix in key_groups.items():
            s1 &= sup[:, ix].sum(axis=1) <= k1[c]
        s2 = sup.sum(axis=1) <= k2
        s3 = good.sum(axis=1) >= k3
        s4 = err.sum(axis=1) <= k4
        s5 = (inv & valid).sum(axis=1) <= k5
        passed += int((s1 & s2 & s3 & s4 & s5).sum())
    return passed / sims


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


def expected_rows(pilot: dict, dropped=()) -> list[dict]:
    """The cards of the design with every pool pair equally often (the key draws pairs
    uniformly): per condition and variant, n_datasets cards spread evenly over the pool, with the
    truth, allowed outcomes and establishability of pilot.json. A stand-in for the realized cards
    before the key."""
    rows = []
    est = pilot.get("establishable", {})
    for c in P.CONDITIONS:
        size = int(pilot["pool_size"][c.background])
        for variant, n in c.variants:
            if P._dropped(c.name, variant, dropped):
                continue
            for i in range(n):
                pair = i % size
                ok = P.allowed(c, variant, pair, pilot)
                truth = P.metric_truth(c, pair, pilot)
                for _ in range(c.cards):
                    rows.append(dict(condition=c.name, key=c.key and variant == c.key_variant,
                                     sup_error_possible=P.SUPPORTED not in ok,
                                     invalid_error_possible=P.NS_INVALID not in ok,
                                     nde_error_possible=P.NDE not in ok, valid=truth == "valid",
                                     establishable=c.oracle and bool(est.get(f"{c.name}:{variant}:{pair}", {})
                                                                      .get("establishable")),
                                     nominal=nominal_error(ok)))
    return rows


def show_error(name: str, n: int, p0: float):
    k, ps, pd = error_rule(n, p0)
    hi = clopper_pearson(k, n)[1] if n else float("nan")
    ok = "meets" if ps >= P_PASS_SOUND and pd <= P_PASS_DOUBLED else "FAILS"
    print(f"  {name}: n = {n}, pass with at most {k} errors (upper 95% CP bound <= {hi:.4f}); "
          f"P(pass | {p0:.4f}) = {ps:.3f}, P(pass | {min(1, 2 * p0):.4f}) = {pd:.3f} — {ok} the principle")


def show_joint(name: str, n: int, p0: float, m: int):
    k, pall, pd, pone = joint_error_rule(n, p0, m)
    hi = clopper_pearson(k, n)[1]
    ok = "meets" if pall >= P_PASS_SOUND and pd <= P_PASS_DOUBLED else "FAILS"
    print(f"  {name}: {m} conditions x n = {n}, each passes with at most {k}/{n} false SUPPORTED (upper "
          f"95% CP bound <= {hi:.4f}); P(one passes | {p0:.3f}) = {pone:.4f}, P(all {m} pass | "
          f"{p0:.3f}) = {pall:.4f}, P(a condition passes | {2 * p0:.3f} there) = {pd:.4f} — {ok} the principle")


def show_decisive(name: str, n: int, d0: float):
    k, ps, pd = decisiveness_rule(n, d0)
    lo = clopper_pearson(k, n)[0] if n else float("nan")
    ok = "meets" if ps >= P_PASS_SOUND and pd <= P_PASS_DOUBLED else "FAILS"
    d1 = 1 - 2 * (1 - d0)
    print(f"  {name}: n = {n}, pass with at least {k} correct definite outcomes (lower 95% CP bound "
          f">= {lo:.4f}); P(pass | {d0:.2f}) = {ps:.3f}, P(pass | {d1:.2f}) = {pd:.3f} — {ok} the principle")


def show_design(title: str, rows: list[dict]):
    """The thresholds and operating characteristics of S1-S5 on a set of cards."""
    print(f"## {title}")
    n = len(rows)
    key = {}
    for r in rows:
        if r["key"]:
            key[r["condition"]] = key.get(r["condition"], 0) + 1
    n_sup = sum(r["sup_error_possible"] for r in rows)
    n_valid = sum(r["valid"] for r in rows)
    n_est = sum(r["establishable"] for r in rows)
    e_bar = float(np.mean([r["nominal"] for r in rows]))
    print(f"  cards {n}: SUPPORTED an error on {n_sup}, valid metric on {n_valid}, establishable {n_est}; "
          f"mean nominal error {e_bar:.4f}")
    for c, m in key.items():
        show_error(f"S1 {c}", m, E_SUPPORTED)
    show_error("S2 false SUPPORTED where it is an error", n_sup, E_SUPPORTED)
    show_decisive("S3 correct definite outcomes, establishable cards", n_est, D_NOMINAL)
    show_error("S4 outside the allowed outcomes, all cards (mean nominal)", n, e_bar)
    show_error("S5 false metric invalid, valid-metric cards", n_valid, E_INVALID)
    print(f"  P(S1-S5 all pass | sound), simulated: {joint_pass_probability(rows):.3f}")
    print()


def scenario_pilot(blind_low: bool, blind_medium_share: float, establishable_share: float) -> dict:
    """A stand-in pilot for the operating characteristics before the pilot: 8 pairs per level
    on both backgrounds, every high pair valid, the low pairs blind (or valid), a share of the
    medium pairs blind; Δ* above the SESOI at the key dose; a share of every case establishable."""
    k = P.MAX_PAIRS_PER_LEVEL
    pool = [dict(index=i, level=P.LEVELS[i // k], pair=["a", "b"]) for i in range(3 * k)]
    truth = {}
    for i, pe in enumerate(pool):
        blind = ((pe["level"] == "low" and blind_low)
                 or (pe["level"] == "medium" and (i % k) < round(blind_medium_share * k)))
        truth[str(i)] = {"class": "blind" if blind else "valid"}
    pilot = dict(pool={b: pool for b in P.BACKGROUNDS}, pool_size={b: len(pool) for b in P.BACKGROUNDS},
                 truth={b: truth for b in P.BACKGROUNDS}, sesoi={b: {lv: 0.1 for lv in P.LEVELS} for b in P.BACKGROUNDS},
                 delta={str(i): {f"{f:g}": dict(value=0.125 * f) for f in (0.25, 0.5, 1.0, 1.5)}
                        for i in range(len(pool))}, establishable={})
    for c in P.CONDITIONS:
        for v, _ in c.variants:
            for i in range(len(pool)):
                pilot["establishable"][f"{c.name}:{v}:{i}"] = dict(establishable=(i % 20) < 20 * establishable_share)
    return pilot


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--pilot", help="pilot.json: the operating characteristics on its expected cards")
    args = p.parse_args(argv)
    print("# operating characteristics of the v1 success criteria (two-sided 95% Clopper-Pearson)")
    print(f"# principle: P(pass | sound) >= {P_PASS_SOUND}, P(pass | doubled nominal error) <= {P_PASS_DOUBLED};")
    print("# S1, a conjunction over the key null conditions, meets it as a whole")
    print()
    if args.pilot:
        pilot = json.loads(open(args.pilot).read())
        show_design("The expected cards of pilot.json (every pool pair equally often)",
                    expected_rows(pilot, pilot.get("dropped", ())))
        return
    m = sum(c.key for c in P.CONDITIONS)
    print("## Nominal rates (errors by cause, decided 2026-10-08)")
    print(f"  false SUPPORTED {E_SUPPORTED} (alpha/2: the effect test is two-sided at alpha, the direction checked")
    print(f"  afterwards); false NOT SUPPORTED (metric invalid, GATE 4/5) on a valid metric {E_INVALID} (2 alpha:")
    print("  GATE 4's upper bound below delta_min alpha/2, GATE 5's negative control alpha, a silent positive")
    print(f"  control shown blind alpha/2); false NO DETECTABLE EFFECT {E_NDE} (alpha, the TOST's size). A card's nominal")
    print("  error is the sum over the routes open on it (panel.allowed); S4 uses the cards' mean. Decisiveness:")
    print(f"  {D_NOMINAL} correct definite outcomes on establishable cards (the oracle reaches >= 0.9), doubled shortfall 0.70.")
    print()
    print(f"## S1 as a whole: the number of datasets per key condition ({m} conditions)")
    first = next(n for n in range(400, 2001) if joint_error_rule(n, E_SUPPORTED, m)[1] >= P_PASS_SOUND)
    fails = [n for n in range(first, 1001) if joint_error_rule(n, E_SUPPORTED, m)[1] < P_PASS_SOUND]
    print(f"  the rule first holds at n = {first}; k* moves in steps, so it fails again at n = "
          f"{_ranges(fails)} (up to 1000); the plan's {P.KEY_N} meets it")
    show_joint("S1 false SUPPORTED, key null conditions", P.KEY_N, E_SUPPORTED, m)
    print()
    for title, args_ in (("Scenario A: the low level blind, half the medium pairs blind, half the cases establishable",
                          (True, 0.5, 0.5)),
                         ("Scenario B: every pair valid, every case establishable", (False, 0.0, 1.0)),
                         ("Scenario C: the low level blind, every medium pair blind, a fifth of the cases establishable",
                          (True, 1.0, 0.2))):
        show_design(title, expected_rows(scenario_pilot(*args_)))
    print("## Precision of a rate from n cards (half-width of the 95% CP interval)")
    for n in (100, 200, 400, 790, 1000, 4000):
        row = ", ".join(f"p={q:.3f}: ±{(clopper_pearson(round(q * n), n)[1] - clopper_pearson(round(q * n), n)[0]) / 2:.3f}"
                        for q in (0.025, 0.1, 0.5))
        print(f"  n={n:4d}: {row}")


if __name__ == "__main__":
    main()
