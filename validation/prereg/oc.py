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
P_JOINT_SOUND = 0.90          # P(S1-S5 together | sound) (decided 2026-10-08, third round)
ALPHA = 0.05
E_SUPPORTED = ALPHA / 2       # false SUPPORTED: the two-sided effect test at alpha, one direction
E_INVALID = 2 * ALPHA         # false "metric invalid": GATE 4 (alpha / 2) + GATE 5's negative control
                              # (alpha) + a silent positive control shown blind (alpha / 2)
E_GATE5 = 1.5 * ALPHA         # the GATE 5 share of it (negative control alpha, silent positive control alpha/2)
E_NDE = ALPHA                 # false NO DETECTABLE EFFECT (or "explained by depth"): the TOST's size
E_OPPOSITE = ALPHA / 2        # a real effect detected against its direction: one tail of the test
D_NOMINAL = 0.85              # correct definite outcomes on establishable cards: where no pilot measures it,
                              # and the most S3's thresholds ask for (a measured rate near 1 leaves no doubled
                              # shortfall to tell apart, which would leave the stratum unjudged)
STRATA = ("effect", "null", "invalid")  # S3's strata (decided after the first review)


def nominal_error(ok: frozenset) -> float:
    """The designed sizes of the routes to an outcome outside the allowed set `ok`, where the pilot
    measures nothing (N4, scenarios): SUPPORTED alpha/2, NOT SUPPORTED (metric invalid) 2 alpha, NO
    DETECTABLE EFFECT alpha, NOT SUPPORTED (opposite direction) alpha/2."""
    return ((E_SUPPORTED if P.SUPPORTED not in ok else 0.0) + (E_INVALID if P.NS_INVALID not in ok else 0.0)
            + (E_NDE if P.NDE not in ok else 0.0) + (E_OPPOSITE if P.NS_OPPOSITE not in ok else 0.0))


def stratum(cond: P.Condition, truth: str) -> str | None:
    """S3's stratum of a card: real effects with a valid metric (the correct definite outcome is
    SUPPORTED, or NO DETECTABLE EFFECT below the SESOI); null data with a valid metric (NO
    DETECTABLE EFFECT or NOT SUPPORTED); a blind or useless metric (metric invalid, or DEGENERATE
    for the constant). Ambiguous metrics belong to none: either kind of verdict is correct there."""
    if truth in ("blind", "useless", "constant"):
        return "invalid"
    if truth == "valid":
        return "effect" if cond.data == "effect" else "null"
    return None


def _shares(counts: dict) -> tuple[dict, int]:
    n = sum(counts.values())
    return ({k: v / n for k, v in counts.items()} if n else {}), n


def sound_model(cond: P.Condition, variant: str, pair: int, pilot: dict) -> dict:
    """What a sound validator — one that follows the engine's pre-registered rules correctly —
    gets on a card, from the pilot's measurements of its case (decided after the first review:
    the nominals follow every route the rules leave open, at its measured size):

    * a valid metric: GATE 4 FAILs with the case's measured odds (pilot.json "gate4", from the
      metric's response on the case's datasets) and GATE 5 with its designed sizes (1.5 alpha):
      a false "metric invalid"; otherwise, after a GATE 4 PASS, the effect's analysis decides as
      the oracle's did on the case's datasets (with the metric's validity taken as known);
    * a blind, ambiguous or useless metric: the engine's rules on the case's datasets, where the
      oracle ran GATE 4's rule — "metric invalid" where it FAILed, the effect's outcome where it
      PASSed, INCONCLUSIVE where it was UNTESTED;
    * the constant: DEGENERATE METRIC.

    Returns the probabilities of a false SUPPORTED (`p_sup`), a false "metric invalid" (`p_inv`),
    any outcome outside the allowed set (`p_err`, both included) and a correct definite outcome
    (`decisive`); `measured` says whether the pilot measured them (else the designed sizes of
    `nominal_error`, and D_NOMINAL where the errors leave room for it)."""
    ok = P.allowed(cond, variant, pair, pilot)
    good = P.definite(cond, variant, pair, pilot)
    truth = P.metric_truth(cond, variant, pair, pilot)
    est = (pilot.get("establishable") or {}).get(f"{cond.name}:{variant}:{pair}") or {}
    if truth == "constant":
        return dict(p_sup=0.0, p_inv=0.0, p_err=0.0, decisive=1.0, measured=True)
    if truth == "valid" and est.get("effect_outcomes"):
        odds = P.truth_record(cond, variant, pair, pilot).get("gate4")
        if odds is not None:
            f, _ = _shares(est["effect_outcomes"])
            reach = float(odds["p_pass"]) * (1 - E_GATE5)  # the effect is analysed only after both gates
            p_inv = float(odds["p_fail"]) + E_GATE5
            p_sup = reach * (f.get(P.SUPPORTED, 0.0) if P.SUPPORTED not in ok else 0.0)
            p_eff = reach * sum(v for o, v in f.items() if o not in ok)
            return dict(p_sup=p_sup, p_inv=p_inv, p_err=p_inv + p_eff,
                        decisive=reach * sum(v for o, v in f.items() if o in good), measured=True)
    if truth in ("blind", "ambiguous", "useless") and est.get("engine_outcomes"):
        g, _ = _shares(est["engine_outcomes"])
        return dict(p_sup=g.get(P.SUPPORTED, 0.0) if P.SUPPORTED not in ok else 0.0, p_inv=0.0,
                    p_err=sum(v for o, v in g.items() if o not in ok),
                    decisive=sum(v for o, v in g.items() if o in good), measured=True)
    err = nominal_error(ok)
    return dict(p_sup=E_SUPPORTED if P.SUPPORTED not in ok else 0.0,
                p_inv=E_INVALID if (truth == "valid" and P.NS_INVALID not in ok) else 0.0,
                p_err=err, decisive=min(D_NOMINAL, 1.0 - err), measured=False)


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
    P(pass | p0), P(pass | 2 p0) for the rule 'at most k* errors' (with p0 = 0, where no doubled
    rate can be told apart, k* = 0: any such error fails). For a mixture of cards with
    different nominal rates, p0 is their mean; the count is then Poisson-binomial, more
    concentrated than the binomial (Hoeffding 1956), so the binomial rule is conservative on
    both sides."""
    p1 = min(1.0, 2 * p0)
    if n <= 0:
        return -1, 1.0, 1.0
    if p0 <= 0:  # the sound validator cannot make this error on these cards: none is allowed
        return 0, 1.0, 1.0
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


def decisiveness_applies(n: int, d0: float) -> bool:
    """Whether a decisiveness criterion on n cards can meet the principle at the nominal d0: some
    threshold passes d0 with >= 0.90 and the doubled shortfall with <= 0.05. Below that a stratum
    is reported, not judged (decided after the first review)."""
    if n <= 0 or not 0 < d0 < 1:
        return False
    _, ps, pd = decisiveness_rule(n, d0)
    return ps >= P_PASS_SOUND and pd <= P_PASS_DOUBLED


def criteria_rules(rows: list[dict], sims: int = 20_000, seed: int = 0) -> dict:
    """The thresholds of S1-S5 on a set of cards (as ``score.card_rows`` and ``expected_rows`` give
    them): S1 per key condition at alpha/2; S2 at the mean measured rate of a false SUPPORTED on its
    cards; S3 per stratum at the mean measured decisiveness of its establishable cards, at most
    D_NOMINAL (0.85, doubled shortfall 0.70); S4 at the cards' mean measured error; S5 at the mean
    measured rate of a false "metric invalid" on the valid-metric cards. S3's strata are judged in the order of STRATA where the principle is
    attainable on the stratum (`decisiveness_applies`) and judging it keeps P(S1-S5 together |
    sound) >= 0.90 (decided in the third round; simulated, `_simulate_passes`); the others are
    reported only. `joint` is P(S1-S5 together | sound) with the strata judged."""
    out = dict(S1={}, S3={})
    for c in P.CONDITIONS:
        if c.key:
            n = sum(1 for r in rows if r["key"] and r["condition"] == c.name)
            out["S1"][c.name] = dict(n=n, max_allowed=error_rule(n, E_SUPPORTED)[0] if n else -1, nominal=E_SUPPORTED)
    s2 = [r for r in rows if r["sup_error_possible"]]
    p2 = float(np.mean([r["p_sup"] for r in s2])) if s2 else 0.0
    out["S2"] = dict(n=len(s2), nominal=p2, max_allowed=error_rule(len(s2), p2)[0] if s2 else -1)
    for st in STRATA:
        rs = [r for r in rows if r["establishable"] and r["stratum"] == st]
        measured = float(np.mean([r["decisive"] for r in rs])) if rs else float("nan")
        d0 = min(measured, D_NOMINAL) if rs else float("nan")
        applies = decisiveness_applies(len(rs), d0)
        out["S3"][st] = dict(n=len(rs), nominal=d0, measured=measured, applies=applies, judged=False,
                             min_required=decisiveness_rule(len(rs), d0)[0] if applies else None)
    e4 = float(np.mean([r["nominal"] for r in rows])) if rows else 0.0
    out["S4"] = dict(n=len(rows), nominal=e4, max_allowed=error_rule(len(rows), e4)[0] if rows else -1)
    s5 = [r for r in rows if r["valid"]]
    p5 = float(np.mean([r["p_inv"] for r in s5])) if s5 else 0.0
    out["S5"] = dict(n=len(s5), nominal=p5, max_allowed=error_rule(len(s5), p5)[0] if s5 else -1)
    if not rows:
        out["joint"] = float("nan")
        return out
    passes = _simulate_passes(rows, out, sims, seed)
    base = passes["S1"] & passes["S2"] & passes["S4"] & passes["S5"]
    for st in STRATA:
        if out["S3"][st]["applies"] and float(np.mean(base & passes[f"S3:{st}"])) >= P_JOINT_SOUND:
            out["S3"][st]["judged"] = True
            base = base & passes[f"S3:{st}"]
    out["joint"] = float(np.mean(base))
    return out


def _simulate_passes(rows: list[dict], rules: dict, sims: int, seed: int) -> dict:
    """For a sound validator, by simulation on the cards `rows` with each card's model
    (``sound_model``): on every card a false SUPPORTED with p_sup, a false "metric invalid" with
    p_inv, another outcome outside the allowed set with p_err - p_sup - p_inv, and on an
    establishable card otherwise a correct definite outcome with `decisive`. The criteria share
    cards, so they are simulated on the same draws: returns, per criterion (S1 as a whole, S2, each
    attainable S3 stratum, S4, S5), whether it passes in each simulation."""
    rng = np.random.default_rng(seed)
    n = len(rows)
    p_sup = np.array([r["p_sup"] for r in rows])
    p_inv = np.array([r["p_inv"] for r in rows])
    p_err = np.maximum(np.array([r["p_err"] for r in rows]), p_sup + p_inv)
    p_good = np.minimum(np.array([r["decisive"] if r["establishable"] else 0.0 for r in rows]), 1 - p_err)
    valid = np.array([r["valid"] for r in rows])
    sup_possible = np.array([r["sup_error_possible"] for r in rows])
    key_groups = {}
    for i, r in enumerate(rows):
        if r["key"]:
            key_groups.setdefault(r["condition"], []).append(i)
    strata = {st: np.array([r["establishable"] and r["stratum"] == st for r in rows]) for st in STRATA
              if rules["S3"][st]["applies"]}
    out = {name: [] for name in ("S1", "S2", "S4", "S5", *(f"S3:{st}" for st in strata))}
    chunk = 500
    for start in range(0, sims, chunk):
        m = min(chunk, sims - start)
        u = rng.random((m, n))
        sup = u < p_sup
        inv = (u >= p_sup) & (u < p_sup + p_inv)
        err = u < p_err
        good = (u >= p_err) & (u < p_err + p_good)
        s1 = np.ones(m, bool)
        for c, ix in key_groups.items():
            s1 &= sup[:, ix].sum(axis=1) <= rules["S1"][c]["max_allowed"]
        out["S1"].append(s1)
        out["S2"].append((sup & sup_possible).sum(axis=1) <= rules["S2"]["max_allowed"] if rules["S2"]["n"]
                         else np.ones(m, bool))
        for st, member in strata.items():
            out[f"S3:{st}"].append((good & member).sum(axis=1) >= rules["S3"][st]["min_required"])
        out["S4"].append(err.sum(axis=1) <= rules["S4"]["max_allowed"])
        out["S5"].append((inv & valid).sum(axis=1) <= rules["S5"]["max_allowed"] if rules["S5"]["n"]
                         else np.ones(m, bool))
    return {k: np.concatenate(v) for k, v in out.items()}


def joint_pass_probability(rows: list[dict], sims: int = 20_000, seed: int = 0) -> float:
    """P(S1-S5 all pass) for a sound validator on the cards `rows`, with S3's judged strata
    (``criteria_rules``)."""
    return criteria_rules(rows, sims, seed)["joint"]


def error_rule_deff(n: int, p0: float, deff: float):
    """`error_rule`'s threshold k* (of the independent design) when the count of errors is
    overdispersed by a design effect `deff` (datasets sharing donors): a beta-binomial with the
    binomial's mean and deff times its variance (intra-class correlation (deff - 1) / (n - 1)).
    Returns k*, P(pass | p0) and P(pass | 2 p0)."""
    k, ps, pd = error_rule(n, p0)
    if deff <= 1.0 or n <= 1:
        return k, ps, pd
    from scipy import stats
    rho = min((deff - 1.0) / (n - 1.0), 1 - 1e-12)

    def cdf(p):
        a, b = p * (1 / rho - 1), (1 - p) * (1 / rho - 1)
        return float(stats.betabinom.cdf(k, n, a, b))
    return k, cdf(p0), cdf(min(1.0, 2 * p0))


def joint_error_rule_deff(n: int, p0: float, m: int, deff: float):
    """S1 as a whole when each condition's count of false SUPPORTED is overdispersed by `deff`
    (`error_rule_deff`), at the thresholds of the independent design. Returns P(all m pass | p0)
    and P(a condition passes | 2 p0 there)."""
    _, ps, pd = error_rule_deff(n, p0, deff)
    return ps ** m, pd


def deff_limit(holds, lo: float = 1.0, hi: float = 2.0) -> float:
    """The largest design effect in [lo, hi] at which `holds(deff)` is still true (bisection;
    `holds` true at lo and turning false once)."""
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if holds(mid) else (lo, mid)
    return lo


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
                truth = P.metric_truth(c, variant, pair, pilot)
                model = sound_model(c, variant, pair, pilot)
                for _ in range(c.cards):
                    rows.append(dict(condition=c.name, key=c.key and variant == c.key_variant,
                                     sup_error_possible=P.SUPPORTED not in ok, valid=truth == "valid",
                                     stratum=stratum(c, truth),
                                     establishable=c.oracle and bool(est.get(f"{c.name}:{variant}:{pair}", {})
                                                                      .get("establishable")),
                                     p_sup=model["p_sup"], p_inv=model["p_inv"], p_err=model["p_err"],
                                     decisive=model["decisive"], nominal=model["p_err"],
                                     measured=model["measured"]))
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
    rules = criteria_rules(rows)
    n_est = sum(r["establishable"] for r in rows)
    measured = sum(r.get("measured", False) for r in rows)
    print(f"  cards {len(rows)}: SUPPORTED an error on {rules['S2']['n']}, valid metric on {rules['S5']['n']}, "
          f"establishable {n_est}; mean nominal error {rules['S4']['nominal']:.4f}"
          + (f"; nominals measured by the pilot on {measured}" if measured else ""))
    for c, rule in rules["S1"].items():
        show_error(f"S1 {c}", rule["n"], E_SUPPORTED)
    show_error("S2 false SUPPORTED where it is an error (mean nominal)", rules["S2"]["n"], rules["S2"]["nominal"])
    for st in STRATA:
        r3 = rules["S3"][st]
        if r3["judged"]:
            show_decisive(f"S3 {st}: correct definite outcomes, establishable cards (nominal: the measured "
                          f"{r3['measured']:.3f}, at most {D_NOMINAL})", r3["n"], r3["nominal"])
        elif r3["applies"]:
            show_decisive(f"S3 {st} (reported only: judging it would bring P(S1-S5 together | sound) below "
                          f"{P_JOINT_SOUND}; measured {r3['measured']:.3f})", r3["n"], r3["nominal"])
        else:
            print(f"  S3 {st}: n = {r3['n']} establishable cards, nominal {r3['nominal']:.3f} (measured "
                  f"{r3['measured']:.3f}): the principle is not attainable, reported only")
    show_error("S4 outside the allowed outcomes, all cards (mean nominal)", rules["S4"]["n"], rules["S4"]["nominal"])
    show_error("S5 false metric invalid, valid-metric cards (mean nominal)", rules["S5"]["n"], rules["S5"]["nominal"])
    print(f"  P(S1-S5 all pass | sound), simulated, S3 judged on {', '.join(st for st in STRATA if rules['S3'][st]['judged']) or 'no stratum'}: "
          f"{rules['joint']:.3f}")
    print()


def scenario_pilot(blind_low: bool, blind_medium_share: float, establishable_share: float) -> dict:
    """A stand-in pilot for the operating characteristics before the pilot: 8 pairs per level
    on both backgrounds, every high pair valid, the low pairs blind (or valid), a share of the
    medium pairs blind; Δ* above the SESOI at the key dose; the cases of a share of every level's
    pairs establishable (every second pair for a half, every fourth for a quarter)."""
    k = P.MAX_PAIRS_PER_LEVEL
    pool = [dict(index=i, level=P.LEVELS[i // k], pair=["a", "b"]) for i in range(3 * k)]
    step = round(1 / establishable_share)  # every step-th pair of each level: the same share per level,
                                           # spread over the level's blind and valid pairs
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
                pilot["establishable"][f"{c.name}:{v}:{i}"] = dict(establishable=(i % k) % step == 0)
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
    print("## Nominal rates (errors by cause, decided 2026-10-08; a sound validator's model after the first review)")
    print(f"  Designed sizes: false SUPPORTED {E_SUPPORTED} (alpha/2: the effect test is two-sided at alpha, the")
    print(f"  direction checked afterwards); false NOT SUPPORTED (metric invalid) on a valid metric {E_INVALID} (2 alpha:")
    print("  GATE 4's upper bound below delta_min alpha/2, GATE 5's negative control alpha, a silent positive")
    print(f"  control shown blind alpha/2); false NO DETECTABLE EFFECT or explained by depth {E_NDE} (alpha, the TOST's")
    print(f"  size); NOT SUPPORTED against a real effect's direction {E_OPPOSITE} (one tail). With pilot.json a card's")
    print("  nominals are a sound validator's on its case (sound_model): GATE 4 at the odds measured on the case,")
    print(f"  GATE 5 at {E_GATE5:g}, then the effect's outcomes as the oracle found them on the case's datasets; a blind,")
    print("  ambiguous or useless metric as the engine's rules gave on them. Without a measurement (the scenarios")
    print(f"  below): the designed sizes of the routes open on the card and decisiveness {D_NOMINAL}. S4 uses the cards'")
    print("  mean nominal error; S3 is judged per stratum (real effects, nulls, blind or useless metrics) where the")
    print("  principle is attainable on its n.")
    print()
    print(f"## S1 as a whole: the number of datasets per key condition ({m} conditions)")
    first = next(n for n in range(400, 2001) if joint_error_rule(n, E_SUPPORTED, m)[1] >= P_PASS_SOUND)
    fails = [n for n in range(first, 1001) if joint_error_rule(n, E_SUPPORTED, m)[1] < P_PASS_SOUND]
    print(f"  the rule first holds at n = {first}; k* moves in steps, so it fails again at n = "
          f"{_ranges(fails)} (up to 1000); the plan's {P.KEY_N} meets it")
    show_joint("S1 false SUPPORTED, key null conditions", P.KEY_N, E_SUPPORTED, m)
    print("  with shared donors the counts are overdispersed; at the same thresholds (beta-binomial, design effect")
    print("  deff; the run reports the realized deff of every primary criterion, above 1.5 a limitation):")
    for deff in (1.0, 1.1, 1.2, 1.5):
        pall, pd = joint_error_rule_deff(P.KEY_N, E_SUPPORTED, m, deff)
        print(f"    deff {deff:.1f}: P(all {m} pass | {E_SUPPORTED:.3f}) = {pall:.3f}, "
              f"P(a condition passes | {2 * E_SUPPORTED:.3f} there) = {pd:.3f}")
    sound = deff_limit(lambda d: joint_error_rule_deff(P.KEY_N, E_SUPPORTED, m, d)[0] >= P_PASS_SOUND)
    doubled = deff_limit(lambda d: joint_error_rule_deff(P.KEY_N, E_SUPPORTED, m, d)[1] <= P_PASS_DOUBLED)
    print(f"    the principle holds up to deff {min(sound, doubled):.3f} (the doubled side up to {doubled:.3f}, "
          f"the sound side up to {sound:.3f})")
    print()
    for title, args_ in (("Scenario A: the low level blind, half the medium pairs blind, half the cases establishable",
                          (True, 0.5, 0.5)),
                         ("Scenario B: every pair valid, every case establishable", (False, 0.0, 1.0)),
                         ("Scenario C: the low level blind, every medium pair blind, a quarter of the cases establishable",
                          (True, 1.0, 0.25))):
        show_design(title, expected_rows(scenario_pilot(*args_)))
    print("## Precision of a rate from n cards (half-width of the 95% CP interval)")
    for n in (100, 200, 400, 790, 1000, 4000):
        row = ", ".join(f"p={q:.3f}: ±{(clopper_pearson(round(q * n), n)[1] - clopper_pearson(round(q * n), n)[0]) / 2:.3f}"
                        for q in (0.025, 0.1, 0.5))
        print(f"  n={n:4d}: {row}")


if __name__ == "__main__":
    main()
