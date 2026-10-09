"""Operating characteristics of the success criteria in validation/prereg/v1.md.

One principle for every criterion (decided by the project owner, 2026-10-08): a sound validator
passes the criterion with probability >= 0.90, and a validator whose error is double the
nominal passes with probability <= 0.05. A criterion is a bound on binomial rates: "at most k*
errors of n" (equivalently: the two-sided 95% Clopper-Pearson upper bound of the rate is <= T),
or for decisiveness "at least k* correct definite outcomes". k* is the most lenient threshold at
which the doubled-error validator still passes with probability <= 0.05. S1 is a conjunction over
the key null conditions, so the principle applies to S1 as a whole; together the criteria pass a
sound validator with probability >= 0.90 (the third round).

Errors are judged by the cause of a verdict, not its label (decided 2026-10-08): every card has
its allowed outcomes (``panel.allowed``). A card's error rates are a sound validator's on its
case (``sound_model``), and a criterion's nominal is, per card, the larger of that rate and the
error's designed size (``rule_nominals``, the second review): the rules are designed to allow that
much. Every kind of wrong outcome is counted where it can occur, per stratum or condition group
(``error_cells``, the third review: errors pooled over unlike cards were absorbed). The n of S2-S7
depend on the truth about the metric on each pair (pilot.json) and on the key's assignment, so
``score.py`` computes every threshold with the rules here from the realized cards; S3's strata are
judged at the tiers fixed before the key (``--write-judged``). This file
shows the rules on the design and on scenarios for the truth, and, given pilot.json, on its
expected cards.

    python validation/prereg/oc.py > validation/prereg/oc.log
    python validation/prereg/oc.py --pilot pilot.json --write-judged > oc_pilot.log
"""
from __future__ import annotations

import argparse
import json

import numpy as np

import panel as P

P_PASS_SOUND = 0.90
P_PASS_DOUBLED = 0.05
P_JOINT_SOUND = 0.90          # P(S1-S7 together | sound) (decided 2026-10-08, third round)
ALPHA = 0.05
E_SUPPORTED = ALPHA / 2       # false SUPPORTED: the two-sided effect test at alpha, one direction
E_INVALID = 2 * ALPHA         # false "metric invalid": GATE 4 (alpha / 2) + GATE 5's negative control
                              # (alpha) + a silent positive control shown blind (alpha / 2)
E_GATE5 = 1.5 * ALPHA         # the GATE 5 share of it (negative control alpha, silent positive control alpha/2)
E_NDE = ALPHA                 # false NO DETECTABLE EFFECT: the TOST's size
E_DEPTH = ALPHA               # false NOT SUPPORTED (explained by depth): the TOST's size as well (with a SESOI
                              # the engine calls depth only where the corrected effect is shown smaller than it)
E_OPPOSITE = ALPHA / 2        # NOT SUPPORTED (opposite direction) where it is an error: one tail of the test
E_CRASH = 0.001               # S7b: engine exceptions and missing reports (the owner's fifth round: a crash gives
                              # no verdict and does not decide validity alone); a sound engine raises none
D_NOMINAL = 0.85              # correct definite outcomes on establishable cards: where no pilot measures it,
                              # and the most S3's thresholds ask for (a measured rate near 1 leaves no doubled
                              # shortfall to tell apart, which would leave the stratum unjudged)
STRATA = ("effect", "null", "invalid")  # S3's strata (decided after the first review)
VALID_STRATA = ("effect", "null")       # S3's strata with a valid metric
GATE0_BOUND = 0.20            # a usefulness bound set in advance (the owner's fifth round): on S3's strata with a
                              # valid metric the sound validator's model takes GATE 0's refusal share at most this,
                              # and a stratum where the frozen engine refuses more often before the key fails S3 (the
                              # tool is silent too often on a good metric); the run goes on. On the blind and useless
                              # metrics the measured share (a refusal there protects from a false conclusion)
TIERS = ("principle", "floor")          # S3's thresholds, strictest first (the second review)
# N3 and N8: the oracle undoes their planted nuisance with the true correction, which the engine
# cannot have (it may refuse there, an allowed outcome), so their decisiveness is reported, not
# judged (the second review)
OUTSIDE_STRATA = ("N3", "N8")
MAX_RAISED = 0.5              # the largest nominal an error criterion is raised to (``error_criterion``)
EPS_CELL = 0.0005             # a sound validator fails an error cell with probability <= this (the 16 cells
                              # together <= 0.008 of the 0.016 that S1 at its design leaves of the joint 0.90)
DATA_STRATA = ("effect", "null", "invalid")  # the strata of the error criteria: every card, by the truth
S2_GROUPS = ("N2", "N3", "N4", "N6", "N7", "E")  # S2's groups: the cards outside S1, by condition
# the error criteria's cells: (criterion, outcome counted, the strata or groups it is counted in);
# every wrong outcome but an unexpected verdict (S7a) is counted in exactly one cell of a card, and NOT
# SUPPORTED against the direction on null data, allowed there, in its own (the fourth review)
CAUSE_FIELD = {P.SUPPORTED: "sup", P.NS_INVALID: "inv", P.NDE: "nde", P.NS_OPPOSITE: "opp", P.NS_DEPTH: "depth"}
FIELDS = (*CAUSE_FIELD.values(), "oppnull")  # "oppnull": NOT SUPPORTED against the direction on null data
DESIGNED = {P.SUPPORTED: E_SUPPORTED, P.NS_INVALID: E_INVALID, P.NDE: E_NDE, P.NS_OPPOSITE: E_OPPOSITE,
            P.NS_DEPTH: E_DEPTH}
ERROR_CRITERIA = ("S2", "S4", "S5", "S6")
CRITERIA = ("S1", "S2", "S3", "S4", "S5", "S6", "S7a", "S7b")
# S7 (split in the fifth round). S7a, the verdicts the engine's deterministic rules cannot give here:
# these (a label that does not match its cause or a cause the engine does not have, a cause the panel
# does not expect, UNIDENTIFIABLE), DEGENERATE METRIC on a metric that varies, and an effect verdict
# where its rules exclude one (panel.excluded): none allowed. S7b, crashes (P.CRASH: an engine
# exception or no report): at most the principle's threshold at E_CRASH
S7_OUTCOMES = (P.ERROR, P.OTHER, P.UNIDENTIFIABLE)
# the flag of a card's model that says where a wrong outcome can occur (``card_model``)
POSSIBLE = {P.SUPPORTED: "sup_error_possible", P.NS_INVALID: "invalid_error_possible", P.NDE: "nde_error_possible",
            P.NS_OPPOSITE: "opp_error_possible", P.NS_DEPTH: "depth_error_possible"}


def nominal_error(ok: frozenset) -> float:
    """The designed sizes of the routes to an outcome outside the allowed set `ok`: SUPPORTED
    alpha/2, NOT SUPPORTED (metric invalid) 2 alpha, NO DETECTABLE EFFECT alpha, NOT SUPPORTED
    (opposite direction) alpha/2, NOT SUPPORTED (explained by depth) alpha."""
    return sum(size for o, size in DESIGNED.items() if o not in ok)


def data_stratum(cond: P.Condition, truth: str) -> str:
    """The stratum of a card for the error criteria (every card, unlike S3's): a blind or useless
    metric; else real effects or nulls, by the data."""
    if truth in ("blind", "useless", "constant"):
        return "invalid"
    return "effect" if cond.data == "effect" else "null"


def s2_group(cond: P.Condition, variant: str) -> str | None:
    """S2's group of a card outside S1 (the third review: false SUPPORTED concentrated in one
    condition, N4's pseudoreplication or N7, passed the pooled S2): N2's steps, N3, N4, the useless
    metrics N6a and N6b, N7, and the real effects (whose blind or useless cards can err). None for
    S1's cards. (N4's card without the replicate unit and N6b's constant give no effect verdict:
    theirs is in S7a, ``card_kind``.)"""
    if cond.key and variant == cond.key_variant:
        return None
    if cond.data == "effect":
        return "E"
    return "N6" if cond.name in ("N6a", "N6b") else cond.name


def error_cells(rows: list[dict]) -> dict:
    """The cells of S2, S4, S5 and S6 on a set of cards: per cell, the outcome it counts, the field
    of its rates (``FIELDS``) and its member cards (indices). S2: a false SUPPORTED per group of the
    cards outside S1; S4: a false NOT SUPPORTED for another cause than the metric (opposite
    direction; explained by depth) per stratum, and on null data, where NOT SUPPORTED against the
    direction is allowed (the claim is false), its rate as a false detection's (the fourth review: a
    validator that said it on 16% of the null cards passed); S5: a false "metric invalid" on a valid
    metric per stratum (real effects, nulls); S6: a false NO DETECTABLE EFFECT per stratum (real
    effects at or above the SESOI, blind or useless metrics). An outcome the engine's rules exclude
    on a card is in none (``card_kind``: S7a)."""
    cells = {}

    def add(name, outcome, field, member):
        cells[name] = dict(outcome=outcome, field=field, members=[i for i, r in enumerate(rows) if member(r)])
    for g in S2_GROUPS:
        add(f"S2:{g}", P.SUPPORTED, "sup", lambda r, g=g: r["sup_error_possible"] and r["s2_group"] == g)
    for st in DATA_STRATA:
        if st == "null":  # allowed there, its rate bounded at a false detection's
            add("S4:null:opposite", P.NS_OPPOSITE, "oppnull", lambda r: r["opp_null_counted"])
        else:
            add(f"S4:{st}:opposite", P.NS_OPPOSITE, "opp",
                lambda r, st=st: r["opp_error_possible"] and r["data_stratum"] == st)
        add(f"S4:{st}:depth", P.NS_DEPTH, "depth",
            lambda r, st=st: r["depth_error_possible"] and r["data_stratum"] == st)
    for st in ("effect", "null"):
        add(f"S5:{st}", P.NS_INVALID, "inv", lambda r, st=st: r["valid"] and r["data_stratum"] == st)
    for st in ("effect", "invalid"):
        add(f"S6:{st}", P.NDE, "nde", lambda r, st=st: r["nde_error_possible"] and r["data_stratum"] == st)
    return cells


def card_kind(row: dict, cond: P.Condition, k: int) -> dict:
    """The row of a dataset's k-th claim card from its case's (``card_model``). The engine's rules
    exclude an effect verdict on some cards (``panel.excluded``: N4's first card, which has no
    replicate unit, and the constant metric): there a sound validator never gives one, so its rates
    are 0, and one given is an unexpected verdict, counted in S7a (none allowed), not in the cells
    (the fourth review: a validator that said SUPPORTED on 10% of N4's cards without the replicate
    unit passed the N4 cell, whose small n needed a raised nominal). The outcomes left are in the
    sound validator's other outcomes (INCONCLUSIVE: insufficient replication)."""
    gone = P.excluded(cond, k)
    if not gone:
        return dict(row)
    out = dict(row, opp_null_counted=False, p_oppnull=0.0, r_oppnull=0.0)
    removed = 0.0
    for o in gone:
        f = CAUSE_FIELD[o]
        removed += out[f"p_{f}"]
        out.update({f"p_{f}": 0.0, f"r_{f}": 0.0, POSSIBLE[o]: False})
    out["p_err"] = max(0.0, out["p_err"] - removed)
    return out


def stratum(cond: P.Condition, truth: str) -> str | None:
    """S3's stratum of a card: real effects with a valid metric (the correct definite outcome is
    SUPPORTED, or NO DETECTABLE EFFECT below the SESOI); null data with a valid metric (NO
    DETECTABLE EFFECT or NOT SUPPORTED); a blind or useless metric (metric invalid, or DEGENERATE
    for the constant). Ambiguous metrics belong to none (either kind of verdict is correct there),
    nor do N3 and N8 (``OUTSIDE_STRATA``)."""
    if cond.name in OUTSIDE_STRATA:
        return None
    if truth in ("blind", "useless", "constant"):
        return "invalid"
    if truth == "valid":
        return "effect" if cond.data == "effect" else "null"
    return None


def _shares(counts: dict) -> tuple[dict, int]:
    n = sum(counts.values())
    return ({k: v / n for k, v in counts.items()} if n else {}), n


def outcome_distribution(cond: P.Condition, variant: str, pair: int, pilot: dict) -> dict | None:
    """The outcomes of a sound validator — one that follows the engine's pre-registered rules
    correctly — on a card, as probabilities, from the pilot's measurements of its case (None where
    the pilot measured nothing). The engine's verdict order (report.decide_cause): GATE 5's failure,
    then GATE 4's, is "metric invalid"; then an effect explained by depth, before the metric's
    validity is asked for; then the effect's outcome where GATE 4 PASSed, INCONCLUSIVE where it was
    UNTESTED (``oracle.engine_outcome``).

    * a valid metric: GATE 5 FAILs with its designed 1.5 alpha, GATE 4 with the case's measured odds
      (pilot.json "gate4", from the metric's response on the case's datasets); after a PASS the
      effect's outcomes as the oracle found them on the case's datasets (the metric's validity taken
      as known), after an UNTESTED its "explained by depth" share, else INCONCLUSIVE;
    * a blind or ambiguous metric: GATE 5 FAILs with 1.5 alpha ("metric invalid", allowed there),
      otherwise the engine's order on the case's datasets, where the oracle ran GATE 4's rule;
    * a useless metric: the engine's order on its datasets (GATE 5's controls are pairs: not modelled);
    * the constant: DEGENERATE METRIC.

    GATE 0 is not modelled here: its refusal is allowed on every card and never an error, and the
    engine decides it first, so ``sound_model`` takes it as a share of the card's S3 stratum
    (``refusal_share``)."""
    truth = P.metric_truth(cond, variant, pair, pilot)
    est = (pilot.get("establishable") or {}).get(f"{cond.name}:{variant}:{pair}") or {}
    if truth == "constant":
        return {P.DEGENERATE: 1.0}
    if truth == "valid" and est.get("effect_outcomes"):
        odds = P.truth_record(cond, variant, pair, pilot).get("gate4")
        if odds is None:
            return None
        f, _ = _shares(est["effect_outcomes"])
        g5, p_pass, p_fail = E_GATE5, float(odds["p_pass"]), float(odds["p_fail"])
        p_untested = max(0.0, 1.0 - p_pass - p_fail)
        dist = {P.NS_INVALID: g5 + (1 - g5) * p_fail}
        for o, v in f.items():
            dist[o] = dist.get(o, 0.0) + (1 - g5) * p_pass * v
        depth = f.get(P.NS_DEPTH, 0.0)
        dist[P.NS_DEPTH] = dist.get(P.NS_DEPTH, 0.0) + (1 - g5) * p_untested * depth
        dist[P.INCONCLUSIVE] = dist.get(P.INCONCLUSIVE, 0.0) + (1 - g5) * p_untested * (1 - depth)
        return dist
    if truth in ("blind", "ambiguous", "useless") and est.get("engine_outcomes"):
        g, _ = _shares(est["engine_outcomes"])
        g5 = E_GATE5 if cond.metric == "norm_pearson" else 0.0
        dist = {o: (1 - g5) * v for o, v in g.items()}
        dist[P.NS_INVALID] = dist.get(P.NS_INVALID, 0.0) + g5
        return dist
    return None


def refusal_share(cond: P.Condition, truth: str, pilot: dict) -> float:
    """GATE 0's refusal share on the card's S3 stratum, measured before the key with the frozen
    engine (pilot.json "gate0_refusals", ``timing.py --refusals``: the fourth review found a sound
    engine failing S3 at a 15% share the model left out), at most GATE0_BOUND on the strata with a
    valid metric (the fifth round: a GATE 0 that refused half the time halved the nominal and still
    passed S3; above the bound the stratum fails, ``refusal_bound_failures``); 0 outside the strata
    and where nothing was measured."""
    st = stratum(cond, truth)
    rec = (pilot.get("gate0_refusals") or {}).get(st) if st else None
    share = float(rec["share"]) if rec and rec.get("n") else 0.0
    return min(share, GATE0_BOUND) if st in VALID_STRATA else share


def refusal_bound_failures(pilot: dict) -> list[str]:
    """S3's strata with a valid metric on which the frozen engine's GATE 0 refused more often than
    GATE0_BOUND before the key (pilot.json "gate0_refusals"): S3 fails on each, whatever the run
    gives (decided before the key; the run goes on, its other measurements are needed)."""
    recs = pilot.get("gate0_refusals") or {}
    return [st for st in VALID_STRATA if (recs.get(st) or {}).get("n") and float(recs[st]["share"]) > GATE0_BOUND]


def sound_model(cond: P.Condition, variant: str, pair: int, pilot: dict) -> dict:
    """A sound validator's rates on a card (``outcome_distribution``), per wrong outcome where it is
    wrong: a false SUPPORTED (`p_sup`), "metric invalid" on a valid metric (`p_inv`), NO DETECTABLE
    EFFECT (`p_nde`), NOT SUPPORTED against the direction (`p_opp`) and explained by depth
    (`p_depth`); any outcome outside the allowed set (`p_err`, all of them included) and a correct
    definite outcome (`decisive`); `measured` says whether the pilot measured them. On null data,
    where NOT SUPPORTED against the direction is allowed, its rate too (`p_oppnull`, a false
    detection in the other tail: the fourth review). Without a measurement: the designed sizes of the
    routes the card's allowed outcomes leave open (``DESIGNED``, ``nominal_error``), alpha/2 for
    `p_oppnull`, and decisiveness D_NOMINAL, or what the errors leave of it. GATE 0 refuses first,
    with its stratum's measured share (``refusal_share``, `refusal`; at most GATE0_BOUND where the
    metric is valid), and every other outcome takes the rest; the constant metric's DEGENERATE
    METRIC comes before it."""
    ok = P.allowed(cond, variant, pair, pilot)
    good = P.definite(cond, variant, pair, pilot)
    truth = P.metric_truth(cond, variant, pair, pilot)
    valid = truth == "valid"
    opp_null = data_stratum(cond, truth) == "null" and P.NS_OPPOSITE in ok
    r0 = 0.0 if truth == "constant" else refusal_share(cond, truth, pilot)
    dist = outcome_distribution(cond, variant, pair, pilot)
    if dist is not None:
        dist = {o: (1 - r0) * v for o, v in dist.items()}
        if r0:
            dist[P.REFUSAL] = dist.get(P.REFUSAL, 0.0) + r0
        rates = {f"p_{f}": (dist.get(o, 0.0) if (o not in ok and (o != P.NS_INVALID or valid)) else 0.0)
                 for o, f in CAUSE_FIELD.items()}
        return dict(rates, opp_null_counted=opp_null, p_oppnull=dist.get(P.NS_OPPOSITE, 0.0) if opp_null else 0.0,
                    p_err=sum(v for o, v in dist.items() if o not in ok),
                    decisive=sum(v for o, v in dist.items() if o in good), measured=True, refusal=r0)
    err = nominal_error(ok)
    rates = {f"p_{f}": (1 - r0) * (DESIGNED[o] if (o not in ok and (o != P.NS_INVALID or valid)) else 0.0)
             for o, f in CAUSE_FIELD.items()}
    return dict(rates, opp_null_counted=opp_null, p_oppnull=(1 - r0) * E_OPPOSITE if opp_null else 0.0,
                p_err=(1 - r0) * err, decisive=(1 - r0) * min(D_NOMINAL, 1.0 - err), measured=False, refusal=r0)


def rule_nominals(model: dict, ok: frozenset, valid: bool) -> dict:
    """A card's nominals for the thresholds of S2, S4, S5 and S6: per wrong outcome where it is
    wrong, the larger of the sound validator's rate (``sound_model``) and the outcome's designed
    size (``DESIGNED``: a false SUPPORTED alpha/2, "metric invalid" on a valid metric 2 alpha, NO
    DETECTABLE EFFECT alpha, NOT SUPPORTED against the direction alpha/2, explained by depth alpha).
    The second review: a measured rate far below the designed size can leave too few expected errors
    to tell a doubled rate apart; the design allows the designed size, and the validation asks no
    more of it. Each outcome has its own nominal (the third review: one nominal for any error let a
    validator trade one kind of error for another); NOT SUPPORTED against the direction on null data
    at least alpha/2 (`r_oppnull`)."""
    out = {f"r_{f}": (max(DESIGNED[o], model[f"p_{f}"]) if (o not in ok and (o != P.NS_INVALID or valid)) else 0.0)
           for o, f in CAUSE_FIELD.items()}
    out["r_oppnull"] = max(E_OPPOSITE, model["p_oppnull"]) if model["opp_null_counted"] else 0.0
    return out


def card_model(cond: P.Condition, variant: str, pair: int, pilot: dict) -> dict:
    """What the criteria need of a card before its outcome is known: where each error can occur,
    its stratum and establishability, the sound validator's rates and the rule nominals."""
    ok = P.allowed(cond, variant, pair, pilot)
    truth = P.metric_truth(cond, variant, pair, pilot)
    est = (pilot.get("establishable") or {}).get(f"{cond.name}:{variant}:{pair}") or {}
    model = sound_model(cond, variant, pair, pilot)
    return dict(condition=cond.name, key=cond.key and variant == cond.key_variant,
                sup_error_possible=P.SUPPORTED not in ok, invalid_error_possible=P.NS_INVALID not in ok,
                nde_error_possible=P.NDE not in ok, opp_error_possible=P.NS_OPPOSITE not in ok,
                depth_error_possible=P.NS_DEPTH not in ok, degenerate_allowed=P.DEGENERATE in ok,
                valid=truth == "valid", stratum=stratum(cond, truth), data_stratum=data_stratum(cond, truth),
                s2_group=s2_group(cond, variant), establishable=cond.oracle and bool(est.get("establishable")),
                **model, **rule_nominals(model, ok, truth == "valid"))


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


def binom_sf(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p)."""
    return 1.0 - binom_cdf(k - 1, n, p)


def meets(p_sound: float, p_doubled: float) -> bool:
    return p_sound >= P_PASS_SOUND and p_doubled <= P_PASS_DOUBLED


def error_rule(n: int, p0: float):
    """For an error-rate criterion: the largest k* with P(X <= k* | 2 p0) <= 0.05, and
    P(pass | p0), P(pass | 2 p0) for the rule 'at most k* errors' (with p0 = 0, where no doubled
    rate can be told apart, k* = 0: any such error fails; k* = -1 where even no error is too likely
    for the doubled rate). For a mixture of cards with different nominal rates, p0 is their mean;
    the count is then Poisson-binomial, more concentrated than the binomial (Hoeffding 1956), so
    the binomial rule is conservative on both sides."""
    p1 = min(1.0, 2 * p0)
    if n <= 0:
        return -1, 1.0, 1.0
    if p0 <= 0:  # the sound validator cannot make this error on these cards: none is allowed
        return 0, 1.0, 1.0
    from scipy import stats
    under = np.nonzero(stats.binom.cdf(np.arange(n + 1), n, p1) <= P_PASS_DOUBLED)[0]
    k = int(under[-1]) if len(under) else -1  # the cdf rises with k: the last k still at or under 0.05
    return k, binom_cdf(k, n, p0), binom_cdf(k, n, p1)


def error_criterion(n: int, nominal: float, sound: float | None = None, eps: float = EPS_CELL) -> dict:
    """An error criterion (a cell) on n cards at its nominal (the cards' mean rule nominal): 'at
    most k* errors' (``error_rule``). The nominal is raised, to the first rate of the grid
    nominal x 1.01^j, where the principle is not attainable at it on n (too few expected errors for a
    doubled rate to be told apart: the second review found P(pass | sound) = 0.51 at 0.001 on 5,600
    cards) or where a sound validator — its cards' mean rate `sound` — would fail the cell with
    probability above `eps` (the third review's cells must leave the joint requirement room); where
    no rate up to MAX_RAISED does, the cell is reported only. Without cards it holds (no card on
    which its error can occur)."""
    out = dict(n=int(n), nominal=float(nominal), rule_nominal=float(nominal), raised=False,
               sound=None if sound is None else float(sound))
    if n <= 0:
        return dict(out, judged=False, max_allowed=None, p_pass_sound=1.0, p_pass_doubled=1.0)
    p = float(nominal)
    while True:
        k, ps, pd = error_rule(n, p)
        sound_ok = sound is None or binom_cdf(k, n, sound) >= 1 - eps
        if p <= 0 or (meets(ps, pd) and sound_ok):
            return dict(out, rule_nominal=p, raised=p > nominal, judged=True, max_allowed=k,
                        p_pass_sound=ps, p_pass_doubled=pd,
                        p_pass_measured=None if sound is None else binom_cdf(k, n, sound))
        if p * 1.01 > MAX_RAISED:
            return dict(out, judged=False, max_allowed=None, p_pass_sound=ps, p_pass_doubled=pd)
        p *= 1.01


def joint_error_rule(n: int, p0: float, m: int):
    """S1 as a whole over m conditions of n datasets each: the per-condition k* (as in
    `error_rule`), P(all m pass | p0 everywhere), P(one condition passes | 2 p0 there)."""
    k, ps, pd = error_rule(n, p0)
    return k, ps ** m, pd, ps


def tier_alternative(d0: float, tier: str) -> float:
    """The decisiveness S3's tier must fail with probability >= 0.95: the principle's doubled
    shortfall 1 - 2 (1 - d0), but at least half the nominal (below d0 = 2/3 doubling the shortfall
    leaves less than half, and below 0.5 nothing: the second review found the rule failing there);
    the floor: half the nominal (a validator half as decisive as a sound one)."""
    return max(2 * d0 - 1, d0 / 2) if tier == "principle" else d0 / 2


def decisiveness_rule(n: int, d0: float, alt: float | None = None):
    """For decisiveness (a share that should be high) at the nominal d0 against the alternative
    `alt` (the principle's tier, ``tier_alternative``, by default): the smallest k* with
    P(X >= k* | alt) <= 0.05, and P(pass | d0), P(pass | alt) for the rule 'at least k* correct
    definite outcomes'. (None, 0, 0) where there is no card or no alternative to tell apart."""
    if alt is None and n > 0 and np.isfinite(d0):
        alt = tier_alternative(d0, "principle")
    if n <= 0 or not np.isfinite(d0) or alt is None or not 0 < alt < 1:
        return None, 0.0, 0.0
    from scipy import stats
    tail = stats.binom.sf(np.arange(-1, n + 1), n, alt)  # tail[k] = P(X >= k), k = 0 .. n + 1
    k = int(np.nonzero(tail <= P_PASS_DOUBLED)[0][0])     # P(X >= n + 1) = 0: there is one
    return k, binom_sf(k, n, d0), binom_sf(k, n, alt)


def decisiveness_tiers(n: int, d0: float) -> dict:
    """S3's two tiers on a stratum of n establishable cards at the nominal d0: each tier's
    threshold, its operating characteristics and whether the principle is attainable with it."""
    out = {}
    for tier in TIERS:
        alt = tier_alternative(d0, tier) if np.isfinite(d0) else float("nan")
        k, ps, pa = decisiveness_rule(n, d0, alt) if np.isfinite(d0) else (None, 0.0, 0.0)
        out[tier] = dict(alternative=alt, min_required=k, p_pass_sound=ps, p_pass_alternative=pa,
                         attainable=k is not None and meets(ps, pa))
    return out


def _mean(rows: list[dict], field: str) -> float:
    return float(np.mean([r[field] for r in rows])) if rows else 0.0


def criteria_rules(rows: list[dict], sims: int = 20_000, seed: int = 0, s3_tiers: dict | None = None) -> dict:
    """The thresholds of S1-S7 on a set of cards (as ``score.card_rows`` and ``expected_rows`` give
    them):

    * S1 per key condition at alpha/2 (790 datasets each, at most 29: the design);
    * S2, S4, S5 and S6 per cell (``error_cells``: a false SUPPORTED per group of the cards outside
      S1; a false NOT SUPPORTED against the direction or explained by depth, "metric invalid" on a
      valid metric and NO DETECTABLE EFFECT, each per stratum), every cell at the mean of its cards'
      rule nominals for its outcome (``rule_nominals``; ``error_criterion``); a criterion holds if
      every cell holds;
    * S7a (a verdict the engine's deterministic rules cannot give, all cards): none allowed — a
      sound validator makes none, and any rate at which one makes them is told apart;
    * S7b (crashes: an engine exception or a missing report, all cards): at most the principle's
      threshold at the nominal E_CRASH (raised as a cell's where the principle is not attainable at
      it; where it is not attainable at all, none allowed); a sound engine crashes on none;
    * S3 per stratum at the mean measured decisiveness of its establishable cards, at most D_NOMINAL,
      at one of two tiers (``decisiveness_tiers``): the principle's, or the floor's (a validator half
      as decisive fails with probability >= 0.95). Every stratum is judged (the second review: a
      validator that never says SUPPORTED passed where the real effects were reported only): at the
      floor where every stratum's floor is attainable and P(S1-S7 together | sound) stays >= 0.90,
      then at the principle's tier in the order of STRATA where it is attainable and the joint
      requirement still holds. Where not every stratum fits, the strata that fit are judged in that
      order and S3 fails: decisiveness is not shown. `s3_tiers` (pilot.json "s3_rules", fixed before
      the key) gives the tiers instead.

    The joint requirement is checked with S1 at its design (a false SUPPORTED at alpha/2 on every key
    card, P(S1 | sound) = 0.915): S1's thresholds are the design's, and a measured rate on the key
    cards that differs is a finding about the engine, not a reason to judge S3 differently. `joint`
    is P(S1-S7 together | sound) so, `joint_measured` the same with every rate as measured."""
    out = dict(S1={}, S3={}, cells={})
    for c in P.CONDITIONS:
        if c.key:
            n = sum(1 for r in rows if r["key"] and r["condition"] == c.name)
            out["S1"][c.name] = dict(n=n, max_allowed=error_rule(n, E_SUPPORTED)[0] if n else None,
                                     nominal=E_SUPPORTED)
    for name, cell in error_cells(rows).items():
        field = cell["field"]
        members = [rows[i] for i in cell["members"]]
        rule = error_criterion(len(members), _mean(members, f"r_{field}"), _mean(members, f"p_{field}"))
        out["cells"][name] = dict(rule, outcome=cell["outcome"], field=field, members=cell["members"])
    out["S7a"] = dict(n=len(rows), nominal=0.0, max_allowed=0, judged=bool(rows))
    crash = error_criterion(len(rows), E_CRASH, 0.0)
    out["S7b"] = crash if crash["judged"] or not rows else dict(crash, max_allowed=0)
    for st in STRATA:
        rs = [r for r in rows if r["establishable"] and r["stratum"] == st]
        measured = _mean(rs, "decisive") if rs else float("nan")
        d0 = min(measured, D_NOMINAL) if rs else float("nan")
        out["S3"][st] = dict(n=len(rs), measured=measured, nominal=d0, tiers=decisiveness_tiers(len(rs), d0),
                             tier=None, min_required=None)
    if not rows:
        out.update(joint=float("nan"), joint_measured=float("nan"), s3_feasible=False)
        return out
    design = _simulate_passes(_at_design(rows), out, sims, seed)
    tiers = dict(s3_tiers) if s3_tiers is not None else _select_tiers(design, out)
    for st in STRATA:
        t = tiers.get(st)
        out["S3"][st]["tier"] = t
        out["S3"][st]["min_required"] = out["S3"][st]["tiers"][t]["min_required"] if t else None
    out["s3_feasible"] = all(out["S3"][st]["min_required"] is not None for st in STRATA)
    out["joint"] = float(np.mean(_together(design, out)))
    measured = _simulate_passes(rows, out, sims, seed) if any(r.get("measured") for r in rows) else design
    out["joint_measured"] = float(np.mean(_together(measured, out)))
    return out


def _at_design(rows: list[dict]) -> list[dict]:
    """The cards with S1 at its design: a false SUPPORTED at alpha/2 on every key card."""
    out = []
    for r in rows:
        if r["key"] and r["sup_error_possible"]:
            others = sum(r[f"p_{f}"] for f in CAUSE_FIELD.values() if f != "sup")
            p_err = max(r["p_err"] - r["p_sup"] + E_SUPPORTED, E_SUPPORTED + others)
            r = dict(r, p_sup=E_SUPPORTED, p_err=p_err, decisive=min(r["decisive"], 1.0 - p_err))
        out.append(r)
    return out


def _base(passes: dict) -> np.ndarray:
    """S1, every cell of S2 and S4-S6, S7a and S7b, per simulation."""
    ok = passes["S1"] & passes["S7a"] & passes["S7b"]
    for name, v in passes.items():
        if name.split(":")[0] in ERROR_CRITERIA:
            ok = ok & v
    return ok


def _together(passes: dict, rules: dict) -> np.ndarray:
    """Whether S1, S2, S4-S7 and S3's judged strata all pass, per simulation (a stratum without a
    tier fails S3)."""
    ok = _base(passes)
    for st in STRATA:
        t = rules["S3"][st]["tier"]
        if t is None or rules["S3"][st]["min_required"] is None:
            return np.zeros_like(ok)
        ok = ok & passes[f"S3:{st}:{t}"]
    return ok


def _select_tiers(passes: dict, rules: dict) -> dict:
    """S3's tiers by the rule of ``criteria_rules``, on simulated passes of a sound validator."""
    base = _base(passes)
    tier_ok = {(st, t): rules["S3"][st]["tiers"][t]["attainable"] for st in STRATA for t in TIERS}
    if all(tier_ok[(st, "floor")] for st in STRATA):
        chosen = {st: "floor" for st in STRATA}
        if np.mean(base & np.logical_and.reduce([passes[f"S3:{st}:floor"] for st in STRATA])) >= P_JOINT_SOUND:
            for st in STRATA:
                if not tier_ok[(st, "principle")]:
                    continue
                trial = dict(chosen, **{st: "principle"})
                together = base & np.logical_and.reduce([passes[f"S3:{s}:{t}"] for s, t in trial.items()])
                if np.mean(together) >= P_JOINT_SOUND:
                    chosen = trial
            return chosen
    chosen, acc = {st: None for st in STRATA}, base
    for st in STRATA:
        for t in TIERS:
            if tier_ok[(st, t)] and np.mean(acc & passes[f"S3:{st}:{t}"]) >= P_JOINT_SOUND:
                chosen[st], acc = t, acc & passes[f"S3:{st}:{t}"]
                break
    return chosen


def _simulate_passes(rows: list[dict], rules: dict, sims: int, seed: int) -> dict:
    """For a sound validator, by simulation on the cards `rows` with each card's model
    (``sound_model``): on every card each wrong outcome where it is wrong with its rate (p_sup,
    p_inv, p_nde, p_opp, p_depth), another outcome outside the allowed set with the rest of p_err,
    and on an establishable card otherwise a correct definite outcome with `decisive`, of which NOT
    SUPPORTED against the direction on null data with p_oppnull (allowed, its rate bounded); never
    a rule violation or a crash (S7a, S7b). The criteria share cards, so they are simulated on the
    same draws: returns, per criterion (S1 as a whole, every cell of S2 and S4-S6, S7a, S7b, every S3
    stratum at every tier with a threshold), whether it passes in each simulation."""
    rng = np.random.default_rng(seed)
    n = len(rows)
    fields = list(CAUSE_FIELD.values())
    rates = np.array([[r[f"p_{f}"] for f in fields] for r in rows]).reshape(n, len(fields))
    upper = np.cumsum(rates, axis=1)          # outcome j of a card: u in [upper[j - 1], upper[j])
    lower = upper - rates
    p_err = np.maximum(np.array([r["p_err"] for r in rows]), upper[:, -1] if n else 0.0)
    p_good = np.minimum(np.array([r["decisive"] if r["establishable"] else 0.0 for r in rows]), 1 - p_err)
    # NOT SUPPORTED against the direction on null data: allowed and definite, so drawn first in the
    # correct definite outcomes' share where the card is establishable, else after the errors
    p_oppnull = np.array([r.get("p_oppnull", 0.0) for r in rows])
    key_groups = {}
    for i, r in enumerate(rows):
        if r["key"]:
            key_groups.setdefault(r["condition"], []).append(i)
    s3 = {f"S3:{st}:{t}": (np.array([r["establishable"] and r["stratum"] == st for r in rows]), v["min_required"])
          for st in STRATA for t, v in rules["S3"][st]["tiers"].items() if v["min_required"] is not None}
    cells = rules.get("cells", {})
    out = {name: [] for name in ("S1", "S7a", "S7b", *cells, *s3)}
    col = {f: j for j, f in enumerate(fields)}
    chunk = 500
    for start in range(0, sims, chunk):
        m = min(chunk, sims - start)
        u = rng.random((m, n))
        hit = {f: (u >= lower[:, j]) & (u < upper[:, j]) for f, j in col.items()}
        good = (u >= p_err) & (u < p_err + p_good)
        hit["oppnull"] = (u >= p_err) & (u < p_err + p_oppnull)
        s1 = np.ones(m, bool)
        for c, ix in key_groups.items():
            s1 &= hit["sup"][:, ix].sum(axis=1) <= rules["S1"][c]["max_allowed"]
        out["S1"].append(s1)
        for name, rule in cells.items():
            ix = rule["members"]
            if rule.get("judged") and ix:
                out[name].append(hit[rule["field"]][:, ix].sum(axis=1) <= rule["max_allowed"])
            else:
                out[name].append(np.ones(m, bool))
        out["S7a"].append(np.ones(m, bool))
        out["S7b"].append(np.ones(m, bool))
        for name, (member, k) in s3.items():
            out[name].append((good & member).sum(axis=1) >= k)
    return {k: np.concatenate(v) for k, v in out.items()}


def joint_pass_probability(rows: list[dict], sims: int = 20_000, seed: int = 0) -> float:
    """P(S1-S7 all pass) for a sound validator on the cards `rows`, S1 at its design
    (``criteria_rules``)."""
    return criteria_rules(rows, sims, seed)["joint"]


def error_rule_deff(n: int, p0: float, deff: float):
    """`error_rule`'s threshold k* (of the independent design) when the count of errors is
    overdispersed by a design effect `deff`: a beta-binomial with the binomial's mean and deff times
    its variance (intra-class correlation (deff - 1) / (n - 1)). Given the background the datasets
    are independent draws and the count is binomial (the third review); this is the count read
    beyond the background's donors, as if they were a random sample of donors of the same kind.
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
    truth, allowed outcomes and establishability of pilot.json (``card_model``). A stand-in for the
    realized cards before the key. Conditions dropped, or whose background is (no pool), have none."""
    rows = []
    for c in P.CONDITIONS:
        size = (pilot.get("pool_size") or {}).get(c.background)
        if not size:
            continue
        for variant, n in c.variants:
            if P._dropped(c.name, variant, dropped):
                continue
            for i in range(n):
                row = card_model(c, variant, i % int(size), pilot)
                rows += [card_kind(row, c, k) for k in range(c.cards)]
    return rows


def show_error(name: str, n: int, p0: float, k: int | None = None):
    k_, ps, pd = error_rule(n, p0)
    k = k_ if k is None else k
    hi = clopper_pearson(k, n)[1] if n and k >= 0 else float("nan")
    print(f"  {name}: n = {n}, pass with at most {k} errors (upper 95% CP bound <= {hi:.4f}); "
          f"P(pass | {p0:.4f}) = {ps:.3f}, P(pass | {min(1, 2 * p0):.4f}) = {pd:.3f} — "
          f"{'meets' if meets(ps, pd) else 'FAILS'} the principle")


def detectable_rate(k: int, n: int, power: float = 0.95) -> float:
    """A cell's resolution: the smallest error rate at which 'at most k errors of n' fails with
    probability >= `power` (bisection on binom_cdf; the fourth review: a small cell, whose nominal is
    raised, passes targeted error rates well above its nominal)."""
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if binom_cdf(k, n, mid) <= 1 - power:
            hi = mid
        else:
            lo = mid
    return hi


def show_criterion(name: str, rule: dict):
    if rule["n"] == 0:
        print(f"  {name}: no card on which the error can occur — holds")
        return
    if not rule["judged"]:
        print(f"  {name}: n = {rule['n']}, nominal {rule['nominal']:.4f}: the principle is not attainable up to "
              f"{MAX_RAISED} — reported only")
        return
    raised = (f" (the cards' mean rule nominal {rule['nominal']:.4f}, raised to {rule['rule_nominal']:.4f}: there "
              f"the principle is attainable and a sound validator fails the cell with probability <= {EPS_CELL})"
              if rule["raised"] else "")
    show_error(name + raised, rule["n"], rule["rule_nominal"], rule["max_allowed"])
    if rule.get("sound") is not None and rule.get("p_pass_measured") is not None:
        print(f"    a sound validator (the cards' mean rate {rule['sound']:.4f}) passes it with "
              f"{rule['p_pass_measured']:.4f}; it fails an error rate of "
              f"{detectable_rate(rule['max_allowed'], rule['n']):.4f} or more with probability >= 0.95")


def show_joint(name: str, n: int, p0: float, m: int):
    k, pall, pd, pone = joint_error_rule(n, p0, m)
    hi = clopper_pearson(k, n)[1]
    ok = "meets" if pall >= P_PASS_SOUND and pd <= P_PASS_DOUBLED else "FAILS"
    print(f"  {name}: {m} conditions x n = {n}, each passes with at most {k}/{n} false SUPPORTED (upper "
          f"95% CP bound <= {hi:.4f}); P(one passes | {p0:.3f}) = {pone:.4f}, P(all {m} pass | "
          f"{p0:.3f}) = {pall:.4f}, P(a condition passes | {2 * p0:.3f} there) = {pd:.4f} — {ok} the principle")


def show_decisive(name: str, n: int, d0: float, alt: float | None = None):
    k, ps, pa = decisiveness_rule(n, d0, alt)
    alt = tier_alternative(d0, "principle") if alt is None else alt
    if k is None:
        print(f"  {name}: n = {n}: no threshold")
        return
    lo = clopper_pearson(k, n)[0] if n else float("nan")
    print(f"  {name}: n = {n}, pass with at least {k} correct definite outcomes (lower 95% CP bound "
          f">= {lo:.4f}); P(pass | {d0:.3f}) = {ps:.3f}, P(pass | {alt:.3f}) = {pa:.3f} — "
          f"{'meets' if meets(ps, pa) else 'FAILS'} the principle")


def show_design(title: str, rows: list[dict], sims: int = 20_000, with_measured: bool = False) -> dict:
    """The thresholds and operating characteristics of S1-S7 on a set of cards (with the joint
    probability at the measured rates too, for pilot.json's cards)."""
    print(f"## {title}")
    rules = criteria_rules(rows, sims)
    n_est = sum(r["establishable"] for r in rows)
    measured = sum(r.get("measured", False) for r in rows)
    cells = rules["cells"]
    print(f"  cards {len(rows)}: valid metric on {sum(r['valid'] for r in rows)}, establishable {n_est}"
          + (f"; rates measured by the pilot on {measured}" if measured else ""))
    for c, rule in rules["S1"].items():
        show_error(f"S1 {c}", rule["n"], E_SUPPORTED)
    for name, rule in cells.items():
        if name.startswith("S2:"):
            show_criterion(f"S2 {name[3:]}: false SUPPORTED where it is an error, the cards outside S1", rule)
    for st in STRATA:
        r3 = rules["S3"][st]
        if not r3["n"]:
            print(f"  S3 {st}: no establishable card — not judged")
            continue
        head = f"S3 {st}: correct definite outcomes, establishable cards (the measured {r3['measured']:.3f}, at most {D_NOMINAL})"
        for t in TIERS:
            v = r3["tiers"][t]
            mark = ("JUDGED" if r3["tier"] == t else
                    "attainable, not judged" if v["attainable"] else "not attainable")
            show_decisive(f"{head}, {t} tier [{mark}]", r3["n"], r3["nominal"], v["alternative"])
    words = {"S4": "false NOT SUPPORTED", "S5": "false metric invalid", "S6": "false NO DETECTABLE EFFECT"}
    for name, rule in cells.items():
        crit, *where = name.split(":")
        if name == "S4:null:opposite":
            show_criterion("S4 null: NOT SUPPORTED against the direction on null data (allowed: its rate bounded "
                           "as a false detection's)", rule)
        elif crit in words:
            cause = {"opposite": " against the direction", "depth": " explained by depth"}.get(where[-1], "")
            show_criterion(f"{crit} {where[0]}: {words[crit]}{cause} where it is an error", rule)
    print(f"  S7a verdicts the engine's deterministic rules cannot give (a label that does not match its cause, an "
          f"unexpected cause, UNIDENTIFIABLE, DEGENERATE METRIC on a varying metric, an effect verdict where the rules "
          f"exclude one: N4 without the replicate unit, the constant metric): n = {rules['S7a']['n']}, none allowed")
    show_criterion(f"S7b crashes (an engine exception or no report; a sound engine has none) at the nominal {E_CRASH:g}",
                   rules["S7b"])
    judged = ", ".join(f"{st} ({rules['S3'][st]['tier']})" for st in STRATA if rules["S3"][st]["tier"])
    print(f"  P(S1-S7 all pass | sound), simulated, S1 at its design, S3 judged on {judged or 'no stratum'}: "
          f"{rules['joint']:.3f}" + ("" if rules["s3_feasible"] else
                                     " — S3 CANNOT BE JUDGED ON EVERY STRATUM: it fails for any validator"))
    if with_measured:
        print(f"  the same with every rate as measured (S1 included): {rules['joint_measured']:.3f}")
    print()
    return rules


def scenario_pilot(blind_low: bool, blind_medium_share: float, establishable_share: float) -> dict:
    """A stand-in pilot for the operating characteristics before the pilot: 8 pairs per level
    on both backgrounds, every high pair valid, the low pairs blind (or valid), a share of the
    medium pairs blind; Δ* above the SESOI at the key dose (E2's and E3's at their depth too); the
    cases of a share of every level's pairs establishable (every second pair for a half, every
    fourth for a quarter)."""
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
                 delta={str(i): {**{f"{f:g}": dict(value=0.125 * f) for f in (0.25, 0.5, 1.0, 1.5)},
                                 f"capture={P.E_CAPTURE:g}": dict(value=0.125)} for i in range(len(pool))},
                 establishable={})
    for c in P.CONDITIONS:
        for v, _ in c.variants:
            for i in range(len(pool)):
                pilot["establishable"][f"{c.name}:{v}:{i}"] = dict(establishable=(i % k) % step == 0)
    return pilot


SCENARIOS = (("Scenario A: the low level blind, half the medium pairs blind, half the cases establishable",
              (True, 0.5, 0.5)),
             ("Scenario B: every pair valid, every case establishable", (False, 0.0, 1.0)),
             ("Scenario C: the low level blind, every medium pair blind, a quarter of the cases establishable",
              (True, 1.0, 0.25)))


def s3_rules_record(rules: dict) -> dict:
    """What pilot.json keeps of S3 before the key (``--write-judged``): every stratum's tier
    (None where it cannot be judged), whether every stratum is judged, and the joint probability."""
    return dict(tiers={st: rules["S3"][st]["tier"] for st in STRATA}, feasible=rules["s3_feasible"],
                joint=rules["joint"], joint_measured=rules["joint_measured"],
                expected={st: dict(n=rules["S3"][st]["n"], nominal=rules["S3"][st]["nominal"]) for st in STRATA})


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--pilot", help="pilot.json: the operating characteristics on its expected cards")
    p.add_argument("--write-judged", action="store_true",
                   help="write S3's tiers (fixed before the key) into the pilot.json given")
    args = p.parse_args(argv)
    print("# operating characteristics of the v1 success criteria (two-sided 95% Clopper-Pearson)")
    print(f"# principle: P(pass | sound) >= {P_PASS_SOUND}, P(pass | doubled nominal error) <= {P_PASS_DOUBLED};")
    print(f"# S1, a conjunction over the key null conditions, meets it as a whole; P(S1-S7 together | sound) >= {P_JOINT_SOUND}")
    print()
    if args.pilot:
        pilot = json.loads(open(args.pilot).read())
        rules = show_design("The expected cards of pilot.json (every pool pair equally often)",
                            expected_rows(pilot, pilot.get("dropped", ())), with_measured=True)
        over = refusal_bound_failures(pilot)
        for st, rec in sorted((pilot.get("gate0_refusals") or {}).items()):
            bound = (f"at most {GATE0_BOUND:g} in the model" + ("; ABOVE THE BOUND: S3 FAILS on this stratum, "
                     "decided before the key (the run goes on)" if st in over else "")) if st in VALID_STRATA \
                else "as measured (a blind or useless metric)"
            print(f"# GATE 0's refusal share before the key, S3 stratum {st}: {rec['refused']} of {rec['n']} "
                  f"({rec['share']:.3f}); {bound}")
        if args.write_judged:
            pilot["s3_rules"] = dict(s3_rules_record(rules), refusal_bound_failed=over)
            with open(args.pilot, "w") as fh:
                fh.write(json.dumps(pilot, indent=1))
            print(f"# S3's tiers written to {args.pilot}: {pilot['s3_rules']['tiers']}")
        return
    m = sum(c.key for c in P.CONDITIONS)
    print("## Nominal rates (errors by cause, decided 2026-10-08; a sound validator's model after the first review,")
    print("## the rule nominals after the second, the cells after the third, GATE 0 and the null's other tail after")
    print("## the fourth)")
    print(f"  Designed sizes: false SUPPORTED {E_SUPPORTED} (alpha/2: the effect test is two-sided at alpha, the")
    print(f"  direction checked afterwards); false NOT SUPPORTED (metric invalid) on a valid metric {E_INVALID} (2 alpha:")
    print("  GATE 4's upper bound below delta_min alpha/2, GATE 5's negative control alpha, a silent positive")
    print(f"  control shown blind alpha/2); false NO DETECTABLE EFFECT {E_NDE} and explained by depth {E_DEPTH} (alpha, the")
    print(f"  TOST's size); NOT SUPPORTED against the direction where it is an error {E_OPPOSITE} (one tail). With pilot.json a card's")
    print("  rates are a sound validator's on its case, in the engine's verdict order (sound_model): GATE 5 at")
    print(f"  {E_GATE5:g}, GATE 4 at the odds measured on the case, then the effect's outcomes as the oracle found them")
    print("  on the case's datasets (after an UNTESTED GATE 4 its explained-by-depth share); a blind, ambiguous or")
    print("  useless metric in the engine's order on them; GATE 0 refuses first with its stratum's share measured")
    print(f"  with the engine before the key (timing.py --refusals), at most {GATE0_BOUND:g} on the strata with a valid")
    print("  metric (a usefulness bound set in advance: a stratum above it fails S3, decided before the key). Without")
    print("  a measurement (the scenarios below): the")
    print(f"  designed sizes of the routes open on the card and decisiveness {D_NOMINAL}. Every wrong outcome is counted")
    print("  in its own cell: S2 a false SUPPORTED per group of the cards outside S1 (N2's steps, N3, N4 with the")
    print("  replicate unit, N6a, N7, the real effects), S4 a false NOT SUPPORTED against the direction or explained")
    print("  by depth, and against the direction on null data (allowed there, bounded as a false detection), S5")
    print("  metric invalid on a valid metric and S6 NO DETECTABLE EFFECT, each per stratum (real effects, nulls,")
    print("  blind or useless metrics); a cell's threshold uses, per card, the larger of the rate and the designed")
    print("  size, and its resolution is the error rate it fails with probability >= 0.95. S3 is judged on every")
    print("  stratum (not N3 and N8), at the principle's tier or the floor's (a validator half as decisive fails);")
    print("  S7a allows no verdict the engine's deterministic rules cannot give (a label that does not match its")
    print("  cause, an unexpected cause, UNIDENTIFIABLE, DEGENERATE METRIC on a varying metric, an effect verdict")
    print("  where the rules give none: N4 without the replicate unit, the constant metric); S7b counts crashes (an")
    print(f"  engine exception or no report, each listed with its traceback) at the nominal {E_CRASH:g} by the principle,")
    print("  a sound engine crashing on none. The joint requirement is checked with S1 at its design.")
    print()
    print(f"## S1 as a whole: the number of datasets per key condition ({m} conditions)")
    first = next(n for n in range(400, 2001) if joint_error_rule(n, E_SUPPORTED, m)[1] >= P_PASS_SOUND)
    fails = [n for n in range(first, 1001) if joint_error_rule(n, E_SUPPORTED, m)[1] < P_PASS_SOUND]
    print(f"  the rule first holds at n = {first}; k* moves in steps, so it fails again at n = "
          f"{_ranges(fails)} (up to 1000); the plan's {P.KEY_N} meets it")
    show_joint("S1 false SUPPORTED, key null conditions", P.KEY_N, E_SUPPORTED, m)
    print("  given the background every dataset is an independent draw, so the counts are binomial and the above is")
    print("  exact; read beyond the background's donors (their shared donors as a random sample: beta-binomial, design")
    print("  effect deff; the run reports the realized deff of every primary criterion, above 1.5 a limitation):")
    for deff in (1.0, 1.1, 1.2, 1.5):
        pall, pd = joint_error_rule_deff(P.KEY_N, E_SUPPORTED, m, deff)
        print(f"    deff {deff:.1f}: P(all {m} pass | {E_SUPPORTED:.3f}) = {pall:.3f}, "
              f"P(a condition passes | {2 * E_SUPPORTED:.3f} there) = {pd:.3f}")
    sound = deff_limit(lambda d: joint_error_rule_deff(P.KEY_N, E_SUPPORTED, m, d)[0] >= P_PASS_SOUND)
    doubled = deff_limit(lambda d: joint_error_rule_deff(P.KEY_N, E_SUPPORTED, m, d)[1] <= P_PASS_DOUBLED)
    print(f"    the principle holds up to deff {min(sound, doubled):.3f} (the doubled side up to {doubled:.3f}, "
          f"the sound side up to {sound:.3f})")
    print()
    budget = P_JOINT_SOUND / joint_error_rule(P.KEY_N, E_SUPPORTED, m)[1]
    print("## S3's tiers on a stratum of n establishable cards")
    for d0 in (0.85, 0.80, 0.70):
        need = {t: next(n for n in range(1, 1000) if all(decisiveness_tiers(k, d0)[t]["p_pass_sound"] >= budget
                                                          for k in range(n, n + 60))) for t in TIERS}
        print(f"  at d = {d0:.2f}: the floor's tier passes a sound validator with >= {budget:.4f} (what S1 at its "
              f"design leaves of {P_JOINT_SOUND}) from n = {need['floor']}, the principle's from n = {need['principle']}")
    print()
    for title, args_ in SCENARIOS:
        show_design(title, expected_rows(scenario_pilot(*args_)))
    print("## Resolution of an error cell by its number of cards (the fourth review): at the nominal alpha/2 and a")
    print("## sound validator at it, the cell's raised nominal and threshold, and the error rate it fails with >= 0.95")
    for n in (100, 200, 300, 500, 790, 1000, 1500, 2000):
        r = error_criterion(n, E_SUPPORTED, E_SUPPORTED)
        print(f"  n={n:4d}: nominal {r['rule_nominal']:.4f}, at most {r['max_allowed']}, resolution "
              f"{detectable_rate(r['max_allowed'], n):.4f}")
    print()
    print("## Precision of a rate from n cards (half-width of the 95% CP interval)")
    for n in (100, 200, 400, 790, 1000, 4000):
        row = ", ".join(f"p={q:.3f}: ±{(clopper_pearson(round(q * n), n)[1] - clopper_pearson(round(q * n), n)[0]) / 2:.3f}"
                        for q in (0.025, 0.1, 0.5))
        print(f"  n={n:4d}: {row}")


if __name__ == "__main__":
    main()
