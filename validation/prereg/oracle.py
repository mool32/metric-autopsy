"""The oracle and the pilot of the confirmatory validation (validation/prereg/v1.md, 3.2 and 8).

The pilot (before the key) fixes, per background and expression level: the SESOI (the smallest
grid value at which the oracle reaches a correct definite outcome on the background's null, N1
on B1 and N7 on B2, in at least 95% of the datasets, counted exactly: one rule, computed on each
background); the response of the metric to the GATE 4 injection over the dose grid (injected
minus sham on whole datasets, as GATE 4 measures it) and the saturation dose, the largest grid
dose without saturation (the smallest dose whose response reaches 95% of the grid maximum); and
per pool pair the truth about the metric: its population response at that dose against
delta_min = 0.5 x SESOI, valid at >= 1.2 delta_min, blind at <= 0.8 delta_min, ambiguous
between - on the background's null, and for the conditions whose data change the pair's counts
(N2-N5, N8, E1-E3: panel.TRUTH_CASE_CONDITIONS) on the condition's own datasets ("truth_case").
On B1 it also fixes, per level, the key dose (the smallest grid dose at which the oracle detects
E1 with power >= 0.9 and the level's Δ* >= 1.25 x SESOI; where there is none, the real effects
use the saturation dose), Δ* per pair (the population difference of the metric under the
injected coupling, signal side minus sham, from a large paired simulation; value and standard
error), and the oracle's power to detect N2's raw difference at c = 0.5. Finally the
establishable (condition, variant, pair) cases - those where the oracle reaches a correct
definite outcome with power >= 0.9 - with what a sound validator's model needs (oc.sound_model):
the oracle's outcomes, the effect's outcomes with the metric's validity known, the outcomes of
the engine's rules where GATE 4 decides, and the odds of GATE 4 on the case's response.

The oracle is told what the engine has to find out. For a valid (or ambiguous) metric it takes
the metric's validity as known, removes a planted nuisance by applying it to the other side as
well (capture loss with the true factor, dropout with the true fraction, N8's per-cell capture
with fresh Beta(2, 2) draws), and tests the per-replicate metric values with the test the graded
replicate rule prescribes: an exact (or Monte Carlo) permutation over replicates at >= 4 per
group, Welch's t at 3; paired designs by sign flips; equivalence by TOST, the (1 - 2 alpha)
interval inside ±SESOI. Its outcome follows the engine's order (explained by depth - with a
SESOI only where the corrected effect is shown smaller than it, as the engine since journal D7 -,
reversed sign, detected in or against the declared direction, equivalent, else INCONCLUSIVE). For a blind
or useless metric the correct definite outcome is "metric invalid", which a validator can only
reach by showing blindness: the oracle runs GATE 4's rule on the dataset (GATE4_REPS injections
against shams, the upper 95% bound below delta_min). The oracle's power is the most a validator
could reach, and only it decides which cases are establishable: the engine never does.

Independence: numpy and the t distribution of scipy.stats only; this file never imports
``metric_autopsy``. The pilot runs on datasets drawn with seeds outside the panel's (they come
from the public PILOT_SEED, the panel's from the beacon's key).

    python validation/prereg/oracle.py --backgrounds backgrounds.json --data-dir DATA --out pilot.json --workers 4
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import frozen  # standard library only

frozen.pin_numerics()  # one numerical path, as run_panel.py: set before numpy loads

import numpy as np  # noqa: E402

import panel as P  # noqa: E402

ALPHA = 0.05
POWER = 0.90          # establishable: oracle power >= 0.90
SESOI_TARGET = 0.95   # the SESOI: the oracle reaches a correct definite outcome on the null this often
PILOT_SEED = 20261008
PILOT_DATASETS = 100  # per establishability case (condition, variant, pair) and per level for the SESOI
DELTA_DRAWS = 2000    # paired donor draws per pair and dose for Δ*
DELTA_MARGIN = 1.25   # the key dose needs the level's Δ* >= 1.25 x SESOI
SESOI_GRID = tuple(round(0.005 * k, 3) for k in range(1, 201))
# The coupling dose enters a logistic keep probability; far above 4 the thinning becomes almost
# all-or-none, the co-detected cells are those that kept everything, and the coupling a
# co-detection metric sees fades again. The grid stays in the range where the response rises.
DOSE_GRID = tuple(round(0.25 * k, 2) for k in range(1, 17))
SATURATION_SHARE = 0.95  # the saturation dose: the smallest grid dose reaching 95% of the maximum
CURVE_DATASETS, CURVE_REPS = 8, 4    # per pair and dose: the level's response curve
TRUTH_DATASETS, TRUTH_REPS = 60, 20  # per pair at the saturation dose: its population response
CASE_DATASETS, CASE_REPS = 40, 10    # per pair on each condition and variant that changes its counts
GATE4_REPS = 200      # injections per dataset in GATE 4's rule (the engine's gates.GATE4_N_REP)
# The sizes of the pilot's simulations; the dry run's smoke test shrinks them (--smoke), the real
# pilot uses these.
SIZES = dict(curve=(CURVE_DATASETS, CURVE_REPS), truth=(TRUTH_DATASETS, TRUTH_REPS),
             case=(CASE_DATASETS, CASE_REPS), gate4=GATE4_REPS, sesoi=PILOT_DATASETS)
SMOKE_SIZES = dict(curve=(2, 2), truth=(4, 4), case=(3, 3), gate4=20, sesoi=8)
SMOKE_DATASETS, SMOKE_DRAWS = 3, 40  # the dry run's smoke test (--smoke): the same code, tiny sizes
E1_FACTORS = (0.25, 0.5, 1.0, 1.5)
MODULE_FOLD, MODULE_FRAC = 2.0, 0.3  # the engine's GATE 4 module injection in the claim cards


# --------------------------------------------------------------------------- #
# the metrics, re-implemented from their definitions
# --------------------------------------------------------------------------- #
def _npc(a: np.ndarray, b: np.ndarray, tot: np.ndarray) -> float:
    """CP10k + log1p Pearson of two genes over the cells where both are detected, with the
    cells' totals given: 0 with fewer than 3 such cells or a constant gene (as the engine)."""
    both = (a > 0) & (b > 0)
    if both.sum() < 3:
        return 0.0
    tot = np.where(tot > 0, tot, 1.0)
    la, lb = np.log1p(a / tot * 1e4)[both], np.log1p(b / tot * 1e4)[both]
    if la.std() == 0 or lb.std() == 0:
        return 0.0
    return float(np.corrcoef(la, lb)[0, 1])


def norm_pearson(X: np.ndarray, ia: int, ib: int) -> float:
    """CP10k + log1p Pearson of genes a and b over the cells where both are detected."""
    return _npc(X[:, ia], X[:, ib], X.sum(axis=1))


def module_score(X: np.ndarray, cols: list) -> float:
    tot = X.sum(axis=1, keepdims=True)
    return float(np.log1p(X[:, cols] / np.where(tot > 0, tot, 1.0) * 1e4).mean())


# --------------------------------------------------------------------------- #
# GATE 4's rule, re-implemented: the response to an injected coupling against its sham
# --------------------------------------------------------------------------- #
def coupling_deltas(X: np.ndarray, ia: int, ib: int, dose: float, rng, reps: int) -> np.ndarray:
    """`reps` responses of log-normalized Pearson on the whole dataset: the pair coupled by one
    shared keep probability 1 / (1 + e^(-dose z)) minus the pair thinned with independent ones
    (``panel.inject_coupling`` / ``sham_coupling`` on the two columns; the other genes' totals
    are unchanged)."""
    a0, b0 = X[:, ia], X[:, ib]
    rest = X.sum(axis=1) - a0 - b0
    n = len(a0)
    out = np.empty(reps)
    for r in range(reps):
        p = P._keep(rng.standard_normal(n), dose)
        a, b = P.thin(a0, p, rng), P.thin(b0, p, rng)
        sa, sb = (P.thin(a0, P._keep(rng.standard_normal(n), dose), rng),
                  P.thin(b0, P._keep(rng.standard_normal(n), dose), rng))
        out[r] = _npc(a, b, rest + a + b) - _npc(sa, sb, rest + sa + sb)
    return out


def module_deltas(X: np.ndarray, score_cols: list, module_cols: list, rng, reps: int) -> np.ndarray:
    """`reps` responses of the random-gene score to the G2M module injection (every gene but the
    module thinned to 1/fold in a random 30% of cells) against its sham (every gene thinned in a
    random 30% of cells). The other genes enter only through each cell's total; the thinned total
    of independent binomial thinnings with one probability is the binomial of the total, so the
    simulation is exact in distribution."""
    S = X[:, score_cols]
    mod = X[:, module_cols].sum(axis=1)
    rest = X.sum(axis=1) - S.sum(axis=1) - mod
    n = X.shape[0]
    q = 1.0 / MODULE_FOLD

    def score(Sx, tot):
        return float(np.log1p(Sx / np.where(tot > 0, tot, 1.0)[:, None] * 1e4).mean())

    out = np.empty(reps)
    for r in range(reps):
        c = rng.random(n) < MODULE_FRAC
        Si, resti = S.copy(), rest.copy()
        Si[c] = P.thin(S[c], q, rng)
        resti[c] = P.thin(rest[c], q, rng)
        inj = score(Si, Si.sum(axis=1) + mod + resti)
        c = rng.random(n) < MODULE_FRAC
        Ss, mods, rests = S.copy(), mod.copy(), rest.copy()
        Ss[c], mods[c], rests[c] = P.thin(S[c], q, rng), P.thin(mod[c], q, rng), P.thin(rest[c], q, rng)
        out[r] = inj - score(Ss, Ss.sum(axis=1) + mods + rests)
    return out


def random_deltas(rng, reps: int) -> np.ndarray:
    """A metric that returns a random number: each response is the difference of two draws."""
    return rng.standard_normal(reps) - rng.standard_normal(reps)


def shows_blind(deltas: np.ndarray, delta_min: float, alpha: float = ALPHA) -> bool:
    """GATE 4's FAIL: the upper bound of the two-sided (1 - alpha) t interval of the mean
    response, in the declared direction (an increase), is below delta_min."""
    return gate4_outcome(deltas, delta_min, alpha) == "FAIL"


def gate4_outcome(deltas: np.ndarray, delta_min: float, alpha: float = ALPHA) -> str:
    """GATE 4's rule on the responses: FAIL if the upper bound of the two-sided (1 - alpha) t
    interval is below delta_min (checked first), PASS if its lower bound is above 0, else UNTESTED."""
    d = np.asarray(deltas, float)
    half = t_ppf(1 - alpha / 2, len(d) - 1) * d.std(ddof=1) / math.sqrt(len(d))
    m = float(d.mean())
    if m + half < delta_min:
        return "FAIL"
    return "PASS" if m - half > 0 else "UNTESTED"


def gate4_odds(mean: float, between_sd: float, within_sd: float, delta_min: float,
               reps: int = GATE4_REPS, alpha: float = ALPHA) -> dict:
    """The probabilities of GATE 4's outcomes on a dataset of this case, from the case's response
    statistics: a dataset's mean response varies around `mean` with `between_sd`, and GATE 4's
    mean of `reps` injections around it with within_sd / sqrt(reps); its interval has the
    half-width t x within_sd / sqrt(reps) (a normal approximation of the rule)."""
    from scipy import stats
    se = within_sd / math.sqrt(reps)
    half = t_ppf(1 - alpha / 2, reps - 1) * se
    sd = math.sqrt(max(between_sd, 0.0) ** 2 + se ** 2)
    if not sd > 0:
        fail = float(mean + half < delta_min)
        passed = float(not fail and mean - half > 0)
    else:
        fail = float(stats.norm.cdf((delta_min - half - mean) / sd))
        passed = float(stats.norm.sf((max(half, delta_min - half) - mean) / sd))
    return dict(p_pass=passed, p_fail=fail, p_untested=max(0.0, 1.0 - passed - fail))


# --------------------------------------------------------------------------- #
# tests on replicate values
# --------------------------------------------------------------------------- #
def t_ppf(q: float, df: float) -> float:
    from scipy import stats  # the oracle's only dependency outside numpy
    return float(stats.t.ppf(q, df))


def _t_cdf(t: float, df: float) -> float:
    from scipy import stats
    return float(stats.t.cdf(t, df))


def welch(a: np.ndarray, b: np.ndarray):
    """Difference mean(a) - mean(b), its standard error and Welch's degrees of freedom."""
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    se = math.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1)) if va + vb > 0 else 1.0
    return float(a.mean() - b.mean()), se, df


_COMBOS: dict = {}


def _combos(n: int, k: int) -> np.ndarray:
    if (n, k) not in _COMBOS:
        _COMBOS[(n, k)] = np.array(list(itertools.combinations(range(n), k)), dtype=np.int64)
    return _COMBOS[(n, k)]


def permutation_p(a: np.ndarray, b: np.ndarray, max_exact: int = 20000, rng=None) -> float:
    """Two-sided permutation p of the difference in means: exact when the number of
    assignments is at most max_exact, else 4999 Monte Carlo draws, (1 + hits) / (1 + draws)."""
    y = np.concatenate([a, b]).astype(float)
    n, k = len(y), len(a)
    obs = abs(a.mean() - b.mean())
    total = math.comb(n, k)
    eps = 1e-12 * max(1.0, obs)
    if total <= max_exact:
        sa = y[_combos(n, k)].sum(axis=1)
        diff = np.abs(sa / k - (y.sum() - sa) / (n - k))
        return float(np.mean(diff >= obs - eps))
    rng = rng if rng is not None else np.random.default_rng(0)
    draws = 4999
    hits = 0
    for _ in range(draws):
        p = rng.permutation(n)
        hits += abs(y[p[:k]].mean() - y[p[k:]].mean()) >= obs - eps
    return (1 + hits) / (1 + draws)


def signflip_p(d: np.ndarray) -> float:
    """Two-sided exact sign-flip p of the mean paired difference."""
    obs = abs(d.mean())
    signs = np.array(list(itertools.product((1.0, -1.0), repeat=len(d))))
    return float(np.mean(np.abs(signs @ d / len(d)) >= obs - 1e-12 * max(1.0, obs)))


def detects(a: np.ndarray, b: np.ndarray, alpha: float = ALPHA) -> tuple[bool, float]:
    """The graded rule's two-sided test of a difference: permutation at >= 4 per group, Welch at 3."""
    est = float(a.mean() - b.mean())
    if min(len(a), len(b)) >= 4:
        return permutation_p(a, b) < alpha, est
    est, se, df = welch(a, b)
    return (se > 0 and 2 * (1 - _t_cdf(abs(est) / se, df)) < alpha), est


def detects_paired(d: np.ndarray, alpha: float = ALPHA) -> tuple[bool, float]:
    return signflip_p(d) < alpha, float(d.mean())


def tost_width(a: np.ndarray, b: np.ndarray, alpha: float = ALPHA) -> float:
    """max |bound| of the (1 - 2 alpha) Welch interval: equivalent within ±SESOI iff < SESOI."""
    est, se, df = welch(a, b)
    half = t_ppf(1 - alpha, df) * se
    return max(abs(est - half), abs(est + half))


def paired_tost_width(d: np.ndarray, alpha: float = ALPHA) -> float:
    est, se, df = float(d.mean()), float(d.std(ddof=1) / math.sqrt(len(d))), len(d) - 1
    half = t_ppf(1 - alpha, df) * se
    return max(abs(est - half), abs(est + half))


# --------------------------------------------------------------------------- #
# per-replicate values of one dataset
# --------------------------------------------------------------------------- #
def per_donor(X, obs, fn) -> tuple[np.ndarray, np.ndarray]:
    """Metric of every donor's cells, split by group (unpaired designs)."""
    vals = {"A": [], "B": []}
    donor, group = np.asarray(obs["donor"]), np.asarray(obs["group"])
    for d in dict.fromkeys(donor):
        m = donor == d
        vals[str(group[m][0])].append(fn(X[m]))
    return np.asarray(vals["A"]), np.asarray(vals["B"])


def per_donor_paired(X, obs, fn) -> np.ndarray:
    out = []
    donor, group = np.asarray(obs["donor"]), np.asarray(obs["group"])
    for d in dict.fromkeys(donor):
        m = donor == d
        out.append(fn(X[m & (group == "A")]) - fn(X[m & (group == "B")]))
    return np.asarray(out)


EQUALIZE_DRAWS = 5  # the engine's equalizations per analysis (effect.estimate_effect's n_equalize)


def equalized_sets(X, obs, rng, draws: int = EQUALIZE_DRAWS) -> list:
    """The engine's composition correction where nothing is planted (equalize.thin_to_match, as
    effect.estimate_effect runs it; the third review found the oracle analysing such data raw): the
    group with the higher mean total thinned by one common ratio to the other's, `draws` times."""
    group = np.asarray(obs["group"])
    tot = X.sum(axis=1)
    a, b = group == "A", group == "B"
    hi, lo = (a, b) if tot[a].mean() >= tot[b].mean() else (b, a)
    p = float(np.clip(tot[lo].mean() / tot[hi].mean(), 0.0, 1.0)) if tot[hi].mean() > 0 else 1.0
    out = []
    for _ in range(draws):
        Y = X.copy()
        Y[hi] = rng.binomial(np.round(Y[hi]).astype(np.int64), p).astype(np.float64)
        out.append(Y)
    return out


def corrected_values(X, obs, entry: dict, fn, rng, paired: bool = False):
    """The replicate values the oracle tests: after the true correction where a nuisance is
    planted (`true_correction`), else after the engine's equalization (`equalized_sets`, the mean
    of the replicate values over its draws, as the engine averages them)."""
    per = per_donor_paired if paired else per_donor
    fixed = true_correction(X, obs, entry, rng)
    if fixed is not None:
        return per(fixed, obs, fn)
    vals = [per(Y, obs, fn) for Y in equalized_sets(X, obs, rng)]
    if paired:
        return np.mean(vals, axis=0)
    return np.mean([v[0] for v in vals], axis=0), np.mean([v[1] for v in vals], axis=0)


def true_correction(X, obs, entry: dict, rng):
    """The planted nuisance applied to the side that lacks it (None if nothing was planted)."""
    name, variant, side = entry["condition"], entry["variant"], entry["side"]
    group = np.asarray(obs["group"])
    other = "B" if side == "A" else "A"
    if name == "N2":
        target, fn = group != side, (lambda x: P.thin(x, float(variant.split("=")[1]), rng))
    elif name == "N3":
        target, fn = group != side, (lambda x: P.drop_out(x, float(variant.split("=")[1]), rng))
    elif name == "N8":
        target, fn = group != side, (lambda x: P.variable_capture(x, rng))
    elif name == "E2":   # capture loss on the signal side: thin the other side too
        target, fn = group == other, (lambda x: P.thin(x, 0.5, rng))
    elif name == "E3":   # capture loss on the other side: thin the signal side too
        target, fn = group == side, (lambda x: P.thin(x, 0.5, rng))
    else:
        return None
    X = X.copy()
    X[target] = fn(X[target])
    return X


# --------------------------------------------------------------------------- #
# the oracle's outcome on one dataset of known truth
# --------------------------------------------------------------------------- #
def outcome_from_values(raw, corr, corrected: bool, sesoi: float, direction: str) -> str:
    """The engine's decision order on replicate values, for an oracle that knows the metric is
    valid: a detected raw difference that the true correction removes (< half retained) is
    explained by depth if the corrected effect is equivalent to zero within the SESOI (TOST; the
    engine's rule with a SESOI, journal D7); one whose sign the correction reverses is
    INCONCLUSIVE; a detected effect is SUPPORTED in the declared direction and NOT SUPPORTED
    (opposite direction) against it; otherwise NO DETECTABLE EFFECT if the TOST establishes
    equivalence, else INCONCLUSIVE. `raw` and `corr` are (A values, B values); the effect is
    A - B, so > 0 is a decrease."""
    hit, est = detects(*corr)
    if corrected:
        raw_hit, raw_est = detects(*raw)
        retained = est / raw_est if raw_est else float("nan")
        if (raw_hit and not hit and np.isfinite(retained) and abs(retained) < 0.5
                and tost_width(*corr) < sesoi):
            return P.NS_DEPTH
        if raw_hit and hit and np.sign(est) != np.sign(raw_est):
            return P.INCONCLUSIVE
    if hit:
        observed = "decrease" if est > 0 else "increase"
        return P.SUPPORTED if observed == direction or direction == "two-sided" else P.NS_OPPOSITE
    return P.NDE if tost_width(*corr) < sesoi else P.INCONCLUSIVE


def paired_outcome(raw, corr, sesoi: float, direction: str) -> str:
    """`outcome_from_values` for the paired design (N5): per-donor differences A - B, raw and
    corrected, tested by sign flips and the paired TOST, in the engine's order."""
    hit, est = detects_paired(corr)
    raw_hit, raw_est = detects_paired(raw)
    retained = est / raw_est if raw_est else float("nan")
    if raw_hit and not hit and np.isfinite(retained) and abs(retained) < 0.5 and paired_tost_width(corr) < sesoi:
        return P.NS_DEPTH
    if raw_hit and hit and np.sign(est) != np.sign(raw_est):
        return P.INCONCLUSIVE
    if hit:
        observed = "decrease" if est > 0 else "increase"
        return P.SUPPORTED if observed == direction or direction == "two-sided" else P.NS_OPPOSITE
    return P.NDE if paired_tost_width(corr) < sesoi else P.INCONCLUSIVE


def effect_outcome(entry: dict, X, obs, card: dict, ia: int, ib: int, rng) -> str:
    """The oracle's outcome on the effect, with the metric's validity known: the engine's analysis
    (its correction where nothing is planted, the true one where a nuisance is), in its order."""
    sesoi = float(card["prereg"]["sesoi"])
    direction = card["prereg"]["direction"]
    fn = (lambda x: norm_pearson(x, ia, ib))
    if entry["condition"] == "N5":
        return paired_outcome(per_donor_paired(X, obs, fn), corrected_values(X, obs, entry, fn, rng, paired=True),
                              sesoi, direction)
    return outcome_from_values(per_donor(X, obs, fn), corrected_values(X, obs, entry, fn, rng), True, sesoi,
                               direction)


def _values_outcome(raw, card: dict, corr=None) -> str:
    """The engine's order on replicate values of null data without a planted nuisance: raw, and
    after the engine's equalization where it is given (`corr`)."""
    return outcome_from_values(raw, raw if corr is None else corr, corr is not None,
                               float(card["prereg"]["sesoi"]), card["prereg"]["direction"])


def oracle_outcomes(entry: dict, X, obs, card: dict, bgs: dict, pilot: dict, rng) -> dict:
    """On one dataset of known truth: ``best``, the outcome the oracle reaches (it knows the
    truth: for a valid metric its validity, for a blind or useless one only GATE 4's rule can show
    it, for an ambiguous one either); ``effect``, the outcome of the effect's analysis taken as
    valid (None where there is none); ``gate4``, GATE 4's outcome on the dataset where the oracle
    ran its rule (blind, ambiguous and useless metrics; None for a valid one, whose odds come from
    the case's response statistics). A sound validator that follows the engine's rules reaches
    what `engine_outcome` gives (oc.sound_model)."""
    cond = P.conditions()[entry["condition"]]
    bg = bgs[cond.background]
    truth = P.metric_truth(cond, entry["variant"], entry["pair"], pilot)
    dmin = float(card["prereg"]["delta_min"])
    reps = int(_sizes()["gate4"])
    if truth == "constant":
        return dict(best=P.DEGENERATE, effect=None, gate4=None)
    if truth == "useless":
        if cond.metric == "random":
            g4 = gate4_outcome(random_deltas(rng, reps), dmin)
            values = (rng.standard_normal(P.DONORS_PER_GROUP), rng.standard_normal(P.DONORS_PER_GROUP))
        else:
            genes = list(bg.genes)
            cols = [genes.index(g) for g in card["score_genes"]]
            g4 = gate4_outcome(module_deltas(X, cols, [genes.index(g) for g in bg.plan["g2m"]], rng, reps), dmin)
            score = (lambda x: module_score(x, cols))
            return dict(best=P.NS_INVALID if g4 == "FAIL" else P.INCONCLUSIVE, gate4=g4,
                        effect=_values_outcome(per_donor(X, obs, score), card,
                                               corrected_values(X, obs, entry, score, rng)))
        return dict(best=P.NS_INVALID if g4 == "FAIL" else P.INCONCLUSIVE, effect=_values_outcome(values, card),
                    gate4=g4)
    pe = P.pool_entry(bg, entry)
    genes = list(bg.genes)
    ia, ib = genes.index(pe["pair"][0]), genes.index(pe["pair"][1])
    effect = effect_outcome(entry, X, obs, card, ia, ib, rng) if cond.oracle else None
    if truth == "valid":
        return dict(best=effect if effect is not None else P.INCONCLUSIVE, effect=effect, gate4=None)
    g4 = gate4_outcome(coupling_deltas(X, ia, ib, float(card["signal_test"]["strength"]), rng, reps), dmin)
    if truth == "blind" or effect is None:
        best = P.NS_INVALID if g4 == "FAIL" else P.INCONCLUSIVE
    else:  # ambiguous: either verdict is correct
        best = effect if effect != P.INCONCLUSIVE else (P.NS_INVALID if g4 == "FAIL" else P.INCONCLUSIVE)
    return dict(best=best, effect=effect, gate4=g4)


def engine_outcome(effect: str | None, gate4: str | None) -> str:
    """What a sound validator gets by the engine's verdict order (report.decide_cause) where GATE 4's
    rule ran: "metric invalid" where it FAILs; else an effect explained by depth, which the engine
    checks before the metric's validity (the second review); else the effect's outcome where GATE 4
    PASSes and INCONCLUSIVE where it is UNTESTED (an untested metric's absence or detection is
    inconclusive)."""
    if gate4 == "FAIL":
        return P.NS_INVALID
    if effect == P.NS_DEPTH:
        return P.NS_DEPTH
    if gate4 == "PASS":
        return effect or P.INCONCLUSIVE
    return P.INCONCLUSIVE


def oracle_outcome(entry: dict, X, obs, card: dict, bgs: dict, pilot: dict, rng) -> str:
    """The outcome the oracle reaches on this dataset (it knows the truth)."""
    return oracle_outcomes(entry, X, obs, card, bgs, pilot, rng)["best"]


# --------------------------------------------------------------------------- #
# Δ*: the population difference of the metric under the injected coupling
# --------------------------------------------------------------------------- #
def delta_star(bg: P.Background, pair_index: int, dose: float, draws: int = DELTA_DRAWS,
               seed: int = PILOT_SEED) -> dict:
    """E[metric | injected coupling] - E[metric | sham] per donor, on CELLS_PER_DONOR cells of a
    random donor: each draw takes one donor's cells and applies both to the same cells (paired),
    so the mean difference estimates Δ* with a small standard error."""
    pe = bg.plan["pool"][pair_index]
    genes = list(bg.genes)
    ia, ib = genes.index(pe["pair"][0]), genes.index(pe["pair"][1])
    rng = np.random.default_rng([seed, pair_index, int(round(dose * 1000))])
    donors = bg.donors
    diffs = np.empty(draws)
    for t in range(draws):
        d = donors[int(rng.integers(len(donors)))]
        idx = np.where(bg.donor == d)[0]
        rows = rng.choice(idx, size=min(P.CELLS_PER_DONOR, len(idx)), replace=False)
        Xd = np.asarray(bg.X[rows], dtype=np.float64)
        diffs[t] = coupling_deltas(Xd, ia, ib, dose, rng, 1)[0]
    return dict(value=float(diffs.mean()), se=float(diffs.std(ddof=1) / math.sqrt(draws)), draws=draws,
                dose=float(dose))


def delta_level(bg: P.Background, level: str, dose: float, draws: int = DELTA_DRAWS) -> dict:
    """The level's Δ*: the mean over its pairs (the key draws a level's pairs uniformly)."""
    idx = [k for k, pe in enumerate(bg.plan["pool"]) if pe["level"] == level]
    per = [delta_star(bg, k, dose, max(2, draws // len(idx))) for k in idx]
    return dict(value=float(np.mean([p["value"] for p in per])),
                se=float(math.sqrt(sum(p["se"] ** 2 for p in per)) / len(per)), pairs=idx, dose=float(dose),
                per_pair={str(k): p["value"] for k, p in zip(idx, per)})


# --------------------------------------------------------------------------- #
# pilot datasets (public seeds) and their parallel map
# --------------------------------------------------------------------------- #
def pilot_entries(cond: P.Condition, variant: str, pair: int, n: int, salt: int = 0) -> list[dict]:
    """Public-seed datasets of one condition, variant and pool pair."""
    vi = [v for v, _ in cond.variants].index(variant)
    ss = np.random.SeedSequence([PILOT_SEED, salt, P.CONDITIONS.index(cond), vi, int(pair)])
    out = []
    for i, child in enumerate(ss.spawn(n)):
        rng = np.random.default_rng(child)
        out.append(dict(id=f"P{i}", condition=cond.name, variant=variant, index=i,
                        side=("A", "B")[int(rng.integers(2))], pair=int(pair),
                        seed=int(child.generate_state(1, np.uint64)[0])))
    return out


def level_pairs(bg: P.Background, level: str) -> list[int]:
    return [k for k, pe in enumerate(bg.plan["pool"]) if pe["level"] == level]


_STATE: dict = {}  # backgrounds, inherited by forked workers


def _sizes() -> dict:
    return _STATE.get("sizes") or SIZES


def _map(fn, jobs: list, workers: int) -> list:
    """`fn` over `jobs` in forked processes, in order; a worker the system kills (out of memory)
    breaks the pool at once rather than at the job's timeout (as run_panel.fork_map)."""
    if workers <= 1 or len(jobs) < 2:
        return [fn(j) for j in jobs]
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(workers, mp_context=mp.get_context("fork")) as ex:
        return list(ex.map(fn, jobs, chunksize=max(1, len(jobs) // (8 * workers))))


def _base_condition(background: str) -> P.Condition:
    """The null design whose datasets define a background's SESOI and truth: N1 on B1, N7 on B2."""
    return P.conditions()["N1" if background == "B1" else "N7"]


def _stand_in(pilot: dict, background: str, level: str, **over) -> dict:
    """A pilot with the fields `panel.build` reads filled in for one background and level."""
    out = json.loads(json.dumps(pilot))
    for key, val in dict(sesoi=1.0, saturation_dose=1.0).items():
        out.setdefault(key, {}).setdefault(background, {}).setdefault(level, val)
    out.setdefault("e_dose", {}).setdefault(level, 1.0)
    for k, v in over.items():
        out[k] = v
    return out


def _response_job(job):
    """Injected-minus-sham responses of one pair on one public-seed dataset of a condition and
    variant (default: the background's null): their mean and variance."""
    background, pair, dose, i, reps, salt, name, variant = job
    bgs, pilot = _STATE["bgs"], _STATE["pilot"]
    bg = bgs[background]
    cond = P.conditions()[name] if name else _base_condition(background)
    variant = variant or cond.variants[0][0]
    e = pilot_entries(cond, variant, pair, i + 1, salt)[i]
    X, _, genes, _ = P.build(e, bgs, _stand_in(pilot, background, bg.plan["pool"][pair]["level"]))
    pe = bg.plan["pool"][pair]
    rng = np.random.default_rng([PILOT_SEED, salt, pair, i, int(round(dose * 100))])
    d = coupling_deltas(X, genes.index(pe["pair"][0]), genes.index(pe["pair"][1]), dose, rng, reps)
    return float(d.mean()), float(d.var(ddof=1)) if len(d) > 1 else 0.0


def response_curve(bgs: dict, background: str, level: str, pilot: dict, workers: int = 1) -> dict:
    """The level's response over the dose grid: per dose the mean over its pairs of the
    injected-minus-sham response on CURVE_DATASETS null datasets x CURVE_REPS injections."""
    _STATE.update(bgs=bgs, pilot=pilot)
    n_ds, reps = _sizes()["curve"]
    pairs = level_pairs(bgs[background], level)
    jobs = [(background, k, d, i, reps, 1, None, None) for d in DOSE_GRID for k in pairs for i in range(n_ds)]
    vals = _map(_response_job, jobs, workers)
    by = {}
    for (_, _, d, *_), (v, _) in zip(jobs, vals):
        by.setdefault(d, []).append(v)
    return {f"{d:g}": dict(response=float(np.mean(v)), se=float(np.std(v, ddof=1) / math.sqrt(len(v))))
            for d, v in by.items()}


def saturation_dose(curve: dict) -> float:
    """The largest grid dose without saturation: the smallest dose whose response reaches 95% of
    the grid's maximum (the largest grid dose where the maximum is not above 0)."""
    resp = {float(d): v["response"] for d, v in curve.items()}
    top = max(resp.values())
    if top <= 0:
        return max(resp)
    return min(d for d, r in resp.items() if r >= SATURATION_SHARE * top)


def response_record(means: list, variances: list, reps: int, dose: float, delta_min: float) -> dict:
    """A pair's response on a case's datasets: the population mean response and its standard
    error, the between-dataset SD (the dataset means' variance less the injections' share), the
    within-dataset SD of one injection, the class against delta_min, and GATE 4's odds."""
    m, v = np.asarray(means, float), np.asarray(variances, float)
    within = float(math.sqrt(max(v.mean(), 0.0)))
    between = float(math.sqrt(max(m.var(ddof=1) - within ** 2 / reps, 0.0))) if len(m) > 1 else 0.0
    r = float(m.mean())
    return {"response": r, "se": float(m.std(ddof=1) / math.sqrt(len(m))) if len(m) > 1 else float("nan"),
            "between_dataset_sd": between, "within_sd": within, "datasets": len(m), "reps": reps,
            "dose": float(dose), "delta_min": float(delta_min), "class": P.classify_response(r, delta_min),
            "gate4": gate4_odds(r, between, within, delta_min, int(_sizes()["gate4"]))}


def pair_truth(bgs: dict, background: str, pilot: dict, workers: int = 1) -> dict:
    """Per pool pair: its population response at its level's saturation dose on the background's
    null (TRUTH_DATASETS datasets x TRUTH_REPS injections), and the class it gives against
    delta_min (`response_record`)."""
    _STATE.update(bgs=bgs, pilot=pilot)
    n_ds, reps = _sizes()["truth"]
    bg = bgs[background]
    jobs = [(background, k, float(pilot["saturation_dose"][background][pe["level"]]), i, reps, 2, None, None)
            for k, pe in enumerate(bg.plan["pool"]) for i in range(n_ds)]
    vals = _map(_response_job, jobs, workers)
    out = {}
    for k, pe in enumerate(bg.plan["pool"]):
        got = [x for job, x in zip(jobs, vals) if job[1] == k]
        dmin = P.DELTA_MIN_FRACTION * float(pilot["sesoi"][background][pe["level"]])
        out[str(k)] = response_record([g[0] for g in got], [g[1] for g in got], reps, jobs[k * n_ds][2], dmin)
    return out


def case_truth(bgs: dict, pilot: dict, workers: int = 1) -> dict:
    """The truth about the metric on the data of every condition and variant that changes the
    pair's counts or the dataset's size (panel.TRUTH_CASE_CONDITIONS: capture loss, dropout, N4's and N5's designs, variable capture, the real
    effects' coupling): per pool pair its population response to GATE 4's injection on
    CASE_DATASETS datasets of the case x CASE_REPS injections, at the card's dose (decided after
    the first review: the truth measured on N1 alone does not hold where the data thin the pair)."""
    _STATE.update(bgs=bgs, pilot=pilot)
    n_ds, reps = _sizes()["case"]
    out, jobs = {}, []
    for cond in P.CONDITIONS:
        if cond.name not in P.TRUTH_CASE_CONDITIONS:
            continue
        bg = bgs[cond.background]
        for variant, _ in cond.variants:
            for k, pe in enumerate(bg.plan["pool"]):
                dose = float(pilot["saturation_dose"][cond.background][pe["level"]])
                jobs += [(cond.background, k, dose, i, reps, 7, cond.name, variant) for i in range(n_ds)]
    vals = _map(_response_job, jobs, workers)
    by = {}
    for job, v in zip(jobs, vals):
        by.setdefault((job[0], job[6], job[7], job[1]), []).append((job[2], v))
    for (background, name, variant, k), got in by.items():
        level = bgs[background].plan["pool"][k]["level"]
        dmin = P.DELTA_MIN_FRACTION * float(pilot["sesoi"][background][level])
        out.setdefault(background, {}).setdefault(f"{name}:{variant}", {})[str(k)] = response_record(
            [g[1][0] for g in got], [g[1][1] for g in got], reps, got[0][0], dmin)
    return out


def _null_job(job):
    """On one public-seed null dataset of the background (N1 / N7), for the SESOI: whether the
    effect test detects a difference, whether against the card's declared direction, and the
    TOST's width."""
    background, pair, i = job
    bgs, pilot = _STATE["bgs"], _STATE["pilot"]
    bg = bgs[background]
    cond = _base_condition(background)
    e = pilot_entries(cond, cond.variants[0][0], pair, i + 1, 3)[i]
    X, obs, genes, cards = P.build(e, bgs, _stand_in(pilot, background, bg.plan["pool"][pair]["level"]))
    pe = bg.plan["pool"][pair]
    ia, ib = genes.index(pe["pair"][0]), genes.index(pe["pair"][1])
    # the engine's analysis of the null: after its equalization of depth (the third review)
    a, b = corrected_values(X, obs, e, lambda x: norm_pearson(x, ia, ib), np.random.default_rng([e["seed"], 2]))
    hit, est = detects(a, b)
    against = hit and ("decrease" if est > 0 else "increase") != cards[-1]["prereg"]["direction"]
    return bool(hit), bool(against), float(tost_width(a, b))


def choose_sesoi(bgs: dict, background: str, level: str, pilot: dict, n: int = PILOT_DATASETS,
                 workers: int = 1) -> tuple[float, bool]:
    """The smallest SESOI on the grid at which the oracle reaches a correct definite outcome on the
    background's null (N1 on B1: 2 x 8 donors; N7 on B2: at most N7_MAX_MICE mice split in two; pairs of this
    level, drawn in turn) in at least SESOI_TARGET of n datasets: NO DETECTABLE EFFECT (not
    detected, TOST within ±SESOI) or NOT SUPPORTED against the declared direction. Counted
    exactly on the grid (decided after the first review: the null must be establishable). Returns
    the SESOI and whether one met the target; where none did, the grid's largest, recorded as not
    found (pilot.json "sesoi_found")."""
    _STATE.update(bgs=bgs, pilot=pilot)
    pairs = level_pairs(bgs[background], level)
    got = _map(_null_job, [(background, pairs[i % len(pairs)], i) for i in range(n)], workers)
    need = math.ceil(SESOI_TARGET * len(got))
    for sesoi in SESOI_GRID:
        if sum(against or (not hit and width < sesoi) for hit, against, width in got) >= need:
            return sesoi, True
    return SESOI_GRID[-1], False


def _outcome_job(job):
    """The oracle's outcomes on one public-seed dataset of a case (`oracle_outcomes`)."""
    name, variant, pair, i, salt = job
    bgs, pilot = _STATE["bgs"], _STATE["pilot"]
    cond = P.conditions()[name]
    e = pilot_entries(cond, variant, pair, i + 1, salt)[i]
    X, obs, _, cards = P.build(e, bgs, pilot)
    return oracle_outcomes(e, X, obs, cards[-1], bgs, pilot, np.random.default_rng([e["seed"], 1]))


def choose_dose(bgs: dict, level: str, pilot: dict, n: int = PILOT_DATASETS, draws: int = DELTA_DRAWS,
                workers: int = 1) -> tuple[float | None, list]:
    """On B1: the smallest grid dose at which the oracle (the metric's validity taken as known)
    detects E1's injected coupling, in the direction of each pair's Δ* at that dose, with power
    >= 0.90, and the level's Δ* >= 1.25 x SESOI; None if no dose meets both."""
    sesoi = float(pilot["sesoi"]["B1"][level])
    pairs = level_pairs(bgs["B1"], level)
    scan = []
    for dose in DOSE_GRID:  # ascending: the smallest dose that meets both
        dl = delta_level(bgs["B1"], level, dose, draws)
        trial = _stand_in(pilot, "B1", level, e_dose={**pilot.get("e_dose", {}), level: dose})
        trial["truth"] = {"B1": {str(k): {"class": "valid"} for k in pairs}}
        trial["truth_case"] = {}
        trial["delta"] = {k: {"1": dict(value=v)} for k, v in dl["per_pair"].items()}  # the cards' direction
        _STATE.update(bgs=bgs, pilot=trial)
        jobs = [("E1", "dose=key", pairs[i % len(pairs)], i, 4) for i in range(n)]
        outs = [o["best"] for o in _map(_outcome_job, jobs, workers)]
        power = sum(o == P.SUPPORTED for o in outs) / n
        scan.append(dict(dose=dose, power=power, delta=dl["value"], delta_se=dl["se"]))
        if power >= POWER and dl["value"] >= DELTA_MARGIN * sesoi:
            return dose, scan
    return None, scan


def _raw_job(job):
    pair, i = job
    bgs, pilot = _STATE["bgs"], _STATE["pilot"]
    bg = bgs["B1"]
    e = pilot_entries(P.conditions()["N2"], "c=0.5", pair, i + 1, 5)[i]
    X, obs, genes, _ = P.build(e, bgs, pilot)
    pe = bg.plan["pool"][pair]
    ia, ib = genes.index(pe["pair"][0]), genes.index(pe["pair"][1])
    return bool(detects(*per_donor(X, obs, lambda x: norm_pearson(x, ia, ib)))[0])


def raw_difference_power(bgs: dict, level: str, pilot: dict, n: int = PILOT_DATASETS, workers: int = 1) -> float:
    """The oracle's power to detect N2's raw (uncorrected) difference at c = 0.5 on this level."""
    _STATE.update(bgs=bgs, pilot=pilot)
    pairs = level_pairs(bgs["B1"], level)
    return float(np.mean(_map(_raw_job, [(pairs[i % len(pairs)], i) for i in range(n)], workers)))


def establishability(bgs: dict, pilot: dict, n: int = PILOT_DATASETS, workers: int = 1) -> dict:
    """Per condition, variant and pool pair: the share of pilot datasets on which the oracle's
    outcome is a correct definite one (establishable where it is >= 0.90), and the measurements
    the sound-validator model needs (oc.sound_model): the outcomes of the effect's analysis taken
    as valid and, where the oracle ran GATE 4's rule, GATE 4's outcomes. N4 is never establishable
    (no oracle); the useless metrics of N6 do not depend on the pair, only on its level's
    delta_min, and are computed once per level."""
    _STATE.update(bgs=bgs, pilot=pilot)
    out, jobs, cases = {}, [], []
    for cond in P.CONDITIONS:
        if cond.background not in bgs:  # dropped with its background (pilot.json "dropped")
            continue
        pool = bgs[cond.background].plan["pool"]
        for variant, _ in cond.variants:
            if not cond.oracle:
                for k in range(len(pool)):
                    out[f"{cond.name}:{variant}:{k}"] = dict(power=0.0, establishable=False, outcomes={})
                continue
            if cond.metric != "norm_pearson":
                reps = {lv: level_pairs(bgs[cond.background], lv)[0] for lv in P.LEVELS}
                for k, pe in enumerate(pool):
                    cases.append((f"{cond.name}:{variant}:{k}", cond, variant, reps[pe["level"]]))
            else:
                for k in range(len(pool)):
                    cases.append((f"{cond.name}:{variant}:{k}", cond, variant, k))
    todo = sorted({(c.name, v, k) for _, c, v, k in cases})
    for name, variant, k in todo:
        jobs += [(name, variant, k, i, 6) for i in range(n)]
    outs = _map(_outcome_job, jobs, workers)
    by = {}
    for (name, variant, k, _, _), o in zip(jobs, outs):
        by.setdefault((name, variant, k), []).append(o)

    def counts(xs):
        xs = [x for x in xs if x is not None]
        return {x: xs.count(x) for x in sorted(set(xs))}
    for key, cond, variant, k in cases:
        got = by[(cond.name, variant, k)]
        good = P.definite(cond, variant, int(key.split(":")[-1]), pilot)
        hits = sum(o["best"] in good for o in got)
        # a sound validator that follows the engine's verdict order, on the datasets where GATE 4's
        # rule ran (engine_outcome)
        engine = [engine_outcome(o["effect"], o["gate4"]) for o in got if o["gate4"] is not None]
        out[key] = dict(power=hits / len(got), establishable=hits >= POWER * len(got), n=len(got),
                        outcomes=counts([o["best"] for o in got]),
                        effect_outcomes=counts([o["effect"] for o in got]),
                        gate4_outcomes=counts([o["gate4"] for o in got]),
                        engine_outcomes=counts(engine), computed_on_pair=k)
    return out


# --------------------------------------------------------------------------- #
# the pilot
# --------------------------------------------------------------------------- #
def run_pilot(bgs: dict, n: int = PILOT_DATASETS, draws: int = DELTA_DRAWS, workers: int = 1,
              sizes: dict | None = None) -> dict:
    """The pilot.json of the backgrounds `bgs` (module docstring); `sizes` overrides SIZES (the
    smoke test's), for this call only."""
    sizes = dict(SIZES, **(sizes or {}))
    _STATE["sizes"] = sizes
    try:
        return _run_pilot(bgs, n, draws, workers, sizes)
    finally:
        _STATE.pop("sizes", None)  # a later call in the same process gets the default sizes


def _run_pilot(bgs: dict, n: int, draws: int, workers: int, sizes: dict) -> dict:
    dropped = P.check_dropped(P.background_drops(bgs))  # a background without a candidate, with its cases
    pilot = dict(pilot_seed=PILOT_SEED, datasets_per_case=n, delta_draws=draws, alpha=ALPHA,
                 power_threshold=POWER, sesoi_target=SESOI_TARGET, delta_margin=DELTA_MARGIN,
                 delta_min_fraction=P.DELTA_MIN_FRACTION, truth_band=P.TRUTH_BAND,
                 saturation_share=SATURATION_SHARE, gate4_reps=sizes["gate4"], sizes=sizes,
                 dose_grid=list(DOSE_GRID), dropped=dropped,
                 pool={b: [dict(index=pe["index"], level=pe["level"], pair=pe["pair"]) for pe in bg.plan["pool"]]
                       for b, bg in bgs.items()},
                 pool_size={b: len(bg.plan["pool"]) for b, bg in bgs.items()},
                 sesoi={b: {} for b in bgs}, sesoi_found={b: {} for b in bgs}, response_curve={b: {} for b in bgs},
                 saturation_dose={b: {} for b in bgs}, truth={})
    for b in bgs:
        for level in P.LEVELS:
            pilot["sesoi"][b][level], pilot["sesoi_found"][b][level] = choose_sesoi(
                bgs, b, level, pilot, int(sizes["sesoi"]), workers)
    for b in bgs:
        for level in P.LEVELS:
            curve = response_curve(bgs, b, level, pilot, workers)
            pilot["response_curve"][b][level] = curve
            pilot["saturation_dose"][b][level] = saturation_dose(curve)
        pilot["truth"][b] = pair_truth(bgs, b, pilot, workers)
    pilot.update(key_dose={}, key_dose_found={}, e_dose={}, dose_scan={}, delta={}, delta_level={},
                 n2_raw_power={})
    for level in P.LEVELS:
        dose, scan = choose_dose(bgs, level, pilot, n, draws, workers)
        pilot["key_dose"][level], pilot["key_dose_found"][level] = dose, dose is not None
        pilot["e_dose"][level] = dose if dose is not None else pilot["saturation_dose"]["B1"][level]
        pilot["dose_scan"][level] = scan
    bg = bgs["B1"]
    for k, pe in enumerate(bg.plan["pool"]):
        pilot["delta"][str(k)] = {f"{f:g}": delta_star(bg, k, f * pilot["e_dose"][pe["level"]], draws)
                                  for f in E1_FACTORS}
    for level in P.LEVELS:
        idx = [str(k) for k in level_pairs(bg, level)]
        pilot["delta_level"][level] = {
            f"{f:g}": dict(value=float(np.mean([pilot["delta"][k][f"{f:g}"]["value"] for k in idx])),
                           se=float(math.sqrt(sum(pilot["delta"][k][f"{f:g}"]["se"] ** 2 for k in idx)) / len(idx)))
            for f in E1_FACTORS}
        pilot["n2_raw_power"][level] = raw_difference_power(bgs, level, pilot, n, workers)
    # reported only (v1.md section 3.2): the levels where N2's capture loss is visible in the raw
    # difference with the oracle's power, so that N2 there tests the engine's correction of depth
    pilot["informative_levels"] = [lv for lv in P.LEVELS if pilot["n2_raw_power"][lv] >= POWER]
    pilot["truth_case"] = case_truth(bgs, pilot, workers)
    pilot["establishable"] = establishability(bgs, pilot, n, workers)
    pilot["backgrounds"] = {k: P.background_sha256(b) for k, b in bgs.items()}
    return pilot


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--backgrounds", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--datasets", type=int,
                   help=f"per establishability case (default {PILOT_DATASETS}; {SMOKE_DATASETS} with --smoke)")
    p.add_argument("--draws", type=int, help=f"donor draws for Δ* (default {DELTA_DRAWS}; {SMOKE_DRAWS} with --smoke)")
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--data-dir", help="the downloaded files of backgrounds.json (default: next to it)")
    p.add_argument("--smoke", action="store_true",
                   help="the dry run's smoke test: the same code with tiny simulations (not a pilot)")
    args = p.parse_args(argv)
    n = args.datasets or (SMOKE_DATASETS if args.smoke else PILOT_DATASETS)
    draws = args.draws or (SMOKE_DRAWS if args.smoke else DELTA_DRAWS)
    pilot = run_pilot(P.load_backgrounds(args.backgrounds, args.data_dir), n, draws, args.workers,
                      SMOKE_SIZES if args.smoke else None)
    if args.smoke:
        pilot["smoke"] = True
    Path(args.out).write_text(json.dumps(pilot, indent=1))
    print(json.dumps({k: v for k, v in pilot.items() if k not in ("establishable", "response_curve", "pool")},
                     indent=1))


if __name__ == "__main__":
    main()
