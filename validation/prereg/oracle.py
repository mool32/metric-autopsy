"""The oracle and the pilot of the confirmatory validation (validation/prereg/v1.md, 3.2 and 8).

The oracle is told what the engine has to find out. It takes the metric's validity as known,
removes a planted nuisance by applying it to the other side as well (capture loss with the true
factor, dropout with the true fraction, N8's per-cell capture with fresh Beta(2, 2) draws), and
tests the per-replicate metric values with the test the graded replicate rule prescribes: an
exact (or Monte Carlo) permutation over replicates at >= 4 per group, Welch's t at 3; paired
designs by sign flips; equivalence by TOST, the (1 - 2 alpha) interval inside ±SESOI. Its
verdict follows the engine's order (explained by depth, reversed sign, detected in or against
the declared direction, equivalent, else INCONCLUSIVE). Its power is the most a validator could
reach, and only it decides which cases are establishable: the engine never does.

The pilot (before the key) fixes, per expression level of the pair pool: the SESOI (the smallest
grid value at which the oracle's TOST establishes N1 with power >= 0.9); Δ*, the population
difference of the metric under the injected coupling (signal side minus sham, a large paired
simulation; value and standard error); the key dose (the smallest grid dose at which the oracle
detects E1 with power >= 0.9 and the level's Δ* >= 1.25 x SESOI); the oracle's power to detect
N2's raw difference at c = 0.5 (N2 and N3 are informative only on a level where it is >= 0.9);
and the establishable (condition, variant, level) cases - those where the oracle reaches a
correct definite verdict (in the allowed set, not INCONCLUSIVE) with power >= 0.9.

Independence: numpy and the t distribution of scipy.stats only; this file never imports
``metric_autopsy``. The pilot runs on datasets drawn with seeds outside the panel's (they come
from the public PILOT_SEED, the panel's from the owner's key).

    python validation/prereg/oracle.py --backgrounds backgrounds.json --data-dir DATA --out pilot.json
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

import panel as P

ALPHA = 0.05
POWER = 0.90          # establishable: oracle power >= 0.90
PILOT_SEED = 20261008
PILOT_DATASETS = 100  # per condition, variant and level
DELTA_DRAWS = 2000    # paired donor draws per pair and dose for Δ*
DELTA_MARGIN = 1.25   # the key dose needs the level's Δ* >= 1.25 x SESOI
SESOI_GRID = tuple(round(0.005 * k, 3) for k in range(1, 201))
# The coupling dose enters a logistic keep probability; far above 4 the thinning becomes almost
# all-or-none, the co-detected cells are those that kept everything, and the coupling a
# co-detection metric sees fades again. The grid stays in the range where power rises.
DOSE_GRID = tuple(round(0.25 * k, 2) for k in range(1, 17))
E1_FACTORS = (0.25, 0.5, 1.0, 1.5)
MODULE_FOLD, MODULE_FRAC = 2.0, 0.3  # the engine's GATE 4 module injection in the claim cards


# --------------------------------------------------------------------------- #
# the metrics, re-implemented from their definitions
# --------------------------------------------------------------------------- #
def norm_pearson(X: np.ndarray, ia: int, ib: int) -> float:
    """CP10k + log1p Pearson of genes a and b over the cells where both are detected."""
    a, b = X[:, ia], X[:, ib]
    both = (a > 0) & (b > 0)
    if both.sum() < 3:
        return 0.0
    tot = X.sum(axis=1)
    tot = np.where(tot > 0, tot, 1.0)
    la, lb = np.log1p(a / tot * 1e4)[both], np.log1p(b / tot * 1e4)[both]
    if la.std() == 0 or lb.std() == 0:
        return 0.0
    return float(np.corrcoef(la, lb)[0, 1])


def module_score(X: np.ndarray, cols: list) -> float:
    tot = X.sum(axis=1, keepdims=True)
    return float(np.log1p(X[:, cols] / np.where(tot > 0, tot, 1.0) * 1e4).mean())


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
# the oracle's verdict on one dataset of known truth
# --------------------------------------------------------------------------- #
def label_from_values(raw, corr, corrected: bool, sesoi: float, direction: str) -> str:
    """The engine's decision order on replicate values, for an oracle that knows the metric is
    valid: a detected raw difference that the true correction removes (< half retained) is
    explained (NOT SUPPORTED); one whose sign the correction reverses is INCONCLUSIVE; a
    detected effect is SUPPORTED in the declared direction and NOT SUPPORTED against it;
    otherwise NO DETECTABLE EFFECT if the TOST establishes equivalence, else INCONCLUSIVE.
    `raw` and `corr` are (A values, B values); the effect is A - B, so > 0 is a decrease."""
    hit, est = detects(*corr)
    if corrected:
        raw_hit, raw_est = detects(*raw)
        retained = est / raw_est if raw_est else float("nan")
        if raw_hit and not hit and np.isfinite(retained) and abs(retained) < 0.5:
            return "NOT SUPPORTED"
        if raw_hit and hit and np.sign(est) != np.sign(raw_est):
            return "INCONCLUSIVE"
    if hit:
        observed = "decrease" if est > 0 else "increase"
        return "SUPPORTED" if observed == direction or direction == "two-sided" else "NOT SUPPORTED"
    return "NO DETECTABLE EFFECT" if tost_width(*corr) < sesoi else "INCONCLUSIVE"


def oracle_label(entry: dict, X, obs, card: dict, bg: P.Background, pilot: dict, rng) -> str:
    """The verdict the oracle reaches on this dataset (it knows the truth)."""
    cond = P.conditions()[entry["condition"]]
    pe = P.pool_entry(bg, entry)
    genes = list(bg.genes)
    ia, ib = genes.index(pe["pair"][0]), genes.index(pe["pair"][1])
    sesoi = float(pilot["sesoi"][pe["level"]])
    direction = card["prereg"]["direction"]
    if cond.name == "N6b":
        return "DEGENERATE METRIC"
    if cond.name in ("N6a", "N6c"):
        dose = float(pilot["key_dose"][pe["level"]])
        return ("NOT SUPPORTED" if valid_metric_responds(X, bg, pe, "coupling" if cond.name == "N6a" else "module",
                                                         dose, rng) else "INCONCLUSIVE")
    if not cond.oracle:
        return "INCONCLUSIVE"
    fn = (lambda x: norm_pearson(x, ia, ib))
    if cond.name == "N5":
        d = per_donor_paired(X, obs, fn)
        hit, est = detects_paired(d)
        if hit:
            observed = "decrease" if est > 0 else "increase"
            return "SUPPORTED" if observed == direction else "NOT SUPPORTED"
        return "NO DETECTABLE EFFECT" if paired_tost_width(d) < sesoi else "INCONCLUSIVE"
    raw = per_donor(X, obs, fn)
    fixed = true_correction(X, obs, entry, rng)
    corr = raw if fixed is None else per_donor(fixed, obs, fn)
    return label_from_values(raw, corr, fixed is not None, sesoi, direction)


def valid_metric_responds(X, bg: P.Background, pe: dict, kind: str, dose: float, rng, k: int = 10) -> bool:
    """For a useless metric (N6) the oracle is the valid metric facing the same injected signal:
    it must respond, 10 injections against 10 shams, by more than 4 standard errors."""
    genes = list(bg.genes)
    if kind == "coupling":
        ia, ib = genes.index(pe["pair"][0]), genes.index(pe["pair"][1])
        inj = [norm_pearson(P.inject_coupling(X, ia, ib, dose, rng), ia, ib) for _ in range(k)]
        sham = [norm_pearson(P.sham_coupling(X, ia, ib, dose, rng), ia, ib) for _ in range(k)]
    else:
        cols = [genes.index(g) for g in bg.plan["g2m"]]
        inj = [module_score(_inject_module(X, cols, rng), cols) for _ in range(k)]
        sham = [module_score(X, cols) for _ in range(k)]
    inj, sham = np.asarray(inj), np.asarray(sham)
    se = math.sqrt(inj.var(ddof=1) / k + sham.var(ddof=1) / k)
    return (inj.mean() - sham.mean()) > 4 * max(se, 1e-12)


def _inject_module(X, cols, rng):
    """Module up-regulation in a share of cells: every other gene thinned by 1/fold there."""
    out = X.copy()
    cells = rng.random(X.shape[0]) < MODULE_FRAC
    rest = np.setdiff1d(np.arange(X.shape[1]), cols)
    out[np.ix_(cells, rest)] = P.thin(X[np.ix_(cells, rest)], 1.0 / MODULE_FOLD, rng)
    return out


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
        diffs[t] = (norm_pearson(P.inject_coupling(Xd, ia, ib, dose, rng), ia, ib)
                    - norm_pearson(P.sham_coupling(Xd, ia, ib, dose, rng), ia, ib))
    return dict(value=float(diffs.mean()), se=float(diffs.std(ddof=1) / math.sqrt(draws)), draws=draws,
                dose=float(dose))


def delta_level(bg: P.Background, level: str, dose: float, draws: int = DELTA_DRAWS) -> dict:
    """The level's Δ*: the mean over its pairs (the key draws a level's pairs uniformly)."""
    idx = [k for k, pe in enumerate(bg.plan["pool"]) if pe["level"] == level]
    per = [delta_star(bg, k, dose, max(2, draws // len(idx))) for k in idx]
    return dict(value=float(np.mean([p["value"] for p in per])),
                se=float(math.sqrt(sum(p["se"] ** 2 for p in per)) / len(per)), pairs=idx, dose=float(dose))


# --------------------------------------------------------------------------- #
# the pilot
# --------------------------------------------------------------------------- #
def pilot_entries(cond: P.Condition, variant: str, level: str, n: int) -> list[dict]:
    """Public-seed datasets of one condition, variant and level (the pair drawn uniformly
    within the level, as the key draws pairs uniformly over the pool)."""
    vi = [v for v, _ in cond.variants].index(variant)
    ss = np.random.SeedSequence([PILOT_SEED, P.CONDITIONS.index(cond), vi, P.LEVELS.index(level)])
    out = []
    for i, child in enumerate(ss.spawn(n)):
        rng = np.random.default_rng(child)
        out.append(dict(id=f"P{i}", condition=cond.name, variant=variant, index=i,
                        side=("A", "B")[int(rng.integers(2))],
                        pair=P.LEVELS.index(level) * P.PAIRS_PER_LEVEL + int(rng.integers(P.PAIRS_PER_LEVEL)),
                        seed=int(child.generate_state(1)[0])))
    return out


def _labels(bgs: dict, cond: P.Condition, variant: str, level: str, pilot: dict, n: int) -> list[str]:
    out = []
    for e in pilot_entries(cond, variant, level, n):
        X, obs, _, cards = P.build(e, bgs, pilot)
        out.append(oracle_label(e, X, obs, cards[-1], bgs[cond.background], pilot,
                                np.random.default_rng([e["seed"], 1])))
    return out


def choose_sesoi(bgs: dict, level: str, n: int = PILOT_DATASETS) -> float:
    """The smallest SESOI on the grid at which the oracle's TOST establishes equivalence on N1
    (2 x 8 donors, no biology, pairs of this level) with power >= 0.90."""
    cond = P.conditions()["N1"]
    bg = bgs["B1"]
    widths = []
    for e in pilot_entries(cond, "null", level, n):
        X, obs, _, _ = P.build(e, bgs, dict(sesoi={lv: 1.0 for lv in P.LEVELS}, key_dose={}))
        pe = P.pool_entry(bg, e)
        ia, ib = bg.genes.index(pe["pair"][0]), bg.genes.index(pe["pair"][1])
        widths.append(tost_width(*per_donor(X, obs, lambda x: norm_pearson(x, ia, ib))))
    q = float(np.quantile(widths, POWER))
    return next((s for s in SESOI_GRID if s > q), SESOI_GRID[-1])


def choose_dose(bgs: dict, level: str, sesoi: float, n: int = PILOT_DATASETS,
                draws: int = DELTA_DRAWS) -> tuple[float, bool, list]:
    """The smallest dose on the grid at which the oracle detects E1's injected coupling, in the
    direction of Δ*, with power >= 0.90, and the level's Δ* >= 1.25 x SESOI. If no dose meets
    both, the largest grid dose is used and the level has no key dose (flag False)."""
    cond = P.conditions()["E1"]
    scan = []
    for dose in DOSE_GRID:  # ascending: the smallest dose that meets both
        dl = delta_level(bgs["B1"], level, dose, draws)
        pilot = dict(sesoi={lv: sesoi for lv in P.LEVELS}, key_dose={lv: dose for lv in P.LEVELS})
        labels = _labels(bgs, cond, "dose=key", level, pilot, n)
        power = sum(lab == "SUPPORTED" for lab in labels) / n
        scan.append(dict(dose=dose, power=power, delta=dl["value"], delta_se=dl["se"]))
        if power >= POWER and dl["value"] >= DELTA_MARGIN * sesoi:
            return dose, True, scan
    return DOSE_GRID[-1], False, scan


def raw_difference_power(bgs: dict, level: str, pilot: dict, n: int = PILOT_DATASETS) -> float:
    """The oracle's power to detect N2's raw (uncorrected) difference at c = 0.5 on this level."""
    cond = P.conditions()["N2"]
    bg = bgs["B1"]
    hits = 0
    for e in pilot_entries(cond, "c=0.5", level, n):
        X, obs, _, _ = P.build(e, bgs, pilot)
        pe = P.pool_entry(bg, e)
        ia, ib = bg.genes.index(pe["pair"][0]), bg.genes.index(pe["pair"][1])
        hits += detects(*per_donor(X, obs, lambda x: norm_pearson(x, ia, ib)))[0]
    return hits / n


def establishability(bgs: dict, pilot: dict, n: int = PILOT_DATASETS) -> dict:
    """Per condition, variant and level: the share of pilot datasets on which the oracle's
    verdict is a correct definite one; establishable where it is >= 0.90."""
    out = {}
    for cond in P.CONDITIONS:
        for variant, _ in cond.variants:
            for level in P.LEVELS:
                key = f"{cond.name}:{variant}:{level}"
                if not cond.oracle:
                    out[key] = dict(power=0.0, establishable=False, labels={})
                    continue
                entries = pilot_entries(cond, variant, level, n)
                labels = _labels(bgs, cond, variant, level, pilot, n)
                hits = sum(lab in P.definite(cond, variant, e["pair"], pilot) for e, lab in zip(entries, labels))
                out[key] = dict(power=hits / n, establishable=hits >= POWER * n,
                                labels={lab: labels.count(lab) for lab in sorted(set(labels))},
                                definite=sorted({v for e in entries for v in P.definite(cond, variant, e["pair"], pilot)}))
    return out


def run_pilot(bgs: dict, n: int = PILOT_DATASETS, draws: int = DELTA_DRAWS) -> dict:
    bg = bgs["B1"]
    pilot = dict(pilot_seed=PILOT_SEED, datasets_per_case=n, delta_draws=draws, alpha=ALPHA,
                 power_threshold=POWER, delta_margin=DELTA_MARGIN, sesoi={}, key_dose={},
                 key_dose_found={}, dose_scan={}, delta={}, delta_level={}, n2_raw_power={})
    for level in P.LEVELS:
        pilot["sesoi"][level] = choose_sesoi(bgs, level, n)
    for level in P.LEVELS:
        dose, found, scan = choose_dose(bgs, level, pilot["sesoi"][level], n, draws)
        pilot["key_dose"][level], pilot["key_dose_found"][level], pilot["dose_scan"][level] = dose, found, scan
    for k, pe in enumerate(bg.plan["pool"]):
        pilot["delta"][str(k)] = {f"{f:g}": delta_star(bg, k, f * pilot["key_dose"][pe["level"]], draws)
                                  for f in E1_FACTORS}
    for level in P.LEVELS:
        idx = [str(k) for k, pe in enumerate(bg.plan["pool"]) if pe["level"] == level]
        pilot["delta_level"][level] = {
            f"{f:g}": dict(value=float(np.mean([pilot["delta"][k][f"{f:g}"]["value"] for k in idx])),
                           se=float(math.sqrt(sum(pilot["delta"][k][f"{f:g}"]["se"] ** 2 for k in idx)) / len(idx)))
            for f in E1_FACTORS}
        pilot["n2_raw_power"][level] = raw_difference_power(bgs, level, pilot, n)
    pilot["informative_levels"] = [lv for lv in P.LEVELS if pilot["n2_raw_power"][lv] >= POWER]
    pilot["n2n3_informative"] = bool(pilot["informative_levels"])
    pilot["establishable"] = establishability(bgs, pilot, n)
    pilot["backgrounds"] = {k: P.background_sha256(b) for k, b in bgs.items()}
    return pilot


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--backgrounds", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--datasets", type=int, default=PILOT_DATASETS)
    p.add_argument("--draws", type=int, default=DELTA_DRAWS)
    p.add_argument("--data-dir", help="the downloaded files of backgrounds.json (default: next to it)")
    args = p.parse_args(argv)
    pilot = run_pilot(P.load_backgrounds(args.backgrounds, args.data_dir), args.datasets, args.draws)
    Path(args.out).write_text(json.dumps(pilot, indent=1))
    print(json.dumps({k: v for k, v in pilot.items() if k != "establishable"}, indent=1))


if __name__ == "__main__":
    main()
