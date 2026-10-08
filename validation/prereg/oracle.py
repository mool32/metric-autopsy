"""The oracle and the pilot of the confirmatory validation (validation/prereg/v1.md, 3.2 and 8).

The oracle is told what the engine has to find out. It takes the metric's validity as known,
removes a planted capture loss with the true factor (the other side is thinned by the same
factor) and a planted dropout by applying the same dropout to the other side, and tests the
per-replicate metric values with the test the graded replicate rule prescribes: an exact (or
Monte Carlo) permutation over replicates at >= 4 per group, Welch's t at 3; paired designs by
sign flips or the paired t; equivalence by TOST, the (1 - 2 alpha) interval inside ±SESOI. Its
power is the most a validator could reach, and only it decides which conditions are
establishable: the engine never does.

Independence: numpy and the t quantile of scipy.stats only; this file never imports
``metric_autopsy``. The pilot runs on datasets drawn with seeds outside the panel's (they come
from the public PILOT_SEED, the panel's from the owner's key seed).

    python validation/prereg/oracle.py --backgrounds backgrounds.json --out pilot.json
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
PILOT_DATASETS = 100  # per condition and variant
SESOI_GRID = tuple(round(0.005 * k, 3) for k in range(1, 201))
# The coupling dose enters a logistic keep probability; far above 4 the thinning becomes almost
# all-or-none, the co-detected cells are those that kept everything, and the coupling a
# co-detection metric sees fades again. The grid stays in the range where power rises.
DOSE_GRID = tuple(round(0.25 * k, 2) for k in range(1, 17))
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
    hits = 0
    for signs in itertools.product((1, -1), repeat=len(d)):
        hits += abs((d * np.asarray(signs)).mean()) >= obs - 1e-12 * max(1.0, obs)
    return hits / 2 ** len(d)


def detects(a: np.ndarray, b: np.ndarray, alpha: float = ALPHA) -> tuple[bool, float]:
    """The graded rule's test of a difference: permutation at >= 4 per group, Welch at 3."""
    if min(len(a), len(b)) >= 4:
        return permutation_p(a, b) < alpha, float(a.mean() - b.mean())
    est, se, df = welch(a, b)
    return (se > 0 and 2 * (1 - _t_cdf(abs(est) / se, df)) < alpha), est


def _t_cdf(t: float, df: float) -> float:
    from scipy import stats
    return float(stats.t.cdf(t, df))


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
    for d in dict.fromkeys(obs["donor"]):
        m = np.asarray(obs["donor"]) == d
        g = str(np.asarray(obs["group"])[m][0])
        vals[g].append(fn(X[m]))
    return np.asarray(vals["A"]), np.asarray(vals["B"])


def per_donor_paired(X, obs, fn) -> np.ndarray:
    out = []
    donor, group = np.asarray(obs["donor"]), np.asarray(obs["group"])
    for d in dict.fromkeys(donor):
        m = donor == d
        out.append(fn(X[m & (group == "A")]) - fn(X[m & (group == "B")]))
    return np.asarray(out)


def equalize(X, obs, side: str, kind: str, amount: float, rng) -> np.ndarray:
    """Apply the planted nuisance to the other side as well (the true correction)."""
    X = X.copy()
    other = np.asarray(obs["group"]) != side
    if kind == "capture":
        X[other] = P.thin(X[other], amount, rng)
    elif kind == "dropout":
        X[other] = P.drop_out(X[other], amount, rng)
    return X


# --------------------------------------------------------------------------- #
# the oracle's verdict on one dataset of known truth
# --------------------------------------------------------------------------- #
def oracle_reaches(entry: dict, X, obs, bg: P.Background, pilot: dict, rng) -> bool:
    """Whether the oracle reaches the condition's definite verdict on this dataset."""
    cond = P.conditions()[entry["condition"]]
    genes = bg.plan["genes"]
    ia, ib = genes.index(bg.plan["pair"][0]), genes.index(bg.plan["pair"][1])
    fn = (lambda x: norm_pearson(x, ia, ib))
    sesoi = float(pilot["sesoi"])
    name, variant, side = cond.name, entry["variant"], entry["side"]
    if name in ("N1", "N7"):
        a, b = per_donor(X, obs, fn)
        return tost_width(a, b) < sesoi
    if name == "N5":
        return paired_tost_width(per_donor_paired(X, obs, fn)) < sesoi
    if name in ("N2", "N3"):
        kind, amount = ("capture", float(variant.split("=")[1])) if name == "N2" else \
                       ("dropout", float(variant.split("=")[1]))
        raw = detects(*per_donor(X, obs, fn))[0]
        fixed = equalize(X, obs, side, kind, amount, rng)
        return raw and not detects(*per_donor(fixed, obs, fn))[0]
    if name == "N6b":
        return True
    if name in ("N6a", "N6c"):
        return valid_metric_responds(X, bg, kind="coupling" if name == "N6a" else "module",
                                     dose=float(pilot["key_dose"]), rng=rng)
    if name.startswith("E"):
        if name in ("E2", "E3"):
            X = equalize(X, obs, side if name == "E2" else ("B" if side == "A" else "A"), "capture", 0.5, rng)
        a, b = per_donor(X, obs, fn)
        hit, est = detects(a, b)
        higher_is_a = est > 0
        return hit and (higher_is_a == (side == "A"))
    return False  # N4: not establishable


def valid_metric_responds(X, bg: P.Background, kind: str, dose: float, rng, k: int = 10) -> bool:
    """For a useless metric (N6) the oracle is the valid metric facing the same injected signal:
    it must respond, 10 injections against 10 shams, by more than 4 standard errors."""
    genes = bg.plan["genes"]
    if kind == "coupling":
        ia, ib = genes.index(bg.plan["pair"][0]), genes.index(bg.plan["pair"][1])
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
# the pilot
# --------------------------------------------------------------------------- #
def pilot_entries(cond: P.Condition, variant: str, n: int) -> list[dict]:
    ss = np.random.SeedSequence([PILOT_SEED, P.CONDITIONS.index(cond), cond.variants.index(
        next(v for v in cond.variants if v[0] == variant))])
    out = []
    for i, child in enumerate(ss.spawn(n)):
        rng = np.random.default_rng(child)
        out.append(dict(id=f"P{i}", condition=cond.name, variant=variant, index=i,
                        side=("A", "B")[int(rng.integers(2))], seed=int(child.generate_state(1)[0])))
    return out


def choose_sesoi(bgs: dict, n: int = PILOT_DATASETS) -> float:
    """The smallest SESOI on the grid at which the oracle's TOST establishes equivalence on
    N1 (2 x 8 donors, no biology) with power >= 0.90."""
    cond = P.conditions()["N1"]
    widths = []
    for e in pilot_entries(cond, "null", n):
        X, obs, _, _ = P.build(e, bgs, dict(sesoi=1.0, key_dose=1.0))
        bg = bgs["B1"]
        genes = bg.plan["genes"]
        ia, ib = genes.index(bg.plan["pair"][0]), genes.index(bg.plan["pair"][1])
        widths.append(tost_width(*per_donor(X, obs, lambda x: norm_pearson(x, ia, ib))))
    q = float(np.quantile(widths, POWER))
    return next(s for s in SESOI_GRID if s > q)


def choose_dose(bgs: dict, sesoi: float, n: int = PILOT_DATASETS) -> float:
    """The smallest dose on the grid at which the oracle detects E1's injected coupling, in the
    right direction, with power >= 0.90."""
    cond = P.conditions()["E1"]
    entries = pilot_entries(cond, "dose=key", n)

    def power(dose):
        hits = 0
        for e in entries:
            X, obs, _, _ = P.build(e, bgs, dict(sesoi=sesoi, key_dose=dose))
            hits += oracle_reaches(e, X, obs, bgs["B1"], dict(sesoi=sesoi, key_dose=dose),
                                   np.random.default_rng(e["seed"] + 1))
        return hits / n

    for dose in DOSE_GRID:  # ascending: the smallest dose that reaches the power
        if power(dose) >= POWER:
            return dose
    raise RuntimeError("no dose on the grid reaches the oracle's power 0.90")


def establishability(bgs: dict, pilot: dict, n: int = PILOT_DATASETS) -> dict:
    """Oracle power of every condition and variant; establishable where it is >= 0.90."""
    out = {}
    for cond in P.CONDITIONS:
        for variant, _ in cond.variants:
            if cond.definite is None:
                out[f"{cond.name}:{variant}"] = dict(power=0.0, establishable=False)
                continue
            hits = 0
            for e in pilot_entries(cond, variant, n):
                X, obs, _, _ = P.build(e, bgs, pilot)
                hits += oracle_reaches(e, X, obs, bgs[cond.background], pilot,
                                       np.random.default_rng(e["seed"] + 1))
            out[f"{cond.name}:{variant}"] = dict(power=hits / n, establishable=hits >= POWER * n)
    return out


def run_pilot(bgs: dict, n: int = PILOT_DATASETS) -> dict:
    sesoi = choose_sesoi(bgs, n)
    dose = choose_dose(bgs, sesoi, n)
    pilot = dict(sesoi=sesoi, key_dose=dose, pilot_seed=PILOT_SEED, datasets_per_condition=n)
    pilot["establishable"] = establishability(bgs, pilot, n)
    return pilot


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--backgrounds", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--datasets", type=int, default=PILOT_DATASETS)
    args = p.parse_args(argv)
    spec = json.loads(Path(args.backgrounds).read_text())
    bgs = {k: P.load_background(k, v) for k, v in spec.items()}
    for bg in bgs.values():
        P.plan_background(bg)
    pilot = run_pilot(bgs, args.datasets)
    Path(args.out).write_text(json.dumps(pilot, indent=1))
    print(json.dumps(pilot, indent=1))


if __name__ == "__main__":
    main()
