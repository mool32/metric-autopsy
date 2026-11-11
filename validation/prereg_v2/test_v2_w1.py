"""Tests of validation v2's first step (v2.md section 11, W1): `num.py`, `exchange.py` and
`comparators.py`.

Those files never import the engine (checked below); this file does, to check that their copies of the
engine's rules give the engine's numbers: the t distribution, the reference metrics, the graded test over
replicates with its TOST and floor, and, end to end, the raw replicate-level test of a `run_autopsy`
report. The exchange is checked against the formula and the numbers of
`validation/exploratory/v1_correction_rule.log` (section 4, as quoted in
`validation/plans/v2_requirements.md` section 1) and against the bootstrap of its script; the design
effects against `validation/prereg/score.py` and `validation/exploratory/v1_pair_deff.py`.
"""
from __future__ import annotations

import ast
import math
import sys
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import comparators as C  # noqa: E402
import exchange as E  # noqa: E402
import num as N  # noqa: E402

from metric_autopsy import SimpleData, metrics, run_autopsy  # noqa: E402
from metric_autopsy import effect as engine_effect  # noqa: E402
from metric_autopsy import report as engine_report  # noqa: E402
from metric_autopsy import stats as engine_stats  # noqa: E402

OWN_FILES = ("num.py", "exchange.py", "comparators.py")
ALLOWED_IMPORTS = {"__future__", "math", "itertools", "dataclasses", "numpy", "num", "scipy", "scipy.stats"}


# --------------------------------------------------------------------------- #
# the files of W1 do not import the engine
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name", OWN_FILES)
def test_w1_files_import_neither_the_engine_nor_anything_else(name):
    tree = ast.parse((HERE / name).read_text())
    seen = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            seen |= {a.name for a in node.names}
        elif isinstance(node, ast.ImportFrom):
            seen.add(node.module)
    assert seen <= ALLOWED_IMPORTS, seen - ALLOWED_IMPORTS
    assert not any("metric_autopsy" in s for s in seen)


# --------------------------------------------------------------------------- #
# num: the engine's t and beta functions, number for number
# --------------------------------------------------------------------------- #
def test_t_and_beta_functions_are_the_engines():
    for df in (1.0, 2.0, 2.7, 3.0, 5.5, 14.0, 60.0):
        for t in (-4.0, -1.3, 0.0, 0.4, 2.1, 7.0):
            assert N.t_cdf(t, df) == engine_stats.t_cdf(t, df)
        for q in (0.025, 0.05, 0.5, 0.9, 0.95, 0.975):
            assert N.t_ppf(q, df) == engine_stats.t_ppf(q, df)
    for a, b in ((0.5, 3.0), (2.0, 7.0), (10.0, 1.5)):
        for x in (0.01, 0.3, 0.77, 0.99):
            assert N.betainc(a, b, x) == engine_stats._betainc(a, b, x)
        for q in (0.025, 0.5, 0.975):
            assert N.beta_ppf(q, a, b) == engine_stats._beta_ppf(q, a, b)


def test_clopper_pearson_matches_the_engines_and_takes_an_effective_size():
    for k, n in ((0, 10), (3, 10), (10, 10), (11, 253), (40, 2445)):
        lo, hi = N.clopper_pearson(k, n)
        elo, ehi = engine_stats.clopper_pearson(k, n)
        assert lo == pytest.approx(elo, abs=1e-9) and hi == pytest.approx(ehi, abs=1e-9)
    lo, hi = N.clopper_pearson(5.5, 40.25)  # Korn-Graubard: a non-integer effective k and n
    assert 0 < lo < 5.5 / 40.25 < hi < 1


# --------------------------------------------------------------------------- #
# exchange: the calls, the scoring, pi* and its intervals
# --------------------------------------------------------------------------- #
def test_every_engine_cause_makes_one_call():
    calls = {}
    for cause in engine_report.CAUSES:
        label = {"replicated": "SUPPORTED — replicated", "provisional": "SUPPORTED (provisional until replicated)",
                 "opposite_direction": "NOT SUPPORTED — the effect is in the direction opposite",
                 "no_detectable_effect": "NO DETECTABLE EFFECT — effect +0.01",
                 "effect_below_sesoi": "EFFECT BELOW SESOI — effect +0.01"}.get(cause, "INCONCLUSIVE — x")
        calls[cause] = E.engine_call(label, cause)
    assert {c for c, v in calls.items() if v == "effect"} == {"replicated", "provisional"}
    assert {c for c, v in calls.items() if v == "against"} == {"opposite_direction"}
    assert {c for c, v in calls.items() if v == "absent"} == {"no_detectable_effect"}
    assert {c for c, v in calls.items() if v == "below"} == {"effect_below_sesoi"}
    assert set(calls.values()) <= set(E.CALLS)
    with pytest.raises(ValueError):
        E.engine_call("INCONCLUSIVE — x", "replicated")  # a label that does not match its cause


def test_the_scoring_of_every_call():
    null_risk, null_blind = E.Truth("null", True), E.Truth("null", False)
    big, small, blind_big = E.Truth("effect", True, True), E.Truth("effect", True, False), E.Truth("effect", False, True)
    # null data: an effect in either direction is false; an absence with a blind metric is unwarranted
    assert E.kinds("effect", null_risk) == (True, False, False)
    assert E.kinds("against", null_blind) == (True, False, False)
    assert E.kinds("absent", null_risk) == (False, False, False)
    assert E.kinds("absent", null_blind) == (False, False, True)
    assert E.kinds("below", null_risk) == (False, False, False)  # change 14: no false effect
    # real effects
    assert E.kinds("against", small) == (True, False, False)
    assert E.kinds("absent", big) == (False, True, False)
    assert E.kinds("absent", small) == (False, False, False)
    assert E.kinds("absent", E.Truth("effect", False, False)) == (False, True, False)  # blind metric
    assert E.kinds("below", big) == (False, True, False)
    assert E.kinds("below", small) == (False, False, False)
    assert E.kinds("effect", blind_big) == (False, False, False)
    assert E.detection("effect", small) and not E.detection("below", small) and not E.detection("effect", null_risk)
    with pytest.raises(ValueError):
        E.Truth("null", True, True)


def test_pi_star_gives_v1s_numbers():
    """v2_requirements.md section 1 (v1_correction_rule.log section 4): at risk against (b), 40 of 2,445
    false conclusions removed on null cards, none on effect cards, 96 of 515 detections lost: 0.26; under
    change 6's rule 58 and 99: 0.33; blind or useless, 525 of 6,295, 21 of 1,285 and 136: 0.89."""
    assert round(E.pi_star(40 / 2445, 0 / 515, 96 / 515), 2) == 0.26
    assert round(E.pi_star(58 / 2445, 0 / 515, 99 / 515), 2) == 0.33
    assert round(E.pi_star(525 / 6295, 21 / 1285, 136 / 1285), 2) == 0.89
    assert E.pi_star(0.0, 0.1, 0.0) == 1.0   # nothing removed on null cards, nothing lost: pays at every share
    assert E.pi_star(0.0, 0.0, 0.4) == 0.0   # only lost detections: at no share
    assert math.isnan(E.pi_star(-0.1, 0.2, 0.0))


def _v1_section4(rs_null, dE, dS, units_of, seed, n_boot):
    """`v1_correction_rule.py` section 4's bootstrap, copied verbatim but for its inputs."""
    rng = np.random.default_rng(seed)
    is_null = np.asarray(rs_null)
    units = sorted(set(units_of))
    uidx = np.array([units.index(u) for u in units_of])
    n0, n1 = int(is_null.sum()), int((~is_null).sum())
    i0, i1 = np.flatnonzero(is_null), np.flatnonzero(~is_null)
    card = []
    for _ in range(n_boot):
        j0, j1 = rng.choice(i0, n0), rng.choice(i1, n1)
        card.append(E.pi_star(dE[j0].sum() / n0, dE[j1].sum() / n1, dS[j1].sum() / n1))
    agg = np.zeros((len(units), 5))
    np.add.at(agg, uidx, np.column_stack([is_null, ~is_null, dE * is_null, dE * ~is_null, dS]))
    pair = []
    for _ in range(n_boot):
        t = np.bincount(rng.integers(0, len(units), len(units)), minlength=len(units)) @ agg
        pair.append(E.pi_star(t[2] / t[0], t[3] / t[1], t[4] / t[1]) if t[0] and t[1] else float("nan"))
    return card, pair


def test_the_exchange_draws_as_v1s_script():
    rng = np.random.default_rng(7)
    n = 400
    truths, ref, an, units = [], [], [], []
    for i in range(n):
        null = i % 3 != 0
        t = E.Truth("null", bool(rng.random() < 0.6)) if null else E.Truth("effect", bool(rng.random() < 0.7),
                                                                           bool(rng.random() < 0.5))
        truths.append(t)
        ref.append(str(rng.choice(E.CALLS)))
        an.append(str(rng.choice(E.CALLS)))
        units.append(f"bg:{int(rng.integers(0, 37))}")
    out = E.exchange(ref, an, truths, units, rng=np.random.default_rng(20261012), n_boot=300)
    dE = np.array([float(E.wrong(a, t)) - float(E.wrong(b, t)) for a, b, t in zip(ref, an, truths)])
    dS = np.array([float(E.detection(a, t)) - float(E.detection(b, t)) for a, b, t in zip(ref, an, truths)])
    card, pair = _v1_section4([t.data == "null" for t in truths], dE, dS, units, 20261012, 300)
    assert out["ci_cards"]["lo"] == np.nanpercentile(card, 2.5)
    assert out["ci_clusters"]["hi"] == np.nanpercentile(pair, 97.5)
    n0 = sum(t.data == "null" for t in truths)
    assert out["null_cards"] == n0 and out["effect_cards"] == n - n0 and out["clusters"] == len(set(units))
    assert out["pi_star"] == E.pi_star(out["removed_null"] / n0, out["removed_effect"] / (n - n0),
                                       out["lost"] / (n - n0))


def test_the_design_effects_are_v1s():
    pytest.importorskip("scipy")
    sys.path.insert(0, str(ROOT / "validation" / "prereg"))
    sys.path.insert(0, str(ROOT / "validation" / "exploratory"))
    import score as v1_score  # noqa: E402
    import v1_pair_deff  # noqa: E402
    rng = np.random.default_rng(3)
    donors = [f"d{i}" for i in range(30)]
    for trial in range(4):
        n = 120
        y = (rng.random(n) < 0.2).astype(float)
        dsets = [list(rng.choice(donors, size=8, replace=False)) for _ in range(n)]
        units = [f"p{int(rng.integers(0, 9))}" for _ in range(n)]
        mine = E.overlap_deff(y, dsets)
        assert mine["deff"] == pytest.approx(v1_pair_deff.overlap_deff(y, dsets)["deff"], rel=1e-12)
        assert mine["deff"] == pytest.approx(v1_score.overlap_interval(y, dsets)["deff"], rel=1e-5)
        assert E.joint_deff(y, dsets, units) == pytest.approx(v1_pair_deff.joint_deff(y, dsets, units), rel=1e-12)
        r = E.rate(y, donor_sets=dsets, units=units)
        assert r["n_eff"] == pytest.approx(n / max(d["deff"] for d in r["deff"].values()))
        assert r["ci95_deff"][0] <= r["rate"] <= r["ci95_deff"][1]


# --------------------------------------------------------------------------- #
# comparators: the metrics, number for number with the engine's
# --------------------------------------------------------------------------- #
def _counts(rng, n_cells, n_genes, coupled=(0, 1), coupling=0.8, depth=1.0, pair_mean=None):
    """Poisson counts with a shared log-normal factor of strength `coupling` on the `coupled` genes
    (whose mean is `pair_mean` where given) and log-normal library sizes."""
    z = rng.normal(size=n_cells)
    base = np.exp(rng.normal(0.0, 1.0, n_genes))
    if pair_mean is not None:
        base[list(coupled)] = pair_mean
    size = np.exp(rng.normal(0.0, 0.3, n_cells)) * depth
    lam = base[None, :] * size[:, None] * np.ones((n_cells, n_genes))
    for j in coupled:
        lam[:, j] *= np.exp(coupling * z)
    return rng.poisson(lam).astype(float)


@pytest.mark.parametrize("sparse", [False, True])
def test_the_metrics_are_the_engines(sparse):
    rng = np.random.default_rng(11)
    names = [f"g{i}" for i in range(12)]
    for trial in range(6):
        n_cells = int(rng.integers(5, 300))
        X = _counts(rng, n_cells, len(names), depth=float(rng.choice([0.2, 1.0, 4.0])))
        if sparse:
            sp = pytest.importorskip("scipy.sparse")
            Xs = sp.csr_matrix(X)
        else:
            Xs = X
        data = SimpleData(Xs, pd.DataFrame(index=range(n_cells)), names)
        pair = C.PerCell.from_matrix(Xs, names, ["g0", "g1"])
        for name in ("norm_pearson", "mi_3bin"):
            engine = getattr(metrics, name)(data, gene_a="g0", gene_b="g1")
            mine = C.Metric(name)(pair)
            assert (math.isnan(engine) and math.isnan(mine)) or engine == mine, (name, engine, mine)
        genes = ["g0", "g1", "g2", "g3"]
        engine = metrics.set_norm_pearson(data, genes=genes)
        mine = C.Metric("set_norm_pearson")(C.PerCell.from_matrix(Xs, names, genes))
        assert (math.isnan(engine) and math.isnan(mine)) or engine == mine
    one = C.PerCell(np.array([[0.0], [2.0], [5.0]]), np.array([10.0, 20.0, 0.0]))
    assert C.Metric("marker_mean_counts")(one) == pytest.approx(7 / 3)
    assert C.Metric("marker_log1p_cp10k")(one) == pytest.approx(np.mean(np.log1p([0.0, 1000.0, 5e4])))
    assert C.Metric("constant")(one) == 0.0


# --------------------------------------------------------------------------- #
# comparators: (b), the engine's graded test coded apart
# --------------------------------------------------------------------------- #
def _same(a, b) -> bool:
    """Equal, or both NaN."""
    return (isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b)) or a == b


def _engine_infer(values, design, alpha=0.05, n_perm=1000, seed=0):
    d = dict(kind=design.kind, reps=list(design.reps),
             membership={r: ({"g0"} if r in design.first else {"g1"}) for r in design.reps},
             strata={r: design.stratum.get(r, "") for r in design.reps})
    return engine_effect._infer(values, d, ("g0", "g1"), "permutation", alpha, n_perm,
                                np.random.default_rng(seed), 20000, 3)


@pytest.mark.parametrize("n_first,n_second,n_strata", [(3, 3, 1), (4, 4, 1), (5, 3, 1), (4, 5, 2), (6, 6, 2),
                                                       (13, 12, 1), (12, 12, 2), (2, 6, 1)])
def test_the_graded_test_over_donors_is_the_engines(n_first, n_second, n_strata):
    rng = np.random.default_rng(n_first * 100 + n_second * 10 + n_strata)
    for trial in range(5):
        reps = [f"r{i}" for i in range(n_first + n_second)]
        first = set(reps[:n_first])
        stratum = {r: f"s{i % n_strata}" for i, r in enumerate(reps)}
        shift = float(rng.choice([0.0, 0.8, 2.5]))
        values = {r: float(rng.normal(shift if r in first else 0.0, 1.0)) for r in reps}
        if trial == 4:
            values[reps[0]] = float("nan")  # an unmeasured replicate is left out
        design = C.Design("nested", reps, first=first, stratum=stratum)
        mine = C.donor_test(values, design, direction="decrease", margin=0.5, rng=np.random.default_rng(trial))
        eng = _engine_infer(values, design, seed=trial)
        for key, theirs in (("est", eng["est"]), ("se", eng["se"]), ("df", eng["df"]), ("p_t", eng["p_t"]),
                            ("p_perm", eng["perm"].p), ("min_attainable_p", eng["perm"].min_attainable_p)):
            assert _same(float(mine[key]), float(theirs)), (key, mine[key], theirs)
        assert mine["perm_method"] == eng["perm"].method
        assert mine["tier"] == eng["tier"] and mine["detected"] == eng["detected"]
        for a, b in zip(mine["ci_2alpha"], eng["ci_2alpha"]):
            assert _same(float(a), float(b))


@pytest.mark.parametrize("n", [3, 5, 9, 16])
def test_the_paired_test_over_donors_is_the_engines(n):
    rng = np.random.default_rng(n)
    for trial in range(4):
        reps = [f"d{i}" for i in range(n)]
        values = {r: float(rng.normal(float(rng.choice([0.0, 1.0])), 1.0)) for r in reps}
        design = C.Design("paired", reps)
        mine = C.donor_test(values, design, direction="increase", margin=0.3, rng=np.random.default_rng(trial))
        eng = _engine_infer(values, design, seed=trial)
        for key, theirs in (("est", eng["est"]), ("se", eng["se"]), ("p_t", eng["p_t"]), ("p_perm", eng["perm"].p)):
            assert _same(float(mine[key]), float(theirs)), (key, mine[key], theirs)
        assert mine["tier"] == eng["tier"] and mine["detected"] == eng["detected"]


def test_the_floor_is_the_engines():
    for se, df in ((0.1, 3.0), (0.05, 7.4), (0.2, 30.0), (0.0, 4.0), (float("nan"), 3.0)):
        a, b = C.design_floor(se, df, 0.05), engine_effect.design_floor(se, df, 0.05)
        assert (math.isnan(a) and math.isnan(b)) or a == b


def test_the_calls_of_the_test_over_donors():
    reps = [f"r{i}" for i in range(10)]
    first = set(reps[:5])
    design = C.Design("nested", reps, first=first)
    big = {r: (3.0 if r in first else 0.0) + 0.01 * i for i, r in enumerate(reps)}
    assert C.donor_test(big, design, direction="decrease")["call"] == "effect"     # groups[0] higher: a decrease
    assert C.donor_test(big, design, direction="increase")["call"] == "against"
    assert C.donor_test(big, design, direction="two-sided")["call"] == "effect"
    assert C.donor_test(big, design, direction="decrease", margin=100.0)["call"] == "below"
    assert C.donor_test(big, design, direction="decrease", margin=100.0, margin_source="floor")["call"] == "effect"
    flat = {r: 0.001 * ((-1) ** i) for i, r in enumerate(reps)}
    assert C.donor_test(flat, design, direction="decrease", margin=0.5)["call"] == "absent"
    assert C.donor_test(flat, design, direction="decrease")["margin_source"] == "floor"
    noisy = {r: float(v) for r, v in zip(reps, np.random.default_rng(0).normal(0.0, 5.0, 10))}
    assert C.donor_test(noisy, design, direction="decrease", margin=0.01)["call"] == "none"
    two = C.Design("nested", reps[:4], first=set(reps[:2]))
    assert C.donor_test(big, two, direction="decrease")["call"] == "none"          # two a group: no verdict


def test_the_test_over_cells():
    """A coupled pair (40 genes, so that the pair is a small share of every cell's total) against an
    uncoupled one: detected, an effect for "decrease" and one against "increase"; a null: no claim."""
    rng = np.random.default_rng(5)
    names = [f"g{i}" for i in range(40)]
    X = np.vstack([_counts(rng, 300, 40, coupling=1.5, pair_mean=6.0),
                   _counts(rng, 300, 40, coupling=0.0, pair_mean=6.0)])
    cells = C.PerCell.from_matrix(X, names, ["g0", "g1"])
    first = np.r_[np.ones(300, bool), np.zeros(300, bool)]
    out = C.cell_test(C.Metric("norm_pearson"), cells, first, direction="decrease", n_perm=199,
                      rng=np.random.default_rng(1))
    assert out["call"] == "effect" and out["est"] > 0 and out["p"] < 0.05
    out = C.cell_test(C.Metric("norm_pearson"), cells, first, direction="increase", n_perm=199,
                      rng=np.random.default_rng(1))
    assert out["call"] == "against"
    null = np.vstack([_counts(rng, 300, 40, coupling=0.5, pair_mean=6.0) for _ in range(2)])
    out = C.cell_test(C.Metric("norm_pearson"), C.PerCell.from_matrix(null, names, ["g0", "g1"]), first,
                      direction="decrease", n_perm=199, rng=np.random.default_rng(2))
    assert out["call"] == "none" and out["p"] >= 0.05


# --------------------------------------------------------------------------- #
# end to end: (b) on a card's data is the engine's raw replicate-level test
# --------------------------------------------------------------------------- #
def _same_value(mine: float, engine) -> bool:
    """A replicate's value: the engine reports an unmeasured one as None, this code as NaN."""
    return (engine is None and not np.isfinite(mine)) or mine == engine


def _nested_data(seed, n_per_group=4, cells=120, coupling=(1.0, 0.3), n_genes=30):
    rng = np.random.default_rng(seed)
    Xs, rows = [], []
    for gi, grp in enumerate(("g0", "g1")):
        for m in range(n_per_group):
            Xs.append(_counts(rng, cells, n_genes, coupled=(0, 1), coupling=coupling[gi], pair_mean=6.0))
            rows += [dict(group=grp, donor=f"{grp}_{m}", sex="f" if m % 2 else "m")] * cells
    names = ["A", "B"] + [f"G{i}" for i in range(n_genes - 2)]
    return SimpleData(np.vstack(Xs), pd.DataFrame(rows), names)


def _paired_data(seed, n_donors=6, cells=100, coupling=(1.0, 0.3), n_genes=30):
    rng = np.random.default_rng(seed)
    Xs, rows = [], []
    for d in range(n_donors):
        for gi, grp in enumerate(("g0", "g1")):
            Xs.append(_counts(rng, cells, n_genes, coupled=(0, 1), coupling=coupling[gi]))
            rows += [dict(group=grp, donor=f"d{d}")] * cells
    names = ["A", "B"] + [f"G{i}" for i in range(n_genes - 2)]
    return SimpleData(np.vstack(Xs), pd.DataFrame(rows), names)


@pytest.mark.parametrize("seed,coupling,within,sesoi", [
    (1, (1.0, 0.3), (), 0.1),        # a real decrease of coupling, a declared SESOI on the metric's scale
    (2, (0.6, 0.6), (), 0.2),        # none
    (3, (0.6, 0.6), (), None),       # none, no SESOI: the floor is the margin
    (4, (1.0, 0.3), ("sex",), 0.1),  # strata
])
def test_b_reproduces_the_engines_raw_test_end_to_end(seed, coupling, within, sesoi):
    data = _nested_data(seed, coupling=coupling)
    prereg = dict(estimand="composition", direction="decrease")
    if sesoi is not None:
        prereg.update(sesoi=sesoi, sesoi_scale="observed")
    rep = run_autopsy(partial(metrics.norm_pearson, gene_a="A", gene_b="B"), data, group_col="group",
                      groups=("g0", "g1"), replicate_col="donor", within=list(within), prereg=prereg, seed=seed)
    det = rep.effect.detail
    rv = det["replicate_values"]
    reps = [r["replicate"] for r in rv]
    first = {r["replicate"] for r in rv if r["group"] == "g0"}
    stratum = {r["replicate"]: "\x1f".join(r["stratum"]) if r["stratum"] else "" for r in rv}
    cells = C.PerCell.from_matrix(data.X, data.var_names, ["A", "B"])
    values = C.donor_values(C.Metric("norm_pearson"), cells, data.obs["donor"], data.obs["group"] == "g0",
                            C.Design("nested", reps, first=first, stratum=stratum))
    assert all(_same_value(values[r["replicate"]], r["raw"]) for r in rv)
    assert det["unmeasured"] == []  # every donor measured here (the paired test leaves some out)
    margin = det["effective_sesoi_raw"]
    b = C.donor_test(values, C.Design("nested", reps, first=first, stratum=stratum), direction="decrease",
                     margin=margin if sesoi is not None else None)
    assert b["est"] == det["raw_effect"]
    assert b["detected"] == det["raw_detected"]
    assert b["p_perm" if b["tier"] == "permutation" else "p_t"] == det["p_raw"]
    assert b["ci_2alpha"] == pytest.approx(det["ci_2alpha_raw"], abs=0, rel=0)
    assert b["equivalent"] == det["tost_equivalent_raw"]
    if sesoi is None:
        assert b["margin"] == margin and b["margin_source"] == "floor"


def test_b_reproduces_the_engines_raw_test_on_a_paired_design():
    data = _paired_data(5)
    prereg = dict(estimand="composition", direction="decrease", sesoi=0.1, sesoi_scale="observed")
    rep = run_autopsy(partial(metrics.norm_pearson, gene_a="A", gene_b="B"), data, group_col="group",
                      groups=("g0", "g1"), replicate_col="donor", prereg=prereg, seed=5)
    det = rep.effect.detail
    assert det["design_kind"] == "paired"
    reps = [r["replicate"] for r in det["replicate_values"]]
    cells = C.PerCell.from_matrix(data.X, data.var_names, ["A", "B"])
    design = C.Design("paired", reps)
    values = C.donor_values(C.Metric("norm_pearson"), cells, data.obs["donor"], data.obs["group"] == "g0", design)
    assert all(_same_value(values[r["replicate"]], r["raw"]) for r in det["replicate_values"])
    assert len(det["unmeasured"]) == sum(1 for v in values.values() if not np.isfinite(v)) > 0
    b = C.donor_test(values, design, direction="decrease", margin=det["effective_sesoi_raw"])
    assert b["est"] == det["raw_effect"] and b["detected"] == det["raw_detected"]
    assert b["equivalent"] == det["tost_equivalent_raw"]
