"""The v0.3 verdict scheme: graded replicate rule, designs, estimands, power, provenance,
the injected-signal gate, the single decision rule, and API / CLI / MCP parity.

Every test builds data whose truth is known by construction (cells are iid within a
replicate unless stated), so the assertion is the verdict a correct validator must give.
"""
from __future__ import annotations

import json
from functools import partial

import numpy as np
import pandas as pd
import pytest

from metric_autopsy import (
    Assessment, GateStatus, SimpleData, __version__, decide, gate0_independence,
    gate4_signal_response, gate5_controls, injected_signal, metrics, run_autopsy,
)
from metric_autopsy import cli, mcp_server
from metric_autopsy import provenance as prov
from metric_autopsy.core import as_dense, gene_column
from metric_autopsy.equalize import thin_to_match
from metric_autopsy.qc import per_cell_qc
from metric_autopsy.report import Autopsy
from metric_autopsy.stats import (empirical_two_sided_p, extend_null, permutation_test_nested,
                                  permutation_test_paired, t_cdf)

from test_gates import MI, NPR, _assemble, _block, add_mice, make_clean, make_confounded

COMPOSITION = {"estimand": "composition"}


# --------------------------------------------------------------------------- #
# builders
# --------------------------------------------------------------------------- #
def nested_mice(n_per_group, c_young=2.6, c_old=0.3, cells=150, mouse_sd=0.0, seed=0):
    """Mice nested in age (one sex); per-mouse coupling ~ N(c_age, mouse_sd); QC matched."""
    rng = np.random.default_rng(seed)
    blocks = []
    for age, c in (("young", c_young), ("old", c_old)):
        for m in range(n_per_group):
            cm = max(rng.normal(c, mouse_sd) if mouse_sd > 0 else c, 0.05)
            X, o = _block(rng, cells, coupling=cm, efficiency=1.0, age=age, sex="female")
            o["mouse"] = f"{age}{m}"
            blocks.append((X, o))
    return _assemble(blocks)


def paired_donors(n_donors, c_treated=2.6, c_control=0.3, cells=150, seed=0):
    """Every donor contributes treated and control cells (a paired design)."""
    rng = np.random.default_rng(seed)
    blocks = []
    for d in range(n_donors):
        for cond, c in (("treated", c_treated), ("control", c_control)):
            X, o = _block(rng, cells, coupling=c, efficiency=1.0, age=cond, sex="female")
            o["donor"] = f"d{d}"
            blocks.append((X, o))
    data = _assemble(blocks)
    data.obs["cond"] = data.obs["age"]
    return data


def with_ercc(content=(1.0, 1.0), capture=(1.0, 1.0), n_mice=4, cells=150, n_genes=40,
              n_ercc=8, spikes=True, seed=0):
    """Endogenous counts ∝ content x capture; ERCC spike-ins ∝ capture only."""
    rng = np.random.default_rng(seed)
    mu = np.exp(rng.normal(0.5, 0.8, n_genes))
    Xs, rows = [], []
    for gi, age in enumerate(("young", "old")):
        for m in range(n_mice):
            size = np.exp(rng.normal(0, 0.2, cells))[:, None]
            endo = rng.poisson(mu[None, :] * content[gi] * capture[gi] * size)
            ercc = rng.poisson(np.full((cells, n_ercc), 3.0) * capture[gi])
            Xs.append(np.column_stack([endo, ercc]) if spikes else endo)
            rows += [dict(age=age, mouse=f"{age}{m}")] * cells
    names = [f"g{i}" for i in range(n_genes)] + ([f"ERCC-{i:05d}" for i in range(n_ercc)] if spikes else [])
    return SimpleData(np.vstack(Xs).astype(float), pd.DataFrame(rows), names)


def null_control_strata(n_strata, cells, seed):
    """Controls behave identically everywhere: Actb-Gapdh coupled, Gene0-Gene1 independent."""
    rng = np.random.default_rng(seed)
    blocks = []
    for k in range(n_strata):
        X, o = _block(rng, cells, coupling=1.5, efficiency=1.0, age="young", sex="x")
        o["stratum"] = f"s{k}"
        blocks.append((X, o))
    return _assemble(blocks)


def log_total_endogenous(data):
    """A content metric: mean log total endogenous counts per cell."""
    X = as_dense(data.X)
    endo = [j for j, g in enumerate(data.var_names) if not str(g).startswith("ERCC-")]
    return float(np.mean(np.log1p(X[:, endo].sum(axis=1))))


def _run(data, metric=NPR, **kw):
    kw.setdefault("group_col", "age")
    kw.setdefault("groups", ("young", "old"))
    kw.setdefault("prereg", COMPOSITION)
    return run_autopsy(metric, data, **kw)


# --------------------------------------------------------------------------- #
# graded replicate rule (decided 2026-10-07)
# --------------------------------------------------------------------------- #
def test_four_or_more_per_group_use_the_exact_permutation_test():
    a = _run(nested_mice(4), replicate_col="mouse")
    assert a.effect.status == "DETECTED"
    assert a.effect.detail["tier"] == "permutation"
    assert a.effect.detail["test"] == "exact permutation over replicates (n=70)"
    assert a.effect.detail["p_perm"] == pytest.approx(2 / 70)
    assert "PARAMETRIC_ONLY" not in a.effect.flags


def test_three_per_group_is_parametric_only_and_says_alpha_is_unattainable():
    a = _run(nested_mice(3), replicate_col="mouse")
    assert a.effect.status == "DETECTED"
    assert a.effect.detail["tier"] == "parametric"
    assert "t interval" in a.effect.detail["test"]
    assert "PARAMETRIC_ONLY" in a.effect.flags and "PARAMETRIC_ONLY" in a.design_adequacy.flags
    assert a.effect.detail["min_attainable_p"] == pytest.approx(0.1)  # 2 / C(6, 3)
    assert any("unattainable by construction" in n for n in a.design_adequacy.detail["notes"])


@pytest.mark.parametrize("n", [1, 2])
def test_two_or_fewer_per_group_give_no_effect_verdict(n):
    a = _run(nested_mice(n), replicate_col="mouse")
    assert a.design_adequacy.status == "INSUFFICIENT_REPLICATION"
    assert a.effect.status != "DETECTED"
    assert a.verdict.startswith("INCONCLUSIVE — insufficient replication")


def test_preregistered_minimum_replicates_is_enforced():
    a = _run(nested_mice(4), replicate_col="mouse", prereg={**COMPOSITION, "min_replicates": 5})
    assert a.design_adequacy.status == "INSUFFICIENT_REPLICATION"


def test_no_replicate_column_means_no_effect_verdict():
    a = _run(nested_mice(6))
    assert a.design_adequacy.status == "INSUFFICIENT_REPLICATION"
    assert a.effect.status == "INCONCLUSIVE"
    assert "replicate_col" in a.design_adequacy.reason


# --------------------------------------------------------------------------- #
# designs
# --------------------------------------------------------------------------- #
def test_paired_design_uses_sign_flips_over_donors():
    a = _run(paired_donors(6), group_col="cond", groups=("treated", "control"), replicate_col="donor")
    assert a.effect.detail["design_kind"] == "paired"
    assert a.effect.detail["test"] == "exact permutation over replicates (n=64)"
    assert a.effect.status == "DETECTED"


def test_paired_design_too_small_for_permutation_falls_back_to_parametric():
    a = _run(paired_donors(5), group_col="cond", groups=("treated", "control"), replicate_col="donor")
    assert a.effect.detail["tier"] == "parametric"  # 2 / 2**5 = 0.0625 > alpha
    assert "PARAMETRIC_ONLY" in a.effect.flags
    notes = a.design_adequacy.detail["notes"]
    assert len(notes) == 1 and "unattainable" in notes[0]


def test_partially_crossed_replicates_are_unidentifiable():
    d = paired_donors(6)
    obs = d.obs.copy()
    obs.loc[(obs["donor"] == "d0") & (obs["cond"] == "control"), "donor"] = "only-control"
    a = _run(SimpleData(d.X, obs, d.var_names), group_col="cond", groups=("treated", "control"),
             replicate_col="donor")
    assert a.design_adequacy.status == "UNIDENTIFIABLE"
    assert a.verdict.startswith("UNIDENTIFIABLE")


def test_no_estimand_is_unidentifiable():
    a = _run(nested_mice(4), replicate_col="mouse", prereg={})
    assert a.design_adequacy.status == "UNIDENTIFIABLE"
    assert "estimand" in a.verdict


# --------------------------------------------------------------------------- #
# estimands: composition (depth thinning) vs content (spike-in capture)
# --------------------------------------------------------------------------- #
def test_content_estimand_without_spikeins_is_unidentifiable():
    a = _run(with_ercc(capture=(1.0, 0.5), spikes=False), metric=log_total_endogenous,
             replicate_col="mouse", prereg={"estimand": "content"})
    assert a.design_adequacy.status == "UNIDENTIFIABLE"
    assert "spike-in" in a.design_adequacy.reason
    assert a.verdict.startswith("UNIDENTIFIABLE")


def test_content_estimand_removes_a_capture_difference():
    """Same RNA content, old capture halved: thinning to equal spike-in capture removes it."""
    a = _run(with_ercc(capture=(1.0, 0.5)), metric=log_total_endogenous, replicate_col="mouse",
             prereg={"estimand": "content"})
    assert a.design_adequacy.detail["correction"] == "capture_thinning"
    assert a.effect.detail["explained_by_depth"] is True
    assert abs(a.effect.detail["retained"]) < 0.1
    assert a.verdict.startswith("NOT SUPPORTED — the raw difference")
    assert "explained by capture" in a.verdict


def test_content_estimand_keeps_a_real_content_difference():
    """Young cells carry twice the RNA at equal capture: the difference must survive."""
    a = _run(with_ercc(content=(2.0, 1.0)), metric=log_total_endogenous, replicate_col="mouse",
             prereg={"estimand": "content"})
    assert a.design_adequacy.detail["correction"] == "capture_thinning"
    assert a.effect.status == "DETECTED"
    assert a.effect.detail["retained"] > 0.9


def test_depth_thinning_equalizes_depth_within_strata_and_leaves_the_other_group():
    d = make_confounded()
    totals = np.asarray(per_cell_qc(d)["total_counts"], dtype=float)
    eq, info = thin_to_match(d, "age", ("young", "old"), totals=totals, within=["sex"],
                             rng=np.random.default_rng(0))
    sex, age = np.asarray(d.obs["sex"]), np.asarray(d.obs["age"])
    eq_tot = eq.X.sum(1)
    ym, om = (sex == "male") & (age == "young"), (sex == "male") & (age == "old")
    assert np.median(eq_tot[ym]) == pytest.approx(np.median(eq_tot[om]), rel=0.1)
    assert np.array_equal(eq.X[om], d.X[om])  # the shallower group is untouched
    male = next(r for r in info["strata"] if r["stratum"] == {"sex": "male"})
    assert male["thinned_group"] == "young" and male["depth_ratio"] < 0.5


def pair_dominated(n=20000, eff_shallow=0.3, seed=0):
    """Identical biology in two groups; the coupled pair is a large share of each cell's total."""
    rng = np.random.default_rng(seed)
    parts = []
    for g, eff in (("deep", 1.0), ("shallow", eff_shallow)):
        lat = rng.normal(0, 1, n)
        smad = np.exp(1.2 + 0.9 * lat + rng.normal(0, 0.3, n))
        col = np.exp(1.2 + 0.9 * lat + rng.normal(0, 0.3, n))
        filler = np.exp(0.4 + rng.normal(0, 0.5, (n, 20)))
        parts.append((rng.poisson(np.column_stack([smad, col, filler]) * eff).astype(float),
                      pd.DataFrame({"grp": [g] * n})))
    return SimpleData(np.vstack([x for x, _ in parts]), pd.concat([o for _, o in parts], ignore_index=True),
                      ["Smad3", "Col1a1"] + [f"g{i}" for i in range(20)])


@pytest.mark.parametrize("seed", [0, 1])
def test_depth_thinning_keeps_the_joint_distribution_of_the_genes(seed):
    """Thinning the deeper group to the other's depth must reproduce the other group's MI when
    biology is identical. Matching each cell to depth *quantiles* failed this (+0.0145 here; on
    the p11b design 18% of a pure artifact survived, 7 SE): a cell's keep-probability then
    depends on its own total, which includes the genes the metric reads. One common ratio does not."""
    d = pair_dominated(seed=seed)
    tot = np.asarray(per_cell_qc(d)["total_counts"], dtype=float)
    lab = np.asarray(d.obs["grp"])
    eq, _ = thin_to_match(d, "grp", ("deep", "shallow"), totals=tot, rng=np.random.default_rng(seed + 1))
    mi = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")
    assert abs(mi(eq[lab == "deep"]) - mi(d[lab == "shallow"])) < 0.008


# --------------------------------------------------------------------------- #
# GATE 0 attenuation -> design adequacy (power), not metric validity
# --------------------------------------------------------------------------- #
def test_attenuation_is_a_power_check_in_design_adequacy():
    a = _run(cli.demo_data(), metric=MI, within=["sex"], replicate_col="mouse",
             gene_pair=("Smad3", "Col1a1"), prereg={**COMPOSITION, "sesoi": 0.05},
             stop_on_first_fail=False)
    assert "ATTENUATION" not in a.metric_validity.flags
    assert "attenuation" not in a.metric_validity.detail
    assert "ATTENUATION" in a.design_adequacy.flags
    assert a.design_adequacy.detail["attenuation"]["depth_downsample"] > 0.1
    power = a.design_adequacy.detail["power"]
    assert power["lambda_"] < 1 and power["effective_sesoi"] == pytest.approx(0.05 * power["lambda_"])


def test_underpowered_design_is_flagged_and_named_in_the_verdict():
    a = _run(nested_mice(4, c_young=1.5, c_old=1.4, mouse_sd=0.4, seed=3), replicate_col="mouse",
             prereg={**COMPOSITION, "sesoi": 0.02})
    assert "UNDERPOWERED" in a.design_adequacy.flags
    assert a.design_adequacy.detail["power"]["mde"] > a.design_adequacy.detail["power"]["effective_sesoi"]
    assert "underpowered relative to the SESOI" in a.verdict


def test_level_metric_is_not_failed_for_an_immaterial_null_signal():
    """A mean log total is moved slightly by the gene shuffle (cell-size covariance is lost).
    That immaterial 'signal' must not become the unit in which depth shifts are judged."""
    d = with_ercc(spikes=False)
    assert gate0_independence(log_total_endogenous, d, effect_scale=0.69).status != GateStatus.FAIL


# --------------------------------------------------------------------------- #
# GATE 4: response to an injected signal
# --------------------------------------------------------------------------- #
def test_injected_signal_separates_a_responsive_metric_from_useless_ones():
    d = make_clean()
    inject = injected_signal.coupling("Smad3", "Col1a1")
    rng = np.random.default_rng(0)
    assert gate4_signal_response(NPR, d, inject).status == GateStatus.PASS
    assert gate4_signal_response(lambda data: float(rng.normal()), d, inject).status == GateStatus.FAIL
    assert gate4_signal_response(lambda data: 1.0, d, inject).status == GateStatus.FAIL
    wrong_genes = injected_signal.coupling("Gene0", "Gene1")
    assert gate4_signal_response(NPR, d, wrong_genes).status == GateStatus.FAIL


def test_injected_signal_can_establish_metric_validity_without_controls():
    a = _run(add_mice(make_clean()), within=["sex"], replicate_col="mouse",
             gene_pair=("Smad3", "Col1a1"), signal_test=injected_signal.coupling("Smad3", "Col1a1"),
             prereg={**COMPOSITION, "judgment_pending": False})
    assert a.metric_validity.status == "PASS"
    assert "responds to coupling" in a.metric_validity.reason
    assert a.verdict == "SUPPORTED (provisional until replicated)"


def test_injected_module_raises_relative_expression_and_keeps_counts():
    d = make_clean()
    out = injected_signal.module(["Smad3", "Col1a1"], fold=2.0, frac=1.0)(d, np.random.default_rng(0))
    assert np.allclose(out.X, np.round(out.X)) and out.X.min() >= 0
    frac = lambda X: X[:, :2].sum() / X.sum()  # noqa: E731
    assert frac(out.X) > 1.5 * frac(d.X)


# --------------------------------------------------------------------------- #
# GATE 5: calibrated in small strata; a silent positive control is not a failure
# --------------------------------------------------------------------------- #
CTRL = (("Actb", "Gapdh"), ("Gene0", "Gene1"))


@pytest.mark.parametrize("n_strata, cells", [(1, 10), (4, 30)])
def test_gate5_does_not_fail_null_controls_in_small_strata(n_strata, cells):
    """Truth: both controls behave in every stratum. Before the fix, 10-cell strata failed
    100% of null datasets (a one-cell depth bin cannot be shuffled, so the positive control's
    null equalled the data) and 4 x 30 cells failed 97% (a silent positive control in an
    underpowered stratum counted as FAIL). FAIL must now stay near alpha."""
    results = [gate5_controls(metrics.norm_pearson, null_control_strata(n_strata, cells, s), *CTRL,
                              within=["stratum"]) for s in range(20)]
    assert sum(r.status == GateStatus.FAIL for r in results) <= 3
    assert all(r.status != GateStatus.FAIL or any(not row["neg_ok"] for row in r.detail["rows"])
               for r in results)  # only the negative control can fail the gate


def test_gate5_skewed_metric_is_calibrated_with_many_strata():
    """mi_3bin's null is right-skewed. With 16 strata alpha/K = 0.003 is below the resolution of
    200 draws; a normal tail there failed 17-33% of null datasets, and re-drawing the same few
    unrelated pairs gave a false resolution (10-23%). Both are gone."""
    results = [gate5_controls(metrics.mi_3bin, null_control_strata(16, 30, s), *CTRL,
                              within=["stratum"], n_power=10) for s in range(10)]
    assert sum(not row["neg_ok"] for r in results for row in r.detail["rows"]) <= 1


def test_matched_null_pairs_are_distinct_and_unrelated():
    from metric_autopsy.gates import _matched_pool
    d = make_clean()
    pool = _matched_pool(d, ("Gene0", "Gene1"), {"Gene0", "Gene1", "Actb", "Gapdh"}, np.random.default_rng(0))
    assert len(pool) == len(set(pool)) and all(a < b for a, b in pool)
    assert not ({g for pr in pool for g in pr} & {"Gene0", "Gene1", "Actb", "Gapdh"})


def test_monte_carlo_p_never_assumes_a_parametric_tail():
    null = np.random.default_rng(0).exponential(size=100)  # skewed
    p, method = empirical_two_sided_p(50.0, null, alpha=0.001)
    assert p == pytest.approx(1 / 101) and "unattainable" in method
    draws = iter([np.random.default_rng(1).exponential(size=200), np.random.default_rng(2).exponential(size=800)])
    null2, p2, method2 = extend_null(50.0, lambda k: next(draws)[:k], 200, alpha=0.0025)
    assert len(null2) == 800 and p2 == pytest.approx(1 / 801) and method2 == "empirical"


def test_positive_control_null_is_not_degenerate_in_a_tiny_stratum():
    res = gate5_controls(metrics.norm_pearson, null_control_strata(1, 10, 0), *CTRL)
    row = res.detail["rows"][0]
    assert row["pos_null_center"] != pytest.approx(row["pos"])


def blind_pair_metric(data, *, gene_a, gene_b):
    """A 'coupling' metric that ignores gene b entirely: blind to any coupling."""
    return float(np.mean(gene_column(data, gene_a)))


def test_controls_count_only_for_the_metric_they_test():
    """GATE 5 runs pair_metric on the control pairs. The API took on trust that pair_metric bound
    to gene_pair is the judged metric: at 677d1b8 a metric blind to gene b given norm_pearson's
    controls was certified and SUPPORTED (make_clean seeds 0-2, 3 of 3), and a random-number
    metric got metric validity PASS. The controls now count only when
    pair_metric(data, *gene_pair) == metric(data)."""
    d = add_mice(make_clean())
    kw = dict(within=["sex"], replicate_col="mouse", gene_pair=("Smad3", "Col1a1"),
              pair_metric=metrics.norm_pearson, pos_pair=("Actb", "Gapdh"), neg_pair=("Gene0", "Gene1"),
              prereg={**COMPOSITION, "judgment_pending": False}, log_path="off")
    blind = _run(d, metric=partial(blind_pair_metric, gene_a="Smad3", gene_b="Col1a1"), **kw)
    g5 = next(r for r in blind.results if r.gate == 5)
    assert g5.status == GateStatus.SKIP and g5.detail["controls_tied"] is False
    assert blind.metric_validity.status == "UNTESTED" and not blind.verdict.startswith("SUPPORTED")
    rng = np.random.default_rng(0)
    assert _run(d, metric=lambda data: float(rng.normal()), **kw).metric_validity.status == "UNTESTED"
    same = _run(d, **kw)  # NPR is norm_pearson on Smad3-Col1a1: its own controls count
    g5 = next(r for r in same.results if r.gate == 5)
    assert g5.status != GateStatus.SKIP and g5.detail["controls_tied"] is True
    assert same.metric_validity.status == "PASS"
    no_pair = _run(d, **{**kw, "gene_pair": None})  # nothing to tie the controls to
    assert next(r for r in no_pair.results if r.gate == 5).status == GateStatus.SKIP


def test_blind_metric_fails_where_the_design_has_power_and_is_untested_where_not():
    """Decided 2026-10-07: a silent positive control is UNTESTED only when the stratum lacked the
    power to show it; with power >= 0.8 for an injected coupling of the pre-registered dose it is
    evidence that the metric is insensitive (FAIL). Otherwise UNTESTED would shelter blind
    metrics forever. The power is the design's (a reference detector), not the metric's."""
    for s in range(3):
        powered = gate5_controls(blind_pair_metric, null_control_strata(1, 600, s), *CTRL)
        row = powered.detail["rows"][0]
        assert powered.status == GateStatus.FAIL and row["pos_insensitive"] and row["pos_power"] >= 0.8
        assert "insensitive" in powered.message
    tiny = gate5_controls(blind_pair_metric, null_control_strata(1, 10, 0), *CTRL)
    assert tiny.status == GateStatus.WARN and tiny.detail["pos_demonstrated"] is False
    assert tiny.detail["rows"][0]["pos_power"] < 0.8


def test_uncoupled_positive_control_in_a_powered_design_invalidates_the_metric_check():
    """Two independent genes offered as the positive control: in a stratum that could detect the
    dose, its silence FAILs GATE 5 (insensitive metric or a control that is not coupled), and
    the metric is not certified."""
    d = add_mice(make_clean())
    res = gate5_controls(metrics.norm_pearson, d, ("Gene2", "Gene3"), ("Gene0", "Gene1"), within=["sex"])
    assert res.status == GateStatus.FAIL and res.detail["pos_demonstrated"] is False
    a = _run(d, within=["sex"], replicate_col="mouse", gene_pair=("Smad3", "Col1a1"),
             pair_metric=metrics.norm_pearson, pos_pair=("Gene2", "Gene3"), neg_pair=("Gene0", "Gene1"))
    assert a.metric_validity.status == "FAIL" and not a.verdict.startswith("SUPPORTED")


def test_positive_control_power_settings_come_from_the_preregistration():
    a = _run(add_mice(make_clean()), within=["sex"], replicate_col="mouse", gene_pair=("Smad3", "Col1a1"),
             pair_metric=metrics.norm_pearson, pos_pair=("Gene2", "Gene3"), neg_pair=("Gene0", "Gene1"),
             prereg={**COMPOSITION, "positive_control_dose": 0.5, "positive_control_power": 0.95})
    g5 = next(r for r in a.results if r.gate == 5)
    assert g5.detail["pos_dose"] == 0.5 and g5.detail["pos_power_min"] == 0.95
    assert g5.status == GateStatus.WARN  # a weak dose cannot be established here: untested, not invalid


def depth_only_pair(n=600, n_genes=300, seed=0):
    """Genes A and B independent given cell size; size varies widely (log-SD 0.8)."""
    rng = np.random.default_rng(seed)
    size = np.exp(rng.normal(0, 0.8, n))
    mu = np.exp(rng.normal(0.5, 0.6, n_genes))
    X = rng.poisson(size[:, None] * mu[None, :] * 3).astype(float)
    obs = pd.DataFrame({"age": ["young"] * (n // 2) + ["old"] * (n - n // 2)})
    return SimpleData(X, obs, ["A", "B"] + [f"g{i}" for i in range(n_genes - 2)])


def test_pair_coupled_only_through_depth_is_not_a_positive_control():
    """The positive control's self-null must keep the dependence on depth. A pair coupled only
    through cell size correlates strongly on raw counts (a plain shuffle rejects it every time),
    but must not beat the depth-matched null more often than alpha. A shuffle within depth bins
    did not achieve this: residual within-bin depth variation let such a pair pass 5/5 times."""
    fired = 0
    for s in range(10):
        d = depth_only_pair(seed=s)
        res = gate5_controls(metrics.pearson, d, ("A", "B"), ("g0", "g1"))
        fired += res.detail["rows"][0]["pos_fires"]
        a, b = gene_column(d, "A"), gene_column(d, "B")
        rng = np.random.default_rng(s)
        plain = np.array([np.corrcoef(a, rng.permutation(b))[0, 1] for _ in range(200)])
        assert np.corrcoef(a, b)[0, 1] > plain.max()  # coupled on raw counts, through depth
    assert fired <= 2


def capture_confound_mice(n_per_group=6, cells=150, eff_old=0.5, seed=0):
    """Identical biology; old cells have half the capture (a pure depth artifact)."""
    rng = np.random.default_rng(seed)
    blocks = []
    for age, eff in (("young", 1.0), ("old", eff_old)):
        for m in range(n_per_group):
            X, o = _block(rng, cells, coupling=1.5, efficiency=eff, age=age, sex="female")
            o["mouse"] = f"{age}{m}"
            blocks.append((X, o))
    return _assemble(blocks)


def mean_lognorm_gene5(data):
    """A level metric: mean log1p(CP10k) of a low-expression gene (its zeros move with depth)."""
    X = as_dense(data.X)
    tot = X.sum(axis=1)
    tot[tot == 0] = 1
    return float(np.mean(np.log1p(X[:, list(data.var_names).index("Gene5")] / tot * 1e4)))


def test_level_metric_driven_by_depth_is_caught_by_the_effect_field():
    """A level metric (no gene-gene structure: GATE 0 can only report its level shift) on a pure
    depth artifact, with a demonstrated response to signal (metric PASS) and resolved judgment.
    The effect field must catch it after equalization: never SUPPORTED."""
    verdicts, explained = [], 0
    for s in range(10):
        a = run_autopsy(mean_lognorm_gene5, capture_confound_mice(seed=s), group_col="age",
                        groups=("young", "old"), replicate_col="mouse",
                        signal_test=injected_signal.module(["Gene5"], fold=2.0, frac=0.3),
                        prereg={**COMPOSITION, "judgment_pending": False})
        assert a.metric_validity.status == "PASS" and "LEVEL_SHIFT" in a.metric_validity.flags
        verdicts.append(a.verdict)
        explained += bool(a.effect.detail["explained_by_depth"])
    assert not any(v.startswith("SUPPORTED") for v in verdicts), verdicts
    assert explained >= 8


def test_effect_reversed_by_the_correction_is_inconclusive():
    """Seed 3 of the capture confound: raw +0.64 (detected) becomes -0.33 (detected) at equal
    depth. Before the rule this was a false SUPPORTED."""
    a = run_autopsy(mean_lognorm_gene5, capture_confound_mice(seed=3), group_col="age",
                    groups=("young", "old"), replicate_col="mouse",
                    signal_test=injected_signal.module(["Gene5"], fold=2.0, frac=0.3),
                    prereg={**COMPOSITION, "judgment_pending": False})
    assert a.effect.detail["reversed_by_correction"] is True
    assert a.effect.status == "INCONCLUSIVE" and a.verdict.startswith("INCONCLUSIVE")
    assert "reversed the sign" in a.verdict


def test_legacy_fixed_band_still_fails_a_silent_positive_control():
    res = gate5_controls(metrics.norm_pearson, make_clean(), ("Gene2", "Gene3"), ("Gene0", "Gene1"),
                         pos_min=0.3, neg_max=0.2)
    assert res.status == GateStatus.FAIL


# --------------------------------------------------------------------------- #
# decide(): the single rule
# --------------------------------------------------------------------------- #
def _autopsy(mv="PASS", da="ADEQUATE", ef="DETECTED", rp="NOT_RUN", pending=False, ef_flags=(),
             da_flags=(), ef_detail=None):
    return Autopsy("m", [], {"judgment_pending": pending},
                   Assessment(mv, "mv"), Assessment(da, "da", list(da_flags)),
                   Assessment(ef, "ef", list(ef_flags), dict(ef_detail or {})), Assessment(rp, "rp"))


@pytest.mark.parametrize("fields, prefix", [
    (dict(mv="DEGENERATE", da="UNIDENTIFIABLE"), "DEGENERATE METRIC"),
    (dict(mv="FAIL"), "NOT SUPPORTED — metric invalid"),
    (dict(da="UNIDENTIFIABLE"), "UNIDENTIFIABLE"),
    (dict(ef="INCONCLUSIVE", ef_detail={"explained_by_depth": True}), "NOT SUPPORTED"),
    (dict(da="INSUFFICIENT_REPLICATION", ef="INCONCLUSIVE"), "INCONCLUSIVE — insufficient replication"),
    (dict(ef="NO_DETECTABLE_EFFECT"), "NO DETECTABLE EFFECT"),
    (dict(mv="UNTESTED", ef="NO_DETECTABLE_EFFECT"), "INCONCLUSIVE"),
    (dict(mv="UNTESTED"), "INCONCLUSIVE — effect detected, but the metric's response"),
    (dict(rp="NOT_REPLICATED"), "NOT SUPPORTED — the effect did not replicate"),
    (dict(pending=True), "INCONCLUSIVE — effect detected; judgment gates"),
    (dict(rp="REPLICATED"), "SUPPORTED — replicated"),
    (dict(), "SUPPORTED (provisional until replicated)"),
])
def test_decide_order(fields, prefix):
    assert decide(_autopsy(**fields)).startswith(prefix)


def test_decide_names_the_assumption_the_verdict_rests_on():
    v = decide(_autopsy(ef_flags=["PARAMETRIC_ONLY"], da_flags=["UNDERPOWERED"]))
    assert v.startswith("SUPPORTED") and "parametric only" in v and "underpowered" in v


# --------------------------------------------------------------------------- #
# provenance: JSON report, hashes, run log
# --------------------------------------------------------------------------- #
def test_json_report_is_strict_and_carries_hashes():
    d = add_mice(make_clean())
    a = _run(d, within=["sex"], replicate_col="mouse", log_path="off")

    def no_nan(token):
        raise AssertionError(f"non-strict JSON token {token}")

    rep = json.loads(a.to_json(), parse_constant=no_nan)
    assert rep["schema"] == prov.SCHEMA and rep["verdict"] == a.verdict
    assert set(rep["fields"]) == {"metric_validity", "design_adequacy", "effect", "replication"}
    p = rep["provenance"]
    assert len(p["data_sha256"]) == len(p["prereg_sha256"]) == len(p["claim_id"]) == 64
    assert p["environment"]["metric_autopsy"] == __version__
    assert rep["params"]["replicate_col"] == "mouse"


def test_hashes_identify_data_and_preregistration():
    d = add_mice(make_clean())
    base = _run(d, within=["sex"], replicate_col="mouse", log_path="off").provenance
    again = _run(d, within=["sex"], replicate_col="mouse", log_path="off").provenance
    assert (base["data_sha256"], base["prereg_sha256"], base["claim_id"]) == \
           (again["data_sha256"], again["prereg_sha256"], again["claim_id"])
    X = d.X.copy()
    X[0, 0] += 1
    changed = _run(SimpleData(X, d.obs, d.var_names), within=["sex"], replicate_col="mouse",
                   log_path="off").provenance
    assert changed["data_sha256"] != base["data_sha256"]
    other_prereg = _run(d, within=["sex"], replicate_col="mouse", log_path="off",
                        prereg={**COMPOSITION, "sesoi": 0.1}).provenance
    assert other_prereg["prereg_sha256"] != base["prereg_sha256"]
    assert other_prereg["claim_id"] != base["claim_id"]


def test_run_log_counts_attempts_per_claim(tmp_path):
    log = tmp_path / "runs.jsonl"
    d = nested_mice(4)
    first = _run(d, replicate_col="mouse", log_path=log)
    second = _run(d, replicate_col="mouse", log_path=log, seed=1)
    other = _run(d, replicate_col="mouse", log_path=log, prereg={**COMPOSITION, "sesoi": 0.2})
    assert first.provenance["log"]["attempt"] == 1
    assert second.provenance["log"]["attempt"] == 2
    assert second.provenance["log"]["previous"][0]["verdict"] == first.verdict
    assert other.provenance["log"]["attempt"] == 1  # a different pre-registration is a new claim
    assert len(prov.read_log(log)) == 3


def test_run_log_is_on_with_a_preregistration_in_the_python_api(tmp_path, monkeypatch):
    """Agents work through the Python API, so a pre-registered run is logged there by default
    (decided 2026-10-07). Without a pre-registration there is no verdict to hunt: no log."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("METRIC_AUTOPSY_LOG", raising=False)
    d = nested_mice(4)
    bare = run_autopsy(NPR, d, group_col="age", groups=("young", "old"), replicate_col="mouse")
    assert bare.provenance["log"] is None and not (tmp_path / prov.DEFAULT_LOG).exists()
    first = _run(d, replicate_col="mouse")
    second = _run(d, replicate_col="mouse")
    assert (tmp_path / prov.DEFAULT_LOG).exists()
    assert (first.provenance["log"]["attempt"], second.provenance["log"]["attempt"]) == (1, 2)
    assert _run(d, replicate_col="mouse", log_path="off").provenance["log"] is None
    monkeypatch.setenv("METRIC_AUTOPSY_LOG", str(tmp_path / "env.jsonl"))
    via_env = _run(d, replicate_col="mouse")
    assert via_env.provenance["log"]["path"].endswith("env.jsonl")


def test_cli_demo_is_not_logged_unless_asked(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("METRIC_AUTOPSY_LOG", raising=False)
    cli.main(["--demo", "--no-stop"])
    assert not (tmp_path / prov.DEFAULT_LOG).exists()
    cli.main(["--demo", "--no-stop", "--log", str(tmp_path / "demo.jsonl")])
    assert len(prov.read_log(tmp_path / "demo.jsonl")) == 1
    capsys.readouterr()


# --------------------------------------------------------------------------- #
# stats
# --------------------------------------------------------------------------- #
def test_exact_permutation_p_values():
    y = np.array([5.0, 6, 7, 8, 1, 2, 3, 4])
    g = np.array([1, 1, 1, 1, 0, 0, 0, 0])
    r = permutation_test_nested(y, g, np.zeros(8, int))
    assert r.method == "exact" and r.n == 70 and r.p == pytest.approx(2 / 70)
    assert r.min_attainable_p == pytest.approx(2 / 70)
    rp = permutation_test_paired(np.array([0.5, 0.4, 0.6, 0.3, 0.7, 0.2]))
    assert rp.method == "exact" and rp.n == 64 and rp.p == pytest.approx(2 / 64)


def test_t_distribution_matches_scipy():
    st = pytest.importorskip("scipy.stats")
    for t, df in ((0.5, 3), (2.1, 4.7), (-3.0, 10), (1.96, 1000)):
        assert t_cdf(t, df) == pytest.approx(st.t.cdf(t, df), abs=1e-8)


# --------------------------------------------------------------------------- #
# one rule for the API, the CLI and the MCP server
# --------------------------------------------------------------------------- #
def _verdict_line(markdown: str) -> str:
    lines = markdown.splitlines()
    i = lines.index("## Verdict")
    return next(line for line in lines[i + 1:] if line.strip()).strip("* ")


def test_demo_verdict_is_identical_through_api_cli_and_mcp(capsys, tmp_path):
    out_json = tmp_path / "cli.json"
    cli.main(["--demo", "--no-stop", "--no-log", "--json", str(out_json)])
    cli_md = capsys.readouterr().out
    mcp_md = mcp_server.demo_report(stop_on_first_fail=False)
    api = run_autopsy(
        MI, cli.demo_data(), group_col="age", groups=("young", "old"), within=["sex"],
        gene_pair=("Smad3", "Col1a1"), replicate_col="mouse", pair_metric=metrics.mi_3bin,
        pos_pair=("Actb", "Gapdh"), neg_pair=("Gene0", "Gene1"), prereg=COMPOSITION,
        stop_on_first_fail=False, log_path="off",
    )
    assert _verdict_line(cli_md) == _verdict_line(mcp_md) == api.verdict
    rep = json.loads(out_json.read_text())
    assert rep["verdict"] == api.verdict
    assert {k: v["status"] for k, v in rep["fields"].items()} == {k: v.status for k, v in api.fields().items()}
    for key in ("data_sha256", "prereg_sha256", "claim_id"):
        assert rep["provenance"][key] == api.provenance[key]


def test_h5ad_verdict_hashes_and_log_are_shared_by_api_cli_and_mcp(capsys, tmp_path):
    anndata = pytest.importorskip("anndata")
    sparse = pytest.importorskip("scipy.sparse")
    sd = cli.demo_data()
    path = tmp_path / "demo.h5ad"
    ad = anndata.AnnData(X=sparse.csr_matrix(sd.X), obs=sd.obs.copy(),
                         var=pd.DataFrame(index=list(sd.var_names)))
    ad.obs_names = [f"cell_{i}" for i in range(sd.n_obs)]
    ad.write_h5ad(path)
    log = tmp_path / "runs.jsonl"
    common = ["--metric", "mi_3bin", "--gene-a", "Smad3", "--gene-b", "Col1a1", "--group-col", "age",
              "--groups", "young", "old", "--within", "sex", "--replicate-col", "mouse",
              "--estimand", "composition", "--pos-pair", "Actb", "Gapdh", "--neg-pair", "Gene0", "Gene1",
              "--no-stop"]
    cli.main(["--h5ad", str(path), *common, "--json", str(tmp_path / "cli.json"), "--log", str(log)])
    cli_verdict = _verdict_line(capsys.readouterr().out)
    mcp_md = mcp_server.autopsy_report(
        str(path), metric="mi_3bin", gene_a="Smad3", gene_b="Col1a1", group_col="age",
        groups=["young", "old"], within=["sex"], replicate_col="mouse", estimand="composition",
        pos_pair=["Actb", "Gapdh"], neg_pair=["Gene0", "Gene1"], stop_on_first_fail=False,
        json_path=str(tmp_path / "mcp.json"), log_path=str(log))
    api = run_autopsy(
        MI, anndata.read_h5ad(path), group_col="age", groups=("young", "old"), within=["sex"],
        gene_pair=("Smad3", "Col1a1"), replicate_col="mouse", pair_metric=metrics.mi_3bin,
        pos_pair=("Actb", "Gapdh"), neg_pair=("Gene0", "Gene1"), prereg=COMPOSITION,
        stop_on_first_fail=False, log_path=log)
    assert cli_verdict == _verdict_line(mcp_md) == api.verdict
    reports = [json.loads((tmp_path / f).read_text()) for f in ("cli.json", "mcp.json")]
    for rep in reports:
        for key in ("data_sha256", "prereg_sha256", "claim_id"):
            assert rep["provenance"][key] == api.provenance[key]
    # one claim, three entry points: the log counts three attempts
    assert [r["provenance"]["log"]["attempt"] for r in reports] == [1, 2]
    assert api.provenance["log"]["attempt"] == 3


def test_mcp_qc_report_uses_the_diagnostic_vocabulary(tmp_path):
    anndata = pytest.importorskip("anndata")
    sd = make_confounded()
    path = tmp_path / "confounded.h5ad"
    anndata.AnnData(X=sd.X, obs=sd.obs.copy(), var=pd.DataFrame(index=list(sd.var_names))).write_h5ad(path)
    md = mcp_server.qc_parity_report(str(path), group_col="age", groups=["young", "old"], within=["sex"])
    assert "**WARN" in md and "FAIL" not in md


def _tool_text(result) -> str:
    """Text of an MCP call_tool result (mcp 1.x returns a tuple or a list, 2.x a result object)."""
    if isinstance(result, tuple):
        result = result[0]
    content = getattr(result, "content", result)
    return "".join(getattr(c, "text", "") for c in content)


def test_mcp_server_registers_the_tools_and_shares_the_verdict(capsys):
    """The MCP front door itself (mcp 1.x FastMCP or 2.x MCPServer): tools registered with the
    v0.3 parameters, and the demo tool returns the CLI's verdict."""
    pytest.importorskip("mcp")
    import asyncio
    server = mcp_server.build_server()
    tools = {t.name: t for t in asyncio.run(server.list_tools())}
    assert set(tools) == {"autopsy_report", "qc_parity_report", "list_metrics", "demo_report"}
    t = tools["autopsy_report"]
    props = (getattr(t, "inputSchema", None) or getattr(t, "input_schema"))["properties"]
    assert {"replicate_col", "estimand", "sesoi", "prereg_path", "json_path"} <= set(props)
    md = _tool_text(asyncio.run(server.call_tool("demo_report", {"stop_on_first_fail": False})))
    cli.main(["--demo", "--no-stop", "--no-log"])
    assert _verdict_line(md) == _verdict_line(capsys.readouterr().out)
