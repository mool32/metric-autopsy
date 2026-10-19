"""v0.4, the owner's second decisions of 2026-10-19 (validation/plans/v0.4_spec.md, changes 17', 18 and 19):

- change 17': for a claim about a pair, after equal depth, the claim genes' detection is equalized in every cell
  set where the effect is estimated (the group that detects a gene more often has its counts of it thinned to
  the other's share); the effect, GATE 4's response s_m and the line of norm_pearson are on the data the effect
  is judged on, and GATE 4 runs on one draw of the correction for every claim; the dossier carries the
  detection before and after and the effect at equal depth alone; no new cause;
- change 18: one direction for every difference of the dossier, groups[0] - groups[1], named in its header;
- change 19: the correction's line, with a warning where equal depth widened the groups' gap in genes per cell.

Change 17 (the differential sensitivity) was reverted by the same decisions; its acceptance on synthetic data
is in `validation/exploratory/differential_sensitivity.log`, the evidence of its rejection. The acceptance of
17' (20 seeds per scenario) is `validation/exploratory/detection_equalization.py`; the tests here are its
mechanics."""
from __future__ import annotations

import functools
from functools import partial

import numpy as np
import pandas as pd
import pytest

from detection_data import METRICS, PREREG, SCENARIOS, deeper_less_complex
from metric_autopsy import SimpleData, injected_signal, metrics, run_autopsy
from metric_autopsy.equalize import detection_keep_probability, equalize_detection


@functools.lru_cache(maxsize=None)
def claim(scenario: str, metric: str):
    fn, mu = METRICS[metric]
    d = deeper_less_complex(mu_claim=mu, **SCENARIOS[scenario])
    return run_autopsy(fn, d, group_col="age", groups=("young", "old"), replicate_col="mouse", gene_pair=("A", "B"),
                       signal_test=injected_signal.coupling("A", "B"), prereg=dict(PREREG), log_path="off",
                       specificity=False)


# --------------------------------------------------------------------------- #
# the equalization
# --------------------------------------------------------------------------- #
def test_the_keep_probability_reaches_the_other_groups_share():
    rng = np.random.default_rng(1)
    x = rng.poisson(rng.gamma(1.0, 1.5, 20000))
    now = (x > 0).mean()
    for target in (0.9 * now, 0.5 * now, 0.05):
        p = detection_keep_probability(x, target)
        realized = np.mean([(rng.binomial(x, p) > 0).mean() for _ in range(10)])
        assert realized == pytest.approx(target, abs=0.004)
    assert detection_keep_probability(x, now) == 1.0 and detection_keep_probability(x, 1.0) == 1.0
    assert detection_keep_probability(x, 0.0) == 0.0


def test_equalize_detection_thins_the_more_detecting_group_in_each_cell_set():
    rng = np.random.default_rng(2)
    n = 4000
    obs = pd.DataFrame(dict(age=np.repeat(["young", "old"], n), sex=np.tile(np.repeat(["f", "m"], n // 2), 2)))
    lam = np.where(obs["age"] == "young", 2.0, 0.6)[:, None] * np.ones((1, 3))
    lam[(obs["sex"] == "m").to_numpy(), 1] = 0.8  # in males gene B is detected alike in both groups
    X = rng.poisson(lam).astype(float)
    d = SimpleData(X.copy(), obs, ["A", "B", "C"])
    recs = equalize_detection(d, ["A", "B"], "age", ("young", "old"), within=["sex"], rng=np.random.default_rng(3))
    by = {(r["cells"]["sex"], r["gene"]): r for r in recs}
    assert set(by) == {("f", "A"), ("f", "B"), ("m", "A"), ("m", "B")}
    for key in (("f", "A"), ("f", "B"), ("m", "A")):
        r = by[key]
        assert r["thinned_group"] == "young" and r["keep_p"] < 1
        assert r["after"]["young"] == pytest.approx(r["before"]["old"], abs=0.03)
        assert r["after"]["old"] == r["before"]["old"]  # the other group is left as it was
    assert by[("m", "B")]["keep_p"] > 0.9  # detected alike: thinned slightly or not at all
    np.testing.assert_array_equal(d.X[:, 2], X[:, 2])  # a gene outside the claim is left as it was
    old = (obs["age"] == "old").to_numpy()
    np.testing.assert_array_equal(d.X[old, 0], X[old, 0])


def test_without_counts_nothing_is_equalized():
    d = deeper_less_complex(**SCENARIOS["gap, equal coupling"], n_mice=4, cells=100)
    lognorm = SimpleData(np.log1p(d.X / d.X.sum(1, keepdims=True) * 1e4), d.obs.copy(), list(d.var_names))
    a = run_autopsy(partial(metrics.norm_pearson, gene_a="A", gene_b="B"), lognorm, group_col="age",
                    groups=("young", "old"), replicate_col="mouse", gene_pair=("A", "B"),
                    signal_test=injected_signal.coupling("A", "B"), prereg=dict(PREREG), log_path="off",
                    specificity=False)
    g4 = next(r for r in a.results if r.gate == 4)
    assert "as they are" in g4.detail["data"]
    assert not (a.effect.detail.get("detection") or {}).get("records")


# --------------------------------------------------------------------------- #
# the verdict's data: the effect, GATE 4 and the reference line
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("metric", ["mi_3bin", "norm_pearson"])
def test_a_pair_is_judged_at_equal_depth_and_detection(metric):
    a = claim("gap, equal coupling", metric)
    e = a.effect.detail
    recs = e["detection"]["records"]
    assert {r["gene"] for r in recs} == {"A", "B"}
    for r in recs:
        assert r["raw"]["young"] > r["raw"]["old"] + 0.05  # the gap of §4
        assert r["thinned_group"] == "young" and r["keep_p"] < 0.5
        assert abs(r["after"]["young"] - r["after"]["old"]) < 0.03
    # the gap's artifact at equal depth alone is larger than what is left once the detection is equal
    assert e["equal_depth_effect"]["effect"] > abs(e["effect"])
    assert e["corrected_for"] == "depth and the claim genes' detection"
    g4 = next(r for r in a.results if r.gate == 4)
    assert g4.detail["data"] == "equal depth, the claim genes' detection equalized"
    assert {r["gene"] for r in g4.detail["equalization"]["detection"]} == {"A", "B"}
    # the reference line is on the effect's own data and draws: its metric is the effect
    ref = a.reference
    assert ref["data"].startswith("the effect's: at equal depth, the pair's detection equalized")
    assert ref["rows"][-1]["metric"]["difference"] == pytest.approx(e["effect"])
    md = a.to_markdown()
    assert "## The pair's detection, equalized after equal depth" in md and "the effect at equal depth alone" in md
    assert "Measured on the data the effect is judged on: equal depth, the claim genes' detection equalized" in md


def test_at_an_equal_detection_the_equalization_barely_moves_the_effect():
    a = claim("equal detection, decrease", "norm_pearson")
    e = a.effect.detail
    assert all(r["keep_p"] > 0.8 for r in e["detection"]["records"])
    assert e["effect"] == pytest.approx(e["equal_depth_effect"]["effect"], rel=0.1)
    assert a.verdict.startswith("SUPPORTED")


def test_gate4_runs_at_equal_depth_for_a_metric_of_level_without_equalizing_detection():
    d = deeper_less_complex(**SCENARIOS["gap, equal coupling"], n_mice=4, cells=120)
    a = run_autopsy(partial(metrics.mi_3bin, gene_a="A", gene_b="B"), d, group_col="age", groups=("young", "old"),
                    replicate_col="mouse", signal_test=injected_signal.module(["A"]), prereg=dict(PREREG),
                    log_path="off", specificity=False)
    g4 = next(r for r in a.results if r.gate == 4)
    assert g4.detail["data"] == "equal depth" and g4.detail["equalization"]["detection"] is None
    assert "detection" not in a.effect.detail


# --------------------------------------------------------------------------- #
# change 18: one direction in every line of the dossier
# --------------------------------------------------------------------------- #
def test_every_difference_of_the_dossier_is_groups0_minus_groups1():
    a = claim("gap, equal coupling", "mi_3bin")
    e = a.effect.detail
    assert a.reference["rows"][-1]["metric"]["difference"] == pytest.approx(e["effect"])
    tech = {(r["cells"], r["covariate"]): r for r in a.technical["rows"]}
    gpc = tech[("all cells", "genes per cell")]
    assert gpc["difference"] == pytest.approx(gpc["first"] - gpc["second"]) and gpc["difference"] > 0
    md = a.to_markdown()
    assert "Every difference below is young − old (groups[0] − groups[1]" in md
    assert "young − old |" in md  # the reference table's columns
    assert a.to_dict()["differences"] == "young - old"


def test_the_specificity_line_gives_t_as_groups0_minus_groups1():
    from metric_autopsy.specificity import describe
    line = dict(applicable=True, statistic="the t", data="uncorrected", side="in the claimed direction (increase "
                "from young to old)", groups=["young", "old"], sign=-1.0, t_claim=-2.5, t_claim_signed=2.5,
                in_direction=dict(k=1, n=10, p=0.2, share_ci=(0, 0.4), outcome="undecided"),
                absolute=dict(k=1, n=10, p=0.2, outcome="undecided"),
                nearest=[dict(genes=["X", "Y"], t=-1.0, t_signed=1.0, at_least_as_far=False)],
                shift=dict(median_t=0.5, ci=(0.1, 0.9), ci_level=0.95, n_relabelings=70, p_one_sided=0.02))
    out = describe(line)
    assert "young − old: t = -2.5" in out[0] and "a negative t" in out[0]
    assert "X–Y -1" in out[2] and "median t -0.5" in out[3] and "[-0.9, -0.1]" in out[3]


# --------------------------------------------------------------------------- #
# change 19: the correction's line
# --------------------------------------------------------------------------- #
def test_the_correction_line_warns_where_equal_depth_widens_the_gap_in_detection():
    from metric_autopsy import report as R
    a = claim("gap, equal coupling", "mi_3bin")
    cap = a.capture
    row = cap["rows"][0]
    assert row["cells"] == "all cells" and row["genes_raw"] > 0
    assert row["widened"] == (abs(row["genes_equal_depth"]) > abs(row["genes_raw"])) and cap["widened"]
    md = a.to_markdown()
    assert "## The correction and the groups' detection" in md and "WARNING: the correction widened" in md
    lines = R.describe_capture(dict(groups=["young", "old"], rows=[
        dict(cells="all cells", genes_raw=0.24, genes_equal_depth=0.29, depth_raw=-0.54, depth_equal_depth=-0.06,
             widened=True),
        dict(cells="the domain's cells", genes_raw=0.17, genes_equal_depth=0.12, depth_raw=-0.6, depth_equal_depth=-0.1,
             widened=False)]))
    assert "+24.0% uncorrected, +29.0% at equal depth" in lines[0] and "WARNING" in lines[0]
    assert "WARNING" not in lines[1]
