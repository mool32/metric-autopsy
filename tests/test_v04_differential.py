"""v0.4, the owner's decisions of 2026-10-19 (validation/plans/v0.4_spec.md, change 17), each with its test:

- the differential sensitivity of a metric of structure: GATE 4's probe, the same change of the construct,
  into each group; where the groups differ in detection, a difference the groups' unequal responses would make
  at an equal coupling, on the claim's side and beyond bias_tolerance x SESOI, is subtracted at its bound before
  the verdict (cause ``differential_sensitivity``); the line of the dossier always;
- one direction for every difference of the dossier, groups[0] - groups[1], named in its header;
- the correction's line: how equal depth changed the groups' gap in genes per cell, with a warning where it grew.

The acceptance is on synthetic data (the owner: not on real data), mice nested in age whose old group is
deeper but less complex, as on preprint §4 (a block of five genes takes most of its counts, so the other genes
get fewer: more counts and fewer genes per cell):
(i) an equal coupling: not SUPPORTED, cause differential_sensitivity (SUPPORTED without the check);
(ii) an equal detection: the verdicts do not change;
(iii) a real decrease of 3 x SESOI at an equal detection: SUPPORTED;
(iv) the same decrease with a gap of detection: SUPPORTED where the effect less the bound passes.
"""
from __future__ import annotations

import functools
from functools import partial

import numpy as np
import pandas as pd
import pytest

from metric_autopsy import SimpleData, injected_signal, metrics, run_autopsy
from metric_autopsy import report as R
from metric_autopsy.core import Assessment

N_FILLER = 150
GENES = ["A", "B"] + [f"D{i}" for i in range(5)] + [f"F{i}" for i in range(N_FILLER)]
PREREG = {"estimand": "composition", "direction": "decrease", "sesoi": 0.1, "judgment_pending": False}


def deeper_less_complex(rho=(0.5, 0.5), dom=(0.1, 0.75), depth=(3000.0, 4200.0), mu_claim=-3.0, n_mice=8,
                        cells=300, var_bio=0.8, seed=0):
    """Mice nested in age (young, old). The pair A, B: log expected expression mu_claim + a shared factor whose
    share of the variance `var_bio` is the coupling rho (a correlation of log expected expression, the
    coupling probe's construct; the genes' variance stays `var_bio`); 150 fillers of lower expression. Per
    group a depth and the share `dom` of the counts that five dominant genes take: the old group deeper but
    less complex by default (more counts, fewer genes per cell, the pair detected less), as on §4."""
    rng = np.random.default_rng(seed)
    fill_mu = np.random.default_rng(123).normal(-3.0, 1.5, N_FILLER)
    blocks, rows = [], []
    for gi, age in enumerate(("young", "old")):
        for m in range(n_mice):
            load, noise = np.sqrt(rho[gi] * var_bio), np.sqrt((1 - rho[gi]) * var_bio)
            lat = rng.normal(size=cells)
            la = mu_claim + load * lat + noise * rng.normal(size=cells)
            lb = mu_claim + load * lat + noise * rng.normal(size=cells)
            lf = fill_mu[None, :] + rng.normal(0, 0.5, (cells, N_FILLER))
            rel = np.column_stack([np.exp(la), np.exp(lb), np.exp(lf)])
            rel /= rel.sum(axis=1, keepdims=True)
            D = depth[gi] * np.exp(rng.normal(0, 0.3, cells))
            rest = (1 - dom[gi]) * D[:, None] * rel
            lam = np.column_stack([rest[:, :2], dom[gi] * D[:, None] * np.full((cells, 5), 0.2), rest[:, 2:]])
            blocks.append(rng.poisson(lam).astype(float))
            rows += [dict(age=age, mouse=f"{age}{m}")] * cells
    return SimpleData(np.vstack(blocks), pd.DataFrame(rows), GENES)


SCENARIOS = {
    "gap, equal coupling": dict(rho=(0.5, 0.5)),
    "equal detection, equal coupling": dict(rho=(0.5, 0.5), dom=(0.1, 0.1), depth=(3000.0, 3000.0)),
    "equal detection, decrease": dict(rho=(0.6, 0.3), dom=(0.1, 0.1), depth=(3000.0, 3000.0)),
    "gap, decrease": dict(rho=(0.6, 0.3)),
    "smaller gap, decrease": dict(rho=(0.6, 0.3), dom=(0.1, 0.6)),
}
METRICS = {"mi_3bin": (partial(metrics.mi_3bin, gene_a="A", gene_b="B"), -3.0),
           "norm_pearson": (partial(metrics.norm_pearson, gene_a="A", gene_b="B"), -2.5)}


@functools.lru_cache(maxsize=None)
def claim(scenario: str, metric: str, check: bool = True):
    """The claim "the pair's coupling falls with age" on a scenario, with the differential sensitivity or (to
    compare) without it."""
    fn, mu = METRICS[metric]
    d = deeper_less_complex(mu_claim=mu, **SCENARIOS[scenario])
    saved = R.STRUCTURE_PROBES
    R.STRUCTURE_PROBES = saved if check else ()
    try:
        return run_autopsy(fn, d, group_col="age", groups=("young", "old"), replicate_col="mouse", gene_pair=("A", "B"),
                           signal_test=injected_signal.coupling("A", "B"), prereg=dict(PREREG), log_path="off",
                           specificity=False)
    finally:
        R.STRUCTURE_PROBES = saved


def ds(a) -> dict:
    return a.effect.detail.get("differential_sensitivity") or {}


# --------------------------------------------------------------------------- #
# the acceptance on synthetic data
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("metric", ["mi_3bin", "norm_pearson"])
def test_i_an_equal_coupling_with_a_gap_of_detection_is_not_supported(metric):
    a, without = claim("gap, equal coupling", metric), claim("gap, equal coupling", metric, check=False)
    assert without.verdict.startswith("SUPPORTED")  # the artifact the check is for
    assert a.cause == "differential_sensitivity" and a.verdict.startswith("NOT SUPPORTED — differential sensitivity")
    rec = ds(a)
    assert rec["detection"]["differs"] and rec["applies"] and rec["rich"] == "young"
    assert 0 < rec["ratio"] < 1 and rec["artifact"] > 0  # the poorer group responds less: young - old > 0
    assert rec["offset"] == pytest.approx(max(rec["artifact_ci"]))
    assert not rec["adjusted"]["passes"]
    md = a.to_markdown()
    assert "## Differential sensitivity" in md and "they differ in detection" in md


@pytest.mark.parametrize("scenario, metric", [("equal detection, equal coupling", "mi_3bin"),
                                              ("equal detection, decrease", "mi_3bin"),
                                              ("equal detection, decrease", "norm_pearson")])
def test_ii_an_equal_detection_leaves_the_verdicts_as_they_were(scenario, metric):
    a, without = claim(scenario, metric), claim(scenario, metric, check=False)
    assert (a.verdict.split(" [")[0], a.cause) == (without.verdict.split(" [")[0], without.cause)
    assert ds(a)["responses"]  # the line of the dossier is there all the same


@pytest.mark.parametrize("metric", ["mi_3bin", "norm_pearson"])
def test_iii_a_real_decrease_at_an_equal_detection_is_supported(metric):
    a = claim("equal detection, decrease", metric)
    assert a.verdict.startswith("SUPPORTED") and a.cause == "provisional"


def test_iv_a_real_decrease_with_a_gap_is_supported_where_the_effect_less_the_bound_passes():
    a = claim("smaller gap, decrease", "norm_pearson")
    rec = ds(a)
    assert rec["applies"] and rec["adjusted"]["passes"]
    assert a.verdict.startswith("SUPPORTED") and "judged on the effect less the differential sensitivity's bound" in a.verdict
    # the adjusted test is the effect's against the bound: the same estimate, shifted by it
    assert rec["adjusted"]["effect"] == pytest.approx(a.effect.detail["effect"] - rec["offset"])


@pytest.mark.parametrize("scenario", ["gap, decrease", "smaller gap, decrease"])
@pytest.mark.parametrize("metric", ["mi_3bin", "norm_pearson"])
def test_iv_with_a_gap_a_real_decrease_is_supported_exactly_where_the_effect_less_the_bound_passes(scenario, metric):
    a, without = claim(scenario, metric), claim(scenario, metric, check=False)
    assert without.verdict.startswith("SUPPORTED")  # the decrease itself is detected
    rec = ds(a)
    assert rec["detection"]["differs"] and rec["applies"]
    assert a.verdict.startswith("SUPPORTED") == rec["adjusted"]["passes"]
    if not rec["adjusted"]["passes"]:
        assert a.cause == "differential_sensitivity"


def test_iv_a_gap_as_on_section_4_removes_a_decrease_of_three_sesoi_from_mi():
    """mi_3bin responds more where the coupling is larger, so with a real decrease its ratio k_old / k_young
    also takes the slope of the curve (a known issue, conservative): the effect less the bound fails here."""
    a, without = claim("gap, decrease", "mi_3bin"), claim("gap, decrease", "mi_3bin", check=False)
    assert without.verdict.startswith("SUPPORTED")
    assert a.cause == "differential_sensitivity" and not ds(a)["adjusted"]["passes"]
    assert ds(a)["ratio"] < ds(claim("gap, equal coupling", "mi_3bin"))["ratio"]


# --------------------------------------------------------------------------- #
# the rule (report.differential_sensitivity), on records
# --------------------------------------------------------------------------- #
def _rec(artifact=0.03, ci=(0.01, 0.05), **kw):
    return dict(artifact=artifact, artifact_ci=ci, rich="young", poor="old", k_rich=0.01, k_poor=0.005, ratio=0.5,
                **kw)


OPEN, CLOSED = dict(differs=True, rows=[]), dict(differs=False, rows=[])


def _effect(e=0.05):
    return Assessment("DETECTED", "x", [], dict(effect=e))


def test_the_rule_applies_on_the_claims_side_beyond_the_tolerance_where_detection_differs():
    pre = dict(direction="decrease", bias_tolerance=0.5)
    r = R.differential_sensitivity(_rec(), _effect(), pre, 0.02, OPEN)
    assert r["applies"] and r["offset"] == pytest.approx(0.05) and r["tolerance"] == pytest.approx(0.01)
    assert not R.differential_sensitivity(_rec(), _effect(), pre, 0.02, CLOSED)["applies"]  # no gap: reported only
    assert not R.differential_sensitivity(_rec(-0.03, (-0.05, -0.01)), _effect(), pre, 0.02, OPEN)["applies"]
    assert not R.differential_sensitivity(_rec(0.004, (0.001, 0.008)), _effect(), pre, 0.02, OPEN)["applies"]
    # the claim's side by its point: a bound beyond the tolerance applies though the interval reaches 0
    assert R.differential_sensitivity(_rec(0.01, (-0.01, 0.03)), _effect(), pre, 0.02, OPEN)["offset"] == pytest.approx(0.03)
    # an increase claim: the claim's side is negative, the offset too
    r = R.differential_sensitivity(_rec(-0.03, (-0.05, -0.01)), _effect(-0.05), dict(pre, direction="increase"), 0.02, OPEN)
    assert r["applies"] and r["offset"] == pytest.approx(-0.05)
    # a two-sided claim takes the observed effect's side
    assert R.differential_sensitivity(_rec(), _effect(0.05), dict(pre, direction="two-sided"), 0.02, OPEN)["applies"]
    assert not R.differential_sensitivity(_rec(), _effect(-0.05), dict(pre, direction="two-sided"), 0.02, OPEN)["applies"]
    # without a SESOI the tolerance is 0
    assert R.differential_sensitivity(_rec(0.004, (0.001, 0.008)), _effect(), pre, None, OPEN)["applies"]


def test_an_artifact_that_cannot_be_sized_is_unsized_only_where_detection_differs():
    pre = dict(direction="decrease")
    assert R.differential_sensitivity(dict(error="x"), _effect(), pre, 0.02, OPEN)["unsized"]
    assert not R.differential_sensitivity(dict(error="x"), _effect(), pre, 0.02, CLOSED)["unsized"]
    nan = _rec(0.03, (float("nan"), float("nan")))
    assert R.differential_sensitivity(nan, _effect(), pre, 0.02, OPEN)["unsized"]


def test_the_detection_gap_reads_genes_per_cell_on_the_claims_cells():
    tech = dict(alpha=0.05, excluded=["A", "B"], rows=[
        dict(cells="all cells", covariate="genes per cell", stratum="all", difference=10.0, p=0.01),
        dict(cells="all cells", covariate="counts per cell", stratum="all", difference=-1e5, p=0.001),
        dict(cells="the domain's cells", covariate="genes per cell", stratum="all", difference=1.0, p=0.4)])
    assert R.detection_gap(tech, "all cells")["differs"]
    assert not R.detection_gap(tech, "the domain's cells")["differs"]
    assert R.detection_gap(dict(error="no pool"), "all cells")["differs"]  # untested: the premise stays open


def test_a_metric_of_level_gets_no_differential_sensitivity():
    from test_v04_dossier import mice
    d = mice(n_per_group=4, cells=120)
    a = run_autopsy(partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1"), d, group_col="age",
                    groups=("young", "old"), replicate_col="mouse", signal_test=injected_signal.module(["Smad3"]),
                    prereg=dict(PREREG, direction="decrease"), log_path="off", specificity=False)
    assert "differential" not in next(g for g in a.results if g.gate == 4).detail
    assert not ds(a)


# --------------------------------------------------------------------------- #
# one direction in every line of the dossier
# --------------------------------------------------------------------------- #
def test_every_difference_of_the_dossier_is_groups0_minus_groups1():
    a = claim("gap, equal coupling", "mi_3bin")
    e = a.effect.detail
    # the reference row's metric is the effect's own uncorrected difference, in the same direction
    assert a.reference["rows"][-1]["metric"]["difference"] == pytest.approx(e["raw_effect"])
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
# the correction's line
# --------------------------------------------------------------------------- #
def test_the_correction_line_warns_where_equal_depth_widens_the_gap_in_detection():
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
