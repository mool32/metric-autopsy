"""v0.4, the owner's second decisions of 2026-10-18 (validation/plans/v0.4_spec.md), each with its regression test:

- norm_pearson gives NaN, not 0.0, where it measures nothing: fewer than MIN_CODETECTED (20) cells carry both
  genes, or one gene does not vary among them; the replicates where the metric gives no value are left out of
  both tests and listed, and where too few are left the comparison is UNIDENTIFIABLE (cause ``unmeasured``);
- the full dossier: after a FAIL of metric validity the effect, the domain's, the mixture's and every report
  line are estimated, marked as not deciding, and the verdict follows the first FAIL;
- change 12: GATE 0 blocks only a bias it cannot show to be equal between the groups (the metric on GATE 0's
  null, estimated between the groups as the effect is);
- the domain of a claim about a pair: the types where both genes are detected in each group; the types where
  they reach the threshold in one group only are a line of the dossier;
- the groups' technical covariates per replicate (counts and genes per cell), with permutation p.
"""
from __future__ import annotations

from functools import partial

import numpy as np
import pandas as pd
import pytest

from metric_autopsy import SimpleData, gates, injected_signal, metrics, run_autopsy
from metric_autopsy import gates as _gates
from metric_autopsy.core import as_dense
from metric_autopsy.effect import estimate_effect
from metric_autopsy.report import CAUSES, technical_difference

RESOLVED = {"estimand": "composition", "judgment_pending": False}
NPR = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
MI = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")


def mice(n_per_group=6, cells=200, n_filler=30, claim_keep=(1.0, 1.0), efficiency=(1.0, 1.0),
         coupling=(1.5, 1.5), sexes=("female",), absent_in=(), seed=0):
    """Mice nested in age (young, old). Per group: a capture efficiency of every gene (`efficiency`),
    and a share of the claim genes' counts kept (`claim_keep`: a dropout of the claim genes alone,
    which equal depth does not remove). Col1a1 has no count in the mice named in `absent_in`."""
    rng = np.random.default_rng(seed)
    rows, blocks = [], []
    for gi, age in enumerate(("young", "old")):
        for sex in sexes:
            for m in range(n_per_group):
                name = f"{age}{m}" if len(sexes) == 1 else f"{age}-{sex}-{m}"
                lat = rng.normal(0, 1, cells)
                pair = np.exp(1.2 + coupling[gi] * 0.6 * lat[:, None] + rng.normal(0, 0.3, (cells, 2)))
                filler = np.exp(0.4 + rng.normal(0, 0.5, (cells, n_filler)))
                X = rng.poisson(np.column_stack([pair, filler]) * efficiency[gi]).astype(float)
                if claim_keep[gi] < 1.0:
                    X[:, :2] = rng.binomial(X[:, :2].astype(np.int64), claim_keep[gi])
                if name in absent_in:
                    X[:, 1] = 0.0
                blocks.append(X)
                rows += [dict(age=age, sex=sex, mouse=name)] * cells
    names = ["Smad3", "Col1a1"] + [f"G{i}" for i in range(n_filler)]
    return SimpleData(np.vstack(blocks), pd.DataFrame(rows), names)


def coupling_and_detection(data):
    """A metric of structure (the pair's normalized Pearson) plus the share of cells with a count of Smad3:
    the second term is no biology of the pair, and dropout moves it on GATE 0's null too."""
    X = as_dense(data.X)
    return NPR(data) + float(np.mean(X[:, list(data.var_names).index("Smad3")] > 0))


def big(**kw):
    """Ten mice per group of 600 cells: enough for the groups' technical difference to be shown within 0.1."""
    return mice(10, cells=600, **kw)


def _only_dropout(data, kind, rng, protect=frozenset()):
    """GATE 0 with every perturbation but the extra dropout replaced by the identity."""
    if kind == "extra_dropout":
        return _ORIG_PERTURB(data, kind, rng, protect)
    return _gates.SimpleData(data.X.copy(), data.obs, data.var_names), {}


_ORIG_PERTURB = _gates._perturb


# --------------------------------------------------------------------------- #
# norm_pearson: NaN where it measures nothing (the owner's item 4)
# --------------------------------------------------------------------------- #
def test_norm_pearson_gives_no_value_where_it_measures_nothing():
    rng = np.random.default_rng(0)
    n = 40
    X = np.zeros((n, 4))
    X[:, 2:] = rng.poisson(6, (n, 2)) + 1
    X[:19, 0] = rng.poisson(3, 19) + 1
    X[:19, 1] = rng.poisson(3, 19) + 1
    d = SimpleData(X, pd.DataFrame({"i": range(n)}), ["a", "b", "c", "e"])
    assert metrics.MIN_CODETECTED == 20
    assert np.isnan(metrics.norm_pearson(d, gene_a="a", gene_b="b"))  # 19 cells carry both genes
    X[19, :2] = (2, 5)
    v = metrics.norm_pearson(d, gene_a="a", gene_b="b")  # 20
    assert np.isfinite(v) and -1 <= v <= 1
    assert np.isnan(metrics.norm_pearson(d, gene_a="a", gene_b="b", min_codetected=21))
    two = SimpleData(X[17:], pd.DataFrame({"i": range(n - 17)}), ["a", "b", "c", "e"])  # 3 cells carry both
    assert np.isfinite(metrics.norm_pearson(two, gene_a="a", gene_b="b", min_codetected=0))
    assert np.isnan(metrics.norm_pearson(SimpleData(X[18:], pd.DataFrame({"i": range(n - 18)}), ["a", "b", "c", "e"]),
                                         gene_a="a", gene_b="b", min_codetected=0))  # never fewer than 3
    flat = SimpleData(np.full((30, 4), 2.0), pd.DataFrame({"i": range(30)}), ["a", "b", "c", "e"])
    assert np.isnan(metrics.norm_pearson(flat, gene_a="a", gene_b="b"))  # no variance among them


def test_replicates_the_metric_gives_no_value_on_are_left_out_and_listed():
    d = mice(6, absent_in=("young0",))
    eff, _ = estimate_effect(NPR, d, group_col="age", groups=("young", "old"), replicate_col="mouse",
                             estimand="composition", seed=0)
    assert [u["replicate"] for u in eff.detail["unmeasured"]] == ["young0"]
    assert eff.detail["unmeasured"][0]["raw"] is None and eff.detail["n_replicates"] == {"young": 5, "old": 6}
    assert np.isfinite(eff.detail["effect"]) and eff.status != "NOT_ESTIMABLE"
    row = next(r for r in eff.detail["replicate_values"] if r["replicate"] == "young0")
    assert row["raw"] is None and row["corrected"] is None  # never a value of 0


def test_too_few_measured_replicates_make_the_comparison_unidentifiable():
    assert CAUSES.index("unidentifiable") < CAUSES.index("unmeasured") < CAUSES.index("effect_not_evaluated")
    d = mice(4, absent_in=("young0", "young1"))
    a = run_autopsy(NPR, d, group_col="age", groups=("young", "old"), replicate_col="mouse",
                    prereg={**RESOLVED, "direction": "decrease"}, log_path="off")
    assert a.effect.status == "NOT_ESTIMABLE" and a.effect.detail["unmeasured_insufficient"]
    assert a.cause == "unmeasured", a.verdict
    assert a.verdict.startswith("UNIDENTIFIABLE — the metric gives no value (not measured) on 2 of the 8 replicates")
    assert [u["replicate"] for u in a.effect.detail["unmeasured"]] == ["young0", "young1"]
    md = a.to_markdown()  # listed in the report, not only in its JSON
    assert "## Left out: the metric gives no value (not measured)" in md
    assert "- replicate young0: uncorrected no value, at equal depth no value" in md


def test_units_of_the_domain_the_metric_gives_no_value_on_are_left_out():
    from test_v04_domain import typed_mice
    d = typed_mice(n_mice=5, seed=1)
    X = np.array(d.X)
    young0_a = (np.asarray(d.obs["mouse"]) == "young0") & (np.asarray(d.obs["cell_type"]) == "A")
    X[young0_a, 1] = 0.0
    X[np.where(young0_a)[0][:5], 1] = 1.0  # Col1a1 in 5 of the unit's cells: fewer than 20 carry both genes
    d = SimpleData(X, d.obs, d.var_names)
    eff, _ = estimate_effect(NPR, d, group_col="age", groups=("young", "old"), replicate_col="mouse",
                             estimand="composition", seed=0, domain=_domain(d))
    dom = eff.detail["domain"]
    assert [(u["replicate"], u["type"]) for u in dom["unmeasured"] if u["in_domain"]] == [("young0", "A")]
    assert dom["n_before"] == 5 and dom["n_after"] == 4 and np.isfinite(eff.detail["effect"])
    # outside the domain B carries no Col1a1: its units give no value either, and are listed
    assert {u["type"] for u in dom["unmeasured"] if not u["in_domain"]} == {"B"}


def _domain(d):
    from metric_autopsy.domain import claim_domain
    return claim_domain(d, type_col="cell_type", group_col="age", groups=("young", "old"), replicate_col="mouse",
                        genes=("Smad3", "Col1a1"), kind="pair")


# --------------------------------------------------------------------------- #
# the full dossier after the first FAIL (the owner's item 1)
# --------------------------------------------------------------------------- #
def test_after_a_validity_fail_the_domain_the_mixture_and_the_lines_are_estimated_as_the_dossier():
    from test_v04_domain import typed_mice
    d = typed_mice(col_in_b=(0.05, 0.0), seed=3)
    a = run_autopsy(lambda data: -MI(data), d, group_col="age", groups=("young", "old"), replicate_col="mouse",
                    gene_pair=("Smad3", "Col1a1"), composition=["cell_type"],
                    signal_test=injected_signal.coupling("Smad3", "Col1a1"),
                    prereg={**RESOLVED, "direction": "decrease", "sesoi": 0.1}, log_path="off")
    assert a.metric_validity.status == "FAIL" and a.cause == "metric_invalid_gate4"
    assert a.verdict.startswith("NOT SUPPORTED — metric invalid")
    ef = a.effect
    assert "DOSSIER_ONLY" in ef.flags and ef.detail["dossier_only"] is True
    assert ef.detail["domain"]["types"] == ["A"] and ef.detail["mixture"]["raw"]["effect"] is not None
    assert a.composition["columns"]["cell_type"]["tvd"] is not None and a.technical["rows"]
    md = a.to_markdown()
    assert "## The domain (the verdict is judged within it)" in md and "(the dossier: estimated after" in md
    assert "## Technical covariates per replicate (reported, not a gate)" in md
    # the verdict is the same as the first FAIL's alone: the dossier never decides
    from metric_autopsy.report import decide_cause
    ef.detail.pop("domain")
    assert decide_cause(a)[1] == "metric_invalid_gate4"


# --------------------------------------------------------------------------- #
# change 12: GATE 0 blocks only a bias not shown equal between the groups (the owner's item 2)
# --------------------------------------------------------------------------- #
def test_the_groups_technical_difference_is_the_metric_on_gate0s_null_between_the_groups():
    kw = dict(group_col="age", groups=("young", "old"), replicate_col="mouse", estimand="composition", margin=0.1,
              seed=0)
    equal = technical_difference(coupling_and_detection, big(), ("Smad3", "Col1a1"), **kw)
    assert equal["equivalent"] and abs(equal["effect"]) < 0.05 and equal["draws"] == 3
    lo, hi = equal["ci_2alpha"]
    assert -0.1 < lo <= equal["effect"] <= hi < 0.1
    # a dropout of the claim genes alone in the old: equal depth does not remove it, nor does the null
    unequal = technical_difference(coupling_and_detection, big(claim_keep=(1.0, 0.3)), ("Smad3", "Col1a1"), **kw)
    assert not unequal["equivalent"] and unequal["ci_2alpha"][0] > 0.1  # young − old


def test_a_bias_equal_between_the_groups_is_reported_and_one_that_differs_blocks(monkeypatch):
    monkeypatch.setattr(_gates, "_perturb", _only_dropout)
    claim = {**RESOLVED, "direction": "two-sided", "sesoi": 0.2, "sesoi_scale": "observed"}  # tolerance 0.1
    kw = dict(group_col="age", groups=("young", "old"), replicate_col="mouse", gene_pair=("Smad3", "Col1a1"),
              prereg=claim, log_path="off", specificity=False)
    same = run_autopsy(coupling_and_detection, big(), **kw)
    g0 = next(r for r in same.results if r.gate == 0)
    r = g0.detail["responses"]["extra_dropout"]
    assert r["classification"] == "bias" and r["bias_kind"] == "null moves"
    assert r["bias_ci_low"] > r["bias_tolerance"] == pytest.approx(0.1)
    assert r["handling"] == "equal between groups" and g0.status.value == "WARN"
    assert g0.detail["between_groups"]["equivalent"] and "shown within it" in g0.message
    assert same.metric_validity.status != "FAIL"
    differ = run_autopsy(coupling_and_detection, big(claim_keep=(1.0, 0.3)), **kw)
    g0d = next(r for r in differ.results if r.gate == 0)
    assert g0d.status.value == "FAIL" and g0d.detail["responses"]["extra_dropout"]["handling"] == "blocking"
    assert not g0d.detail["between_groups"]["equivalent"] and "not shown within" in g0d.message
    assert differ.cause == "metric_invalid_gate0"
    # the dossier still estimates the effect after GATE 0's FAIL
    assert "DOSSIER_ONLY" in differ.effect.flags and differ.effect.detail["effect"] is not None


def test_without_replicates_a_blocking_bias_still_blocks(monkeypatch):
    monkeypatch.setattr(_gates, "_perturb", _only_dropout)
    a = run_autopsy(coupling_and_detection, big(), group_col="age", groups=("young", "old"),
                    gene_pair=("Smad3", "Col1a1"),
                    prereg={**RESOLVED, "direction": "two-sided", "sesoi": 0.2, "sesoi_scale": "observed"},
                    log_path="off", specificity=False)
    g0 = next(r for r in a.results if r.gate == 0)
    assert g0.status.value == "FAIL" and "between_groups" not in g0.detail


def test_gate0_alone_keeps_its_result_and_change_12_is_applied_by_the_autopsy(monkeypatch):
    """`gate0_independence` called alone has no groups to compare: change 12 is applied by `run_autopsy`
    (``gates.gate0_between_groups``), so a direct call still FAILs a blocking bias."""
    monkeypatch.setattr(_gates, "_perturb", _only_dropout)
    g0 = gates.gate0_independence(coupling_and_detection, big(), sesoi=0.2)
    assert g0.status.value == "FAIL" and "between_groups" not in g0.detail
    between = dict(equivalent=True, margin=0.1, alpha=0.05, effect=0.001, ci_2alpha=[-0.01, 0.012],
                   estimate_over="the replicates")
    g1 = gates.gate0_between_groups(g0, between)
    assert g1.status.value == "WARN" and g0.status.value == "FAIL"  # the original is not changed
    assert g1.detail["bias_handling"]["extra_dropout"] == "equal between groups"
    assert gates.gate0_between_groups(g0, dict(between, equivalent=False)).status.value == "FAIL"


# --------------------------------------------------------------------------- #
# the domain of a pair: both genes in each group (the owner's item 3)
# --------------------------------------------------------------------------- #
def test_a_type_where_the_pair_appears_in_one_group_is_outside_the_domain_and_listed():
    from metric_autopsy.domain import claim_domain, describe
    from test_v04_domain import typed_mice
    induced = typed_mice(col_in_a=(0.5, 0.03))  # Col1a1 in about 40% of the young A cells, 3% of the old
    dom = claim_domain(induced, type_col="cell_type", group_col="age", groups=("young", "old"),
                       replicate_col="mouse", genes=("Smad3", "Col1a1"), kind="pair")
    a = next(r for r in dom["strata"][0]["types"] if r["type"] == "A")
    assert a["detection"]["old"]["Col1a1"] < 0.10 <= a["detection"]["young"]["Col1a1"]
    assert not a["in_domain"] and a["detected_in"] == ["young"] and dom["domain"] == {"": []}
    assert [r["type"] for r in dom["one_group"]] == ["A"] and dom["rule"]["groups_needed"] == "each"
    assert a["index"]["young"] > a["index"]["old"] > 0
    assert any("reaches the threshold in young only" in ln for ln in describe(dom))
    assert dom["rule"]["text"].startswith("types where, in each group, both genes are detected in at least 10%")
    # a claim about a set keeps the type: a change of expression may be the expression's appearance
    as_set = claim_domain(induced, type_col="cell_type", group_col="age", groups=("young", "old"),
                          replicate_col="mouse", genes=("Col1a1",), kind="set")
    assert "A" in as_set["domain"][""] and as_set["one_group"] == []


def test_a_claim_about_a_pair_that_appears_is_not_a_claim_about_coupling():
    from test_v04_domain import _typed, typed_mice
    a = _typed(typed_mice(col_in_a=(0.5, 0.03), col_in_b=(0.0, 0.0)))
    assert a.cause == "empty_domain" and a.verdict.startswith("INCONCLUSIVE — the claim's domain is empty")
    assert "reaches the threshold in young only" in a.to_markdown()


# --------------------------------------------------------------------------- #
# the groups' technical covariates per replicate (the owner's item 9)
# --------------------------------------------------------------------------- #
def test_the_technical_covariates_are_compared_per_replicate_with_their_permutation_p():
    from metric_autopsy.composition import describe_technical, technical_covariates
    d = mice(4, efficiency=(1.0, 0.5), sexes=("female", "male"))
    tech = technical_covariates(d, "age", ("young", "old"), "mouse", ["sex"])
    rows = {(r["covariate"], r["stratum"]): r for r in tech["rows"]}
    assert set(rows) == {(c, s) for c in ("counts per cell", "genes per cell") for s in ("female", "male", "combined")}
    for s in ("female", "male"):
        r = rows[("counts per cell", s)]
        assert r["difference"] > 0 and r["method"] == "exact" and r["n_relabelings"] == 70
        assert r["p"] == pytest.approx(2 / 70)  # every young replicate deeper than every old one
        assert r["replicates"] == {"young": 4, "old": 4}
    comb = rows[("counts per cell", "combined")]
    assert comb["n_relabelings"] == 70 * 70 and comb["p"] < 0.05
    assert {"cells": "all cells", "covariate": "counts per cell", "stratum": "combined"} in tech["differs"]
    assert rows[("counts per cell", "female")]["first"] == pytest.approx(np.mean([
        np.median(as_dense(d[np.asarray(d.obs["mouse"]) == f"young-female-{m}"].X).sum(axis=1)) for m in range(4)]))
    lines = describe_technical(tech)
    assert any(ln.startswith("all cells, genes per cell, male: young") for ln in lines)
    same = technical_covariates(mice(4, seed=5), "age", ("young", "old"), "mouse")
    assert all(r["stratum"] == "all" for r in same["rows"]) and same["source"] == "data"


def test_the_technical_covariates_use_every_gene_of_the_pool_and_never_decide():
    from metric_autopsy.composition import technical_covariates
    d = mice(4, efficiency=(1.0, 0.5))
    small = SimpleData(as_dense(d.X)[:, :6], d.obs, list(d.var_names)[:6])  # the data's cut of the genes
    full = technical_covariates(small, "age", ("young", "old"), "mouse", pool=d)
    cut = technical_covariates(small, "age", ("young", "old"), "mouse")
    assert full["source"] == "pool" and cut["source"] == "data"
    f = next(r for r in full["rows"] if r["covariate"] == "genes per cell")
    c = next(r for r in cut["rows"] if r["covariate"] == "genes per cell")
    assert f["first"] > c["first"] and c["first"] <= 6
    a = run_autopsy(NPR, d, group_col="age", groups=("young", "old"), replicate_col="mouse",
                    prereg={**RESOLVED, "direction": "decrease"}, log_path="off")
    b = run_autopsy(NPR, d, group_col="age", groups=("young", "old"), replicate_col="mouse",
                    prereg={**RESOLVED, "direction": "decrease"}, log_path="off", specificity_pool=d)
    assert a.technical["rows"] and a.verdict == b.verdict and a.cause == b.cause
    assert a.to_dict()["technical"]["rows"][0]["covariate"] == "counts per cell"
