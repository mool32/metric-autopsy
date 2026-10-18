"""v0.4, the owner's decisions of 2026-10-18 (validation/plans/v0.4_spec.md), each with its regression test:

- changes 3 and 4: the SESOI on the construct's scale, translated by GATE 4's probe into the metric's
  units (s_m, its interval over the replicates); GATE 4 decides only whether the metric is blind; s_m in
  GATE 0's tolerance, the TOST and EFFECT BELOW SESOI; the design's floor beside it, and in its place
  without a SESOI;
- the domain: with cell-type labels the claim is judged within the types where its genes are detected,
  the types weighted equally, GATE 4 on the domain's cells, the whole mixture reported beside;
- change 13: the check of specificity, a line of the report, not a condition of SUPPORTED; GATE 5's
  negative control reported, not failed;
- GATE 0's known issue: its rate of false FAIL over seeds against a target fixed in advance.
"""
from __future__ import annotations

from functools import partial

import numpy as np
import pandas as pd
import pytest

from metric_autopsy import SimpleData, gate0_independence, gates, injected_signal, metrics, run_autopsy
from metric_autopsy.core import as_dense

RESOLVED = {"estimand": "composition", "judgment_pending": False}
MI = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")
NPR = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")


def typed_mice(n_mice=4, cells=(120, 120), coupling=(1.5, 1.5), col_in_b=(0.0, 0.0), col_in_a=None, seed=0,
               n_filler=40, filler_sd=0.5):
    """Two cell types in every mouse, young and old. Type A expresses the pair, Smad3 and Col1a1, with a
    per-group coupling; type B expresses Smad3 and holds only `col_in_b` (per group: young, old) of
    Col1a1 per cell, uncoupled. `col_in_a`, if given, replaces Col1a1's mean in A per group (detection
    then differs between the groups)."""
    rng = np.random.default_rng(seed)
    rows, blocks = [], []
    for gi, age in enumerate(("young", "old")):
        for m in range(n_mice):
            for t, n in zip(("A", "B"), cells):
                lat = rng.normal(0, 1, n)
                smad = np.exp(1.2 + (coupling[gi] * 0.6 * lat if t == "A" else 0) + rng.normal(0, 0.3, n))
                if t == "A":
                    base = 1.2 if col_in_a is None else np.log(col_in_a[gi])
                    col = np.exp(base + coupling[gi] * 0.6 * lat + rng.normal(0, 0.3, n))
                else:
                    col = np.full(n, col_in_b[gi])
                filler = np.exp(0.4 + rng.normal(0, filler_sd, (n, n_filler)))
                X = rng.poisson(np.column_stack([smad, col, filler])).astype(float)
                blocks.append(X)
                rows += [dict(age=age, mouse=f"{age}{m}", cell_type=t, sex="female")] * n
    names = ["Smad3", "Col1a1"] + [f"G{i}" for i in range(n_filler)]
    return SimpleData(np.vstack(blocks), pd.DataFrame(rows), names)


def _typed(d, **kw):
    kw.setdefault("prereg", {**RESOLVED, "direction": "decrease", "sesoi": 0.1})
    kw.setdefault("signal_test", injected_signal.coupling("Smad3", "Col1a1"))
    return run_autopsy(MI, d, group_col="age", groups=("young", "old"), replicate_col="mouse",
                       gene_pair=("Smad3", "Col1a1"), composition=["cell_type"], log_path="off", **kw)


def _g4(a):
    return next(r for r in a.results if r.gate == 4)


# --------------------------------------------------------------------------- #
# GATE 4 decides only whether the metric is blind (change 3)
# --------------------------------------------------------------------------- #
def test_gate4_shows_a_response_however_small_and_blindness_only_at_or_below_zero():
    iv = gates.response_interval([0.010, 0.012, 0.011, 0.013, 0.009, 0.011])
    assert gates.judge_response(iv) == "PASS" and gates.judge_response(iv, 0.05) == "PASS"  # delta_min unused
    assert gates.judge_response(gates.response_interval([-0.02, 0.03, 0.01, -0.01])) == "UNTESTED"
    assert gates.judge_response(gates.response_interval([-0.020, -0.012, -0.015, -0.011])) == "FAIL"
    assert gates.judge_response(gates.response_interval([0.0, 0.0, 0.0])) == "FAIL"  # a metric that cannot move


def test_the_module_probe_is_a_change_of_the_sesoi_on_the_log_scale():
    from test_v04 import mean_lognorm_xist, sex_marker
    a = run_autopsy(mean_lognorm_xist, sex_marker(), group_col="sex", groups=("female", "male"), replicate_col="mouse",
                    signal_test=injected_signal.module(["Xist"]),
                    prereg={**RESOLVED, "direction": "decrease", "sesoi": 0.5}, log_path="off")
    g = _g4(a)
    assert g.detail["dose"]["fold"] == pytest.approx(np.exp(0.5)) and g.detail["construct_change"] == 0.5
    s_m = g.detail["sesoi_metric"]
    assert s_m > 0 and a.effect.detail["effective_sesoi"] == pytest.approx(s_m)
    assert "s_m" in a.design_adequacy.detail["power"]["sesoi_label"]
    # the dose-response curve beside the judgment: 0.5x, 1x and 2x the SESOI, rising with the dose
    curve = g.detail["dose_response"]
    assert [p["construct_change"] for p in curve] == [0.25, 0.5, 1.0]
    sign = -1.0  # a decrease
    assert sign * curve[0]["response"] < sign * curve[1]["response"] < sign * curve[2]["response"]


def test_the_coupling_dose_raises_the_correlation_of_log_expected_expression_by_the_sesoi():
    from metric_autopsy.injected_signal import coupling_strength_for, log_sigmoid_variance
    for S, va, vb in ((0.1, 0.3, 0.5), (0.2, 1.0, 1.0), (0.05, 0.02, 2.0)):
        info = coupling_strength_for(S, va, vb)
        v = log_sigmoid_variance(info["strength"])
        assert v / np.sqrt((va + v) * (vb + v)) == pytest.approx(S, rel=1e-6)
    assert coupling_strength_for(0.1, 0.0, 0.5)["floored"] == ["a"]
    with pytest.raises(ValueError):
        coupling_strength_for(1.0, 0.3, 0.3)


def test_a_construct_sesoi_without_a_probe_is_not_translated_and_the_floor_stands_in():
    from test_v03_verdict import nested_mice
    a = run_autopsy(metrics.norm_pearson.__wrapped__ if hasattr(metrics.norm_pearson, "__wrapped__") else NPR,
                    nested_mice(4, c_young=1.5, c_old=1.5), group_col="age", groups=("young", "old"),
                    replicate_col="mouse", prereg={**RESOLVED, "direction": "two-sided", "sesoi": 0.1}, log_path="off")
    assert "SESOI_NOT_TRANSLATED" in a.design_adequacy.flags
    assert a.effect.detail["sesoi_source"] == "floor"
    assert a.effect.detail["effective_sesoi"] == pytest.approx(a.effect.detail["floor"])


def test_without_a_sesoi_absence_is_judged_at_the_designs_floor():
    """Change 4 (decided 2026-10-11): without a declared SESOI the floor is used, and the report says
    that absence is confirmed only for effects of at least it."""
    from test_v03_verdict import nested_mice
    from metric_autopsy.effect import design_floor, estimate_effect
    eff, design = estimate_effect(NPR, nested_mice(6, c_young=1.5, c_old=1.5, seed=2), group_col="age",
                                  groups=("young", "old"), replicate_col="mouse", estimand="composition", seed=0)
    floor = design["power"]["floor"]
    assert floor == pytest.approx(design_floor(eff.detail["se"], eff.detail["df"], 0.05))
    assert design["power"]["sesoi_source"] == "floor" and eff.detail["effective_sesoi"] == pytest.approx(floor)
    if eff.status == "NO_DETECTABLE_EFFECT":
        assert "NO_SESOI_FLOOR" in eff.flags and "confirmed only for effects of at least" in eff.reason
    else:  # the corrected interval reached past the floor: never NO DETECTABLE EFFECT then
        assert eff.status in ("INCONCLUSIVE", "DETECTED")


def test_a_sesoi_below_the_floor_is_flagged():
    from test_v03_verdict import nested_mice
    a = run_autopsy(NPR, nested_mice(4, c_young=1.5, c_old=1.5), group_col="age", groups=("young", "old"),
                    replicate_col="mouse", prereg={**RESOLVED, "direction": "two-sided", "sesoi": 1e-4,
                                                   "sesoi_scale": "observed"}, log_path="off")
    assert "SESOI_BELOW_FLOOR" in a.design_adequacy.flags
    assert "cannot show the absence of a change that matters" in a.design_adequacy.reason


# --------------------------------------------------------------------------- #
# the domain
# --------------------------------------------------------------------------- #
def test_the_domain_holds_the_types_where_both_genes_are_detected_in_each_group():
    from metric_autopsy.domain import claim_domain
    d = typed_mice(col_in_b=(0.05, 0.0))  # Col1a1 in about 5% of the young B cells: below 10%
    dom = claim_domain(d, type_col="cell_type", group_col="age", groups=("young", "old"), replicate_col="mouse",
                       genes=("Smad3", "Col1a1"), kind="pair")
    assert dom["domain"] == {"": ["A"]} and dom["outside"] == {"": ["B"]}
    b = next(r for r in dom["strata"][0]["types"] if r["type"] == "B")
    assert b["detection"]["young"]["Col1a1"] < 0.10 and b["detection"]["young"]["Smad3"] > 0.9
    # in each group (the owner's decision of 2026-10-18, by the definition of coupling; the rule was "in at
    # least one group" before): a pair that appears with the change (Col1a1 from 3% to 40% of A's cells) is
    # a claim about expression, outside the domain (tests/test_v04_dossier.py has the line it gets)
    induced = typed_mice(col_in_a=(0.5, 0.03))
    dom2 = claim_domain(induced, type_col="cell_type", group_col="age", groups=("young", "old"),
                        replicate_col="mouse", genes=("Smad3", "Col1a1"), kind="pair")
    a = next(r for r in dom2["strata"][0]["types"] if r["type"] == "A")
    assert a["detection"]["old"]["Col1a1"] < 0.10 <= a["detection"]["young"]["Col1a1"] and not a["in_domain"]


def test_the_claim_is_judged_within_its_domain_and_the_mixture_is_reported_beside():
    """Within A the pair's coupling is the same in both groups; the young B cells carry stray Col1a1
    counts (5% of cells, below the domain's 10%), which move MI over the whole mixture."""
    d = typed_mice(col_in_b=(0.05, 0.0), seed=3)
    a = _typed(d)
    dom = a.effect.detail["domain"]
    assert dom["types"] == ["A"] and a.effect.detail["estimate_over"].startswith("the domain")
    # the domain's raw estimate is A's own (one type: the replicate-level difference of A's cells)
    from metric_autopsy.effect import _infer, _replicate_design, _rep_values
    sa = d[np.asarray(d.obs["cell_type"]) == "A"]
    design = _replicate_design(sa.obs, "age", ("young", "old"), "mouse", [])
    inf = _infer(_rep_values(MI, sa, "age", ("young", "old"), "mouse", design), design, ("young", "old"),
                 "permutation", 0.05, 1000, np.random.default_rng(0))
    assert a.effect.detail["raw_effect"] == pytest.approx(inf["est"], abs=1e-12)
    # beside it: the whole mixture, the types outside the domain, every type on its own
    assert a.effect.detail["mixture"]["raw"]["effect"] is not None
    assert dom["outside"]["types"] == ["B"] and {r["type"] for r in dom["per_type"]} == {"A", "B"}
    md = a.to_markdown()
    assert "## The domain (the verdict is judged within it)" in md and "the whole mixture" in md
    assert "judged within the domain: A" in a.verdict
    # GATE 4 ran on the domain's cells
    g = _g4(a)
    assert g.detail["units"] == dict(n=8, types=["A"], domain=True)


def test_a_claim_whose_genes_are_detected_nowhere_has_an_empty_domain():
    d = typed_mice(col_in_b=(0.0, 0.0), col_in_a=(0.02, 0.02))
    a = _typed(d)
    assert a.cause == "empty_domain" and a.verdict.startswith("INCONCLUSIVE — the claim's domain is empty")


def test_without_cell_types_nothing_changes():
    from test_v03_verdict import nested_mice
    d = nested_mice(4)
    a = run_autopsy(NPR, d, group_col="age", groups=("young", "old"), replicate_col="mouse",
                    prereg={**RESOLVED, "direction": "decrease"}, log_path="off")
    assert "domain" not in a.effect.detail and a.effect.detail["estimate_over"] == "the replicates"
    assert a.effect.detail["p_perm"] == pytest.approx(2 / 70)


# --------------------------------------------------------------------------- #
# change 13: the check of specificity, a line of the report
# --------------------------------------------------------------------------- #
def many_genes(n_mice=4, cells=150, n_genes=160, coupled=(1.5, 0.2), shared=(0.0, 0.0), seed=0):
    """The pair Smad3-Col1a1 coupled per group (young, old); every other gene independent of it, with
    means spread over two orders; `shared` adds a factor common to every gene, per group."""
    rng = np.random.default_rng(seed)
    base = rng.normal(0.3, 0.9, n_genes)
    rows, blocks = [], []
    for gi, age in enumerate(("young", "old")):
        for m in range(n_mice):
            lat = rng.normal(0, 1, cells)
            common = rng.normal(0, 1, cells)[:, None] * shared[gi]
            pair = 1.2 + coupled[gi] * 0.6 * lat[:, None] + rng.normal(0, 0.3, (cells, 2))
            other = base[None, :] + common + rng.normal(0, 0.5, (cells, n_genes))
            blocks.append(rng.poisson(np.exp(np.column_stack([pair, other]))).astype(float))
            rows += [dict(age=age, mouse=f"{age}{m}")] * cells
    names = ["Smad3", "Col1a1"] + [f"G{i}" for i in range(n_genes)]
    return SimpleData(np.vstack(blocks), pd.DataFrame(rows), names)


def _spec(d, **kw):
    return run_autopsy(NPR, d, group_col="age", groups=("young", "old"), replicate_col="mouse",
                       gene_pair=("Smad3", "Col1a1"), pair_metric=metrics.norm_pearson,
                       prereg={**RESOLVED, "direction": "decrease"}, log_path="off", **kw)


def test_a_change_of_one_pair_is_specific_among_matched_controls():
    a = _spec(many_genes())
    sp = a.specificity
    assert sp["applicable"] and sp["n_controls"] > 100 and sp["data"] == "uncorrected"
    assert sp["in_direction"]["k"] <= 0.05 * sp["n_controls"]
    assert len(sp["nearest"]) == 10 and all("genes" in c for c in sp["nearest"])
    assert not sp["shift"]["significant"] and sp["shift"]["warning"] is None
    assert sp["matched"]["Smad3"]["detection"] > 0.9
    md = a.to_markdown()
    assert "## Specificity (reported, not a gate)" in md and "matched controls at least as far" in md


def test_a_shift_shared_by_every_gene_is_flagged_by_the_controls():
    a = _spec(many_genes(coupled=(1.5, 1.5), shared=(0.6, 0.0), seed=1))
    sh = a.specificity["shift"]
    assert sh["significant"] and sh["median_t"] > 0 and "shared by many genes" in sh["warning"]


def test_specificity_is_reported_and_never_decides():
    d = many_genes(seed=2)
    on, off = _spec(d), _spec(d, specificity=False)
    assert on.verdict == off.verdict and on.cause == off.cause and off.specificity == {}
    glob = run_autopsy(lambda data: float(np.mean(as_dense(data.X))), d, group_col="age", groups=("young", "old"),
                       replicate_col="mouse", prereg={**RESOLVED, "direction": "two-sided"}, log_path="off")
    assert glob.specificity["applicable"] is False and "unspecific by construction" in glob.specificity["reason"]


# --------------------------------------------------------------------------- #
# GATE 5: the negative control is reported, not failed (change 13)
# --------------------------------------------------------------------------- #
def test_a_negative_control_outside_its_null_no_longer_fails_the_metric():
    from test_gates import make_clean
    d = make_clean()
    # "negative" control: a pair that is coupled (Actb, Gapdh share a factor); the positive control is real
    g = gates.gate5_controls(metrics.norm_pearson, d, ("Smad3", "Col1a1"), ("Actb", "Gapdh"))
    assert g.status.value != "FAIL"
    assert g.detail["negative_control_outside_null"] and "reported, not a failure" in g.message


# --------------------------------------------------------------------------- #
# GATE 0's known issue: false FAIL over seeds against a target fixed in advance
# --------------------------------------------------------------------------- #
GATE0_FALSE_FAIL_TARGET = 0.05  # at most alpha of the seeds, fixed before measuring (owner, 2026-10-18)
GATE0_SEEDS = 40


_GATE0_RUNS: dict = {}


def _gate0_over_seeds():
    """GATE 0 on the known issue's data over the seeds, computed once for the two tests below."""
    if not _GATE0_RUNS:
        from test_v03_verdict import nested_mice
        from metric_autopsy.effect import _cell_effect
        d = nested_mice(4, c_young=1.5, c_old=1.4, mouse_sd=0.4, seed=3)
        scale = max(abs(_cell_effect(NPR, d, "age", ("young", "old"))), 0.02)
        _GATE0_RUNS.update(data=d, results=[gate0_independence(
            NPR, d, seed=s, effect_scale=scale, sesoi=0.02, corrected=("depth_downsample",), bias_tolerance=0.5)
            for s in range(GATE0_SEEDS)])
    return _GATE0_RUNS["data"], _GATE0_RUNS["results"]


@pytest.mark.xfail(strict=True, reason=(
    "known issue (CHANGELOG 0.4.0.dev0): on these data GATE 0 calls 20% extra dropout a bias of 0.0195-0.0213 "
    "against a tolerance of 0.01 (0.5 x SESOI 0.02) and FAILs in 4 of 40 seeds (8, 11, 36, 37); the groups "
    "do not differ in dropout. GATE 0 alone has no groups: change 12 is applied by run_autopsy (next test)"))
def test_gate0_false_fail_rate_over_seeds_is_within_its_target():
    _, results = _gate0_over_seeds()
    fails = sum(g.status.value == "FAIL" for g in results)
    assert fails / GATE0_SEEDS <= GATE0_FALSE_FAIL_TARGET, f"{fails} of {GATE0_SEEDS} seeds FAIL"


@pytest.mark.xfail(strict=True, reason=(
    "known issue (CHANGELOG 0.4.0.dev0): change 12 cannot show the groups' technical difference within the "
    "tolerance of 0.01 on 4 mice per group, its 90% interval spanning about 0.25 (seeds 8, 11, 36, 37: "
    "-0.081 [-0.210, +0.049], -0.090 [-0.221, +0.041], -0.072 [-0.200, +0.056], -0.073 [-0.187, +0.040]), so "
    "the 4 FAILs stand"))
def test_gate0_false_fail_rate_with_change_12_is_within_its_target():
    from metric_autopsy.report import technical_difference
    d, results = _gate0_over_seeds()
    fails = 0
    for s, g0 in enumerate(results):
        if g0.status.value == "FAIL":
            between = technical_difference(NPR, d, ("Smad3", "Col1a1"), group_col="age", groups=("young", "old"),
                                           replicate_col="mouse", estimand="composition",
                                           margin=g0.detail["bias_tolerance_abs"], seed=s)
            fails += gates.gate0_between_groups(g0, between).status.value == "FAIL"
    assert fails / GATE0_SEEDS <= GATE0_FALSE_FAIL_TARGET, f"{fails} of {GATE0_SEEDS} seeds FAIL with change 12"
