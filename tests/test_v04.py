"""v0.4: the changes of validation/plans/v0.4_spec.md, each with its regression test.

- change 1: GATE 4 measures the response with the claim's own statistic, injected into the group
  where the probed genes have more counts per cell, its interval over the replicates; every random
  draw keyed by the data and the claim's canonical form;
- change 2: the module probe lowers the module by thinning its own genes, without a sham;
- changes 6, 7, 8, 14, 15 and 16 (the decisions of 2026-10-16 and 2026-10-17); those of 2026-10-18 are
  in test_v04_domain.py.
"""
from __future__ import annotations

from functools import partial

import numpy as np
import pandas as pd
import pytest

from metric_autopsy import SimpleData, gate4_signal_response, injected_signal, metrics, provenance, run_autopsy
from metric_autopsy.core import as_dense

RESOLVED = {"estimand": "composition", "judgment_pending": False}


def sex_marker(n_mice=4, cells=150, genes=30, xist_female=8.0, seed=0, constant_marker=None):
    """Females carry a marker (Xist) in every cell, males none; the other genes are alike."""
    rng = np.random.default_rng(seed)
    means = np.exp(rng.normal(0.5, 1.0, genes))
    rows, blocks = [], []
    for sex in ("female", "male"):
        for m in range(n_mice):
            depth = rng.lognormal(0, 0.2, (cells, 1))
            X = rng.poisson(means[None, :] * depth)
            marker = (np.full(cells, constant_marker) if constant_marker is not None else
                      rng.poisson(xist_female * depth[:, 0]) if sex == "female" else np.zeros(cells))
            blocks.append(np.column_stack([X, marker]))
            rows += [dict(sex=sex, mouse=f"{sex[0]}{m}", age="young")] * cells
    names = [f"G{i}" for i in range(genes)] + ["Xist"]
    return SimpleData(np.vstack(blocks).astype(float), pd.DataFrame(rows), names)


def mean_lognorm_xist(data):
    X = as_dense(data.X)
    tot = X.sum(axis=1)
    tot[tot == 0] = 1
    return float(np.mean(np.log1p(X[:, list(data.var_names).index("Xist")] / tot * 1e4)))


def _claim(data, groups, direction, **kw):
    kw.setdefault("signal_test", injected_signal.module(["Xist"]))
    kw.setdefault("replicate_col", "mouse")
    return run_autopsy(mean_lognorm_xist, data, group_col="sex", groups=groups,
                       prereg={**RESOLVED, "direction": direction, "sesoi": 0.5}, log_path="off", **kw)


def _g4(a):
    return next(r for r in a.results if r.gate == 4)


# --------------------------------------------------------------------------- #
# change 1: GATE 4 on the claim's statistic, in the group with more probed counts
# --------------------------------------------------------------------------- #
def test_gate4_injects_where_the_marker_is_whatever_the_wording():
    """R1's lesson (D1): "Xist is higher in females" and "Xist is lower in males" are one claim.
    The module goes into the females, who carry Xist, in both wordings; the response is the change
    of the replicate-level difference, its interval over the 4 female mice."""
    d = sex_marker()
    a = _claim(d, ("female", "male"), "decrease")
    b = _claim(d, ("male", "female"), "increase")
    ga, gb = _g4(a), _g4(b)
    for g in (ga, gb):
        assert g.detail["mode"] == "claim statistic" and g.detail["injected_group"] == "female"
        assert g.detail["n_rep"] == 4 and sorted(g.detail["replicate_responses"]) == [f"f{i}" for i in range(4)]
        assert g.detail["direction"] == "decrease" and g.status.value == "PASS"
        assert g.detail["probed_counts_per_cell"]["male"] == 0
    # the same draws: the same response in the females, the difference moving in opposite senses
    assert ga.detail["mean_response"] == pytest.approx(gb.detail["mean_response"], abs=1e-12)
    assert ga.detail["mean_response"] < -0.3
    assert ga.detail["statistic_change"] == pytest.approx(-gb.detail["statistic_change"], abs=1e-12)
    assert ga.detail["statistic_change"] < 0  # female - male falls when the females' Xist falls
    # the probe's own direction (a decrease), not the claim's, is what the response is judged by
    assert "expected decrease" in ga.message and "injected into female" in ga.message
    # the module removes only its share of the counts, and says so
    assert 0 < ga.detail["depth_loss"] < 0.2 and "removes" in ga.message
    # v0.3's all-cell response stays as a diagnostic
    assert ga.detail["all_cells"]["n_rep"] == 200 and ga.detail["all_cells"]["mean_response"] < 0


def test_a_mirrored_wording_gets_the_same_draws_and_the_same_verdict():
    d = sex_marker(seed=1)
    a = _claim(d, ("female", "male"), "decrease")
    b = _claim(d, ("male", "female"), "increase")
    assert a.params["claim"] == b.params["claim"]
    assert a.verdict == b.verdict and a.cause == b.cause
    assert {k: v.status for k, v in a.fields().items()} == {k: v.status for k, v in b.fields().items()}
    assert a.effect.detail["effect"] == pytest.approx(-b.effect.detail["effect"], abs=1e-12)
    assert a.effect.detail["p_perm"] == b.effect.detail["p_perm"]
    for ra, rb in zip(sorted(a.results, key=lambda r: r.gate), sorted(b.results, key=lambda r: r.gate)):
        assert ra.gate == rb.gate and ra.status == rb.status


def test_two_claims_on_the_same_data_get_different_draws():
    c1 = provenance.canonical_claim("m", ["Xist"], ("female", "male"), "decrease")
    c2 = provenance.canonical_claim("m", ["Xist"], ("male", "female"), "increase")
    c3 = provenance.canonical_claim("m", ["Tsix"], ("female", "male"), "decrease")
    c4 = provenance.canonical_claim("m", ["Xist"], ("female", "male"), "two-sided")
    assert c1 == c2 and c1 != c3 and c1 != c4
    seeds = {provenance.claim_seed("d" * 64, c, 0, "gate4") for c in (c1, c3, c4)}
    assert len(seeds) == 3
    assert provenance.claim_seed("d" * 64, c1, 0, "gate4") != provenance.claim_seed("d" * 64, c1, 0, "effect")
    assert provenance.claim_seed("d" * 64, c1, 0, "gate4") != provenance.claim_seed("d" * 64, c1, 1, "gate4")


def test_on_an_exact_tie_both_groups_are_injected_and_the_weaker_is_judged():
    """A marker with the same counts in every cell of both groups: the probe runs in each group. The rule is
    GATE 4's own, so it is tested on the gate: since the owner's second decisions of 2026-10-19 the autopsy
    gives GATE 4 one draw of the correction (equal depth), whose thinning breaks a tie of the raw counts."""
    d = sex_marker(constant_marker=3)
    probe = injected_signal.module(["Xist"])
    g = gate4_signal_response(mean_lognorm_xist, d, probe, direction=getattr(probe, "direction", None), alpha=0.05,
                              seed=0, group_col="sex", groups=("female", "male"), replicate_col="mouse",
                              construct_change=0.5)
    assert g.detail["tie"] is True and len(g.detail["runs"]) == 2
    weaker = min(g.detail["runs"], key=lambda r: r["interval"]["signed_lower"])
    assert g.detail["injected_group"] == weaker["group"]
    assert _g4(_claim(d, ("female", "male"), "two-sided")).detail["data"] == "equal depth"


def test_without_a_replicate_unit_gate4_judges_the_all_cell_response():
    d = sex_marker()
    g = _g4(_claim(d, ("female", "male"), "decrease", replicate_col=None))
    assert g.detail["mode"] == "all cells" and g.detail["n_rep"] == 200
    assert "over all cells" in g.message


def test_gate4_on_a_coupling_claim_uses_the_replicates_of_the_group_with_more_counts():
    from test_gates import _assemble, _block
    rng = np.random.default_rng(4)
    blocks = []
    for age, eff in (("young", 1.0), ("old", 0.6)):  # the young carry more counts of the pair
        for m in range(4):
            X, o = _block(rng, 150, coupling=1.5, efficiency=eff, age=age, sex="female")
            o["mouse"] = f"{age}{m}"
            blocks.append((X, o))
    d = _assemble(blocks)
    npr = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
    a = run_autopsy(npr, d, group_col="age", groups=("young", "old"), replicate_col="mouse",
                    gene_pair=("Smad3", "Col1a1"), signal_test=injected_signal.coupling("Smad3", "Col1a1"),
                    prereg={**RESOLVED, "direction": "decrease"}, log_path="off")
    g = _g4(a)
    assert g.detail["injected_group"] == "young" and g.detail["n_rep"] == 4
    assert g.detail["reference"] == "sham" and g.detail["direction"] == "increase"
    assert g.status.value == "PASS" and g.detail["statistic_change"] > 0


def test_gate4_called_alone_keeps_the_all_cell_rule():
    """Without the claim's groups (a direct call) GATE 4 is v0.3's all-cell gate."""
    from test_gates import make_clean
    g = gate4_signal_response(partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1"), make_clean(),
                              injected_signal.coupling("Smad3", "Col1a1"), delta_min=0.05)
    assert g.detail["mode"] == "all cells" and g.status.value == "PASS"


def test_the_data_hash_does_not_depend_on_the_storage():
    sparse = pytest.importorskip("scipy.sparse")
    d = sex_marker(n_mice=2, cells=40)

    class Sparse:
        pass
    s = Sparse()
    s.X, s.obs, s.var_names = sparse.csr_matrix(d.X.astype(np.float32)), d.obs, d.var_names
    assert provenance.sha256_data(d, ["sex"]) == provenance.sha256_data(s, ["sex"])
    e = SimpleData(d.X.copy(), d.obs, d.var_names)
    e.X[3, 4] += 1
    assert provenance.sha256_data(d, ["sex"]) != provenance.sha256_data(e, ["sex"])


# --------------------------------------------------------------------------- #
# change 6: a correction can only remove a conclusion, never create one
# --------------------------------------------------------------------------- #
def created_by_the_correction(n_mice=6, cells=300, seed=0):
    """Group a is twice as deep, but X has the same absolute counts in both groups (its share of
    a's counts is half): X's detection is equal on the raw data and falls in a at equal depth."""
    rng = np.random.default_rng(seed)
    means = np.exp(rng.normal(0.5, 1.0, 30))
    rows, blocks = [], []
    for grp, depth in (("a", 2.0), ("b", 1.0)):
        for m in range(n_mice):
            X = rng.poisson(means[None, :] * depth, (cells, 30))
            blocks.append(np.column_stack([X, rng.poisson(1.0, cells)]))
            rows += [dict(grp=grp, mouse=f"{grp}{m}")] * cells
    return SimpleData(np.vstack(blocks).astype(float), pd.DataFrame(rows), [f"G{i}" for i in range(30)] + ["X"])


def detection_of_x(data):
    return float(np.mean(np.asarray(data.X)[:, -1] > 0))


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_an_effect_only_the_correction_shows_is_not_detected(seed):
    from metric_autopsy.effect import estimate_effect
    eff, design = estimate_effect(detection_of_x, created_by_the_correction(seed=seed), group_col="grp",
                                  groups=("a", "b"), replicate_col="mouse", estimand="composition",
                                  sesoi=0.1, seed=seed)
    d = eff.detail
    assert d["detected_only_after_correction"] is True and not d["raw_detected"]
    assert d["effect"] < -0.15 and d["p_perm"] < 0.05  # the corrected test alone would detect it
    assert eff.status == "INCONCLUSIVE" and "can only remove a conclusion" in eff.reason
    # the capture model's check is in the report: equal depth after the thinning
    cc = design["equalization"]["capture_check"]
    assert abs(cc["equalized"]["relative_difference"]["depth"]) < 0.02
    assert cc["raw"]["relative_difference"]["depth"] > 0.9


def test_the_verdict_names_an_effect_created_by_the_correction():
    a = run_autopsy(detection_of_x, created_by_the_correction(), group_col="grp", groups=("a", "b"),
                    replicate_col="mouse", signal_test=injected_signal.module(["X"]),
                    prereg={**RESOLVED, "direction": "increase", "sesoi": 0.1}, log_path="off")
    assert a.cause == "detected_only_after_correction" and a.verdict.startswith("INCONCLUSIVE")
    g2 = next(r for r in a.results if r.gate == 2)
    assert g2.status.value == "WARN" and "only the corrected effect is detected" in g2.message


# --------------------------------------------------------------------------- #
# change 14: a detected effect smaller than the SESOI is an outcome of its own
# --------------------------------------------------------------------------- #
def test_a_detected_effect_within_the_sesoi_is_effect_below_sesoi():
    d = sex_marker(seed=2)
    small = _claim(d, ("female", "male"), "decrease")  # SESOI 0.5 is far below the effect
    assert small.cause in ("provisional", "replicated") and small.verdict.startswith("SUPPORTED")
    big = run_autopsy(mean_lognorm_xist, d, group_col="sex", groups=("female", "male"), replicate_col="mouse",
                      signal_test=injected_signal.module(["Xist"]),
                      prereg={**RESOLVED, "direction": "decrease", "sesoi": 50.0, "sesoi_scale": "observed"},
                      log_path="off")  # on the metric's scale: GATE 4 is not what this test is about
    assert big.effect.status == "DETECTED" and big.effect.detail["tost_equivalent"] is True
    assert big.cause == "effect_below_sesoi" and big.verdict.startswith("EFFECT BELOW SESOI")
    assert "smaller than the smallest effect that matters" in big.verdict


def test_effect_below_sesoi_in_the_order_of_the_rule():
    from metric_autopsy.core import Assessment
    from metric_autopsy.report import Autopsy, decide_cause

    def autopsy(ef_detail, rp="NOT_RUN", pending=False, direction="two-sided"):
        return Autopsy("m", [], {"judgment_pending": pending, "direction": direction},
                       Assessment("PASS", "mv"), Assessment("ADEQUATE", "da"),
                       Assessment("DETECTED", "effect +1 [+0.5, +1.5], p=0.01", [], ef_detail), Assessment(rp, "rp"))
    below = dict(effect=1.0, tost_equivalent=True, ci_2alpha=(0.6, 1.4), effective_sesoi=2.0)
    assert decide_cause(autopsy(below))[1] == "effect_below_sesoi"
    assert decide_cause(autopsy(below, pending=True))[1] == "effect_below_sesoi"
    assert decide_cause(autopsy(below, rp="NOT_REPLICATED"))[1] == "not_replicated"
    assert decide_cause(autopsy(below, direction="increase"))[1] == "opposite_direction"  # +1: a decrease
    assert decide_cause(autopsy(dict(below, tost_equivalent=False)))[1] == "provisional"
    created = Autopsy("m", [], {"judgment_pending": False, "direction": "two-sided"}, Assessment("PASS", "mv"),
                      Assessment("ADEQUATE", "da"),
                      Assessment("INCONCLUSIVE", "x", [], dict(detected_only_after_correction=True)),
                      Assessment("NOT_RUN", "rp"))
    assert decide_cause(created)[1] == "detected_only_after_correction"


def test_an_absence_only_the_correction_shows_is_flagged():
    """The TOST runs on the uncorrected data too (decided 2026-10-12): where only the corrected
    data are equivalent, NO DETECTABLE EFFECT stands but carries ABSENCE_ONLY_AFTER_CORRECTION."""
    from metric_autopsy.effect import estimate_effect
    rng = np.random.default_rng(5)
    rows, blocks = [], []
    # the counts of X scale with depth; a's mice span depths 1 to 3, so the raw difference is wide
    for grp, depths in (("a", (1.0, 1.0, 3.0, 3.0, 1.0, 3.0)), ("b", (1.0,) * 6)):
        for m, dep in enumerate(depths):
            X = rng.poisson(2.0 * dep, (400, 20))
            blocks.append(np.column_stack([X, rng.poisson(0.5 * dep, 400)]))
            rows += [dict(grp=grp, mouse=f"{grp}{m}")] * 400
    d = SimpleData(np.vstack(blocks).astype(float), pd.DataFrame(rows), [f"G{i}" for i in range(20)] + ["X"])

    def mean_x(data):
        return float(np.asarray(data.X)[:, -1].mean())
    eff, _ = estimate_effect(mean_x, d, group_col="grp", groups=("a", "b"), replicate_col="mouse",
                             estimand="composition", sesoi=0.3, sesoi_scale="observed", seed=0)
    assert eff.detail["tost_equivalent_raw"] is False
    if eff.status == "NO_DETECTABLE_EFFECT":
        assert "ABSENCE_ONLY_AFTER_CORRECTION" in eff.flags and "only after the depth correction" in eff.reason
    else:  # the corrected interval must then reach the SESOI: never flagged without an absence
        assert "ABSENCE_ONLY_AFTER_CORRECTION" not in eff.flags


# --------------------------------------------------------------------------- #
# change 7: composition, reported with the claim; the claim per type, not combined
# --------------------------------------------------------------------------- #
def two_types(n_mice=4, cells=200, share_b=(0.2, 0.6), seed=0, only_young_c=False):
    """Type B expresses X twice as much as type A; the old mice hold more B. Within each type the
    groups are alike, so the mixture's difference of X is composition alone."""
    rng = np.random.default_rng(seed)
    means = np.exp(rng.normal(0.5, 1.0, 20))
    rows, blocks = [], []
    for age, sb in zip(("young", "old"), share_b):
        for m in range(n_mice):
            n_b = int(round(cells * sb))
            types = ["A"] * (cells - n_b) + ["B"] * n_b
            if only_young_c and age == "young":
                types = types[:-30] + ["C"] * 30
            for t in types:
                rows.append(dict(age=age, mouse=f"{age}{m}", cell_type=t, tissue="lung" if t != "C" else "spleen"))
            x_rate = np.array([4.0 if t == "B" else 2.0 for t in types])
            X = rng.poisson(means[None, :], (cells, 20))
            blocks.append(np.column_stack([X, rng.poisson(x_rate)]))
    return SimpleData(np.vstack(blocks).astype(float), pd.DataFrame(rows), [f"G{i}" for i in range(20)] + ["X"])


def mean_lognorm_x(data):
    X = as_dense(data.X)
    tot = X.sum(axis=1)
    tot[tot == 0] = 1
    return float(np.mean(np.log1p(X[:, -1] / tot * 1e4)))


def _mixture(d, **kw):
    return run_autopsy(mean_lognorm_x, d, group_col="age", groups=("young", "old"), replicate_col="mouse",
                       signal_test=injected_signal.module(["X"]), composition=["cell_type", "tissue"],
                       prereg={**RESOLVED, "direction": "increase", "sesoi": 0.05}, log_path="off", **kw)


def test_the_composition_is_compared_at_the_replicate_level():
    from metric_autopsy.composition import composition_difference
    d = two_types()
    c = composition_difference(d.obs, "cell_type", "age", ("young", "old"), "mouse")
    assert c["method"] == "exact" and c["n_relabelings"] == 70 and c["design"] == "nested"
    b = next(r for r in c["levels"] if r["level"] == "B")
    assert b["first"] == pytest.approx(0.2) and b["second"] == pytest.approx(0.6) and b["difference"] == pytest.approx(0.4)
    assert c["tvd"] == pytest.approx(0.4) and c["tvd_p"] == pytest.approx(2 / 70) and c["differs"]
    same = composition_difference(two_types(share_b=(0.4, 0.4)).obs, "cell_type", "age", ("young", "old"), "mouse")
    assert same["tvd"] == pytest.approx(0.0) and not same["differs"]


def test_a_claim_about_a_mixture_carries_its_difference_of_composition():
    """Since 2026-10-18 the claim is judged within its domain, the types weighted equally: within each
    type the groups are alike, so the domain finds no rise, while the whole mixture, reported beside
    with the difference of composition, rises."""
    a = _mixture(two_types())
    comp = a.composition["columns"]["cell_type"]
    assert comp["differs"] and "composition in cell_type differs" in a.verdict
    dom = a.effect.detail["domain"]
    assert dom["types"] == ["A", "B"] and not a.verdict.startswith("SUPPORTED")
    assert {r["type"]: r["units"] for r in dom["per_type"]} == {"A": 8, "B": 8}
    mix = a.effect.detail["mixture"]["raw"]
    assert mix["effect"] < 0 and mix["detected"]  # young - old: the mixture rises with the share of B
    assert abs(a.effect.detail["raw_effect"]) < abs(mix["effect"]) / 3
    md = a.to_markdown()
    assert "## Composition (reported, not a gate)" in md and "## The domain (the verdict is judged within it)" in md
    assert a.to_dict()["composition"]["columns"]["cell_type"]["tvd"] == pytest.approx(0.4)


def test_a_type_absent_from_one_group_is_unidentifiable_and_small_types_are_not_judged():
    a = _mixture(two_types(only_young_c=True))
    rows = {r["type"]: r for r in a.composition["domain"]["strata"][0]["types"]}
    assert rows["C"]["absent_from"] == "old" and not rows["C"]["qualifies"]
    assert a.composition["columns"]["cell_type"]["absent"] == {"C": "old"}
    assert "C" not in a.effect.detail["domain"]["types"]
    few = _mixture(two_types(cells=60, share_b=(0.2, 0.25)))  # 12 and 15 cells of B per mouse
    rows = {r["type"]: r for r in few.composition["domain"]["strata"][0]["types"]}
    assert not rows["B"]["qualifies"] and few.effect.detail["domain"]["types"] == ["A"]


# --------------------------------------------------------------------------- #
# change 15: the engine's norm_pearson beside a metric of a pair's coupling
# --------------------------------------------------------------------------- #
def test_norm_pearson_is_reported_beside_a_coupling_metric_by_stratum():
    from test_gates import _assemble, _block
    rng = np.random.default_rng(6)
    blocks = []
    for sex in ("female", "male"):
        for age, c in (("young", 2.0), ("old", 0.4)):
            for m in range(4):
                X, o = _block(rng, 120, coupling=c, efficiency=1.0, age=age, sex=sex)
                o["mouse"] = f"{sex}{age}{m}"
                blocks.append((X, o))
    d = _assemble(blocks)
    mi = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")
    a = run_autopsy(mi, d, group_col="age", groups=("young", "old"), within=["sex"], replicate_col="mouse",
                    gene_pair=("Smad3", "Col1a1"), prereg={**RESOLVED, "direction": "decrease"}, log_path="off")
    ref = a.reference
    assert ref["measure"].startswith("norm_pearson") and [r["stratum"] for r in ref["rows"]] == [
        "female", "male", "over the strata"]
    for r in ref["rows"]:
        # young − old, the effect test's convention since 2026-10-19: both fall with age, so both are positive
        assert r["metric"]["difference"] > 0 and r["reference"]["difference"] > 0
    assert "beside the metric (reported, not a gate)" in a.to_markdown()


# --------------------------------------------------------------------------- #
# change 16: units whose ids name several animals
# --------------------------------------------------------------------------- #
def test_units_that_name_several_animals_are_flagged():
    d = sex_marker()
    mouse = d.obs["mouse"].to_numpy().copy()
    mouse[(mouse == "m0")] = "m0/m1"  # a composite unit beside the unit m1
    d2 = SimpleData(d.X, d.obs.assign(mouse=mouse), d.var_names)
    a = _claim(d2, ("female", "male"), "decrease")
    assert "COMPOSITE_UNITS" in a.design_adequacy.flags
    cu = a.design_adequacy.detail["composite_units"]
    assert cu == [dict(unit="m0/m1", names=["m0", "m1"], overlaps=["m1"])]
    assert "name several animals" in " ".join(a.design_adequacy.detail["notes"])


def test_the_composition_reaches_the_same_report_through_the_api_the_cli_and_mcp(tmp_path, capsys):
    import json
    anndata = pytest.importorskip("anndata")
    from metric_autopsy import cli, mcp_server
    d = two_types()
    path = tmp_path / "types.h5ad"
    ad = anndata.AnnData(X=d.X, obs=d.obs.copy(), var=pd.DataFrame(index=list(d.var_names)))
    ad.obs_names = [f"cell_{i}" for i in range(d.n_obs)]
    ad.write_h5ad(path)
    args = ["--metric", "norm_pearson", "--gene-a", "G0", "--gene-b", "X", "--group-col", "age", "--groups",
            "young", "old", "--within", "--replicate-col", "mouse", "--estimand", "composition",
            "--direction", "increase", "--composition", "cell_type", "tissue", "--no-log"]
    cli.main(["--h5ad", str(path), *args, "--json", str(tmp_path / "cli.json")])
    capsys.readouterr()
    mcp_server.autopsy_report(str(path), metric="norm_pearson", gene_a="G0", gene_b="X", group_col="age",
                              groups=["young", "old"], replicate_col="mouse", estimand="composition",
                              direction="increase", composition=["cell_type", "tissue"],
                              json_path=str(tmp_path / "mcp.json"), log_path="off")
    api = run_autopsy(partial(metrics.norm_pearson, gene_a="G0", gene_b="X"), anndata.read_h5ad(path),
                      group_col="age", groups=("young", "old"), gene_pair=("G0", "X"), replicate_col="mouse",
                      prereg={"estimand": "composition", "direction": "increase"}, composition=["cell_type", "tissue"],
                      log_path="off")
    reports = [json.loads((tmp_path / f).read_text()) for f in ("cli.json", "mcp.json")]
    for rep in reports:
        assert rep["verdict"] == api.verdict
        assert rep["composition"]["columns"]["cell_type"]["tvd"] == pytest.approx(api.composition["columns"]["cell_type"]["tvd"])
        assert ({r["type"]: r["raw"]["effect"] for r in rep["fields"]["effect"]["detail"]["domain"]["per_type"]}
                == pytest.approx({r["type"]: r["raw"]["effect"] for r in api.effect.detail["domain"]["per_type"]}))
        assert rep["reference"]["rows"][0]["reference"]["difference"] == pytest.approx(
            api.reference["rows"][0]["reference"]["difference"])


# --------------------------------------------------------------------------- #
# change 8: no number without a script in the docstrings
# --------------------------------------------------------------------------- #
def test_the_qc_docstring_cites_no_number_without_a_script():
    from metric_autopsy import qc
    assert "2.4x" not in (qc.__doc__ or "") and "REPORT.md" in qc.__doc__
