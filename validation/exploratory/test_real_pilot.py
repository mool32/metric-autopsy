"""Tests of the pilot's rules and kernels (real_pilot_select.py, real_pilot.py), on fabricated tables
and random counts only. Not collected by the repository's pytest (validation/exploratory is not a
test path); the pilot's workflow runs them in its dry run, with the engine of the tag:

    python -m pytest -q validation/exploratory/test_real_pilot.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import real_pilot as RP  # noqa: E402
import real_pilot_select as RS  # noqa: E402
from metric_autopsy import SimpleData, metrics  # noqa: E402


def _tms(rows):
    """A TMS table from (donor, sex, stage, cell_type, cells, codetected) rows."""
    out = []
    for donor, sex, stage, ct, n, both in rows:
        for k in range(n):
            out.append(dict(soma_joinid=len(out), dataset_id="d", donor_id=donor, sex=sex, development_stage=stage,
                            cell_type=ct, tissue="t", assay="Smart-seq2", suspension_type="cell",
                            g_Smad3=int(k < both), g_Col1a1=int(k < both), g_Actb=5, g_Gapdh=5, g_Tst=0, g_Lrrc42=0,
                            g_Xist=int(sex == "female"), g_Ddx3y=int(sex == "male"), g_Eif2s3y=0, g_Kdm5d=0, g_Uty=0))
    return pd.DataFrame(out)


Y, O = RS.YOUNG, RS.OLD


def test_case1_rule_counts_qualifying_mice_and_codetection():
    rows = []
    for i in range(3):  # type A: 3 young + 3 old males of >= 20 cells; few co-detected
        rows += [(f"y{i}_M", "male", Y, "A", 30, 2), (f"o{i}_M", "male", O, "A", 30, 2)]
        rows += [(f"y{i}_M", "male", Y, "B", 25, 20), (f"o{i}_M", "male", O, "B", 25, 20)]  # type B: more co-detected
    rows += [("y0_M/y1_M", "male", Y, "C", 400, 400)]  # composite: not a mouse
    rows += [(f"y{i}_F", "female", Y, "C", 30, 30) for i in range(3)] + [("o0_F", "female", O, "C", 30, 30)]
    choice, table = RS.case1_rule(_tms(rows))
    assert list(table["cell_type"]) == ["B", "A"]  # C: only one old female, and the composite id is excluded
    assert choice["cell_type"] == "B" and choice["sexes"] == ["male"] and len(choice["units"]) == 6
    assert choice["codetected"] == 6 * 20


def test_case1_rule_mice_below_20_cells_do_not_count():
    rows = [(f"y{i}_M", "male", Y, "A", 30, 1) for i in range(3)] + [(f"o{i}_M", "male", O, "A", 19, 1) for i in range(3)]
    choice, table = RS.case1_rule(_tms(rows))
    assert choice is None and table.empty


def test_case2_structure_calls_and_exclusions():
    rows = [(f"f{i}", "female", Y, "A", 30, 0) for i in range(4)] + [(f"m{i}", "male", Y, "A", 30, 0) for i in range(4)]
    t = _tms(rows)
    t.loc[t["donor_id"] == "f3", "g_Ddx3y"] = 1   # a female label with Y genes in every cell: contradicts
    st = RS.case2_structure(t).set_index("mouse")
    assert st.loc["f0", "call"] == "female" and st.loc["m0", "call"] == "male"
    assert not st.loc["f3", "kept"] and st.loc["f3", "call"] == "neither"
    t.loc[t["donor_id"] == "f3", "g_Xist"] = 0
    st = RS.case2_structure(t).set_index("mouse")
    assert st.loc["f3", "call"] == "male" and "contradicts" in st.loc["f3", "excluded_because"]


def test_case2_rule_every_mouse_then_fallback():
    rows = []
    for i in range(3):
        rows += [(f"f{i}", "female", Y, "A", 25, 0), (f"m{i}", "male", Y, "A", 25, 0)]
        rows += [(f"f{i}", "female", Y, "B", 40, 0), (f"m{i}", "male", Y, "B", 40 if i else 10, 0)]
    t = _tms(rows)
    st = RS.case2_structure(t)
    choice, table, why = RS.case2_rule(t, st)
    assert choice["cell_type"] == "A" and "every" in choice["rule"]  # B misses one mouse though it has more cells
    t2 = t[~((t["cell_type"] == "A") & (t["donor_id"] == "m0"))]
    rows2 = rows + [("m3", "male", Y, "B", 40, 0), ("m3", "male", Y, "A", 5, 0)]
    t3 = _tms(rows2)
    t3 = t3[~((t3["cell_type"] == "A") & (t3["donor_id"] == "m0"))]
    choice, _, _ = RS.case2_rule(t3, RS.case2_structure(t3))
    assert choice["cell_type"] == "B" and "most mice" in choice["rule"] and len(choice["units"]) == 6
    assert len(t2)


def test_case2_does_not_run_with_two_mice_of_a_sex():
    rows = [(f"f{i}", "female", Y, "A", 30, 0) for i in range(4)] + [(f"m{i}", "male", Y, "A", 30, 0) for i in range(2)]
    t = _tms(rows)
    choice, _, why = RS.case2_rule(t, RS.case2_structure(t))
    assert choice is None and "fewer than 3" in why


def _counts3(rows):
    return pd.DataFrame(rows, columns=["organism", "dataset_id", "donor_id", "assay", "suspension_type", "tissue",
                                       "cell_type", "sex", "development_stage", "cells"])


def _datasets3():
    return pd.DataFrame(dict(
        dataset_id=["d1", "d2", "k1", "t1"], collection_id=["c", "c", "k", "t"],
        collection_name=["Good", "Good", "Multiplexed droplet single-cell RNA-sequencing using natural genetic variation",
                         "Tabula Muris Senis"],
        collection_doi=["10.1/x", "10.1/x", "10.1038/nbt.4042", "10.1/t"], dataset_title=["a", "b", "k", "t"],
        citation=["", "", "", ""], soma_joinid=[0, 1, 2, 3]))


def test_case3_candidates_rules():
    rows = []
    for k in range(5):
        sex = "male" if k else "female"
        rows.append(("hs", "d1", f"D{k}", "10x 3' v3", "cell", "lv", "fib", sex, "50y", 60))
        rows.append(("hs", "d2", f"D{k}", "10x 3' v3", "nucleus", "lv", "fib", sex if k != 4 else "female", "50y", 70))
    rows.append(("hs", "d2", "D9", "10x 3' v3", "nucleus", "lv", "fib", "male", "50y", 49))  # below 50 cells
    for k in range(8):
        rows += [("hs", "k1", f"K{k}", "10x 3' v1", "cell", "blood", "mono", "male", "a", 60),
                 ("hs", "k1", f"K{k}", "10x 3' v2", "cell", "blood", "mono", "male", "a", 60)]
    cand, donors = RS.case3_candidates(_counts3(rows), _datasets3())
    assert list(cand["collection_id"]) == ["k", "c"]  # Kang's set ranks first by donors but is excluded
    assert "DOI 10.1038/nbt.4042" in cand.loc[0, "excluded"] and cand.loc[1, "excluded"] == ""
    assert int(cand.loc[1, "donors"]) == 4  # D4's sex differs between the conditions; D9 is below 50 cells
    assert cand.loc[1, "condition_1"] == "10x 3' v3 | cell" and cand.loc[1, "condition_2"] == "10x 3' v3 | nucleus"


def test_case3_candidates_need_four_donors():
    rows = []
    for k in range(3):
        rows += [("hs", "d1", f"D{k}", "10x 3' v3", "cell", "lv", "fib", "male", "50y", 60),
                 ("hs", "d2", f"D{k}", "10x 3' v3", "nucleus", "lv", "fib", "male", "50y", 60)]
    cand, _ = RS.case3_candidates(_counts3(rows), _datasets3())
    assert cand.empty


def test_sample_ids_caps_and_is_deterministic():
    ids = {"a": np.arange(1000, 1700), "b": np.arange(10, 40)}
    s1, s2 = RS.sample_ids(ids, 500), RS.sample_ids(ids, 500)
    assert np.array_equal(s1, s2) and len(s1) == 530 and np.all(np.isin(np.arange(10, 40), s1))
    assert np.all(np.diff(s1) > 0)


def test_cut_genes_is_d1s_rule():
    rng = np.random.default_rng(0)
    X = rng.poisson(rng.gamma(0.5, 2.0, 2600), size=(50, 2600))
    names = [f"g{j}" for j in range(2600)]
    names[7] = names[8] = "dup"
    dense, kept, missing = RS.cut_genes(X, names, named=("g2599", "absent"))
    assert dense.shape == (50, 2000) and missing == ["absent"] and "g2599" in kept
    assert "dup" in kept or "dup-1" in kept or all(n not in kept for n in ("dup", "dup-1"))
    mean = X.mean(axis=0)
    rest = [j for j in np.argsort(-mean, kind="mergesort") if j != 2599][:1999]
    assert sorted(rest + [2599]) == [int(k) for k in np.where(np.isin(np.arange(2600), rest + [2599]))[0]]
    assert dense.dtype == np.int32


def _random_data(n=300, g=40, seed=1, level=0.0, spread=1.2):
    rng = np.random.default_rng(seed)
    lam = np.outer(np.exp(rng.normal(0, 0.4, n)), np.exp(rng.normal(level, spread, g)))
    z = rng.normal(0, 1, n)
    lam[:, 0] *= np.exp(0.7 * z)
    lam[:, 1] *= np.exp(0.7 * z)
    X = rng.poisson(lam).astype(float)
    genes = [f"G{j}" for j in range(g)]
    return X, genes


@pytest.mark.parametrize("metric", ["mi_3bin", "norm_pearson", "log_cp10k_mean"])
def test_kernels_equal_the_engines_metrics(metric):
    X, genes = _random_data()
    claim = dict(metric=metric, genes=["G0", "G1"] if metric != "log_cp10k_mean" else ["G0", "G3"])
    kern = RP.Kernel(X, genes, claim)
    fn, _, _ = RP.bind(claim, genes)
    data = SimpleData(X, pd.DataFrame(dict(k=np.arange(len(X)) % 3)), genes)
    rng = np.random.default_rng(2)
    masks = rng.random((25, len(X))) < 0.4
    want = np.array([fn(data[m]) for m in masks])
    assert np.allclose(kern.values(masks), want, rtol=0, atol=1e-10)
    assert all(abs(kern.value(m) - w) <= 1e-12 for m, w in zip(masks, want))


def test_pooled_pairs_use_norm_pearson():
    X, genes = _random_data(n=400, g=30, seed=3, level=3.0, spread=0.3)  # well detected: no pair of few cells
    pairs = RS.pooled_pairs(X, genes, list(range(30)), 5)
    data = SimpleData(X, pd.DataFrame(index=range(len(X))), genes)
    used = [g for a, b, _ in pairs for g in (a, b)]
    assert len(set(used)) == len(used) == 10
    for a, b, r in pairs:
        assert abs(metrics.norm_pearson(data, gene_a=a, gene_b=b) - r) < 1e-9
    assert pairs[0][:2] in (("G0", "G1"), ("G1", "G0"))  # the coupled pair is the most correlated


def test_pooled_pairs_of_rare_genes_reach_r_one():
    """norm_pearson's ratio artifact: two genes with one count in each of a few co-detected cells
    are both log1p(1e4 / total), so r is 1 and such pairs rank first (case 3's ranks 501-600)."""
    X = np.zeros((50, 4))
    X[:, 2] = 20 + np.arange(50) % 7
    X[:, 3] = 30 + (np.arange(50) * 3) % 11
    X[:5, 0] = X[:5, 1] = 1
    X[5:8, 0] = 3
    pairs = RS.pooled_pairs(X, ["a", "b", "c", "d"], [0, 1, 2, 3], 1)
    assert pairs[0][:2] == ("a", "b") and abs(pairs[0][2] - 1.0) < 1e-12
    assert [r for _, _, r in pairs] == sorted((r for _, _, r in pairs), reverse=True)


def test_cells_and_donors_tests_on_null_and_effect():
    rng = np.random.default_rng(4)
    n_donor, per = 8, 60
    donor = np.repeat([f"m{k}" for k in range(n_donor)], per)
    group = np.where(np.isin(donor, [f"m{k}" for k in range(4)]), "young", "old")
    strata = np.array([""] * len(donor))
    X = rng.poisson(5.0, size=(len(donor), 20)).astype(float)
    X[group == "old", 0] = rng.poisson(15.0, size=int((group == "old").sum()))
    genes = [f"G{j}" for j in range(20)]
    obs = pd.DataFrame(dict(age=group, mouse=donor))
    data = SimpleData(X, obs, genes)
    claim = dict(metric="log_cp10k_mean", genes=["G0"])
    kern = RP.Kernel(X, genes, claim)
    fn, _, _ = RP.bind(claim, genes)
    a = RP.cells_test(kern, group, strata, ("young", "old"), "increase", np.random.default_rng(0))
    b = RP.donors_test(fn, data, group, donor, strata, ("young", "old"), "increase", 0.5, np.random.default_rng(0))
    assert a["conclusion"] == "effect" and a["p"] == 1 / 1000
    assert b["design"] == "nested" and b["perm_method"] == "exact" and b["perm_n"] == 70
    assert b["conclusion"] == "effect" and b["delta"] > 0
    null = dict(metric="log_cp10k_mean", genes=["G5"])
    fn0, _, _ = RP.bind(null, genes)
    b0 = RP.donors_test(fn0, data, group, donor, strata, ("young", "old"), "increase", 0.5, np.random.default_rng(0))
    assert b0["conclusion"] == "absence" and b0["tost"]  # 90% CI well inside +-0.5


def test_donors_test_paired():
    rng = np.random.default_rng(5)
    donor = np.repeat([f"d{k}" for k in range(6)], 80)
    group = np.tile(np.repeat(["A", "B"], 40), 6)
    X = rng.poisson(4.0, size=(len(donor), 10)).astype(float)
    genes = [f"G{j}" for j in range(10)]
    data = SimpleData(X, pd.DataFrame(dict(group=group, donor=donor)), genes)
    fn, _, _ = RP.bind(dict(metric="log_cp10k_mean", genes=["G1"]), genes)
    b = RP.donors_test(fn, data, group, donor, np.array([""] * len(donor)), ("A", "B"), "two-sided", 0.5,
                       np.random.default_rng(0))
    assert b["design"] == "paired" and b["perm_n"] == 64 and b["conclusion"] in ("absence", "none")


def test_d_conclusions_and_labels():
    assert RP.label_of("SUPPORTED (provisional until replicated) [x]") == "SUPPORTED"
    assert RP.label_of("NOT SUPPORTED — metric invalid: y") == "NOT SUPPORTED"
    assert RP.d_conclusion("SUPPORTED", "provisional") == "effect"
    assert RP.d_conclusion("NOT SUPPORTED", "opposite_direction") == "against"
    assert RP.d_conclusion("NOT SUPPORTED", "metric_invalid_gate0") == "none"
    assert RP.d_conclusion("NO DETECTABLE EFFECT", "no_detectable_effect") == "absence"
    assert RP.call(0.3, 0.01, "decrease") == "against" and RP.call(-0.3, 0.01, "two-sided") == "effect"
    assert RP.call(0.3, 0.2, "increase") == "none"


# --------------------------------------------------------------------------- #
# the second run (the plan's section "Deviation: a second run of cases 1 and 2")
# --------------------------------------------------------------------------- #
def test_sex_call_by_contrast():
    assert RS.sex_call(0.4213, 0.0154, "threshold") == "neither"  # 3_38_F in the first run
    assert RS.sex_call(0.4213, 0.0154, "contrast") == "female"    # and by the rule fixed for the second
    assert RS.sex_call(0.0137, 0.7022, "contrast") == "male"
    assert RS.sex_call(0.30, 0.06, "contrast") == "neither"       # a Y gene in more than 5% of cells
    assert RS.sex_call(0.30, 0.04, "contrast") == "neither"       # Xist in less than ten times that share
    assert RS.sex_call(0.04, 0.0, "contrast") == "neither"        # Xist in less than 5% of cells
    assert RS.sex_call(0.05, 0.005, "contrast") == "female" and RS.sex_call(0.005, 0.05, "contrast") == "male"
    with pytest.raises(ValueError):
        RS.sex_call(0.5, 0.0, "other")


def test_case2_structure_by_contrast_keeps_a_low_xist_female():
    rows = [(f"f{i}", "female", Y, "A", 100, 0) for i in range(3)] + [(f"m{i}", "male", Y, "A", 100, 0) for i in range(3)]
    t = _tms(rows)
    low = (t["donor_id"] == "f0").to_numpy()
    t.loc[low, "g_Xist"] = (np.arange(int(low.sum())) < 42).astype(int)  # Xist in 42% of f0's cells
    t.loc[low, "g_Ddx3y"] = (np.arange(int(low.sum())) < 2).astype(int)  # a Y gene in 2%
    by_threshold = RS.case2_structure(t).set_index("mouse")
    by_contrast = RS.case2_structure(t, "contrast").set_index("mouse")
    assert not by_threshold.loc["f0", "kept"] and by_contrast.loc["f0", "kept"]
    assert by_contrast["kept"].all() and by_contrast.loc["f0", "call"] == "female"


def test_case1_audit_takes_every_cell_and_composite_units():
    rows = [("3_1_M", "male", Y, "A", 30, 3), ("3_1_M", "male", Y, "B", 10, 1), ("3_1_M/3_2_M", "male", Y, "A", 5, 0),
            ("3_2_M", "male", Y, "B", 7, 0), ("24_1_M", "male", O, "A", 40, 4), ("21_1_F", "female", O, "C", 12, 0),
            ("3_3_F", "female", Y, "C", 9, 0), ("18_1_M", "male", "18-month-old stage", "A", 50, 5)]
    t, per, check = RS.case1_audit(_tms(rows))
    assert len(t) == 30 + 10 + 5 + 7 + 40 + 12 + 9  # every cell of the two stages, no type rule, no cap
    assert dict(zip(per["donor_id"], per["cells"])) == {"3_3_F": 9, "21_1_F": 12, "3_1_M": 40, "3_1_M/3_2_M": 5,
                                                         "3_2_M": 7, "24_1_M": 40}
    assert per.set_index("donor_id").loc["3_1_M/3_2_M", "composite"] and check is None
    audit = [dict(donor_id=d, sex=sx, age=a, n_cells=n) for d, sx, a, n in
             (("3_1_M", "male", "3m", 40), ("3_1_M/3_2_M", "male", "3m", 5), ("3_2_M", "male", "3m", 7),
              ("24_1_M", "male", "20m", 40), ("21_1_F", "female", "20m", 12), ("3_3_F", "female", "3m", 9))]
    assert RS.case1_audit(_tms(rows), audit)[2]["same"]
    audit[0]["n_cells"] = 41
    check = RS.case1_audit(_tms(rows), audit + [dict(donor_id="x", sex="male", age="3m", n_cells=1)])[2]
    assert not check["same"] and [d["donor"] for d in check["differences"]] == ["x", "3_1_M"]


def test_design_floor_is_the_tost_at_power_095():
    """At the floor computed from the same SE, the TOST confirms absence exactly where the t test is
    not significant at 5%: |Delta| / SE < t_0.975,df, which under no effect has probability 0.95."""
    te = RP.stats.TEstimate(0.0, 0.004, 13.0, "pooled-stratified")
    assert abs(RP.design_floor(te) - 0.004 * (RP.stats.t_ppf(0.95, 13) + RP.stats.t_ppf(0.975, 13))) < 1e-15
    assert np.isnan(RP.design_floor(RP.stats.TEstimate(0.0, float("nan"), 0.0, "welch")))
    rng = np.random.default_rng(6)
    seen = set()
    for k in range(40):
        donor = np.repeat([f"m{i}" for i in range(9)], 30)
        group = np.where(np.isin(donor, [f"m{i}" for i in range(5)]), "young", "old")
        sex = np.where(np.isin(donor, ["m0", "m1", "m5"]), "female", "male")
        X = rng.poisson(rng.gamma(2.0, 2.0, size=(len(donor), 1)) * np.ones((1, 6))).astype(float)
        X[:, 0] += np.repeat(rng.normal(0, 3, 9).clip(-2, None), 30).round().clip(0, None)
        genes = [f"G{j}" for j in range(6)]
        data = SimpleData(X, pd.DataFrame(dict(age=group, mouse=donor, sex=sex)), genes)
        fn, _, _ = RP.bind(dict(metric="log_cp10k_mean", genes=["G0"]), genes)
        b = RP.donors_test(fn, data, group, donor, sex, ("young", "old"), "decrease", "floor", np.random.default_rng(k))
        assert b["sesoi_used"] == b["floor"] and set(b["by_stratum"]) == {"female", "male"}
        assert b["tost"] == (abs(b["delta"]) / b["se"] < RP.stats.t_ppf(0.975, b["df"]))
        seen.add(b["tost"])
    assert seen == {True, False}


def test_stratum_tests_are_the_audits_exact_permutation():
    import itertools
    young = [0.022069283318606848, 0.016562863987525704, 0.026035638404892772, 0.031368831094429686, 0.04692125067323071]
    old = [0.006265352010074675, 0.01112434189219214, 0.01074778295129472, 0.0072395731320595205]
    y = np.array(young + old)
    g = np.array([1] * 5 + [0] * 4)
    r = RP.stratum_tests(y, g, np.zeros(9, int), {"male": 0}, ("young", "old"), np.random.default_rng(0))["male"]
    obs = abs(np.mean(young) - np.mean(old))
    splits = [abs(y[list(c)].mean() - np.delete(y, list(c)).mean()) for c in itertools.combinations(range(9), 5)]
    assert r["method"] == "exact" and r["n_relabelings"] == 126
    assert abs(r["p"] - np.mean(np.array(splits) >= obs - 1e-12)) < 1e-15 and abs(r["p"] - 2 / 126) < 1e-15
    assert r["delta"] < 0 and r["n"] == {"young": 5, "old": 4}


def test_against_audit_compares_per_mouse_values():
    audit = [dict(donor_id="a", mi_3bin=0.01), dict(donor_id="b", mi_3bin=0.02), dict(donor_id="c", mi_3bin=0.03)]
    v = RP.against_audit({"a": 0.01, "b": 0.0200001, "d": 0.5}, audit, ["Smad3", "Col1a1"])
    assert v["compared"] == 2 and v["only_here"] == ["d"] and v["only_audit"] == ["c"] and abs(v["max_abs_diff"] - 1e-7) < 1e-12
    assert "note" in RP.against_audit({}, audit, ["Actb", "Gapdh"])


def test_merge_expects_only_the_parts_named(tmp_path):
    import json
    for part in ("case1", "case2"):
        d = tmp_path / part
        (d / "reports").mkdir(parents=True)
        (d / f"{part}.json").write_text(json.dumps(dict(synthetic=True, part=part, complete=True, seconds=1.0,
                                                        cases={part: dict(status="run", claims=[])})))
        (d / f"{part}.log").write_text(f"# {part}\n")
    RP.merge([tmp_path / "case1", tmp_path / "case2"], tmp_path / "m2", expect=["case1", "case2"])
    text = (tmp_path / "m2" / "real_pilot.log").read_text()
    assert "MISSING" not in text and "merged 2 of 2 parts" in text
    RP.merge([tmp_path / "case1", tmp_path / "case2"], tmp_path / "m5")
    text = (tmp_path / "m5" / "real_pilot.log").read_text()
    assert text.count("MISSING") == 3 and "merged 2 of 5 parts" in text


def test_a_floor_claim_reaches_the_engine_with_its_floor(tmp_path):
    """1a-floor: (a) and (b) are 1a's (the same permutations), and (d) runs with (b)'s floor as its
    SESOI on the observed scale (delta_min 0.5 x the floor, the engine's default)."""
    import hashlib
    import json
    rng = np.random.default_rng(8)
    mice = [(f"y{i}", "young", "male" if i < 3 else "female") for i in range(5)] + \
           [(f"o{i}", "old", "male" if i < 2 else "female") for i in range(4)]
    rows, blocks = [], []
    for m, age, sex in mice:
        n = 60
        z = rng.normal(0, 1, n)
        lam = np.exp(rng.normal(1.0, 0.5, (n, 12)))
        lam[:, 0] *= np.exp(0.6 * z)
        lam[:, 1] *= np.exp(0.6 * z)
        blocks.append(rng.poisson(lam))
        rows += [dict(mouse=m, age=age, sex=sex)] * n
    X = np.vstack(blocks).astype(np.int32)
    genes = ["Smad3", "Col1a1", "Actb", "Gapdh", "Tst", "Lrrc42"] + [f"G{j}" for j in range(6)]
    RS.P.save_npz(tmp_path / "case1.npz", X, pd.DataFrame(rows), genes)
    sha = hashlib.sha256((tmp_path / "case1.npz").read_bytes()).hexdigest()
    base = dict(metric="mi_3bin", genes=["Smad3", "Col1a1"], direction="decrease",
                signal=dict(kind="coupling", strength=1.0), controls=dict(pos=["Actb", "Gapdh"], neg=["Tst", "Lrrc42"]))
    case = dict(data=dict(file="case1.npz", sha256=sha), group_col="age", groups=["young", "old"], within=["sex"],
                replicate="mouse", claims=[dict(base, id="1a", sesoi=0.1),
                                           dict(base, id="1a-floor", sesoi="floor", sesoi_scale="observed", seed_id="1a")])
    (tmp_path / "out" / "reports").mkdir(parents=True)
    res = RP.run_case("case1", case, tmp_path, tmp_path / "out", lambda *m: None)
    one, fl = res["claims"]
    assert json.dumps(one["a"], sort_keys=True, default=str).replace(str(one["a"]["seconds"]), "") == \
        json.dumps(fl["a"], sort_keys=True, default=str).replace(str(fl["a"]["seconds"]), "")
    assert one["b"]["p"] == fl["b"]["p"] and one["b"]["floor"] == fl["b"]["floor"] == fl["sesoi_used"]
    assert "error" not in fl["d"], fl["d"].get("traceback")
    assert fl["d"]["prereg"]["sesoi"] == fl["b"]["floor"] and fl["d"]["prereg"]["sesoi_scale"] == "observed"
    assert one["d"]["prereg"]["sesoi"] == 0.1 and "sesoi_scale" not in one["d"]["prereg"]
    report = json.loads((tmp_path / "out" / "reports" / "1a-floor.json").read_text())
    assert abs(report["prereg"]["sesoi"] - fl["b"]["floor"]) < 1e-15
