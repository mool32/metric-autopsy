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
