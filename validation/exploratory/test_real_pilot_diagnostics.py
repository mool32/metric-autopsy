"""Tests of the diagnostics' rules (real_pilot_diagnostics.py), on fabricated tables and random counts
only. Not collected by the repository's pytest; the workflow real-pilot-3.yml runs them in its dry run,
with the engine of the tag:

    python -m pytest -q validation/exploratory/test_real_pilot_diagnostics.py
"""
from __future__ import annotations

import gzip
import sys
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import real_pilot as RP  # noqa: E402
import real_pilot_diagnostics as D  # noqa: E402
from anchors import log_cp10k_mean  # noqa: E402
from metric_autopsy import SimpleData, metrics  # noqa: E402


def _counts(n=300, g=12, seed=0):
    rng = np.random.default_rng(seed)
    lam = np.exp(rng.normal(0.0, 1.0, (n, g)))
    return rng.poisson(lam).astype(np.int32), [f"G{i}" for i in range(g)]


def _case(mice=(("y1", "young", "male"), ("y2", "young", "male"), ("o1", "old", "male"), ("o2", "old", "male"),
                ("y3", "young", "female"), ("o3", "old", "female")), cells=40, g=12, seed=0):
    X, genes = _counts(len(mice) * cells, g, seed)
    obs = pd.DataFrame([dict(mouse=m, age=a, sex=s) for m, a, s in mice for _ in range(cells)])
    case = dict(group_col="age", groups=["young", "old"], replicate="mouse", within=["sex"])
    return X, genes, obs, case


def _claim(**kw):
    c = dict(id="t", case="case1", metric="mi_3bin", genes=["G0", "G1"], direction="decrease", males=False,
             sesoi=0.1, other="norm_pearson")
    c.update(kw)
    return c


@pytest.mark.parametrize("kind", ["mi_3bin", "norm_pearson", "log_cp10k_mean"])
def test_reduced_matrix_gives_the_full_matrix_value(kind):
    X, genes = _counts()
    X[:50, 0] = 0
    tot = X.sum(axis=1, dtype=np.int64)
    names = ["G0", "G1"] if kind != "log_cp10k_mean" else ["G0", "G3", "G5"]
    full = SimpleData(X.astype(float), pd.DataFrame(index=np.arange(len(X))), genes)
    if kind == "log_cp10k_mean":
        want = log_cp10k_mean(full, cols=[genes.index(g) for g in names])
    else:
        want = getattr(metrics, kind)(full, gene_a=names[0], gene_b=names[1])
    got = D.measure_fn(kind, names)(D.as_data(D.reduced(X, genes, names, tot), names))
    assert got == want  # exactly: the same columns, the same row sums
    rows = np.arange(100, 200)
    sub = SimpleData(X[rows].astype(float), pd.DataFrame(index=np.arange(len(rows))), genes)
    want_sub = (log_cp10k_mean(sub, cols=[genes.index(g) for g in names]) if kind == "log_cp10k_mean"
                else getattr(metrics, kind)(sub, gene_a=names[0], gene_b=names[1]))
    assert D.measure_fn(kind, names)(D.as_data(D.reduced(X, genes, names, tot, rows), names)) == want_sub


def test_detection_measure_is_the_share_of_cells_with_the_genes():
    M = np.array([[0, 0, 5], [1, 0, 5], [0, 2, 5], [0, 0, 5]], float)
    assert D.measure_fn("detection", ["A", "B"])(D.as_data(M, ["A", "B"])) == 0.5


def test_design_delta_is_the_second_runs_b():
    X, genes, obs, case = _case()
    d = D.Design(obs, case, _claim())
    vals = D.per_mouse("mi_3bin", X, genes, ["G0", "G1"], X.sum(axis=1, dtype=np.int64), d)
    delta, te = d.delta(vals)
    data = SimpleData(X.astype(float), obs.astype(str), genes)
    b = RP.donors_test(partial(metrics.mi_3bin, gene_a="G0", gene_b="G1"), data, obs["age"].to_numpy(),
                       obs["mouse"].to_numpy(), obs["sex"].to_numpy(), ("young", "old"), "decrease", 0.1,
                       np.random.default_rng(0))
    assert delta == pytest.approx(b["delta"], abs=1e-15) and te.se == pytest.approx(b["se"])
    assert vals == pytest.approx(b["values"])


def test_design_males_only_and_nested_only():
    X, genes, obs, case = _case()
    d = D.Design(obs, case, _claim(males=True))
    assert d.within == [] and d.mice == ["o1", "o2", "y1", "y2"] and list(d.labels) == [""]
    obs2 = obs.copy()
    obs2.loc[obs2.index[:5], "age"] = "old"  # y1 now in both groups
    with pytest.raises(ValueError, match="both groups"):
        D.Design(obs2, case, _claim())


def test_nearest_ties_by_file_order():
    det = np.array([0.5, 0.4, 0.6, 0.4, 0.6, 0.1])
    assert D.nearest(det, 0, [1, 2, 3, 4, 5], 4) == [1, 2, 3, 4]
    assert D.nearest(det, 0, [5, 4, 3, 2, 1], 2) == [1, 2]


def _string(tmp_path, genes, links, unmapped=()):
    ids = {g: f"10090.P{i:05d}" for i, g in enumerate(genes) if g not in unmapped}
    with gzip.open(tmp_path / D.STRING_FILES[0], "wt") as fh:
        fh.write("#string_protein_id\tpreferred_name\tprotein_size\tannotation\n")
        for g, p in ids.items():
            fh.write(f"{p}\t{g}\t100\tx\n")
    with gzip.open(tmp_path / D.STRING_FILES[1], "wt") as fh:
        fh.write("protein1 protein2 combined_score\n")
        for a, b, s in links:
            fh.write(f"{ids[a]} {ids[b]} {s}\n{ids[b]} {ids[a]} {s}\n")
    return D.string_names(tmp_path / D.STRING_FILES[0]), tmp_path / D.STRING_FILES[1]


def test_control_pairs_follow_the_rule(tmp_path):
    genes = ["Smad3", "Col1a1", "Actb", "A", "B", "C", "D", "E", "Xist", "Ddx3y"]
    det = np.array([0.30, 0.50, 0.30, 0.31, 0.49, 0.29, 0.51, 0.30, 0.30, 0.50])
    names, links = _string(tmp_path, genes, [("A", "B", 900), ("C", "D", 399)], unmapped=("E",))
    sel = D.control_pairs(det, genes, ["Smad3", "Col1a1"], names, links, n=100, k=3)
    got = {frozenset(p) for p in sel["controls"]}
    # candidates: A, B, C, D (Actb a control of section 4, Xist and Ddx3y sex genes, E not in STRING)
    assert sel["candidates"] == 4 and sel["linked_excluded"] == 1  # A-B linked >= 400; C-D (399) kept
    assert frozenset(("A", "B")) not in got and frozenset(("C", "D")) in got
    assert all(len(p) == 2 for p in got) and len(got) == len(sel["controls"]) == sel["eligible"]
    assert D.control_pairs(det, genes, ["Smad3", "Col1a1"], names, links, n=2, k=3)["controls"] == \
        D.control_pairs(det, genes, ["Smad3", "Col1a1"], names, links, n=2, k=3)["controls"]


def test_string_headers_are_checked(tmp_path):
    with gzip.open(tmp_path / "info.gz", "wt") as fh:
        fh.write("id\tname\n")
    with pytest.raises(ValueError, match="header"):
        D.string_names(tmp_path / "info.gz")
    with gzip.open(tmp_path / "links.gz", "wt") as fh:
        fh.write("a b c\n")
    with pytest.raises(ValueError, match="header"):
        D.string_links(tmp_path / "links.gz", set())


def test_control_genes_and_sets():
    rng = np.random.default_rng(1)
    genes = ["Xist", "Tsix", *D.Y_GENES, "Ddx3x"] + [f"G{i}" for i in range(300)]
    det = rng.random(len(genes))
    sel = D.control_genes(det, genes, "Xist", n=100, k=200)
    flat = [c[0] for c in sel["controls"]]
    assert len(flat) == len(set(flat)) == 100 and not set(flat) & set(D.SEX_GENES)
    sets = D.control_sets(det, genes, list(D.Y_GENES), n=100, k=50)
    assert len(sets["controls"]) == 100 and all(len(set(s)) == 4 for s in sets["controls"])
    assert len({frozenset(s) for s in sets["controls"]}) == 100
    assert not {g for s in sets["controls"] for g in s} & set(D.SEX_GENES)


def test_position_in_the_claimed_direction():
    c = list(np.linspace(-1, 1, 100))
    p = D.position(-0.99, c, "decrease")  # beyond 99 of the 100
    assert p["at_least_as_far"] == 1 and p["p_one_sided"] == pytest.approx(2 / 101) and p["outcome"] == D.PASSES
    p = D.position(-0.9, c, "decrease")
    assert p["at_least_as_far"] == 5 and p["outcome"] == D.FLAGS  # (1 + 5) / 101 >= 0.05
    p = D.position(0.95, c, "decrease")  # against the claim: flags
    assert p["outcome"] == D.FLAGS
    p = D.position(1.5, c, "increase")
    assert p["at_least_as_far"] == 0 and p["p_one_sided"] == pytest.approx(1 / 101) and p["outcome"] == D.PASSES


def test_outcome_of_check_i():
    assert D.outcome_i({})[0] == D.NA
    base = dict(raw_effect=0.02, effect=0.019, p_raw=0.01, raw_detected=True, explained_by_depth=False,
                reversed_by_correction=False, status="DETECTED")
    assert D.outcome_i(base)[0] == D.PASSES
    assert D.outcome_i(dict(base, raw_detected=False, p_raw=0.2))[0] == D.NA
    assert D.outcome_i(dict(base, status="INCONCLUSIVE", effect=0.004))[0] == D.FLAGS
    assert D.outcome_i(dict(base, status="NO_DETECTABLE_EFFECT"))[0] == D.FLAGS
    assert D.outcome_i(dict(base, explained_by_depth=True, status="INCONCLUSIVE"))[1].startswith("explained")
    assert D.outcome_i(dict(base, reversed_by_correction=True, status="INCONCLUSIVE"))[1].startswith("the correction")


def _per(values, removal, doses):
    return {m: dict(value=v, removal=removal[m], dose=doses[m]) for m, v in values.items()}


def test_check_iii_parts_a_and_b():
    X, genes, obs, case = _case()
    d = D.Design(obs, case, _claim(males=True))  # young y1, y2 (0.03); old o1, o2 (0.01): |Delta| 0.02
    vals = {"y1": 0.03, "y2": 0.03, "o1": 0.01, "o2": 0.01}
    ok = _per(vals, {"y1": -0.029, "y2": -0.029, "o1": -0.009, "o2": -0.009},
              {m: {"1.0": 0.001, "50.0": 0.03} for m in vals})
    r = D.check_iii(_claim(males=True), d, ok)
    s = r["strata"]["all"]
    assert r["outcome"] == D.PASSES and s["higher"] == "young" and s["max_dose"] == "50.0"
    small = _per(vals, {m: ok[m]["removal"] for m in vals}, {m: {"1.0": 0.001, "50.0": 0.015} for m in vals})
    assert D.check_iii(_claim(males=True), d, small)["failed"] == ["all: (b)"]
    weak = _per(vals, {m: -0.01 for m in vals}, {m: {"1.0": 0.03} for m in vals})
    assert D.check_iii(_claim(males=True), d, weak)["failed"] == ["all: (a)"]


def test_check_iv_signs_and_detection():
    X, genes, obs, case = _case(mice=tuple((f"y{i}", "young", "male") for i in range(4))
                                + tuple((f"o{i}", "old", "male") for i in range(4)))
    d = D.Design(obs, case, _claim(males=True))
    up = {m: (1.0 if m.startswith("o") else 0.0) + 0.01 * i for i, m in enumerate(d.mice)}
    down = {m: -v for m, v in up.items()}
    flat = {m: 0.01 * i for i, m in enumerate(d.mice)}
    assert D.check_iv(_claim(), d, up, down)["outcome"] == D.FLAGS  # opposite and detected (p = 1/70 * 2)
    assert D.check_iv(_claim(), d, up, up)["outcome"] == D.PASSES
    noisy = {m: v * (-1 if i % 2 else 1) for i, (m, v) in enumerate(flat.items())}
    r = D.check_iv(_claim(), d, noisy, {m: -v for m, v in noisy.items()})
    assert all(not s["own"]["detected"] and not s["other"]["detected"] for s in r["strata"].values())
    assert r["outcome"] == D.PASSES  # opposite signs, neither detected


def test_responses_extremes():
    rng = np.random.default_rng(3)
    n = 400
    x = rng.poisson(3.0, n)
    M = np.column_stack([x, x, rng.poisson(50.0, n)]).astype(float)  # a perfectly coupled pair
    r = D.responses("mi_3bin", ["A", "B"], M, (1, 2, 3))
    assert r["value"] > 0.3 and r["removal"] < -0.25 and set(r["dose"]) == {str(s) for s in D.STRENGTHS}
    L = np.column_stack([np.zeros(n), rng.poisson(50.0, n)]).astype(float)  # a gene with no counts
    r = D.responses("log_cp10k_mean", ["Xist"], L, (1, 2, 3))
    assert r["value"] == 0 and r["removal"] == 0 and all(v == 0 for v in r["dose"].values())
    L[:, 0] = rng.poisson(2.0, n)
    r = D.responses("log_cp10k_mean", ["Xist"], L, (1, 2, 3))
    assert r["removal"] == pytest.approx(-r["value"]) and r["dose"]["100.0"] > r["dose"]["2.0"] > 0


def test_separates_and_the_prediction():
    F, P, N = D.FLAGS, D.PASSES, D.NA
    assert D.separates({"1a": F, "1a-males": F, "2a": P, "2b": P}) == "separates"
    assert D.separates({"1a": N, "1a-males": F, "2a": P, "2b": P}) == "partly"
    assert D.separates({"1a": F, "1a-males": F, "2a": F, "2b": P}) == "does not separate"
    assert D.separates({"1a": P, "1a-males": N, "2a": P, "2b": P}) == "does not separate"
    assert D.separates({"1a": "error", "1a-males": F, "2a": P, "2b": P}).startswith("undetermined")
    table = {c: {"1a": F, "1a-males": F, "1b": F, "2a": P, "2b": P} for c in D.CHECKS}
    assert D.prediction(table)["holds"]
    table["iv"]["2b"] = F
    assert not D.prediction(table)["holds"]
