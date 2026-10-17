"""Tests of the second diagnostics' rules (real_pilot_4_controls.py, real_pilot_4.py), on fabricated
tables and random counts only. The workflow real-pilot-4.yml runs them in its dry run, with the
engine of the tag:

    python -m pytest -q validation/exploratory/test_real_pilot_4.py
"""
from __future__ import annotations

import gzip
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import real_pilot_4 as R  # noqa: E402
import real_pilot_4_controls as C  # noqa: E402
import real_pilot_diagnostics as D  # noqa: E402
from metric_autopsy import SimpleData, metrics, stats  # noqa: E402


def test_the_lists_are_the_first_diagnostics():
    assert C.SEX_GENES == D.SEX_GENES and C.PAIR == D.PAIR and C.Y_GENES == D.Y_GENES and C.S4_CONTROLS == D.S4_CONTROLS


def test_gene_stats_sparse_and_dense_agree():
    rng = np.random.default_rng(0)
    X = rng.poisson(0.4, (200, 30)).astype(np.int32)
    m1, d1 = C.gene_stats(X)
    m2, d2 = C.gene_stats(sparse.csr_matrix(X))
    assert np.allclose(m1, m2) and np.allclose(d1, d2)
    assert np.allclose(m1, X.mean(axis=0)) and np.allclose(d1, (X > 0).mean(axis=0))


def test_pool_rules():
    symbols = ["A", "B", "B", "C", "Xist", "D", "E"]
    mean = np.array([1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0])
    det = np.full(7, 0.5)
    assert C.pool(symbols, mean, det, exclude={"Xist"}) == [0, 5, 6]  # B twice, C without counts, Xist excluded
    assert C.pool(symbols, mean, det, exclude={"Xist"}, mapped={"A", "E"}) == [0, 6]


def test_nearest_on_mean_and_detection_ties_by_order():
    mean = np.array([1.0, 1.0, 10.0, 1.0, 0.1])
    det = np.array([0.5, 0.5, 0.5, 0.9, 0.5])
    near = C.nearest(mean, det, [1, 2, 3, 4], 0, 2)
    assert near[0] == 1  # identical
    d = C.distances(mean, det, [1, 2, 3, 4], 0)
    assert d[0] == 0 and d[1] == pytest.approx(d[3])  # 10x and 0.1x are equally far on log10
    assert C.nearest(mean, det, [1, 2, 3, 4], 0, 3)[1:] == sorted(C.nearest(mean, det, [1, 2, 3, 4], 0, 3)[1:])


def _string(tmp_path, genes, links):
    ids = {g: f"10090.P{i:05d}" for i, g in enumerate(genes)}
    with gzip.open(tmp_path / C.STRING_FILES[0], "wt") as fh:
        fh.write("#string_protein_id\tpreferred_name\tprotein_size\tannotation\n")
        for g, p in ids.items():
            fh.write(f"{p}\t{g}\t100\tx\n")
    with gzip.open(tmp_path / C.STRING_FILES[1], "wt") as fh:
        fh.write("protein1 protein2 combined_score\n")
        for a, b, s in links:
            fh.write(f"{ids[a]} {ids[b]} {s}\n{ids[b]} {ids[a]} {s}\n")
    return C.string_names(tmp_path / C.STRING_FILES[0]), tmp_path / C.STRING_FILES[1]


def test_draws(tmp_path):
    rng = np.random.default_rng(1)
    genes = ["Smad3", "Col1a1", "Actb", "Xist", *C.Y_GENES] + [f"G{i}" for i in range(400)]
    mean = np.exp(rng.normal(0, 1, len(genes)))
    det = rng.uniform(0.05, 0.9, len(genes))
    names, links = _string(tmp_path, genes, [("G1", "G2", 900), ("G3", "G4", 300)])
    s = C.draw_pairs(genes, mean, det, list(C.PAIR), names, links, k=50, nearest_k=30)
    pairs = [frozenset(p) for p in s["controls"]]
    assert len(pairs) == len(set(pairs)) == 50 and all(len(p) == 2 for p in pairs)
    assert not {g for p in pairs for g in p} & ({"Smad3", "Col1a1", "Actb", "Xist"} | set(C.Y_GENES))
    assert s["controls"] == C.draw_pairs(genes, mean, det, list(C.PAIR), names, links, k=50, nearest_k=30)["controls"]
    g = C.draw_genes(genes, mean, det, "Xist", k=40, nearest_k=80)
    assert len(g["controls"]) == len({c[0] for c in g["controls"]}) == 40
    st = C.draw_sets(genes, mean, det, list(C.Y_GENES), k=60, nearest_k=30)
    assert len(st["controls"]) == 60 and all(len(set(x)) == 4 for x in st["controls"])
    assert len({frozenset(x) for x in st["controls"]}) == 60


def test_relabelings_and_exact_correlation():
    R4 = R.relabelings(5, 2)
    assert R4.shape == (10, 5) and (R4.sum(axis=1) == 2).all() and len({tuple(r) for r in R4}) == 10
    x, y = [1.0, 2.0, 3.0, 4.0], [2.0, 1.0, 4.0, 3.0]
    c = R.corr_perm(x, y)
    r_obs = np.corrcoef(x, y)[0, 1]
    brute = np.mean([abs(np.corrcoef(x, np.array(y)[list(p)])[0, 1]) >= abs(r_obs) - 1e-12
                     for p in itertools.permutations(range(4))])
    assert c["pearson"]["r"] == pytest.approx(r_obs) and c["pearson"]["p"] == pytest.approx(brute) and c["orderings"] == 24


def test_t_of_is_second_minus_first_over_its_se():
    vals = {"a": 1.0, "b": 1.2, "c": 0.9, "d": 2.0, "e": 2.4, "f": 2.1}
    d, se, t = R.t_of(vals, sorted(vals), second={"d", "e", "f"})
    te = stats.nested_t(np.array([1.0, 1.2, 0.9, 2.0, 2.4, 2.1]), np.array([1, 1, 1, 0, 0, 0]), np.zeros(6, int))
    assert d == pytest.approx(-te.est) and d > 0 and se == pytest.approx(te.se) and t == pytest.approx(d / se)


def test_outcome3():
    assert R.outcome3(0, 1000)[0] == "specific"
    assert R.outcome3(200, 1000)[0] == "not specific"
    assert R.outcome3(45, 1000)[0] == "undecided"  # the interval reaches across 0.05


def test_reduced_values_are_the_full_matrix_values():
    rng = np.random.default_rng(2)
    X = rng.poisson(1.0, (120, 8)).astype(np.int32)
    genes = [f"G{i}" for i in range(8)]
    tot = X.sum(axis=1, dtype=np.int64)
    rows = np.arange(30, 90)
    full = SimpleData(X[rows].astype(float), pd.DataFrame(index=np.arange(len(rows))), genes)
    assert R.value("mi_3bin", X, genes, ["G1", "G2"], tot, rows) == metrics.mi_3bin(full, gene_a="G1", gene_b="G2")
    assert R.value("norm_pearson", X, genes, ["G1", "G2"], tot, rows) == metrics.norm_pearson(full, gene_a="G1", gene_b="G2")


def _cells():
    """Two types in 4 young and 4 old male units; type B only in two old units (does not qualify)."""
    rng = np.random.default_rng(3)
    rows = []
    for u, age in [(f"y{i}", "young") for i in range(4)] + [(f"o{i}", "old") for i in range(4)]:
        for ct, n in (("fibroblast of lung", 30), ("B cell", 25 if u in ("o0", "o1") else 5), ("T cell", 22)):
            rows += [dict(mouse=u, age=age, sex="male", cell_type=ct, tissue="lung" if ct != "T cell" else "spleen")] * n
    obs = pd.DataFrame(rows)
    genes = ["Smad3", "Col1a1", "G2", "G3"]
    X = rng.poisson(1.5, (len(obs), len(genes))).astype(np.int32)
    return X, genes, obs


def test_composition_rules():
    X, genes, obs = _cells()
    res = R.composition(X, genes, obs, None, lambda *m: None)
    m = res["sexes"]["male"]
    assert set(m["within"]["mi_3bin"]["types"]) == {"fibroblast of lung", "T cell"}  # B cell: >= 20 in 2 old units only
    assert res["fibroblast_like_labels"] == {"fibroblast of lung": 240}
    assert all(0 < v < 1 for v in m["fibroblast_like_share"].values())
    # the combined Delta is the mean of the types' Deltas; its p is over the 70 relabelings
    w = m["within"]["mi_3bin"]
    assert w["combined"]["delta"] == pytest.approx(np.mean([v["delta"] for v in w["types"].values()]))
    assert w["combined"]["relabelings"] == 70 and 0 < w["combined"]["p"] <= 1
    assert m["cell_type"]["relabelings"] == 70 and 0 <= m["cell_type"]["tvd"] <= 1
    assert "prediction" in res and "outcome" in res


def test_the_second_runs_values_by_their_structure():
    """real_pilot.json of the second run: cases -> case -> claims -> [{id, b: {values: {unit: float}}}]."""
    X, genes, obs = _cells()
    tot = X.sum(axis=1, dtype=np.int64)
    unit = obs["mouse"].to_numpy()
    mine = {u: R.value("mi_3bin", X, genes, list(C.PAIR), tot, np.where(unit == u)[0]) for u in sorted(set(unit))}
    run2 = {"cases": {"case1": {"claims": [{"id": "1a-floor", "b": {"values": {}}},
                                           {"id": "1a", "b": {"values": mine, "delta": 0.0}}]},
                      "case2": {"claims": [{"id": "2a", "a": {}}]}}}
    assert R.run2_values(run2, "case1", "1a") == mine
    assert R.run2_values(run2, "case2", "2a") is None and R.run2_values(run2, "case3", "1a") is None
    assert R.run2_values(None, "case1", "1a") is None
    m = R.composition(X, genes, obs, R.run2_values(run2, "case1", "1a"), lambda *a: None)["sexes"]["male"]
    assert m["mi_vs_run2_units"] == 8 and m["mi_vs_run2_max_abs_diff"] == 0.0
