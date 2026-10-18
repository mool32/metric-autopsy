"""Tests of the third diagnostics' rules (real_pilot_5.py, real_pilot_5_extract.py, real_pilot_5_spec.py), on
fabricated tables and random counts only. The workflow real-pilot-5.yml runs them in its dry run, with the
engine of the pull request's head:

    python -m pytest -q -p no:cacheprovider -c /dev/null validation/exploratory/test_real_pilot_5.py
"""
from __future__ import annotations

import gzip
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import real_pilot_5 as R  # noqa: E402
import real_pilot_5_extract as E  # noqa: E402
import real_pilot_5_spec as S  # noqa: E402
import real_pilot_select as RS  # noqa: E402
import real_pilot_diagnostics as D  # noqa: E402


# --------------------------------------------------------------------------- #
# the registered lists
# --------------------------------------------------------------------------- #
def test_the_lists_are_the_plans():
    assert S.FIB_GENES == ("Col1a1", "Col1a2", "Col3a1", "Col5a1", "Dcn", "Fbln1", "Fbln2", "Loxl1", "Lum", "Sparc")
    assert set(S.FIB_GENES) == set(S.MUHL_ECM) | set(S.OWNER_EXAMPLES)
    assert S.INDEX_GENES == ("Col1a2", "Col3a1", "Col5a1", "Dcn", "Fbln1", "Fbln2", "Loxl1", "Lum")
    assert "Col1a1" in S.INDEX_PLUS and "Sparc" not in S.INDEX_PLUS
    assert S.SEX_GENES == D.SEX_GENES and S.PAIR == D.PAIR and S.Y_GENES == D.Y_GENES
    assert len(S.MARSH_SET) == 25 and set(S.MARSH_PAIR) <= set(S.MARSH_SET)
    arms = [a for _, a, _ in S.MARSH_SAMPLES.values()]
    assert all(arms.count(a) == 3 for a in S.ARMS) and len(S.MARSH_SAMPLES) == 12
    assert {c: v[0] for c, v in S.CONTRASTS.items()} == {"M1": ("ENZ-INHIB", "ENZ-NONE"),
                                                      "M2": ("DNC-NONE", "ENZ-NONE"),
                                                      "M3": ("DNC-NONE", "DNC-INHIB")}


def test_the_acceptance_claims_are_the_plans():
    ids = [c["id"] for c in S.ACCEPT]
    assert ids == ["A1", "A2", "A3", "A4", "A5", "A6", "A7", "A5m", "A6m", "A7m"]
    by = {c["id"]: c for c in S.ACCEPT}
    assert by["A3"]["males"] and not by["A4"].get("males") and by["A4"]["within"] == ["sex"]
    assert by["A5"]["metric"] == "norm_pearson" and by["A5m"]["metric"] == "mi_3bin"
    assert by["A1"]["sesoi"] == by["A2"]["sesoi"] == 0.5 and by["A3"]["sesoi"] == by["A5"]["sesoi"] == 0.1
    assert set(S.PARTS) == {"b", "a", "case2", "case1m", "case1", "marsh", "case3a", "case3b", "case3c"}


# --------------------------------------------------------------------------- #
# (b): the exact combined difference and the permutations
# --------------------------------------------------------------------------- #
def test_the_combined_difference_is_the_mean_of_the_types_differences():
    V = np.array([[1.0, 2.0, 3.0, np.nan, 5.0],
                  [0.0, np.nan, 1.0, 1.0, 1.0]])
    young = np.array([True, True, False, False, False])
    d, per = R.combined_delta(V, young)
    assert per[0] == pytest.approx(1.5 - 4.0) and per[1] == pytest.approx(0.0 - 1.0)
    assert d == pytest.approx(np.mean([1.5 - 4.0, -1.0]))


def test_the_exact_p_counts_every_relabeling():
    rng = np.random.default_rng(0)
    V = rng.normal(size=(3, 6))
    young = np.array([True, True, True, False, False, False])
    ex = R.exact_combined(V, young)
    assert ex["relabelings"] == 20
    stats = []
    for pick in itertools.combinations(range(6), 3):
        lab = np.zeros(6, bool)
        lab[list(pick)] = True
        stats.append(R.combined_delta(V, lab)[0])
    assert ex["p"] == pytest.approx(np.mean(np.abs(stats) >= abs(ex["delta"]) - 1e-12))
    V[:, :3] += 10.0  # the young far above: the smallest p, 2 of 20 (both extremes)
    assert R.exact_combined(V, young)["p"] == pytest.approx(0.1)


def test_the_association_is_within_cells_and_its_permutations_stay_within_them():
    rng = np.random.default_rng(1)
    cells = np.repeat(["t1|young", "t1|old", "t2|young", "t2|old"], 6)
    shift = np.repeat([5.0, -5.0, 3.0, -3.0], 6)
    x = rng.normal(size=24) + shift
    y = rng.normal(size=24) - 7 * shift  # the cells' means are anti-aligned: no within-cell association
    r = R.perm_r(x, y, cells, rng, 199)
    assert abs(r["r"]) < 0.6 and r["p"] > 0.01
    y2 = 2 * (x - shift) + rng.normal(scale=0.1, size=24) - 7 * shift  # strongly associated within cells
    r2 = R.perm_r(x, y2, cells, rng, 199)
    assert r2["r"] > 0.9 and r2["p"] == pytest.approx(1 / 200)


def test_fibroblast_like_follows_the_second_diagnostics_label_rule():
    assert R.fib_like("fibroblast of cardiac tissue") and R.fib_like("Mesenchymal stem cell") and R.fib_like("stromal cell")
    assert not R.fib_like("B cell") and not R.fib_like("endothelial cell")


# --------------------------------------------------------------------------- #
# (b4): the metadata's columns and GEO's ids
# --------------------------------------------------------------------------- #
def test_a_batch_separates_the_ages_a_plate_is_nested_and_an_id_is_left_out():
    units = ["y1", "y2", "y3", "o1", "o2", "o3"]
    rows = []
    for k, u in enumerate(units):
        for c in range(4):
            rows.append(dict(unit=u, batch="TM" if u[0] == "y" else "TMS", plate=f"P{k}{c % 2}", mouse_copy=u,
                             method="facs", depth=float(100 + 10 * k + c + (50 if u[0] == "y" else 0))))
    fr = pd.DataFrame(rows)
    sex_of = {u: "male" for u in units}
    age_of = {u: "young" if u[0] == "y" else "old" for u in units}
    out = R.column_check(fr, "unit", sex_of, age_of, set(), lambda *a: None, "test")
    by = {r["column"]: r for r in out}
    assert by["batch"]["outcome"] == "separates the ages (shared within an age)"
    assert by["plate"]["outcome"] == "nested within units"
    assert by["method"]["outcome"] == "one value"
    assert "mouse_copy" not in by  # one value per unit, one unit per value: an id
    assert by["depth"]["kind"] == "numeric" and by["depth"]["young_minus_old"] > 0 and by["depth"]["relabelings"] == 20


def test_unit_ids_match_whole_tokens_with_any_separator_and_either_order():
    pats = R.unit_pattern("3_8_M")
    hit = lambda s: any(p.search(s) for p in pats)  # noqa: E731
    for s in ("mouse 3_8_M plate", "3-8-M", "id: 3.8.m", "3-M-8 lung", "(3_M_8)"):
        assert hit(s), s
    for s in ("13_8_M", "3_8_MX", "3_18_M", "3_8", "23-M-8"):
        assert not hit(s), s


def test_the_statement_on_tabula_muris_is_scored_by_gse109774(tmp_path):
    def soft(path, ids):
        lines = ["^SERIES = GSE109774"]
        for i, m in enumerate(ids):
            lines += [f"^SAMPLE = GSM{i}", f"!Sample_title = plate {i}", f"!Sample_characteristics_ch1 = mouse: {m}"]
        path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))
    (tmp_path / "GSE132042_family.soft.gz").write_bytes(gzip.compress(b"^SERIES = GSE132042\n!Series_summary = x\n"))
    soft(tmp_path / "GSE109774_family.soft.gz", S.YOUNG_UNITS)
    assert R.geo_check(tmp_path, lambda *a: None)["statement_score"].startswith("holds")
    soft(tmp_path / "GSE109774_family.soft.gz", [*S.YOUNG_UNITS, "24_58_M"])
    assert R.geo_check(tmp_path, lambda *a: None)["statement_score"].startswith("does not hold")
    soft(tmp_path / "GSE109774_family.soft.gz", ["10X_P4_0"])
    assert R.geo_check(tmp_path, lambda *a: None)["statement_score"].startswith("undecided")


# --------------------------------------------------------------------------- #
# Marsh et al. 2022: the fabricated GEO files, their checks and their merge
# --------------------------------------------------------------------------- #
def test_the_marsh_files_are_checked_and_merged(tmp_path):
    geo, out = tmp_path / "geo", tmp_path / "out"
    out.mkdir()
    E.fabricate_geo(geo)
    assert E.check_marsh_soft(geo / "GSE152183_family.soft.gz") == []
    args = type("A", (), dict(geo=geo, out=out, synthetic=True))()
    res = {}
    E.run_marsh(args, lambda *a: None, res)
    m = res["marsh"]
    assert m["set_missing"] == [] and m["left_out"] == 2  # the fabricated symbol that names two features
    X, genes, obs = R.load_npz(out / "marsh.npz")
    assert X.shape[0] == m["cells"] and sorted(obs["mouse"].unique()) == sorted(S.MARSH_SAMPLES)
    assert set(obs["cell_type"]) == {S.MARSH_TYPE} and "Dup2" not in genes and {"Fos", "Jun"} <= set(genes)
    data, pool = R.load_marsh(out)
    assert data.X.shape[1] <= R.N_GENES and pool.X.shape == X.shape and {"Fos", "Jun"} <= set(data.var_names)


def test_a_registered_table_that_geo_contradicts_stops_the_extraction(tmp_path):
    geo = tmp_path / "geo"
    E.fabricate_geo(geo)
    text = gzip.decompress((geo / "GSE152183_family.soft.gz").read_bytes()).decode()
    bad = text.replace("batch: Batch 1", "batch: Batch 3", 1)
    (geo / "GSE152183_family.soft.gz").write_bytes(gzip.compress(bad.encode()))
    problems = E.check_marsh_soft(geo / "GSE152183_family.soft.gz")
    assert len(problems) == 1 and "batch" in problems[0]


def test_the_fabricated_census_keeps_real_pilot_selects_counts():
    base, ours = RS.SyntheticCensus(), E.SyntheticCensus()
    n0 = base._X["mus_musculus"].shape[1]
    assert np.array_equal(base._X["mus_musculus"], ours._X["mus_musculus"][:, :n0])
    assert np.array_equal(base._X["homo_sapiens"], ours._X["homo_sapiens"])
    added = list(ours._var["mus_musculus"]["feature_name"])[n0:]
    assert set(added) == set(S.FIB_GENES) - set(RS.TMS_GENES)


def test_the_d1_cut_keeps_the_named_genes_then_the_most_expressed():
    rng = np.random.default_rng(3)
    X = rng.poisson(np.linspace(0.01, 5, 3000), size=(40, 3000))
    genes = [f"g{i}" for i in range(3000)]
    dense, kept = R.d1_cut(X, genes, ("g0", "g1"))
    assert len(kept) == R.N_GENES and kept[:2] == ["g0", "g1"] and "g2999" in kept
    assert kept == sorted(kept, key=genes.index)


# --------------------------------------------------------------------------- #
# the scores
# --------------------------------------------------------------------------- #
def _rec(label, status="DETECTED", domain=None, mix=None, outside=None):
    return dict(label=label, fields=dict(effect=dict(status=status)),
                domain=domain, mixture=dict(corrected=mix) if mix else None)


def test_the_criteria_are_scored_as_the_plan_fixes_them():
    dom = dict(types=["fibroblast"], empty=False, outside=dict(corrected=dict(detected=True, effect=0.02)))
    recs = {"A1": _rec("SUPPORTED"), "A2": _rec("SUPPORTED"), "A4": _rec("INCONCLUSIVE"),
            "A5": _rec("SUPPORTED"), "A7": _rec("NO DETECTABLE EFFECT"),
            "A3": _rec("INCONCLUSIVE", status="INCONCLUSIVE", domain=dom, mix=dict(detected=True, effect=0.019)),
            **{f"A8-3a{i}": _rec("EFFECT BELOW SESOI") for i in range(3)}}
    sc = R.score_claims(recs)
    assert all(sc[k]["met"] for k in ("A1", "A2", "A3", "A4", "A5", "A7"))
    assert sc["A8"]["strict"] and sc["A8"]["lenient"]
    recs["A3"] = _rec("INCONCLUSIVE", status="INCONCLUSIVE", domain=dom, mix=dict(detected=False, effect=0.019))
    recs["A8-3a0"] = _rec("INCONCLUSIVE")
    sc = R.score_claims(recs)
    assert sc["A3"]["met"] is False and sc["A3"]["drop_outside"] is False and sc["A3"]["no_effect_in_domain"]
    assert not sc["A8"]["strict"] and sc["A8"]["lenient"]
    recs["A8-3a1"] = _rec("SUPPORTED")
    recs["A5"] = _rec("INCONCLUSIVE")
    sc = R.score_claims(recs)
    assert not sc["A8"]["lenient"] and sc["A8"]["supported"] == ["A8-3a1"] and sc["A5"]["met"] is False
    empty = dict(types=[], empty=True, outside=None)
    recs["A3"] = _rec("INCONCLUSIVE", status="INCONCLUSIVE", domain=empty, mix=dict(detected=True, effect=0.02))
    assert R.score_claims(recs)["A3"]["met"] is False


def test_the_prediction_on_foreign_counts_needs_both_parts():
    b = dict(b2=dict(F=dict(non_fibroblast=dict(delta=0.3, p=0.016))), b3_i=dict(r=0.4, p=0.003))
    assert R.score_b(b)["score"] == "holds"
    b["b3_i"] = dict(r=0.4, p=0.2)
    assert R.score_b(b)["score"] == "partly holds"
    b["b2"]["F"]["non_fibroblast"] = dict(delta=-0.3, p=0.016)
    assert R.score_b(b)["score"] == "does not hold"


def test_the_verdicts_label():
    assert R.label_of("SUPPORTED (provisional until replicated) [x]") == "SUPPORTED"
    assert R.label_of("SUPPORTED — replicated") == "SUPPORTED"
    assert R.label_of("EFFECT BELOW SESOI — effect +0.1") == "EFFECT BELOW SESOI"
    assert R.label_of("NOT SUPPORTED — metric invalid: x") == "NOT SUPPORTED"
