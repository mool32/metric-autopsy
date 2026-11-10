"""Level B, step S6 (`plan_s6.md`): the controls' rules, on synthetic tables and the fabricated run."""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import controls as C  # noqa: E402


def _counts(rows):
    return pd.DataFrame(rows, columns=["organism", *C.COUNTED, "cells"]).astype({c: str for c in C.COUNTED})


def _unit(org, ds, tissue, ct, donor, assay, susp, cells, sex="female", stage="50-year-old stage", disease="normal"):
    return [org, ds, tissue, ct, donor, assay, susp, disease, sex, stage, cells]


DATASETS = pd.DataFrame([
    dict(dataset_id="a", collection_id="c1", collection_name="One", collection_doi="10.0/one", dataset_title="A"),
    dict(dataset_id="b", collection_id="c1", collection_name="One", collection_doi="10.0/one", dataset_title="B"),
    dict(dataset_id="x", collection_id="c2", collection_name="SEA-AD: atlas", collection_doi="10.0/x", dataset_title="X"),
    dict(dataset_id="y", collection_id="c3", collection_name="Other", collection_doi="10.1038/nbt.4042", dataset_title="Y"),
    dict(dataset_id="37a17b78-4864-4a42-b67b-31c00962795a", collection_id="c4", collection_name="Cohort",
         collection_doi="", dataset_title="MSSM_Cohort"),
])


def test_the_constants_are_the_plans_and_the_engines_domain_rule():
    from metric_autopsy import domain
    assert C.DOMAIN_DETECTION == domain.DOMAIN_DETECTION
    assert (C.MIN_DONORS, C.MIN_CELLS, C.MAX_DONORS, C.MAX_CELLS, C.MAX_TRIED, C.N_NEGATIVE) == (6, 200, 200, 400, 10, 3)
    assert (C.SEED_CELLS, C.SEED_DIRECTIONS, C.SEED_MATCHING) == (20261110, 20261111, 20261112)
    assert len(C.ISGS) == 12 and C.SET_LEVELS == ("high", "medium")


def test_exclusions_by_dataset_name_and_doi():
    ex = C.exclusions(DATASETS)
    assert set(ex) == {"c2", "c3", "c4"}
    assert "SEA-AD" in ex["c2"] and "Kang" in ex["c3"] and "B1" in ex["c4"]


def test_a_candidate_needs_six_donors_of_200_cells_alike_in_both_conditions():
    rows = []
    for k in range(6):  # c1: six donors, cells (dataset a) and nuclei (dataset b)
        rows.append(_unit("homo_sapiens", "a", "heart", "fibroblast", f"D{k}", "10x 3' v3", "cell", 250))
        rows.append(_unit("homo_sapiens", "b", "heart", "fibroblast", f"D{k}", "10x 3' v3", "nucleus", 210 + k))
    rows.append(_unit("homo_sapiens", "b", "heart", "fibroblast", "D0", "10x 3' v3", "nucleus", 30, stage="51-year-old stage"))
    rows.append(_unit("homo_sapiens", "a", "heart", "fibroblast", "D9/D8", "10x 3' v3", "cell", 900))  # composite id
    for k in range(5):  # the same with five donors: no candidate
        rows.append(_unit("homo_sapiens", "a", "heart", "endothelial cell", f"D{k}", "10x 3' v3", "cell", 250))
        rows.append(_unit("homo_sapiens", "b", "heart", "endothelial cell", f"D{k}", "10x 3' v3", "nucleus", 250))
    for k in range(8):  # c2 (excluded): more donors
        rows.append(_unit("homo_sapiens", "x", "brain", "astrocyte", f"S{k}", "10x 3' v3", "nucleus", 300))
        rows.append(_unit("homo_sapiens", "x", "brain", "astrocyte", f"S{k}", "10x multiome", "nucleus", 300))
    cand, donors = C.negative_candidates(_counts(rows), DATASETS)
    assert list(cand["collection_id"]) == ["c2", "c1"]
    top = cand.set_index("collection_id")
    assert top.loc["c2", "excluded"] and top.loc["c2", "rank"] == 0
    assert top.loc["c1", "rank"] == 1 and top.loc["c1", "donors"] == 6 and top.loc["c1", "mismatched"] == 0
    assert top.loc["c1", "condition_a"] == "10x 3' v3 | cell" and top.loc["c1", "condition_b"] == "10x 3' v3 | nucleus"
    assert top.loc["c1", "median_cells"] == np.median([250] * 6 + [210 + k for k in range(6)])
    assert list(C.first_per_collection(cand)["collection_id"]) == ["c1"]
    # a donor whose sex differs between the conditions leaves five: no candidate
    rows2 = [r if not (r[4] == "D1" and r[6] == "nucleus") else r[:8] + ["male"] + r[9:] for r in rows]
    cand2, _ = C.negative_candidates(_counts(rows2), DATASETS)
    assert "c1" not in set(cand2["collection_id"])


def test_the_sample_is_capped_per_donor_and_condition_and_reproducible():
    rng = np.random.default_rng(0)
    cells = pd.DataFrame(dict(soma_joinid=rng.permutation(5000), donor_id=np.repeat(["d1", "d2"], 2500),
                              condition=np.tile(np.repeat(["A", "B"], 1250), 2)))
    ids = C.sample_units(cells, ["d2", "d1"])
    assert len(ids) == 4 * C.MAX_CELLS and np.all(np.diff(ids) > 0)
    assert np.array_equal(ids, C.sample_units(cells, ["d1", "d2"]))


def test_the_claims_take_every_pool_pair_two_sets_and_seeded_directions():
    pool = [dict(level=lv, index=i, pair=[f"G{i}a", f"G{i}b"], r=0.5, neg_pair=["n1", "n2"], pos_pair=["p1", "p2"])
            for i, lv in enumerate(["high"] * 4 + ["medium"] * 4 + ["low"] * 4)]
    claims = C.negative_claims(dict(pool=pool), 1)
    assert [c["kind"] for c in claims] == ["pair"] * 12 + ["set", "set"]
    assert claims[-2]["genes"] == [g for i in range(4) for g in (f"G{i}a", f"G{i}b")] and claims[-1]["level"] == "medium"
    assert claims == C.negative_claims(dict(pool=pool), 1)
    assert [c["higher_in"] for c in claims] != [c["higher_in"] for c in C.negative_claims(dict(pool=pool), 2)]
    assert {c["higher_in"] for c in claims} == {"A", "B"}


def test_the_sex_call():
    assert C.sex_call(0.6, 0.0) == "female" and C.sex_call(0.0, 0.4) == "male"
    assert C.sex_call(0.3, 0.2) == "unclear" and C.sex_call(0.02, 0.01) == "unclear"


def test_kangs_names_and_table():
    assert C.kang_lane_names(["a-1", "b-1"], ["b-1", "c-1"]) == ["b-11", "c-1"]
    raw = gzip.compress(b"tsne1\ttsne2\tind\tstim\tcluster\tcell\tmultiplets\nA-1\t0\t0\t101\tctrl\t1\tB cells\tsinglet\n")
    t = C.kang_cells(raw)
    assert list(t.columns) == ["name", "tsne1", "tsne2", "ind", "stim", "cluster", "cell", "multiplets"]
    assert t.loc[0, "name"] == "A-1" and t.loc[0, "ind"] == "101"
    m = C.read_mtx(gzip.compress(b"%%MatrixMarket matrix coordinate integer general\n%\n3 2 2\n1 1 5\n3 2 7\n"))
    assert m.shape == (3, 2) and m[0, 0] == 5 and m[2, 1] == 7
    with pytest.raises(ValueError):
        C.read_mtx(gzip.compress(b"%%MatrixMarket matrix coordinate integer general\n3 2 3\n1 1 5\n"))


def test_the_dry_run_end_to_end(tmp_path):
    assert C.main(["--out", str(tmp_path), "--synthetic"]) == 0
    sel = json.loads((tmp_path / "selection.json").read_text())
    claims = json.loads((tmp_path / "claims.json").read_text())
    assert sel["synthetic"] and claims["synthetic"]
    neg = sel["negative"]
    assert neg["candidates"] == 4 and neg["excluded"] == 1 and neg["tried"] == 3
    assert [c["collection_id"] for c in neg["chosen"]] == ["c-heart", "c-mouse"]  # the kidney's counts are not integers
    tried = pd.read_csv(tmp_path / "tried.tsv", sep="\t")
    assert list(tried["outcome"]) == ["chosen: N1", "fails: the counts are not non-negative integers", "chosen: N2"]
    n1 = claims["negative"][0]
    assert n1["groups"] == {"A": "10x 3' v3 | cell", "B": "10x 3' v3 | nucleus"} and not n1["sex_conflicts"]
    assert sum(c["kind"] == "pair" for c in n1["claims"]) == 3 * n1["pool"]["pairs_per_level"]
    assert [c["kind"] for c in n1["claims"]][-2:] == ["set", "set"]
    assert len((tmp_path / "joinids" / "N1.txt").read_text().split()) == n1["cells"]
    pos = claims["positive"][0]
    assert pos["passed"] and len(pos["donors_in_both"]) == 8 and "IFIT2" not in pos["eligible"]
    assert len(pos["eligible"]) == 11 and pos["unmatched"] and [c["kind"] for c in pos["claims"]].count("pair") == 5
    assert {g for c in pos["claims"] if c["kind"] == "pair" for g in c["genes"]} | set(pos["unmatched"]) == set(pos["eligible"])
    assert (tmp_path / "data" / "K1.npz").exists() and "SYNTHETIC" in (tmp_path / "controls.log").read_text()
