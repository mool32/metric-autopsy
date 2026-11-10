"""Level B, step S5 (`plan_s5.md`): the data checks on synthetic data."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import data_checks as D  # noqa: E402


def test_rule4_is_the_protocols():
    assert D.rule4("nested", 4, 4)["passes"] and D.rule4("nested", 3, 5)["passes"]
    assert not D.rule4("nested", 3, 3)["passes"] and not D.rule4("nested", 3, 4)["passes"]
    assert D.rule4("paired", 6)["passes"] and not D.rule4("paired", 5)["passes"]
    assert D.rule4("nested", 4, 4)["relabelings"] == 70 and D.rule4("paired", 6)["min_p"] == 2 / 64


def test_counts_donors_genes_and_source():
    X = sparse.csr_matrix(np.array([[0, 2], [1, 0], [3, 1.5]]))
    assert not D.counts_check(X)["raw"] and D.counts_check(np.array([[0, 2], [1, 0]]))["raw"]
    obs = pd.DataFrame(dict(donor=["a", "a", "b", "c", "d/e", "f"], group=["x", "y", "x", "y", "x", "y"],
                            batch=["1", "1", "1", "2", "1", "2"]))
    d = D.donors_check(obs, "donor", "group", ("x", "y"))
    assert d["design"] == "mixed" and d["composite"] == ["d/e"] and not d["rule4"]["passes"]
    nested = obs[obs["donor"] != "a"]
    assert D.donors_check(nested, "donor", "group", ("x", "y"))["design"] == "nested"
    paired = pd.DataFrame(dict(donor=[f"p{k}" for k in range(6)] * 2, group=["x"] * 6 + ["y"] * 6))
    assert D.donors_check(paired, "donor", "group", ("x", "y"))["rule4"]["passes"]
    assert D.genes_check(["ACE2", "Tmprss2"], ["ace2", "TMPRSS2", "XIAP"])["missing"] == ["XIAP"]
    src = D.source_check(obs, "group", ("x", "y"), "batch")
    assert src["verdict"] == "partly" and src["shared"] == ["1"]
    conf = pd.DataFrame(dict(group=["x", "x", "y", "y"], src=["old", "old", "new", "new"]))
    assert D.source_check(conf, "group", ("x", "y"), "src")["verdict"] == "fails"
    bal = pd.DataFrame(dict(group=["x", "y", "x", "y"], src=["s1", "s1", "s2", "s2"]))
    assert D.source_check(bal, "group", ("x", "y"), "src")["verdict"] == "balanced"


def test_the_sex_check_finds_two_individuals_under_one_id_and_a_label_against_the_genes():
    rng = np.random.default_rng(0)
    genes = ["XIST", "RPS4Y1", "DDX3Y", "G1"]
    rows, donors, labels = [], [], []
    for d, kind, label in (("f1", "f", "female"), ("m1", "m", "male"), ("mix", "both", "female"), ("bad", "m", "female")):
        for _ in range(100):
            x = rng.poisson(3) if kind in ("f", "both") and (kind == "f" or rng.random() < 0.5) else 0
            y = rng.poisson(2) if kind in ("m", "both") and (kind == "m" or x == 0) else 0
            rows.append([x, y, 0, 5])
            donors.append(d)
            labels.append(label)
    X = sparse.csr_matrix(np.array(rows, float))
    obs = pd.DataFrame(dict(donor=donors, sex=labels))
    res = D.sex_check(X, genes, obs, "donor", "homo_sapiens", label_col="sex")
    by = {r["donor"]: r for r in res["donors"]}
    assert res["both_sexes"] and by["f1"]["call"] == "female" and by["m1"]["call"] == "male"
    assert by["mix"]["split"] and by["bad"]["against_label"] and sorted(res["fails"]) == ["bad", "mix"]
