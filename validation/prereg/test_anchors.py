"""Anchor R1's loading of a sparse background (DEVIATIONS.md, D1), on synthetic files only.

In v1's run, B2 (stored sparse) reached SimpleData as a scipy matrix and R1 raised ValueError
before the engine ran. The dry runs' backgrounds are dense, so no run took this path. These tests
take it with small synthetic files: never with B2 itself, whose R1 verdict nobody has seen.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import anchors  # noqa: E402
import panel as P  # noqa: E402
from metric_autopsy import SimpleData  # noqa: E402


def _background(tmp: Path, *, dense: bool = False, mice: int = 4, cells: int = 30, fillers: int = 30, seed: int = 0):
    """A B2-like file: mice of both sexes in one age, Xist in female cells, the Y genes in male ones
    and filler genes of increasing expression; written as panel.save_npz writes a selection."""
    rng = np.random.default_rng(seed)
    genes = ["Xist", *anchors.Y_GENES, *[f"G{j}" for j in range(fillers)]]
    rows, obs = [], []
    for sex in ("female", "male"):
        for m in range(mice):
            for _ in range(cells):
                lam = np.r_[3.0 if sex == "female" else 0.0, np.full(len(anchors.Y_GENES), 1.0 if sex == "male" else 0.0),
                            np.linspace(0.1, 6.0, fillers)]
                rows.append(rng.poisson(lam))
                obs.append(dict(sex=sex, development_stage="3-month-old stage", donor_id=f"{sex[0]}{m}"))
    X = np.asarray(rows, dtype=np.int32)
    path = tmp / "B2.npz"
    P.save_npz(path, X if dense else sparse.csr_matrix(X), pd.DataFrame(obs), genes)
    spec = tmp / "backgrounds.json"
    spec.write_text(json.dumps({"B2": {"file": "B2.npz", "sha256": P.sha256(path), "counts": "X", "donor": "donor_id"}}))
    return spec, X, genes


def test_a_sparse_background_failed_before_d1_and_loads_dense_after(tmp_path):
    spec, X, genes = _background(tmp_path)
    with pytest.raises(ValueError):  # the path of v1's run: a scipy matrix into SimpleData
        SimpleData(sparse.csr_matrix(X), pd.DataFrame(index=range(len(X))), genes)
    got, kept, obs = anchors._load(spec, "B2", tmp_path, named=("Xist", *anchors.Y_GENES))
    assert isinstance(got, np.ndarray) and got.shape == X.shape  # fewer genes than N_GENES: all kept
    assert kept == genes and np.array_equal(got, X)
    SimpleData(got, obs, kept)


def test_the_cut_follows_the_panel_rule(tmp_path, monkeypatch):
    """N_GENES in all: the named genes first, then the most expressed by mean count (ties by
    column), in column order."""
    spec, X, genes = _background(tmp_path, fillers=12)
    monkeypatch.setattr(anchors.P, "N_GENES", 8)
    got, kept, _ = anchors._load(spec, "B2", tmp_path, named=("XIST", "Uty"))  # case-insensitive, as _cols
    named = [genes.index("Xist"), genes.index("Uty")]
    by_mean = [j for j in np.argsort(-X.mean(axis=0), kind="mergesort") if j not in named]
    want = sorted(named + by_mean[:8 - len(named)])
    assert kept == [genes[j] for j in want]
    assert np.array_equal(got, X[:, want])


def test_a_dense_background_is_unchanged(tmp_path, monkeypatch):
    spec, X, genes = _background(tmp_path, dense=True)
    monkeypatch.setattr(anchors.P, "N_GENES", 4)  # no cut for a dense file (B3)
    got, kept, _ = anchors._load(spec, "B2", tmp_path, named=("Xist",))
    assert kept == genes and np.array_equal(got, X)


def test_r1_runs_on_a_sparse_background(tmp_path):
    """The whole of R1 with the frozen engine on a synthetic sparse B2: three claims, each with a
    verdict and its allowed set (which verdict is not the point: the data are synthetic)."""
    spec, _, _ = _background(tmp_path, mice=4, cells=25, fillers=20)
    out = anchors.r1(spec, tmp_path)
    assert [r["claim"] for r in out] == ["R1 Xist", "R1 Y genes", "R1 sham (female vs female)"]
    for r in out:
        assert "error" not in r and r["verdict"] and r["label"]
        assert set(r["fields"]) == {"metric_validity", "design_adequacy", "effect", "replication"}
    assert out[0]["allowed"] == ["SUPPORTED"]  # 4 mice of each sex: the replicate rule gives a verdict
