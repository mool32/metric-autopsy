"""Tests of the confirmatory panel's own code: panel.py, oracle.py, score.py, run_panel.py.

They run on simulated backgrounds (no network). The panel and the oracle must not import the
engine; scoring and the oracle need scipy (skipped without it, as in the core-only CI job).
"""
from __future__ import annotations

import ast
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel as P  # noqa: E402

PILOT = dict(sesoi=0.15, key_dose=3.0)


def simulated_background(name="B1", donors=20, cells=220, n_genes=400, seed=0) -> P.Background:
    """Negative-binomial counts with cell-size variation; genes 0-1 and 2-3 coupled pairs."""
    rng = np.random.default_rng(seed)
    mu = np.exp(rng.normal(-0.5, 1.2, n_genes))
    mu[:4] = 4.0
    genes = list(P.G2M_GENES[:20]) + [f"GENE{i}" for i in range(n_genes - 20)]
    genes[:4], genes[20:24] = genes[20:24], genes[:4]  # keep the coupled genes outside the G2M list
    X, donor = [], []
    for d in range(donors):
        size = np.exp(rng.normal(0, 0.35, cells)) * np.exp(rng.normal(0, 0.1))
        m = np.outer(size, mu)
        for a in (0, 2):
            z = rng.normal(0, 1, cells)
            m[:, a] *= np.exp(0.5 * z)
            m[:, a + 1] *= np.exp(0.5 * z)
        X.append(rng.negative_binomial(3.0, 3.0 / (3.0 + m)))
        donor += [f"{name}d{d}"] * cells
    bg = P.Background(np.vstack(X).astype(float), genes, np.asarray(donor), name)
    P.plan_background(bg)
    return bg


@pytest.fixture(scope="module")
def bgs():
    return {"B1": simulated_background("B1"), "B2": simulated_background("B2", donors=12, seed=1)}


# --------------------------------------------------------------------------- #
# the key
# --------------------------------------------------------------------------- #
def test_assign_is_deterministic_complete_and_keyed():
    a, b, c = P.assign(1), P.assign(1), P.assign(2)
    assert a == b and a != c
    assert len(a) == P.n_datasets() == 4350
    assert sum(P.conditions()[e["condition"]].cards for e in a) == 4450  # N4 has two cards
    assert len({e["id"] for e in a}) == len(a) and {e["id"] for e in a} == {e["id"] for e in c}
    counts = Counter((e["condition"], e["variant"]) for e in a)
    assert counts == {(cond.name, v): n for cond in P.CONDITIONS for v, n in cond.variants}
    # the ID says nothing: the same ID holds different conditions under different keys
    assert sum(x["condition"] != y["condition"] for x, y in zip(a, c)) > 0.5 * len(a)
    assert {e["side"] for e in a} == {"A", "B"}


def test_key_null_conditions_have_600_datasets_at_their_key_variant():
    for cond in P.CONDITIONS:
        if cond.key:
            assert dict(cond.variants)[cond.key_variant] == P.KEY_N == 600


# --------------------------------------------------------------------------- #
# datasets and claim cards
# --------------------------------------------------------------------------- #
def _entry(condition, variant, side="A", seed=0, i=0):
    return dict(id=f"T{i}", condition=condition, variant=variant, index=i, side=side, seed=seed)


def test_claim_cards_do_not_reveal_the_condition(bgs):
    """Cards of the norm_pearson conditions differ only in the id and the claimed direction."""
    cards = []
    for cond, variant in (("N1", "null"), ("N2", "c=0.5"), ("N3", "f=0.2"), ("E1", "dose=key"),
                          ("E2", "against"), ("E3", "with")):
        cards.append(P.build(_entry(cond, variant, seed=3), bgs, PILOT)[3][0])
    strip = [{k: v for k, v in c.items() if k not in ("id", "prereg")} for c in cards]
    assert all(s == strip[0] for s in strip)
    assert all({k: v for k, v in c["prereg"].items() if k != "direction"} ==
               {k: v for k, v in cards[0]["prereg"].items() if k != "direction"} for c in cards)
    # nulls get a coin for the direction, effects their true one (the signal side is higher)
    dirs = Counter(P.build(_entry("N1", "null", seed=s, i=s), bgs, PILOT)[3][0]["prereg"]["direction"]
                   for s in range(30))
    assert set(dirs) == {"increase", "decrease"}
    assert P.build(_entry("E1", "dose=key", side="A"), bgs, PILOT)[3][0]["prereg"]["direction"] == "decrease"
    assert P.build(_entry("E1", "dose=key", side="B"), bgs, PILOT)[3][0]["prereg"]["direction"] == "increase"


def test_designs_and_artifacts_are_planted_as_specified(bgs):
    X, obs, genes, _ = P.build(_entry("N1", "null", seed=1), bgs, PILOT)
    per_group = obs.groupby("group")["donor"].nunique()
    assert per_group["A"] == per_group["B"] == P.DONORS_PER_GROUP
    assert X.shape[0] == 2 * P.DONORS_PER_GROUP * P.CELLS_PER_DONOR
    X, obs, _, _ = P.build(_entry("N2", "c=0.5", side="B", seed=1), bgs, PILOT)
    tot = obs.groupby("group")["total_counts"].mean()
    assert tot["B"] / tot["A"] == pytest.approx(0.5, abs=0.08)
    X, obs, _, _ = P.build(_entry("N3", "f=0.4", side="A", seed=1), bgs, PILOT)
    det = (X > 0).mean(axis=1)
    g = np.asarray(obs["group"])
    assert det[g == "A"].mean() / det[g == "B"].mean() == pytest.approx(0.6, abs=0.08)
    X, obs, _, cards = P.build(_entry("N4", "3v3", seed=1), bgs, PILOT)
    assert obs.groupby("group")["donor"].nunique().tolist() == [3, 3]
    assert [c["replicate_col"] for c in cards] == [None, "donor"]
    X, obs, _, _ = P.build(_entry("N5", "sham", seed=1), bgs, PILOT)
    assert (obs.groupby("donor")["group"].nunique() == 2).all()
    X, obs, _, _ = P.build(_entry("N7", "mice", seed=1), bgs, PILOT)
    assert obs["donor"].str.startswith("B2").all()


def test_injected_coupling_raises_the_pair_correlation_on_its_side_only(bgs):
    import oracle as O
    genes = bgs["B1"].plan["genes"]
    ia, ib = genes.index(bgs["B1"].plan["pair"][0]), genes.index(bgs["B1"].plan["pair"][1])
    X, obs, _, _ = P.build(_entry("E1", "dose=key", side="A", seed=2), bgs, dict(PILOT, key_dose=6.0))
    g = np.asarray(obs["group"])
    ra, rb = O.norm_pearson(X[g == "A"], ia, ib), O.norm_pearson(X[g == "B"], ia, ib)
    assert ra > rb + 0.1
    # the sham thins the other side's pair alike, so the means stay comparable
    ma, mb = X[g == "A"][:, [ia, ib]].mean(), X[g == "B"][:, [ia, ib]].mean()
    assert ma / mb == pytest.approx(1.0, abs=0.15)


def test_backgrounds_load_from_their_spec_filtered_and_planned(tmp_path):
    """The real-data path of panel.py, oracle.py and timing.py (`load_backgrounds`): a background
    named in a backgrounds JSON loads as raw counts, keeps only the filtered cells and gets the
    plan the same counts get in memory; counts that are not integers are refused."""
    import pandas as pd
    bg = simulated_background("B1", donors=6, cells=60, n_genes=300)
    cell_type = np.where(np.arange(len(bg.donor)) % 10 == 0, "other", "fibroblast")
    P.save_npz(tmp_path / "b1.npz", bg.X, pd.DataFrame({"donor_id": bg.donor, "cell_type": cell_type}),
               bg.genes)
    spec = {"B1": {"path": str(tmp_path / "b1.npz"), "donor": "donor_id",
                   "filter": {"cell_type": "fibroblast"}}}
    (tmp_path / "backgrounds.json").write_text(json.dumps(spec))
    got = P.load_backgrounds(tmp_path / "backgrounds.json")["B1"]
    keep = cell_type == "fibroblast"
    assert np.array_equal(got.X, bg.X[keep]) and list(got.donor) == list(bg.donor[keep])
    ref = P.Background(bg.X[keep], bg.genes, bg.donor[keep], "B1")
    P.plan_background(ref)
    assert got.plan == ref.plan and got.plan_cols == ref.plan_cols and got.plan["pair"]
    np.savez_compressed(tmp_path / "b1.npz", X=bg.X + 0.5, genes=np.asarray(bg.genes),
                        obs_donor_id=bg.donor, obs_cell_type=cell_type)
    with pytest.raises(ValueError, match="raw counts"):
        P.load_backgrounds(tmp_path / "backgrounds.json")


def test_the_panel_and_the_oracle_never_import_the_engine():
    for name in ("panel.py", "oracle.py", "score.py", "oc.py"):
        tree = ast.parse((HERE / name).read_text())
        mods = {a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
        mods |= {n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module}
        assert not any(m.split(".")[0] == "metric_autopsy" for m in mods), (name, mods)


# --------------------------------------------------------------------------- #
# the oracle
# --------------------------------------------------------------------------- #
def test_oracle_tests_on_replicate_values():
    pytest.importorskip("scipy")
    import oracle as O
    hi, lo = np.array([5.0, 6, 7, 8]), np.array([1.0, 2, 3, 4])
    assert O.permutation_p(hi, lo) == pytest.approx(2 / 70)
    assert O.permutation_p(np.ones(4), np.ones(4)) == 1.0
    assert O.detects(hi, lo)[0] and O.detects(hi, lo)[1] == pytest.approx(4.0)
    rng = np.random.default_rng(0)
    a, b = rng.normal(0, 0.01, 8), rng.normal(0, 0.01, 8)
    assert O.tost_width(a, b) < 0.05 and O.tost_width(a + 1, b) > 0.9
    assert O.signflip_p(np.array([1.0, 2, 3, 4, 5])) == pytest.approx(2 / 32)


def test_the_pilot_sets_sesoi_and_dose_and_finds_establishable_conditions(bgs):
    pytest.importorskip("scipy")
    import oracle as O
    sesoi = O.choose_sesoi(bgs, n=20)
    assert sesoi in O.SESOI_GRID
    dose = O.choose_dose(bgs, sesoi, n=20)
    assert dose in O.DOSE_GRID
    pilot = dict(sesoi=sesoi, key_dose=dose)
    est = O.establishability(bgs, pilot, n=20)
    assert est["E1:dose=key"]["establishable"] and est["N6b:constant"]["establishable"]
    assert est["N4:3v3"]["establishable"] is False
    assert est["N1:null"]["power"] >= 0.8  # the SESOI was chosen on these pilot seeds for 0.9
    assert est["E1:dose=0.25"]["power"] <= est["E1:dose=key"]["power"]


# --------------------------------------------------------------------------- #
# scoring
# --------------------------------------------------------------------------- #
def _verdicts(entries, fn):
    out = {}
    for e in entries:
        cond = P.conditions()[e["condition"]]
        ids = [e["id"]] if cond.cards == 1 else [f"{e['id']}a", f"{e['id']}b"]
        for i in ids:
            out[i] = fn(cond)
    return out


def _all_establishable():
    return dict(establishable={f"{c.name}:{v}": dict(power=1.0, establishable=c.definite is not None)
                               for c in P.CONDITIONS for v, _ in c.variants})


def test_score_labels():
    pytest.importorskip("scipy")
    import score as S
    assert S.label("SUPPORTED (provisional until replicated) [non-directional claim]") == "SUPPORTED"
    assert S.label("SUPPORTED — replicated") == "SUPPORTED"
    assert S.label("NOT SUPPORTED — metric invalid: x") == "NOT SUPPORTED"
    assert S.label("NO DETECTABLE EFFECT — x [parametric only]") == "NO DETECTABLE EFFECT"
    assert S.label("DEGENERATE METRIC — x") == "DEGENERATE METRIC"
    assert S.label(None) == S.label("garbage") == "ERROR"


def test_score_passes_a_perfect_engine_and_fails_the_two_failure_modes():
    pytest.importorskip("scipy")
    import score as S
    entries, pilot = P.assign(3), _all_establishable()
    perfect = S.score(entries, _verdicts(entries, lambda c: f"{c.definite or c.allowed[0]} — x"), pilot)
    assert perfect["passed"] and all(perfect["criteria"][s]["passed"] for s in ("S1", "S2", "S3", "S4"))
    shy = S.score(entries, _verdicts(entries, lambda c: "INCONCLUSIVE — x"), pilot)
    assert not shy["criteria"]["S3"]["passed"] and shy["criteria"]["S1"]["passed"]
    credulous = S.score(entries, _verdicts(entries, lambda c: "SUPPORTED (provisional until replicated)"), pilot)
    assert not credulous["criteria"]["S1"]["passed"] and not credulous["criteria"]["S2"]["passed"]
    assert not credulous["criteria"]["S4"]["passed"]
    missing = S.score(entries, {}, pilot)
    assert missing["criteria"]["S4"]["rate"] == 1.0


# --------------------------------------------------------------------------- #
# the runner, end to end on the real engine
# --------------------------------------------------------------------------- #
def test_runner_runs_the_engine_on_claim_cards(bgs, tmp_path):
    entries = P.assign(5)
    pick = [next(e["id"] for e in entries if e["condition"] == c) for c in ("N1", "E1", "N6b")]
    P.write_panel(5, bgs, dict(PILOT, key_dose=6.0), tmp_path / "panel", only=pick)
    import run_panel as R
    summary = R.run(tmp_path / "panel", tmp_path / "reports", workers=1)
    assert summary["run"] == 3 and summary["errors"] == 0 and not summary["logged_twice"]
    assert summary["machine"]["cpus"] >= 1
    for cid in pick:
        rep = json.loads((tmp_path / "reports" / f"{cid}.json").read_text())
        assert rep["verdict"] and rep["elapsed_seconds"] > 0
    assert R.run(tmp_path / "panel", tmp_path / "reports", workers=1)["skipped"] == 3  # one attempt per card
    log = (tmp_path / "reports" / "runlog.jsonl").read_text().splitlines()
    assert len(log) == 3


def test_runner_pins_one_blas_thread_per_worker_before_numpy_loads():
    """The thread variables count only if they are set before numpy loads its BLAS; a forked
    worker inherits the pool numpy was loaded with."""
    import os
    import subprocess
    env = {k: v for k, v in os.environ.items() if k not in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
    code = "import json, run_panel; print(json.dumps(run_panel.machine()))"
    out = subprocess.run([sys.executable, "-c", code], cwd=HERE, env=env, capture_output=True,
                         text=True, check=True).stdout
    m = json.loads(out)
    assert set(m["blas_env"].values()) == {"1"}
    if m["blas_threads_in_use"] is not None:  # read with threadpoolctl where it is installed
        assert m["blas_threads_in_use"] == [1]
