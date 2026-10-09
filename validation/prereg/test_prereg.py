"""Tests of the confirmatory panel's own code: panel.py, oracle.py, score.py, oc.py,
run_panel.py, blind.py and beacon.py.

They run on simulated backgrounds (simulate.py; no network). The panel, the oracle and the
scoring must not import the engine; scoring and the oracle need scipy, the beacon py_ecc
(skipped without them, as in the core-only CI job). This file itself may import the engine, to
check that the panel's copies of the engine's rules and constants agree with it.
"""
from __future__ import annotations

import ast
import hashlib
import json
import math
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel as P  # noqa: E402
from simulate import simulated_background  # noqa: E402

KEY = "00112233445566778899aabbccddeeff" * 2
OTHER_KEY = "ffeeddccbbaa99887766554433221100" * 2


def _pilot(bgs, sesoi=0.15, dose=3.0, delta=0.3, truth=None):
    """A pilot with every field the panel reads; `truth` maps a level to the metric's truth on
    its pairs (default: valid everywhere)."""
    truth = truth or {}
    out = dict(pool={}, pool_size={}, sesoi={}, saturation_dose={}, truth={}, e_dose={lv: dose for lv in P.LEVELS},
               establishable={}, dropped=[])
    for name, bg in bgs.items():
        pool = bg.plan["pool"]
        out["pool"][name] = [dict(index=pe["index"], level=pe["level"], pair=pe["pair"]) for pe in pool]
        out["pool_size"][name] = len(pool)
        out["sesoi"][name] = {lv: sesoi for lv in P.LEVELS}
        out["saturation_dose"][name] = {lv: dose for lv in P.LEVELS}
        out["truth"][name] = {str(k): {"class": truth.get(pe["level"], "valid")} for k, pe in enumerate(pool)}
    out["delta"] = {str(k): {**{f"{f:g}": dict(value=delta * f, se=0.01) for f in (0.25, 0.5, 1.0, 1.5)},
                             f"capture={P.E_CAPTURE:g}": dict(value=0.6 * delta, se=0.01)}  # E2, E3 at their depth
                    for k in range(len(bgs["B1"].plan["pool"]))}
    return out


@pytest.fixture(scope="module")
def bgs():
    return {"B1": simulated_background("B1"), "B2": simulated_background("B2", donors=12, seed=1)}


def _entry(condition, variant, side="A", seed=0, i=0, pair=0):
    return dict(id=f"T{i}", condition=condition, variant=variant, index=i, side=side, seed=seed, pair=pair)


def _sizes(bgs):
    return {k: len(b.plan["pool"]) for k, b in bgs.items()}


# --------------------------------------------------------------------------- #
# the key and the assignment
# --------------------------------------------------------------------------- #
def test_the_key_is_a_drand_randomness_of_64_hex():
    for bad in (1, "0" * 32, "0" * 63, "0" * 65, "G" * 64, "AB" * 32, " " + "0" * 64):
        with pytest.raises(ValueError):
            P.check_key(bad)
    assert P.check_key(KEY) == KEY
    assert not hasattr(P, "key_commitment")  # no secret key, no commitment (decided 2026-10-08)


def test_assign_is_deterministic_complete_and_keyed():
    a, b, c = P.assign(KEY), P.assign(KEY), P.assign(OTHER_KEY)
    assert a == b and a != c
    assert len(a) == P.n_datasets() == 6400
    assert P.n_cards() == 6700 and P.n_cards(data="null") == 5800  # N4 has two cards
    assert len({e["id"] for e in a}) == len(a) and {e["id"] for e in a} == {e["id"] for e in c}
    counts = Counter((e["condition"], e["variant"]) for e in a)
    assert counts == {(cond.name, v): n for cond in P.CONDITIONS for v, n in cond.variants}
    assert sum(x["condition"] != y["condition"] for x, y in zip(a, c)) > 0.5 * len(a)
    assert {e["side"] for e in a} == {"A", "B"}
    assert {e["pair"] for e in a} == set(range(len(P.LEVELS) * P.MAX_PAIRS_PER_LEVEL))


def test_the_pair_is_the_pair_draw_modulo_the_pool_size():
    a = P.assign(KEY, pool_sizes={"B1": 15, "B2": 12})
    for e in a:
        size = 12 if P.conditions()[e["condition"]].background == "B2" else 15
        assert e["pair"] == e["pair_draw"] % size
    assert {e["pair"] for e in a if e["condition"] == "N7"} == set(range(12))
    # the pool size changes the pair, nothing else
    full = P.assign(KEY)
    assert [(x["id"], x["condition"], x["variant"], x["side"], x["seed"]) for x in a] == \
           [(x["id"], x["condition"], x["variant"], x["side"], x["seed"]) for x in full]


def test_key_null_conditions_have_790_datasets_and_n8_is_one_of_them():
    keyed = {c.name for c in P.CONDITIONS if c.key}
    assert keyed == {"N1", "N2", "N5", "N6c", "N8"}
    for cond in P.CONDITIONS:
        if cond.key:
            assert dict(cond.variants)[cond.key_variant] == P.KEY_N == 790


def test_the_gene_pair_is_drawn_independently_of_the_condition():
    """The pair (and so its level and its truth) does not depend on the condition: the pair counts
    of N1 and of the real effects are both uniform over the pool."""
    stats = pytest.importorskip("scipy.stats")
    a = P.assign(KEY)
    k = len(P.LEVELS) * P.MAX_PAIRS_PER_LEVEL
    for group in ({"N1"}, {"E1", "E2", "E3"}, {"N8"}):
        counts = np.bincount([e["pair"] for e in a if e["condition"] in group], minlength=k)
        assert stats.chisquare(counts).pvalue > 1e-3
    table = np.array([np.bincount([e["pair"] for e in a if e["condition"] == c], minlength=k) for c in ("N1", "N8")])
    assert stats.chi2_contingency(table)[1] > 1e-3


def test_drops_follow_the_pre_registered_order_and_never_touch_key_conditions():
    full = P.n_datasets()
    assert P.n_datasets(["N3:f=0.1", "N3:f=0.2", "N3:f=0.4"]) == full - 300
    assert len(P.assign(KEY, ["N3:f=0.1"])) == full - 100
    for bad in (["N3:f=0.2"], ["N1:null"], ["E1:dose=key"], ["N3:f=0.1", "N2:c=0.5"]):
        with pytest.raises(ValueError):
            P.assign(KEY, bad)
    keyed = {f"{c.name}:{c.key_variant}" for c in P.CONDITIONS if c.key}
    assert not keyed & {f"{n}:{v}" for n, v in P.DROP_ORDER}
    assert not {f"E{i}" for i in (1, 2, 3)} & {f"{n}:{v}" for n, v in P.DROP_ORDER if v in ("dose=key",)}


# --------------------------------------------------------------------------- #
# backgrounds: the pair pool by level and its controls
# --------------------------------------------------------------------------- #
def test_the_plan_finds_eight_disjoint_pairs_per_level_with_their_controls(bgs):
    plan = bgs["B1"].plan
    pool = plan["pool"]
    k = plan["pairs_per_level"]
    assert k == P.MAX_PAIRS_PER_LEVEL == 8 and not plan["fewer_pairs_because"]
    assert [pe["level"] for pe in pool] == [lv for lv in P.LEVELS for _ in range(k)]
    assert [pe["index"] for pe in pool] == list(range(len(pool)))
    genes = [g for pe in pool for g in (*pe["pair"], *pe["neg_pair"])] + plan["positive_control"]["pair"]
    assert len(genes) == len(set(genes))  # pairs, negative controls and the positive control disjoint
    for pe in pool:
        assert abs(pe["neg_r"]) < P.NEG_MAX_R
        assert pe["pos_pair"] == plan["positive_control"]["pair"]
        for m, d in zip(pe["mean"], pe["detection"]):
            assert P.level_of(m, d) == pe["level"]
        for m in pe["neg_mean"]:
            assert m == pytest.approx(np.mean(pe["mean"]), rel=1.0)
    by_level = {lv: [pe for pe in pool if pe["level"] == lv] for lv in P.LEVELS}
    assert min(pe["r"] for pe in by_level["high"]) > 0.1  # the planted couplings are found
    assert set(plan["genes"]) >= set(genes) | set(plan["g2m"]) | set(plan["random_genes"])
    assert bgs["B1"].X.dtype == np.int32 and bgs["B1"].X.shape[1] == len(plan["genes"])
    assert P.level_of_pair(0, len(pool)) == "high" and P.level_of_pair(len(pool) - 1, len(pool)) == "low"


def test_the_pool_takes_the_largest_common_number_of_pairs_down_to_four():
    """Where a level cannot give 8 pairs with their negative controls, every level gets the
    largest number it can, down to 4; below 4 the background does not qualify."""
    small = simulated_background("B1", n_genes=96, coupled=6)  # 32 genes per level: room for 6 pairs
    k = small.plan["pairs_per_level"]
    assert 4 <= k < 8 and small.plan["fewer_pairs_because"]
    assert [pe["level"] for pe in small.plan["pool"]] == [lv for lv in P.LEVELS for _ in range(k)]
    with pytest.raises(ValueError, match="no pool"):
        simulated_background("B1", n_genes=48, coupled=3)


def test_level_rule():
    assert P.level_of(5.0, 0.95) == "high" and P.level_of(1.2, 0.6) == "medium"
    assert P.level_of(0.3, 0.25) == "low" and P.level_of(3.0, 0.4) == "low"
    assert P.level_of(0.05, 0.05) is None


def test_a_compact_background_round_trips(bgs, tmp_path):
    P.save_compact(bgs["B1"], tmp_path / "B1.npz")
    back = P.load_compact(tmp_path / "B1.npz")
    assert P.background_sha256(back) == P.background_sha256(bgs["B1"])
    assert back.plan == bgs["B1"].plan


def test_backgrounds_load_from_their_spec_filtered_sparse_and_planned(tmp_path):
    """The real-data path of panel.py, oracle.py, timing.py and blind.py (`load_backgrounds`):
    a background named in a backgrounds JSON (dense or sparse .npz) loads as raw counts, keeps
    only the filtered cells, and gets the plan the same counts get in memory; counts that are not
    integers are refused."""
    sparse = pytest.importorskip("scipy.sparse")
    bg = simulated_background("B1", donors=6, cells=240, plan=False)
    cell_type = np.where(np.arange(len(bg.donor)) % 12 == 0, "other", "fibroblast")
    obs = pd.DataFrame({"donor_id": bg.donor, "cell_type": cell_type})
    keep = cell_type == "fibroblast"
    for name, X in (("dense", bg.X), ("sparse", sparse.csr_matrix(bg.X))):
        P.save_npz(tmp_path / f"{name}.npz", X, obs, bg.genes)
        spec = {"B1": {"path": f"{name}.npz", "donor": "donor_id", "filter": {"cell_type": "fibroblast"}}}
        (tmp_path / "backgrounds.json").write_text(json.dumps(spec))
        got = P.load_backgrounds(tmp_path / "backgrounds.json")["B1"]
        ref = P.Background(bg.X[keep], bg.genes, bg.donor[keep], "B1")
        P.plan_background(ref)
        assert got.plan == ref.plan and np.array_equal(got.X, ref.X) and list(got.donor) == list(ref.donor)
    # the form backgrounds.json takes: a file name in the data directory, checked against its sha256
    data = tmp_path / "data"
    data.mkdir()
    P.save_npz(data / "b1.npz", bg.X, obs, bg.genes)
    spec = {"B1": {"file": "b1.npz", "sha256": P.sha256(data / "b1.npz"), "donor": "donor_id",
                   "filter": {"cell_type": "fibroblast"}, "url": "https://example.org/b1"}}
    (tmp_path / "backgrounds.json").write_text(json.dumps(spec))
    assert P.load_backgrounds(tmp_path / "backgrounds.json", data)["B1"].plan == ref.plan
    spec["B1"]["sha256"] = "0" * 64
    (tmp_path / "backgrounds.json").write_text(json.dumps(spec))
    with pytest.raises(ValueError, match="sha256"):
        P.load_backgrounds(tmp_path / "backgrounds.json", data)
    np.savez_compressed(tmp_path / "dense.npz", X=bg.X + 0.5, genes=np.asarray(bg.genes),
                        obs_donor_id=bg.donor, obs_cell_type=cell_type)
    spec = {"B1": {"path": "dense.npz", "donor": "donor_id"}}
    (tmp_path / "backgrounds.json").write_text(json.dumps(spec))
    with pytest.raises(ValueError, match="raw counts"):
        P.load_backgrounds(tmp_path / "backgrounds.json")


def test_an_h5ad_background_uses_raw_counts_and_gene_symbols(tmp_path):
    """CELLxGENE files hold normalized values in X, raw counts in raw.X, Ensembl IDs as
    var_names and symbols in var['feature_name']; the cells are filtered before loading."""
    anndata = pytest.importorskip("anndata")
    sparse = pytest.importorskip("scipy.sparse")
    bg = simulated_background("B1", donors=6, cells=240, plan=False)
    ens = [f"ENSG{i:011d}" for i in range(len(bg.genes))]
    var = pd.DataFrame({"feature_name": bg.genes}, index=ens)
    obs = pd.DataFrame({"donor_id": bg.donor, "cell_type": "fibroblast"}, index=[f"c{i}" for i in range(len(bg.donor))])
    raw = anndata.AnnData(sparse.csr_matrix(bg.X), obs=obs, var=var)
    ad = anndata.AnnData(sparse.csr_matrix(np.log1p(bg.X)), obs=obs, var=var)
    ad.raw = raw
    ad.write_h5ad(tmp_path / "b1.h5ad")
    got = P.load_background("B1", {"path": str(tmp_path / "b1.h5ad"), "donor": "donor_id",
                                   "filter": {"cell_type": "fibroblast"}})
    assert got.genes == bg.genes
    assert np.array_equal(got.X.toarray(), bg.X)


# --------------------------------------------------------------------------- #
# datasets and claim cards
# --------------------------------------------------------------------------- #
def test_claim_cards_do_not_reveal_the_condition(bgs):
    """For a given pair, cards of the norm_pearson conditions on B1 — every condition and variant
    (the fourth review: N5 and E1's other doses were not checked) — differ only in the id and the
    claimed direction, and N4's first card in its missing replicate unit (another design, by
    construction); the pair, its level's SESOI and delta_min, the injection's strength and the
    controls come from the pair."""
    pilot = _pilot(bgs)
    cases = [(c.name, v) for c in P.CONDITIONS if c.background == "B1" and c.metric == "norm_pearson"
             for v, _ in c.variants]
    assert {"N5", "N4"} <= {c for c, _ in cases} and len([c for c in cases if c[0] == "E1"]) == 4
    for pair in (0, 9, 20):
        cards = []
        for c, v in cases:
            built = P.build(_entry(c, v, seed=3, pair=pair), bgs, pilot)[3]
            if c == "N4":
                assert {k: x for k, x in built[0].items() if k not in ("id", "replicate_col")} == \
                    {k: x for k, x in built[1].items() if k not in ("id", "replicate_col")}
                assert built[0]["replicate_col"] is None and built[1]["replicate_col"] == "donor"
            cards.append(built[-1])
        strip = [{k: v for k, v in c.items() if k not in ("id", "prereg")} for c in cards]
        assert all(s == strip[0] for s in strip)
        assert all({k: v for k, v in c["prereg"].items() if k != "direction"} ==
                   {k: v for k, v in cards[0]["prereg"].items() if k != "direction"} for c in cards)
        assert cards[0]["gene_pair"] == bgs["B1"].plan["pool"][pair]["pair"]
        assert cards[0]["signal_test"]["strength"] == 3.0
        assert cards[0]["prereg"]["delta_min"] == pytest.approx(P.DELTA_MIN_FRACTION * 0.15)
    dirs = Counter(P.build(_entry("N1", "null", seed=s, i=s), bgs, pilot)[3][0]["prereg"]["direction"]
                   for s in range(30))
    assert set(dirs) == {"increase", "decrease"}
    assert P.build(_entry("E1", "dose=key", side="A"), bgs, pilot)[3][0]["prereg"]["direction"] == "decrease"
    assert P.build(_entry("E1", "dose=key", side="B"), bgs, pilot)[3][0]["prereg"]["direction"] == "increase"
    # a negative Δ* makes the signal side the lower one
    neg = _pilot(bgs, delta=-0.3)
    assert P.build(_entry("E1", "dose=key", side="A"), bgs, neg)[3][0]["prereg"]["direction"] == "increase"
    # B2's cards carry B2's own SESOI and saturation dose
    p2 = _pilot(bgs)
    p2["sesoi"]["B2"] = {lv: 0.3 for lv in P.LEVELS}
    p2["saturation_dose"]["B2"] = {lv: 1.5 for lv in P.LEVELS}
    card = P.build(_entry("N7", "mice", pair=1), bgs, p2)[3][0]
    assert card["prereg"]["sesoi"] == 0.3 and card["signal_test"]["strength"] == 1.5


def test_designs_and_artifacts_are_planted_as_specified(bgs):
    pilot = _pilot(bgs)
    X, obs, genes, _ = P.build(_entry("N1", "null", seed=1), bgs, pilot)
    per_group = obs.groupby("group")["donor"].nunique()
    assert per_group["A"] == per_group["B"] == P.DONORS_PER_GROUP
    assert X.shape[0] == 2 * P.DONORS_PER_GROUP * P.CELLS_PER_DONOR
    X, obs, _, _ = P.build(_entry("N2", "c=0.5", side="B", seed=1), bgs, pilot)
    tot = obs.groupby("group")["total_counts"].mean()
    assert tot["B"] / tot["A"] == pytest.approx(0.5, abs=0.08)
    X, obs, _, _ = P.build(_entry("N3", "f=0.4", side="A", seed=1), bgs, pilot)
    det = (X > 0).mean(axis=1)
    g = np.asarray(obs["group"])
    assert det[g == "A"].mean() / det[g == "B"].mean() == pytest.approx(0.6, abs=0.08)
    # N8: per-cell capture ~ Beta(2, 2): the same mean loss as c = 0.5, a wider spread of depth
    X8, obs8, _, _ = P.build(_entry("N8", "beta(2,2)", side="B", seed=1), bgs, pilot)
    X2, obs2, _, _ = P.build(_entry("N2", "c=0.5", side="B", seed=1), bgs, pilot)
    t8 = obs8.groupby("group")["total_counts"]
    assert t8.mean()["B"] / t8.mean()["A"] == pytest.approx(0.5, abs=0.08)
    cv = lambda s: float(np.std(np.log(s + 1)))  # noqa: E731
    assert cv(obs8.loc[obs8.group == "B", "total_counts"]) > cv(obs2.loc[obs2.group == "B", "total_counts"]) + 0.1
    X, obs, _, cards = P.build(_entry("N4", "3v3", seed=1), bgs, pilot)
    assert obs.groupby("group")["donor"].nunique().tolist() == [3, 3]
    assert [c["replicate_col"] for c in cards] == [None, "donor"]
    X, obs, _, _ = P.build(_entry("N5", "sham", seed=1), bgs, pilot)
    assert (obs.groupby("donor")["group"].nunique() == 2).all()
    X, obs, _, _ = P.build(_entry("N7", "mice", seed=1), bgs, pilot)
    assert obs["donor"].str.startswith("B2").all()
    assert obs["donor"].nunique() == len(bgs["B2"].donors) == 12  # fewer than N7_MAX_MICE: every mouse
    with pytest.raises(ValueError, match="outside"):
        P.build(_entry("N1", "null", pair=99), bgs, pilot)
    # the third review: a real effect's direction is the sign of its Δ*, never a default
    no_delta = dict(pilot, delta={k: v for k, v in pilot["delta"].items() if k != "1"})
    with pytest.raises(ValueError, match="Δ\\* of pair 1"):
        P.build(_entry("E1", "dose=key", pair=1), bgs, no_delta)
    P.build(_entry("E1", "dose=key", pair=0), bgs, no_delta)  # another pair's is known


def test_n7_takes_at_most_its_cap_of_mice_so_its_memory_is_bounded(bgs):
    """The third review's B4: N7 took every qualifying mouse of B2, whose number has no bound, and
    a worker killed for memory would leave no results after the key. At most N7_MAX_MICE mice per
    dataset, drawn per dataset, split into two halves (with an odd number B has one more)."""
    big = {"B1": bgs["B1"], "B2": simulated_background("B2", donors=26, cells=210, n_genes=420, seed=3)}
    pilot = _pilot(big)
    seen = set()
    for seed in range(3):
        X, obs, _, _ = P.build(_entry("N7", "mice", seed=seed), big, pilot)
        per = obs.groupby("group")["donor"].nunique()
        assert obs["donor"].nunique() == P.N7_MAX_MICE and per.to_dict() == {"A": 12, "B": 12}
        assert X.shape[0] == P.N7_MAX_MICE * P.CELLS_PER_DONOR
        seen |= set(obs["donor"])
    assert len(seen) > P.N7_MAX_MICE  # the mice are drawn per dataset
    odd = {"B1": bgs["B1"], "B2": simulated_background("B2", donors=13, cells=210, n_genes=420, seed=4)}
    _, obs, _, _ = P.build(_entry("N7", "mice"), odd, _pilot(odd))
    assert obs.groupby("group")["donor"].nunique().to_dict() == {"A": 6, "B": 7}


def test_the_memory_rule_stops_the_pilot_before_the_key_and_a_killed_worker_fails_fast():
    """The third review's B4: the drops weighed time only, and multiprocessing.Pool waits for a
    worker the system killed until the job's timeout. The timing pilot checks the run's workers at
    a worker's measured peak against the runner's memory and stops before the key; the pools break
    at once when a worker dies."""
    import time
    import timing as T
    import run_panel as R
    assert "4 workers = 2000 MB" in T.check_memory(500.0, 16000.0)
    with pytest.raises(SystemExit, match="stops before the key"):
        T.check_memory(3500.0, 16000.0)
    with pytest.raises(SystemExit, match="unknown"):
        T.check_memory(500.0, None)
    assert "14000 MB" in T.check_memory(3500.0, 16000.0, enforce=False)  # timing.log reports it, the pilot enforces it
    assert "unknown" in T.check_memory(500.0, None, enforce=False)
    assert T.runner_memory_mb() > 0
    from concurrent.futures.process import BrokenProcessPool
    t0 = time.time()
    with pytest.raises(BrokenProcessPool):
        R.fork_map(_die_on_three, list(range(8)), 2)
    assert time.time() - t0 < 60
    assert R.fork_map(abs, [-1, -2, 3], 2) == [1, 2, 3]


def _die_on_three(x):
    if x == 3:
        os._exit(9)  # as the kernel's out-of-memory killer ends a process: no exception, no cleanup
    return x


def test_injected_coupling_raises_the_pair_correlation_on_its_side_only(bgs):
    import oracle as O
    bg = bgs["B1"]
    pe = bg.plan["pool"][0]
    ia, ib = bg.genes.index(pe["pair"][0]), bg.genes.index(pe["pair"][1])
    X, obs, _, _ = P.build(_entry("E1", "dose=key", side="A", seed=2, pair=0), bgs, _pilot(bgs, dose=3.0))
    g = np.asarray(obs["group"])
    ra, rb = O.norm_pearson(X[g == "A"], ia, ib), O.norm_pearson(X[g == "B"], ia, ib)
    assert ra > rb + 0.1
    ma, mb = X[g == "A"][:, [ia, ib]].mean(), X[g == "B"][:, [ia, ib]].mean()
    assert ma / mb == pytest.approx(1.0, abs=0.15)  # the sham thins the other side's pair alike


def test_dataset_and_card_hashes_are_canonical(bgs):
    pilot = _pilot(bgs)
    e = _entry("N2", "c=0.5", seed=7, pair=4)
    X, obs, genes, cards = P.build(e, bgs, pilot)
    X2, obs2, genes2, cards2 = P.build(e, bgs, pilot)
    assert P.dataset_sha256(X, obs, genes) == P.dataset_sha256(X2, obs2, genes2)
    assert P.card_sha256(cards[0]) == P.card_sha256(json.loads(json.dumps(cards2[0])))
    X2[0, 0] += 1
    assert P.dataset_sha256(X, obs, genes) != P.dataset_sha256(X2, obs, genes)


def test_the_panel_the_oracle_the_scoring_and_the_beacon_never_import_the_engine():
    for name in ("panel.py", "oracle.py", "score.py", "oc.py", "simulate.py", "beacon.py", "frozen.py"):
        tree = ast.parse((HERE / name).read_text())
        mods = {a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
        mods |= {n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module}
        assert not any(m.split(".")[0] == "metric_autopsy" for m in mods), (name, mods)


def test_the_panels_copies_of_the_engines_rules_agree_with_it():
    """The panel never imports the engine, so the rules it shares with it are copies: delta_min's
    default, GATE 4's number of injections, the module injection, every cause the engine can
    give, and GATE 4's FAIL rule."""
    pytest.importorskip("scipy")
    import oracle as O
    import run_panel as R
    from metric_autopsy import gates as G
    from metric_autopsy.report import CAUSES
    assert P.DELTA_MIN_FRACTION == G.DELTA_MIN_FRACTION and O.GATE4_REPS == G.GATE4_N_REP
    assert (O.MODULE_FOLD, O.MODULE_FRAC) == (R.MODULE_FOLD, R.MODULE_FRAC)
    assert set(CAUSES) == set(P.CAUSE_OUTCOME)
    rng = np.random.default_rng(0)
    for _ in range(300):
        d = rng.normal(rng.uniform(-0.05, 0.1), rng.uniform(0.001, 0.2), size=int(rng.integers(5, 60)))
        dmin = float(rng.uniform(0.01, 0.1))
        assert O.shows_blind(d, dmin) == (G.judge_response(G.response_interval(d), dmin) == "FAIL")
    # the engine's verdict order where GATE 4 ran (oracle.engine_outcome against report.decide_cause):
    # "explained by depth" comes before the metric's validity, which the second review found unmodelled
    from metric_autopsy.core import Assessment
    from metric_autopsy.report import decide_cause
    from types import SimpleNamespace
    effects = {P.NS_DEPTH: Assessment("INCONCLUSIVE", "x", [], dict(explained_by_depth=True)),
               P.NDE: Assessment("NO_DETECTABLE_EFFECT", "x", [], {}),
               P.SUPPORTED: Assessment("DETECTED", "x", [], dict(effect=-1.0)),
               P.INCONCLUSIVE: Assessment("INCONCLUSIVE", "x", [], {})}
    validity = {"FAIL": Assessment("FAIL", "x", [], dict(failed_gate=4)), "PASS": Assessment("PASS", "x"),
                "UNTESTED": Assessment("UNTESTED", "x")}
    for g4, mv in validity.items():
        for outcome, ef in effects.items():
            a = SimpleNamespace(metric_validity=mv, design_adequacy=Assessment("ADEQUATE", "x"), effect=ef,
                                replication=Assessment("NOT_RUN", "x"), params=dict(groups=["A", "B"]),
                                prereg=dict(direction="increase"), judgment_pending=False)
            label, cause = decide_cause(a)
            assert P.outcome(label, cause) == O.engine_outcome(outcome, g4), (g4, outcome, cause)
    # the third review: the negative control is chosen against GATE 5's matched null, and the
    # oracle's analysis equalizes depth as the engine's effect estimate does
    import inspect
    from metric_autopsy import effect as EF
    assert P.MATCHED_NEIGHBOURS_MIN == G.MATCHED_NEIGHBOURS_MIN
    assert P.MATCHED_NEIGHBOURS_FRAC == inspect.signature(G._neighbourhood).parameters["frac"].default
    assert P.NEG_NULL_PAIRS == inspect.signature(G.gate5_controls).parameters["n_null"].default
    assert O.EQUALIZE_DRAWS == inspect.signature(EF.estimate_effect).parameters["n_equalize"].default
    # the fourth review: the oracle's Monte Carlo permutation test with the engine's numbers (run_panel
    # sets neither, so run_autopsy's n_perm and estimate_effect's max_exact apply)
    from metric_autopsy.report import run_autopsy
    assert O.N_PERM == inspect.signature(run_autopsy).parameters["n_perm"].default
    assert O.MAX_EXACT == inspect.signature(EF.estimate_effect).parameters["max_exact"].default
    src = inspect.getsource(R.make_run_args)
    assert "n_perm" not in src and "max_exact" not in src
    from metric_autopsy import SimpleData, metrics  # and the metric the controls are chosen by
    X = rng.poisson(rng.gamma(0.5, 2.0, size=(300, 6))).astype(float)
    data = SimpleData(X, pd.DataFrame(index=range(300)), [f"g{i}" for i in range(6)])
    for a, b in ((0, 1), (2, 3), (4, 5)):
        assert P.norm_pearson_cols(X[:, a], X[:, b], X.sum(axis=1)) == pytest.approx(
            metrics.norm_pearson(data, gene_a=f"g{a}", gene_b=f"g{b}"), abs=1e-12)


def test_the_negative_control_is_the_candidate_most_typical_of_gate5s_null(monkeypatch):
    """The third review's B3: controls chosen by partial correlation alone stood out of GATE 5's
    matched null on some pairs. Of the NEG_CANDIDATES qualifying pairs closest in mean the plan
    takes the one with the smallest |z| against that null (ties: the closest), and records it."""
    bg = simulated_background("B1", plan=False)
    names = list(bg.genes)
    calls = []

    def score(value, means, kept, pair, exclude, rng):  # a fixed pseudo-score per candidate
        c, d = pair
        calls.append((frozenset(exclude) - {c, d}, (c, d), ((c * 7919 + d * 104729) % 997) / 100))
        return calls[-1][2]
    monkeypatch.setattr(P, "_gate5_typicality", score)
    plan = P.plan_background(bg)
    assert not plan["fewer_pairs_because"]  # one pass, so the calls are this pool's
    for pe in plan["pool"]:
        a, b = names.index(pe["pair"][0]), names.index(pe["pair"][1])
        cand = [(z, cd) for owner, cd, z in calls if {a, b} <= owner]
        assert 1 < len(cand) <= P.NEG_CANDIDATES
        z, (c, d) = min(cand, key=lambda t: t[0])
        assert pe["neg_pair"] == [names[c], names[d]] and pe["neg_gate5_z"] == z


def test_gate5_on_the_panels_null_datasets_fails_no_more_than_a_sound_validator(bgs):
    """The sound validator's model (oc.sound_model) lets GATE 5 FAIL a null dataset with 1.5 alpha;
    the third review found the panel's own negative controls making it FAIL far more often on
    some pairs. The engine's GATE 5, as run_autopsy calls it for a card, on one N1 and one N7
    dataset per pool pair: its failures are within what 1.5 alpha gives at the 1% level."""
    import oc
    from metric_autopsy import SimpleData, gates as G, metrics
    pilot = _pilot(bgs)
    fails = runs = 0
    for name in ("N1", "N7"):
        cond = P.conditions()[name]
        for k in range(len(bgs[cond.background].plan["pool"])):
            X, obs, genes, cards = P.build(_entry(name, cond.variants[0][0], seed=1000 * k, pair=k), bgs, pilot)
            card = cards[0]
            g5 = G.gate5_controls(metrics.norm_pearson, SimpleData(X, obs, genes), tuple(card["pos_pair"]),
                                  tuple(card["neg_pair"]), exclude=card["gene_pair"], seed=k,
                                  delta_min=card["prereg"]["delta_min"])
            fails += g5.status.value == "FAIL"
            runs += 1
    import math  # the binomial tail without scipy, which the core-only CI jobs lack
    p = 1.5 * oc.ALPHA
    below = sum(math.comb(runs, j) * p ** j * (1 - p) ** (runs - j) for j in range(fails))  # P(X <= fails - 1)
    assert runs >= 24 and below < 0.99, (fails, runs)


def test_the_oracles_gate4_response_matches_the_engines_injection(bgs):
    """The oracle re-implements GATE 4's measurement (injected minus sham on the whole dataset)
    without copying the matrix; on the same dataset its mean response agrees with the engine's."""
    pytest.importorskip("scipy")
    from functools import partial

    import oracle as O
    from metric_autopsy import SimpleData, gates as G, injected_signal, metrics
    pilot = _pilot(bgs)
    pe = bgs["B1"].plan["pool"][1]
    X, obs, genes, _ = P.build(_entry("N1", "null", seed=4, pair=1), bgs, pilot)
    ia, ib = genes.index(pe["pair"][0]), genes.index(pe["pair"][1])
    ours = O.coupling_deltas(X, ia, ib, 2.0, np.random.default_rng(1), 200)
    theirs = np.asarray(G.injected_deltas(partial(metrics.norm_pearson, gene_a=pe["pair"][0], gene_b=pe["pair"][1]),
                                          SimpleData(X, obs, genes), injected_signal.coupling(*pe["pair"], strength=2.0),
                                          200, np.random.default_rng(2)))
    se = np.sqrt(ours.var(ddof=1) / len(ours) + theirs.var(ddof=1) / len(theirs))
    assert abs(ours.mean() - theirs.mean()) < 4 * se and ours.mean() > 0.05


def test_the_oracles_module_response_matches_the_engines_injection(bgs):
    """The first review's C16: the module injection is checked against the engine itself, not only
    between the panel's copies: the fold and share of run_panel and the oracle are the engine's
    defaults, and on the same N6c dataset the oracle's response of the random-gene score to the
    G2M module (module_deltas) agrees with the engine's (injected_signal.module)."""
    pytest.importorskip("scipy")
    import inspect

    import oracle as O
    import run_panel as R
    from metric_autopsy import SimpleData, gates as G, injected_signal
    sig = inspect.signature(injected_signal.module).parameters
    assert (O.MODULE_FOLD, O.MODULE_FRAC) == (R.MODULE_FOLD, R.MODULE_FRAC) == (sig["fold"].default, sig["frac"].default)
    pilot = _pilot(bgs)
    X, obs, genes, _ = P.build(_entry("N6c", "random-genes", seed=4, pair=1), bgs, pilot)
    plan = bgs["B1"].plan
    score_cols = [genes.index(g) for g in plan["random_genes"]]
    module_cols = [genes.index(g) for g in plan["g2m"]]
    ours = O.module_deltas(np.asarray(X, float), score_cols, module_cols, np.random.default_rng(1), 200)
    theirs = np.asarray(G.injected_deltas(lambda d: O.module_score(np.asarray(d.X, float), score_cols),
                                          SimpleData(np.asarray(X, float), obs, genes),
                                          injected_signal.module(plan["g2m"], fold=R.MODULE_FOLD, frac=R.MODULE_FRAC),
                                          200, np.random.default_rng(2)))
    se = np.sqrt(ours.var(ddof=1) / len(ours) + theirs.var(ddof=1) / len(theirs))
    assert abs(ours.mean() - theirs.mean()) < 4 * se and abs(theirs.mean()) > 4 * theirs.std(ddof=1) / np.sqrt(200)


# --------------------------------------------------------------------------- #
# outcomes, the truth and the allowed (label, cause) pairs
# --------------------------------------------------------------------------- #
def test_outcome_reads_the_label_and_the_cause():
    assert P.outcome("NOT SUPPORTED — metric invalid: injected signal: blind", "metric_invalid_gate4") == P.NS_INVALID
    assert P.outcome("NOT SUPPORTED — metric invalid: controls: x", "metric_invalid_gate5") == P.NS_INVALID
    assert P.outcome("NOT SUPPORTED — metric invalid: 'extra_dropout' biases", "metric_invalid_gate0") == P.REFUSAL
    assert P.outcome("NOT SUPPORTED — the raw difference ...", "explained_by_depth") == P.NS_DEPTH
    assert P.outcome("NOT SUPPORTED — the effect is in the direction opposite", "opposite_direction") == P.NS_OPPOSITE
    assert P.outcome("SUPPORTED (provisional until replicated) [non-directional claim]", "provisional") == P.SUPPORTED
    assert P.outcome("NO DETECTABLE EFFECT — x [underpowered relative to the SESOI]", "no_detectable_effect") == P.NDE
    assert P.outcome("INCONCLUSIVE — x", "detected_untested_metric") == P.INCONCLUSIVE
    assert P.outcome("DEGENERATE METRIC — x", "degenerate_metric") == P.DEGENERATE
    # a label that does not match its cause, an unknown cause, no report
    assert P.outcome("SUPPORTED — replicated", "no_detectable_effect") == P.ERROR
    assert P.outcome("INCONCLUSIVE — x", "something_new") == P.ERROR
    assert P.outcome(None, None) == P.ERROR and P.outcome("garbage", "provisional") == P.ERROR
    assert P.outcome("NOT SUPPORTED — the effect did not replicate", "not_replicated") == P.OTHER


def test_the_truth_about_the_metric_has_a_band_around_delta_min():
    assert P.classify_response(0.06, 0.05) == "valid" and P.classify_response(0.0599, 0.05) == "ambiguous"
    assert P.classify_response(0.04, 0.05) == "blind" and P.classify_response(0.0401, 0.05) == "ambiguous"
    assert P.classify_response(-0.3, 0.05) == "blind"


def test_the_allowed_outcomes_follow_the_owners_table_in_every_condition(bgs):
    """Decided 2026-10-08, the same for N1-N8 and E1-E3: a blind or useless metric allows metric
    invalid and INCONCLUSIVE (the constant also DEGENERATE METRIC); a valid metric on null data NO
    DETECTABLE EFFECT, INCONCLUSIVE, NOT SUPPORTED for the opposite direction, and explained by
    depth only where an artifact is planted; on a real effect SUPPORTED and INCONCLUSIVE, and NO
    DETECTABLE EFFECT only below the SESOI; an ambiguous metric the union; a refusal by GATE 0
    everywhere; metric invalid on a valid metric is never allowed."""
    conds = P.conditions()
    pilot = _pilot(bgs, sesoi=0.15, delta=0.3, truth={"high": "valid", "medium": "ambiguous", "low": "blind"})
    hi, med, lo = 0, 9, 20
    blind = {P.NS_INVALID, P.INCONCLUSIVE, P.REFUSAL}
    null = {P.NDE, P.INCONCLUSIVE, P.NS_OPPOSITE, P.REFUSAL}
    for c in P.CONDITIONS:
        for v, _ in c.variants:
            pair = hi
            got = P.allowed(c, v, pair, pilot)
            assert P.REFUSAL in got and P.REFUSAL not in P.definite(c, v, pair, pilot)
            if c.metric == "constant":
                assert got == blind | {P.DEGENERATE}
            elif c.metric != "norm_pearson":
                assert got == blind
            elif c.data == "null":
                assert got == (null | {P.NS_DEPTH} if c.artifact else null)
                assert P.NS_INVALID not in got and P.SUPPORTED not in got
            if c.metric == "norm_pearson":
                assert P.allowed(c, v, lo, pilot) == blind
                assert P.allowed(c, v, med, pilot) == P.data_allowed(c, v, med, pilot) | blind
    # real effects: NO DETECTABLE EFFECT only where |Δ*| < SESOI (0.3 x factor against 0.15)
    for name, variant in (("E1", "dose=key"), ("E2", "against"), ("E3", "with")):
        assert P.allowed(conds[name], variant, hi, pilot) == {P.SUPPORTED, P.INCONCLUSIVE, P.REFUSAL}
    assert P.allowed(conds["E1"], "dose=0.25", hi, pilot) == {P.SUPPORTED, P.INCONCLUSIVE, P.REFUSAL, P.NDE}
    assert P.NS_DEPTH not in P.allowed(conds["E3"], "with", hi, pilot)  # an artifact with a real effect
    # below the SESOI, where an artifact is planted, "explained by depth" is correct as NDE is: with a
    # SESOI the engine says it only when the corrected effect is shown smaller (journal D7)
    small = _pilot(bgs, sesoi=0.15, delta=0.1)
    assert P.allowed(conds["E3"], "with", hi, small) == {P.SUPPORTED, P.INCONCLUSIVE, P.REFUSAL, P.NDE, P.NS_DEPTH}
    assert P.NS_DEPTH not in P.allowed(conds["E1"], "dose=0.25", hi, small)  # no artifact planted
    # E2 and E3 by their Δ* at the depth of their analysis (the fourth review): above the SESOI at full
    # depth (0.2 against 0.15), below it with both sides at half capture (0.12), NDE is correct there
    halved = _pilot(bgs, sesoi=0.15, delta=0.2)
    assert P.delta_of("E2", "against", hi, halved) is halved["delta"][str(hi)][f"capture={P.E_CAPTURE:g}"]
    assert P.allowed(conds["E1"], "dose=key", hi, halved) == {P.SUPPORTED, P.INCONCLUSIVE, P.REFUSAL}
    for name, variant in (("E2", "against"), ("E3", "with")):
        assert P.allowed(conds[name], variant, hi, halved) == {P.SUPPORTED, P.INCONCLUSIVE, P.REFUSAL, P.NDE,
                                                               P.NS_DEPTH}


def test_the_truth_follows_the_data_of_the_condition(bgs):
    """The first review: the truth measured on N1 does not hold where a condition thins the pair.
    The conditions that change the pair's counts take the truth measured on their own datasets; the
    fourth review: a pilot that measured its background's cases but not one of them is refused (the
    background's null stood in for it silently), and only a stand-in without case records (oc.py's
    scenarios) takes the null's."""
    conds = P.conditions()
    pilot = _pilot(bgs)
    assert "truth_case" not in pilot and P.metric_truth(conds["E1"], "dose=key", 0, pilot) == "valid"
    pilot["truth_case"] = {"B1": {"E2:against": {str(k): {"class": "blind"} for k in range(24)}}}
    assert P.metric_truth(conds["E2"], "against", 0, pilot) == "blind"
    assert P.allowed(conds["E2"], "against", 0, pilot) == P.INVALID_ALLOWED
    with pytest.raises(ValueError, match="no truth for E1:dose=key"):
        P.metric_truth(conds["E1"], "dose=key", 0, pilot)
    assert P.metric_truth(conds["N1"], "null", 0, pilot) == "valid" and "N1" not in P.TRUTH_CASE_CONDITIONS
    # the third review: N4's and N5's designs are smaller than N1's, so GATE 4's odds are their own
    assert {"N4", "N5"} <= set(P.TRUTH_CASE_CONDITIONS)
    pilot["truth_case"]["B1"]["N5:sham"] = {str(k): {"class": "blind"} for k in range(24)}
    assert P.metric_truth(conds["N5"], "sham", 0, pilot) == "blind"


def test_duplicate_gene_symbols_are_made_unique_and_the_engine_runs(tmp_path):
    """The second review's K2: several Ensembl genes can share a symbol; the engine refuses duplicate
    var_names, so every card of such a background would have been an engine error. The loader makes
    the symbols unique (anndata's convention) before the plan, and a card runs."""
    import run_panel as R
    assert P.unique_names(["A", "B", "A", "A", "A-1"]) == ["A", "B", "A-2", "A-3", "A-1"]
    bg = simulated_background("B1", plan=False)
    genes = list(bg.genes)
    for j in (1, 5, 9, 130, 131):  # a coupled pair's genes among them
        genes[j] = "HMGB2"
    P.save_npz(tmp_path / "B1.npz", bg.X, pd.DataFrame({"donor_id": bg.donor}), genes)
    spec = {"B1": dict(file="B1.npz", sha256=P.sha256(tmp_path / "B1.npz"), donor="donor_id", counts="X")}
    (tmp_path / "backgrounds.json").write_text(json.dumps(spec))
    got = P.load_backgrounds(tmp_path / "backgrounds.json", tmp_path)["B1"]
    assert len(set(got.genes)) == len(got.genes)
    bgs = {"B1": got, "B2": simulated_background("B2", donors=12, seed=1)}
    pilot = _pilot(bgs)
    entries = P.assign(KEY, pool_sizes=_sizes(bgs))
    pick = [next(e for e in entries if e["condition"] == "N1")]
    summary = R.run(pick, bgs, pilot, tmp_path / "out", workers=1)
    assert summary["errors"] == 0


def test_the_panel_loads_only_its_own_backgrounds(tmp_path):
    """The first review's K1: a backgrounds.json with the anchors' B3 (three batches of 96 cells, as
    select_backgrounds.fetch_b3 writes it) must not stop the pilot, the preparation or the timing:
    the panel plans B1 and B2 only; anchors.py loads B3 itself."""
    pytest.importorskip("scipy")
    import anchors as A
    b1 = simulated_background("B1", n_genes=240)
    P.save_npz(tmp_path / "B1.npz", b1.X, pd.DataFrame({"donor_id": b1.donor}), list(b1.genes))
    rng = np.random.default_rng(0)
    genes3 = [f"G{j}" for j in range(30)] + [f"ERCC-{j:05d}" for j in range(5)]
    obs3 = pd.DataFrame({"phase": np.repeat(["G1", "S", "G2M"], 96)})
    obs3["batch"] = obs3["phase"]
    P.save_npz(tmp_path / "B3.npz", rng.poisson(2.0, size=(288, 35)), obs3, genes3)
    spec = {"B1": dict(file="B1.npz", sha256=P.sha256(tmp_path / "B1.npz"), donor="donor_id", counts="X"),
            "B3": dict(file="B3.npz", sha256=P.sha256(tmp_path / "B3.npz"), donor="batch", counts="X")}
    (tmp_path / "backgrounds.json").write_text(json.dumps(spec))
    bgs = P.load_backgrounds(tmp_path / "backgrounds.json", tmp_path)
    assert set(bgs) == {"B1"} and len(bgs["B1"].plan["pool"]) >= 12
    X, genes, obs = A._load(tmp_path / "backgrounds.json", "B3", tmp_path)
    assert X.shape == (288, 35) and set(obs["phase"]) == {"G1", "S", "G2M"}


def test_the_pilot_and_the_run_go_through_their_command_lines_on_files(tmp_path, monkeypatch, capsys):
    """The first review's V7 (the pilot had never run as the workflow runs it) and K1: on simulated
    files whose backgrounds.json also names the anchors' B3 (simulate.py write, the dry run's
    selection), the pilot's commands - oracle.py --smoke, timing.py --write-drops, oc.py --pilot -
    and the run's - blind.py prepare, run, collect, and score.py - go through their command lines
    as the workflow calls them. Only the timing's engine runs are stubbed (run_panel.run has its
    own tests); the oracle's sizes are cut to keep the test short."""
    pytest.importorskip("scipy")
    import blind as B
    import oc
    import oracle as O
    import score as S
    import simulate
    import timing as T
    small = {"B1": simulated_background("B1", plan=False),
             "B2": simulated_background("B2", donors=12, seed=1, plan=False)}
    data, spec, pilot_path = tmp_path / "data", tmp_path / "backgrounds.json", tmp_path / "pilot.json"
    simulate.write_backgrounds(data, spec, tmp_path / "selection.md", backgrounds=small)
    assert set(json.loads(spec.read_text())) == {"B1", "B2", "B3"}
    monkeypatch.setattr(O, "SMOKE_SIZES", dict(curve=(1, 2), truth=(2, 3), case=(1, 2), gate4=8, sesoi=4))
    O.main(["--backgrounds", str(spec), "--data-dir", str(data), "--out", str(pilot_path), "--workers", "2",
            "--smoke", "--datasets", "1", "--draws", "10"])
    pilot = json.loads(pilot_path.read_text())
    assert pilot["smoke"] and set(pilot["backgrounds"]) == {"B1", "B2"} and pilot["datasets_per_case"] == 1

    def fake_run(entries, bgs, pilot, out_dir, workers=1):  # 30 s a card, 90 s for N6c
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        rows = [dict(id=e["id"], cards=[dict(id=c) for c in P.card_ids(e)]) for e in entries]
        secs = {c["id"]: 90.0 if e["condition"] == "N6c" else 30.0 for e, r in zip(entries, rows) for c in r["cards"]}
        (out_dir / "manifest.json").write_text(json.dumps(dict(datasets=rows)))
        (out_dir / "runtime.json").write_text(json.dumps(dict(cards=[dict(id=k, seconds=v, peak_rss_mb=500.0)
                                                                   for k, v in secs.items()])))
        return dict(run=len(secs), wall_seconds=sum(secs.values()) / workers, mean_seconds=float(np.mean(list(secs.values()))),
                    median_seconds=30.0, errors=0)
    with monkeypatch.context() as m:  # run_panel is the module blind.py runs too: stub it for the timing only
        m.setattr(T.R, "run", fake_run)
        T.main(["--workers", "4", "--cards", "13", "--backgrounds", str(spec), "--data-dir", str(data),
                "--pilot", str(pilot_path), "--write-drops", "--refusals", "1"])
    timing_out = capsys.readouterr().out
    assert "gate0_refusals" in json.loads(pilot_path.read_text())  # written before the key (the fourth review)
    # every card at its condition's time: (6700 - 790) x 30 s + 790 x 90 s over 20 shards x 4 workers
    expected = ((P.n_cards() - 790) * 30.0 + 790 * 90.0) / 80 / 3600
    assert f"a shard's expected time {expected:.2f} h" in timing_out and "drops written" in timing_out
    assert "# memory: a worker's peak 500 MB x 4 workers = 2000 MB" in timing_out
    assert json.loads(pilot_path.read_text())["dropped"] == []
    oc.main(["--pilot", str(pilot_path), "--write-judged"])
    assert "P(S1-S7 all pass | sound)" in capsys.readouterr().out
    fixed = json.loads(pilot_path.read_text())["s3_rules"]  # S3's tiers, fixed before the key
    assert set(fixed["tiers"]) == set(oc.STRATA) and 0 <= fixed["joint"] <= 1
    B.main(["prepare", "--backgrounds", str(spec), "--data-dir", str(data), "--out", str(tmp_path / "compact")])
    B.main(["run", "--dry-run", "--compact", str(tmp_path / "compact"), "--pilot", str(pilot_path), "--shard", "0",
            "--shards", "1", "--workers", "2", "--limit", "2", "--out", str(tmp_path / "shards" / "shard-0")])
    B.main(["collect", "--shards-dir", str(tmp_path / "shards"), "--out", str(tmp_path / "results")])
    capsys.readouterr()
    S.main(["--key", simulate.DRY_RUN_KEY, "--results", str(tmp_path / "results"), "--pilot", str(pilot_path),
            "--out", str(tmp_path / "scores.json")])
    scores = json.loads((tmp_path / "scores.json").read_text())
    assert scores["n_cards"] == P.n_cards() and set(scores["criteria"]) == set(oc.CRITERIA)
    assert scores["criteria"]["S3"]["tiers_fixed_before_the_key"]
    assert "validation" in capsys.readouterr().out


# --------------------------------------------------------------------------- #
# the oracle
# --------------------------------------------------------------------------- #
def test_oracle_tests_on_replicate_values():
    pytest.importorskip("scipy")
    import oracle as O
    a, b = np.array([8.0, 7, 6, 5]), np.array([4.0, 3, 2, 1])
    assert O.permutation_p(a, b) == pytest.approx(2 / 70)                      # two-sided, exact
    assert O.permutation_p(np.array([8.0, 7, 6, 4]), np.array([5.0, 3, 2, 1])) == pytest.approx(4 / 70)
    assert O.signflip_p(np.array([1.0, 2, 3, 4, 5])) == pytest.approx(2 / 32)
    assert O.tost_width(np.zeros(8) + 0.01 * np.arange(8), np.zeros(8)) < 0.1


def test_the_oracle_decides_in_the_engine_order():
    pytest.importorskip("scipy")
    import oracle as O
    hi, lo = np.array([1.0, 1.1, 1.2, 1.3, 1.05, 1.15]), np.array([0.0, 0.1, 0.2, 0.05, 0.15, 0.1])
    same = np.array([0.5, 0.52, 0.48, 0.51, 0.49, 0.5])
    assert O.outcome_from_values((hi, lo), (same, same + 0.001), True, 0.5, "decrease") == P.NS_DEPTH
    # the corrected effect not shown smaller than the SESOI: not depth, INCONCLUSIVE (journal D7)
    assert O.outcome_from_values((hi, lo), (same, same + 0.001), True, 1e-6, "decrease") == P.INCONCLUSIVE
    assert O.outcome_from_values((hi, lo), (lo, hi), True, 0.5, "decrease") == P.INCONCLUSIVE
    assert O.outcome_from_values((hi, lo), (hi, lo), False, 0.5, "decrease") == P.SUPPORTED
    assert O.outcome_from_values((hi, lo), (hi, lo), False, 0.5, "increase") == P.NS_OPPOSITE
    assert O.outcome_from_values((same, same), (same, same), False, 0.5, "increase") == P.NDE
    assert O.outcome_from_values((same, same), (same, same), False, 0.001, "increase") == P.INCONCLUSIVE


def test_gate4_odds_and_outcome_follow_the_rule():
    """GATE 4's rule on responses (FAIL first, on the upper bound below delta_min; PASS on the lower
    bound above 0; else UNTESTED) and its odds from a case's response statistics."""
    pytest.importorskip("scipy")
    import oracle as O
    rng = np.random.default_rng(0)
    assert O.gate4_outcome(rng.normal(0.3, 0.05, 200), 0.05) == "PASS"
    assert O.gate4_outcome(rng.normal(0.0, 0.01, 200), 0.05) == "FAIL"
    assert O.gate4_outcome(rng.normal(0.02, 0.01, 200), 0.05) == "FAIL"     # inside (0, delta_min): blind
    assert O.gate4_outcome(rng.normal(0.0, 0.5, 200), 0.05) == "UNTESTED"
    far = O.gate4_odds(0.3, 0.01, 0.05, 0.05)
    assert far["p_pass"] > 0.999 and far["p_fail"] < 1e-6
    zero = O.gate4_odds(0.0, 0.0, 0.01, 0.05)
    assert zero["p_fail"] > 0.999
    near = O.gate4_odds(0.06, 0.02, 0.1, 0.05)                              # valid, but datasets vary
    assert 0.05 < near["p_fail"] < 0.5 and near["p_pass"] + near["p_fail"] + near["p_untested"] == pytest.approx(1)


def test_the_sound_validator_model_takes_every_open_route_at_its_measured_size(bgs):
    """The first review's V4: a card's rates are what a sound validator following the engine's rules
    gets on its case, from the pilot's measurements, in the engine's verdict order: on a valid metric
    GATE 5's designed 1.5 alpha, then GATE 4's measured odds of FAILing, then the effect's measured
    outcomes after a PASS; on a blind metric GATE 5 too, then the engine's rules on the case's
    datasets (a PASS then a NO DETECTABLE EFFECT is an error); without measurements, the designed
    sizes. The second review: after an UNTESTED GATE 4 an effect explained by depth is still called
    (the engine checks it before the metric's validity), and the thresholds use, per card, the larger
    of the rate and the designed size."""
    import oc
    conds = P.conditions()
    pilot = _pilot(bgs, truth={"low": "blind"})
    pilot["truth"]["B1"]["0"]["gate4"] = dict(p_pass=0.9, p_fail=0.08, p_untested=0.02)
    pilot["establishable"] = {
        "N1:null:0": dict(establishable=True, effect_outcomes={P.NDE: 90, P.SUPPORTED: 3, P.NS_OPPOSITE: 2,
                                                               P.INCONCLUSIVE: 5}),
        "N1:null:20": dict(establishable=True, engine_outcomes={P.NS_INVALID: 80, P.NDE: 15, P.INCONCLUSIVE: 5})}
    m = oc.sound_model(conds["N1"], "null", 0, pilot)
    g5 = oc.E_GATE5
    reach = 0.9 * (1 - g5)
    assert m["measured"] and m["p_inv"] == pytest.approx(g5 + (1 - g5) * 0.08)
    assert m["p_sup"] == pytest.approx(reach * 0.03) and m["p_err"] == pytest.approx(m["p_inv"] + reach * 0.03)
    assert m["decisive"] == pytest.approx(reach * 0.92) and m["p_nde"] == 0.0  # NDE is allowed on a null
    row = oc.card_model(conds["N1"], "null", 0, pilot)
    assert row["r_sup"] == oc.E_SUPPORTED and row["r_nde"] == 0.0  # the designed size above the rate
    assert row["r_inv"] == pytest.approx(m["p_inv"])                   # the rate above the designed size
    assert row["r_opp"] == 0.0 and row["r_depth"] == oc.E_DEPTH        # against the direction is right on a null;
    assert row["data_stratum"] == "null" and row["s2_group"] is None   # depth is wrong without an artifact
    # against the direction on a null: allowed, its rate bounded as a false detection's (the fourth review)
    assert m["opp_null_counted"] and m["p_oppnull"] == pytest.approx(reach * 0.02)
    assert row["r_oppnull"] == oc.E_OPPOSITE                           # the designed alpha/2 above the rate
    b = oc.sound_model(conds["N1"], "null", 20, pilot)                     # blind: PASS then NDE is an error
    assert b["p_err"] == pytest.approx((1 - g5) * 0.15) and b["p_nde"] == pytest.approx((1 - g5) * 0.15)
    assert b["p_inv"] == 0 and b["decisive"] == pytest.approx(g5 + (1 - g5) * 0.80)
    assert not b["opp_null_counted"] and b["p_oppnull"] == 0.0  # blind: an error there, counted as one (S4 invalid)
    t = oc.sound_model(conds["N1"], "null", 1, pilot)                      # nothing measured: the designed sizes
    assert not t["measured"] and t["p_err"] == pytest.approx(0.175) and t["decisive"] == pytest.approx(0.825)
    assert t["p_sup"] == oc.E_SUPPORTED and t["p_inv"] == oc.E_INVALID and t["p_depth"] == oc.E_DEPTH
    assert t["p_oppnull"] == oc.E_OPPOSITE
    # an UNTESTED GATE 4 with the effect explained by depth: NOT SUPPORTED (explained by depth), an
    # error on a real effect above the SESOI (E2 at the key dose), not INCONCLUSIVE
    pilot["truth_case"] = {"B1": {"E2:against": {"0": dict(pilot["truth"]["B1"]["0"],
                                                           gate4=dict(p_pass=0.0, p_fail=0.0, p_untested=1.0))}}}
    pilot["establishable"]["E2:against:0"] = dict(establishable=True,
                                                  effect_outcomes={P.NS_DEPTH: 40, P.SUPPORTED: 60})
    e = oc.sound_model(conds["E2"], "against", 0, pilot)
    assert P.NS_DEPTH not in P.allowed(conds["E2"], "against", 0, pilot)
    assert e["p_err"] == pytest.approx(g5 + (1 - g5) * 0.4) and e["decisive"] == 0.0
    assert oc.stratum(conds["E1"], "valid") == "effect" and oc.stratum(conds["N1"], "valid") == "null"
    assert oc.stratum(conds["N8"], "valid") is None and oc.stratum(conds["N3"], "blind") is None  # the oracle's
    assert oc.stratum(conds["N6c"], "useless") == "invalid" and oc.stratum(conds["N1"], "ambiguous") is None


def test_gate0s_measured_refusals_are_in_the_sound_validators_model(bgs, tmp_path, monkeypatch):
    """The fourth review's V3: GATE 0 was left out of the sound validator's model, so at a 15% refusal
    share of a stratum a sound engine failed S3. timing.py --refusals runs the frozen engine before the
    key on establishable cards of every S3 stratum (a public key's datasets; not the constant metric,
    whose DEGENERATE METRIC comes first) and writes the share per stratum into pilot.json; the model
    then refuses first with it, and every other outcome takes the rest."""
    import oc
    import timing as T
    conds = P.conditions()
    pilot = _est_pilot(bgs)
    entries = T.refusal_sample(pilot, 5, np.random.default_rng(1))
    by = Counter(e["stratum"] for e in entries)
    assert set(by) == set(oc.STRATA) and all(v == 5 for v in by.values())
    assert all(conds[e["condition"]].metric != "constant" and conds[e["condition"]].oracle for e in entries)

    def fake_run(es, bgs_, pilot_, out_dir, workers=1):  # GATE 0 refuses every second card
        (Path(out_dir) / "reports").mkdir(parents=True, exist_ok=True)
        for i, e in enumerate(es):
            rep = (dict(verdict="NOT SUPPORTED — metric invalid: a nuisance bias", cause="metric_invalid_gate0")
                   if i % 2 == 0 else dict(verdict="INCONCLUSIVE — the effect is inconclusive", cause="effect_inconclusive"))
            (Path(out_dir) / "reports" / f"{e['id']}.json").write_text(json.dumps(rep))
    monkeypatch.setattr(T.R, "run", fake_run)
    rec = T.measure_refusals(entries, bgs, pilot, tmp_path / "refusals", workers=1)
    assert sum(r["n"] for r in rec.values()) == len(entries) == 15
    assert sum(r["refused"] for r in rec.values()) == 8 and all(r["share"] == r["refused"] / r["n"] for r in rec.values())
    before = {c: oc.sound_model(conds[c], v, 0, pilot) for c, v in (("N1", "null"), ("N6a", "random"))}
    pilot["gate0_refusals"] = {"null": dict(n=10, refused=2, share=0.2), "invalid": dict(n=10, refused=5, share=0.5)}
    after = oc.sound_model(conds["N1"], "null", 0, pilot)
    assert after["refusal"] == 0.2 and after["decisive"] == pytest.approx(0.8 * before["N1"]["decisive"])
    assert after["p_err"] == pytest.approx(0.8 * before["N1"]["p_err"])
    assert oc.sound_model(conds["N6a"], "random", 0, pilot)["decisive"] == pytest.approx(0.5 * before["N6a"]["decisive"])
    assert oc.sound_model(conds["E1"], "dose=key", 0, pilot)["refusal"] == 0.0     # not measured there
    assert oc.sound_model(conds["N6b"], "constant", 0, pilot)["refusal"] == 0.0    # DEGENERATE METRIC first
    rows = oc.expected_rows(pilot)
    null = [r for r in rows if r["stratum"] == "null" and r["establishable"]]
    assert null and all(r["refusal"] == 0.2 for r in null)


def test_the_oracle_applies_the_planted_nuisance_to_the_side_that_lacks_it(bgs):
    import oracle as O
    pilot = _pilot(bgs)
    for name, variant, lacks in (("N2", "c=0.5", "A"), ("N8", "beta(2,2)", "A"), ("E2", "against", "A"),
                                 ("E3", "with", "B")):
        e = _entry(name, variant, side="B", seed=5)
        X, obs, _, _ = P.build(e, bgs, pilot)
        fixed = O.true_correction(X, obs, e, np.random.default_rng(0))
        g = np.asarray(obs["group"])
        for side in ("A", "B"):
            changed = not np.array_equal(fixed[g == side], X[g == side])
            assert changed == (side == lacks), (name, side)
    assert O.true_correction(*P.build(_entry("N1", "null"), bgs, pilot)[:2], _entry("N1", "null"),
                             np.random.default_rng(0)) is None


def test_the_oracle_shows_blindness_where_the_metric_is_blind_and_useless(bgs):
    """For a blind pair or a useless metric the correct definite outcome is metric invalid, which
    takes showing blindness: the low-level pairs of the simulated background (where log-normalized
    Pearson does not respond, probe p14) and the random-gene score facing the G2M module are shown
    blind; a random number is not (its noise keeps the interval wide); the constant is DEGENERATE."""
    pytest.importorskip("scipy")
    import oracle as O
    pilot = _pilot(bgs, sesoi=0.1, dose=2.25, truth={"low": "blind"})
    low = next(pe["index"] for pe in bgs["B1"].plan["pool"] if pe["level"] == "low")
    outs = []
    for s in range(4):
        e = _entry("N1", "null", seed=s, i=s, pair=low)
        X, obs, _, cards = P.build(e, bgs, pilot)
        outs.append(O.oracle_outcome(e, X, obs, cards[0], bgs, pilot, np.random.default_rng(s)))
    assert outs.count(P.NS_INVALID) >= 3
    e = _entry("N6c", "random-genes", seed=1)
    X, obs, _, cards = P.build(e, bgs, pilot)
    assert O.oracle_outcome(e, X, obs, cards[0], bgs, pilot, np.random.default_rng(0)) == P.NS_INVALID
    assert O.module_deltas(X, [bgs["B1"].genes.index(g) for g in cards[0]["score_genes"]],
                           [bgs["B1"].genes.index(g) for g in bgs["B1"].plan["g2m"]], np.random.default_rng(0),
                           50).mean() < 0
    e = _entry("N6a", "random", seed=1)
    X, obs, _, cards = P.build(e, bgs, pilot)
    assert O.oracle_outcome(e, X, obs, cards[0], bgs, pilot, np.random.default_rng(0)) == P.INCONCLUSIVE
    e = _entry("N6b", "constant", seed=1)
    X, obs, _, cards = P.build(e, bgs, pilot)
    assert O.oracle_outcome(e, X, obs, cards[0], bgs, pilot, np.random.default_rng(0)) == P.DEGENERATE


def test_the_saturation_dose_is_the_smallest_dose_reaching_95_percent_of_the_maximum():
    import oracle as O
    curve = {f"{d:g}": dict(response=r) for d, r in ((0.5, 0.05), (1.0, 0.15), (2.0, 0.24), (2.5, 0.255),
                                                      (3.0, 0.256), (4.0, 0.25))}
    assert O.saturation_dose(curve) == 2.5
    assert O.saturation_dose({f"{d:g}": dict(response=-0.01) for d in (0.5, 1.0, 4.0)}) == 4.0


def test_delta_star_is_a_paired_estimate_with_a_small_standard_error(bgs):
    """And E2's and E3's at the depth of their analysis, both sides at half capture (the fourth
    review): smaller than at full depth, the coupling's correlation diluted by the thinning."""
    import oracle as O
    d = O.delta_star(bgs["B1"], 0, 3.0, draws=200)
    assert d["value"] > 0.05 and d["se"] < d["value"] / 5 and d["capture"] == 1.0
    assert O.delta_star(bgs["B1"], 0, 3.0, draws=200) == d  # public seed: reproducible
    half = O.delta_star(bgs["B1"], 0, 3.0, draws=200, capture=P.E_CAPTURE)
    assert half["capture"] == P.E_CAPTURE and 0 < half["value"] < d["value"] - 2 * math.hypot(d["se"], half["se"])


@pytest.fixture(scope="module")
def small_bgs():
    """Backgrounds with 4 pairs per level, for the pilot's run in a test."""
    old = P.MAX_PAIRS_PER_LEVEL
    P.MAX_PAIRS_PER_LEVEL = 4
    try:
        return {"B1": simulated_background("B1", n_genes=300), "B2": simulated_background("B2", donors=12, n_genes=300, seed=1)}
    finally:
        P.MAX_PAIRS_PER_LEVEL = old


def test_the_pilot_fixes_sesoi_saturation_truth_doses_and_establishable_cases(small_bgs, monkeypatch):
    pytest.importorskip("scipy")
    import oracle as O
    monkeypatch.setattr(O, "DOSE_GRID", (1.0, 2.0, 3.0))
    pilot = O.run_pilot(small_bgs, n=3, draws=20, sizes=dict(curve=(1, 2), truth=(3, 4), case=(2, 3), gate4=20, sesoi=6))
    for b in ("B1", "B2"):
        assert pilot["pool_size"][b] == 12 and [pe["level"] for pe in pilot["pool"][b]] == [
            lv for lv in P.LEVELS for _ in range(4)]
        for lv in P.LEVELS:
            assert pilot["sesoi"][b][lv] in O.SESOI_GRID
            assert pilot["sesoi_found"][b][lv] or pilot["sesoi"][b][lv] == O.SESOI_GRID[-1]
            assert pilot["saturation_dose"][b][lv] in O.DOSE_GRID
        for k, t in pilot["truth"][b].items():
            assert t["class"] == P.classify_response(t["response"], t["delta_min"])
            assert t["delta_min"] == pytest.approx(P.DELTA_MIN_FRACTION * pilot["sesoi"][b][pilot["pool"][b][int(k)]["level"]])
    for lv in P.LEVELS:
        if pilot["key_dose_found"][lv]:
            last = pilot["dose_scan"][lv][-1]
            assert pilot["e_dose"][lv] == pilot["key_dose"][lv] == last["dose"]
            assert last["power"] >= O.POWER and last["delta"] >= O.DELTA_MARGIN * pilot["sesoi"]["B1"][lv]
        else:  # no key dose: the real effects use the saturation dose
            assert pilot["key_dose"][lv] is None and pilot["e_dose"][lv] == pilot["saturation_dose"]["B1"][lv]
        assert 0 <= pilot["n2_raw_power"][lv] <= 1
    assert set(pilot["delta"]) == {str(k) for k in range(12)}
    for k, pe in enumerate(pilot["pool"]["B1"]):  # E2's and E3's Δ* at the key dose and half capture
        half = pilot["delta"][str(k)][f"capture={P.E_CAPTURE:g}"]
        assert half["capture"] == P.E_CAPTURE and half["dose"] == pytest.approx(pilot["e_dose"][pe["level"]])
    keys = {f"{c.name}:{v}:{k}" for c in P.CONDITIONS for v, _ in c.variants for k in range(12)}
    assert set(pilot["establishable"]) == keys
    assert not pilot["establishable"]["N4:3v3:0"]["establishable"]
    assert pilot["establishable"]["N6b:constant:5"]["establishable"]  # always DEGENERATE METRIC
    assert pilot["backgrounds"] == {k: P.background_sha256(b) for k, b in small_bgs.items()}
    # the truth on the data of every condition that changes the pair's counts, with GATE 4's odds
    cases = {f"{c.name}:{v}" for c in P.CONDITIONS if c.name in P.TRUTH_CASE_CONDITIONS for v, _ in c.variants}
    assert set(pilot["truth_case"]["B1"]) == cases
    for rec in list(pilot["truth_case"]["B1"]["E2:against"].values()) + list(pilot["truth"]["B1"].values()):
        assert rec["class"] == P.classify_response(rec["response"], rec["delta_min"])
        odds = rec["gate4"]
        assert odds["p_pass"] + odds["p_fail"] + odds["p_untested"] == pytest.approx(1.0)
    # what the sound-validator model reads: the effect's outcomes, and the engine's rules where GATE 4 ran
    case = pilot["establishable"]["N1:null:11"]
    assert sum(case["effect_outcomes"].values()) == 3
    for key, case in pilot["establishable"].items():
        if case.get("gate4_outcomes"):
            assert sum(case["engine_outcomes"].values()) == sum(case["gate4_outcomes"].values())


# --------------------------------------------------------------------------- #
# scoring and the criteria
# --------------------------------------------------------------------------- #
CAUSE_OF = {P.SUPPORTED: ("SUPPORTED (provisional until replicated)", "provisional"),
            P.NDE: ("NO DETECTABLE EFFECT — x", "no_detectable_effect"),
            P.INCONCLUSIVE: ("INCONCLUSIVE — x", "effect_inconclusive"),
            P.NS_INVALID: ("NOT SUPPORTED — metric invalid: x", "metric_invalid_gate4"),
            P.NS_DEPTH: ("NOT SUPPORTED — x", "explained_by_depth"),
            P.NS_OPPOSITE: ("NOT SUPPORTED — x", "opposite_direction"),
            P.DEGENERATE: ("DEGENERATE METRIC — x", "degenerate_metric"),
            P.REFUSAL: ("NOT SUPPORTED — metric invalid: bias", "metric_invalid_gate0")}


def test_a_background_without_a_candidate_is_dropped_with_its_cases(small_bgs, monkeypatch):
    """The second review's V6: v1.md 3.1 drops a background without a qualifying candidate with the
    cases that need it, but without B2 the pilot stopped (KeyError 'B2'). The pilot records 'B2:*',
    which drops N7 from the key's assignment, the expected cards, the timing and the scoring; the
    order's drops come after it. Without B1 there is no validation to run."""
    pytest.importorskip("scipy")
    import oc
    import oracle as O
    import timing as T
    monkeypatch.setattr(O, "DOSE_GRID", (1.0, 2.0, 3.0))
    pilot = O.run_pilot({"B1": small_bgs["B1"]}, n=2, draws=10,
                        sizes=dict(curve=(1, 2), truth=(2, 3), case=(2, 2), gate4=10, sesoi=4))
    assert pilot["dropped"] == ["B2:*"] and "B2" not in pilot["pool_size"]
    assert not any(k.startswith("N7:") for k in pilot["establishable"])
    entries = P.assign(KEY, pilot["dropped"], pilot["pool_size"])
    assert {e["condition"] for e in entries} == {c.name for c in P.CONDITIONS} - {"N7"}
    n7 = sum(n for _, n in P.conditions()["N7"].variants)
    assert len(entries) == P.n_datasets() - n7 == P.n_datasets(pilot["dropped"]) and n7 == 300
    assert not any(r["condition"] == "N7" for r in oc.expected_rows(pilot, pilot["dropped"]))
    slow = {c.name: 160.0 for c in P.CONDITIONS}  # slow enough that some of the order's drops are needed
    order = T.decide_drops(slow, dropped=pilot["dropped"])
    assert order and P.check_dropped(pilot["dropped"] + order) == ["B2:*"] + order
    assert T.shard_hours(slow, ["B2:*"] + order) <= T.BUDGET_HOURS
    with pytest.raises(ValueError):
        P.check_dropped(["N3:f=0.1", "B2:*"])
    with pytest.raises(ValueError):
        P.check_dropped(["B1:*"])
    with pytest.raises(SystemExit):
        P.background_drops({"B2"})


def _est_pilot(bgs, delta=0.3, truth=None):
    pilot = _pilot(bgs, delta=delta, truth=truth)
    pilot["establishable"] = {f"{c.name}:{v}:{k}": dict(establishable=c.oracle)
                              for c in P.CONDITIONS for v, _ in c.variants for k in range(24)}
    return pilot


def _reports(entries, fn):
    """Reports by card id; fn(condition, entry) gives one outcome per dataset, or fn(condition, entry,
    card id) one per card where it takes three arguments."""
    import inspect
    conds = P.conditions()
    per_card = len(inspect.signature(fn).parameters) == 3
    return {cid: CAUSE_OF[fn(conds[e["condition"]], e, cid) if per_card else fn(conds[e["condition"]], e)]
            for e in entries for cid in P.card_ids(e)}


ORDER = (P.NS_INVALID, P.SUPPORTED, P.NS_DEPTH, P.NDE, P.NS_OPPOSITE, P.DEGENERATE)


def _card(e, cid):
    return P.card_ids(e).index(cid)


def _perfect(pilot):
    def fn(c, e, cid):  # a correct definite outcome of the card (what the engine's rules exclude on it aside)
        good = P.card_allowed(c, e["variant"], e["pair"], pilot, _card(e, cid)) & P.DEFINITE
        return next((o for o in ORDER if o in good), P.INCONCLUSIVE)
    return fn


def _bad_validators(pilot):
    """The bad validators of v1.md section 6, with the criteria each must fail (at least)."""
    perfect = _perfect(pilot)

    def no_supported(c, e, cid):  # the first review's: perfect, but never SUPPORTED
        o = perfect(c, e, cid)
        return P.INCONCLUSIVE if o == P.SUPPORTED else o

    def nde_not_supported(c, e, cid):  # the second review's: perfect, but NO DETECTABLE EFFECT for SUPPORTED
        o = perfect(c, e, cid)
        return P.NDE if o == P.SUPPORTED else o

    def blind_or_nde(c, e):  # shows blindness where it can, else NO DETECTABLE EFFECT, never tests the effect
        return P.NS_INVALID if P.NS_INVALID in P.definite(c, e["variant"], e["pair"], pilot) else P.NDE

    def blind_or_unsure(c, e):  # shows blindness where it can, else INCONCLUSIVE
        return P.NS_INVALID if P.NS_INVALID in P.definite(c, e["variant"], e["pair"], pilot) else P.INCONCLUSIVE

    def settled(c, e):
        est = (pilot["establishable"].get(f"{c.name}:{e['variant']}:{e['pair']}") or {}).get("establishable")
        return bool(c.oracle and est)

    def realistic(c, e, cid):  # correct definite where the case is establishable, else INCONCLUSIVE
        return perfect(c, e, cid) if settled(c, e) else P.INCONCLUSIVE

    def nde_on_unestablishable_blind(c, e, cid):  # the second review's V8: NO DETECTABLE EFFECT on a quarter of
        truth = P.metric_truth(c, e["variant"], e["pair"], pilot)  # the blind or useless cards it cannot settle
        if truth in ("blind", "useless") and not settled(c, e) and int(e["id"][1:]) % 4 == 0:
            return P.NDE
        return perfect(c, e, cid)

    def opposite_for_unsure(c, e, cid):  # the third review's X1: NOT SUPPORTED (opposite direction) where a
        o = realistic(c, e, cid)          # sound validator would say INCONCLUSIVE or NO DETECTABLE EFFECT
        return P.NS_OPPOSITE if o in (P.INCONCLUSIVE, P.NDE) else o

    def nde_for_unsure(c, e, cid):  # the third review's X4: NO DETECTABLE EFFECT instead of INCONCLUSIVE
        o = realistic(c, e, cid)
        return P.NDE if o == P.INCONCLUSIVE else o

    def pseudoreplication(c, e, cid):  # the third review's X2b: SUPPORTED on N4's card without the replicate unit
        return P.SUPPORTED if c.name == "N4" and cid.endswith("a") else perfect(c, e, cid)

    def supported_on_n7(c, e, cid):  # the third review's X3: SUPPORTED on a quarter of N7
        return P.SUPPORTED if c.name == "N7" and int(e["id"][1:]) % 4 == 0 else perfect(c, e, cid)

    def opposite_on_nulls(c, e, cid):  # the fourth review's V1: NOT SUPPORTED against the direction on a sixth
        o = perfect(c, e, cid)          # of the null cards where it is allowed (16% in its demonstration)
        allowed = P.card_allowed(c, e["variant"], e["pair"], pilot, _card(e, cid))
        null = c.data == "null" and P.metric_truth(c, e["variant"], e["pair"], pilot) in ("valid", "ambiguous")
        return P.NS_OPPOSITE if null and P.NS_OPPOSITE in allowed and int(e["id"][1:]) % 6 == 0 else o

    def supported_on_n6a(c, e, cid):  # the fourth review's V2: SUPPORTED on an eighth of N6a (12% passed)
        return P.SUPPORTED if c.name == "N6a" and int(e["id"][1:]) % 8 == 0 else perfect(c, e, cid)

    def pseudoreplication_tenth(c, e, cid):  # the fourth review's V2: SUPPORTED on a tenth of N4's cards
        if c.name == "N4" and cid.endswith("a") and int(e["id"][1:]) % 10 == 0:  # without the replicate unit
            return P.SUPPORTED
        return perfect(c, e, cid)
    return (("always invalid", lambda c, e: P.NS_INVALID, {"S3", "S5"}),
            ("always inconclusive", lambda c, e: P.INCONCLUSIVE, {"S3"}),
            ("always refuses", lambda c, e: P.REFUSAL, {"S3"}),
            ("always supported", lambda c, e: P.SUPPORTED, {"S1", "S2"}),
            ("always no effect", lambda c, e: P.NDE, {"S3", "S6"}),
            ("perfect but never SUPPORTED", no_supported, {"S3"}),
            ("perfect but NO DETECTABLE EFFECT for SUPPORTED", nde_not_supported, {"S3"}),
            ("blindness, else NO DETECTABLE EFFECT", blind_or_nde, {"S3"}),
            ("blindness, else INCONCLUSIVE", blind_or_unsure, {"S3"}),
            ("NO DETECTABLE EFFECT on a quarter of the unsettled blind or useless cards",
             nde_on_unestablishable_blind, {"S6"}),
            ("NOT SUPPORTED against the direction where unsure (X1)", opposite_for_unsure, {"S4"}),
            ("NO DETECTABLE EFFECT where unsure (X4)", nde_for_unsure, {"S6"}),
            ("pseudoreplication on N4 (X2b)", pseudoreplication, {"S7"}),
            ("SUPPORTED on a quarter of N7 (X3)", supported_on_n7, {"S2"}),
            ("NOT SUPPORTED against the direction on a sixth of the nulls (V1)", opposite_on_nulls, {"S4"}),
            ("SUPPORTED on an eighth of N6a (V2)", supported_on_n6a, {"S2"}),
            ("SUPPORTED on a tenth of N4 without the replicate unit (V2)", pseudoreplication_tenth, {"S7"}))


UNSETTLED_ONLY = ("NO DETECTABLE EFFECT on a quarter", "NOT SUPPORTED against the direction where unsure",
                  "NO DETECTABLE EFFECT where unsure")


def _failed(res):
    import oc
    return {s for s in oc.CRITERIA if not res["criteria"][s]["passed"]}


def test_score_passes_a_perfect_engine_and_fails_the_bad_validators(bgs):
    """A validator that always gives a correct definite outcome passes; one that always rejects
    the metric, always says INCONCLUSIVE, always refuses, always says SUPPORTED or always NO
    DETECTABLE EFFECT fails, and so do the reviews' perfect-but validators (the reviewer's question:
    no bad validator passes). An engine that fails on cards fails S7, however few."""
    pytest.importorskip("scipy")
    import score as S
    pilot = _est_pilot(bgs, truth={"low": "blind"})
    entries = P.assign(KEY, pool_sizes=_sizes(bgs))
    perfect = _perfect(pilot)
    res = S.score(entries, _reports(entries, perfect), pilot, sims=2000)
    assert res["passed"] and not _failed(res)
    assert res["joint_pass_probability_sound"] >= 0.9
    for name, fn, fails in _bad_validators(pilot):
        if name.startswith(UNSETTLED_ONLY):
            continue  # every case of this pilot is establishable: see the scenario tests
        got = _failed(S.score(entries, _reports(entries, fn), pilot, sims=2000))
        assert fails <= got, (name, got)
    got = S.score(entries, _reports(entries, _bad_validators(pilot)[5][1]), pilot, sims=2000)["criteria"]["S3"]["strata"]
    assert not got["effect"]["passed"] and got["null"]["passed"] and got["invalid"]["passed"]
    nothing = S.score(entries, {}, pilot, sims=2000)                    # no reports: every card an error
    assert nothing["pooled"]["errors"]["rate"] == 1.0 and nothing["criteria"]["S7"]["rate"] == 1.0
    reports = _reports(entries, perfect)
    lost = [e["id"] for e in entries if e["condition"] == "N1"][0]
    reports.pop(lost)                                     # one missing report among 6,700 cards
    res = S.score(entries, reports, pilot, sims=2000)
    assert _failed(res) == {"S7"} and res["criteria"]["S7"]["k"] == 1 and not res["passed"]


@pytest.mark.parametrize("scenario", [(True, 1.0, 0.25), (True, 0.5, 0.5)], ids=["C", "A"])
def test_the_bad_validators_fail_where_the_real_effects_are_few(scenario):
    """The second review's K1: in scenario C of oc.log (the low and medium levels blind, a quarter of
    the cases establishable) only 81 real-effect cards are establishable; with that stratum reported
    only, a validator that never says SUPPORTED, one that shows blindness where it can and otherwise
    says NO DETECTABLE EFFECT, and one that says NO DETECTABLE EFFECT for SUPPORTED passed every
    criterion; and in scenario A one that says NO DETECTABLE EFFECT on a quarter of the blind or
    useless cards it cannot settle, and an engine failing on half of N1's cards (V8). Every stratum is
    judged now (in C the real effects at the strictest tier the joint requirement allows on the key's
    cards), S6 counts a false NO DETECTABLE EFFECT and S7 an engine error: they all fail, and a
    perfect validator passes."""
    pytest.importorskip("scipy")
    import oc
    import score as S
    pilot = oc.scenario_pilot(*scenario)
    entries = P.assign(KEY, pool_sizes=pilot["pool_size"])
    res = S.score(entries, _reports(entries, _perfect(pilot)), pilot, sims=4000)
    assert not _failed(res) and res["joint_pass_probability_sound"] >= 0.90
    assert res["criteria"]["S3"]["strata"]["effect"]["tier"] in oc.TIERS
    for name, fn, fails in _bad_validators(pilot):
        got = _failed(S.score(entries, _reports(entries, fn), pilot, sims=4000))
        assert fails <= got, (name, got)
    half = set(sorted(e["id"] for e in entries if e["condition"] == "N1")[::2])
    reports = {cid: rec for cid, rec in _reports(entries, _perfect(pilot)).items() if cid not in half}
    assert _failed(S.score(entries, reports, pilot, sims=4000)) == {"S7"}  # crashing on half of N1: S4 absorbed it


def test_s3_counts_only_correct_definite_outcomes_and_nde_on_a_large_effect_is_an_error(bgs):
    pytest.importorskip("scipy")
    import score as S
    pilot = _est_pilot(bgs, delta=0.3)  # |Δ*| >= SESOI at the key dose
    entries = P.assign(KEY, pool_sizes=_sizes(bgs))
    res = S.score(entries, _reports(entries, lambda c, e: P.NDE if c.metric == "norm_pearson" else P.NS_INVALID), pilot)
    e1 = res["per_condition"]["E1:dose=key"]
    assert e1["errors"]["rate"] == 1.0 and e1["correct_definite"]["rate"] == 0.0
    assert e1["outcomes"][P.NDE]["rate"] == 1.0  # the dose-verdict curves: every outcome's share
    assert res["per_condition"]["E1:dose=0.25"]["errors"]["rate"] == 0.0  # 0.075 < SESOI 0.15
    e3 = res["criteria"]["S3"]["strata"]["effect"]
    assert e3["definite_any"]["rate"] > e3["rate"] and not res["criteria"]["S3"]["passed"]


def test_s1_fails_when_any_one_key_condition_fails(bgs):
    pytest.importorskip("scipy")
    import score as S
    pilot = _est_pilot(bgs)
    entries = P.assign(KEY, pool_sizes=_sizes(bgs))
    bad = {e["id"] for e in [e for e in entries if e["condition"] == "N8"][:40]}  # 40/790 > 29
    res = S.score(entries, _reports(entries, lambda c, e: P.SUPPORTED if e["id"] in bad else P.INCONCLUSIVE), pilot)
    assert not res["criteria"]["S1"]["conditions"]["N8"]["passed"] and not res["criteria"]["S1"]["passed"]
    assert all(v["passed"] for k, v in res["criteria"]["S1"]["conditions"].items() if k != "N8")


def test_s5_counts_false_invalid_only_where_the_metric_is_valid(bgs):
    pytest.importorskip("scipy")
    import score as S
    pilot = _est_pilot(bgs, truth={"low": "blind", "medium": "ambiguous"})
    entries = P.assign(KEY, pool_sizes=_sizes(bgs))
    res = S.score(entries, _reports(entries, lambda c, e: P.NS_INVALID), pilot)
    s5 = res["criteria"]["S5"]  # per stratum since the third review: real effects, nulls
    n_valid = sum(P.metric_truth(P.conditions()[e["condition"]], e["variant"], e["pair"], pilot) == "valid"
                  for e in entries for _ in P.card_ids(e))
    assert set(s5["cells"]) == {"S5:effect", "S5:null"}
    assert sum(c["n"] for c in s5["cells"].values()) == n_valid == sum(c["k"] for c in s5["cells"].values())
    assert not s5["passed"] and not any(c["passed"] for c in s5["cells"].values() if c["n"])
    assert res["per_truth"]["blind"]["errors"]["rate"] == 0.0 and res["per_truth"]["ambiguous"]["errors"]["rate"] == 0.0


def test_the_frozen_oc_log_is_what_oc_py_prints(capsys):
    """oc.log is frozen with the tag and v1.md quotes it: it must be what the code prints (the fourth
    review checked it by hand)."""
    pytest.importorskip("scipy")
    import oc
    oc.main([])
    assert capsys.readouterr().out == (HERE / "oc.log").read_text()


def test_the_criteria_meet_the_principle_and_s1_as_a_whole():
    pytest.importorskip("scipy")
    import oc
    k, p_all, p_doubled, _ = oc.joint_error_rule(790, 0.025, 5)
    assert k == 29 and p_all >= 0.90 and p_doubled <= 0.05
    assert oc.joint_error_rule(600, 0.025, 5)[1] < 0.90  # 600 per condition fails S1 as a whole
    for n, p0 in ((5100, oc.E_SUPPORTED), (6000, 0.1), (3000, oc.E_INVALID)):
        _, ps, pd = oc.error_rule(n, p0)
        assert ps >= 0.90 and pd <= 0.05
    _, ps, pd = oc.decisiveness_rule(1000, oc.D_NOMINAL)
    assert ps >= 0.90 and pd <= 0.05
    assert oc.nominal_error(frozenset({P.NDE, P.INCONCLUSIVE, P.NS_OPPOSITE, P.REFUSAL})) == pytest.approx(0.175)
    assert oc.nominal_error(P.INVALID_ALLOWED) == pytest.approx(0.15)  # SUPPORTED, NDE, opposite direction, depth
    assert oc.nominal_error(frozenset({P.SUPPORTED, P.INCONCLUSIVE, P.REFUSAL})) == pytest.approx(0.225)  # + invalid
    for scenario in ((True, 0.5, 0.5), (False, 0.0, 1.0), (True, 1.0, 0.25)):
        rules = oc.criteria_rules(oc.expected_rows(oc.scenario_pilot(*scenario)), sims=4000)
        assert rules["joint"] >= 0.90 and rules["s3_feasible"], scenario
        for name, r in rules["cells"].items():  # the second review's V4: every error cell meets the principle,
            if not r["n"]:                       # and (the third's cells) passes a sound validator with >= 0.9995
                continue
            assert r["judged"] and r["p_pass_sound"] >= 0.90 and r["p_pass_doubled"] <= 0.05, (scenario, name)
            assert r["p_pass_measured"] >= 1 - oc.EPS_CELL, (scenario, name)
    # where the designed sizes leave less room than D_NOMINAL, the fallback's decisiveness is what is left
    assert oc.sound_model(P.conditions()["E1"], "dose=key", 0, oc.scenario_pilot(False, 0.0, 1.0))["decisive"] \
        == pytest.approx(1 - 0.225)
    # the design effect of shared donors (review 1, V1): S1 as a whole at the same thresholds
    sens = {d: oc.joint_error_rule_deff(790, 0.025, 5, d) for d in (1.0, 1.1, 1.2, 1.5)}
    assert sens[1.0][0] == pytest.approx(0.915, abs=1e-3) and sens[1.1][0] == pytest.approx(0.892, abs=1e-3)
    assert sens[1.5] == pytest.approx((0.800, 0.085), abs=1e-3)


def oc_s2_group(cond):
    import oc
    return oc.s2_group(cond, cond.variants[0][0])


def _rows(stratum, n, decisive, **kw):
    """Hand-made cards for oc.criteria_rules: establishable, of one stratum, with a sound validator's
    decisiveness and no error."""
    cond = P.conditions()["N6a" if stratum == "invalid" else "N1" if stratum == "null" else "E1"]
    row = dict(condition=cond.name, key=False, sup_error_possible=stratum != "effect",
               invalid_error_possible=stratum != "invalid", nde_error_possible=stratum != "null",
               opp_error_possible=stratum != "null", depth_error_possible=True, degenerate_allowed=False,
               valid=stratum != "invalid", stratum=stratum, data_stratum=stratum,
               s2_group=oc_s2_group(cond), establishable=True, p_sup=0.0, p_inv=0.0, p_nde=0.0, p_opp=0.0,
               p_depth=0.0, p_err=0.0, decisive=decisive, measured=True, r_sup=0.025,
               r_inv=0.1 if stratum != "invalid" else 0.0, r_nde=0.05, r_opp=0.025, r_depth=0.05,
               opp_null_counted=stratum == "null", p_oppnull=0.0, r_oppnull=0.025 if stratum == "null" else 0.0)
    return [dict(row, **kw) for _ in range(n)]


def test_s3_judges_a_stratum_whose_measured_decisiveness_is_one():
    """The dry run's smoke pilot measured the blind and useless cards' decisiveness at 1.000: doubling
    a shortfall of 0 leaves nothing to tell apart, so that stratum went unjudged and a validator that
    never shows blindness would have passed S3. The nominal is the measured rate at most 0.85."""
    pytest.importorskip("scipy")
    import oc
    rows = _rows("invalid", 400, 1.0) + _rows("null", 400, 0.85) + _rows("effect", 400, 0.85)
    rules = oc.criteria_rules(rows, sims=2000)
    r3 = rules["S3"]["invalid"]
    assert r3["measured"] == 1.0 and r3["nominal"] == oc.D_NOMINAL and r3["tier"] == "principle"
    assert 0 < r3["min_required"] <= 0.85 * 400 and rules["joint"] >= 0.90 and rules["s3_feasible"]
    assert oc.error_rule(400, 0.0)[0] == 0  # where the sound validator cannot err, no error is allowed


def test_s3_judges_every_stratum_at_the_strictest_tier_the_joint_requirement_allows():
    """81 establishable real-effect cards at d = 0.831 pass a sound validator at the principle's tier
    with 0.973, which with S1 as a whole (0.915) would bring P(S1-S7 together | sound) below 0.90
    (decided in the third round): that stratum is judged at the floor's tier (a validator half as
    decisive fails), the nulls and the blind cards at the principle's. Before the second review such
    a stratum was reported only. In oc.log's scenarios every stratum now takes the principle's tier:
    scenario C's 81 real-effect cards, at the d = 0.786 their designed errors leave since the third
    review, pass a sound validator there with 0.991."""
    pytest.importorskip("scipy")
    import oc
    key = [r for c in P.CONDITIONS if c.key  # S1's cards, so that the joint holds S1 at its design
           for r in _rows("null", P.KEY_N, 0.85, key=True, condition=c.name, establishable=False)]
    rows = key + _rows("invalid", 400, 0.85) + _rows("null", 400, 0.85) + _rows("effect", 81, 0.831)
    rules = oc.criteria_rules(rows, sims=4000)
    eff = rules["S3"]["effect"]
    assert eff["n"] == 81 and eff["tiers"]["principle"]["attainable"] and eff["tier"] == "floor"
    assert eff["tiers"]["principle"]["p_pass_sound"] < 0.90 / 0.915 < eff["tiers"]["floor"]["p_pass_sound"]
    assert eff["min_required"] == eff["tiers"]["floor"]["min_required"] and eff["tiers"]["floor"]["alternative"] \
        == pytest.approx(eff["nominal"] / 2)
    assert rules["S3"]["null"]["tier"] == rules["S3"]["invalid"]["tier"] == "principle" and rules["joint"] >= 0.90
    rules = oc.criteria_rules(oc.expected_rows(oc.scenario_pilot(True, 1.0, 0.25)), sims=4000)  # scenario C
    eff = rules["S3"]["effect"]
    assert eff["n"] == 81 and eff["tier"] == "principle" and eff["tiers"]["principle"]["p_pass_sound"] > 0.99
    assert all(r["tier"] == "principle" for r in rules["S3"].values()) and rules["joint"] >= 0.90
    rules = oc.criteria_rules(oc.expected_rows(oc.scenario_pilot(True, 0.5, 0.5)), sims=4000)    # scenario A
    assert all(r["tier"] == "principle" for r in rules["S3"].values())


def test_s3_fails_when_a_stratum_cannot_be_judged():
    """Every stratum is judged (the second review's K1): with 5 establishable real-effect cards not even
    the floor's tier is attainable, so S3 fails for any validator, a perfect one included, and the
    joint probability of a sound validator is 0 — visible in oc_pilot.log before the key."""
    pytest.importorskip("scipy")
    import oc
    rows = _rows("invalid", 400, 0.9) + _rows("null", 400, 0.9) + _rows("effect", 5, 0.9)
    rules = oc.criteria_rules(rows, sims=2000)
    assert not rules["S3"]["effect"]["tiers"]["floor"]["attainable"] and rules["S3"]["effect"]["tier"] is None
    assert not rules["s3_feasible"] and rules["joint"] == 0.0
    assert rules["S3"]["null"]["tier"] is not None  # the strata that fit are still judged, and reported


def test_error_criteria_and_decisiveness_meet_the_principle_at_any_rate():
    """The second review's V4 and V5: a mean nominal of 0.001 on 5,600 cards passed a sound validator
    with 0.51 (too few expected errors to tell a doubled rate apart); the nominal is raised to the
    first rate at which the principle holds. Decisiveness below 0.5 crashed the rule (a negative
    doubled rate); the alternative is at least half the nominal."""
    pytest.importorskip("scipy")
    import oc
    assert oc.error_rule(5600, 0.001)[1] < 0.90
    r = oc.error_criterion(5600, 0.001)
    assert r["raised"] and r["judged"] and r["p_pass_sound"] >= 0.90 and r["p_pass_doubled"] <= 0.05
    assert 0.001 < r["rule_nominal"] < 0.004
    assert not oc.error_criterion(5600, 0.025)["raised"]
    assert oc.error_criterion(0, 0.1)["judged"] is False                     # no card: holds
    for d0 in (0.2, 0.45, 0.5, 0.51, 0.6):
        k, ps, pa = oc.decisiveness_rule(400, d0)
        assert k is not None and pa <= 0.05 and oc.tier_alternative(d0, "principle") >= d0 / 2
    assert oc.decisiveness_rule(0, 0.8)[0] is None and oc.decisiveness_rule(100, float("nan"))[0] is None
    assert oc.tier_alternative(0.85, "principle") == pytest.approx(0.70)
    assert oc.tier_alternative(0.85, "floor") == pytest.approx(0.425)


def test_the_tiers_fixed_before_the_key_are_the_ones_scored(bgs, tmp_path):
    """The second review's V3: score.py recomputed which strata were judged from the realized cards,
    so the judged set was not known before the key. `oc.py --pilot ... --write-judged` writes S3's
    tiers into pilot.json; score.py judges at them, with thresholds from the realized n."""
    pytest.importorskip("scipy")
    import oc
    import score as S
    pilot = oc.scenario_pilot(True, 1.0, 0.25)
    path = tmp_path / "pilot.json"
    path.write_text(json.dumps(pilot))
    oc.main(["--pilot", str(path), "--write-judged"])
    fixed = json.loads(path.read_text())["s3_rules"]
    # scenario C: since the third review its 81 real-effect cards are judged at the principle's tier too
    assert fixed["tiers"] == {"effect": "principle", "null": "principle", "invalid": "principle"} and fixed["feasible"]
    entries = P.assign(KEY, pool_sizes=pilot["pool_size"])
    pilot["s3_rules"] = dict(fixed, tiers={"effect": "floor", "null": "floor", "invalid": "principle"})
    res = S.score(entries, _reports(entries, _perfect(pilot)), pilot, sims=2000)["criteria"]["S3"]
    assert res["tiers_fixed_before_the_key"]
    assert {st: v["tier"] for st, v in res["strata"].items()} == pilot["s3_rules"]["tiers"]
    eff = res["strata"]["effect"]
    assert eff["min_required"] == oc.decisiveness_tiers(eff["n"], eff["nominal"])["floor"]["min_required"]


def test_the_shared_donor_interval_widens_only_when_donors_drive_the_outcome():
    pytest.importorskip("scipy")
    import score as S
    rng = np.random.default_rng(0)
    donors = [f"d{i}" for i in range(40)]
    n, cover_naive, cover_adj, reps = 300, 0, 0, 60
    deffs_indep = []
    for r in range(reps):
        sets = [list(rng.choice(donors, 8, replace=False)) for _ in range(n)]
        y = rng.random(n) < 0.1  # no donor effect
        deffs_indep.append(S.overlap_interval(y, sets)["deff"])
        u = {d: rng.normal(0, 1.2) for d in donors}  # strong donor effects
        lin = np.array([np.mean([u[d] for d in s]) for s in sets])
        truth = 0.1
        y = rng.random(n) < 1 / (1 + np.exp(-(np.log(truth / (1 - truth)) + 1.5 * lin)))
        lo, hi = S.cp(int(y.sum()), n)
        ov = S.overlap_interval(y, sets)
        p_true = np.mean(1 / (1 + np.exp(-(np.log(truth / (1 - truth)) + 1.5 * rng.normal(0, 1.2 / np.sqrt(8), 20000)))))
        cover_naive += lo <= p_true <= hi
        cover_adj += ov["ci95"][0] <= p_true <= ov["ci95"][1]
    assert np.median(deffs_indep) < 1.3
    assert cover_adj > cover_naive and cover_adj / reps >= 0.85


def test_a_design_effect_above_1_5_is_reported_as_a_limitation(bgs):
    pytest.importorskip("scipy")
    import score as S
    pilot = _est_pilot(bgs)
    entries = [e for e in P.assign(KEY, pool_sizes=_sizes(bgs)) if e["condition"] in ("N1", "N8")]
    rng = np.random.default_rng(3)
    donors = {e["id"]: [f"d{j}" for j in rng.choice(30, 8, replace=False)] for e in entries}
    bad = {f"d{j}" for j in range(6)}  # the datasets with these donors err: a donor-driven outcome
    reports = _reports(entries, lambda c, e: P.SUPPORTED if bad & set(donors[e["id"]]) else P.INCONCLUSIVE)
    res = S.score(entries, reports, pilot, donors)
    assert res["secondary"]["shared_donors"]["S1:N1"]["deff"] > 1.5
    assert any(lim.startswith("S1:N1") for lim in res["limitations"])


# --------------------------------------------------------------------------- #
# the runner and the blind run
# --------------------------------------------------------------------------- #
def test_runner_runs_the_engine_on_the_fly(bgs, tmp_path):
    import run_panel as R
    pilot = _pilot(bgs, dose=2.0)
    entries = P.assign(KEY, pool_sizes=_sizes(bgs))
    pick = [next(e for e in entries if e["condition"] == c) for c in ("N1", "E1", "N6b")]
    summary = R.run(pick, bgs, pilot, tmp_path / "out", workers=1)
    assert summary["run"] == 3 and summary["errors"] == 0 and not summary["logged_twice"]
    assert summary["machine"]["cpus"] >= 1
    manifest = json.loads((tmp_path / "out" / "manifest.json").read_text())
    for row, e in zip(manifest["datasets"], sorted(pick, key=lambda x: x["id"])):
        X, obs, genes, cards = P.build(e, bgs, pilot)
        assert row["data_sha256"] == P.dataset_sha256(X, obs, genes)  # rebuilt from the key
        assert row["donors"] == sorted(set(obs["donor"]))
        rep = json.loads((tmp_path / "out" / "reports" / f"{e['id']}.json").read_text())
        assert rep["verdict"] and rep["cause"] and "elapsed_seconds" not in rep
        assert not set(R.RUNTIME_FIELDS) & set(rep["provenance"])
    runtime = json.loads((tmp_path / "out" / "runtime.json").read_text())
    assert all(c["seconds"] > 0 and c["timestamp_utc"] for c in runtime["cards"])
    again = R.run(pick, bgs, pilot, tmp_path / "out", workers=1)
    assert again["skipped"] == 3  # one attempt per card
    assert len((tmp_path / "out" / "runlog.jsonl").read_text().splitlines()) == 3
    if pytest.importorskip("scipy"):
        import score as S
        reports, donors = S.read_results(tmp_path / "out")
        assert all(len(r) == 3 and r[2]["gate4"] in ("PASS", "FAIL", "UNTESTED", "SKIP", None)
                   for r in reports.values() if r)
        res = S.score(pick, reports, pilot, donors)
        assert res["n_cards"] == 3 and "gate4" in res["per_gate"] and "S5" in res["criteria"]


DETERMINISM = r"""
import json, sys
from pathlib import Path
sys.path.insert(0, {here!r})
import panel as P, run_panel as R, simulate
bgs = {{"B1": simulate.simulated_background("B1", n_genes=240), "B2": simulate.simulated_background("B2", donors=12, n_genes=240, seed=1)}}
pilot = simulate.dry_pilot(bgs)
entries = P.assign({key!r}, pool_sizes={{k: len(b.plan["pool"]) for k, b in bgs.items()}})
conds = P.conditions()
pick = [next(e for e in entries if e["condition"] == c.name) for c in P.CONDITIONS]  # every condition
cards = sum(conds[e["condition"]].cards for e in pick)
for e in entries:  # then the key's order, up to 20 claim cards
    if cards >= 20:
        break
    if e not in pick and conds[e["condition"]].cards == 1:
        pick.append(e)
        cards += 1
R.run(pick, bgs, pilot, Path({out!r}), workers={workers})
"""


def test_the_same_cards_give_byte_identical_reports(tmp_path):
    """The run is a deterministic function of the key, the backgrounds and the pilot: the same 20
    claim cards, run twice in fresh processes with different hash seeds and numbers of workers,
    give byte-identical reports and manifests (the time and the machine go to runtime.json)."""
    pytest.importorskip("scipy")
    outs = []
    for run, (seed, workers) in enumerate((("1", 1), ("2", 2))):
        out = tmp_path / f"run{run}"
        code = DETERMINISM.format(here=str(HERE), key=KEY, out=str(out), workers=workers)
        env = dict(os.environ, PYTHONHASHSEED=seed)
        subprocess.run([sys.executable, "-c", code], cwd=HERE, env=env, check=True, capture_output=True)
        outs.append(out)
    reports = [sorted((p.name, p.read_bytes()) for p in (o / "reports").glob("*.json")) for o in outs]
    assert len(reports[0]) == 20 and reports[0] == reports[1]
    assert (outs[0] / "manifest.json").read_bytes() == (outs[1] / "manifest.json").read_bytes()
    conds = {json.loads(b)["id"] for _, b in reports[0]}
    assert len(conds) == 20


def test_runner_pins_one_blas_thread_per_worker_before_numpy_loads():
    """The thread and kernel variables count only if they are set before numpy loads its BLAS; a
    forked worker inherits the pool numpy was loaded with. Every script that runs cards sets the
    same numerical environment (frozen.NUMERIC_ENV) unless its caller did."""
    import frozen
    env = {k: v for k, v in os.environ.items() if k not in frozen.NUMERIC_ENV}
    for script in ("run_panel", "blind", "timing"):
        code = (f"import json, {script}, run_panel; "
                "print(json.dumps(run_panel.machine()))")
        out = subprocess.run([sys.executable, "-c", code], cwd=HERE, env=env,
                             capture_output=True, text=True, check=True).stdout
        m = json.loads(out)
        assert m["numeric_env"] == frozen.NUMERIC_ENV, script
        if m["numpy_simd"] is not None:  # and numpy dispatches as the environment says (the third review)
            assert not {"X86_V4", "AVX512_ICL", "AVX512_SPR"} & set(m["numpy_simd"]["dispatch"]), script
        if m["blas_threads_in_use"] is not None:  # read with threadpoolctl where it is installed
            assert m["blas_threads_in_use"] == [1]
    env["OMP_NUM_THREADS"] = "2"  # a value the caller set wins
    out = subprocess.run([sys.executable, "-c", "import json, run_panel; "
                          "print(json.dumps(run_panel.machine()))"], cwd=HERE, env=env,
                         capture_output=True, text=True, check=True).stdout
    assert json.loads(out)["numeric_env"]["OMP_NUM_THREADS"] == "2"


def test_shards_partition_the_key_order():
    import blind as B
    entries = P.assign(KEY)
    parts = [B.shard_entries(entries, k, 20) for k in range(20)]
    ids = [e["id"] for p in parts for e in p]
    assert sorted(ids) == sorted(e["id"] for e in entries) and len(ids) == len(set(ids))
    assert max(len(p) for p in parts) - min(len(p) for p in parts) <= 1
    with pytest.raises(SystemExit):
        B.shard_entries(entries, 20, 20)


def _fake_shard(root: Path, k: int, shards: int, key: str, ids: list):
    d = root / f"shard-{k}"
    (d / "reports").mkdir(parents=True)
    rows = []
    for i in ids:
        text = json.dumps(dict(id=i, verdict="INCONCLUSIVE — x", cause="effect_inconclusive"))
        (d / "reports" / f"{i}.json").write_text(text)
        rows.append(dict(id=i, data_sha256="0" * 64, donors=["a"], cards=[
            dict(id=i, card_sha256="1" * 64, report_sha256=hashlib.sha256(text.encode()).hexdigest())]))
    (d / "manifest.json").write_text(json.dumps(dict(datasets=rows)))
    (d / "runtime.json").write_text(json.dumps(dict(summary=dict(wall_seconds=1.0, workers=4), cards=[])))
    (d / "runlog.jsonl").write_text("".join(json.dumps(dict(claim_id=i)) + "\n" for i in ids))
    meta = dict(shard=k, shards=shards, key=key, beacon=dict(chain="quicknet", round=7), pilot_sha256="p",
                backgrounds={"B1": {}}, entries=len(ids))
    (d / "shard.json").write_text(json.dumps(meta))


def test_collect_merges_shards_and_refuses_inconsistent_ones(tmp_path):
    import blind as B
    _fake_shard(tmp_path / "ok", 0, 2, KEY, ["D1", "D3"])
    _fake_shard(tmp_path / "ok", 1, 2, KEY, ["D2"])
    s = B.collect(tmp_path / "ok", tmp_path / "res", expected_datasets=3, reruns=[dict(shard=1, reason="runner lost")])
    assert s["datasets"] == 3 and s["reruns"] == 1 and (tmp_path / "res" / "SHA256SUMS").exists()
    manifest = json.loads((tmp_path / "res" / "manifest.json").read_text())
    assert [r["id"] for r in manifest["datasets"]] == ["D1", "D2", "D3"] and manifest["key"] == KEY
    with pytest.raises(SystemExit):
        B.collect(tmp_path / "ok", tmp_path / "res2", expected_datasets=4)
    _fake_shard(tmp_path / "keys", 0, 2, KEY, ["D1"])
    _fake_shard(tmp_path / "keys", 1, 2, OTHER_KEY, ["D2"])
    with pytest.raises(SystemExit, match="key"):
        B.collect(tmp_path / "keys", tmp_path / "res3")
    _fake_shard(tmp_path / "missing", 0, 3, KEY, ["D1"])
    with pytest.raises(SystemExit, match="shards present"):
        B.collect(tmp_path / "missing", tmp_path / "res4")
    _fake_shard(tmp_path / "twice", 0, 2, KEY, ["D1"])
    _fake_shard(tmp_path / "twice", 1, 2, KEY, ["D1"])
    with pytest.raises(SystemExit, match="twice"):
        B.collect(tmp_path / "twice", tmp_path / "res5")


def test_the_manifest_keeps_the_round_not_when_or_where_it_was_fetched(tmp_path):
    """The manifest is a function of the code, the backgrounds, the pilot and the round: the key
    record's fetch time and relay answers stay in key.json."""
    import blind as B
    rec = dict(chain="quicknet", round=7, signature="ab" * 48, fetched_utc="2026-10-08T17:35:10+00:00",
               answers=[dict(relay="https://api.drand.sh", verified=True)])
    rec["randomness"] = rec["key"] = hashlib.sha256(bytes.fromhex(rec["signature"])).hexdigest()
    (tmp_path / "key.json").write_text(json.dumps(rec))
    key, beacon = B.read_key(True, tmp_path / "key.json", None)
    assert key == rec["key"] and beacon["round"] == 7 and beacon["randomness"] == key
    assert "fetched_utc" not in beacon and "answers" not in beacon


def test_verify_re_runs_a_finished_run_and_compares_its_reports(tmp_path):
    """`blind.py verify`, the check anyone can run: the first datasets of a finished run, re-run from
    its key, give the published datasets' and reports' sha256; a manifest whose hash was altered, or
    another pilot, does not pass. The second review's V1: between CPU models the last digits of a
    report's numbers can differ, so a report that says the same (verdict, cause, every field; numbers
    within 1e-9) is the same, not different; one with another verdict is different; a dataset whose
    hash differs is reported; and verify refuses a non-empty output (the runner reuses reports there)."""
    pytest.importorskip("scipy")
    import blind as B
    import simulate
    bgs = {"B1": simulate.simulated_background("B1", n_genes=240),
           "B2": simulate.simulated_background("B2", donors=12, n_genes=240, seed=1)}
    pilot = simulate.dry_pilot(bgs)
    rec = dict(chain="quicknet", round=7, signature="ab" * 48)
    rec["randomness"] = rec["key"] = hashlib.sha256(bytes.fromhex(rec["signature"])).hexdigest()
    compact = tmp_path / "compact"
    B.write_compact(bgs, compact)
    B.run_shard(compact, pilot, rec["key"], dict(rec), 0, 1, tmp_path / "shards" / "shard-0", workers=1, limit=3)
    B.collect(tmp_path / "shards", tmp_path / "res")
    (tmp_path / "res" / "key.json").write_text(json.dumps(rec))
    res = B.verify(tmp_path / "res", compact, pilot, tmp_path / "rerun", workers=1, datasets=2)
    assert res["datasets"] == 2 and res["cards"] >= 2 and res["identical"] == res["cards"] and res["key"] == rec["key"]
    assert res["datasets_identical"] == 2 and res["different"] == 0
    with pytest.raises(SystemExit, match="not empty"):
        B.verify(tmp_path / "res", compact, pilot, tmp_path / "rerun", workers=1, datasets=1)
    m = json.loads((tmp_path / "res" / "manifest.json").read_text())
    first, second = m["datasets"][0]["cards"][0], m["datasets"][1]["cards"][0]
    m["datasets"][0]["cards"][0]["report_sha256"] = "0" * 64        # a hash that binds no report
    path = tmp_path / "res" / "reports" / f"{second['id']}.json"   # as if run on another CPU model
    rep = json.loads(path.read_text())
    nums = [k for k, v in rep["fields"]["effect"]["detail"].items() if isinstance(v, float) and v]
    assert nums
    rep["fields"]["effect"]["detail"][nums[0]] *= 1 + 8e-13
    path.write_text(json.dumps(rep, sort_keys=True))
    m["datasets"][1]["cards"][0]["report_sha256"] = P.sha256(path)
    m["datasets"][2]["data_sha256"] = "1" * 64
    (tmp_path / "res" / "manifest.json").write_text(json.dumps(m))
    # the third review: verify checks SHA256SUMS first, so files changed after it are refused...
    with pytest.raises(SystemExit, match="differs from SHA256SUMS"):
        B.verify(tmp_path / "res", compact, pilot, tmp_path / "rerun-sums", workers=1)
    B.write_sums(tmp_path / "res")  # ...and the rest of this test is a run published as it stands
    res = B.verify(tmp_path / "res", compact, pilot, tmp_path / "rerun2", workers=1)
    assert res["files_checked"] >= res["cards"] + 1
    row = {r["card"]: r for r in res["rows"]}
    assert res["datasets"] == 3 and res["identical"] == res["cards"] - 2 and res["datasets_identical"] == 2
    assert not row[first["id"]]["same"] and row[second["id"]]["same"] and not row[second["id"]]["identical"]
    assert res["different"] == 1
    rep["verdict"] = "INCONCLUSIVE — x"                             # another verdict: different
    path.write_text(json.dumps(rep, sort_keys=True))
    m["datasets"][1]["cards"][0]["report_sha256"] = P.sha256(path)
    (tmp_path / "res" / "manifest.json").write_text(json.dumps(m))
    B.write_sums(tmp_path / "res")
    assert not {r["card"]: r for r in B.verify(tmp_path / "res", compact, pilot, tmp_path / "rerun3",
                                               workers=1)["rows"]}[second["id"]]["same"]
    with pytest.raises(SystemExit, match="pilot"):
        B.verify(tmp_path / "res", compact, dict(pilot, dropped=["N3 steps"]), tmp_path / "rerun4", workers=1)
    assert B.same_report({"a": [1.0, "x", True]}, {"a": [1.0 + 1e-12, "x", True]})
    assert not B.same_report({"a": [1.0, "x", True]}, {"a": [1.0 + 1e-6, "x", True]})
    assert not B.same_report({"a": True}, {"a": 1}) and not B.same_report({"a": 1}, {"b": 1})


def test_the_guard_and_the_key_follow_the_frozen_tag_and_the_run_tag(tmp_path, monkeypatch):
    """In a scratch repository: the guard passes when only data files changed after the frozen tag
    and fails when frozen code changed; the key record must be the randomness (sha256 of the
    signature) of the round the run tag's message names."""
    import blind as B
    repo = tmp_path / "repo"
    (repo / "validation" / "prereg").mkdir(parents=True)
    (repo / "src").mkdir()

    def git(*a):
        subprocess.run(["git", *a], cwd=repo, check=True, capture_output=True,
                       env={"GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t", "GIT_COMMITTER_NAME": "t",
                            "GIT_COMMITTER_EMAIL": "t@t", "PATH": "/usr/bin:/bin"})
    git("init", "-q")
    (repo / "src" / "engine.py").write_text("x = 1\n")
    (repo / "validation" / "prereg" / "panel.py").write_text("y = 1\n")
    (repo / "validation" / "prereg" / "v1.md").write_text("the protocol\n")
    git("add", "-A")
    git("commit", "-qm", "frozen")
    git("tag", "-a", "v0.3.0-prereg", "-m", "frozen")
    import frozen as F
    monkeypatch.setattr(B, "HERE", repo)
    monkeypatch.setattr(F, "HERE", repo)
    (repo / "validation" / "prereg" / "pilot.json").write_text("{}")
    git("add", "-A")
    git("commit", "-qm", "data only")
    git("tag", "-a", "panel-v1-run", "-m", "run\n\ndrand quicknet round: 1234567\n")
    assert B.guard("v0.3.0-prereg")["frozen_tag"] == "v0.3.0-prereg"
    assert B.round_from_tag("panel-v1-run") == ("quicknet", 1234567)
    sig = "ab" * 48
    rec = dict(chain="quicknet", round=1234567, signature=sig, randomness=hashlib.sha256(bytes.fromhex(sig)).hexdigest())
    rec["key"] = rec["randomness"]
    (tmp_path / "key.json").write_text(json.dumps(rec))
    key, beacon = B.read_key(False, tmp_path / "key.json", "panel-v1-run")
    assert key == rec["randomness"] and beacon["round"] == 1234567
    (tmp_path / "late.json").write_text(json.dumps(dict(rec, round=1234568)))
    with pytest.raises(SystemExit, match="names"):
        B.read_key(False, tmp_path / "late.json", "panel-v1-run")
    (tmp_path / "forged.json").write_text(json.dumps(dict(rec, key="0" * 64, randomness="0" * 64)))
    with pytest.raises(SystemExit, match="sha256 of its signature"):
        B.read_key(False, tmp_path / "forged.json", "panel-v1-run")
    with pytest.raises(SystemExit, match="key record"):
        B.read_key(False, None, "panel-v1-run")
    (repo / "src" / "engine.py").write_text("x = 2\n")
    git("add", "-A")
    git("commit", "-qm", "engine change")
    with pytest.raises(SystemExit, match="frozen files differ"):
        B.guard("v0.3.0-prereg")
    git("reset", "-q", "--hard", "HEAD~1")  # back to the data-only commit: the protocol is frozen too
    (repo / "validation" / "prereg" / "v1.md").write_text("the protocol, edited after the tag\n")
    git("add", "-A")
    git("commit", "-qm", "protocol change")
    with pytest.raises(SystemExit, match="frozen files differ"):
        B.guard("v0.3.0-prereg")
    # the third review: the workflow's commits add only its records (v1.md section 2's last row)
    git("reset", "-q", "--hard", "HEAD~1")
    (repo / "validation" / "prereg" / "notes on the run.md").write_text("a file the workflow never writes\n")
    git("add", "-A")
    git("commit", "-qm", "another file")
    with pytest.raises(SystemExit, match="other than the workflow's records.*notes on the run.md"):
        B.guard("v0.3.0-prereg")
    # and the real run is one attempt: a job of a re-run workflow run stops
    assert F.one_attempt({}) == 1 and F.one_attempt({"GITHUB_RUN_ATTEMPT": "1"}) == 1
    with pytest.raises(SystemExit, match="attempt 2 of workflow run 99"):
        F.one_attempt({"GITHUB_RUN_ATTEMPT": "2", "GITHUB_RUN_ID": "99"})


# --------------------------------------------------------------------------- #
# the beacon
# --------------------------------------------------------------------------- #
def _beacon_keys():
    """A BLS key pair of quicknet's scheme (signatures on G1, public key on G2) and a signer."""
    pytest.importorskip("py_ecc")
    from py_ecc.bls.hash_to_curve import hash_to_G1
    from py_ecc.bls.point_compression import compress_G1, compress_G2
    from py_ecc.optimized_bls12_381 import G2, multiply
    import beacon as BC
    sk = 987654321987654321
    pk = compress_G2(multiply(G2, sk))
    pk_hex = pk[0].to_bytes(48, "big").hex() + pk[1].to_bytes(48, "big").hex()

    def sign(rnd):
        return compress_G1(multiply(hash_to_G1(BC.message(rnd), BC.DST_G1, hashlib.sha256), sk)).to_bytes(48, "big").hex()
    return pk_hex, sign


def test_the_beacon_verifies_a_round_and_derives_its_randomness():
    import beacon as BC
    pk, sign = _beacon_keys()
    sig = sign(4242)
    assert BC.verify(4242, sig, pk) and not BC.verify(4243, sig, pk)
    assert not BC.verify(4242, "c0" + "00" * 47, pk) and not BC.verify(4242, sig[:-2], pk)
    assert BC.randomness_of(sig) == hashlib.sha256(bytes.fromhex(sig)).hexdigest()
    assert len(BC.randomness_of(sig)) == 64 and P.check_key(BC.randomness_of(sig))


def test_the_chain_constants_verify_a_real_quicknet_round():
    """Round 32892512 of drand quicknet (17:35:00 UTC on 2026-10-08), as four relays served it to
    the dry run of workflow run 37817228721: the constants, the message, the hash to G1 and the
    pairing check reproduce its verification and its randomness offline."""
    pytest.importorskip("py_ecc")
    import beacon as BC
    sig = ("b6c102fe1446e998a199aea2167cc4351b8d9dfcc14c06b4889ceac55493359d"
           "86dcc5d7a8b59cbf4147b5da3243a735")
    assert BC.round_time(32892512) == 1791480900  # 2026-10-08T17:35:00Z
    assert BC.verify(32892512, sig) and not BC.verify(32892511, sig)
    assert BC.randomness_of(sig) == "680404b04626fd169ccfa9050543611b4a6308352d9be06ea9278dd521106f03"


def test_the_beacon_round_lies_at_least_an_hour_after_the_tag():
    import beacon as BC
    g = BC.CHAIN["genesis_time"]
    assert BC.round_time(1) == g and BC.first_round_at(g) == 1 and BC.first_round_at(g + 1) == 2
    for t in (g + 10_000.5, g + 123_456_789.0):
        r = BC.round_for_tag(t)
        assert BC.round_time(r) >= t + BC.MIN_DELAY > BC.round_time(r - 1)


def test_the_beacon_accepts_only_verified_and_agreeing_answers(monkeypatch):
    import beacon as BC
    pk, sign = _beacon_keys()
    chain = dict(BC.CHAIN, public_key=pk)
    good = dict(round=77, signature=sign(77), randomness=BC.randomness_of(sign(77)))
    answers = {"https://a": good, "https://b": dict(good, signature=sign(78)), "https://c": None}

    def fake_get(url, timeout=20.0):
        relay = url.split("/" + chain["hash"])[0]
        if answers[relay] is None:
            raise OSError("down")
        return answers[relay]
    monkeypatch.setattr(BC, "_get", fake_get)
    rec = BC.fetch(77, relays=tuple(answers), chain=chain)
    assert rec["key"] == good["randomness"] and rec["round"] == 77
    assert [a.get("verified") for a in rec["answers"]] == [True, False, None]
    answers["https://a"] = dict(good, signature=sign(79))
    with pytest.raises(RuntimeError, match="no relay"):
        BC.fetch(77, relays=tuple(answers), chain=chain)


def test_the_gap_to_the_round_is_checked_against_the_beacons_newest_round(monkeypatch):
    """The fourth review's K10: the hour between the run tag's push and its round was measured by the
    runner's clock alone; the newest verified round the relays serve after the push bounds it too,
    and the smaller of the two counts."""
    import beacon as BC
    pk, sign = _beacon_keys()
    chain = dict(BC.CHAIN, public_key=pk)
    newest = {"https://a": dict(round=1000, signature=sign(1000)),
              "https://b": dict(round=1300, signature=sign(1)),  # does not verify as round 1300
              "https://c": None}

    def fake_get(url, timeout=20.0):
        relay = url.split("/" + chain["hash"])[0]
        if newest[relay] is None:
            raise OSError("down")
        return newest[relay]
    monkeypatch.setattr(BC, "_get", fake_get)
    assert BC.latest(relays=tuple(newest), chain=chain) == 1000
    pushed = BC.round_time(1000, chain)
    assert BC.gap_after_push(2200, pushed, 1000, chain) == dict(round=2200, newest=1000, by_clock=3600,
                                                               by_beacon=3600, seconds=3600)
    slow = BC.gap_after_push(2200, pushed - 600, 1000, chain)  # the runner's clock ten minutes behind
    assert slow["by_clock"] == 4200 and slow["seconds"] == 3600
    assert BC.gap_after_push(2200, pushed, 1300, chain)["seconds"] == 2700 < BC.MIN_DELAY
    newest["https://a"] = None
    with pytest.raises(RuntimeError, match="no relay"):
        BC.latest(relays=tuple(newest), chain=chain)


# --------------------------------------------------------------------------- #
# the selection of the backgrounds and the anchors
# --------------------------------------------------------------------------- #
def _census_like(rng, groups):
    """An obs table as the Census returns it: (dataset, cell type, donor, cells, sex) blocks."""
    rows = []
    for ds, ct, donor, n, sex in groups:
        rows += [dict(dataset_id=ds, cell_type=ct, tissue="lung", donor_id=donor, sex=sex,
                      development_stage="3-month-old stage") for _ in range(n)]
    obs = pd.DataFrame(rows)
    obs["soma_joinid"] = rng.permutation(len(obs))
    return obs


def test_the_selection_rule_lists_the_candidates_and_picks_by_the_rule():
    import select_backgrounds as SB
    rng = np.random.default_rng(0)
    groups = ([("d1", "fibroblast", f"a{i}", 250, "female") for i in range(30)]          # 30 donors
              + [("d1", "fibroblast", "small", 150, "male")]                            # below 200 cells
              + [("d2", "T cell", f"b{i}", 300, "male") for i in range(30)]               # 30, more cells each
              + [("d2", "T cell", "unknown", 900, "male")]                                # no donor ID
              + [("d3", "B cell", f"c{i}", 400, "male") for i in range(10)])               # too few donors
    obs = _census_like(rng, groups)
    cand = SB.candidates(SB.count_cells([obs], ["dataset_id", "cell_type", "donor_id"]), ["dataset_id", "cell_type"],
                         SB.B1_MIN_DONORS)
    assert list(cand["dataset_id"]) == ["d2", "d1"]  # tie on 30 donors: more cells per donor wins
    assert list(cand["donors"]) == [30, 30] and cand["cells_per_donor"].iloc[0] == 300
    ids = SB.sample_cells(obs, cand.iloc[0].to_dict(), ["dataset_id", "cell_type"], max_donors=20, max_cells=100)
    picked = obs.set_index("soma_joinid").loc[ids]
    assert picked["donor_id"].nunique() == 20 and (picked.groupby("donor_id").size() == 100).all()
    assert "unknown" not in set(picked["donor_id"])
    assert np.array_equal(ids, SB.sample_cells(obs, cand.iloc[0].to_dict(), ["dataset_id", "cell_type"], 20, 100))
    # B2: both sexes among the qualifying mice
    mice = ([("m1", "fibroblast", f"f{i}", 250, "female") for i in range(8)]
            + [("m1", "fibroblast", f"m{i}", 250, "male") for i in range(6)]
            + [("m2", "fibroblast", f"x{i}", 250, "male") for i in range(20)])  # one sex only
    keys2 = ["dataset_id", "tissue", "cell_type"]
    cand2 = SB.candidates(SB.count_cells([_census_like(rng, mice)], keys2 + ["donor_id", "sex"]), keys2, SB.B2_MIN_MICE,
                          both_sexes=True)
    assert list(cand2["dataset_id"]) == ["m1"] and cand2["donors"].iloc[0] == 14


def test_the_b3_tables_are_read_by_their_cell_columns():
    """E-MTAB-2805's tables: gene ID, transcript ID, symbol and gene length (missing on the ERCC rows)
    before the cells; only the columns named as cells are counts, so a missing gene length does not
    make the counts look non-raw."""
    import select_backgrounds as SB
    hdr = "\t".join(["EnsemblGeneID", "EnsemblTranscriptID", "AssociatedGeneName", "GeneLength"]
                    + [f"G1_cell{i}_count" for i in range(1, 4)])
    rows = ["ENSMUSG01\tENSMUST01\tActb\t1800\t5\t0\t3", "ERCC-00002\t\t\t\t10\t12\t9",
            "ENSMUSG02\tENSMUST02\tGapdh\t1300\t7\t1\t0"]
    text = "\n".join([hdr] + rows) + "\n"
    X, genes, obs = SB.parse_buettner({"G1": text, "G2M": text.replace("G1_", "G2M_")})
    assert X.shape == (6, 3) and genes == ["Actb", "ERCC-00002", "Gapdh"] and SB.raw_counts(X)
    assert list(obs["phase"]) == ["G1"] * 3 + ["G2M"] * 3
    assert any("missing 1" in ln for ln in SB.describe_table(text, "B3 G1"))


def test_b4_is_checked_against_its_criteria_on_the_deposits_tables():
    """GSE146773's tables: cells x Ensembl genes, the cells named <well>_<plate>, and one FUCCI row per
    cell. RSEM expected counts fail the raw-counts criterion; a cell without FUCCI intensities fails
    "for every cell"; integer counts with FUCCI for every cell and plate IDs pass."""
    import select_backgrounds as SB
    rng = np.random.default_rng(2)
    cells = [f"{w}{i}_{p}" for p in (355, 356) for w in "ABCDEFGH" for i in range(1, 25)]  # 384 cells
    counts = pd.DataFrame(rng.poisson(3.0, size=(len(cells), 5)).astype(float), index=cells,
                          columns=[f"ENSG0000000000{j}" for j in range(5)])
    fucci = pd.DataFrame({"cell": cells, "raw_green530": 10.0, "raw_red585": 20.0})
    fails, facts = SB.check_b4(counts, fucci)
    assert fails == [] and facts["cells"] == 384 and facts["plates"] == ["355", "356"]
    rsem = counts.copy()
    rsem.iloc[0, 0] = 25.53
    fails, facts = SB.check_b4(rsem, fucci.iloc[1:])
    assert facts["non_integer"] == 1 and len(fails) == 2
    assert fails[0].startswith("not raw counts") and fails[1] == "FUCCI intensities for 383 of 384 cells"
    fails, _ = SB.check_b4(counts.iloc[:100], fucci)
    assert fails == ["100 cells"]


def test_the_census_is_counted_in_chunks_on_its_category_codes():
    """The Census obs arrives as Arrow chunks with dictionary columns (75 million human cells): counting
    chunk by chunk on the codes gives the counts of the whole table; missing donor IDs are not donors;
    a donor's sex is the commonest among its cells; the extraction sample reads categoricals."""
    pa = pytest.importorskip("pyarrow")
    import select_backgrounds as SB
    rng = np.random.default_rng(3)
    groups = ([("d1", "T cell", f"a{i}", 210 + i, "female" if i % 2 else "male") for i in range(26)]
              + [("d1", "T cell", "unknown", 500, "male"), ("d1", "B cell", "a1", 300, "female")])
    obs = _census_like(rng, groups)
    obs.loc[obs.index[:5], "donor_id"] = None                       # missing donor IDs
    obs.loc[(obs["donor_id"] == "a3").to_numpy().nonzero()[0][:30], "sex"] = "unknown"  # a minority label
    cols = ["dataset_id", "cell_type", "donor_id", "sex"]
    whole = SB.count_cells([obs], cols)
    chunks = [pa.Table.from_pandas(obs.iloc[i:i + 1000][cols].astype("category"), preserve_index=False)
              for i in range(0, len(obs), 1000)]
    assert all(pa.types.is_dictionary(c.type) for c in chunks[0].schema)
    assert SB.count_cells(chunks, cols).equals(whole)
    assert whole["cells"].sum() == len(obs) and "nan" in set(whole["donor_id"])
    keys = ["dataset_id", "cell_type"]
    cand = SB.candidates(whole, keys, 24, both_sexes=True)
    assert len(cand) == 1 and cand["donors"].iloc[0] == 26 and cand["sexes"].iloc[0] == "female,male"
    cat = obs.astype({c: "category" for c in cols})
    ids = SB.sample_cells(cat, cand.iloc[0].to_dict(), keys, max_donors=10, max_cells=200)
    assert np.array_equal(ids, SB.sample_cells(obs, cand.iloc[0].to_dict(), keys, max_donors=10, max_cells=200))
    picked = obs.set_index("soma_joinid").loc[ids]
    assert set(picked["cell_type"]) == {"T cell"} and picked["donor_id"].notna().all()
    assert "unknown" not in set(picked["donor_id"]) and len(ids) == 2000


def _fake_census(monkeypatch, obs_by_org, raw_ok=True, no_pool=()):
    """A stand-in for cellxgene_census with the API the selection uses: obs read in Arrow chunks with
    dictionary columns (a TableReadIter with .concat()), and get_anndata returning only the obs columns
    asked for, on a string index. Its counts are simulate.py's (coupled pairs at every level, so that
    the pair pool exists), except for the datasets in `no_pool` (independent Poisson counts)."""
    import sys
    import types

    import simulate
    pa = pytest.importorskip("pyarrow")
    sparse = pytest.importorskip("scipy.sparse")

    class ReadIter:
        def __init__(self, tables):
            self.tables = tables

        def __iter__(self):
            return iter(self.tables)

        def concat(self):
            return pa.concat_tables(self.tables)

    def _filter(df, value_filter):
        if 'dataset_id == "' in value_filter:
            df = df[df["dataset_id"] == value_filter.split('dataset_id == "')[1].split('"')[0]]
        return df

    class Obs:
        def __init__(self, df):
            self.df = df

        def read(self, value_filter, column_names):
            df = _filter(self.df, value_filter)[list(column_names)]
            cat = df.astype({c: "category" for c in column_names if c != "soma_joinid"})
            return ReadIter([pa.Table.from_pandas(cat.iloc[i:i + 700], preserve_index=False)
                             for i in range(0, len(cat), 700)])

    class Census(dict):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    class Datasets:
        def read(self):
            ids = sorted({d for df in obs_by_org.values() for d in df["dataset_id"]})
            return ReadIter([pa.table({"dataset_id": ids, "dataset_title": [f"title of {d}" for d in ids]})])

    class Ad:
        def __init__(self, X, obs, var):
            self.X, self.obs, self.var = X, obs, var
            self.n_obs, self.n_vars = X.shape

        def __getitem__(self, idx):
            return Ad(self.X[idx], self.obs.iloc[idx], self.var)

    def get_anndata(census, organism, X_name, obs_coords, obs_column_names, var_column_names):
        df = obs_by_org[organism].set_index("soma_joinid", drop=False).loc[list(obs_coords)]
        obs = df[list(obs_column_names)].reset_index(drop=True)
        obs.index = obs.index.astype(str)
        donors = sorted(set(df["donor_id"]))
        per = df.groupby("donor_id").cumcount().to_numpy()
        cells = int(per.max()) + 1
        if set(df["dataset_id"]) & set(no_pool):
            X = np.random.default_rng(0).poisson(1.0, size=(len(obs), 420)).astype(np.float32)
            genes = [f"GENE{j}" for j in range(420)]
        else:
            bg = simulate.simulated_background("S", donors=len(donors), cells=cells, n_genes=420, plan=False)
            rows = [donors.index(d) * cells + j for d, j in zip(df["donor_id"], per)]
            X, genes = bg.X[rows].astype(np.float32), bg.genes
        if not raw_ok:
            X = X * 0.5
        return Ad(sparse.csr_matrix(X), obs, pd.DataFrame({"feature_name": genes}))

    census = Census(census_info={"datasets": Datasets()},
                    census_data={o.lower().replace(" ", "_"): types.SimpleNamespace(obs=Obs(df))
                                 for o, df in obs_by_org.items()})
    mod = types.SimpleNamespace(
        get_census_version_description=lambda v: {"release_build": "2099-01-01", "lts": True},
        open_soma=lambda census_version: census, get_anndata=get_anndata)
    monkeypatch.setitem(sys.modules, "cellxgene_census", mod)


def test_the_census_selection_runs_end_to_end_on_a_stand_in(monkeypatch, tmp_path):
    """select_census on a stand-in Census: the rehearsal names no identity and keeps no file; the real
    selection lists the candidates, extracts the chosen one with its sha256, and moves to the next
    rank when a candidate's counts are not raw or have no pair pool (review 1, V5)."""
    import select_backgrounds as SB
    rng = np.random.default_rng(5)
    human = _census_like(rng, [("h1", "T cell", f"a{i}", 210, "female") for i in range(25)]
                         + [("h2", "B cell", f"b{i}", 230, "male") for i in range(25)])
    mouse = _census_like(rng, [("m1", "fibroblast", f"f{i}", 205, "female") for i in range(7)]
                         + [("m1", "fibroblast", f"m{i}", 205, "male") for i in range(6)])
    for df in (human, mouse):
        df["development_stage"] = "adult"
    orgs = {"Homo sapiens": human, "Mus musculus": mouse}
    _fake_census(monkeypatch, orgs)
    spec, lines = SB.select_census(tmp_path, rehearsal=True)
    assert spec == {} and not list(tmp_path.iterdir())
    text = "\n".join(lines)
    assert "h2" not in text and "B cell" not in text and "m1" not in text and "chosen rank 1" in text
    assert "sha256" not in text and "donors" not in text  # nothing that singles the choice out (the third review)
    spec, lines = SB.select_census(tmp_path)
    assert spec["B1"]["source"]["dataset_id"] == "h2" and spec["B1"]["source"]["rank"] == 1
    assert spec["B1"]["source"]["extracted_cells"] == 25 * 230  # tie on 25 donors: more cells per donor
    assert spec["B2"]["source"]["qualifying_donors"] == 13 and spec["B2"]["source"]["rule"].startswith("all donors")
    z = np.load(tmp_path / "B1.npz")
    assert P.sha256(tmp_path / "B1.npz") == spec["B1"]["sha256"] and set(z["obs_cell_type"]) == {"B cell"}
    assert z["X_shape"][0] == spec["B1"]["source"]["extracted_cells"] and len(set(z["obs_donor_id"])) == 25
    bg = P.load_background("B1", P.resolve_background(tmp_path / "backgrounds.json", "B1", spec["B1"], tmp_path))
    P.plan_background(bg)  # the chosen background has the pool the pilot needs
    assert len(bg.plan["pool"]) >= 3 * P.MIN_PAIRS_PER_LEVEL
    _fake_census(monkeypatch, orgs, no_pool=("h2",))  # rank 1 of B1 has no pair pool: rank 2 is taken
    spec, lines = SB.select_census(tmp_path)
    assert spec["B1"]["source"]["dataset_id"] == "h1" and spec["B1"]["source"]["rank"] == 2
    assert sum("rank 1 fails the pool rule" in x for x in lines) == 1 and spec["B2"]["source"]["rank"] == 1
    _fake_census(monkeypatch, orgs, raw_ok=False)
    spec, lines = SB.select_census(tmp_path)
    assert spec == {} and sum("fails the raw-counts criterion" in x for x in lines) == 2 + 1
    assert sum("none of the first" in x for x in lines) == 2


def test_the_anchors_run_on_their_backgrounds(tmp_path):
    """R1 on a B2-like file (sex, age, mouse; Xist in females, Y genes in males) and R2 on a
    B3-like file (phases as single batches, ERCC rows): the claims run and are scored against
    their allowed sets; R3 is dropped with B4."""
    pytest.importorskip("scipy")
    import anchors as A
    rng = np.random.default_rng(1)
    genes = ["Xist", *A.Y_GENES, *[f"g{i}" for i in range(40)]]
    rows, obs = [], []
    for sex in ("female", "male"):
        for m in range(5):
            for _ in range(60):
                base = rng.poisson(3.0, len(genes)).astype(float)
                base[0] = rng.poisson(20) if sex == "female" else 0
                base[1:5] = rng.poisson(5, 4) if sex == "male" else 0
                rows.append(base)
                obs.append(dict(sex=sex, development_stage="3m", donor_id=f"{sex}{m}"))
    P.save_npz(tmp_path / "B2.npz", np.asarray(rows), pd.DataFrame(obs), genes)
    g3 = ["ERCC-00002", "ERCC-00003", "Cdk1", "Top2a", "Mki67", *[f"h{i}" for i in range(30)]]
    rows3, obs3 = [], []
    for phase, up in (("G1", 1.0), ("S", 1.3), ("G2M", 2.0)):
        for _ in range(60):
            x = rng.poisson(4.0, len(g3)).astype(float)
            x[2:5] = rng.poisson(4.0 * up, 3)
            rows3.append(x)
            obs3.append(dict(phase=phase, batch=phase))
    P.save_npz(tmp_path / "B3.npz", np.asarray(rows3), pd.DataFrame(obs3), g3)
    spec = {"B2": {"file": "B2.npz", "donor": "donor_id"}, "B3": {"file": "B3.npz", "donor": "batch"}}
    (tmp_path / "backgrounds.json").write_text(json.dumps(spec))
    A.main(["--backgrounds", str(tmp_path / "backgrounds.json"), "--data-dir", str(tmp_path),
            "--out", str(tmp_path / "anchors.json")])
    res = {r["claim"]: r for r in json.loads((tmp_path / "anchors.json").read_text())}
    assert res["R1 Xist"]["allowed"] == ["SUPPORTED"] and res["R1 Xist"]["cause"]
    assert set(res["R1 sham (female vs female)"]["allowed"]) == {"NO DETECTABLE EFFECT", "INCONCLUSIVE"}
    assert res["R2c total RNA, ERCC removed"]["label"] == "UNIDENTIFIABLE"
    assert res["R2b total RNA, ERCC present"]["label"] == "INCONCLUSIVE"  # one capture batch per phase
    assert res["R1 Xist"]["in_allowed"] and res["R2a G2M score"]["in_allowed"]
    assert res["R3"]["skipped"].startswith("B4 was dropped")
    assert (tmp_path / "anchors.md").exists()
    # R1's allowed outcome follows the engine's replicate rule over all strata (review 1, C11): an
    # age with no female mouse does not stop SUPPORTED when 3 females and 3 males are there in all
    obs = pd.DataFrame(dict(sex=["female"] * 3 + ["male"] * 4, age=["3m", "3m", "3m", "3m", "3m", "24m", "24m"],
                            mouse=[f"m{i}" for i in range(7)]))
    assert A.r1_allowed(obs) == ("SUPPORTED",)
    assert A.r1_allowed(obs[obs["mouse"] != "m0"]) == ("INCONCLUSIVE",)
    # R2a and R2b by the captures per phase (the third review: the branch of >= 3 was missing)
    one = pd.DataFrame(dict(phase=["G1"] * 4 + ["G2M"] * 4, batch=["b1"] * 4 + ["b2"] * 4))
    three = pd.DataFrame(dict(phase=["G1"] * 6 + ["G2M"] * 6, batch=[f"b{i % 3}" for i in range(6)]
                              + [f"c{i % 3}" for i in range(6)]))
    assert A.r2_allowed(one) == ("INCONCLUSIVE",)
    assert A.r2_allowed(three) == ("SUPPORTED", "INCONCLUSIVE")
    assert A.r2_allowed(three[~((three["phase"] == "G2M") & (three["batch"] == "c2"))]) == ("INCONCLUSIVE",)
