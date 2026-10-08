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
    out["delta"] = {str(k): {f"{f:g}": dict(value=delta * f, se=0.01) for f in (0.25, 0.5, 1.0, 1.5)}
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
    assert len(a) == P.n_datasets() == 5900
    assert P.n_cards() == 6000 and P.n_cards(data="null") == 5100  # N4 has two cards
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
    """For a given pair, cards of the norm_pearson conditions on B1 differ only in the id and the
    claimed direction; the pair, its level's SESOI and delta_min, the injection's strength and the
    controls come from the pair."""
    pilot = _pilot(bgs)
    for pair in (0, 9, 20):
        cards = [P.build(_entry(c, v, seed=3, pair=pair), bgs, pilot)[3][0]
                 for c, v in (("N1", "null"), ("N2", "c=0.5"), ("N3", "f=0.2"), ("N8", "beta(2,2)"),
                              ("E1", "dose=key"), ("E2", "against"), ("E3", "with"))]
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
    with pytest.raises(ValueError, match="outside"):
        P.build(_entry("N1", "null", pair=99), bgs, pilot)


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
    assert O.outcome_from_values((hi, lo), (lo, hi), True, 0.5, "decrease") == P.INCONCLUSIVE
    assert O.outcome_from_values((hi, lo), (hi, lo), False, 0.5, "decrease") == P.SUPPORTED
    assert O.outcome_from_values((hi, lo), (hi, lo), False, 0.5, "increase") == P.NS_OPPOSITE
    assert O.outcome_from_values((same, same), (same, same), False, 0.5, "increase") == P.NDE
    assert O.outcome_from_values((same, same), (same, same), False, 0.001, "increase") == P.INCONCLUSIVE


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
    import oracle as O
    d = O.delta_star(bgs["B1"], 0, 3.0, draws=200)
    assert d["value"] > 0.05 and d["se"] < d["value"] / 5
    assert O.delta_star(bgs["B1"], 0, 3.0, draws=200) == d  # public seed: reproducible


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
    for name, val in dict(CURVE_DATASETS=1, CURVE_REPS=2, TRUTH_DATASETS=3, TRUTH_REPS=4, GATE4_REPS=20,
                          DOSE_GRID=(1.0, 2.0, 3.0)).items():
        monkeypatch.setattr(O, name, val)
    pilot = O.run_pilot(small_bgs, n=3, draws=20)
    for b in ("B1", "B2"):
        assert pilot["pool_size"][b] == 12 and [pe["level"] for pe in pilot["pool"][b]] == [
            lv for lv in P.LEVELS for _ in range(4)]
        for lv in P.LEVELS:
            assert pilot["sesoi"][b][lv] in O.SESOI_GRID
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
    keys = {f"{c.name}:{v}:{k}" for c in P.CONDITIONS for v, _ in c.variants for k in range(12)}
    assert set(pilot["establishable"]) == keys
    assert not pilot["establishable"]["N4:3v3:0"]["establishable"]
    assert pilot["establishable"]["N6b:constant:5"]["establishable"]  # always DEGENERATE METRIC
    assert pilot["backgrounds"] == {k: P.background_sha256(b) for k, b in small_bgs.items()}


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


def _est_pilot(bgs, delta=0.3, truth=None):
    pilot = _pilot(bgs, delta=delta, truth=truth)
    pilot["establishable"] = {f"{c.name}:{v}:{k}": dict(establishable=c.oracle)
                              for c in P.CONDITIONS for v, _ in c.variants for k in range(24)}
    return pilot


def _reports(entries, fn):
    conds = P.conditions()
    return {cid: CAUSE_OF[fn(conds[e["condition"]], e)] for e in entries for cid in P.card_ids(e)}


def test_score_passes_a_perfect_engine_and_fails_the_bad_validators(bgs):
    """A validator that always gives a correct definite outcome passes; one that always rejects
    the metric, always says INCONCLUSIVE, always refuses, always says SUPPORTED or always NO
    DETECTABLE EFFECT fails (the reviewer's question: no bad validator passes)."""
    pytest.importorskip("scipy")
    import score as S
    pilot = _est_pilot(bgs, truth={"low": "blind"})
    entries = P.assign(KEY, pool_sizes=_sizes(bgs))
    order = (P.NS_INVALID, P.SUPPORTED, P.NS_DEPTH, P.NDE, P.NS_OPPOSITE, P.DEGENERATE)

    def perfect(c, e):
        good = P.definite(c, e["variant"], e["pair"], pilot)
        return next((o for o in order if o in good), P.INCONCLUSIVE)
    res = S.score(entries, _reports(entries, perfect), pilot)
    assert res["passed"] and all(res["criteria"][s]["passed"] for s in ("S1", "S2", "S3", "S4", "S5"))
    assert res["joint_pass_probability_sound"] >= 0.9
    for name, fn, fails in (
            ("always invalid", lambda c, e: P.NS_INVALID, {"S3", "S4", "S5"}),
            ("always inconclusive", lambda c, e: P.INCONCLUSIVE, {"S3"}),
            ("always refuses", lambda c, e: P.REFUSAL, {"S3"}),
            ("always supported", lambda c, e: P.SUPPORTED, {"S1", "S2", "S4"}),
            ("always no effect", lambda c, e: P.NDE, {"S3", "S4"})):
        got = S.score(entries, _reports(entries, fn), pilot)["criteria"]
        failed = {s for s in ("S1", "S2", "S3", "S4", "S5") if not got[s]["passed"]}
        assert fails <= failed, (name, failed)
    assert S.score(entries, {}, pilot)["criteria"]["S4"]["rate"] == 1.0  # no reports: every card an error


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
    assert res["criteria"]["S3"]["definite_any"]["rate"] > res["criteria"]["S3"]["rate"]


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
    s5 = res["criteria"]["S5"]
    n_valid = sum(P.metric_truth(P.conditions()[e["condition"]], e["pair"], pilot) == "valid"
                  for e in entries for _ in P.card_ids(e))
    assert s5["n"] == n_valid and s5["k"] == n_valid and not s5["passed"]
    assert res["per_truth"]["blind"]["errors"]["rate"] == 0.0 and res["per_truth"]["ambiguous"]["errors"]["rate"] == 0.0


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
    assert oc.nominal_error(frozenset({P.NDE, P.INCONCLUSIVE, P.NS_OPPOSITE, P.REFUSAL})) == pytest.approx(0.125)
    assert oc.nominal_error(P.INVALID_ALLOWED) == pytest.approx(0.075)
    for scenario in ((True, 0.5, 0.5), (False, 0.0, 1.0)):
        assert oc.joint_pass_probability(oc.expected_rows(oc.scenario_pilot(*scenario)), sims=4000) >= 0.90


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
    """The thread variables count only if they are set before numpy loads its BLAS; a forked
    worker inherits the pool numpy was loaded with."""
    env = {k: v for k, v in os.environ.items() if k not in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
    code = "import json, run_panel; print(json.dumps(run_panel.machine()))"
    out = subprocess.run([sys.executable, "-c", code], cwd=HERE, env=env, capture_output=True,
                         text=True, check=True).stdout
    m = json.loads(out)
    assert set(m["blas_env"].values()) == {"1"}
    if m["blas_threads_in_use"] is not None:  # read with threadpoolctl where it is installed
        assert m["blas_threads_in_use"] == [1]


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
    its key, give the published reports' sha256; a manifest whose hash was altered, or another pilot,
    does not pass."""
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
    m = json.loads((tmp_path / "res" / "manifest.json").read_text())
    m["datasets"][0]["cards"][0]["report_sha256"] = "0" * 64
    (tmp_path / "res" / "manifest.json").write_text(json.dumps(m))
    res = B.verify(tmp_path / "res", compact, pilot, tmp_path / "rerun2", workers=1)
    assert res["datasets"] == 3 and res["identical"] == res["cards"] - 1
    with pytest.raises(SystemExit, match="pilot"):
        B.verify(tmp_path / "res", compact, dict(pilot, dropped=["N3 steps"]), tmp_path / "rerun3", workers=1)


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


def _fake_census(monkeypatch, obs_by_org, raw_ok=True):
    """A stand-in for cellxgene_census with the API the selection uses: obs read in Arrow chunks with
    dictionary columns (a TableReadIter with .concat()), and get_anndata returning only the obs columns
    asked for, on a string index."""
    import sys
    import types
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
        rng = np.random.default_rng(0)
        X = rng.poisson(1.0, size=(len(obs), 30)).astype(np.float32)
        if not raw_ok:
            X = X * 0.5
        return Ad(sparse.csr_matrix(X), obs, pd.DataFrame({"feature_name": [f"G{j}" for j in range(30)]}))

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
    rank when a candidate's counts are not raw."""
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
    spec, lines = SB.select_census(tmp_path)
    assert spec["B1"]["source"]["dataset_id"] == "h2" and spec["B1"]["source"]["rank"] == 1
    assert spec["B1"]["source"]["extracted_cells"] == 25 * 230  # tie on 25 donors: more cells per donor
    assert spec["B2"]["source"]["qualifying_donors"] == 13 and spec["B2"]["source"]["rule"].startswith("all donors")
    z = np.load(tmp_path / "B1.npz")
    assert P.sha256(tmp_path / "B1.npz") == spec["B1"]["sha256"] and set(z["obs_cell_type"]) == {"B cell"}
    assert z["X_shape"][0] == spec["B1"]["source"]["extracted_cells"] and len(set(z["obs_donor_id"])) == 25
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
