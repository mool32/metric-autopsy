"""Tests of the confirmatory panel's own code: panel.py, oracle.py, score.py, oc.py,
run_panel.py and blind.py.

They run on simulated backgrounds (simulate.py; no network). The panel, the oracle and the
scoring must not import the engine; scoring and the oracle need scipy (skipped without it, as in
the core-only CI job).
"""
from __future__ import annotations

import ast
import hashlib
import json
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

KEY = "00112233445566778899aabbccddeeff"
OTHER_KEY = "ffeeddccbbaa99887766554433221100"


def _pilot(bgs, sesoi=0.15, dose=3.0, delta=0.3):
    pool = bgs["B1"].plan["pool"]
    return dict(sesoi={lv: sesoi for lv in P.LEVELS}, key_dose={lv: dose for lv in P.LEVELS},
                delta={str(k): {f"{f:g}": dict(value=delta * f, se=0.01) for f in (0.25, 0.5, 1.0, 1.5)}
                       for k in range(len(pool))}, establishable={}, dropped=[])


@pytest.fixture(scope="module")
def bgs():
    return {"B1": simulated_background("B1"), "B2": simulated_background("B2", donors=12, seed=1)}


def _entry(condition, variant, side="A", seed=0, i=0, pair=0):
    return dict(id=f"T{i}", condition=condition, variant=variant, index=i, side=side, seed=seed, pair=pair)


# --------------------------------------------------------------------------- #
# the key
# --------------------------------------------------------------------------- #
def test_the_key_is_128_bits_of_hex_and_its_commitment_is_the_sha256_of_that_string():
    for bad in (1, "0" * 31, "0" * 33, "G" * 32, "AB" * 16, " " + "0" * 32):
        with pytest.raises(ValueError):
            P.check_key(bad)
    assert P.key_commitment(KEY) == hashlib.sha256(KEY.encode("utf-8")).hexdigest()


def test_assign_is_deterministic_complete_and_keyed():
    a, b, c = P.assign(KEY), P.assign(KEY), P.assign(OTHER_KEY)
    assert a == b and a != c
    assert len(a) == P.n_datasets() == 5900
    assert P.n_cards() == 6000 and P.n_cards(null=True) == 5100  # N4 has two cards
    assert len({e["id"] for e in a}) == len(a) and {e["id"] for e in a} == {e["id"] for e in c}
    counts = Counter((e["condition"], e["variant"]) for e in a)
    assert counts == {(cond.name, v): n for cond in P.CONDITIONS for v, n in cond.variants}
    assert sum(x["condition"] != y["condition"] for x, y in zip(a, c)) > 0.5 * len(a)
    assert {e["side"] for e in a} == {"A", "B"}
    assert {e["pair"] for e in a} == set(range(len(P.LEVELS) * P.PAIRS_PER_LEVEL))


def test_key_null_conditions_have_790_datasets_and_n8_is_one_of_them():
    keyed = {c.name for c in P.CONDITIONS if c.key}
    assert keyed == {"N1", "N2", "N5", "N6c", "N8"}
    for cond in P.CONDITIONS:
        if cond.key:
            assert dict(cond.variants)[cond.key_variant] == P.KEY_N == 790


def test_the_gene_pair_is_drawn_independently_of_the_condition():
    """The pair (and so its level) does not depend on the condition: the pair counts of N1 and
    of the real effects are both uniform over the pool."""
    stats = pytest.importorskip("scipy.stats")
    a = P.assign(KEY)
    k = len(P.LEVELS) * P.PAIRS_PER_LEVEL
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
def test_the_plan_finds_disjoint_pairs_at_three_levels_with_their_controls(bgs):
    plan = bgs["B1"].plan
    pool = plan["pool"]
    assert [pe["level"] for pe in pool] == [lv for lv in P.LEVELS for _ in range(P.PAIRS_PER_LEVEL)]
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
    """For a given pair, cards of the norm_pearson conditions differ only in the id and the
    claimed direction; the pair, its level's SESOI and its controls come from the pair."""
    pilot = _pilot(bgs)
    for pair in (0, 5, 10):
        cards = [P.build(_entry(c, v, seed=3, pair=pair), bgs, pilot)[3][0]
                 for c, v in (("N1", "null"), ("N2", "c=0.5"), ("N3", "f=0.2"), ("N8", "beta(2,2)"),
                              ("E1", "dose=key"), ("E2", "against"), ("E3", "with"))]
        strip = [{k: v for k, v in c.items() if k not in ("id", "prereg")} for c in cards]
        assert all(s == strip[0] for s in strip)
        assert all({k: v for k, v in c["prereg"].items() if k != "direction"} ==
                   {k: v for k, v in cards[0]["prereg"].items() if k != "direction"} for c in cards)
        assert cards[0]["gene_pair"] == bgs["B1"].plan["pool"][pair]["pair"]
    dirs = Counter(P.build(_entry("N1", "null", seed=s, i=s), bgs, pilot)[3][0]["prereg"]["direction"]
                   for s in range(30))
    assert set(dirs) == {"increase", "decrease"}
    assert P.build(_entry("E1", "dose=key", side="A"), bgs, pilot)[3][0]["prereg"]["direction"] == "decrease"
    assert P.build(_entry("E1", "dose=key", side="B"), bgs, pilot)[3][0]["prereg"]["direction"] == "increase"
    # a negative Δ* makes the signal side the lower one
    neg = _pilot(bgs, delta=-0.3)
    assert P.build(_entry("E1", "dose=key", side="A"), bgs, neg)[3][0]["prereg"]["direction"] == "increase"


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
    t8, t2 = obs8.groupby("group")["total_counts"], obs2.groupby("group")["total_counts"]
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


def test_the_panel_the_oracle_and_the_scoring_never_import_the_engine():
    for name in ("panel.py", "oracle.py", "score.py", "oc.py", "simulate.py"):
        tree = ast.parse((HERE / name).read_text())
        mods = {a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
        mods |= {n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module}
        assert not any(m.split(".")[0] == "metric_autopsy" for m in mods), (name, mods)


# --------------------------------------------------------------------------- #
# allowed sets and the oracle
# --------------------------------------------------------------------------- #
def test_no_detectable_effect_is_allowed_for_a_real_effect_only_below_the_sesoi(bgs):
    conds = P.conditions()
    small, large = _pilot(bgs, sesoi=0.15, delta=0.1), _pilot(bgs, sesoi=0.15, delta=0.3)
    for name, variant in (("E1", "dose=key"), ("E2", "against"), ("E3", "with")):
        assert "NO DETECTABLE EFFECT" in P.allowed(conds[name], variant, 0, small)
        assert "NO DETECTABLE EFFECT" not in P.allowed(conds[name], variant, 0, large)
        assert "INCONCLUSIVE" not in P.definite(conds[name], variant, 0, large)
    # E1's dose ladder reads Δ* at its own dose: 0.25 x 0.3 = 0.075 < 0.15
    assert "NO DETECTABLE EFFECT" in P.allowed(conds["E1"], "dose=0.25", 0, large)
    assert P.allowed(conds["N8"], "beta(2,2)", 0, large) == P.NULL_ALLOWED


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
    # explained: a detected raw difference that the true correction removes
    assert O.label_from_values((hi, lo), (same, same + 0.001), True, 0.5, "decrease") == "NOT SUPPORTED"
    # reversed by the correction: INCONCLUSIVE
    assert O.label_from_values((hi, lo), (lo, hi), True, 0.5, "decrease") == "INCONCLUSIVE"
    # detected in and against the declared direction (A - B > 0 is a decrease)
    assert O.label_from_values((hi, lo), (hi, lo), False, 0.5, "decrease") == "SUPPORTED"
    assert O.label_from_values((hi, lo), (hi, lo), False, 0.5, "increase") == "NOT SUPPORTED"
    # equivalent within the SESOI, or not established
    assert O.label_from_values((same, same), (same, same), False, 0.5, "increase") == "NO DETECTABLE EFFECT"
    assert O.label_from_values((same, same), (same, same), False, 0.001, "increase") == "INCONCLUSIVE"


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


def test_delta_star_is_a_paired_estimate_with_a_small_standard_error(bgs):
    import oracle as O
    d = O.delta_star(bgs["B1"], 0, 3.0, draws=200)
    assert d["value"] > 0.05 and d["se"] < d["value"] / 5
    assert O.delta_star(bgs["B1"], 0, 3.0, draws=200) == d  # public seed: reproducible


def test_the_pilot_fixes_sesoi_dose_and_delta_per_level_and_finds_establishable_cases(bgs):
    pytest.importorskip("scipy")
    import oracle as O
    pilot = O.run_pilot(bgs, n=6, draws=20)
    assert set(pilot["sesoi"]) == set(P.LEVELS) and all(s in O.SESOI_GRID for s in pilot["sesoi"].values())
    for lv in P.LEVELS:
        assert pilot["key_dose"][lv] in O.DOSE_GRID
        if pilot["key_dose_found"][lv]:
            last = pilot["dose_scan"][lv][-1]
            assert last["power"] >= O.POWER and last["delta"] >= O.DELTA_MARGIN * pilot["sesoi"][lv]
        assert 0 <= pilot["n2_raw_power"][lv] <= 1
    assert set(pilot["delta"]) == {str(k) for k in range(len(bgs["B1"].plan["pool"]))}
    assert pilot["n2n3_informative"] == bool(pilot["informative_levels"])
    keys = {f"{c.name}:{v}:{lv}" for c in P.CONDITIONS for v, _ in c.variants for lv in P.LEVELS}
    assert set(pilot["establishable"]) == keys
    assert not pilot["establishable"]["N4:3v3:high"]["establishable"]
    assert pilot["establishable"]["N6b:constant:low"]["establishable"]  # always DEGENERATE METRIC


# --------------------------------------------------------------------------- #
# scoring and the criteria
# --------------------------------------------------------------------------- #
def test_score_labels():
    pytest.importorskip("scipy")
    import score as S
    assert S.label("SUPPORTED (provisional until replicated) [non-directional claim]") == "SUPPORTED"
    assert S.label("SUPPORTED — replicated") == "SUPPORTED"
    assert S.label("NOT SUPPORTED — metric invalid: x") == "NOT SUPPORTED"
    assert S.label("NO DETECTABLE EFFECT — effect +0.01 [underpowered relative to the SESOI]") == "NO DETECTABLE EFFECT"
    assert S.label(None) == "ERROR" and S.label("garbage") == "ERROR"


def _est_pilot(bgs, delta=0.3):
    pilot = _pilot(bgs, delta=delta)
    pilot["establishable"] = {f"{c.name}:{v}:{lv}": dict(establishable=c.oracle)
                              for c in P.CONDITIONS for v, _ in c.variants for lv in P.LEVELS}
    return pilot


def _verdicts(entries, pilot, fn):
    conds = P.conditions()
    return {cid: fn(conds[e["condition"]], e) for e in entries for cid in P.card_ids(e)}


def test_score_passes_a_perfect_engine_and_fails_the_failure_modes(bgs):
    pytest.importorskip("scipy")
    import score as S
    entries, pilot = P.assign(KEY), _est_pilot(bgs)

    def perfect(c, e):
        good = P.definite(c, e["variant"], e["pair"], pilot)
        return f"{good[0] if good else 'INCONCLUSIVE'} — x"
    res = S.score(entries, _verdicts(entries, pilot, perfect), pilot)
    assert res["passed"] and all(res["criteria"][s]["passed"] for s in ("S1", "S2", "S3", "S4"))
    shy = S.score(entries, _verdicts(entries, pilot, lambda c, e: "INCONCLUSIVE — x"), pilot)
    assert not shy["criteria"]["S3"]["passed"] and shy["criteria"]["S1"]["passed"]
    credulous = S.score(entries, _verdicts(entries, pilot, lambda c, e: "SUPPORTED (provisional until replicated)"), pilot)
    assert not credulous["criteria"]["S1"]["passed"] and not credulous["criteria"]["S2"]["passed"]
    assert not credulous["criteria"]["S4"]["passed"]
    assert S.score(entries, {}, pilot)["criteria"]["S4"]["rate"] == 1.0


def test_s3_counts_only_correct_definite_verdicts_and_nde_on_a_large_effect_is_an_error(bgs):
    pytest.importorskip("scipy")
    import score as S
    entries, pilot = P.assign(KEY), _est_pilot(bgs, delta=0.3)  # |Δ*| >= SESOI at the key dose

    def nde_everywhere(c, e):
        return "NO DETECTABLE EFFECT — x" if not c.name.startswith("N6") else "NOT SUPPORTED — x"
    res = S.score(entries, _verdicts(entries, pilot, nde_everywhere), pilot)
    e1 = res["per_condition"]["E1:dose=key"]
    assert e1["outside"]["rate"] == 1.0 and e1["correct_definite"]["rate"] == 0.0
    assert res["per_condition"]["E1:dose=0.25"]["outside"]["rate"] == 0.0  # 0.075 < SESOI 0.15
    assert res["criteria"]["S3"]["definite_any"]["rate"] > res["criteria"]["S3"]["rate"]


def test_s1_fails_when_any_one_key_condition_fails(bgs):
    pytest.importorskip("scipy")
    import score as S
    entries, pilot = P.assign(KEY), _est_pilot(bgs)
    n8 = [e for e in entries if e["condition"] == "N8"]
    bad = {e["id"] for e in n8[:40]}  # 40/790 > 29
    res = S.score(entries, {cid: ("SUPPORTED (provisional until replicated)" if e["id"] in bad else "INCONCLUSIVE — x")
                            for e in entries for cid in P.card_ids(e)}, pilot)
    assert not res["criteria"]["S1"]["conditions"]["N8"]["passed"] and not res["criteria"]["S1"]["passed"]
    assert all(v["passed"] for k, v in res["criteria"]["S1"]["conditions"].items() if k != "N8")


def test_the_criteria_meet_the_principle_and_s1_as_a_whole():
    pytest.importorskip("scipy")
    import oc
    k, p_all, p_doubled, _ = oc.joint_error_rule(790, 0.025, 5)
    assert k == 29 and p_all >= 0.90 and p_doubled <= 0.05
    assert oc.joint_error_rule(600, 0.025, 5)[1] < 0.90  # the old plan fails S1 as a whole
    counts = oc.design_counts()
    assert counts["key"] == {"N1": 790, "N2": 790, "N5": 790, "N6c": 790, "N8": 790}
    assert counts["all"] == 6000 and counts["null"] == 5100
    for n in (counts["null"], counts["all"]):
        _, ps, pd = oc.error_rule(n, 0.025)
        assert ps >= 0.90 and pd <= 0.05
    _, ps, pd = oc.decisiveness_rule(1000, 0.85)
    assert ps >= 0.90 and pd <= 0.05
    assert oc.joint_pass_probability(counts, 1000, sims=20000) >= 0.90


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
        # the expected rate over donor draws, by simulation of the same model
        lo, hi = S.cp(int(y.sum()), n)
        ov = S.overlap_interval(y, sets)
        p_true = np.mean(1 / (1 + np.exp(-(np.log(truth / (1 - truth)) + 1.5 * rng.normal(0, 1.2 / np.sqrt(8), 20000)))))
        cover_naive += lo <= p_true <= hi
        cover_adj += ov["ci95"][0] <= p_true <= ov["ci95"][1]
    assert np.median(deffs_indep) < 1.3
    assert cover_adj > cover_naive and cover_adj / reps >= 0.85


# --------------------------------------------------------------------------- #
# the runner and the blind run
# --------------------------------------------------------------------------- #
def test_runner_runs_the_engine_on_the_fly(bgs, tmp_path):
    import run_panel as R
    pilot = _pilot(bgs, dose=6.0)
    entries = P.assign(KEY)
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
        assert rep["verdict"] and rep["elapsed_seconds"] > 0
    again = R.run(pick, bgs, pilot, tmp_path / "out", workers=1)
    assert again["skipped"] == 3  # one attempt per card
    assert len((tmp_path / "out" / "runlog.jsonl").read_text().splitlines()) == 3


def test_runner_pins_one_blas_thread_per_worker_before_numpy_loads():
    """The thread variables count only if they are set before numpy loads its BLAS; a forked
    worker inherits the pool numpy was loaded with."""
    import os
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


def _fake_shard(root: Path, k: int, shards: int, key_sha: str, ids: list):
    d = root / f"shard-{k}"
    (d / "reports").mkdir(parents=True)
    rows = []
    for i in ids:
        text = json.dumps(dict(id=i, verdict="INCONCLUSIVE — x"))
        (d / "reports" / f"{i}.json").write_text(text)
        rows.append(dict(id=i, data_sha256="0" * 64, donors=["a"], cards=[
            dict(id=i, card_sha256="1" * 64, report_sha256=hashlib.sha256(text.encode()).hexdigest(), seconds=1.0)]))
    (d / "manifest.json").write_text(json.dumps(dict(datasets=rows)))
    (d / "runlog.jsonl").write_text("".join(json.dumps(dict(claim_id=i)) + "\n" for i in ids))
    meta = dict(shard=k, shards=shards, key_sha256=key_sha, pilot_sha256="p", backgrounds={"B1": {}},
                provenance={}, summary=dict(wall_seconds=1.0, machine=dict(cpus=4), workers=4))
    (d / "shard.json").write_text(json.dumps(meta))


def test_collect_merges_shards_and_refuses_inconsistent_ones(tmp_path):
    import blind as B
    _fake_shard(tmp_path / "ok", 0, 2, "k", ["D1", "D3"])
    _fake_shard(tmp_path / "ok", 1, 2, "k", ["D2"])
    s = B.collect(tmp_path / "ok", tmp_path / "res", expected_datasets=3)
    assert s["datasets"] == 3 and (tmp_path / "res" / "SHA256SUMS").exists()
    assert [r["id"] for r in json.loads((tmp_path / "res" / "manifest.json").read_text())["datasets"]] == ["D1", "D2", "D3"]
    with pytest.raises(SystemExit):
        B.collect(tmp_path / "ok", tmp_path / "res2", expected_datasets=4)
    _fake_shard(tmp_path / "keys", 0, 2, "k", ["D1"])
    _fake_shard(tmp_path / "keys", 1, 2, "other", ["D2"])
    with pytest.raises(SystemExit, match="key_sha256"):
        B.collect(tmp_path / "keys", tmp_path / "res3")
    _fake_shard(tmp_path / "missing", 0, 3, "k", ["D1"])
    with pytest.raises(SystemExit, match="shards present"):
        B.collect(tmp_path / "missing", tmp_path / "res4")
    _fake_shard(tmp_path / "twice", 0, 2, "k", ["D1"])
    _fake_shard(tmp_path / "twice", 1, 2, "k", ["D1"])
    with pytest.raises(SystemExit, match="twice"):
        B.collect(tmp_path / "twice", tmp_path / "res5")


def test_guard_and_key_commitment_follow_the_frozen_tag(tmp_path, monkeypatch):
    """In a scratch repository: the guard passes when only data files changed after the tag and
    fails when frozen code changed; the key must match the commitment in the tag's message."""
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
    git("add", "-A")
    git("commit", "-qm", "frozen")
    git("tag", "-a", "v0.3.0-prereg", "-m", f"frozen\n\nkey sha256: {P.key_commitment(KEY)}\n")
    monkeypatch.setattr(B, "HERE", repo)
    (repo / "validation" / "prereg" / "pilot.json").write_text("{}")
    git("add", "-A")
    git("commit", "-qm", "data only")
    assert B.guard("v0.3.0-prereg")["frozen_tag"] == "v0.3.0-prereg"
    monkeypatch.setenv("KEY_SEED", KEY)
    assert B.read_key(False, "v0.3.0-prereg") == KEY
    monkeypatch.setenv("KEY_SEED", OTHER_KEY)
    with pytest.raises(SystemExit, match="commitment"):
        B.read_key(False, "v0.3.0-prereg")
    (repo / "src" / "engine.py").write_text("x = 2\n")
    git("add", "-A")
    git("commit", "-qm", "engine change")
    with pytest.raises(SystemExit, match="frozen files differ"):
        B.guard("v0.3.0-prereg")
