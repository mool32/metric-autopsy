"""v0.4 change 9 (validation/plans/v0.4_spec.md, section 9; coded on the owner's decision of 2026-11-05): a
coupling probe on a set of genes.

``injected_signal.coupling_set`` thins the genes it couples with one per-cell keep probability (a shared
factor) against a sham that thins each gene of the set with its own; ``couple`` fixes which genes it couples
("all", or "half": a random half at every injection). Its dose follows change 3: the strength at which the
factor raises the mean correlation of log expected expression over the coupled pairs by the SESOI. A set of
two coupled whole is the pair's probe, draw for draw.

The synthetic data: mice nested in age, a set of six genes sharing a latent factor whose share of their
variance is the set's coordination, one value per group."""
from __future__ import annotations

import functools
from functools import partial
from itertools import combinations

import numpy as np
import pandas as pd
import pytest

from metric_autopsy import SimpleData, injected_signal, metrics, run_autopsy
from metric_autopsy.core import as_dense, unique_col_index
from metric_autopsy.gates import injected_deltas
from metric_autopsy.report import decide

SET = [f"S{i}" for i in range(6)]
N_FILLER = 100
PREREG = dict(estimand="composition", direction="decrease", sesoi=0.1, judgment_pending=False)


def set_data(rho=(0.3, 0.3), mu=-2.0, n_mice=6, cells=300, var_bio=0.6, depth=3000.0, seed=0):
    """Mice nested in age (young, old): the set's log expected expression is `mu` + a latent factor shared
    by its genes, whose share of their variance `var_bio` is the set's coordination rho (one per group);
    fillers of lower expression."""
    rng = np.random.default_rng(seed)
    fill = np.random.default_rng(1).normal(-3.0, 1.0, N_FILLER)
    blocks, rows = [], []
    for gi, age in enumerate(("young", "old")):
        for m in range(n_mice):
            lat = rng.normal(size=(cells, 1))
            load, noise = np.sqrt(rho[gi] * var_bio), np.sqrt((1 - rho[gi]) * var_bio)
            ls = mu + load * lat + noise * rng.normal(size=(cells, len(SET)))
            lf = fill[None, :] + rng.normal(0, 0.5, (cells, N_FILLER))
            rel = np.exp(np.hstack([ls, lf]))
            rel /= rel.sum(axis=1, keepdims=True)
            D = depth * np.exp(rng.normal(0, 0.3, cells))
            blocks.append(rng.poisson(D[:, None] * rel).astype(float))
            rows += [dict(age=age, mouse=f"{age}{m}")] * cells
    return SimpleData(np.vstack(blocks), pd.DataFrame(rows), SET + [f"F{i}" for i in range(N_FILLER)])


def mean_norm_pearson(data, *, genes):
    """A metric of a set's coordination: the mean of norm_pearson over its pairs."""
    v = [metrics.norm_pearson(data, gene_a=a, gene_b=b) for a, b in combinations(genes, 2)]
    v = [x for x in v if np.isfinite(x)]
    return float(np.mean(v)) if v else float("nan")


def mean_share_correlation(data, *, genes):
    """A metric that reads what is relative within the set: the mean correlation of the genes' shares of
    the set's counts in each cell."""
    X = as_dense(data.X)[:, [unique_col_index(list(data.var_names), g) for g in genes]]
    t = X.sum(axis=1)
    P = X[t > 0] / t[t > 0, None]
    return float(np.corrcoef(P.T)[np.triu_indices(len(genes), 1)].mean())


def mean_log_level(data, *, genes):
    """A metric of the set's level: the mean log1p CP10k of its genes (a coupling does not move it)."""
    X = as_dense(data.X)
    tot = X.sum(axis=1)
    idx = [unique_col_index(list(data.var_names), g) for g in genes]
    return float(np.log1p(1e4 * X[:, idx] / tot[:, None]).mean())


@functools.lru_cache(maxsize=None)
def claim(rho: tuple, metric: str = "mean_norm_pearson", couple: str = "all"):
    fn = dict(mean_norm_pearson=mean_norm_pearson, mean_log_level=mean_log_level)[metric]
    return run_autopsy(partial(fn, genes=SET), set_data(rho=rho, seed=2), group_col="age", groups=("young", "old"),
                       replicate_col="mouse", signal_test=injected_signal.coupling_set(SET, couple=couple),
                       prereg=dict(PREREG), log_path="off", specificity=False, gene_metric=fn)


def _g4(a):
    return next(r for r in a.results if r.gate == 4)


# --------------------------------------------------------------------------- #
# the probe
# --------------------------------------------------------------------------- #
def test_a_set_of_two_coupled_whole_is_the_pairs_probe_draw_for_draw():
    X = np.random.default_rng(0).poisson(5, (200, 2)).astype(float)
    pair, two = injected_signal.coupling("A", "B", 1.3), injected_signal.coupling_set(["A", "B"], 1.3)
    assert np.array_equal(pair.columns(X, np.random.default_rng(5)), two.columns(X, np.random.default_rng(5)))
    assert np.array_equal(pair.sham.columns(X, np.random.default_rng(5)), two.sham.columns(X, np.random.default_rng(5)))
    assert two.kind == "coupling_set" and two.direction == "increase" and two.genes == ("A", "B")


def test_the_set_probe_refuses_what_it_cannot_couple():
    for bad, kw in ((["A"], {}), (["A", "A", "B"], {}), (SET, dict(couple="some")), (["A", "B", "C"], dict(couple="half"))):
        with pytest.raises(ValueError):
            injected_signal.coupling_set(bad, **kw)


@pytest.mark.parametrize("couple, coupled", [("all", 6), ("half", 3)])
def test_the_factor_couples_the_set_or_one_half_and_the_sham_none(couple, coupled):
    # deep, equal counts: the coupled genes keep the same share of them in each cell, the others do not
    X = np.full((400, 6), 1000.0)
    probe = injected_signal.coupling_set(SET, strength=4.0, couple=couple)
    for cols, n in ((probe.columns(X, np.random.default_rng(1)), coupled), (probe.sham.columns(X, np.random.default_rng(1)), 0)):
        C = np.corrcoef(cols.T)[np.triu_indices(6, 1)]
        assert int((C > 0.9).sum()) == n * (n - 1) // 2 and int((np.abs(C) < 0.2).sum()) == 15 - n * (n - 1) // 2


def test_a_factor_shared_by_the_whole_set_cancels_in_its_shares_and_a_half_moves_them():
    sub = set_data()
    sub = sub[np.asarray(sub.obs["mouse"]) == "young0"]
    resp = {}
    for couple in ("all", "half"):
        d = injected_deltas(partial(mean_share_correlation, genes=SET), sub,
                            injected_signal.coupling_set(SET, 2.0, couple), 40, np.random.default_rng(3))
        resp[couple] = float(np.mean(d))
    assert abs(resp["all"]) < 0.002 and resp["half"] > 0.005
    corr = {c: float(np.mean(injected_deltas(partial(mean_norm_pearson, genes=SET), sub,
                                             injected_signal.coupling_set(SET, 2.0, c), 10, np.random.default_rng(3))))
            for c in ("all", "half")}
    assert corr["all"] > corr["half"] > 0  # a metric of the pairs' correlation reads both


# --------------------------------------------------------------------------- #
# the dose (change 3's rule over the coupled pairs)
# --------------------------------------------------------------------------- #
def test_the_dose_raises_the_mean_correlation_over_the_pairs_by_the_change():
    pair = injected_signal.coupling_strength_for(0.1, 0.3, 0.3)
    for k in (2, 6):  # equal variances: the pair's dose, whatever the set's size
        assert injected_signal.coupling_set_strength_for(0.1, [0.3] * k)["strength"] == pytest.approx(pair["strength"])
    var = [0.05, 0.2, 0.5, 1.0, 0.01, float("nan")]
    d = injected_signal.coupling_set_strength_for(0.15, var)
    v, vs = d["added_variance"], np.array(d["variances"])
    rises = [v / np.sqrt((vs[i] + v) * (vs[j] + v)) for i, j in combinations(range(6), 2)]
    assert np.mean(rises) == pytest.approx(0.15) and d["n_pairs"] == 15 and d["floored"] == [5]
    assert injected_signal.coupling_set_strength_for(0.3, var)["strength"] > d["strength"]
    with pytest.raises(ValueError):
        injected_signal.coupling_set_strength_for(1.0, var)


def test_the_probe_is_redosed_keeping_its_genes_and_what_it_couples():
    probe, info = injected_signal.at_construct_change(injected_signal.coupling_set(SET, couple="half"), 0.1,
                                                      variances=[0.2] * 5 + [float("nan")])
    assert probe.kind == "coupling_set" and probe.couple == "half" and probe.genes == tuple(SET)
    assert info["kind"] == "coupling_set" and info["floored"] == ["S5"] and probe.strength == info["strength"]
    assert "half it couples" in info["scale"]


# --------------------------------------------------------------------------- #
# a claim about a set's coordination, through the autopsy
# --------------------------------------------------------------------------- #
def test_a_decrease_of_a_sets_coordination_is_supported_with_the_probe_dosed_to_the_sesoi():
    a = claim((0.5, 0.15))
    g4 = _g4(a)
    assert g4.status.value == "PASS" and g4.detail["dose"]["kind"] == "coupling_set" and g4.detail["dose"]["n_pairs"] == 15
    assert 0 < g4.detail["sesoi_metric"] < 0.1 and g4.detail.get("dose_response")
    assert a.effect.detail["effect"] > 0 and decide(a).startswith("SUPPORTED")
    assert not any(r.gate == 5 for r in a.results)  # GATE 5's controls are pairs: no control for a set
    assert "## The probe's dose and the metric's response" in a.to_markdown()


def test_an_equal_coordination_shows_no_detectable_effect_within_the_sets_s_m():
    a = claim((0.3, 0.3))
    assert _g4(a).status.value == "PASS" and decide(a).startswith("NO DETECTABLE EFFECT")


def test_a_metric_of_the_sets_level_is_not_shown_to_respond_to_its_coupling():
    assert _g4(claim((0.3, 0.3), "mean_log_level")).status.value != "PASS"
