"""Synthetic tests: plant a known confound and assert the gates catch it.

The generators build single-cell-like count matrices where the *biology* is known by
construction, so we can assert that:
  - a QC confound hidden in a sex x age interaction is found by GATE 1 (a diagnostic),
  - a metric whose effect is pure QC is explained by depth once depth is equalized,
  - mi_3bin's dropout sensitivity is measured by GATE 0 (as attenuation) while a robust
    metric shows none,
  - a clean dataset with matched QC and real signal reaches SUPPORTED only with the full
    evidence (estimand, replicates, demonstrated metric response, resolved judgment).

Revised for the v0.3 verdict scheme: each test keeps the truth it checked in v0.1.1 and
states the v0.3 meaning of the verdict; see the per-test notes.
"""
from functools import partial

import numpy as np
import pandas as pd
import pytest

from metric_autopsy import (
    SimpleData, GateStatus, metrics, run_autopsy,
    gate0_independence, gate1_qc_parity, gate2_ngenes_matching, gate5_controls,
    gate6_replication,
)

GENES = ["Smad3", "Col1a1", "Actb", "Gapdh"] + [f"Gene{i}" for i in range(36)]
MI = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")
NPR = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
CONTROLS = dict(pair_metric=metrics.norm_pearson, pos_pair=("Actb", "Gapdh"), neg_pair=("Gene0", "Gene1"))


def _block(rng, n, coupling, efficiency, age, sex):
    """One (age, sex) block. `coupling` sets Smad3<->Col1a1 covariance (the biology);
    `efficiency` scales capture (the QC nuisance -> lower efficiency = more dropout)."""
    latent = rng.normal(0, 1, n)
    smad = np.exp(1.2 + coupling * 0.6 * latent + rng.normal(0, 0.3, n))
    col = np.exp(1.2 + coupling * 0.6 * latent + rng.normal(0, 0.3, n))
    hk = rng.normal(0, 1, n)  # housekeeping co-regulation, always on
    actb = np.exp(1.6 + 0.9 * hk + rng.normal(0, 0.2, n))
    gapdh = np.exp(1.6 + 0.9 * hk + rng.normal(0, 0.2, n))
    filler = np.exp(0.4 + rng.normal(0, 0.5, (n, 36)))
    lam = np.column_stack([smad, col, actb, gapdh, filler]) * efficiency
    counts = rng.poisson(lam).astype(float)
    obs = pd.DataFrame({"age": [age] * n, "sex": [sex] * n})
    return counts, obs


def _assemble(blocks):
    Xs = [b[0] for b in blocks]
    obs = pd.concat([b[1] for b in blocks], ignore_index=True)
    return SimpleData(np.vstack(Xs), obs, GENES)


def make_confounded(seed=0, n=400):
    """Biology identical everywhere; only male-old cells are QC-degraded (efficiency 0.3)."""
    rng = np.random.default_rng(seed)
    return _assemble([
        _block(rng, n, coupling=1.5, efficiency=1.0, age="young", sex="male"),
        _block(rng, n, coupling=1.5, efficiency=1.0, age="young", sex="female"),
        _block(rng, n, coupling=1.5, efficiency=0.30, age="old", sex="male"),   # <- confound
        _block(rng, n, coupling=1.5, efficiency=1.0, age="old", sex="female"),
    ])


def add_mice(data, per_block=4, cols=("age", "sex"), seed=0):
    """Split every (age, sex) block into `per_block` mice: a replicate column. Cells are iid
    within a block, so the mice add no between-mouse variance."""
    rng = np.random.default_rng(seed)
    obs = data.obs.copy()
    key = obs[list(cols)].astype(str).agg("-".join, axis=1)
    mouse = np.empty(len(obs), dtype=object)
    for k in pd.unique(key):
        idx = np.where(key.values == k)[0]
        mouse[idx] = [f"{k}-m{j}" for j in rng.permutation(np.arange(len(idx)) % per_block)]
    obs["mouse"] = mouse
    return SimpleData(data.X, obs, data.var_names)


def make_clean(seed=1, n=400):
    """Real coupling difference (young strong, old weak), QC matched across all groups."""
    rng = np.random.default_rng(seed)
    return _assemble([
        _block(rng, n, coupling=2.6, efficiency=1.0, age="young", sex="male"),
        _block(rng, n, coupling=2.6, efficiency=1.0, age="young", sex="female"),
        _block(rng, n, coupling=0.3, efficiency=1.0, age="old", sex="male"),
        _block(rng, n, coupling=0.3, efficiency=1.0, age="old", sex="female"),
    ])


# --------------------------------------------------------------------------- #
def test_gate1_catches_factorial_confound():
    """The QC confound hidden in the sex x age interaction is found, in the male stratum only.
    v0.3: GATE 1 is a design diagnostic (WARN), not a blocking FAIL; the imbalance goes to
    the estimand-dependent correction."""
    data = make_confounded()
    # Pooled (no strata) also flags here, but the point is the sex-stratified check.
    res = gate1_qc_parity(data, "age", ("young", "old"), within=["sex"])
    assert res.status == GateStatus.WARN
    flagged_sexes = {r["stratum"]["sex"] for r in res.detail["flagged"]}
    assert flagged_sexes == {"male"}, f"only male stratum should flag, got {flagged_sexes}"


def test_gate2_blocks_in_confounded_stratum():
    """Within the confounded stratum (males), young/old n_genes don't overlap -> STOP.
    This is the 'only 3 cells in overlap' failure; the groups are simply incomparable."""
    data = make_confounded()
    males = data[np.asarray(data.obs["sex"]) == "male"]
    m = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")
    res = gate2_ngenes_matching(m, males, "age", ("young", "old"))
    assert res.blocking  # STOP (no overlap) or FAIL (effect collapses)
    assert res.status == GateStatus.STOP


def test_gate2_pooled_dilutes_the_confound():
    """Pooling old = male+female hides the male-only confound: the apparent effect is ~0,
    so the pooled test 'passes' — the lesson being that you must stratify (GATE 1)."""
    data = make_confounded()
    m = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")
    res = gate2_ngenes_matching(m, data, "age", ("young", "old"))
    assert abs(res.detail["unmatched_effect"]) < 0.02  # confound diluted away by pooling


def test_gate0_separates_mi_sensitivity_from_robust_metric():
    """On QC-matched data GATE 0 isolates intrinsic nuisance sensitivity: dropout and depth
    shrink mi_3bin's signal heavily and norm_pearson's not measurably. v0.3: shrinkage
    toward the null is attenuation, reported with its size and used as a power check; it is
    not bias, so neither metric FAILs (v0.1.1 failed mi_3bin here)."""
    data = make_clean()  # matched QC, so GATE 0 isolates intrinsic nuisance-sensitivity
    mi = gate0_independence(MI, data)
    robust = gate0_independence(NPR, data)
    assert mi.detail["attenuation"]["extra_dropout"] > 0.4
    assert mi.detail["attenuation"]["depth_downsample"] > 0.1
    assert robust.detail["attenuation"] == {}
    assert not any(r["confounded"] for r in mi.detail["responses"].values())
    assert mi.status == robust.status == GateStatus.PASS


def test_gate5_controls_positive_fires_negative_null():
    data = make_confounded()
    res = gate5_controls(
        metrics.mi_3bin, data,
        pos_pair=("Actb", "Gapdh"), neg_pair=("Gene0", "Gene1"), within=["sex"],
    )
    assert res.status == GateStatus.PASS, res.message


def test_full_autopsy_confounded_is_not_supported_at_equal_depth():
    """Biology identical everywhere; only male-old capture is degraded. The apparent MI
    difference must be diagnosed, not believed. v0.3 thins to equal depth within each sex
    (composition estimand) and infers at the mouse level: nothing is detected at equal
    depth and most of the raw difference is gone. When the raw difference is itself
    detected across mice the verdict is NOT SUPPORTED — explained by depth; when it is not
    (4 mice per group, 100 cells each: borderline), INCONCLUSIVE. Never SUPPORTED, for any
    split of the cells into mice. (v0.1.1 said FAIL without saying why.)"""
    for split in range(3):
        data = add_mice(make_confounded(), seed=split)
        autopsy = run_autopsy(
            MI, data, group_col="age", groups=("young", "old"), within=["sex"],
            replicate_col="mouse", prereg={"estimand": "composition"}, stop_on_first_fail=False,
        )
        eff = autopsy.effect
        assert autopsy.design_adequacy.status == "CORRECTED"
        assert eff.status != "DETECTED"
        assert abs(eff.detail["retained"]) < 0.5
        assert eff.detail["explained_by_depth"] == eff.detail["raw_detected"]
        expected = "NOT SUPPORTED" if eff.detail["raw_detected"] else "INCONCLUSIVE"
        assert autopsy.verdict.startswith(expected), autopsy.verdict


def test_full_autopsy_clean_is_supported_only_with_full_evidence():
    """A real coupling difference with matched QC. v0.3 meaning of the positive verdict:
    SUPPORTED needs a declared estimand, replicate-level inference, a metric whose response
    to signal is demonstrated (GATE 5 positive control) and resolved judgment gates. Remove
    any one of them and the verdict must not be SUPPORTED. (v0.1.1 said PASS as soon as no
    automatic gate had failed, on cells, with no control.)"""
    data = add_mice(make_clean())
    full = dict(group_col="age", groups=("young", "old"), within=["sex"],
                gene_pair=("Smad3", "Col1a1"), replicate_col="mouse", **CONTROLS,
                prereg={"estimand": "composition", "direction": "decrease", "judgment_pending": False})
    autopsy = run_autopsy(NPR, data, **full)
    assert autopsy.metric_validity.status == "PASS"
    assert autopsy.design_adequacy.status == "ADEQUATE"
    assert autopsy.effect.status == "DETECTED"
    assert autopsy.verdict == "SUPPORTED (provisional until replicated)", autopsy.verdict

    without = {
        "estimand": {**full, "prereg": {"direction": "decrease", "judgment_pending": False}},
        "direction": {**full, "prereg": {"estimand": "composition", "judgment_pending": False}},
        "replicates": {**full, "replicate_col": None},
        "controls": {k: v for k, v in full.items() if k not in CONTROLS},
        "judgment": {**full, "prereg": {"estimand": "composition", "direction": "decrease"}},
    }
    for missing, kw in without.items():
        verdict = run_autopsy(NPR, data, **kw).verdict
        assert not verdict.startswith("SUPPORTED"), f"SUPPORTED without {missing}: {verdict}"


def test_simpledata_masking_roundtrip():
    data = make_clean(n=50)
    mask = np.asarray(data.obs["age"]) == "young"
    sub = data[mask]
    assert sub.n_obs == mask.sum()
    assert list(sub.var_names) == GENES


def test_gate6_stratified_catches_interaction_confound_on_replication():
    """GATE 6 applies the stratified QC diagnostic to the replication data too: a confound
    hidden in the sex x age interaction (identical biology) must not yield a replication
    claim. v0.3: the imbalance is reported (male stratum) and corrected, and the corrected,
    replicate-level effect decides — it is not detected, so the gate does not PASS.
    (v0.1.1 FAILed on QC parity alone.)"""
    data2 = add_mice(make_confounded())
    res = gate6_replication(NPR, data2, "age", ("young", "old"), within=["sex"],
                            replicate_col="mouse", estimand="composition", sesoi=0.1,
                            primary_effect=0.5)
    assert res.status != GateStatus.PASS
    assert res.detail["replication"] != "REPLICATED"
    assert res.detail["qc_status"] == "WARN"
    assert {r["stratum"]["sex"] for r in res.detail["qc_parity"]["flagged"]} == {"male"}
    assert res.detail["within"] == ["sex"]


def test_gate6_replicates_at_the_replicate_level_only():
    """A clean independent dataset with the real effect replicates — with mice and an
    estimand. v0.3: without a replicate unit the same data are INCONCLUSIVE (v0.1.1 passed
    on cells), and an effect detected with the opposite sign is NOT_REPLICATED."""
    data2 = add_mice(make_clean())
    kw = dict(within=["sex"], estimand="composition")
    res = gate6_replication(NPR, data2, "age", ("young", "old"), replicate_col="mouse",
                            primary_effect=0.5, **kw)
    assert res.status == GateStatus.PASS, res.message
    assert res.detail["replication"] == "REPLICATED"
    on_cells = gate6_replication(NPR, make_clean(), "age", ("young", "old"), primary_effect=0.5, **kw)
    assert on_cells.status == GateStatus.WARN and on_cells.detail["replication"] == "INCONCLUSIVE"
    flipped = gate6_replication(NPR, data2, "age", ("young", "old"), replicate_col="mouse",
                                primary_effect=-0.5, **kw)
    assert flipped.status == GateStatus.FAIL and flipped.detail["replication"] == "NOT_REPLICATED"


def test_whole_matrix_metric_probed_with_gene_subsample():
    """A whole-matrix metric (no bound gene pair) is probed with the gene_subsample
    perturbation via run_autopsy — previously that path was unreachable from any entrypoint."""
    data = make_clean()
    autopsy = run_autopsy(metrics.spectral_entropy, data,
                          group_col="age", groups=("young", "old"), within=["sex"])
    g0 = next(r for r in autopsy.results if r.gate == 0)
    assert "gene_subsample" in g0.detail["responses"]


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
