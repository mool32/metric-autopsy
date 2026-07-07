"""Synthetic tests: plant a known confound and assert the gates catch it.

The generators build single-cell-like count matrices where the *biology* is known by
construction, so we can assert that:
  - a QC confound hidden in a sex x age interaction is caught by GATE 1,
  - a metric whose effect is pure QC collapses under GATE 2 n_genes matching,
  - the sparsity-confounded metric (mi_3bin) fails GATE 0 while a robust metric passes,
  - a clean dataset with matched QC and real signal passes the full autopsy.
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
    data = make_confounded()
    # Pooled (no strata) also flags here, but the point is the sex-stratified check.
    res = gate1_qc_parity(data, "age", ("young", "old"), within=["sex"])
    assert res.status == GateStatus.FAIL
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


def test_gate0_flags_mi_but_not_robust_metric():
    data = make_clean()  # matched QC, so GATE 0 isolates intrinsic nuisance-sensitivity
    mi = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")
    robust = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
    assert gate0_independence(mi, data).status == GateStatus.FAIL
    assert gate0_independence(robust, data).status == GateStatus.PASS


def test_gate5_controls_positive_fires_negative_null():
    data = make_confounded()
    res = gate5_controls(
        metrics.mi_3bin, data,
        pos_pair=("Actb", "Gapdh"), neg_pair=("Gene0", "Gene1"), within=["sex"],
    )
    assert res.status == GateStatus.PASS, res.message


def test_full_autopsy_confounded_fails():
    data = make_confounded()
    m = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")
    autopsy = run_autopsy(
        m, data, group_col="age", groups=("young", "old"),
        within=["sex"], stop_on_first_fail=False,
    )
    assert autopsy.verdict.startswith("FAIL")


def test_full_autopsy_clean_passes():
    data = make_clean()
    m = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
    autopsy = run_autopsy(
        m, data, group_col="age", groups=("young", "old"), within=["sex"],
    )
    assert autopsy.verdict.startswith("PASS"), autopsy.verdict


def test_simpledata_masking_roundtrip():
    data = make_clean(n=50)
    mask = np.asarray(data.obs["age"]) == "young"
    sub = data[mask]
    assert sub.n_obs == mask.sum()
    assert list(sub.var_names) == GENES


def test_gate6_stratified_catches_interaction_confound_on_replication():
    """GATE 6 applies the stratified QC-parity check to the replication data too: a confound
    hidden in the sex x age interaction blocks a clean replication claim instead of passing
    because the pooled effect is diluted to ~0 (the trap that motivated the tool)."""
    data2 = make_confounded()
    m = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
    res = gate6_replication(m, data2, "age", ("young", "old"), within=["sex"])
    assert res.status == GateStatus.FAIL
    assert res.detail["within"] == ["sex"]


def test_gate6_passes_on_clean_replication():
    data2 = make_clean()
    m = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
    res = gate6_replication(m, data2, "age", ("young", "old"), within=["sex"])
    assert res.status == GateStatus.PASS, res.message


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
