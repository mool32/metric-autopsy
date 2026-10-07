"""Dev-set regression tests: one or more per exploratory probe (validation/probes/pNN_*.py).

Every test states the TRUTH, known by construction, and the verdict a correct validator
must give on it. Expected verdicts are frozen here, before any fix. A fix has to satisfy
them without editing them.

Tests marked ``xfail(strict=True)`` documented the known failures of v0.1.1. Every one was
flipped by the v0.3.0.dev0 rework and its marker replaced by a ``# fixed in`` comment that
names the original failure. These probes found the bugs, and the fixes were developed
against them: this is a DEVELOPMENT set with no confirmatory weight (see README.md).

The tests are written against the target (v0.3) API, which v0.1.1 does not have yet:

* ``run_autopsy(..., replicate_col=, signal_test=, seed=)`` returns an ``Autopsy`` with
  four independent assessments ``metric_validity``, ``design_adequacy``, ``effect``,
  ``replication``, each with ``.status`` / ``.reason`` / ``.detail``, plus ``.verdict``.
* metric_validity  ∈ PASS | FAIL | UNTESTED | DEGENERATE
* design_adequacy  ∈ ADEQUATE | CORRECTED | INSUFFICIENT_REPLICATION | UNIDENTIFIABLE
  (flags: UNDERPOWERED, PARAMETRIC_ONLY)
* effect           ∈ DETECTED | NO_DETECTABLE_EFFECT | INCONCLUSIVE | NOT_ESTIMABLE | NOT_RUN
* verdict starts with SUPPORTED | NOT SUPPORTED | NO DETECTABLE EFFECT | INCONCLUSIVE |
  UNIDENTIFIABLE | DEGENERATE METRIC
* prereg keys: ``estimand`` (composition | content), ``sesoi``, ``min_replicates``,
  ``judgment_pending`` (default True)
"""
from __future__ import annotations

import sys
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests"))  # 40-gene generators live in tests/test_gates.py

import probe_sim  # noqa: E402
from test_gates import _assemble, _block, make_clean  # noqa: E402

from metric_autopsy import (  # noqa: E402
    GateStatus, SimpleData, gate0_independence, gate1_qc_parity, gate5_controls, metrics,
    run_autopsy,
)
from metric_autopsy import gates as _gates  # noqa: E402
from metric_autopsy.cli import demo_data  # noqa: E402
from metric_autopsy.core import as_dense, unique_col_index  # noqa: E402

NPR = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
MI = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")
G2M_GENES = [f"M{i}" for i in range(probe_sim.N_G2M)]
COMPOSITION = {"estimand": "composition"}
RESOLVED = {"estimand": "composition", "judgment_pending": False}


# --------------------------------------------------------------------------- #
# data builders (deterministic)
# --------------------------------------------------------------------------- #
def _phases(p_s, p_g2m):
    def f(rng, n):
        u = rng.random(n)
        return np.where(u < p_s, "S", np.where(u < p_s + p_g2m, "G2M", "G1"))
    return f


def _sorted_cell_cycle(g2m_content=2.4, plates=4, cells=250, seed=1):
    """Cells sorted into G1 or G2M (FUCCI/Hoechst-style ground truth), several plates each."""
    content = (1.0, 1.4, g2m_content)
    blocks = [dict(n=cells, obs=dict(sorted_phase=ph, plate=f"{ph}-p{i}"), phase=ph, content=content)
              for ph in ("G1", "G2M") for i in range(plates)]
    return probe_sim.make(blocks, seed=seed)[0]


def _proliferation(content=(1.0, 1.6, 2.4), mice=3, cells=150, seed=2):
    """Young tissue: 35% cycling cells; old: 5%. Cycling cells carry more RNA (biology)."""
    blocks = []
    for sex in ("male", "female"):
        for age, phase in (("young", _phases(0.15, 0.20)), ("old", _phases(0.02, 0.03))):
            for i in range(mice):
                blocks.append(dict(n=cells, obs=dict(age=age, sex=sex, mouse=f"{age}-{sex}-{i}"),
                                   phase=phase, content=content))
    return probe_sim.make(blocks, seed=seed)[0]


def _with_mice(data, per_block, cols, seed=0):
    """Split every block of `cols` into `per_block` mice (a replicate column)."""
    rng = np.random.default_rng(seed)
    obs = data.obs.copy()
    key = obs[cols].astype(str).agg("-".join, axis=1)
    mouse = np.empty(len(obs), dtype=object)
    for k in pd.unique(key):
        idx = np.where(key.values == k)[0]
        mouse[idx] = [f"{k}-m{j}" for j in rng.permutation(np.arange(len(idx)) % per_block)]
    obs["mouse"] = mouse
    return SimpleData(data.X, obs, data.var_names)


def _xist_demo(mice=4):
    """Demo data (male-old capture 30%) plus Xist, expressed in female cells only."""
    d = demo_data()
    rng = np.random.default_rng(7)
    female = np.asarray(d.obs["sex"]) == "female"
    old_male = (~female) & (np.asarray(d.obs["age"]) == "old")
    xist = rng.poisson(np.where(female, 25.0, 0.0) * np.where(old_male, 0.30, 1.0)).astype(float)
    dx = SimpleData(np.column_stack([d.X, xist]), d.obs, list(d.var_names) + ["Xist"])
    return _with_mice(dx, mice, ["sex", "age"])


def _mean_lognorm_xist(data):
    X = as_dense(data.X)
    tot = X.sum(1)
    tot[tot == 0] = 1
    return float(np.mean(np.log1p(X[:, list(data.var_names).index("Xist")] / tot * 1e4)))


def _null_qc_strata(K, n, seed):
    """Identical log-normal n_genes distributions in both groups of every stratum."""
    rng = np.random.default_rng(seed)
    rows = []
    for k in range(K):
        for grp in ("young", "old"):
            ng = np.round(np.exp(rng.normal(np.log(2000), 0.45, n)))
            rows += [dict(age=grp, stratum=f"s{k}", n_genes_by_counts=v) for v in ng]
    obs = pd.DataFrame(rows)
    return SimpleData(np.ones((len(obs), 2)), obs, ["a", "b"])


def _null_controls(K, n, seed):
    """Controls behave identically everywhere: Actb-Gapdh coupled, Gene0-Gene1 independent."""
    rng = np.random.default_rng(seed)
    blocks = []
    for k in range(K):
        X, o = _block(rng, n, coupling=1.5, efficiency=1.0, age="young", sex="x")
        o["stratum"] = f"s{k}"
        blocks.append((X, o))
    return _assemble(blocks)


def _mice(n_mice, c_young=1.5, c_old=1.5, mouse_sd=0.0, cells=300, seed=0):
    """Mice nested in age; per-mouse coupling ~ N(c_age, mouse_sd). No QC difference."""
    rng = np.random.default_rng(seed)
    blocks = []
    for age, c in (("young", c_young), ("old", c_old)):
        for m in range(n_mice):
            cm = max(rng.normal(c, mouse_sd) if mouse_sd > 0 else c, 0.05)
            X, o = _block(rng, cells, coupling=cm, efficiency=1.0, age=age, sex="female")
            o["mouse"] = f"{age}{m}"
            blocks.append((X, o))
    return _assemble(blocks)


def _lognorm_pearson(data, *, gene_a, gene_b):
    """The most common gene-gene correlation in scRNA-seq: Pearson on log1p(CP10k), all cells."""
    X = as_dense(data.X)
    v = list(data.var_names)
    L = probe_sim.lognorm(X)
    a, b = L[:, unique_col_index(v, gene_a)], L[:, unique_col_index(v, gene_b)]
    if a.std() == 0 or b.std() == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def _coupled_pair(target_mean=3.0, n=1000, G=600, seed=11):
    """A truly coupled gene pair at a realistic depth; groups are random halves (depth-matched)."""
    rng = np.random.default_rng(5)
    mu = probe_sim.base_means(rng, G)
    cand = [i for i in np.argsort(mu) if i >= probe_sim.N_S + probe_sim.N_G2M]
    i = min(cand, key=lambda j: abs(mu[j] - target_mean))
    j = min((c for c in cand if c != i), key=lambda c: abs(mu[c] - target_mean * 1.1))
    X = probe_sim.simulate_cells(np.random.default_rng(seed), n, mu, phase="G1",
                                 coupled=(i, j), coupling=0.8)
    names = ([f"S{k}" for k in range(probe_sim.N_S)] + [f"M{k}" for k in range(probe_sim.N_G2M)]
             + [f"g{k}" for k in range(G - probe_sim.N_S - probe_sim.N_G2M)])
    obs = pd.DataFrame({"age": ["young"] * (n // 2) + ["old"] * (n - n // 2)})
    return SimpleData(X, obs, names), names[i], names[j]


def _demo_male_with_mice(mice=4):
    """Demo male stratum only: identical biology, old-male capture 30% (pure QC artifact)."""
    d = demo_data()
    males = d[np.asarray(d.obs["sex"]) == "male"]
    return _with_mice(males, mice, ["age"])


# --------------------------------------------------------------------------- #
# p01 — useless, constant and offset metrics
# --------------------------------------------------------------------------- #
# fixed in v0.3.0.dev0 — v0.1.1 failure (p01): a metric returning random numbers gets PASS
def test_p01_random_metric_is_not_certified():
    """Truth: the metric ignores the data. It must not be certified valid or 'supported'."""
    rng = np.random.default_rng(0)
    a = run_autopsy(lambda data: float(rng.normal()), make_clean(),
                    group_col="age", groups=("young", "old"), within=["sex"])
    assert a.metric_validity.status in ("UNTESTED", "FAIL")
    assert not a.verdict.startswith("SUPPORTED")


# fixed in v0.3.0.dev0 — v0.1.1 failure (p01): a constant metric is diagnosed as a QC artifact
def test_p01_constant_metric_is_reported_degenerate():
    """Truth: the metric is constant. The verdict must say so, not blame QC."""
    a = run_autopsy(lambda data: 0.0, make_clean(), group_col="age", groups=("young", "old"))
    assert a.metric_validity.status == "DEGENERATE"
    assert a.verdict.startswith("DEGENERATE METRIC")
    assert "QC artifact" not in a.to_markdown()


# fixed in v0.3.0.dev0 — v0.1.1 failure (p01): GATE 0 divides the shift by |baseline|
@pytest.mark.parametrize("C", [0.5, 1.0, 10.0])
def test_p01_gate0_is_invariant_to_adding_a_constant(C):
    """Truth: mi_3bin + C is the same metric as mi_3bin. GATE 0 must agree on both."""
    d = make_clean()
    base = gate0_independence(MI, d)
    shifted = gate0_independence(lambda data: MI(data) + C, d)
    assert shifted.status == base.status
    for kind, r in base.detail["responses"].items():
        assert shifted.detail["responses"][kind]["shift_std"] == pytest.approx(r["shift_std"], rel=1e-6)


@pytest.mark.parametrize("k", [100.0, 0.01, -1.0])
def test_p01_gate0_is_invariant_to_rescaling(k):
    """Truth: k·mi_3bin is the same metric. Rescaling or sign-flipping must leave GATE 0's
    status, every classification (bias / attenuation / none) and every |shift| (in scale
    units) unchanged. (Revised for the v0.3 scheme: v0.1.1 compared statuses only.)"""
    d = make_clean()
    base = gate0_independence(MI, d)
    scaled = gate0_independence(lambda data: k * MI(data), d)
    assert scaled.status == base.status
    for kind, r in base.detail["responses"].items():
        s_r = scaled.detail["responses"][kind]
        assert s_r["classification"] == r["classification"]
        assert abs(s_r["shift_std"]) == pytest.approx(abs(r["shift_std"]), rel=1e-6)


# --------------------------------------------------------------------------- #
# p02 — sorted G1 vs G2M: the positive control whose biology moves QC
# --------------------------------------------------------------------------- #
# fixed in v0.3.0.dev0 — v0.1.1 failure (p02): GATE 1 kills a real effect that GATE 2 retains at 99%
def test_p02_sorted_g2m_vs_g1_is_supported():
    """Truth: sorted G2M cells have a far higher G2M score than G1 cells; G2M cells also
    carry 2.4x RNA (biology). With a positive signal test, the claim must be SUPPORTED."""
    from metric_autopsy import injected_signal
    a = run_autopsy(probe_sim.mean_g2m_score, _sorted_cell_cycle(),
                    group_col="sorted_phase", groups=("G2M", "G1"), replicate_col="plate",
                    signal_test=injected_signal.module(G2M_GENES, fold=2.0, frac=0.3), prereg=RESOLVED)
    assert a.effect.status == "DETECTED" and a.effect.detail["effect"] > 0
    assert a.metric_validity.status == "PASS"
    assert a.design_adequacy.status in ("ADEQUATE", "CORRECTED")
    assert a.verdict.startswith("SUPPORTED")


# --------------------------------------------------------------------------- #
# p03 — proliferation decline: n_genes is downstream of the biology
# --------------------------------------------------------------------------- #
# fixed in v0.3.0.dev0 — v0.1.1 failure (p03): n_genes matching removes a real effect
def test_p03_proliferation_decline_is_supported_and_not_removed():
    """Truth: young tissue proliferates (35% cycling) and old does not (5%); cycling cells
    carry more RNA. For a composition estimand the effect must survive depth correction."""
    from metric_autopsy import injected_signal
    a = run_autopsy(probe_sim.mean_g2m_score, _proliferation(),
                    group_col="age", groups=("young", "old"), within=["sex"], replicate_col="mouse",
                    signal_test=injected_signal.module(G2M_GENES, fold=2.0, frac=0.3), prereg=RESOLVED)
    assert a.effect.status == "DETECTED"
    assert a.effect.detail["retained"] >= 0.8
    assert a.metric_validity.status == "PASS"
    assert a.verdict.startswith("SUPPORTED")


# --------------------------------------------------------------------------- #
# p04 — Xist on the demo data: the same-data positive control
# --------------------------------------------------------------------------- #
# fixed in v0.3.0.dev0 — v0.1.1 failure (p04): Xist female>male dies at GATE 1
def test_p04_xist_female_vs_male_is_supported_on_demo_data():
    """Truth: Xist is expressed in female cells only. The old stratum has a 2x QC gap
    (male-old capture 30%). The claim must be SUPPORTED, with the QC gap corrected."""
    from metric_autopsy import injected_signal
    a = run_autopsy(_mean_lognorm_xist, _xist_demo(), group_col="sex", groups=("female", "male"),
                    within=["age"], replicate_col="mouse",
                    signal_test=injected_signal.module(["Xist"], fold=2.0, frac=0.3), prereg=RESOLVED)
    assert a.effect.status == "DETECTED"
    assert a.effect.detail["retained"] >= 0.9
    assert a.design_adequacy.status == "CORRECTED"
    assert a.verdict.startswith("SUPPORTED")


# --------------------------------------------------------------------------- #
# p05 — null strata: the diagnostics must be calibrated
# --------------------------------------------------------------------------- #
# fixed in v0.3.0.dev0 — v0.1.1 failure (p05): GATE 1 flags 82% of null datasets (64 x 20 cells)
def test_p05_gate1_does_not_flag_null_strata():
    """Truth: no stratum has any QC difference. Flag rate must stay <= 20% (30 datasets)."""
    flagged = sum(
        gate1_qc_parity(_null_qc_strata(64, 20, s), "age", ("young", "old"),
                        within=["stratum"]).status.value in ("FAIL", "WARN", "STOP")
        for s in range(30))
    assert flagged <= 6


# fixed in v0.3.0.dev0 — v0.1.1 failure (p05): GATE 5 fails 70% of null datasets (4 x 400 cells)
def test_p05_gate5_calibrated_on_null_controls():
    """Truth: controls behave in every stratum. FAIL rate must stay <= 20% (20 datasets)."""
    fails = sum(
        gate5_controls(metrics.norm_pearson, _null_controls(4, 400, s), ("Actb", "Gapdh"),
                       ("Gene0", "Gene1"), within=["stratum"]).status == GateStatus.FAIL
        for s in range(20))
    assert fails <= 4


# --------------------------------------------------------------------------- #
# p06 — the negative control must be judged against an empirical null
# --------------------------------------------------------------------------- #
# fixed in v0.3.0.dev0 — v0.1.1 failure (p06): closure-inflated negative control fails GATE 5
def test_p06_negative_control_judged_against_empirical_null():
    """Truth: Gene0-Gene1 are independent in every stratum of the demo data; CP10k closure
    makes unrelated genes correlate (~+0.08). GATE 5 must pass and surface that null centre."""
    res = gate5_controls(metrics.norm_pearson, demo_data(), ("Actb", "Gapdh"), ("Gene0", "Gene1"),
                         within=["sex"])
    assert res.status == GateStatus.PASS
    assert any(abs(r["null_center"]) > 0.03 for r in res.detail["rows"])


# --------------------------------------------------------------------------- #
# p07 — pseudoreplication
# --------------------------------------------------------------------------- #
# fixed in v0.3.0.dev0 — v0.1.1 failure (p07): cell-level permutation certifies mouse noise
def test_p07_pseudoreplicated_null_yields_no_effect_claim():
    """Truth: no age effect; mice differ from each other (3 vs 3). Graded replicate rule
    (decided 2026-10-07): without a declared replicate unit there is no effect verdict; with
    3 mice per group the t interval is the test, flagged PARAMETRIC_ONLY. It is calibrated, so
    false detections stay near alpha (<= 3 of 20), and every detection carries the flag."""
    detected = 0
    for s in range(20):
        d = _mice(3, mouse_sd=0.35, seed=s)
        with_rep = run_autopsy(NPR, d, group_col="age", groups=("young", "old"),
                               replicate_col="mouse", prereg=COMPOSITION)
        without_rep = run_autopsy(NPR, d, group_col="age", groups=("young", "old"), prereg=COMPOSITION)
        assert without_rep.effect.status != "DETECTED"
        assert without_rep.design_adequacy.status == "INSUFFICIENT_REPLICATION"
        if with_rep.effect.status == "DETECTED":
            detected += 1
            assert "PARAMETRIC_ONLY" in with_rep.effect.flags
    assert detected <= 3


# fixed in v0.3.0.dev0 — v0.1.1 failure (p07): no replicate-level inference
def test_p07_replicate_level_inference_is_calibrated():
    """Truth: no age effect, 6 vs 6 mice with mouse-level variation. A calibrated test at
    alpha = 0.05 detects an effect in at most ~13% of 30 null datasets."""
    detected = sum(
        run_autopsy(NPR, _mice(6, mouse_sd=0.35, seed=s), group_col="age", groups=("young", "old"),
                    replicate_col="mouse", prereg=COMPOSITION).effect.status == "DETECTED"
        for s in range(30))
    assert detected <= 4


# --------------------------------------------------------------------------- #
# p08 — attenuation is a reliability property, not bias
# --------------------------------------------------------------------------- #
# fixed in v0.3.0.dev0 — v0.1.1 failure (p08): attenuation is scored as confounding
def test_p08_attenuation_is_not_scored_as_bias():
    """Truth: a truly coupled pair; groups are depth-matched random halves. Thinning shrinks
    the correlation toward its null (attenuation). That must be reported as attenuation,
    and must not FAIL the metric as confounded."""
    d, ga, gb = _coupled_pair()
    g0 = gate0_independence(partial(_lognorm_pearson, gene_a=ga, gene_b=gb), d, protect_genes=(ga, gb))
    assert g0.status != GateStatus.FAIL
    assert g0.detail["responses"]["depth_downsample"]["classification"] == "attenuation"


# --------------------------------------------------------------------------- #
# p09 — absence of evidence is not evidence of absence; verdicts must be stable
# --------------------------------------------------------------------------- #
# fixed in v0.3.0.dev0 — v0.1.1 failure (p09): 'the groups simply do not differ' without a test
def test_p09_equivalence_requires_a_sesoi():
    """Truth: a moderate real effect (coupling 1.5 vs 1.2) must never be called 'no
    difference'. On null data, 'no detectable effect' needs a pre-registered SESOI
    (equivalence test); without one the effect is INCONCLUSIVE."""
    for s in range(3):
        a = run_autopsy(NPR, _mice(6, c_old=1.2, seed=s), group_col="age", groups=("young", "old"),
                        replicate_col="mouse", prereg=COMPOSITION)
        assert a.effect.status != "NO_DETECTABLE_EFFECT"
        assert "do not differ" not in a.to_markdown()
    null = _mice(6, seed=0)
    no_sesoi = run_autopsy(NPR, null, group_col="age", groups=("young", "old"),
                           replicate_col="mouse", prereg=COMPOSITION)
    with_sesoi = run_autopsy(NPR, null, group_col="age", groups=("young", "old"),
                             replicate_col="mouse", prereg={**COMPOSITION, "sesoi": 0.15})
    assert no_sesoi.effect.status == "INCONCLUSIVE"
    assert with_sesoi.effect.status == "NO_DETECTABLE_EFFECT"


# fixed in v0.3.0.dev0 — v0.1.1 failure (p09): null floor from 20 permutations flips verdicts
def test_p09_effect_verdict_is_stable_across_tool_seeds():
    """Truth: one fixed dataset with a moderate real effect. The tool's own randomness
    (permutations, thinning draws) must not change the effect verdict."""
    d = _mice(6, c_old=1.2, seed=0)
    statuses = {run_autopsy(NPR, d, group_col="age", groups=("young", "old"), replicate_col="mouse",
                            prereg=COMPOSITION, seed=s).effect.status for s in range(4)}
    assert len(statuses) == 1


# --------------------------------------------------------------------------- #
# p10 — the headline sensitivity of mi_3bin is robust (lock-in, holds in v0.1.1)
# --------------------------------------------------------------------------- #
def test_p10_mi3bin_detection_sensitivity_seen_at_5pct_dropout(monkeypatch):
    """Truth: mi_3bin is driven by detection. Even a 5% extra dropout must move it by
    >25% of its value, reliably (z > 4). Keeps GATE 0's sensitivity from regressing, and
    checks that the sensitivity is now reported as attenuation rather than a FAIL."""
    orig = _gates._perturb

    def mild(data, kind, rng, protect=frozenset()):
        if kind != "extra_dropout":
            return orig(data, kind, rng, protect)
        X = data.X.copy()
        nz = np.argwhere(X > 0)
        pick = nz[rng.choice(len(nz), size=int(0.05 * len(nz)), replace=False)]
        X[pick[:, 0], pick[:, 1]] = 0.0
        return SimpleData(X, data.obs, data.var_names), {}

    monkeypatch.setattr(_gates, "_perturb", mild)
    g0 = _gates.gate0_independence(MI, demo_data())
    r = g0.detail["responses"]["extra_dropout"]
    rel = abs(r["mean_perturbed"] - g0.detail["baseline"]) / abs(g0.detail["baseline"])
    assert rel > 0.25 and r["z"] > 4
    # v0.3 meaning (decided 2026-10-07): this sensitivity is attenuation — reported with its
    # size, not failed; a differential version is caught by the effect field after equalization
    assert r["classification"] == "attenuation" and r["signal_loss"] > 0.25
    assert g0.status != GateStatus.FAIL


# --------------------------------------------------------------------------- #
# p11 — a pure depth artifact must be explained by depth
# --------------------------------------------------------------------------- #
# fixed in v0.3.0.dev0 — v0.1.1 failure (p11): no depth equalization
def test_p11_depth_artifact_is_explained_by_depth():
    """Truth: biology is identical; old-male capture is 30%, so mi_3bin is lower in old
    males. At equal depth the difference must vanish, and the claim is NOT SUPPORTED."""
    a = run_autopsy(MI, _demo_male_with_mice(), group_col="age", groups=("young", "old"),
                    replicate_col="mouse", prereg=COMPOSITION)
    assert a.effect.detail["explained_by_depth"] is True
    assert a.effect.status != "DETECTED"
    assert abs(a.effect.detail["retained"]) <= 0.25
    assert a.verdict.startswith("NOT SUPPORTED")


# --------------------------------------------------------------------------- #
# p12 — atlas-scale memory: warn before densifying
# --------------------------------------------------------------------------- #
# fixed in v0.3.0.dev0 — v0.1.1 failure (p12): silent densification of sparse X
def test_p12_warns_before_densifying_a_large_sparse_matrix(monkeypatch):
    """Truth: densifying this sparse matrix exceeds the configured budget. The engine must
    warn (naming the dense size) before it allocates."""
    anndata = pytest.importorskip("anndata")
    sparse = pytest.importorskip("scipy.sparse")
    monkeypatch.setenv("METRIC_AUTOPSY_DENSE_WARN_GB", "0.001")
    X = sparse.random(2000, 500, density=0.02, format="csr", random_state=0)
    obs = pd.DataFrame({"age": ["young", "old"] * 1000}, index=[f"c{i}" for i in range(2000)])
    var = pd.DataFrame(index=["mt-Co1"] + [f"g{i}" for i in range(499)])
    ad = anndata.AnnData(X=X, obs=obs, var=var)
    with pytest.warns(UserWarning, match="dense"):
        gate1_qc_parity(ad, "age", ("young", "old"))
