"""Regression tests for the code-testable subset of the 24 audit findings (18 of 24 are locked
here; #8, #9, #13, #14, #16, #23 are covered by code review / other suites, not a dedicated test
here). Each test name cites the finding it locks.

Where a test builds data with an explicit ``n_genes_by_counts`` / ``signal`` obs column, it
uses a stub metric that reads ``obs`` directly — this makes the gate's decision logic
deterministic and independent of any particular biological metric, which is exactly what a
metric-agnostic engine should allow.
"""
from functools import partial

import numpy as np
import pandas as pd
import pytest

from metric_autopsy import (
    SimpleData, GateStatus, metrics, run_autopsy,
    gate0_independence, gate1_qc_parity, gate2_ngenes_matching, gate5_controls,
)
from metric_autopsy.qc import _ks_2samp, _dist_overlap
from metric_autopsy.report import _metric_name, _cell
from metric_autopsy import cli

from test_gates import make_confounded, make_clean, GENES


# --- helpers ---------------------------------------------------------------
def _stub_data(rows, n_genes_col=True):
    """rows: list of dicts with age, n_genes, signal. Builds SimpleData whose obs drives
    a stub metric; X is filler so total_counts is defined."""
    obs = pd.DataFrame(rows)
    X = np.ones((len(obs), 4), dtype=float)
    if n_genes_col:
        obs = obs.rename(columns={"n_genes": "n_genes_by_counts"})
    return SimpleData(X, obs, ["g0", "g1", "g2", "g3"])


def _signal_metric(d):
    return float(np.mean(np.asarray(d.obs["signal"], dtype=float)))


# --- #1 / #2 : _dist_overlap and gate1 on comparable/identical groups ------
def test_dist_overlap_peaked_identical_is_one():  # #2
    peaked = np.full(200, 8.0)
    assert _dist_overlap(peaked, peaked) == 1.0


def test_dist_overlap_nested_scores_high():  # #1 / #18
    wide = np.arange(0, 100).astype(float)
    narrow = np.arange(45, 55).astype(float)  # fully inside wide, same location
    assert _dist_overlap(wide, narrow) > 0.9


def test_gate1_passes_identical_groups():  # #1 / #2
    rng = np.random.default_rng(0)
    X = rng.poisson(6, size=(600, 8)).astype(float)
    obs = pd.DataFrame({"age": ["young"] * 300 + ["old"] * 300})
    d = SimpleData(X, obs, [f"g{i}" for i in range(8)])
    assert gate1_qc_parity(d, "age", ("young", "old")).status == GateStatus.PASS


# --- #3 : match_by_ngenes on peaked identical distributions ----------------
def test_gate2_not_stop_on_peaked_identical():
    rng = np.random.default_rng(0)
    X = rng.poisson(6, size=(400, 4)).astype(float)
    d = SimpleData(X, pd.DataFrame({"age": ["young"] * 200 + ["old"] * 200}),
                   ["Smad3", "Col1a1", "Actb", "Gapdh"])
    m = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
    assert gate2_ngenes_matching(m, d, "age", ("young", "old")).status != GateStatus.STOP


# --- #5 : gate1 STOP when groups confounded with the stratifier ------------
def test_gate1_stop_when_no_common_stratum():
    d = make_clean()
    # make sex perfectly predict age -> no stratum has both groups
    d.obs["sex"] = np.where(np.asarray(d.obs["age"]) == "young", "male", "female")
    assert gate1_qc_parity(d, "age", ("young", "old"), within=["sex"]).status == GateStatus.STOP


# --- #6 / #7 : GATE 0 noise-vs-bias and near-zero baseline -----------------
def test_gate0_constant_metric_is_degenerate_not_confounded():
    # #7 locked in that a zero baseline is not a false 'confounded' FAIL. v0.3 keeps that and
    # names the real problem: a metric that never varies is DEGENERATE (v0.1.1 said PASS).
    d = make_clean()
    res = gate0_independence(lambda data: 0.0, d)
    assert res.status == GateStatus.DEGENERATE
    assert not any(r.get("confounded") for r in res.detail["responses"].values())


def test_gate0_scale_invariant_response_not_flagged():
    # library_scale (per-cell scaling) is divided out by norm_pearson -> must not flag.
    d = make_clean()
    m = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
    resp = gate0_independence(m, d).detail["responses"]
    assert resp["library_scale"]["confounded"] is False


# --- #10 : GATE 2 amplification / collapse / no-effect branches ------------
def test_gate2_amplification_fails():
    rows = ([{"age": "young", "n_genes": 100, "signal": 1.0}] * 200
            + [{"age": "young", "n_genes": 10, "signal": -0.9}] * 200
            + [{"age": "old", "n_genes": 100, "signal": 0.0}] * 400)
    d = _stub_data(rows)
    res = gate2_ngenes_matching(_signal_metric, d, "age", ("young", "old"), min_effect=0.0)
    assert res.status == GateStatus.FAIL
    assert "AMPLIF" in res.message.upper()


def test_gate2_collapse_fails():
    rows = ([{"age": "young", "n_genes": 100, "signal": 0.1}] * 200
            + [{"age": "young", "n_genes": 10, "signal": 1.9}] * 200
            + [{"age": "old", "n_genes": 100, "signal": 0.0}] * 400)
    d = _stub_data(rows)
    res = gate2_ngenes_matching(_signal_metric, d, "age", ("young", "old"), min_effect=0.0)
    assert res.status == GateStatus.FAIL
    assert "collapse" in res.message.lower()


def test_gate2_no_effect_is_not_a_false_pass_claim():
    rng = np.random.default_rng(3)
    rows = [{"age": a, "n_genes": 100, "signal": float(rng.normal(0, 0.01))}
            for a in (["young"] * 300 + ["old"] * 300)]
    d = _stub_data(rows)
    res = gate2_ngenes_matching(_signal_metric, d, "age", ("young", "old"))
    assert res.status == GateStatus.PASS
    assert "nothing to preserve" in res.message


# --- #4 / #11 / #17 : GATE 5 controls --------------------------------------
def test_gate5_skip_when_all_strata_too_small():
    d = make_clean(n=3)  # 12 cells, split by sex -> no stratum >= 10
    res = gate5_controls(metrics.mi_3bin, d, ("Actb", "Gapdh"), ("Gene0", "Gene1"),
                         within=["sex"], min_cells=10)
    assert res.status == GateStatus.SKIP


def test_gate5_fails_broken_negative_control():
    d = make_clean()
    # negative control is actually the (coupled) housekeeping pair -> must FAIL
    res = gate5_controls(metrics.mi_3bin, d, ("Actb", "Gapdh"), ("Actb", "Gapdh"))
    assert res.status == GateStatus.FAIL


# --- #15 : KS asymptotic fallback ------------------------------------------
def test_ks_fallback_identical_is_one():
    assert _ks_2samp(np.arange(50.0), np.arange(50.0)) == 1.0


def test_ks_fallback_disjoint_is_small():
    p = _ks_2samp(np.arange(0.0, 50.0), np.arange(1000.0, 1050.0))
    assert p < 0.01


# --- #19 : markdown escaping ------------------------------------------------
def test_markdown_cell_escapes_pipe_and_newline():
    assert _cell("a|b\nc") == "a\\|b<br>c"


# --- #21 : stop_on_first_fail honored for pairwise gates -------------------
def test_stop_halts_before_gate6_after_gate5_fail():
    d = make_clean()
    m = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
    autopsy = run_autopsy(
        m, d, group_col="age", groups=("young", "old"), within=["sex"],
        gene_pair=("Smad3", "Col1a1"),
        pair_metric=metrics.mi_3bin, pos_pair=("Actb", "Gapdh"), neg_pair=("Actb", "Gapdh"),
        data2=d, stop_on_first_fail=True,
    )
    gates_run = {r.gate for r in autopsy.results}
    assert 5 in gates_run and 6 not in gates_run  # GATE 5 blocked -> GATE 6 never ran


# --- #22 : duplicate var_names rejected ------------------------------------
def test_duplicate_var_names_rejected():
    with pytest.raises(ValueError):
        SimpleData(np.ones((3, 3)), pd.DataFrame({"a": [1, 2, 3]}), ["g", "g", "h"])


# --- #24 : partial-bound metric keeps its name -----------------------------
def test_metric_name_unwraps_partial():
    assert _metric_name(partial(metrics.mi_3bin, gene_a="x", gene_b="y")) == "mi_3bin"


def _verdict_line(markdown: str) -> str:
    """The verdict sentence printed under '## Verdict' (bold markers stripped)."""
    lines = markdown.splitlines()
    i = lines.index("## Verdict")
    return next(line for line in lines[i + 1:] if line.strip()).strip("* ")


VERDICT_PREFIXES = ("SUPPORTED", "NOT SUPPORTED", "NO DETECTABLE EFFECT", "INCONCLUSIVE",
                    "UNIDENTIFIABLE", "DEGENERATE METRIC")


# --- #12 : spectral_entropy (no-genes metric) does not crash via CLI -------
def test_spectral_entropy_cli_no_crash(capsys):
    cli.main(["--demo", "--metric", "spectral_entropy", "--no-stop", "--no-log"])
    out = capsys.readouterr().out
    assert "spectral_entropy" in out
    assert _verdict_line(out).startswith(VERDICT_PREFIXES)


# --- #20 : --resolve-judgment is plumbed through, and CLI == API -----------
def test_cli_verdict_matches_api_with_resolve_judgment(capsys):
    """The CLI must print exactly the verdict the Python API gives for the same inputs."""
    cli.main(["--demo", "--metric", "norm_pearson", "--resolve-judgment", "--no-stop", "--no-log"])
    cli_verdict = _verdict_line(capsys.readouterr().out)
    api = run_autopsy(
        partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1"), cli.demo_data(),
        group_col="age", groups=("young", "old"), within=["sex"], gene_pair=("Smad3", "Col1a1"),
        pair_metric=metrics.norm_pearson, pos_pair=("Actb", "Gapdh"), neg_pair=("Gene0", "Gene1"),
        replicate_col="mouse", prereg={"estimand": "composition", "judgment_pending": False},
        stop_on_first_fail=False, log_path="off",
    )
    assert cli_verdict == api.verdict


def test_resolving_judgment_is_necessary_but_not_sufficient_for_supported():
    """v0.3 meaning: resolving the judgment gates is the last step to SUPPORTED, never a way
    around the other three requirements. With a valid metric (positive control), a declared
    estimand and a replicate-level DETECTED effect, pending judgment holds the verdict at
    INCONCLUSIVE and resolving it — only that — gives the provisional SUPPORTED. Without
    replicates or without a control, resolving it changes nothing. (v0.1.1: resolving the
    judgment turned INCONCLUSIVE into PASS on a cell-level comparison with no control.)"""
    from test_gates import CONTROLS, NPR, add_mice
    d = add_mice(make_clean())
    kw = dict(group_col="age", groups=("young", "old"), within=["sex"], gene_pair=("Smad3", "Col1a1"))
    full = dict(kw, replicate_col="mouse", **CONTROLS)
    pending = run_autopsy(NPR, d, prereg={"estimand": "composition"}, **full)
    resolved = run_autopsy(NPR, d, prereg={"estimand": "composition", "judgment_pending": False}, **full)
    assert pending.verdict.startswith("INCONCLUSIVE") and "judgment" in pending.verdict
    assert resolved.verdict == "SUPPORTED (provisional until replicated)"
    assert ({k: v.status for k, v in pending.fields().items()}
            == {k: v.status for k, v in resolved.fields().items()})  # only the judgment differs
    no_reps = run_autopsy(NPR, d, prereg={"estimand": "composition", "judgment_pending": False},
                          **dict(kw, **CONTROLS))
    no_control = run_autopsy(NPR, d, prereg={"estimand": "composition", "judgment_pending": False},
                             **dict(kw, replicate_col="mouse"))
    for a in (no_reps, no_control):
        assert a.verdict.startswith("INCONCLUSIVE"), a.verdict


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
