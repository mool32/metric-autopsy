"""v0.4 change 11 (validation/plans/v0.4_spec.md, section 11; coded on the owner's decision of 2026-11-05): the module
probe in the CLI and the MCP server.

``--inject-signal module --module-genes ... [--module-fold F]`` (the MCP tool's ``inject_signal="module"``,
``module_genes``, ``module_fold``) runs GATE 4 with ``injected_signal.module``; with a SESOI on the construct's scale
(a change of the module's natural-log expression) GATE 4 doses it to exp(SESOI)-fold. The probe is built by one
function for every front end (``injected_signal.named``), and the API, the CLI and the MCP tool give one verdict.

The set's coupling probe (change 9) the same way, by the owner's decision of 2026-11-05: ``--inject-signal
coupling_set --set-genes ... [--couple half]`` (``inject_signal="coupling_set"``, ``set_genes``, ``couple``); with a
SESOI on the construct's scale GATE 4 doses it to that change of the mean correlation over the pairs it couples."""
from __future__ import annotations

import json
from functools import partial

import pandas as pd
import pytest

from metric_autopsy import cli, injected_signal, mcp_server, metrics, run_autopsy

MODULE = ["Smad3", "Col1a1"]
PREREG = dict(estimand="composition", direction="decrease", sesoi=0.2)


def _verdict_line(markdown: str) -> str:
    lines = markdown.splitlines()
    i = lines.index("## Verdict")
    return next(line for line in lines[i + 1:] if line.strip()).strip("* ")


def _gate4(report: dict) -> dict:
    return next(g for g in report["gates"] if g["gate"] == 4)


def test_a_front_end_names_its_probe_and_a_module_needs_its_genes():
    assert injected_signal.named(None) is None and injected_signal.named("") is None
    c = injected_signal.named("coupling", gene_pair=("A", "B"))
    assert c.kind == "coupling" and c.genes == ("A", "B")
    assert injected_signal.named("coupling") is None  # no pair, no coupling (the CLI's behaviour before v0.4)
    m = injected_signal.named("module", module_genes=["A", "B", "C"], module_fold=3.0)
    assert m.kind == "module" and m.genes == ("A", "B", "C") and m.fold == 3.0 and m.direction == "decrease"
    for bad in (dict(kind="module"), dict(kind="module", module_genes=[]), dict(kind="spike")):
        with pytest.raises(ValueError):
            injected_signal.named(bad.pop("kind"), **bad)


def test_the_cli_runs_the_module_probe_dosed_to_the_sesoi_with_the_apis_verdict(capsys, tmp_path):
    out = tmp_path / "cli.json"
    cli.main(["--demo", "--inject-signal", "module", "--module-genes", *MODULE, "--sesoi", "0.2",
              "--direction", "decrease", "--no-log", "--json", str(out)])
    md = capsys.readouterr().out
    rep = json.loads(out.read_text())
    g4 = _gate4(rep)
    assert g4["detail"]["declared_injection"].startswith("module(2 genes")
    assert g4["detail"]["dose"]["kind"] == "module" and g4["detail"]["dose"]["construct_change"] == 0.2
    api = run_autopsy(partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1"), cli.demo_data(), group_col="age",
                      groups=("young", "old"), within=["sex"], gene_pair=("Smad3", "Col1a1"), replicate_col="mouse",
                      pair_metric=metrics.mi_3bin, pos_pair=("Actb", "Gapdh"), neg_pair=("Gene0", "Gene1"),
                      signal_test=injected_signal.module(MODULE), prereg=dict(PREREG), stop_on_first_fail=True,
                      log_path="off")
    assert _verdict_line(md) == rep["verdict"] == api.verdict


def test_the_cli_refuses_a_module_probe_without_its_genes(capsys):
    with pytest.raises(SystemExit) as e:
        cli.main(["--demo", "--inject-signal", "module", "--no-log"])
    assert "module probe needs its genes" in str(e.value)


def test_the_mcp_tool_runs_the_module_probe_with_the_clis_verdict(capsys, tmp_path):
    anndata = pytest.importorskip("anndata")
    sparse = pytest.importorskip("scipy.sparse")
    sd = cli.demo_data()
    path = tmp_path / "demo.h5ad"
    ad = anndata.AnnData(X=sparse.csr_matrix(sd.X), obs=sd.obs.copy(), var=pd.DataFrame(index=list(sd.var_names)))
    ad.obs_names = [f"cell_{i}" for i in range(sd.n_obs)]
    ad.write_h5ad(path)
    common = ["--metric", "norm_pearson", "--gene-a", "Smad3", "--gene-b", "Col1a1", "--group-col", "age",
              "--groups", "young", "old", "--within", "sex", "--replicate-col", "mouse", "--estimand", "composition",
              "--direction", "decrease", "--sesoi", "0.2", "--inject-signal", "module", "--module-genes", *MODULE,
              "--module-fold", "3", "--no-stop", "--no-log"]
    cli.main(["--h5ad", str(path), *common, "--json", str(tmp_path / "cli.json")])
    cli_md = capsys.readouterr().out
    mcp_md = mcp_server.autopsy_report(
        str(path), metric="norm_pearson", gene_a="Smad3", gene_b="Col1a1", group_col="age", groups=["young", "old"],
        within=["sex"], replicate_col="mouse", estimand="composition", direction="decrease", sesoi=0.2,
        inject_signal="module", module_genes=MODULE, module_fold=3.0, stop_on_first_fail=False,
        json_path=str(tmp_path / "mcp.json"), log_path="off")
    assert _verdict_line(cli_md) == _verdict_line(mcp_md)
    reports = [json.loads((tmp_path / f).read_text()) for f in ("cli.json", "mcp.json")]
    for rep in reports:
        assert _gate4(rep)["detail"]["declared_injection"] == "module(2 genes, 1/3 by thinning, frac=1)"
        assert _gate4(rep)["detail"]["dose"]["kind"] == "module"
    assert reports[0]["provenance"]["data_sha256"] == reports[1]["provenance"]["data_sha256"]


def test_the_mcp_tool_declares_the_module_probes_parameters():
    pytest.importorskip("mcp")
    import asyncio
    tools = {t.name: t for t in asyncio.run(mcp_server.build_server().list_tools())}
    t = tools["autopsy_report"]
    props = (getattr(t, "inputSchema", None) or getattr(t, "input_schema"))["properties"]
    assert {"inject_signal", "module_genes", "module_fold", "set_genes", "couple"} <= set(props)


# --------------------------------------------------------------------------- #
# the set's probe (change 9) in the front ends, by the owner's decision of 2026-11-05
# --------------------------------------------------------------------------- #
SET = ["Smad3", "Col1a1", "Actb", "Gapdh"]
SET_PREREG = dict(estimand="composition", direction="decrease", sesoi=0.1)


def test_a_front_end_names_the_sets_probe_and_it_needs_its_genes():
    assert "coupling_set" in injected_signal.FRONT_END_PROBES
    s = injected_signal.named("coupling_set", set_genes=SET)
    assert s.kind == "coupling_set" and s.genes == tuple(SET) and s.couple == "all" and s.direction == "increase"
    assert injected_signal.named("coupling_set", set_genes=SET, couple="half").couple == "half"
    for bad in (dict(), dict(set_genes=[]), dict(set_genes=SET[:3], couple="half"), dict(set_genes=SET, couple="some")):
        with pytest.raises(ValueError):
            injected_signal.named("coupling_set", **bad)


def test_the_cli_runs_the_sets_probe_dosed_to_the_sesoi_with_the_apis_verdict(capsys, tmp_path):
    out = tmp_path / "cli.json"
    cli.main(["--demo", "--inject-signal", "coupling_set", "--set-genes", *SET, "--sesoi", "0.1",
              "--direction", "decrease", "--no-log", "--json", str(out)])
    md = capsys.readouterr().out
    rep = json.loads(out.read_text())
    g4 = _gate4(rep)
    assert g4["detail"]["declared_injection"] == "coupling_set(4 genes, couple=all, strength=1)"
    dose = g4["detail"]["dose"]
    assert dose["kind"] == "coupling_set" and dose["n_pairs"] == 6 and dose["construct_change"] == 0.1
    api = run_autopsy(partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1"), cli.demo_data(), group_col="age",
                      groups=("young", "old"), within=["sex"], gene_pair=("Smad3", "Col1a1"), replicate_col="mouse",
                      pair_metric=metrics.mi_3bin, pos_pair=("Actb", "Gapdh"), neg_pair=("Gene0", "Gene1"),
                      signal_test=injected_signal.coupling_set(SET), prereg=dict(SET_PREREG), stop_on_first_fail=True,
                      log_path="off")
    assert _verdict_line(md) == rep["verdict"] == api.verdict
    assert _gate4(json.loads(api.to_json()))["detail"]["dose"] == dose


def test_the_cli_refuses_a_sets_probe_without_its_genes_or_a_half_of_fewer_than_two(capsys):
    with pytest.raises(SystemExit) as e:
        cli.main(["--demo", "--inject-signal", "coupling_set", "--no-log"])
    assert "set's probe needs its genes" in str(e.value)
    with pytest.raises(SystemExit) as e:
        cli.main(["--demo", "--inject-signal", "coupling_set", "--set-genes", "Smad3", "Col1a1", "Actb",
                  "--couple", "half", "--no-log"])
    assert "at least 4 genes" in str(e.value)


def test_the_mcp_tool_runs_the_sets_probe_with_the_clis_verdict(capsys, tmp_path):
    anndata = pytest.importorskip("anndata")
    sparse = pytest.importorskip("scipy.sparse")
    sd = cli.demo_data()
    path = tmp_path / "demo.h5ad"
    ad = anndata.AnnData(X=sparse.csr_matrix(sd.X), obs=sd.obs.copy(), var=pd.DataFrame(index=list(sd.var_names)))
    ad.obs_names = [f"cell_{i}" for i in range(sd.n_obs)]
    ad.write_h5ad(path)
    common = ["--metric", "norm_pearson", "--gene-a", "Smad3", "--gene-b", "Col1a1", "--group-col", "age",
              "--groups", "young", "old", "--within", "sex", "--replicate-col", "mouse", "--estimand", "composition",
              "--direction", "decrease", "--sesoi", "0.1", "--inject-signal", "coupling_set", "--set-genes", *SET,
              "--couple", "half", "--no-stop", "--no-log"]
    cli.main(["--h5ad", str(path), *common, "--json", str(tmp_path / "cli.json")])
    cli_md = capsys.readouterr().out
    mcp_md = mcp_server.autopsy_report(
        str(path), metric="norm_pearson", gene_a="Smad3", gene_b="Col1a1", group_col="age", groups=["young", "old"],
        within=["sex"], replicate_col="mouse", estimand="composition", direction="decrease", sesoi=0.1,
        inject_signal="coupling_set", set_genes=SET, couple="half", stop_on_first_fail=False,
        json_path=str(tmp_path / "mcp.json"), log_path="off")
    assert _verdict_line(cli_md) == _verdict_line(mcp_md)
    reports = [json.loads((tmp_path / f).read_text()) for f in ("cli.json", "mcp.json")]
    for rep in reports:
        g4 = _gate4(rep)
        assert g4["detail"]["declared_injection"] == "coupling_set(4 genes, couple=half, strength=1)"
        assert g4["detail"]["dose"]["kind"] == "coupling_set" and g4["detail"]["dose"]["couple"] == "half"
    assert reports[0]["provenance"]["data_sha256"] == reports[1]["provenance"]["data_sha256"]
    assert _gate4(reports[0])["detail"]["dose"] == _gate4(reports[1])["detail"]["dose"]
