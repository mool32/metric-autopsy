"""MCP server — expose the gates to any MCP-capable agent (Claude Desktop, Cursor, Cline…).

This is the "use it through an agent" front door: a researcher points their agent at an
`.h5ad` and the agent calls `autopsy` as a tool, getting the gate-by-gate verdict back. It
complements the Claude Code skill (`SKILL.md`) and the CLI.

Run it (after `pip install "metric-autopsy[mcp]"`):

    python -m metric_autopsy.mcp_server          # stdio transport
    # or the console script:
    metric-autopsy-mcp

Then register it in your agent's MCP config, e.g. Claude Desktop `claude_desktop_config.json`:

    {"mcpServers": {"metric-autopsy": {"command": "metric-autopsy-mcp"}}}

The core report functions below are plain and dependency-light (anndata only); the FastMCP
wrappers are added lazily in `build_server()` so importing this module never requires `mcp`.
"""
from __future__ import annotations

from functools import partial
from typing import Sequence

from . import metrics
from . import qc as _qc
from .report import run_autopsy


def _load(path: str):
    try:
        import anndata
    except ImportError as e:  # pragma: no cover
        raise RuntimeError("reading .h5ad needs anndata — `pip install \"metric-autopsy[mcp]\"`") from e
    return anndata.read_h5ad(path)


def _bind(metric_name: str, gene_a: str | None, gene_b: str | None):
    fn = getattr(metrics, metric_name, None)
    if fn is None:
        raise ValueError(f"unknown metric {metric_name!r}; try list_metrics()")
    pairwise = metric_name != "spectral_entropy"
    if pairwise and gene_a and gene_b:
        return partial(fn, gene_a=gene_a, gene_b=gene_b), fn, pairwise
    return fn, fn, pairwise


# --------------------------------------------------------------------------- #
# Plain report functions (testable without `mcp`)
# --------------------------------------------------------------------------- #
def autopsy_report(
    h5ad_path: str,
    metric: str = "mi_3bin",
    gene_a: str | None = None,
    gene_b: str | None = None,
    group_col: str = "age",
    groups: Sequence[str] = ("young", "old"),
    within: Sequence[str] = (),
    pos_pair: Sequence[str] | None = None,
    neg_pair: Sequence[str] | None = None,
    data2_path: str | None = None,
    resolve_judgment: bool = False,
    stop_on_first_fail: bool = True,
) -> str:
    """Run the full metric-autopsy gate sequence on an AnnData `.h5ad` and return a Markdown
    report (gate-by-gate table + a verdict decided by the first blocking gate).

    metric: one of the reference metrics (see list_metrics) or `spectral_entropy` (no genes).
    groups: the two levels of `group_col` to compare, e.g. ["young", "old"].
    within: factorial obs columns to stratify QC/controls by, e.g. ["sex", "tissue"].
    pos_pair/neg_pair: control gene pairs for GATE 5. data2_path: second .h5ad for GATE 6.
    resolve_judgment: mark judgment gates 4 & 7 resolved so a provisional PASS is reachable.
    """
    data = _load(h5ad_path)
    metric_fn, raw_fn, pairwise = _bind(metric, gene_a, gene_b)
    gene_pair = (gene_a, gene_b) if (pairwise and gene_a and gene_b) else None
    do_controls = bool(pairwise and pos_pair and neg_pair)
    data2 = _load(data2_path) if data2_path else None

    autopsy = run_autopsy(
        metric_fn, data,
        group_col=group_col, groups=tuple(groups),
        gene_pair=gene_pair, within=list(within),
        pair_metric=raw_fn if do_controls else None,
        pos_pair=tuple(pos_pair) if do_controls else None,
        neg_pair=tuple(neg_pair) if do_controls else None,
        data2=data2,
        prereg={"judgment_pending": not resolve_judgment},
        stop_on_first_fail=stop_on_first_fail,
    )
    autopsy.metric_name = metric
    return autopsy.to_markdown()


def qc_parity_report(
    h5ad_path: str,
    group_col: str = "age",
    groups: Sequence[str] = ("young", "old"),
    within: Sequence[str] = (),
) -> str:
    """GATE 1 only: per-stratum QC comparison of two groups across factorial `within`
    columns. Returns a Markdown table of median n_genes ratio + overlap per stratum, with
    flagged strata marked. The cheapest, highest-yield check — run it first."""
    data = _load(h5ad_path)
    res = _qc.qc_parity(data, group_col, tuple(groups), within=list(within))
    lines = [f"# QC parity — {group_col} {tuple(groups)} within {list(within) or 'pooled'}", "",
             "| Stratum | n(A) | n(B) | median n_genes A | median n_genes B | ratio | overlap | flagged |",
             "|---|---|---|---|---|---|---|---|"]
    for r in res["table"]:
        lines.append(f"| {r['stratum']} | {r['n_a']} | {r['n_b']} | {r['median_n_genes_a']:.0f} | "
                     f"{r['median_n_genes_b']:.0f} | {r['n_genes_ratio']:.2f} | {r['overlap']:.2f} | "
                     f"{'⚠️ YES' if r['flagged'] else 'no'} |")
    verdict = ("FAIL — some strata are QC-disparate; matching required (GATE 2)"
               if res["flagged"] else "PASS — all strata within QC tolerance")
    lines += ["", f"**{verdict}** (worst ratio {res['worst_ratio']:.2f}× across "
              f"{res['n_compared']} strata)."]
    return "\n".join(lines)


def list_metrics() -> str:
    """List the reference metrics available by name, and how to bring your own."""
    rows = [
        ("mi_3bin", "MI with a zero/low/high discretization — the confounded antihero (pairwise)"),
        ("pearson", "Pearson correlation of two genes across cells (pairwise)"),
        ("codetected_spearman", "Spearman on co-detected cells — partial fix (pairwise)"),
        ("norm_pearson", "CP10k+log1p Pearson on co-detected cells — robust reference (pairwise)"),
        ("spectral_entropy", "Shannon entropy of the gene-gene correlation eigenspectrum (whole-matrix)"),
    ]
    out = ["# Reference metrics", "", "| name | description |", "|---|---|"]
    out += [f"| `{n}` | {d} |" for n, d in rows]
    out += ["", "**Bring your own:** from Python, pass any `metric(data) -> float` to "
            "`metric_autopsy.run_autopsy`. The MCP tools cover the reference set above by name."]
    return "\n".join(out)


def demo_report(stop_on_first_fail: bool = False) -> str:
    """Run the reference `mi_3bin` autopsy on the bundled synthetic confound (no data file
    needed) — the fastest way for an agent to see what an autopsy looks like."""
    from .cli import demo_data
    data = demo_data()
    m = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")
    autopsy = run_autopsy(
        m, data, group_col="age", groups=("young", "old"), within=["sex"],
        gene_pair=("Smad3", "Col1a1"),
        pair_metric=metrics.mi_3bin, pos_pair=("Actb", "Gapdh"), neg_pair=("Gene0", "Gene1"),
        stop_on_first_fail=stop_on_first_fail,
    )
    autopsy.metric_name = "mi_3bin (demo)"
    return autopsy.to_markdown()


# --------------------------------------------------------------------------- #
# MCP wrapper (lazy — importing this module never requires `mcp`)
# --------------------------------------------------------------------------- #
def build_server():
    """Construct the FastMCP server with the tools registered. Requires the `mcp` extra."""
    try:
        from mcp.server.fastmcp import FastMCP
    except ImportError as e:  # pragma: no cover
        raise RuntimeError(
            "the MCP server needs the mcp SDK — `pip install \"metric-autopsy[mcp]\"`"
        ) from e
    server = FastMCP("metric-autopsy")
    server.tool()(autopsy_report)
    server.tool()(qc_parity_report)
    server.tool()(list_metrics)
    server.tool()(demo_report)
    return server


def main():
    build_server().run()


if __name__ == "__main__":
    main()
