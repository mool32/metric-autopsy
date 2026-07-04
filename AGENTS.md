# AGENTS.md — operating guide for AI agents

You (an AI agent) are working with **metric-autopsy**: a metric-agnostic gate system that
red-teams a computed single-cell metric to tell biological signal apart from QC, technical,
and mathematical artifacts. Its thesis is that metric validation should be *agent-driven* —
so driving it is exactly what you're for. Default posture: **state commitments → red-team →
then believe**, never "compute → believe."

## Three ways to drive it

1. **This repo, directly (you're here).** Use the CLI or the Python API.
   - Fastest check: `python scripts/run_gates.py --demo --no-stop` (reference `mi_3bin` failure, no data).
   - Real data: `metric-autopsy --h5ad DATA.h5ad --metric mi_3bin --gene-a A --gene-b B --group-col age --groups young old --within sex --pos-pair Actb Gapdh --neg-pair G1 G2`
   - Python: `from metric_autopsy import run_autopsy, metrics` → bind a metric with `functools.partial` → `run_autopsy(...)` → `.to_markdown()`.
2. **As a Claude Code skill.** `SKILL.md` — elicit a pre-registration first, then run the gates, then emit the autopsy. This is the skill's required behavior; follow it verbatim when invoked.
3. **As an MCP server** (for any MCP agent). `pip install "metric-autopsy[mcp]"` then `metric-autopsy-mcp`. Tools: `autopsy_report`, `qc_parity_report`, `list_metrics`, `demo_report`. See `src/metric_autopsy/mcp_server.py`.

## The one rule that matters

A metric that changes between conditions is **not** a finding. Run it through the gates and
**stop at the first blocking one** (FAIL or STOP) — a later gate cannot rescue an earlier
failure. GATE 1 (factorial QC parity) is the cheapest and catches the most; run it first, and
always stratify by every factor (`--within sex tissue …`) — confounds hide in interactions
that pooling erases.

## Layout (source of truth = `src/metric_autopsy/`)

| Path | What |
|---|---|
| `src/metric_autopsy/gates.py` | `gate0…gate6` as metric-agnostic functions |
| `src/metric_autopsy/qc.py` | factorial QC parity + n_genes matching (GATES 1–2 core) |
| `src/metric_autopsy/metrics.py` | reference metrics (`mi_3bin`, `norm_pearson`, …) |
| `src/metric_autopsy/report.py` | `run_autopsy` orchestration + verdict |
| `src/metric_autopsy/mcp_server.py` | MCP tools |
| `references/` | full gate definitions, red flags, the pre-reg form (load on demand) |
| `examples/mi_coupling_tms/` | worked example + executed notebook |
| `paper/manuscript.md` | the methods write-up |

## Working on the code

- **Install:** `pip install -e ".[dev]"` (core needs only numpy+pandas; scipy/anndata/matplotlib/mcp are optional extras).
- **Test:** `pytest -q` — 32 tests must stay green (`tests/`). Add a regression test for any behavior change.
- **Contract:** the gates take a black-box `metric(data) -> float`; `data` is any object with `.X`, `.obs`, `.var_names`, and `data[mask]` — real `anndata.AnnData` or the bundled `SimpleData`. Keep gates metric-agnostic; keep optional deps lazily imported.
- **Adding a reference metric:** add a keyword-only `fn(data, *, gene_a, gene_b) -> float` to `metrics.py`, then it's usable by name in the CLI and MCP tools.

## Don't

- Don't report a metric as biology without running the gates.
- Don't soften a FAIL ("no rescue language") — record the number that killed it.
- Don't make an optional dependency (scipy/anndata/matplotlib/mcp) a hard import.
