# AGENTS.md — operating guide for AI agents

You (an AI agent) are working with **metric-autopsy**: a metric-agnostic gate system that
red-teams a computed single-cell metric to tell biological signal apart from QC, technical,
and mathematical artifacts. Its thesis is that metric validation should be *agent-driven* —
so driving it is exactly what you're for. Default posture: **state commitments → red-team →
then believe**, never "compute → believe."

> **Status: v0.x under validation.** The validator has known failures in both directions:
> it can pass useless metrics and block real biology. They are listed in
> `validation/probes/README.md`. Never present a verdict from this version as validated.

## Three ways to drive it

1. **This repo, directly (you're here).** Use the CLI or the Python API.
   - Fastest check: `python scripts/run_gates.py --demo --no-stop` (reference `mi_3bin` failure, no data).
   - Real data: `metric-autopsy --h5ad DATA.h5ad --metric mi_3bin --gene-a A --gene-b B --group-col age --groups young old --within sex --replicate-col mouse --estimand composition --direction decrease --sesoi 0.1 --pos-pair Actb Gapdh --neg-pair G1 G2 --json autopsy.json`
   - Python: `from metric_autopsy import run_autopsy, metrics` → bind a metric with `functools.partial` → `run_autopsy(..., replicate_col="mouse", prereg={"estimand": "composition", "direction": "decrease", "sesoi": 0.1})` → `.to_markdown()` / `.save_json(path)`.
2. **As a Claude Code skill.** `SKILL.md` — elicit a pre-registration first, then run the gates, then emit the autopsy. This is the skill's required behavior; follow it verbatim when invoked.
3. **As an MCP server** (for any MCP agent). `pip install "metric-autopsy[mcp] @ git+https://github.com/mool32/metric-autopsy.git"` then `metric-autopsy-mcp` (not on PyPI yet). Tools: `autopsy_report`, `qc_parity_report`, `list_metrics`, `demo_report`. See `src/metric_autopsy/mcp_server.py`.

## The one rule that matters

A metric that changes between conditions is **not** a finding. Declare the commitments first
(estimand, direction, replicate unit, SESOI — `references/prereg_template.md`), then run the gates, then
report the verdict that `report.decide` gives from four fields: **metric validity, design
adequacy, effect, replication**. An invalid metric (nuisance *bias*, a failed negative control,
no response to an injected signal) stops the analysis — a later field cannot rescue it. The
effect is judged across biological replicates at equal depth (or spike-in capture), never
across cells. SUPPORTED needs a metric whose response is demonstrated (positive control or
injected signal), a replicate-level DETECTED effect in the pre-registered direction and resolved
judgment gates. Always stratify
by every factor (`--within sex tissue …`) — confounds hide in interactions that pooling
erases. GATE 1 is the quickest look at them; it is a diagnostic now, not a kill switch.

## Layout (source of truth = `src/metric_autopsy/`)

| Path | What |
|---|---|
| `src/metric_autopsy/gates.py` | `gate0…gate6` as metric-agnostic functions |
| `src/metric_autopsy/qc.py` | factorial QC parity (GATE 1 core) + legacy n_genes matching |
| `src/metric_autopsy/effect.py` | estimand-dependent correction, replicate-level effect, graded rule, power |
| `src/metric_autopsy/equalize.py` | binomial thinning to equal depth / spike-in capture |
| `src/metric_autopsy/stats.py` | t / permutation tests (exact or Monte Carlo), TOST, MDE — no scipy needed |
| `src/metric_autopsy/injected_signal.py` | known construct changes for GATE 4 (`coupling`, `module`) |
| `src/metric_autopsy/provenance.py` | data / pre-registration hashes, versions, run log |
| `src/metric_autopsy/metrics.py` | reference metrics (`mi_3bin`, `norm_pearson`, …) |
| `src/metric_autopsy/report.py` | `run_autopsy` orchestration, the four fields, `decide` (the verdict rule) |
| `src/metric_autopsy/mcp_server.py` | MCP tools |
| `references/` | full gate definitions, red flags, the pre-reg form (load on demand) |
| `examples/mi_coupling_tms/` | worked example + executed notebook |
| `paper/manuscript.md` | the methods write-up |

## Working on the code

- **Install:** `pip install -e ".[dev]"` (core needs only numpy+pandas; scipy/anndata/matplotlib/mcp are optional extras).
- **Test:** `pytest -q` must stay green. It collects `tests/` and the dev set `validation/probes/`. Known validator failures are `xfail(strict=True)`: a fix that makes one pass shows up as an XPASS error until its marker is removed. Add a regression test for any behavior change, and never edit a dev probe's expected verdict to make it pass.
- **Contract:** the gates take a black-box `metric(data) -> float`; `data` is any object with `.X`, `.obs`, `.var_names`, and `data[mask]` — real `anndata.AnnData` or the bundled `SimpleData`. Keep gates metric-agnostic; keep optional deps lazily imported.
- **One rule:** the verdict is `report.decide`, shared by the API, the CLI and MCP; never re-implement it in a front end. Numbers meant for the paper come only from scripts in this repo, with the hashes from the JSON report.
- **Adding a reference metric:** add a keyword-only `fn(data, *, gene_a, gene_b) -> float` to `metrics.py`, then it's usable by name in the CLI and MCP tools.

## Don't

- Don't report a metric as biology without running the gates.
- Don't soften a FAIL ("no rescue language") — record the number that killed it.
- Don't re-run a claim with new settings until it passes: the run log counts the attempts.
- Don't make an optional dependency (scipy/anndata/matplotlib/mcp) a hard import.
