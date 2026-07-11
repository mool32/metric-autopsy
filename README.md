[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21195679.svg)](https://doi.org/10.5281/zenodo.21195679)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Preprint](https://img.shields.io/badge/Preprint-forthcoming-lightgrey)](paper/manuscript.md)
[![CI](https://github.com/mool32/metric-autopsy/actions/workflows/ci.yml/badge.svg)](https://github.com/mool32/metric-autopsy/actions/workflows/ci.yml)

# metric-autopsy: a metric-agnostic gate system for separating biological signal from QC and technical artifacts in single-cell metrics

**Red-team a computed metric before you believe it — invert the default from "compute → believe" to "state your commitments → red-team → then believe."**

Theodor Spiro | [ORCID 0009-0004-5382-9346](https://orcid.org/0009-0004-5382-9346) | tspiro@vaika.org

📄 **Preprint:** [`paper/main.pdf`](paper/main.pdf) — arXiv-ready (`q-bio.QM`); see [`paper/ARXIV_SUBMISSION.md`](paper/ARXIV_SUBMISSION.md)
🧮 **Run the gates:** [`scripts/run_gates.py`](scripts/run_gates.py) · CLI `metric-autopsy --demo`
📦 **Archived release (Zenodo DOI):** [10.5281/zenodo.21195679](https://doi.org/10.5281/zenodo.21195679)
📊 **Worked-example notebook:** [`examples/mi_coupling_tms/notebook.ipynb`](examples/mi_coupling_tms/notebook.ipynb)

---

## Brief Summary

A metric that changes between conditions is not a finding — it might be dropout, library size, a batch effect, a factorial-interaction confound, or plain mathematics wearing a lab coat. `metric-autopsy` runs a single-cell metric through eight **gates**, each built to catch one way a number fakes biology, and stops at the first one it fails. It ships as both a Claude Code skill and a pip package.

1. **Born from three real failures.** Entropy anticorrelation (ρ = −0.54, vanished on 10x, *reversed* at low depth), cardiac β (a conduction-geometry constant read as biology), and SMAD→ECM mutual information (a detection-rate confound hiding in a sex×age interaction, male-old cells detecting 2.4× fewer genes) — each survived weeks before a 45-second QC check killed it.
2. **Eight gates, six automatic, five that can kill.** Mathematical independence, factorial QC parity, n_genes matching, stratified controls, and cross-dataset replication run from the data *and can block* a metric; raw-data visibility (GATE 3) also runs automatically but only exports the scatter and a dropout hint for you to read — it never blocks on its own; GATE 4 (does it measure what you think) and GATE 7 (is the effect size meaningful) are judgment gates the skill *elicits*, not scripts.
3. **The reference metric dies 0/N.** On the worked example, `mi_3bin` fails GATE 0 (expectation shifts 61%, z = 44.8, under simulated dropout) and GATE 1 (male stratum 1.94× QC ratio, 0.00 n_genes overlap); a library-normalized reference passes — the confound is *avoidable*, not universal.
4. **Metric-as-plugin.** You pass `metric(data) -> float` and your factorial `obs` column names; the gates treat the metric as a black box and probe the data and its response to controlled perturbations. Metric-agnostic, domain-locked to scRNA-seq QC.
5. **Necessary, not sufficient (honest limit).** Passing all eight gates removes only the artifacts these gates know about; no correlation metric is fully depth-invariant under dropout. The engine ships with 35 tests and was itself put through an adversarial code + docs audit (24 + 19 findings fixed) plus a follow-up integrity pass.

Two front doors, one engine: the **skill** catches the audience inside the Claude ecosystem; the **pip package** catches everyone outside it.

---

## Using it through an agent

The tool is built to be driven by an agent — a researcher points their agent at an `.h5ad` and
the agent runs the autopsy. Three entry points, one engine:

- **Claude Code skill** — [`SKILL.md`](SKILL.md). The skill elicits a pre-registration, then runs the gates.
- **MCP server** — for any MCP-capable agent (Claude Desktop, Cursor, Cline…):
  ```bash
  pip install "metric-autopsy[mcp]"
  metric-autopsy-mcp          # stdio transport
  ```
  Register it in your agent (e.g. Claude Desktop `claude_desktop_config.json`):
  ```json
  {"mcpServers": {"metric-autopsy": {"command": "metric-autopsy-mcp"}}}
  ```
  Tools exposed: `autopsy_report` (full gate sequence on an `.h5ad`), `qc_parity_report`
  (GATE 1 only), `list_metrics`, and `demo_report` (runs on bundled synthetic data — no file
  needed). See [`src/metric_autopsy/mcp_server.py`](src/metric_autopsy/mcp_server.py).
- **CLI** — `metric-autopsy --demo` (or `--h5ad …`), scriptable from any agent shell.

Agents landing in the repo should read [`AGENTS.md`](AGENTS.md); [`llms.txt`](llms.txt) is a
curated doc map for LLMs.

---

## Datasets

| Dataset | Source | N | Notes |
|---|---|---|---|
| Synthetic confound | bundled (`metric_autopsy.cli.demo_data`) | 1,600 cells × 40 genes | planted sex×age QC confound; no download; powers `--demo` and the tests |
| Tabula Muris Senis (FACS) | [CELLxGENE Census](https://cellxgene.cziscience.com/) / figshare | 110K cells, 23 tissues | worked example (real-data path); mouse young 3mo vs old 24mo |
| Human skin fibroblasts (10x) | [CELLxGENE](https://cellxgene.cziscience.com/) | 84K cells, 179 donors | GATE 6 cross-platform replication |

---

## Repository structure

```
├── SKILL.md                  # the skill: elicit pre-reg → run gates → autopsy
├── AGENTS.md · llms.txt      # agent operating guide + LLM-friendly doc map
├── references/               # progressive-disclosure depth (gates, red flags, prereg form)
├── src/metric_autopsy/       # engine — single source of truth (gates, qc, metrics, report, mcp_server)
├── scripts/run_gates.py      # thin CLI wrapper the skill invokes
├── examples/mi_coupling_tms/ # worked example: TMS → MI → 0/N → "not biology" (+ notebook, figures)
├── paper/                    # manuscript + figures  (CC-BY-4.0)
└── tests/                    # 35 tests: synthetic gates + audit regressions + AnnData compat
```

---

## Reproducing the analysis

```bash
# 1. Install (only numpy + pandas required; scipy/anndata/matplotlib are optional extras)
pip install -e ".[dev]"

# 2. Run the test suite (32 tests)
pytest -q

# 3. Reproduce the reference failure end-to-end, no downloads (~2 s)
metric-autopsy --demo --no-stop

# 4. Re-execute the worked-example notebook (~30 s)
jupyter nbconvert --to notebook --execute --inplace examples/mi_coupling_tms/notebook.ipynb
```

To run the worked example on the real data instead of the synthetic confound, run
`python examples/mi_coupling_tms/download_data.py --which all` and point the notebook's data
cell at the downloaded `.h5ad` — every downstream cell is unchanged, because the gates accept
any AnnData.

---

## Citation

If you use this work, please cite:

> Spiro, T. (2026). *Metric autopsy: a metric-agnostic gate system for separating biological signal from QC and technical artifacts in single-cell metrics.* Preprint.

And the archived software release (when applicable):

> Spiro, T. (2026). *metric-autopsy* (v0.1.1). Zenodo. https://doi.org/10.5281/zenodo.21195679

Citation metadata is in [`CITATION.cff`](CITATION.cff).

---

## Contact

Theodor Spiro — tspiro@vaika.org

## License

- **Code** (`src/`, `scripts/`, `tests/`, `examples/*.py`): **MIT** — see [LICENSE](LICENSE).
- **Manuscript and figures** (`paper/`, `examples/**/figures/`): **CC-BY-4.0** — see [`paper/LICENSE-CC-BY-4.0.md`](paper/LICENSE-CC-BY-4.0.md).
