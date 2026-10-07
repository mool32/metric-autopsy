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

> **Status: v0.x under validation. Do not treat a verdict as validated.** Probing the
> validator itself (exploratory dev set: [`validation/probes/`](validation/probes/))
> showed that the released v0.1.1 errs in both directions:
> - **It passes useless metrics.** A metric that returns random numbers receives
>   "PASS — cleared 3 auto gates".
> - **GATE 1 blocks real biology.** "Xist is higher in female cells" dies at GATE 1
>   on the demo data because one stratum has a QC gap, even though GATE 2 retains 100% of
>   the effect. A sorted G1-vs-G2M cell-cycle control and a proliferation shift die the
>   same way once cycling cells carry ~2× more RNA.
> - **Permutations run over cells, not biological replicates.** With 3 vs 3 mice and no
>   age effect, 22 of 40 null runs report a QC-robust "effect".
>
> This branch (v0.3.0.dev0) reworks the engine: a four-field verdict (metric validity,
> design adequacy, effect, replication), replicate-level inference with a graded rule,
> depth or spike-in correction chosen by the pre-registered estimand, empirical nulls for
> the controls, and data/pre-registration hashes in every report. All 18 failures found by
> the probes now pass their regression tests — but the fixes were developed against those
> probes, so that is a development result with no confirmatory weight. Confirmatory
> validation will run on a frozen tag and a new, blind panel.

---

## Brief Summary

A metric that changes between conditions is not a finding — it might be dropout, library size, a batch effect, a factorial-interaction confound, or plain mathematics wearing a lab coat. `metric-autopsy` runs a single-cell metric through eight **gates**, each built to catch one way a number fakes biology, and decides one verdict from four fields: is the metric valid, is the design adequate, is there an effect at equal depth across biological replicates, does it replicate. It ships as a Claude Code skill, a pip package and an MCP server.

1. **Born from three real failures.** Entropy anticorrelation (ρ = −0.54, vanished on 10x, *reversed* at low depth), cardiac β (a conduction-geometry constant read as biology), and SMAD→ECM mutual information (a detection-rate confound hiding in a sex×age interaction, male-old cells detecting 2.4× fewer genes) — each survived weeks before a 45-second QC check killed it.
2. **Eight gates, four fields, one rule.** GATE 0 separates nuisance *bias* (FAIL) from *attenuation* (reported, and used as a power check); GATE 4 checks that the metric responds to an injected signal; GATE 5 judges positive and negative controls against empirical nulls — together they decide **metric validity**, which needs a demonstrated response, not just invariance. GATE 1 (stratified QC parity) is a diagnostic; GATE 2 removes the technical difference the way the pre-registered **estimand** allows (thinning to equal depth for composition, to equal spike-in capture for content, otherwise UNIDENTIFIABLE) and infers the **effect** across biological replicates (exact permutation at ≥ 4 per group, *parametric only* at 3, no verdict at ≤ 2). GATE 6 re-estimates it on independent data (**replication**). GATE 3 exports the raw scatter; GATES 4 (alternative explanations) and 7 (effect size, declared as the SESOI) are judgment the skill elicits.
3. **The reference metric is diagnosed, not just killed.** On the bundled demo (`mi_3bin`, biology identical, male-old capture degraded) GATE 0 measures attenuation (dropout −61%, depth halving −24%), GATE 1 flags the male stratum (1.94× QC ratio, 0.00 n_genes overlap), and at equal depth the raw MI difference shrinks to −13% of itself and is not detected across 16 mice: INCONCLUSIVE, not supported. The real-data run of the preprint (§4) used v0.1.1 and is under audit ([`validation/flagship_audit/`](validation/flagship_audit/)).
4. **Metric-as-plugin.** You pass `metric(data) -> float` and your factorial `obs` column names; the gates treat the metric as a black box and probe the data and its response to controlled perturbations. Metric-agnostic, domain-locked to scRNA-seq (RNA only in v1).
5. **Necessary, not sufficient (honest limit).** Passing the gates removes only the artifacts these gates know about; no correlation metric is fully depth-invariant under dropout. Every report carries the data and pre-registration hashes, and a run log counts repeated attempts at the same claim. The validator's own operating characteristics are not established yet (see Status).

Two front doors, one engine: the **skill** catches the audience inside the Claude ecosystem; the **pip package** catches everyone outside it.

---

## Using it through an agent

The tool is built to be driven by an agent — a researcher points their agent at an `.h5ad` and
the agent runs the autopsy. Three entry points, one engine:

- **Claude Code skill** — [`SKILL.md`](SKILL.md). The skill elicits a pre-registration, then runs the gates.
- **MCP server** — for any MCP-capable agent (Claude Desktop, Cursor, Cline…):
  ```bash
  # not on PyPI yet — install from the repository
  pip install "metric-autopsy[mcp] @ git+https://github.com/mool32/metric-autopsy.git"
  metric-autopsy-mcp          # stdio transport
  ```
  Register it in your agent (e.g. Claude Desktop `claude_desktop_config.json`):
  ```json
  {"mcpServers": {"metric-autopsy": {"command": "metric-autopsy-mcp"}}}
  ```
  Tools exposed: `autopsy_report` (full gate sequence on an `.h5ad`, four-field verdict,
  optional JSON report), `qc_parity_report` (GATE 1 only), `list_metrics`, and `demo_report`
  (runs on bundled synthetic data — no file needed). The CLI, the MCP server and the Python
  API share one verdict rule (`report.decide`). See
  [`src/metric_autopsy/mcp_server.py`](src/metric_autopsy/mcp_server.py).
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
├── tests/                    # synthetic gate tests, v0.3 verdict tests, audit regressions, AnnData compat
└── validation/
    ├── probes/               # exploratory dev set: probes of the validator + regression tests
    └── flagship_audit/       # audit of the preprint §4 real-data run (pinned to v0.1.1)
```

---

## Reproducing the analysis

```bash
# 1. Install (only numpy + pandas required; scipy/anndata/matplotlib are optional extras)
pip install -e ".[dev]"

# 2. Run the test suite (tests/ + the validation/probes dev set)
pytest -q

# 3. Run the reference autopsy end-to-end on the bundled demo, no downloads (~5 s);
#    --json writes the full report with data and pre-registration hashes
metric-autopsy --demo --no-stop --json demo_autopsy.json

# 4. Re-execute the worked-example notebook (~30 s). Its committed outputs and narrative
#    were produced with v0.1.1; re-executing on this branch shows the v0.3 fields.
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
