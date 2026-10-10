[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21195679.svg)](https://doi.org/10.5281/zenodo.21195679)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Preprint](https://img.shields.io/badge/Preprint-forthcoming-lightgrey)](paper/manuscript.md)
[![CI](https://github.com/mool32/metric-autopsy/actions/workflows/ci.yml/badge.svg)](https://github.com/mool32/metric-autopsy/actions/workflows/ci.yml)

# metric-autopsy: a metric-agnostic gate system for separating biological signal from QC and technical artifacts in single-cell metrics

**Red-team a computed metric before you believe it — invert the default from "compute → believe" to "state your commitments → red-team → then believe."**

Theodor Spiro | [ORCID 0009-0004-5382-9346](https://orcid.org/0009-0004-5382-9346) | tspiro@vaika.org

📄 **Preprint (v0.1.1, superseded; do not upload):** [`paper/main.pdf`](paper/main.pdf); see [`paper/ARXIV_SUBMISSION.md`](paper/ARXIV_SUBMISSION.md). A full draft of the new manuscript on v0.3, [`paper/manuscript_v03.md`](paper/manuscript_v03.md), is under the owner's review; it is not a preprint yet.
🧮 **Run the gates:** [`scripts/run_gates.py`](scripts/run_gates.py) · CLI `metric-autopsy --demo`
📦 **Archived release (Zenodo DOI):** [10.5281/zenodo.21195679](https://doi.org/10.5281/zenodo.21195679)
📊 **Worked-example notebook:** [`examples/mi_coupling_tms/notebook.ipynb`](examples/mi_coupling_tms/notebook.ipynb)

> **Status (2026-10-10): v0.3.0-prereg passed the pre-registered blind validation v1, within
> its scope only.**
> - **Result.** Every criterion S1–S7b passed ([`scores.txt`](https://github.com/mool32/metric-autopsy/blob/results/panel-v1/scores.txt)).
>   - Protocol: [`validation/prereg/v1.md`](validation/prereg/v1.md).
>   - Frozen tag `v0.3.0-prereg`; one attempt, [workflow run 37990113304](https://github.com/mool32/metric-autopsy/actions/runs/37990113304).
>   - Key: drand quicknet round 32929913.
>
>   The pre-registered prediction that v0.3 fails S1 on N8 was not confirmed: 11 false
>   SUPPORTED of 790 datasets, where S1 allows 29.
> - **Scope.**
>   - Data: scRNA-seq counts only.
>   - Metric: one, the log-normalized Pearson correlation of a gene pair (`norm_pearson`).
>   - Backgrounds, from the CELLxGENE Census
>     ([metadata](validation/exploratory/census_metadata.py), workflow run 38033204255):
>     - B1: human oligodendrocytes of the MSSM cohort; nuclei, 10x 3' v3.
>     - B2: mouse pancreatic islet beta cells; cells, 10x 3' v2 and v3; 50 donor ids (samples),
>       6 of them pools of both sexes; all single females are the NOD strain.
>   - Artifacts: the binomial-thinning family, plus N8's per-cell variable capture.
>
>   Everything else is not validated: other metrics (`mi_3bin` included), other data and
>   other artifacts.
> - **Limitations.**
>   - **v0.3 does not correct capture that varies from cell to cell.** Its depth correction
>     brings whole groups to one depth, and per-cell variable capture is outside that family.
>   - On these backgrounds the metric is valid only at high expression (numbers below).
>   - The design effects reach 5.62 (datasets share donors), so on other donors the error rates
>     are known less precisely.
>   - **Anchor R1, the only real positive control, did not give the known sex difference.** It
>     ran once, after the results and outside v1 ([D1](validation/prereg/DEVIATIONS.md);
>     [`results/panel-v1-r1`](https://github.com/mool32/metric-autopsy/tree/results/panel-v1-r1);
>     [`b2_sex_structure.log`](validation/exploratory/b2_sex_structure.log)). Two of its three
>     claims are outside their allowed sets, for different causes:
>     - Y genes (and the design of the Xist claim): UNIDENTIFIABLE. That reads B2 correctly:
>       within a development stage the sexes meet only inside the 6 pooled ids. The protocol
>       erred: R1's allowed set counted the ids of each sex over all stages, and checked neither
>       that an id is one animal nor that both sexes share a stratum.
>     - Xist: "metric invalid", which is false: a defect of GATE 4 (next point).
>     - The sham, female against female: NO DETECTABLE EFFECT, which is correct.
>
>     v0.1.1 did not confirm Xist because of GATE 1, on TMS. v0.3 did not either, because of
>     GATE 4 and the design, on B2.
>   - **GATE 4's module probe measures the response over all cells.** So for a marker present in
>     only some cells it can call a valid metric invalid. In R1, 74.2% of B2's cells have no
>     Xist. The probe's 2-fold Xist signal moved the metric by +0.163 over all cells, below
>     delta_min 0.25, and by +0.631 over the cells with Xist. This is outside v1's scope (v1's
>     panel used the module probe only on the useless metric N6c) and a target of v0.4.
>   - **With the module probe's default, a valid mean-expression metric of a gene does not pass
>     GATE 4 at a SESOI of 0.5.** `injected_signal.module` raises the marker 2-fold in a random 30%
>     of cells, so the response over all cells is about 0.3 × the share of cells with the marker ×
>     the response in those cells, which for a well-detected marker is at most log 2 ≈ 0.69. At the
>     depth of B2's cells it stayed below delta_min 0.25 even with the marker in every cell
>     (largest upper bound +0.206 on synthetic counts, +0.191 on B2's Xist cells;
>     [`gate4_module_dilution.log`](validation/exploratory/gate4_module_dilution.log)). The module
>     probe is reached through the Python API only: the CLI and the MCP server offer the coupling
>     probe.
> - **The error rates by expression level** (exploratory, not a criterion;
>   [`v1_by_level.log`](validation/exploratory/v1_by_level.log), from `scores.json`):
>   - All 2,458 cards with a valid metric are at the high expression level. At the medium and
>     low levels the metric is blind (55 medium-level cards are ambiguous), and no card there
>     was SUPPORTED. So only 253 of N8's 790 cards were at risk.
>   - N8 at the high level:
>     - false SUPPORTED on 11/253 = 4.3% (95% CI 2.2–7.6%), at a nominal 2.5%;
>     - SUPPORTED or against the direction on 22/253 = 8.7%, at a nominal 5%.
>   - N3 f = 0.4 at the high level: 7/83 = 8.4% and 10/83 = 12.0%.
>   - For comparison, at the high level:
>     - N1, N5 and N7 gave SUPPORTED or against the direction on 5.9%, 3.8% and 4.0%. N7
>       permutes the group labels over B2's donor ids: a null at the level of the id, not
>       of the mouse;
>     - the pure nulls N1, N5 and N7 together gave a false SUPPORTED on 18/812 = 2.2%.
>   - Had N8 kept its high-level rate of 4.3% on all 790 cards, about 34 false SUPPORTED would
>     be expected, against S1's 29: S1 would pass with probability ≈ 0.20. The development
>     probe p16 gave 11.3%.
>
>   A correction for per-cell capture is planned for v0.4, with a new pre-registration.
> - **Next.** Step 4 of the validation plan, external verdicts, is ahead.
> - **What a verdict is.** A statement about one dataset under the declared assumptions, not
>   about biology.
>
> The released v0.1.1 errs in both directions on its own probes (exploratory dev set:
> [`validation/probes/`](validation/probes/)):
> - **It passes useless metrics.** A metric that returns random numbers receives
>   "PASS — cleared 3 auto gates".
> - **GATE 1 blocks real biology.** "Xist is higher in female cells" dies at GATE 1 on the
>   demo data because one stratum has a QC gap, even though GATE 2 retains 100% of the effect.
> - **Permutations run over cells, not biological replicates.** With 3 vs 3 mice and no
>   age effect, 22 of 40 null runs report a QC-robust "effect".

---

## Brief Summary

A metric that changes between conditions is not a finding — it might be dropout, library size, a batch effect, a factorial-interaction confound, or plain mathematics wearing a lab coat. `metric-autopsy` runs a single-cell metric through eight **gates**, each built to catch one way a number fakes biology, and decides one verdict from four fields: is the metric valid, is the design adequate, is there an effect at equal depth across biological replicates, does it replicate. It ships as a Claude Code skill, a pip package and an MCP server.

1. **Born from three real failures.** Entropy anticorrelation (ρ = −0.54, vanished on 10x, *reversed* at low depth), cardiac β (a conduction-geometry constant read as biology), and SMAD→ECM mutual information (a detection-rate confound hiding in a sex×age interaction, male-old cells detecting 2.4× fewer genes) — each survived weeks before a 45-second QC check killed it.
2. **Eight gates, four fields, one rule.** GATE 0 separates nuisance *bias* (FAIL) from *attenuation* (reported, and used as a power check); GATE 4 checks that the metric responds to an injected signal; GATE 5 judges positive and negative controls against empirical nulls — together they decide **metric validity**, which needs a demonstrated response, not just invariance. GATE 1 (stratified QC parity) is a diagnostic; GATE 2 removes the technical difference the way the pre-registered **estimand** allows (thinning to equal depth for composition, to equal spike-in capture for content, otherwise UNIDENTIFIABLE) and infers the **effect** across biological replicates (exact permutation at ≥ 4 per group, *parametric only* at 3, no verdict at ≤ 2). GATE 6 re-estimates it on independent data (**replication**). GATE 3 exports the raw scatter; GATES 4 (alternative explanations) and 7 (effect size, declared as the SESOI) are judgment the skill elicits.
3. **The reference metric is diagnosed, not just killed.** On the bundled demo (`mi_3bin`, biology identical, male-old capture degraded) GATE 0 measures attenuation (dropout −61%, depth halving −24%), GATE 1 flags the male stratum (1.94× QC ratio, 0.00 n_genes overlap), and at equal depth the raw MI difference shrinks to −13% of itself and is not detected across 16 mice: INCONCLUSIVE, not supported. The real-data run of the preprint (§4) used v0.1.1. Its audit ([`validation/flagship_audit/REPORT.md`](validation/flagship_audit/REPORT.md)) found that v0.1.1's verdict reproduces, but the design §4 describes does not hold.
4. **Metric-as-plugin.** You pass `metric(data) -> float` and your factorial `obs` column names; the gates treat the metric as a black box and probe the data and its response to controlled perturbations. Metric-agnostic, domain-locked to scRNA-seq (RNA only in v1).
5. **Necessary, not sufficient (honest limit).** Passing the gates removes only the artifacts these gates know about; no correlation metric is fully depth-invariant under dropout. Every report carries the data and pre-registration hashes, and a run log counts repeated attempts at the same claim. The validator's operating characteristics are established only within the scope of validation v1 (see Status): one metric, two backgrounds, the binomial-thinning family and N8.

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
