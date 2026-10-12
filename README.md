[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21195679.svg)](https://doi.org/10.5281/zenodo.21195679)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CI](https://github.com/mool32/metric-autopsy/actions/workflows/ci.yml/badge.svg)](https://github.com/mool32/metric-autopsy/actions/workflows/ci.yml)

# metric-autopsy: a metric-agnostic gate system for separating biological signal from QC and technical artifacts in single-cell metrics

**Red-team a computed metric before you believe it — invert the default from "compute → believe" to "state your commitments → red-team → then believe."**

Theodor Spiro | [ORCID 0009-0004-5382-9346](https://orcid.org/0009-0004-5382-9346) | tspiro@vaika.org

📄 **No preprint yet.** The v0.1.1 manuscript ([`paper/main.pdf`](paper/main.pdf)) is superseded and will not be posted. The draft of a new manuscript on v0.3, [`paper/manuscript_v03.md`](paper/manuscript_v03.md), is paused (2026-10-11): it is not for posting and will be rewritten after v0.4, validation v2 and the external verdicts.
🧮 **Run the gates:** [`scripts/run_gates.py`](scripts/run_gates.py) · CLI `metric-autopsy --demo`
📦 **Archived release (Zenodo DOI):** [10.5281/zenodo.21195679](https://doi.org/10.5281/zenodo.21195679)
📊 **Worked-example notebook:** [`examples/mi_coupling_tms/notebook.ipynb`](examples/mi_coupling_tms/notebook.ipynb)

> **Status (2026-10-11): v0.3.0-prereg passed the pre-registered blind validation v1, within
> its scope only. The pass excludes only gross excesses of false SUPPORTED on the gene pairs where
> the metric works.**
> - **Result.** Every criterion S1–S7b passed ([`scores.txt`](https://github.com/mool32/metric-autopsy/blob/results/panel-v1/scores.txt)).
>   - Protocol: [`validation/prereg/v1.md`](validation/prereg/v1.md).
>   - Frozen tag `v0.3.0-prereg`; one attempt, [workflow run 37990113304](https://github.com/mool32/metric-autopsy/actions/runs/37990113304).
>   - Key: drand quicknet round 32929913.
>   - No useless metric was SUPPORTED: 0 of 1,630 cards.
>   - The pre-registered prediction that v0.3 fails S1 on N8 was not confirmed: 11 false
>     SUPPORTED of 790 datasets, where S1 allows 29. Its S4 part was not confirmed either (S4's
>     null cell: 58 of 2,165, limit 102, 11 of them N8's).
> - **What the pass excludes** (computed after the results;
>   [`v1_null_cells.log`](validation/exploratory/v1_null_cells.log)). A false SUPPORTED needs a
>   metric that responds, and in each key null of the pair metric (N1, N2, N5, N8) only 253–271
>   of the 790 cards had a valid or ambiguous metric. On those cards S1's cell fails with
>   probability 0.95 only at a
>   false-SUPPORTED rate of 14.3–15.3%; a validator wrong on 7.5% of them passes with probability
>   0.98–0.99. The pilot showed this mix before the key
>   (`oc_pilot.log`: a valid metric on 2,453 of the 8,910 expected pair-metric cards). At the
>   rates the pilot measured, a sound validator, with the oracle's true correction on N3 and N8,
>   passes the criteria together with probability 1.000: a property of the criteria, not a
>   prediction for the engine.
> - **Scope.**
>   - Data: scRNA-seq counts only.
>   - Metric: one, the log-normalized Pearson correlation of a gene pair (`norm_pearson`), on
>     the 8 most correlated pairs of each expression level.
>   - Backgrounds, from the CELLxGENE Census
>     ([`census_metadata.log`](validation/exploratory/census_metadata.log)):
>     - B1: human oligodendrocytes of the MSSM cohort; nuclei, 10x 3' v3. It carried every
>       condition but N7: 9,750 of the 10,540 cards.
>     - B2: mouse pancreatic islet beta cells; cells, 10x 3' v2 (39 donor ids) and v3 (11). Its
>       50 donor ids are samples, 6 of them pools of both sexes; all single females are the NOD
>       strain. It carried only the null N7 and the anchor R1.
>   - Artifacts: N2, N3 and N8 (capture loss and entry-level dropout, of v1's binomial-thinning
>     family, and per-cell capture, outside the engine's correction family). N2 (capture loss in
>     one group) did not test the correction of depth: the oracle's power to detect its raw difference was 0.21,
>     0.14 and 0.16 at the three levels, never the 0.9 that makes it informative
>     ([`v1_pilot_table.log`](validation/exploratory/v1_pilot_table.log)).
>   - The real effects ran at the saturation dose, since no level had a key dose. Of their 1,800
>     cards 206 had a valid metric, 189 of them establishable and 121 at or above the SESOI; all
>     1,285 with a blind metric were establishable too, with "metric invalid" as the correct
>     outcome. At the high level their Δ* (0.258) is about the SESOI (0.255).
>   - Not tested: GATE 6 (no card had a second dataset, so every SUPPORTED was provisional) and
>     the content estimand (only in the anchor R2).
>
>   Everything else is not validated: other metrics (`mi_3bin` included), other data and
>   other artifacts.
> - **Limitations.**
>   - **The SESOI and GATE 4's dose came from the pilot**, that is, from the oracle, on every
>     card. A user sets both without a pilot.
>   - **The truth about validity is circular.** The oracle calls a metric valid when it responds
>     to the injection that GATE 4 itself uses, so v1's S5 measures GATE 4's sampling error and
>     could not find a probe built the wrong way (R1's was one). A next validation needs a truth
>     about validity that does not come from GATE 4's mechanism.
>   - **v0.3 does not correct capture that varies from cell to cell, nor entry-level dropout.**
>     Its depth correction brings whole groups to one depth.
>   - **Design effects** (datasets share donors) exceed 1.5 only for S3's null and invalid strata
>     (4.74, 5.62), S5's null cell (5.28) and S3's E1 condition (1.54). In the false-SUPPORTED
>     cells they are at most 1.055; at N5's, S1's cell passes a doubled error with probability
>     0.051, just above the principle's 0.05.
>   - **Anchor R1, the one real positive control expected to succeed, failed: two of its three
>     claims were outside their allowed sets.** It ran once, after the results and outside v1
>     ([D1](validation/prereg/DEVIATIONS.md);
>     [`results/panel-v1-r1`](https://github.com/mool32/metric-autopsy/tree/results/panel-v1-r1)).
>     Read after the outcome
>     ([`b2_sex_structure.log`](validation/exploratory/b2_sex_structure.log)): the Y genes'
>     UNIDENTIFIABLE reads B2 correctly (within a stage the sexes meet only inside the pools, so
>     the protocol's allowed set was wrong); Xist's "metric invalid" by GATE 4 is false (next
>     point); the sham's NO DETECTABLE EFFECT is correct. v0.3 confirmed no real biological
>     difference.
>   - **GATE 4's module probe measures the response over all cells.** For a metric that is a mean
>     over cells that is the share of cells with the marker times the response in them, which is
>     below ln 2 for a well-detected marker. With Xist's response of 0.631 in B2's Xist cells, R1's
>     probe needed at least 39.6% of cells with Xist (0.25 / 0.631), and B2 has 25.8%: with
>     anchors.py's constants the failure was predictable from B2's composition. With the module's default (2-fold in 30% of cells) a valid
>     mean-expression metric of a well-detected gene cannot pass GATE 4 at a SESOI of 0.5: at most
>     0.3 × ln 2 = 0.208 against delta_min 0.25
>     ([`gate4_module_analytic.log`](validation/exploratory/gate4_module_analytic.log); simulated
>     in [`gate4_module_dilution.log`](validation/exploratory/gate4_module_dilution.log)). At a
>     depth of 2,000 counts it passes a marker of single counts for the wrong reason: of a
>     response of 0.360 the sham's loss of the single counts makes 0.178, the 2-fold change
>     0.182, since the sham thins the marker and the injection does not. The module probe is
>     reached through the Python API only: the CLI and the MCP server offer the coupling probe.
>     v1's panel used it only on the useless metric N6c.
> - **Error rates by level and truth** (pre-registered secondary outcomes, v1.md sections 5 and
>   7; [`v1_by_level.log`](validation/exploratory/v1_by_level.log)):
>   - All 2,458 cards with a valid metric are at the high expression level. At the medium and
>     low levels the metric is blind (55 medium-level cards of N7 are ambiguous), and no card
>     there was SUPPORTED (0 of 6,988).
>   - In the twelve high-level null cells, after Holm's correction over the twelve (decided after
>     the results) no cell differs from the nominal 2.5%; without it, N3 at f = 0.4's interval
>     excludes it. The largest are N3 at f = 0.4, 7/83 = 8.4% (95% CI 3.5–16.6%; Holm 0.058),
>     and N8, 11/253 = 4.3% (2.2–7.6%; Holm 0.61). Both are hypotheses for v0.4 with a
>     mechanism, not measured defects. Found afterwards (exploratory,
>     [`v1_analyzers.log`](validation/exploratory/v1_analyzers.log)): N3's 7 were all made by the
>     depth correction, which enlarged raw differences the replicate-level test had not detected;
>     a test over donors without the correction gives none of the 83.
>   - GATE 0 refused 31 cards, all on N7 (26/276 = 9.4% at the high level, 5 at the medium).
>     GATE 0 refuses where extra dropout would bias the metric beyond its tolerance, as it is built
>     to, and a refusal is allowed, never definite; but on these null data the refusals cost
>     decisiveness (S3's N7 cell: 204 of 237). Decided after the results: 29 of the 31 fell on
>     three of B2's gene pairs, where dropout shifted the correlation beyond the tolerance; eight
>     tests found no association with B2's chemistry or depth (one at p = 0.049)
>     ([`b2_depth_by_chemistry.log`](validation/exploratory/b2_depth_by_chemistry.log)).
>   - Decided after the results: N8 on its cards with a valid metric 8/221 = 3.6%. The
>     development probe p16 (11.3% on simulated data) was 8.1 times N8's real rate over all
>     levels, and at most 3.7 times on the datasets at risk; it also ran with a SESOI of 0.1 and
>     a dose of 2.0, where B1's high level had 0.255 and 2.5, so simulation and settings are
>     confounded.
>
>   Planned for v0.4 (a draft specification, [`validation/plans/v0.4_spec.md`](validation/plans/v0.4_spec.md)):
>   GATE 4's response measured with the claim's own replicate-level test, the signal injected in
>   the group where the probed genes have more counts per cell, whatever the claim's wording; a
>   module probe that creates no response by losing counts; and a correction that can only remove
>   a conclusion, never create one. A correction for per-cell capture is postponed to v0.5. Then
>   a new pre-registration, and the test on published claims
>   ([`validation/plans/real_claims.md`](validation/plans/real_claims.md)).
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
>   age effect, 22 of 40 null runs report a QC-robust "effect" (an exact test over the mice
>   cannot reach p < 0.05 there: its smallest p is 0.10).

---

## Brief Summary

A metric that changes between conditions is not a finding — it might be dropout, library size, a batch effect, a factorial-interaction confound, or plain mathematics wearing a lab coat. `metric-autopsy` runs a single-cell metric through eight **gates**, each built to catch one way a number fakes biology, and decides one verdict from four fields: is the metric valid, is the design adequate, is there an effect at equal depth across biological replicates, does it replicate. It ships as a Claude Code skill, a pip package and an MCP server.

1. **Born from three real failures.** Entropy anticorrelation (it vanished on 10x and *reversed* at low depth), cardiac β (a conduction-geometry constant read as biology), and SMAD→ECM mutual information (a fall in old male cells that followed their lower detection rate; see the audit, [`validation/flagship_audit/REPORT.md`](validation/flagship_audit/REPORT.md)) — each survived weeks before a QC check killed it.
2. **Eight gates, four fields, one rule.** GATE 0 separates nuisance *bias* (FAIL) from *attenuation* (reported, and used as a power check); GATE 4 checks that the metric responds to an injected signal; GATE 5 judges positive and negative controls against empirical nulls — together they decide **metric validity**, which needs a demonstrated response, not just invariance. GATE 1 (stratified QC parity) is a diagnostic; GATE 2 removes the technical difference the way the pre-registered **estimand** allows (thinning to equal depth for composition, to equal spike-in capture for content, otherwise UNIDENTIFIABLE) and infers the **effect** across biological replicates (exact permutation at ≥ 4 per group, *parametric only* at 3, no verdict at ≤ 2). GATE 6 re-estimates it on independent data (**replication**). GATE 3 exports the raw scatter; GATES 4 (alternative explanations) and 7 (effect size, declared as the SESOI) are judgment the skill elicits.
3. **The reference metric is diagnosed, not just killed.** On the bundled demo (`mi_3bin`, biology identical, male-old capture degraded) GATE 0 measures attenuation (dropout −61%, depth halving −24%), GATE 1 flags the male stratum (1.94× QC ratio, 0.00 n_genes overlap), and at equal depth the raw MI difference shrinks to −13% of itself and is not detected across 16 mice: INCONCLUSIVE, not supported. The real-data run of the preprint (§4) used v0.1.1. Its audit ([`validation/flagship_audit/REPORT.md`](validation/flagship_audit/REPORT.md)) found that v0.1.1's verdict reproduces, but the design §4 describes does not hold.
4. **Metric-as-plugin.** You pass `metric(data) -> float` and your factorial `obs` column names; the gates treat the metric as a black box and probe the data and its response to controlled perturbations. Metric-agnostic, domain-locked to scRNA-seq (RNA only in v1).
5. **Necessary, not sufficient (honest limit).** Passing the gates removes only the artifacts these gates know about; no correlation metric is fully depth-invariant under dropout. Every report carries the data and pre-registration hashes, and a run log counts repeated attempts at the same claim. The validator's operating characteristics are established only within the scope of validation v1 (see Status): one metric, background B1 (B2 only for one null), the artifacts N2, N3 and N8; and there the pass excludes only gross excesses of false SUPPORTED on the pairs where the metric works.

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
| Tabula Muris Senis (FACS) | [CELLxGENE Census](https://cellxgene.cziscience.com/) / figshare | 110K cells, 23 tissues | worked example (real-data path); young (3 months) against the Census's "20 months" stage (24-month males, 21-month females; see the audit) |
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

If you use metric-autopsy, please cite the version of the software you used:

> Spiro, T. (2026). *metric-autopsy* (the version you used). Zenodo. https://doi.org/10.5281/zenodo.21195679

Citation metadata is in [`CITATION.cff`](CITATION.cff).

---

## Contact

Theodor Spiro — tspiro@vaika.org

## License

- **Code** (`src/`, `scripts/`, `tests/`, `examples/*.py`): **MIT** — see [LICENSE](LICENSE).
- **Manuscript and figures** (`paper/`, `examples/**/figures/`): **CC-BY-4.0** — see [`paper/LICENSE-CC-BY-4.0.md`](paper/LICENSE-CC-BY-4.0.md).
