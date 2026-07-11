# Changelog

All notable changes to `metric-autopsy` are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[SemVer](https://semver.org/). Downstream papers should cite a fixed version for
comparability.

## [Unreleased]

## [0.1.1] — 2026-07-07

### Fixed
- **`download_data.py` resolves the Tabula Muris Senis FACS dataset by name** from the
  CELLxGENE Census (the old query used a non-existent slug and a non-`obs` column and fetched
  nothing); it now filters to Smart-seq2 + primary cells and maps the Census QC columns so the
  pull is gate-ready — enabling the first real-data autopsy (preprint §4: `mi_3bin` dies at
  GATE 0 on 110,824 real TMS FACS cells, reproducing the synthetic verdict).
- **GATE 6 now stratifies.** `gate6_replication` re-runs GATE 1 with the same `within` factors
  as the primary analysis and reports replication only if QC parity does not fail/STOP on the
  independent data (a confound hidden in an interaction there is no longer silently passed);
  previously it ran QC parity pooled and decided on GATE 2 alone. `run_autopsy` threads
  `within` into GATE 6. Two new tests exercise it (previously GATE 6 had zero executing tests).
- **Removed the `variance_inflation` perturbation** from GATE 0. The only whole-matrix
  reference metric (`spectral_entropy`) is computed from the *correlation* matrix and is
  therefore scale-invariant, so inflating per-gene variance moved it by a measured 0.0% — a
  false probe. GATE 0's whole-matrix probe is now `gene_subsample`, auto-enabled by
  `run_autopsy` for metrics with no bound gene pair (the matrix-perturbation path was
  previously unreachable from every entrypoint; a test now covers it). `spectral_entropy`'s
  docstring is corrected: it typically *passes* GATE 0 — it is an honest whole-matrix
  reference, not an antihero.
- **Hardened `_safe_call`** to swallow any exception a black-box metric throws on a perturbed
  input (not just four types), matching its documented "never crash the gate" contract.
- Documented that `per_cell_qc` trusts precomputed `obs` QC columns (a stale-column footgun).

## [0.1.0] — 2026-07-04

First release. Turns the internal "metric validation checklist" (written after the ACP,
cardiac-β, and MI-coupling failures) into runnable behavior with two front doors.

### Hardened
- Engine passed a 5-dimension adversarial audit; **24 findings fixed**. GATE 0 separates
  nuisance *bias* from estimator *noise* via a bootstrap null + z-test (no more failing
  noisy-but-unbiased or near-null metrics); GATE 1's distribution overlap is matchability-based
  and reports *why* a stratum flagged (ratio vs. overlap); GATE 2 guards the no-effect
  false-pass and runaway-amplification cases and checks matched-subset balance; GATE 5 anchors
  the null band to the positive control, not the negative one under test; markdown output is
  escaped; the KS fallback is correct at D=0; duplicate `var_names` are rejected. SimpleData ↔
  AnnData interchangeability (sparse X, string index) is tested. 32 tests total.
- Documentation (README, SKILL, references, manuscript, PROJECT) passed an adversarial
  docs↔code consistency review; **19 drift findings fixed** so every documented claim,
  sample output, default, and dependency matches the implementation.

### Added
- **Three agent front doors.** The **Claude Code skill** (`SKILL.md`), the **pip package**
  (`metric-autopsy` CLI / `run_autopsy` API), and an **MCP server**
  (`pip install "metric-autopsy[mcp]"`, console script `metric-autopsy-mcp`) exposing
  `autopsy_report`, `qc_parity_report`, `list_metrics`, and `demo_report` to any MCP agent
  (Claude Desktop, Cursor, Cline…). The `mcp` import is lazy; core install stays numpy+pandas.
  Agent-oriented docs: `AGENTS.md` and `llms.txt`.
- **Metric-agnostic gate engine** (`metric_autopsy.gates`): `gate0_independence`,
  `gate1_qc_parity`, `gate2_ngenes_matching`, `gate3_raw_visibility`, `gate5_controls`,
  `gate6_replication`. Gates take a black-box `metric(data) -> float`; judgment gates 4 and 7
  are elicited, not scripted.
- **Data contract** (`metric_autopsy.SimpleData`) — a numpy+pandas AnnData stand-in; real
  `anndata.AnnData` satisfies the same duck type, so either works.
- **Reference metrics** (`metric_autopsy.metrics`): `mi_3bin` (the confounded antihero),
  `pearson`, `codetected_spearman`, `norm_pearson` (library-normalized, robust),
  `spectral_entropy`.
- **QC core** (`metric_autopsy.qc`): per-cell QC, factorial strata tables, QC-parity ratios,
  n_genes matching with an overlap/KS balance check (bundled KS fallback; scipy optional).
- **Autopsy report** (`metric_autopsy.report`): gate-by-gate table + a verdict decided by the
  first blocking gate — no rescue language.
- **Two front doors**: the `metric-autopsy` Claude Code skill (`SKILL.md`) and the
  `metric-autopsy` CLI / `run_autopsy` Python API (`pip install metric-autopsy`).
- **Progressive-disclosure references** (`references/`): full gate definitions, the red-flag
  table, and the pre-registration/elicitation form.
- **Worked example** (`examples/mi_coupling_tms/`) and a synthetic `--demo` that reproduces
  the reference failure with no downloads.
- **Tests**: synthetic planted-confound suite proving the gates separate a confounded metric
  from a clean one on the same data.

[Unreleased]: https://github.com/mool32/metric-autopsy/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/mool32/metric-autopsy/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/mool32/metric-autopsy/releases/tag/v0.1.0
