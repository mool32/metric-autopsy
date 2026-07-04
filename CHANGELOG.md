# Changelog

All notable changes to `metric-autopsy` are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[SemVer](https://semver.org/). Downstream papers should cite a fixed version for
comparability.

## [Unreleased]

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

[Unreleased]: https://github.com/mool32/metric-autopsy/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/mool32/metric-autopsy/releases/tag/v0.1.0
