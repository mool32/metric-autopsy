# Changelog

All notable changes to `metric-autopsy` are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[SemVer](https://semver.org/). Downstream papers should cite a fixed version for
comparability.

## [Unreleased]

Version `0.3.0.dev0` — the verdict logic reworked against the exploratory probes in
`validation/probes/` (a development set with no confirmatory weight). **Breaking:** the
verdict vocabulary, GATE 1/2/5 semantics and several result fields changed.

### Changed
- **Four-field verdict.** `Autopsy` carries `metric_validity`, `design_adequacy`, `effect`
  and `replication` (`core.Assessment`: status, reason, flags, detail); `report.decide` is
  the single verdict rule used by the API, the CLI and the MCP server. Verdicts: SUPPORTED —
  replicated / SUPPORTED (provisional until replicated) / NOT SUPPORTED / NO DETECTABLE EFFECT
  / INCONCLUSIVE / UNIDENTIFIABLE / DEGENERATE METRIC, qualified by *parametric only* and
  *underpowered*. Metric validity PASS now needs a demonstrated response (a positive control
  that beats its null, or an injected signal); a random-number metric is no longer certified.
- **GATE 0** classifies each nuisance response as *bias* (the nuisance moves the null,
  reverses or inflates the signal), *attenuation* (the signal shrinks toward the null:
  reported, passed to design adequacy as a power check) or *level shift* (no material
  structure to classify it: reported). Rules are invariant to affine re-expressions of the
  metric (v0.1 divided by |baseline|). A metric that never varies is DEGENERATE. A detectable
  but immaterial null signal (< `tol` × the effect scale) no longer becomes the unit in which
  level metrics are judged. Decided 2026-10-08:
  - the null gives every gene, independently, a depth-matched draw from neighbouring cells (as
    GATE 5's positive control); the shuffle within depth bins made a pair coupled only
    through cell size look coupled above the null;
  - a bias blocks only where it matters for the claim: a depth bias that the declared
    correction removes between groups is reported, not failed; any other bias FAILs only if it
    and the lower bound of its 95% interval exceed `bias_tolerance` × SESOI (pre-registered,
    default 0.5); without a SESOI it is reported as unsized, SUPPORTED is withheld and the
    message says to declare a SESOI. At 800 cells per donor the earlier rule failed
    log-normalized Pearson on a truly coupled pair in 4 of 4 datasets (probe p13).
- **GATE 4 by the interval of its response (decided 2026-10-08).** 200 injections, each against
  its sham; PASS when the lower 95% bound of the mean response is above 0 in the declared
  direction, FAIL only when the upper bound is below `delta_min` (pre-registered; default 0.5 ×
  SESOI; `--delta-min` in the CLI, `delta_min` in MCP), UNTESTED (WARN) otherwise: "metric invalid"
  is a proof of blindness. An interval inside (0, `delta_min`) is FAIL. Where GATE 4 is UNTESTED a
  positive control on another pair no longer makes the metric valid for the claim. The former rule
  (z ≥ 3 over 10 injections) failed the valid metric where its response is real but weak (probe
  p14, now a dev-set regression). A coupling injection rewrites the two genes of one private copy
  in place (the same draws and values, about 20 times faster).
- **Verdict causes.** `report.decide_cause` returns the verdict and its cause (`report.CAUSES`;
  `cause` in the JSON report and the Markdown): which gate made the metric invalid, explained by
  depth, opposite direction, and so on.
- **Directional claims (decided 2026-10-08).** The pre-registration states the claimed
  `direction` of the change from `groups[0]` to `groups[1]` (`increase`, `decrease`, or
  `two-sided`). SUPPORTED needs the effect in that direction (two-sided test at alpha, so a null
  gives a false SUPPORTED at alpha/2; `test_a_directional_claim_is_tested_two_sided_at_alpha`
  pins it); a detected effect the other way is NOT SUPPORTED; without a direction the verdict
  stays INCONCLUSIVE; a non-directional claim is marked.
- **GATE 1 is a diagnostic.** Bootstrap intervals, Bonferroni across assessable strata and
  `min_cells`; a confident imbalance is WARN (handled by the correction), STOP only when no
  stratum contains both groups. (v0.1 flagged 82% of null datasets at 64 × 20 cells.)
- **GATE 2 is the estimand-dependent correction** (`effect.estimate_effect`): binomial
  thinning by one common ratio per stratum, to equal mean depth for a *composition* estimand
  and to equal spike-in capture for a *content* estimand; UNIDENTIFIABLE without spike-ins or
  without a declared estimand. (Per-cell quantile matching was tried first and biased the
  joint distribution of the genes a metric reads.)
  `gate2_ngenes_matching` is deprecated (n_genes is downstream of biology) and unused by
  `run_autopsy`.
- **Replicate-level inference with a graded rule** (`replicate_col`): ≥ 4 replicates per group
  — exact (or Monte Carlo, with MC error) permutation over replicates is the test, the t
  interval the estimate; 3 — the t interval is the test, flagged PARAMETRIC_ONLY, the
  permutation p reported as unable to reach alpha; ≤ 2 or no replicate unit —
  INSUFFICIENT_REPLICATION, no effect verdict. Nested and paired designs; partially crossed
  replicates are UNIDENTIFIABLE. NO DETECTABLE EFFECT needs a pre-registered SESOI (TOST) and a
  valid metric. A raw difference detected across replicates that vanishes after the correction
  is NOT SUPPORTED — explained by depth (or capture) — and, with a SESOI, only where the TOST
  shows the corrected effect smaller than the SESOI (else INCONCLUSIVE; journal D7); one whose
  sign the correction reverses (both detected) is INCONCLUSIVE.
- **Power.** GATE 0's attenuation under depth halving gives λ, the fraction of a construct-scale
  difference that survives at the analysed depth; the design is UNDERPOWERED when the MDE
  exceeds λ × SESOI. A pre-registration may state its SESOI on the metric's observed scale
  (`sesoi_scale: "observed"`, CLI `--sesoi-scale`, MCP `sesoi_scale`): λ = 1 for the TOST and the
  power check (journal D7; the default `construct` is the behaviour so far).
- **GATE 5** judges controls against empirical nulls: the negative control against
  expression-matched unrelated pairs, the positive control against depth-matched draws of
  itself; Bonferroni across strata. Only the negative control can FAIL the gate; a positive
  control that does not beat its null is WARN (absence of evidence), and the metric stays
  UNTESTED if it fires nowhere — unless the metric is shown blind: a coupling of the
  pre-registered dose (`positive_control_dose`, default 2.0) is injected into the silent
  control's genes, with their own coupling removed, and the control FAILs only if the upper
  bound of the metric's response is below `delta_min` (decided 2026-10-08, probe p15; the rule
  of 2026-10-07, power ≥ 0.8 of a reference detector, failed a valid metric whose control was
  coupled, but weakly, in 17 of 20 datasets). The positive control's self-null replaces gene b with a
  depth-matched, thinned draw from neighbouring cells (a shuffle within depth bins let a pair
  coupled only through depth pass). p values are rank-based Monte Carlo p values with a
  two-stage extension to the resolution alpha/K needs (no normal tail; null pairs drawn without
  replacement). `pos_min`/`neg_max` select the legacy band (not the default).
- **GATE 6** re-estimates the effect on the second dataset with the same estimand, correction,
  replicate rule and strata: REPLICATED / NOT_REPLICATED (equivalent to zero, or opposite sign)
  / INCONCLUSIVE.
- The CLI demo data gain a `mouse` column (4 per sex × age block, drawn from a separate RNG, so
  X is unchanged); the demo now runs with `--replicate-col mouse --estimand composition`.

### Added
- `gate4_signal_response` and `injected_signal.coupling` / `injected_signal.module`: response
  to a known construct change planted by binomial thinning, measured against a matched sham
  (`inject.sham`: the same thinning without the signal).
- `stats.clopper_pearson` / `fmt_rate`: exact binomial intervals for every reported rate.
- `validation/prereg/`: the confirmatory validation's protocol (`v1.md`, third round) and its
  code — `select_backgrounds.py` (the mechanical choice of the backgrounds from the CELLxGENE
  Census and the named B3/B4 sources), `panel.py` (design, the assignment by a 256-bit key, a pool
  of up to 8 gene pairs per expression level drawn by the key independently of the condition,
  truth generators including N8's per-cell variable capture, claim cards, the allowed (label,
  cause) outcomes by the truth about the metric and the data, canonical sha256 of datasets and
  cards; never imports the engine), `oracle.py` (per background and level: SESOI, the response
  curve to GATE 4's injection and its saturation dose, the truth about the metric on every pair;
  on B1 the key dose and Δ\*; establishability per condition, variant and pair), `beacon.py`
  (the key: a drand quicknet round named in the run tag before it exists, BLS-verified),
  `run_panel.py` (datasets built on the fly, one engine run per claim card, deterministic
  reports), `blind.py` (guard, shards, collection, and `verify`: anyone re-runs datasets of a
  finished run, checks the datasets' sha256 and compares the reports by what they say),
  `score.py` (S1-S7, per-level and per-truth rates, per-gate outcomes, design effects), `oc.py`
  (thresholds by one principle with nominals by cause; every error criterion a set of cells — S2
  per group of conditions, S4-S6 per stratum and cause — each of which must pass; S3 judged on
  every stratum at the principle's tier or the floor's, fixed before the key), `anchors.py`,
  `simulate.py`, `timing.py` (with the drop rule and the memory rule), `frozen.py` (the frozen
  files, the workflow's records only, one attempt) and `test_prereg.py`; every environment pinned
  as its whole closure; `.github/workflows/validation.yml` runs everything from the tag
  `v0.3.0-prereg` to the scores (results committed to `results/panel-v1`), and
  `.github/workflows/audit.yml` the flagship audit on v0.1.1.
- `validation/probes/p14_gate4_by_expression_level.py` (now a dev-set regression of GATE 4's
  interval rule), `p15_gate5_weak_positive_control.py` (GATE 5's power rule against the
  response rule on weak but valid controls) and `p16_variable_capture_n8.py` (the whole engine on
  N8's per-cell variable capture against N1).
- `validation/probes/verdicts_v03.py` scores errors (verdict outside the set allowed by the
  design, false SUPPORTED) and decisiveness (definite verdicts where the design makes the truth
  establishable) per case, with Clopper-Pearson intervals. `validation/probes/JOURNAL.md`
  records every change to a frozen dev-set expectation (D1: p11 split into p11a / p11b by
  design; D2: lock-in tests).
- `stats` (t and permutation tests without scipy, TOST, MDE), `equalize` (thinning),
  `provenance` (hashes, versions, run log), `effect`.
- **Provenance.** Every report carries `data_sha256`, `prereg_sha256`, a `claim_id`,
  versions and the seed; `Autopsy.to_json()` / `save_json()` (strict JSON). A JSON-lines run
  log counts attempts per claim: every run with a pre-registration is logged in every interface,
  the Python API included (`log_path=`, else `$METRIC_AUTOPSY_LOG`, else
  `metric_autopsy_runs.jsonl`; `"off"` disables); runs without one and the demo are not.
- CLI flags `--replicate-col`, `--estimand`, `--direction`, `--sesoi`, `--bias-tolerance`,
  `--min-replicates`, `--prereg`,
  `--inject-signal coupling`, `--seed`, `--json`, `--log`, `--no-log`; the same parameters on the
  MCP `autopsy_report` tool.
- `DenseMemoryWarning` (threshold `$METRIC_AUTOPSY_DENSE_WARN_GB`, default 2).
- Tests: `tests/test_v03_verdict.py` (graded rule, designs, estimands, power, provenance,
  GATE 4/5, `decide`, API/CLI/MCP parity); the dev-set tests in `validation/probes/` are
  collected by `pytest`. The v0.1.1 tests were revised to the new semantics, each keeping the
  truth it checked.

### Fixed
- **The module injection under pandas' copy-on-write.** `injected_signal.module` and its sham
  wrote the thinned library sizes into the array of the `total_counts` column, which pandas 3
  (copy-on-write) hands out read-only: GATE 4 with a module signal failed on data that carry
  `total_counts` (CI, pandas 3; the panel's environment pins pandas 2.3.3). They write into a copy.
- **GATE 5's controls count only for the metric they test.** `run_autopsy` took on trust that
  `pair_metric` bound to `gene_pair` is the judged metric: a metric blind to gene b given
  `norm_pearson`'s controls was certified and SUPPORTED (3 of 3 dev datasets), and a
  random-number metric passed metric validity. The controls now run only when
  `pair_metric(data, *gene_pair)` equals `metric(data)`; otherwise GATE 5 is SKIP with the two
  values. The CLI and the MCP server always bind both from the same function.
- **The MCP server starts with mcp 2.x.** `pip install "metric-autopsy[mcp]"` now resolves
  to mcp 2.x, where `mcp.server.fastmcp.FastMCP` was renamed to `mcp.server.mcpserver.MCPServer`;
  `metric-autopsy-mcp` failed at start-up for every fresh install. `build_server` supports both.
  The `dev` extra now includes `mcp` (Python ≥ 3.10), so CI exercises the MCP front door.
- `per_cell_qc` and the gates warn before densifying a large sparse matrix (probe p12).

### Known limitations
- Development results only: the fixes were developed against the probes that found the bugs
  (errors and decisiveness on the dev cases: `validation/probes/verdicts_v0.3.0.dev0.log`).
- With 4 replicates per group the exact permutation reaches p < 0.05 only at complete
  separation, so a pure depth artifact is often INCONCLUSIVE rather than "explained by depth"
  (p11a: NOT SUPPORTED in 3 of 10 splits of its cells into mice); with 20 mice per group it is
  explained in 40 of 40 datasets (p11b; dev-set journal D1).
- GATE 0's depth-matched null can sit slightly above the data for leverage-heavy raw metrics
  in small samples, because the deepest cells have few deeper neighbours; for a depth-only pair
  the signal above the null stays below 0.01 (`test_gate0_null_keeps_depth_so_a_depth_only_pair_shows_no_structure`).
- GATE 4 judges the response on the dataset at hand: its interval covers the injections, not
  the sampling of the dataset, so near `delta_min` the verdict varies between datasets (p14
  reports each pair's spread between datasets). A random-number metric is never shown blind by
  it: its noise keeps the interval wide, so it stays UNTESTED.
- "Explained by depth" is a diagnosis the design must be able to establish: in p11b's
  conditions it needs about 20 mice per group; on typical designs the tool protects against a
  false SUPPORTED but does not prove the artifact (`paper/drafts/limitations.md`).

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
