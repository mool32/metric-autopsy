# metric-autopsy — status card

> The single source of truth for where this project is. Update it at every stage transition.
> Lifecycle and rules: ../_meta/RESEARCH_FLOW.md

**Stage:** 0 Seed · 1 Pre-reg · **2 Execute** · 3 Verdict · 4 Write-up · 5 Publish · 6 Archive  ← current
*(Moved back from 5 Publish on 2026-10-07: the validator is under validation. See the validation plan below.)*
**One-liner:** A gate system — shipped as a Claude Code skill, a pip package, *and* an MCP server — that red-teams a computed single-cell metric to tell biological signal apart from QC/technical/mathematical artifacts.
**Started:** 2026-07-04   **Last update:** 2026-10-07

## Links
- GitHub: https://github.com/mool32/metric-autopsy (public) · CI green (3.9–3.12)
- Preprint: none yet — manuscript ready in `paper/manuscript.md`; bioRxiv posting pending
- Zenodo DOI: **10.5281/zenodo.21195679** (concept, resolves to latest) · v0.1.1 pending release
- Portfolio entry: pending — `mool32.github.io/_data/publications.yml` + `papers.bib`

## Origin
Direct descendant of `perceptual_modules/paper1/oscilatory/docs/metric_validation_checklist.md`
— the "work over the errors" written after ACP, β-heart, and MI-coupling each failed. This repo
turns that checklist into runnable behavior.

## Seed (stage 0)
- **Question:** Can the "compute → believe" default be replaced by a reusable "state commitments
  → red-team → then believe" workflow that catches QC/technical/mathematical artifacts *before* a
  claim is made — packaged so it reaches both the Claude-skill audience and the pip audience?
- **Why it matters:** The three failures cost weeks each and would have produced false papers.
  The checks cost ~1 hour. A tool that enforces them, with dual distribution, is high-leverage.
- **What would falsify it:** If the gates cannot separate a known-confounded metric (mi_3bin on a
  planted QC confound) from a known-clean one on the same data — i.e. if they neither catch the
  artifact nor pass the real signal. *(Status 2026-10-07, v0.1.1: **not met at the verdict
  level.** On the planted demo data `norm_pearson` also gets FAIL (GATE 1 is data-level; GATE
  5's negative control trips on a closure artifact). The separation holds only for GATE 0 in
  isolation. See `validation/probes/`. v0.3.0.dev0 (this branch, dev set only): on the demo,
  where biology is identical, neither metric is supported (both INCONCLUSIVE after depth
  thinning; `norm_pearson` metric validity PASS); on QC-matched data with a real coupling
  difference `norm_pearson` reaches the provisional SUPPORTED only with replicates, an estimand,
  a positive control and resolved judgment. Whether this holds outside the dev set is step 3.)*

## Design decisions (locked for v1)
- **Metric as a plugin.** Gates take `metric(data) -> float` as a black box; they know scRNA-seq
  QC, not your metric. Keeps focus on single-cell while accepting any metric.
- **Triple purpose, one engine.** Skill (`SKILL.md`) + pip package (`metric-autopsy`) + worked
  example (`examples/`). `src/metric_autopsy/` is the single source of truth.
- **Scripted in v1:** GATES 0, 1, 2, 3, 5 (1 & 2 are the crown jewels). GATE 6 runs when a second
  dataset is supplied; it is fully demonstrated only in the worked example. GATES 4 & 7 are
  judgment — the skill elicits them, they are not scripted.

## Self-test (this is a tool, not a hypothesis test)
- **Result (v0.1.1):** the regression suite passes (synthetic gate tests + 24-finding
  adversarial-audit regressions + SimpleData↔AnnData compatibility), and engine + docs passed an
  adversarial multi-agent review. A green suite shows the code runs as specified. It does
  **not** show the validator is right.
- **Self-probe (2026-10-07):** the exploratory probes in `validation/probes/` show the validator
  errs in both directions. A random-number metric gets PASS. Real biology that moves QC (Xist
  female>male, a sorted cell-cycle control, a proliferation shift) dies at GATE 1/2. Cell-level
  permutations certify mouse-to-mouse noise. Each failure is a strict-xfail test. These probes
  are a dev set with no confirmatory weight.
- **Rework (v0.3.0.dev0, step 2):** four-field verdict, replicate-level inference with the graded
  rule, estimand-dependent correction, empirical-null controls, hashes and a run log. All 18
  dev-set failures pass their frozen tests. Further miscalibrations were found and fixed while
  checking robustness (GATE 0 failed level metrics on an immaterial null signal; GATE 5 failed
  small strata on a degenerate null and on silent positive controls, and skewed metrics at many
  strata on a normal tail and on repeated null pairs). One dev-set expectation asked for more
  than its design can establish: p11 (4 vs 4 mice) is split by design into p11a (NOT SUPPORTED
  or INCONCLUSIVE) and p11b (20 vs 20 mice, explained by depth in 40/40), journal D1.
- **Step 2b (2026-10-07):** a silent positive control FAILs where the design had the power to
  show it; the positive control's null is depth-matched; the depth correction thins by one
  common ratio; a sign reversed by the correction is INCONCLUSIVE; GATE 4 is contrasted with a
  sham; every pre-registered run is logged. On the dev cases (`verdicts_v0.3.0.dev0.log`, at the
  current engine): false SUPPORTED 1/150 (the effect test's alpha), outside the allowed set 1/190,
  decisiveness 130/150 (147/150 before the third round's GATE 4/5 rules, which call a metric
  invalid only on a proof of blindness). Development results only.
- **Step 2c (decisions of 2026-10-08):** directional claims (SUPPORTED only in the
  pre-registered direction); GATE 0's null is depth-matched; a depth bias the declared
  correction removes no longer blocks, other biases block only beyond `bias_tolerance` x SESOI
  (probe p13 is the regression); the limitation of "explained by depth" is drafted in
  `paper/drafts/limitations.md`. Pre-registration of step 3: `validation/prereg/v1.md`.
- **Step 2d (second round of 2026-10-08; the first final draft was not approved):** the
  unsized-bias message says to declare a SESOI; a test pins that a directional claim is tested
  two-sided at alpha. The panel: gene pairs at three expression levels drawn by the key, N8
  (per-cell variable capture) as a fifth key null condition, allowed sets of the real effects by
  the true effect Δ\*, S1 as a whole (790 datasets per key condition), S3 on correct definite
  verdicts, a 128-bit key with a sha256 commitment, and the blind run in GitHub Actions from the
  tag with the key in a secret and no dataset stored (`blind.py`, `.github/workflows/panel.yml`).
  Probe p14 found that GATE 4 fails the valid metric where its response is weak or absent.
- **Step 2e (third round of 2026-10-08; the second final draft was not approved):** GATE 4 by
  the interval of its response (PASS on the lower 95% bound above 0, FAIL only on the upper bound
  below delta_min = 0.5 x SESOI, else UNTESTED; p14 a dev-set regression, journal D5); GATE 5's
  silent positive control by the same rule (p15: the power rule failed a weak but valid control in
  17/20 datasets; D6); verdicts carry their cause. The panel scores (label, cause) pairs against
  the truth about the metric on every pair (the oracle's population response to GATE 4's
  injection), with a new criterion S5 (false "metric invalid") and nominals by cause; the key is a
  drand round named in the run tag before it exists; the run is deterministic and automatic from
  the tag to the scores in GitHub Actions (`.github/workflows/validation.yml`), the backgrounds
  selected from the CELLxGENE Census by the rule of v1.md 3.1, the results and scores committed to
  `results/panel-v1`; the flagship audit on v0.1.1 runs in `audit.yml` on the same tag.
- **Step 2f (after the first independent review):** the review found two critical defects — the
  anchors' B3 stopped the pilot, and a validator that never says SUPPORTED passed S3 — and
  important ones (the engine's equivalence bounds narrower than the oracle's, the truth measured on
  the null applied to every condition, error routes without nominals, the pilot never run through
  its command lines, the run tied to the work branch). Fixed: the panel reads only B1/B2; S3 by
  stratum within the joint requirement; `sesoi_scale` and a TOST-backed "explained by depth" in
  the engine (journal D7); the truth per case; a sound validator's nominals measured by the pilot;
  the pool rule per candidate; the scripts pin the numerical environment; the dry run rehearses
  the whole pilot; the run's own branches. v1.md section 10 lists every change.
- **Step 2g (after the second independent review):** the review found two critical defects —
  bad validators passed where S3's real-effect stratum was only reported (scenario C), and
  duplicate gene symbols would have made every card of a background an engine error — and
  important ones (cross-CPU last-digit differences against a byte-for-byte `verify`, S3's judged
  strata decided after the key, error criteria outside the principle at small nominals, the
  decisiveness rule failing below 0.5, a missing B2 stopping the pilot, the sound validator's model
  off the engine's verdict order, targeted errors absorbed by S4). Fixed: S3 judges every stratum
  at the principle's tier or the floor's, fixed in pilot.json before the key; unique gene symbols;
  `verify` by what the reports say; rule nominals at least the designed sizes; B2 dropped with N7;
  the model in the engine's order with N3 and N8 outside S3's strata; new criteria S6 (false NO
  DETECTABLE EFFECT) and S7 (no engine error). The review predicts S1 failing on N8 (the engine's
  ~8% false SUPPORTED there in simulation), an open point for the owner.
- **Step 2h (after the third independent review):** the review found one critical defect — the
  pooled error criteria let errors concentrated on the real effects, on N4 or on N7 through — and
  important ones (one attempt enforced only at the start of the run, the sound validator's model
  off the engine on GATE 5's controls, N4/N5's GATE 4 odds and the oracle's equalization, N7's
  unbounded memory). Fixed: every error criterion is a set of cells (S2 per group of conditions,
  S4-S6 per stratum and cause), each of which must pass; one attempt per job with the attempt
  recorded and the results pushed without force; negative controls typical of GATE 5's null, N4/N5
  measured on their own datasets, the oracle equalizing as the engine does; N7 at most 24 mice and a
  memory rule before the key; every environment pinned as its whole closure. The development probe
  p16 measures N8 through the whole engine on simulated backgrounds: a false SUPPORTED in 18 of 160
  datasets (11%; N1 2 of 160 then, 3 of 160 at the fourth review's fixes), where S1 allows 3.7% —
  the open point for the owner stands.
- **Step 2i (after the fourth independent review):** no critical defect; important ones — NOT
  SUPPORTED against the direction on null data never bounded (a validator saying it on 16% of the
  null cards passed), S2's small groups passing targeted false SUPPORTED (6-12%), GATE 0's refusals
  outside the sound validator's model (a sound engine could fail S3), and E2/E3's Δ* taken at full
  depth. Fixed: a cell for the null's other tail; effect verdicts the engine's rules exclude (N4
  without the replicate unit, the constant) in S7; every S2 group at least 300 cards (N4, N6a, N7 at
  300 datasets: 6,400 datasets, 6,700 claim cards) with every cell's resolution printed; GATE 0's
  refusal share measured with the frozen engine before the key and refusing first in the model;
  Δ* of E2/E3 at half capture; thirteen smaller items (v1.md section 10). p16 re-run at the fixes:
  N8 still 18 of 160 false SUPPORTED, and 19 of 160 against the direction, which the new null cell
  counts. The fixes were not reviewed again (no critical finding); then the owner's «утверждаю».
  No tag yet.

## Data
- See DATASETS.md. All public (Tabula Muris Senis, human skin CELLxGENE). Nothing irreplaceable.

## Open threads / next
- [x] Worked-example notebook `examples/mi_coupling_tms/notebook.ipynb` — executed on synthetic
  `demo_data`, with committed outputs and figures. (Optional: re-run on real TMS after data pull.)
- [x] `references/gates.md` reworked into the preprint methods (`paper/manuscript.md`).
- [x] Public repo + `gh` metadata (description, homepage→DOI, topics incl. `tool`/`single-cell`).
- [x] Release v0.1.0 → Zenodo DOI (concept + version) → DOI badge; CI green across Python 3.9–3.12.
- [ ] Post preprint to bioRxiv (manuscript ready); then fill preprint DOI in README/CITATION/manuscript.
- [ ] Portfolio: add to `mool32.github.io/_data/publications.yml` + `papers.bib`.
- [ ] v1.1: turnkey GATE 6 second-platform replication, more example datasets.
- [ ] Validation plan (2026-10-07):
  0. preserve the dev set and set the status;
  1. audit the flagship TMS data;
  2. rework the verdict logic;
  3. confirmatory validation on frozen tag `v0.3.0-prereg` with a new, blind panel;
  4. external verdicts;
  5. paper untouched until step 3 is done.
