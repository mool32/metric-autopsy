# metric-autopsy — status card

> The single source of truth for where this project is. Update it at every stage transition.
> Lifecycle and rules: ../_meta/RESEARCH_FLOW.md

**Stage:** 0 Seed · 1 Pre-reg · 2 Execute · **3 Verdict** · 4 Write-up · 5 Publish · 6 Archive  ← current
*(Moved back from 5 Publish on 2026-10-07 to validate the validator. On 2026-10-10 v0.3.0-prereg passed the pre-registered validation v1, within its scope only, and the pass excludes only gross excesses of false SUPPORTED on the pairs where the metric works; see the validation plan below. The draft of the new manuscript, `paper/manuscript_v03.md`, is paused since 2026-10-11: not for posting, to be rewritten after v0.4, validation v2 and the external verdicts.)*
**One-liner:** A gate system — shipped as a Claude Code skill, a pip package, *and* an MCP server — that red-teams a computed single-cell metric to tell biological signal apart from QC/technical/mathematical artifacts.
**Started:** 2026-07-04   **Last update:** 2026-10-11

## Links
- GitHub: https://github.com/mool32/metric-autopsy (public) · CI green (3.9–3.12)
- Preprint: none yet. `paper/manuscript.md` is the superseded v0.1.1 text (not to be uploaded); the v0.3 draft, `paper/manuscript_v03.md`, is paused (not for posting; to be rewritten after v0.4, validation v2 and the external verdicts)
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
  counts. The fixes were not reviewed again (no critical finding).
- **Step 2j (the owner's fifth round, 2026-10-09; the fourth review's version was not approved):**
  the seven open items decided — N8 unchanged, with the prediction from p16 (S1 fails on N8) and the
  wording of both outcomes recorded in v1.md section 6, and v0.4's correction planned with a v2
  pre-registration; S7 split into rule violations (S7a, none) and crashes (S7b, 0.1% by the
  principle, each with its traceback); an unjudgeable S3 stratum still fails; E1–E3 doubled; explained
  by depth allowed on a blind metric where an artifact is planted; the five S2 groups at 790 (9,750
  datasets, 10,540 claim cards); GATE 0's refusal share in the model at most 0.20 on a valid metric,
  a stratum above it failing S3. The full pilot was run before the tag on simulated backgrounds of
  the panel's sizes on the Actions runner (`pilot-rehearsal.yml`, `validation/prereg/pilot_rehearsal.log`):
  in one job about 240 minutes, over the bound of 170, so it is split over parallel jobs (base, four
  parts of establishability, merge), with the same pilot.json; split, every job within the bound
  (45 min, 27-47 min, 63 min; the largest process 6.5 GB), and the select job's full extraction 456 s.
  The fifth and last independent review (at e519e2b; posted unedited on the PR) found no critical
  defect and the full pilot within its limits; its four important findings were fixed — an
  ambiguous metric's allowed set with the blind row's explained by depth (E2, E3), the reported
  rates' numerators on cards whose effect verdicts the engine excludes, S3 also per condition (a
  validator undecided on a whole condition passed its stratum), the beacon's newest round asked
  before the round is named and again after the push — and its six cosmetic ones (the joint
  probability over 100,000 simulations: 0.912, 0.913, 0.912; GATE 0's refusal pilot keeping only
  refusals and crashes). The fixes were not reviewed again (no critical finding).
- **Approved (2026-10-09):** the owner wrote «Утверждаю протокол v1.md на коммите 90a7718.» The
  annotated tag `v0.3.0-prereg` goes on 90a7718f77a46440ca11cd90150733a67f6b138d (tree
  f063c3ccc8b79b8e6be4b132921f7b613cdbf2c5); its message, in `validation/prereg/v0.3.0-prereg-tag.txt`,
  quotes the approval with its date and lists the sha256 of the 35 frozen files (`frozen.py`
  FROZEN_PATHS; v1.md's own: 41c6d8b55c6dc7573fe387eedad3601cf79a601f0e22fcd9f053bd17baa0886a). The
  agent's session cannot push tags (HTTP 403 from its proxy), so the owner pushes it: that push starts
  `validation.yml` and `audit.yml`, and from then on nobody decides anything.
- **Step 1 done (the flagship audit, run 2026-10-09):** `audit.yml` ran once on the tag (workflow
  run 37990113298; outputs on `results/flagship-audit`). The step-1 report,
  `validation/flagship_audit/REPORT.md`, finds:
  - v0.1.1's verdict on §4 reproduces: `mi_3bin` dies at GATE 0.
  - The design §4 describes does not hold. Its "old (20 months)" group is 24-month males and
    21-month females, all from the mammary gland. There are no spike-ins. The female arm and the
    sex-by-age interaction are not established per mouse.
  - v0.1.1 kills the same-data positive control (Xist, Y genes) at GATE 1.
  The manuscript stays untouched until step 3; the report's edits are the author's to decide.
- **Step 3 done (validation v1, run 2026-10-09/10):** `validation.yml` ran once from the tag
  (workflow run 37990113304, attempt 1; the key is drand quicknet round 32929913). The results,
  scores and anchors are on `results/panel-v1` (commits 3cb3c21, adc7d65, 445c0dd).
  - `scores.txt`: every criterion S1–S7b passes ("validation PASSES"). 9,750 datasets and
    10,540 cards, with 0 crashes and 0 rule violations.
  - The pre-registered prediction for N8 was not confirmed: 11 false SUPPORTED of 790, where S1
    allows 29.
  - Scope: scRNA-seq counts; one metric, `norm_pearson` on a gene pair, on the 8 most correlated
    pairs of each expression level; two backgrounds from the CELLxGENE Census
    (`validation/exploratory/census_metadata.log`, workflow run 38122733079):
    - B1: human oligodendrocytes, MSSM cohort, nuclei, 10x 3' v3. It carried every condition but
      N7: 9,750 of the 10,540 cards.
    - B2: mouse islet beta cells, cells, 10x 3' v2 (39 donor ids) and v3 (11). Its 50 donor ids
      are samples, 6 of them pools of both sexes; all single females are the NOD strain. It
      carried only the null N7 and the anchor R1.
    - Artifacts: N2, N3 and N8 (capture loss and entry-level dropout, of v1's binomial-thinning
      family, and per-cell capture, outside the engine's correction family). N2 did not test the
      correction of depth:
      the oracle's power to detect its raw difference was 0.21, 0.14 and 0.16 at the three levels,
      never the 0.9 that makes it informative (`validation/exploratory/v1_pilot_table.log`).
    - Real effects: no level had a key dose, so E1–E3 ran at the saturation dose (B1: 2.5, 2.75
      and 3.0). Of their 1,800 cards 206 had a valid metric, 189 of them establishable and 121
      at or above the SESOI; all 1,285 with a blind metric were establishable too (S3's invalid
      stratum: E1 668, E2 358, E3 259). At the high level Δ* (0.258) is about the SESOI (0.255).
    - Not tested: GATE 6 (no card had a second dataset) and the content estimand (only in R2).

    Nothing else is validated, `mi_3bin` included.
  - What the pass excludes (computed after the results; `validation/exploratory/v1_null_cells.log`).
    A false SUPPORTED needs a metric that responds, and in each key null of the pair metric only
    253–271 of the 790 cards had a valid or ambiguous metric. There S1's cell fails with
    probability 0.95 only at a false-SUPPORTED rate of 14.3–15.3%, and a validator wrong on 7.5% of
    them passes with probability 0.98–0.99. The pilot showed this mix before the key
    (`oc_pilot.log`). The pass excludes only gross excesses on the pairs where the metric works.
  - The truth about validity is circular. The oracle calls a metric valid when it responds to the
    injection GATE 4 itself uses, so S5 measures GATE 4's sampling error and could not find a
    probe built the wrong way; R1's was one.
  - By level and truth (pre-registered secondary outcomes, v1.md sections 5 and 7;
    `validation/exploratory/v1_by_level.log`):
    - All 2,458 cards with a valid metric are at the high level. At the medium and low levels
      the metric is blind (55 medium-level cards of N7 ambiguous), and no card there was
      SUPPORTED (0 of 6,988).
    - In the twelve high-level null cells, after Holm's correction over the twelve (decided after
      the results) no cell differs from the nominal 2.5%; without it, N3 f = 0.4's interval
      excludes it. The largest are N3 f = 0.4, 7/83 = 8.4% (95% CI 3.5–16.6%; Holm 0.058), and
      N8, 11/253 = 4.3% (2.2–7.6%; Holm 0.61). Both are hypotheses for v0.4, not measured
      defects.
    - Decided after the results: N8 on its cards with a valid metric 8/221 = 3.6%; the pure nulls
      N1, N5 and N7 pooled 18/812 = 2.2%; SUPPORTED or against the direction 22/253 = 8.7% on N8
      and 10/83 = 12.0% on N3 f = 0.4, at a nominal 5%. N7 permutes the labels over B2's donor
      ids: a null at the level of the id, not of the mouse.
    - The S4 part of the N8 prediction was not confirmed either: S4's null cell held 58 of 2,165
      (limit 102), 11 of them N8's.
  - Design effects above 1.5: S3's null and invalid strata (4.74, 5.62), S5's null cell (5.28)
    and S3's E1 condition (1.54). In the false-SUPPORTED cells at most 1.055, at S1's N5, where a
    doubled error passes with probability 0.051, just above the principle's 0.05.
  - GATE 0 refused 31 cards, all on N7 (26 at the high level). GATE 0 refuses where extra
    dropout would bias the metric beyond its tolerance, as it is built to, and a refusal is
    allowed, never definite; but on these null data the refusals cost decisiveness (S3's N7 cell
    204 of 237, against 250/271, 346/364 and 186/203 on N1, N2 and N5). 29 of the 31 fell on
    three of B2's gene pairs, where dropout shifted the correlation beyond the tolerance; eight
    tests found no association with B2's chemistry or depth (one at p = 0.049; decided after the
    results; `validation/exploratory/b2_depth_by_chemistry.log`).
  - Anchors: R1 did not run. `anchors.py` passes B2's sparse matrix to `SimpleData`, which
    accepts only dense arrays (ValueError). R2a–R2c are in their allowed sets; R3 was dropped
    with B4.
  - R1, once more, after the results and outside v1 (DEVIATIONS.md D1: the loader, and with it
    the cut to 2,000 genes, which changed the denominator of CP10k; workflow run 38033821978;
    `results/panel-v1-r1` at 36a04bf). It ran, and two of its three claims are outside their
    allowed sets. The causes, read after the outcome (`validation/exploratory/b2_sex_structure.log`):
    - Y genes (and the design of the Xist claim): UNIDENTIFIABLE, a correct reading of B2. Its 50
      donor ids are 34 single males, 10 single females (all NOD) and 6 pools of both sexes, and
      within a development stage the sexes meet only inside the pools (2 weeks; 20 months and
      over). The protocol erred: `r1_allowed` counted the ids of each sex over all stages, and
      checked neither that an id is one animal of one sex nor that both sexes share a stratum.
    - Xist: NOT SUPPORTED, "metric invalid" by GATE 4, which is false: an engine defect. GATE
      4's module probe raises Xist 2-fold in every cell and measures the metric over all cells.
      74.2% of B2's cells have no Xist, so the response is +0.163 over all cells, below
      delta_min 0.25, and +0.631 over the cells with Xist. With these constants the failure was
      predictable from B2's composition: for a mean over cells the response over all cells is
      the share of cells with the marker times the response in them, below ln 2 for a
      well-detected marker, so at Xist's 0.631 the probe needed at least 0.25 / 0.631 = 39.6% of
      cells with Xist (`validation/exploratory/gate4_module_analytic.log`). v1 never tested this path: its panel
      used the module probe only on the useless metric N6c. In the anchors it ran on four valid
      metrics and passed on three: R2a's G2M score, R1's Y-gene score, and R1's sham, which is
      the same Xist metric on the female cells only.
    - The sham, female against female: NO DETECTABLE EFFECT, which is correct.

    v0.1.1 did not confirm Xist because of GATE 1, on TMS; v0.3, whose GATE 1 is a diagnostic,
    was not run on TMS, and on B2 it did not confirm Xist either, because of GATE 4 and the
    design. v0.3 confirmed no real biological difference. Nothing is changed without the owner's
    decision.
  - Checked after the run, on another CPU:
    - `score.py` reproduces every criterion.
    - `blind.py verify` on the first 20 datasets: 20 of 20 datasets identical; all 22 reports
      the same (16 byte-identical, 6 within 1e-9).
  - One secondary number is fragile. `oc.error_rule_deff` switches to a beta-binomial for a
    design effect even 2e-16 above 1, and there scipy overflows and reads 1.0 / 1.0. In the run's
    `scores.json` (`secondary.s1_at_realized_deff`) that hit S1's N2 (design effect 1 + 1.6e-15);
    re-scored on another CPU, N8 too. At a design effect of 1 both are 0.982 / 0.046, N1's.

## Data
- See DATASETS.md. All public (Tabula Muris Senis, human skin CELLxGENE). Nothing irreplaceable.

## Open threads / next
- [x] Worked-example notebook `examples/mi_coupling_tms/notebook.ipynb` — executed on synthetic
  `demo_data`, with committed outputs and figures. (Optional: re-run on real TMS after data pull.)
- [x] `references/gates.md` reworked into the preprint methods (`paper/manuscript.md`).
- [x] Public repo + `gh` metadata (description, homepage→DOI, topics incl. `tool`/`single-cell`).
- [x] Release v0.1.0 → Zenodo DOI (concept + version) → DOI badge; CI green across Python 3.9–3.12.
- [ ] ~~Post preprint to bioRxiv (manuscript ready)~~ — superseded (2026-10-10): the v0.1.1
  manuscript is not to be uploaded (`paper/ARXIV_SUBMISSION.md`). The draft on v0.3,
  `paper/manuscript_v03.md`, is paused (2026-10-11): the tool first has to show itself on real
  tasks and against a donor-level test.
- [ ] Portfolio: add to `mool32.github.io/_data/publications.yml` + `papers.bib`.
- [ ] v1.1: turnkey GATE 6 second-platform replication, more example datasets.
- [ ] Validation plan (2026-10-07):
  0. preserve the dev set and set the status;
  1. audit the flagship TMS data;
  2. rework the verdict logic;
  3. confirmatory validation on frozen tag `v0.3.0-prereg` with a new, blind panel;
  4. external verdicts;
  5. paper untouched until step 3 is done.

  Status, 2026-10-11: steps 0–3 are done (step 3: every criterion of v1 passed, within its
  scope); step 4 is next. The main track is now the test on real claims
  (`validation/plans/real_claims.md`): if the tool shows no use on real tasks, it is not needed
  (the owner's decision of 2026-10-11).
- [ ] Lessons for v2, from v1's result and R1's (2026-10-11; computed after the results, from
  `validation/exploratory/v1_null_cells.log`, `v1_by_level.log`, `v1_pilot_table.log`,
  `b2_sex_structure.log` and `gate4_module_analytic.log`):
  - **Size the criteria on the cards at risk.** A false SUPPORTED needs a metric that responds,
    and on the real backgrounds the pair metric was blind at the medium and low levels. S1's 790
    cards per key null held only 253–271 at risk, so S1 could catch only a rate of about 15% among
    them with probability 0.95. Count the criteria of false SUPPORTED on the cards whose metric
    is valid or ambiguous, or per level; the pilot that showed the mix before the key could have
    resized them.
  - **A truth about validity that does not come from GATE 4.** v1's oracle called a metric valid
    when it responded to GATE 4's own injection, so the panel measured GATE 4's sampling error,
    not whether its probe is built right. R1 showed a probe built wrong.
  - **N2 must test the correction of depth, and the real effects must clear the SESOI.** N2's raw
    difference had power 0.14–0.21, never 0.9. No level had a key dose, so the real effects ran at
    the saturation dose, with Δ* about the SESOI at the high level (0.258 against 0.255).
  - **The targets of v0.4: N8, N3 f = 0.4 and GATE 4** (the owner's decisions of 2026-10-11 are in
    `validation/plans/v0.4_spec.md`: N8 gets only a diagnostic, its correction postponed to v0.5;
    N3's errors are met by a GATE 2 by which a correction can only remove a conclusion). Hypotheses
    with a mechanism, not measured defects (Holm 0.61 and 0.058 over the twelve high-level null
    cells):
    - Per-cell variable capture (N8): at the high level, 11/253 = 4.3% false SUPPORTED at a
      nominal 2.5%.
    - Entry-level dropout (N3 f = 0.4): 7/83 = 8.4%, and 10/83 = 12.0% with the answers against
      the direction. Found after the results by the decisive experiment (`v1_analyzers.log`,
      section 5): all 7 were made by the depth correction, which enlarged a raw difference that the
      replicate-level test did not detect; a test over donors without the correction gives 0 of 83.
    - GATE 4 must measure the response where the signal is claimed (within a group, or over the
      replicates), not over all cells.
    - Housekeeping with the next engine edits (src/ is frozen at the tag): the docstring of
      `src/metric_autopsy/qc.py` cites the internal checklist's "2.4x", a number without a
      script; remove it, with an entry in CHANGELOG.md.
    - The module probe's sham must create no response by losing counts. At a depth of 2,000
      counts the sham thins a marker of single counts that the injection keeps, so of a response
      of 0.360 the sham's loss makes 0.178 and the 2-fold change 0.182
      (`gate4_module_analytic.log`): on shallow data the probe can pass for the wrong reason.
  - **The development data erred in both directions.** The development probe on simulated data
    (p16, 11.3%) was 8.1 times N8's real rate over all levels and at most 3.7 times on the
    datasets at risk, and it ran with a SESOI of 0.1 and a dose of 2.0, where B1's high level had
    0.255 and 2.5, so simulation and settings are confounded. The
    development data did not show GATE 4's failure on a marker present in only some cells. So
    v2's development backgrounds must be real. Under v1.md section 6, v1's datasets are now
    development data.
  - **GATE 4's module probe measures the response over all cells.** For a mean over cells that
    is, by identity, the share of cells with the marker times the response in them, below ln 2 for a well-detected marker. The default (2-fold in 30% of cells) cannot pass a
    well-detected mean-expression metric at a SESOI of 0.5 (at most 0.3 × ln 2 = 0.208 against
    delta_min 0.25); R1's probe (frac 1.0) needed at least 39.6% of cells with Xist at its
    response of 0.631 in B2's Xist cells (`gate4_module_analytic.log`), and B2 has 25.8%. At a low depth
    single-count markers pass through the sham's loss as much as through the 2-fold change (a
    target of v0.4 above). Panel v2 includes valid metrics with the module probe, among them
    markers present in only some cells; in v1 it ran only on the useless metric N6c.
  - **The probe's dose and the SESOI.** In v1 both came from the pilot, that is, from the oracle.
    A user sets them, so v2 must test the engine with the defaults, or with the rule by which a
    user picks them.
  - **A positive control is checked before it is pre-registered, by code independent of the
    engine.** Each id is one animal of one sex, checked by Xist and the Y genes, and both sexes
    share a stratum; in B2, 6 of the 50 ids are pools of both sexes. The allowed set, and the
    checks of the data's structure behind it, are computed by independent code, as the oracle
    is: computed by the engine's own design rules (`effect._replicate_design`, GATE 1's
    stratum rule), the anchor would not test those rules. By those rules R1's design on B2 was
    UNIDENTIFIABLE, with the pools (ids partially crossed with sex) and without them (no stage
    holds both sexes).

## Requirements for v2
Collected on 2026-10-11 from the two reviews of the v0.3 draft, the lessons of v1 (above, with
their logs) and D1 (`validation/prereg/DEVIATIONS.md`). Review 1 is the review of the full
draft by an independent agent that had not seen the work on the text, and Review 2 its re-check of
0ed504c, both 2026-10-11, both relayed by the owner; they are not reviews by human experts. Until
their texts are in `paper/reviews/` with their provenance, they are cited by the owner's item
numbers in the decisions of that day. The
operational draft is `validation/plans/v2_requirements.md`.

**The criteria**
1. **Size the criteria of false SUPPORTED on the cards at risk** (a valid or ambiguous metric),
   per condition and level, and print for every cell the rate it catches with probability 0.95.
   *Source:* v1 lesson (S1's 790 cards held 253–271 at risk); Review 1 A1; owner, 2026-10-11.
2. **Pre-register every analysis a statement rests on**: the rates by level and truth, pooled
   nulls, the test of each high-level null cell against its nominal with the correction over
   the cells. Everything else is marked post hoc. *Source:* Review 1 A2, A3; Review 2 A3.
3. **Test each prediction with a cell that can confirm it**, sized to detect the excess v1
   measured (N8 4.3% and N3 at f = 0.4 8.4% at the high level, nominal 2.5%), and name every
   criterion a prediction binds (N8's S4 part). *Source:* Review 1 A3 and C; v1 lesson (the
   targets of v0.4).
4. **Design effects by donor and by gene pair**; every cell keeps the principle (a doubled error
   passes with probability ≤ 0.05) at its realized design effect. *Source:* Review 1 B11 (N5:
   0.051); owner, 2026-10-11 (the design effect by pair).
5. **Operating characteristics for a validator sound within the engine's correction family**,
   not with the oracle's true correction, so that what the criteria predict for an
   out-of-family artifact is explicit. *Source:* Review 2 B9.

**The truth**
6. **A truth about validity independent of GATE 4's mechanism**, among the candidates the same
   Tabula Muris Senis mice measured by FACS (Smart-seq2) and by droplet (10x), the deeper method
   as the truth, once code independent of the engine has shown the mice are really shared.
   *Source:* Review 1 B7; v1 lesson; owner, 2026-10-11.
7. **N2 tests the correction of depth, and the real effects clear the SESOI**: the run needs a
   level where the oracle detects N2's raw difference with power ≥ 0.9, and a key dose with Δ*
   clearly above the SESOI, or it stops before the key. *Source:* Review 1 A5; v1 lesson.
8. **No condition on one background only**, and each background's chemistry recorded and
   declared as a stratum where it is mixed. *Source:* Review 1 B6 (B2 carried only N7 and R1)
   and C (B2's 10x v2 and v3).

**The engine under test**
9. **The engine as a user runs it**: the SESOI, GATE 4's dose and delta_min from a rule a user
   can apply (v0.4), not from the oracle's pilot. *Source:* v1 lesson (the probe's dose and the
   SESOI); Review 1 (its abstract: "user-chosen settings"); owner, 2026-10-11.
10. **GATE 6 and the content estimand in the panel**: cards with a second dataset, and data with
    spike-ins. *Source:* Review 1 B12.
11. **Pairs beyond the 8 most correlated of each level**, among them random and weakly
    correlated pairs. *Source:* Review 1 B12; owner, 2026-10-11.
12. **Metrics of the types the real claims use** (`validation/plans/real_claims.md`), each with
    a probe v0.4 has: coordination (a pair correlation; a gene-set measure such as GCL with the
    set coupling probe of v0.4 change 9) and mean expression of markers present in only some
    cells or absent in one group, through the module probe. Noise and entropy, which need probes
    v0.4 does not have, move to the next version. *Source:* owner, 2026-10-11 (narrowed the same
    day); v1 lesson (the module probe ran only on N6c); D1 (R1's Xist).

**Controls and anchors**
13. **Real positive controls checked by independent code before registration**: each id one
    animal of one sex (by Xist and the Y genes), both sexes in one stratum, and the allowed set
    computed by code independent of the engine's design rules, as the oracle is. *Source:* v1
    lesson; Review 1 B8 and B14; D1 (the causes); owner, 2026-10-11.
14. **Every anchor runs end to end on its real data before the tag**, with the run's loader,
    gene cut and library size. *Source:* D1 (the sparse loader; the cut to 2,000 genes changed
    CP10k's denominator).
15. **Each probe checked by arithmetic against the data's composition before registration**
    (for the module probe: the share of cells with the marker times its response there, against
    delta_min). *Source:* Review 1 B8 and B10; D1.

**Development data**
16. **Real development backgrounds, at the confirmatory run's settings**, and every comparison
    with v1 on the same denominators. *Source:* v1 lesson (p16 erred in both directions);
    Review 1 A4; Review 2 A5.
17. **No outcome fixed by a test's discreteness**: the smallest attainable p of each design
    checked before registration (p07: an exact test at 3 against 3 cannot reach 0.05).
    *Source:* Review 1 B13.

**Comparators and reporting**
18. **The comparison with a donor-level test is v2's primary outcome**: every card also run by
    a test over donors without gates (with a TOST for absence) and by the engine without GATE 0,
    2, 4 and 5; false definite conclusions and sensitivity compared per stratum, with an exchange
    fixed in advance (proposed: the gates are worth it where the false conclusions they remove are
    at least a quarter of the true detections they lose). On v1's cards the exchange fails where
    the metric works and holds on blind and useless metrics (`v1_correction_rule.log`). *Source:*
    owner, 2026-10-11.
19. **Report the establishable cards per stratum, and the cards at risk next to every rate.**
    *Source:* Review 2 A1; Review 1 A1.
20. **Count the errors a correction makes**: on every artifact condition, the cards whose corrected
    effect is detected while the raw one is not. *Source:* the decisive experiment
    (`validation/exploratory/v1_analyzers.log`, section 5: 15 of the 21 false effects the gates
    added on v1's null cards at risk; all 7 of N3's at f = 0.4 on the high level).

The operational drafts: `validation/plans/v2_requirements.md` (v2) and
`validation/plans/v0.4_spec.md` (the engine changes, each with its v2 test; no code before the
owner approves it).
