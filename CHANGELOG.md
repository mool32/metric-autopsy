# Changelog

All notable changes to `metric-autopsy` are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[SemVer](https://semver.org/). Downstream papers should cite a fixed version for
comparability.

## [0.4.0.dev0] — in development, not released

The decided changes of `validation/plans/v0.4_spec.md` (the owner, 2026-10-11 to 2026-10-19):
1, 2, 3, 4, 6, 7, 8, 12, 13 (as a line of the report), 14, 15, 16, 17, 18 and 19, and the claim's domain
(the decisions of 2026-10-18, after the second diagnostics on the pilot's data,
`validation/exploratory/real_pilot_4.log`, the second decisions of 2026-10-18, after the third,
`validation/exploratory/real_pilot_5_summary.log`, and the decisions of 2026-10-19, after the fourth run,
`validation/exploratory/real_pilot_6_summary.md`). Validation v1's scope stays the tag
`v0.3.0-prereg`'s: nothing here is validated.

### Changed
- **The differential sensitivity of a metric of structure (change 17, the owner's decisions of
  2026-10-19).** GATE 0's premise that an attenuation cannot create an effect holds only where it is the
  same in both groups; on preprint §4 extra dropout removed 64% of `mi_3bin`'s signal while the young
  units detected 1,166 more genes per cell. For a coupling probe (`report.STRUCTURE_PROBES`; a metric of
  level is not checked) GATE 4 now injects the same change of the construct into the other group as
  well, dosed by its own biological variances, from a random stream of its own (GATE 4's numbers do not
  change); `gates._differential` gives each group's response k, their ratio k_poor / k_rich, the
  metric minus GATE 0's null of the claim's genes in the more sensitive group (M, 3 draws per unit), and
  A = M × (1 − k_poor / k_rich), the difference at an equal coupling (groups[0] − groups[1]), with a
  delta-method interval over the replicates (`detail["differential"]` of GATE 4). Where the groups differ
  in detection (`report.detection_gap`: genes per cell, the claim's genes left out, permutation p below
  alpha on the claim's cells), A lies on the claim's side and its bound there exceeds bias_tolerance ×
  SESOI, the effect is judged less that bound — the same estimate and test on shifted values
  (`report._shifted`) — and a verdict that does not pass is NOT SUPPORTED, cause
  `differential_sensitivity` (`report.differential_sensitivity`, `effect.detail["differential_sensitivity"]`);
  where A cannot be sized there, a detected effect is INCONCLUSIVE with the same cause. A line of the
  dossier always. GATE 0's docstring and `references/gates.md` corrected.
- **One direction for every difference of the dossier (change 18).** The line of `norm_pearson` beside
  the metric gave groups[1] − groups[0]; it now gives groups[0] − groups[1], as every other line, and the
  dossier's header says so (`differences` in the JSON report). The specificity line gives its t and its
  controls' in the same direction and names the claim's side.
- **The full dossier after the first FAIL (the owner's second decisions of 2026-10-18).** The verdict
  follows the first FAIL as before, but the effect — within the domain, the mixture's, outside the
  domain, per type — GATE 6, the composition, the technical covariates and every other report line
  are estimated after it too. The effect then carries the flag `DOSSIER_ONLY` (and
  `detail["dossier_only"]`), its reason opens with "(the dossier: estimated after metric validity's
  FAIL, which decides the verdict)", and `decide` never reads it before the FAIL. `stop_on_first_fail`
  (the CLI's `--no-stop`, the MCP tools' argument) is kept for compatibility and recorded; it no
  longer stops anything.
- **GATE 0 blocks only a bias it cannot show to be equal between the groups (change 12, decided
  2026-10-12, coded on the owner's second decisions of 2026-10-18).** Where GATE 0 FAILs on a blocking
  bias and the data have a replicate unit, `run_autopsy` estimates the groups' technical difference
  (`report.technical_difference`): the claim's metric on GATE 0's null of the claim's genes (each
  gene's counts redrawn independently from cells at least as deep within the unit, `gates.claim_null`:
  no gene-gene biology, each group's own technology), averaged over `TECHNICAL_DRAWS` = 3 draws per
  unit, estimated between the groups as the effect is (the same units, domain, correction and test).
  Where its 90% interval lies within ±GATE 0's tolerance (bias_tolerance × the SESOI in the metric's
  units, the TOST) the bias is reported, handling "equal between groups", GATE 0 WARN; otherwise it
  still blocks and the test is recorded in GATE 0's message and `detail["between_groups"]`
  (`gates.gate0_between_groups`). A direct call of `gate0_independence` has no groups and is
  unchanged. Without a replicate unit nothing changes.
- **`norm_pearson` gives NaN, not 0.0, where it measures nothing (the owner's second decisions of
  2026-10-18).** Below `metrics.MIN_CODETECTED` = 20 cells with counts of both genes (never fewer than
  3), or where one gene takes one value among them (tested by equality: the standard deviation of
  equal floats can round to 1e-16, and their correlation to 1.0). The minimum, fixed from principle
  before any run, is the 20 cells a unit needs to count (`composition.MIN_CELLS`), applied to the
  cells the correlation is computed on (Fisher's z over 20 cells has a standard error of 0.24). v0.3
  returned 0.0, a value not measured reported as no correlation.
- **Replicates the metric gives no value on are left out and listed.** A replicate (or a unit
  (replicate, type) of the domain) whose metric is not finite, uncorrected or at equal depth, leaves
  both tests and is listed in `detail["unmeasured"]` (and the domain's `unmeasured`) and in the
  Markdown report's section "Left out: the metric gives no value (not measured)". Where that
  leaves fewer replicates per group than the pre-registered minimum, which there were before, the
  effect is NOT_ESTIMABLE and the verdict UNIDENTIFIABLE, cause `unmeasured`.
- **The SESOI on the construct's scale, brought into the metric's units by GATE 4 (change 3,
  decided 2026-10-18).** `prereg["sesoi"]` is by default a change of the construct
  (`sesoi_scale="construct"`): for a coupling probe a change of the correlation of the pair's log
  expected expression, in (0, 1); for a module probe a change of the module's natural-log
  expression. GATE 4's probe is dosed to that change (`injected_signal.at_construct_change`): a
  module lowered exp(SESOI)-fold; a coupling at the strength that raises the correlation by the
  SESOI given the two genes' biological variances in the probed cells (`biological_log_variance`, a
  Poisson-lognormal method of moments, floored at `BIO_VAR_FLOOR` = 0.01; `coupling_strength_for`).
  Its response s_m, with a 95% interval over the replicates (the injection repeated within each,
  a cluster-robust mean), is the SESOI in the metric's units for GATE 0's tolerance, the TOST and
  EFFECT BELOW SESOI. A SESOI on the metric's own scale is declared with `sesoi_scale="observed"`
  and used as declared; a SESOI on the construct's scale that no probe translates is not used
  (SESOI_NOT_TRANSLATED), and the design's floor stands in.
- **GATE 4 decides only blindness (change 3).** PASS where the lower 95% bound of s_m lies above
  zero in the injection's direction, FAIL (blind) where the upper bound lies at or below zero,
  UNTESTED (WARN) otherwise. `delta_min` is recorded, not used (the CLI's `--delta-min` and the MCP
  tool's `delta_min` are accepted and ignored). The detail gives the dose–response at 0.5×, 1× and
  2× the SESOI (`GATE4_CURVE`; `GATE4_N_REP_CURVE` = 20 injections per replicate at the extra
  doses). The protection that delta_min gave is carried by the checks that explain a difference
  away (depth, change 6; the domain and composition, change 7).
- **The design's floor (change 4).** se × (t(1 − α, df) + t(0.975, df)) of the replicate-level
  estimate (`effect.design_floor`): the smallest SESOI at which this design confirms absence by
  the TOST with power 0.95 at a true difference of 0 (v1's rule). A SESOI in the metric's units
  below it is flagged SESOI_BELOW_FLOOR; without a SESOI the floor is the margin of absence, flagged
  NO_SESOI_FLOOR, with "absence is confirmed only for effects of at least X". EFFECT BELOW SESOI and
  the TOST that explains a difference by depth need a pre-registered SESOI.
- **The claim's domain (decided 2026-10-18).** With `composition` (the cell type first), the claim
  is judged in the cell types where its genes are expressed: for a pair, both genes detected in at
  least 10% of the type's cells **in each group** (`domain.DOMAIN_DETECTION`, the lower bound of
  v1's low level; the owner's second decisions of 2026-10-18, by the definition of coupling: a pair's
  coupling is not defined in a group where a gene is not expressed — the third diagnostics ran with
  "in at least one group"); for a set, any of its genes in at least one group (a change of expression
  may be its appearance); for a claim without genes, every type; each type
  with at least 20 cells in at least 3 replicates of each group, per stratum
  (`domain.claim_domain`). The verdict follows the combined estimate within the domain: the metric
  of each (replicate, type) unit, at equal depth within its (stratum, type), the types weighted
  equally within a stratum and the strata by precision, with a cluster-robust interval over the
  replicates (Satterthwaite df) and relabelings of the replicates within strata (`combined.py`).
  The whole mixture's estimate, the estimate outside the domain and each type's are reported beside
  the difference of composition. An empty domain gives INCONCLUSIVE, cause `empty_domain`.
- **GATE 4 and GATE 0 run on the domain's cells** (decided 2026-10-18), not on the whole mixture.
- **GATE 5's negative control is reported, not failed.** An unrelated pair outside its null of
  unrelated pairs is reported (`negative_control_outside_null`); GATE 5 fails only a positive
  control that GATE 4's rule shows blind, or one in the legacy band. The 10 pairs with an
  equivalence or Holm rule of the spec's proposal are not built (change 13, decided 2026-10-17).
- The verdicts per type of change 7 (`report.per_type_verdicts`) give way to the domain's estimates
  per type; a type absent from every replicate of one group cannot be estimated and is listed with
  the group it is absent from (`absent_from`).
- **GATE 4 on the claim's own statistic (change 1).** With a replicate unit the signal goes into
  the group whose replicates have more of the probed genes' counts per cell, whatever the claim's
  wording (on an exact tie into each group, the weaker response judged); the response is the
  change of the replicate-level difference that the effect test uses, with its interval over the
  injected group's replicates (`GATE4_N_REP_REPLICATE` = 50 injections within each), judged
  against the injection's own direction (`inject.direction`). v0.3's all-cell response stays in
  the detail (`all_cells`) and is the judged response without a replicate unit; a direct call
  without the claim's groups is v0.3's gate.
- **Keyed draws (change 1).** Every gate's random stream is seeded by the data's hash, the claim's
  canonical form (`provenance.canonical_claim`: the metric, the genes, the groups sorted and the
  direction for that order) and `seed` (`provenance.claim_seed`): a mirrored wording gets the
  same draws and verdict, two claims on the same data different draws. The tie rules of GATE 1's
  bootstrap and of the thinning follow the groups' labels, not their order. Numbers that depend on
  draws differ from v0.3's for the same `seed`.
- **The data hash** (`provenance.sha256_data`) reads the values whatever their storage (a canonical
  CSR of float64 values): a SimpleData and an AnnData with a sparse X of the same data have one
  hash. Hashes differ from v0.3's, so the run log counts attempts afresh.
- **The module probe lowers the module by thinning its own genes (change 2)**, `fold`-times in a
  random `frac` of cells (all by default), against the untouched data and without a sham; it
  names its genes, its direction (a decrease) and the share of counts it removes
  (`inject.depth_loss`). v0.3 raised the module by thinning every other gene against a sham.
- **A correction can only remove a conclusion (change 6).** An effect counts as detected only if
  the raw replicate-level test detects it too, in the same direction; otherwise INCONCLUSIVE, cause
  `detected_only_after_correction`, and GATE 2 warns. The TOST also runs on the uncorrected data:
  a NO DETECTABLE EFFECT that only the corrected data show carries ABSENCE_ONLY_AFTER_CORRECTION.
  The report checks the data against the capture model (`equalization.capture_check`).
- **EFFECT BELOW SESOI (change 14):** a detected effect in the claimed direction whose (1 − 2α)
  interval lies within the effective SESOI; cause `effect_below_sesoi`. It comes after the
  replication check and before the judgment gates.

### Added
- **Composition (change 7, the decided minimum):** `run_autopsy(composition=[...])`, the CLI's
  `--composition` and the MCP tool's `composition` name obs columns of an annotation (cell type
  first): their shares per replicate are compared between the groups (the total variation
  distance and each level's difference, with permutation p over the replicates;
  `composition.composition_difference`), the claim is judged within each cell type with at least
  20 cells in at least 3 replicates of each group, the verdicts not combined
  (`report.per_type_verdicts`; a type absent from every replicate of one group is UNIDENTIFIABLE),
  and a difference of composition is noted in the verdict. Composition removes no conclusion:
  by which rule it would is not decided (on section 4 the drop of MI stayed within types).
- **The engine's `norm_pearson` beside a pair's metric (change 15)**, by stratum, on the
  uncorrected data (`Autopsy.reference`, `report.coupling_reference`); reported, not a gate.
- **Units whose ids name several animals (change 16)** are flagged COMPOSITE_UNITS, with the
  units they overlap (`effect.composite_units`); a warning, not an exclusion.
- `stats.nested_relabelings` and `stats.relabeling_p`: the effect test's relabelings for a
  statistic of any shape.
- **The check of specificity, a line of the report (change 13, decided 2026-10-18).** For a claim
  about given genes (a pair or a set), K = 1,000 controls drawn from every gene of
  `specificity_pool` (the same cells with all their genes; by default the data), matched on the log10
  mean and the detection, without the pairs of `prereg["known_interactions"]` and the genes of
  `prereg["known_to_differ"]`; each control's donor-level t within the domain, on the uncorrected
  data, as the claim's. The line gives how many of the K lie at least as far as the claim in its
  direction (with the share's Clopper–Pearson interval and an outcome: specific, not specific,
  undecided where the interval crosses α), the 10 nearest controls and every control at least as far,
  and the controls' common shift (their median t, with an interval by relabelings of the replicates
  under a shift model), with a warning where the shift is significant in the claim's direction. It
  is not a condition of SUPPORTED (`specificity.check`; `run_autopsy(gene_metric=, specificity=,
  specificity_pool=)`).
- `combined.py` (the units, the cluster-robust estimator and its relabelings), `domain.py` (the
  domain's rule), `specificity.py`.
- **A test of GATE 0's false FAIL rate** over 40 seeds against a target of 5% fixed before measuring
  (`tests/test_v04_domain.py`), xfail(strict=True) while the known issue below stands; and the same
  rate with change 12 applied, xfail(strict=True) too.
- **The groups' technical covariates per replicate, a line of the dossier (the owner's second
  decisions of 2026-10-18).** Per cell its total count and its number of genes with a count, over
  every gene of `specificity_pool` (else of the data); per replicate the median over its cells; the
  difference groups[0] − groups[1] of the replicates' means per stratum and over the strata, with
  the effect test's relabelings (exact up to 20,000), on all the comparison's cells and on the
  domain's (`composition.technical_covariates`, `Autopsy.technical`). GATE 1's ratio of median genes
  passed §4 in males at 1.32 while the young units detected 1,298 more genes per cell (exact p 0.016).
- **The claim genes' counts per type and group** (`index` of each type of the domain: the mean over
  replicates of log1p(10^4 × the genes' counts / the cell's total)), in the domain's table, and a
  line for every type outside a pair's domain where the genes reach the threshold in one group
  only (`domain["one_group"]`): counts from elsewhere, or the expression's appearance; a line of the
  dossier, not part of the rule.
- `tests/test_v04_dossier.py`: the regression tests of the second decisions of 2026-10-18.
- **The correction's line (change 19).** How equal depth changed the groups' relative gap in genes
  detected per cell, on the domain's cells and on all cells, with a WARNING where the gap grew
  (`Autopsy.capture`, `report.capture_summary`); on §4 it grew from 24.4% to 29.1%.
- `composition.technical_covariates(exclude=)`: genes left out of both counts.
- `tests/test_v04_differential.py`: the acceptance of the differential sensitivity on synthetic data
  (an equal coupling with a gap of detection not SUPPORTED; an equal detection unchanged; a real decrease
  SUPPORTED, with a gap where the effect less the bound passes) and the tests of the dossier's direction
  and of the correction's line; `validation/exploratory/differential_sensitivity.py` and its log.

### Known issues
- **The differential sensitivity over-corrects a metric that is not linear in the coupling.** `mi_3bin`
  responds more where the coupling is already larger, so with a real difference of coupling its ratio
  of responses also takes the slope of its curve and A takes a part of the true difference: on the
  synthetic data a decrease of 3, 5 or 7 × SESOI with a gap of detection as on §4 does not pass
  (`validation/exploratory/differential_sensitivity.log`, sections 2 and 4). Conservative: power is lost,
  no false SUPPORTED is added. Measured at a common coupling (GATE 0's null) the confound goes, but
  `mi_3bin`'s response there is of second order and too noisy.
- **The differential sensitivity costs power wherever the groups differ in detection**: the effect is
  judged less the bound of A's interval, whose width grows with the metric's distance from its null. With
  the gap as on §4 `norm_pearson`'s real decrease of 3 × SESOI does not pass either; 5 × does.
- **The differential sensitivity is measured on the uncorrected data**, as GATE 4 is, while the effect is
  judged at equal depth, where the gap in detection can be larger (the correction's line warns).
- **The differential sensitivity's condition reads a difference in detection by its p, not its size**: at
  an equal detection it opens at alpha's rate, and the groups' responses then differ by their biology
  and the metric's curve as well.
- **GATE 0 falsely fails a valid metric in 4 of 40 seeds, with change 12 too.** On the nested-mice
  test data (`norm_pearson`, groups without a difference in dropout, SESOI 0.02 on the metric's
  scale), the `extra_dropout` perturbation moves the metric by 0.0195 to 0.0213 (about 0.021) against
  a tolerance of 0.01 (0.5 × SESOI), and GATE 0 FAILs in seeds 8, 11, 36 and 37. The power test that
  met it is isolated from GATE 0's dropout (`test_underpowered_design_is_flagged_and_named_in_the_verdict`);
  the rate's own test records the failure. Change 12 does not lift it: with 4 mice per group the
  groups' technical difference cannot be shown within ±0.01, its 90% interval spanning about 0.25
  (seed 8: −0.081 [−0.210, +0.049]); the second rate test records it.
- **`codetected_spearman` and `pearson` still return 0.0 where they measure nothing** (fewer than 3
  co-detected cells, or no variance), and `mi_3bin` 0.0 on no cells: the correctness error decided
  for `norm_pearson`, not changed for them in this round.
- **v1's panel copies the tag's rules.** `validation/prereg/test_prereg.py` checks the panel's and the
  oracle's copies against the tag's engine; on this branch it checks that v0.4 differs from them in
  exactly the decided ways (the causes of changes 6, 14, the domain and the unmeasured replicates;
  GATE 4's FAIL only at or below zero; the module probe's default share and its direction).

### Fixed
- **GATE 5's note of the null's centre** printed "nan..nan" (with an All-NaN warning) where a metric gives
  no value in a stratum's null; the strata without a centre are left out of the range.
- **Two expectations of the dev set that encoded v0.3's mechanics** (`validation/probes/JOURNAL.md`,
  D8 and D9; the owner's second decisions of 2026-10-18). p09's null without a SESOI: INCONCLUSIVE was
  v0.3's way of claiming no absence without a margin; it is now NO_DETECTABLE_EFFECT at the design's
  floor, the floor named, the verdict INCONCLUSIVE while the metric's response is untested. p14's low
  level: "FAIL in at least 7 of 8" rested on v0.3's threshold 0.8 × delta_min; by v0.4's definition of
  blindness two of its pairs are blind and one responds, and GATE 4 PASSes none of the blind (0 of 10
  datasets) and FAILs none of the one that responds (0 of 5; `validation/probes/d8_p14_low_level.log`).
  Both tests are updated, their xfails removed.
- **`norm_pearson`'s 0.0 where it measured nothing** (above, Changed). On the third diagnostics'
  Marsh primary contrast (`validation/exploratory/real_pilot_5_summary.log`) all three ENZ-INHIB mice
  took this value against 0.336 to 0.394 in ENZ-NONE, so that contrast's SUPPORTED compared measured
  correlations with the fallback; by the domain's rule in each group it is outside the domain, not a
  claim about coupling, and Marsh et al. 2022 is no longer a control for coupling.

### Removed
- **The internal checklist's "2.4x" in `qc.py`'s docstring (change 8)**, a number without a
  script; the docstring points to the flagship audit's report instead.

## [Unreleased]

The record of validation v1; the release is postponed until the check on real data (the owner's
decision of 2026-10-15; draft notes in `RELEASE_NOTES_v0.3.0.md`). The code is the tag
`v0.3.0-prereg`'s and reports its version as `0.3.0.dev0`: the verdict logic reworked against the
exploratory probes in `validation/probes/` (a development set with no confirmatory weight).
**Breaking:** the verdict vocabulary, GATE 1/2/5 semantics and several result fields changed.

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
- `validation/prereg/`: the confirmatory validation's protocol (`v1.md`, fifth round) and its
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
  `score.py` (S1-S7, per-level and per-truth rates, per-gate outcomes, design effects, every
  crash with its traceback), `oc.py` (thresholds by one principle with nominals by cause; every
  error criterion a set of cells — S2 per group of conditions, 790 cards each, S4-S6 per stratum
  and cause, NOT SUPPORTED against the direction on null data bounded as a false detection — each
  of which must pass, with its resolution printed; S7a, no violation of the engine's deterministic
  rules (an effect verdict where they give none among them); S7b, crashes at a nominal of 0.1%; S3
  judged on every stratum at the principle's tier or the floor's, fixed before the key, and on
  every condition of a stratum (a validator undecided on a whole condition fails), with GATE
  0's refusals in the sound validator's model, at most 0.20 on a valid metric, above which the
  stratum fails; the joint probability over 100,000 simulations with its standard error), `anchors.py`, `simulate.py` (also backgrounds of the panel's sizes for the
  pilot's rehearsal), `timing.py` (with the drop rule, the memory rule and GATE 0's refusal shares
  measured before the key, only the refusals and the crashes kept), `frozen.py` (the frozen
  files, the workflow's records only, one attempt) and `test_prereg.py`; every environment pinned
  as its whole closure; `.github/workflows/validation.yml` runs everything from the tag
  `v0.3.0-prereg` to the scores (results committed to `results/panel-v1`; the pilot split over
  parallel jobs, `oracle.py --stage base | part | merge`, with the same pilot.json),
  `.github/workflows/pilot-rehearsal.yml` the full pilot before the tag on simulated backgrounds of
  the panel's sizes and the selection's full extraction, timed (its run's lines in
  `validation/prereg/pilot_rehearsal.log`: split, every job within 63 minutes, under the bound of
  170), and `.github/workflows/audit.yml`
  the flagship audit on v0.1.1.
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
- `validation/flagship_audit/REPORT.md`: the step-1 report on the flagship audit's one run on
  real data (workflow run 37990113298, outputs on `results/flagship-audit`). It decides which
  claims of preprint §4 hold and lists the edits for the author; the manuscript is not edited.

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
