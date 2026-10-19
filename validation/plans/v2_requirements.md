# Validation v2: requirements (draft for the owner)

A draft, not a pre-registration. It makes operational the requirements collected in `PROJECT.md`
("Requirements for v2", items 1–21, each with its source) and the owner's decisions of 2026-10-11
and 2026-10-12. Review 1 and Review 2 are the owner's paraphrase of two reviews received from an AI
assistant on 2026-10-11; their texts are not archived.
The engine v2 tests is v0.4 (`validation/plans/v0.4_spec.md`); each of its changes names its test
here. Numbers come from the logs named, each from a script in this repository. v1's datasets are
development data now (v1.md section 6), so v2's panel is drawn afresh, from backgrounds v1 did not
use.

v2 is the simulation track. Whether the tool is needed at all is decided on real claims
(`validation/plans/real_claims.md`, the main track, the owner's decision of 2026-10-11); v2 measures
what the engine does where the truth is known, and the real track says how often that matters.

## 1. The primary comparison: the engine against a donor-level test

The main question is what the gates add over a careful test at the level of the donors, so the
comparison is v2's primary outcome, not a secondary one (the owner, 2026-10-11; requirement 18).
Every card of v2 is run through (b) a test over donors without gates, here with a TOST at the SESOI
for its claims of absence and the outcome EFFECT BELOW SESOI (below), so that it can make the same
errors as the engine; (c) the engine with GATE 0, 2, 4 and 5 off; and the engine, v0.4. (a), the
test over cells, is reported, not compared.

**Errors:** false definite conclusions, counted per null card and per effect card apart: a claim of
an effect, in either direction, on null data or against the true direction on a real effect; a claim
of absence (NO DETECTABLE EFFECT) on a real effect with |Δ\*| ≥ SESOI or on a blind or useless
metric. Counted separately on the cards at risk (a valid or ambiguous metric) and on blind and
useless metrics. A claim of absence on null data with a blind or useless metric is true of the data
but not shown by the metric; it is reported next to the errors, not among them.

**EFFECT BELOW SESOI** (v0.4 change 14; the owner, 2026-10-17): a detected effect whose 90%
interval lies inside ± the (effective) SESOI. It is decided after change 6's rule: an effect that
rule does not support is INCONCLUSIVE, not EFFECT BELOW SESOI. Scored on null data as no false
effect; on a real effect with |Δ\*| ≥ SESOI as a false conclusion; on a real effect with
|Δ\*| < SESOI as correct. **(b) gets the same rule:** its TOST, and EFFECT BELOW SESOI where it
detects an effect in the claimed direction whose 90% interval lies inside ±SESOI, so that the
comparison measures what the gates add, not a difference in the logic of the verdict.

**Sensitivity:** effects found in the true direction per effect card at risk, and on the valid
cards with |Δ\*| ≥ SESOI.

**The exchange** (the owner, 2026-10-12): one false conclusion is worth four lost detections, the
convention α = 0.05 against β = 0.20. Errors removed per 1,000 cards of a stratum against
detections lost per 1,000 effect cards mix two denominators, so such a comparison depends on the
shares of nulls and effects in the panel, which is a decision of design. Instead, with a share π
of true effects among the claims, the gates pay off against a comparator where
(1 − π)·r + π·r₁ ≥ π·l/4: r the false conclusions they remove per null card, r₁ per effect card,
l the true detections they lose per effect card. (r₁, accepted by the owner on 2026-10-12, counts
the false conclusions a comparator draws on effect cards, (b)'s effects in the wrong direction and
(c)'s claims of absence; on v1 it is 0 against (b) where the metric works.)

**Primary outcome: the threshold π\* = r / (r + l/4 − r₁)**, the largest share of true effects at
which the gates pay off, per stratum, against (b) and against (c), with its 95% interval over cards
and the design effects by donor and by gene pair. On real data the same threshold is measured on
controls with the truth by construction (π\*_real) and compared with the share of true effects
among the real claims, estimated from their independent verdicts (`real_claims.md`); v2's π\* is
compared with π\*_real, as a check of the simulation.

**On v1's cards** (an illustration after the results, with (c) standing in for (b) with a TOST;
`v1_correction_rule.log`, section 4; intervals by bootstrap over cards, then over gene pairs, the
donors that datasets share not resampled):
- At risk (2,445 null and 515 effect cards), against (b): the engine removes 40 false conclusions
  (r = 40/2,445) and loses 96 detections (l = 96/515), π\* = 0.26 (0.14–0.37; 0.00–0.56); under
  the rule of v0.4 change 6, 58 and 99, π\* = 0.33 (0.24–0.41; 0.02–0.61). At v1's own share of
  effects (515 of 2,960, 0.17) the exchange holds: 40 ≥ 96/4 = 24.
- On blind and useless metrics (6,295 null and 1,285 effect cards), against (b): 525 false
  conclusions removed on null cards and 21 on effect cards, 136 detections lost, π\* = 0.89
  (0.82–0.98; 0.78–0.98). Against (c), with 596 false conclusions on the effect cards (575 of them
  claims of absence), the gates pay off at every share.
- The interval over gene pairs rests on 18 pairs at risk: π\* is decisive only with many more pairs
  (section 4).

## 2. The criteria count the cards at risk

v1's criteria of false SUPPORTED were counted over all of a condition's cards, and most of them
had a blind metric: in each key null 253–271 of 790 were at risk, so S1's cell failed with
probability 0.95 only at a false-SUPPORTED rate of 14.3–15.3% among them (`v1_null_cells.log`).
Requirement 1.

- Every criterion of false SUPPORTED is a cell of the cards at risk, per condition and level, sized
  so that it fails with probability ≥ 0.95 at a rate the protocol names in advance (for instance
  twice the nominal) and passes a sound validator with probability ≥ 0.95, at the design effect
  expected by donor and by gene pair (section 4).
- The blind and useless cards form their own cells (a false SUPPORTED there needs GATE 4 or 5 to
  pass the metric first), so that dilution cannot make a criterion easy.
- The pilot that measures the mix of truths before the key resizes the cells (v1's pilot showed
  the mix and did not: `oc_pilot.log`, a valid metric on 2,453 of 8,910 pair-metric cards).
- Every analysis a statement rests on is pre-registered, the test of each high-level null cell
  against its nominal with its correction over the cells included; everything else is marked post
  hoc (requirement 2). A prediction gets a cell that can confirm it, and names every criterion it
  binds (requirement 3).

## 3. A truth about validity independent of GATE 4's mechanism

v1's oracle called a metric valid when it responded to the injection GATE 4 itself uses, so S5
measured GATE 4's sampling error and could not find a probe built the wrong way; R1's module
probe was one (Review 1 B7; D1). The real effects, the truth about validity and GATE 4's probe were
all one injection, so that the engine's refusal on a blind metric was right holds only for metrics
blind to that injection. Requirement 6.

Candidates, for the owner's choice (they can be combined):
- **Two methods on the same animals:** Tabula Muris Senis mice measured both by FACS (Smart-seq2)
  and by droplet (10x), with the deeper method as the truth about whether a difference or a
  coupling exists that the shallower data's metric must see. First check, by code independent of
  the engine, that the mice and tissues are really shared between the two methods (the mouse ids,
  their ages and sexes, the tissues each method covers). The 3-month FACS data are by every sign in the metadata the
  original Tabula Muris's (requirement 7's check of the source), so the comparison stays within one
  age, or within the later time points.
- **A generative model other than GATE 4's:** the coupling planted by a shared latent program in a
  simulator fitted to the real background (marginals and a copula), with the truth the metric's
  population response to the program's strength; GATE 4 keeps its thinning injection.
- **External biology:** pairs and modules whose co-regulation is known outside the data (members of
  one protein complex, or co-expressed across many independent datasets), against pairs known not
  to be, with the truth from that knowledge, not from any injection.
- **Physical standards:** mixtures with known proportions (cell lines, spike-ins) for composition
  and content claims.

The truth's mechanism and its thresholds are fixed before the pilot, and the oracle that applies
them never imports the engine, as in v1.

**A real control for coupling: searched for, none verified (2026-10-18).** The owner removed Marsh et
al. 2022 as a control for coupling (its primary contrast is the appearance of Fos's expression, 23.0%
of the cells against 0.2%) and asked for a control where both genes are detected in at least 10% of
the cells in both groups and the difference of coupling follows from a known cause, not Kang et al.
2018 or Levy et al. 2020. A web search from this session (the journals' and GEO's pages are not
reachable from it, so nothing below has been read in full or checked on data) found no benchmark of
that form; the summary of a 2024 methods paper ("Leveraging gene correlations in single cell
transcriptomic data", BMC Bioinformatics, doi:10.1186/s12859-024-05926-z) names the absence of ground
truth as what hinders such assessments. Two candidates, each to be checked by code
independent of the engine before it could be registered (requirement 7):
- **Mosaic loss of chromosome Y in blood (LOY).** In a man with LOY in a share f of his leukocytes the
  Y genes are absent from the same cells, so the correlation of two Y genes across the donor's cells
  rises with f(1 − f); in men without LOY it is only the cells' common variation. Both genes stay
  detected in at least 10% of the cells while f is moderate, the donors are the replicates, and the
  cause is a physical one. LOY is scored from single-cell expression in published work (Cell Mol Life
  Sci 2021, doi:10.1007/s00018-021-03822-w; a 2025 medRxiv preprint on OneK1K's 982 donors,
  doi:10.1101/2025.06.01.25328624). The circularity: where
  the donors' LOY is called from the same cells' Y genes, the groups are defined by the data the
  metric reads; the calling genes must exclude the pair, and a DNA-based LOY call (a genotyping array)
  would remove it. It is a coupling of co-absence, created by a mixture of cells within the donor:
  a metric validated on it is validated for that mechanism.
- **Promoter–reporter lines (scCoBench, bioRxiv 2025, doi:10.1101/2025.05.26.656221; Arabidopsis
  roots).** A reporter driven by a gene's promoter is coupled with the native gene by construction;
  between two lines with different promoters the reporter's coupling with one native gene differs by
  a known cause. Whether the four datasets hold replicates of a line, and both genes at 10% in both,
  is not known; without replicates it measures a metric's sensitivity within a sample, not the engine's
  verdict between groups.

**If neither passes its check, the sensitivity to coupling is measured in v2 by planted causes**, each
labelled with its mechanism: (i) a latent program in a generative model fitted to the real background
(the candidate above), with the truth the metric's population response to the program's strength;
(ii) a mixture of real cells of one type in two states, mixed in donor-specific proportions, with the
truth the coupling the mixture creates; (iii) the natural or constructed causes above, where checked.
**The caveat of circularity:** each truth defines coupling by its mechanism — the probe's, the model's,
the mixture's — so a metric's sensitivity is measured to that mechanism and not to regulatory
coordination at large; where the truth is GATE 4's own injection the measure is only the probe's
sampling error (v1's S5), so at least one mechanism must be independent of GATE 4's.

## 4. Pairs beyond the most correlated, and enough of them

v1 drew the 8 most correlated pairs of each expression level (Review 1 B12). The design effect of
cards sharing a gene pair, which v1 did not compute, is large where the outcome depends on the pair:
3.91 for S1's N8 cell (6 of its 11 false SUPPORTED on one pair), 29.8 for S5's null cell, 5.93 for
S3's null stratum (`v1_pair_deff.log`); and v1's sensitivity rests on 3 pairs (118 of the 121 valid
real-effect cards with |Δ\*| ≥ SESOI SUPPORTED, all on 3 pairs). Requirements 4 and 11.

- Pairs drawn at random within each expression level, next to the most correlated and to weakly
  correlated ones, with the stratum of each pair recorded.
- Enough pairs that every cell's effective sample size, after the design effect by pair, reaches
  its size in section 2; the design effect is computed by donor, by pair and by both
  (`v1_pair_deff.py`'s estimator) for every cell.
- **For the primary outcome, at least 200 gene pairs at risk** (with a valid or ambiguous metric),
  against v1's 18, so that the 95% interval of π\* over pairs is within ±0.1 (the owner,
  2026-10-12). Resampling v1's pairs K at a time, the half-width falls to 0.1 at K = 200 against
  (b), 150 under the rule of v0.4 change 6, and 125 against (c); the square-root rule from v1's own
  interval gives 116 to 145 (`v1_correction_rule.log`, section 5). Against (b) v1's interval
  reaches 0, where π\* is clipped, so there the square-root rule understates the need (145 against
  200).

## 5. The engine as a user runs it

In v1 the SESOI, GATE 4's dose and delta_min came from the pilot, that is, from the oracle
(requirement 9). v2 runs v0.4 with its own rules (v0.4 changes 3 and 4) and with the settings a
user would choose, never with the oracle's. Where the rule leaves a choice to the user, v2 includes
cards with defaults and cards with a typical user's choice, as separate conditions.

## 6. Metrics of the types the real claims use

v1 tested one metric, the log-normalized Pearson correlation of a gene pair, and three useless
metrics. v2 is narrowed to the types of metric in the list of real claims
(`validation/plans/real_claims.md`), each with a probe v0.4 has (requirement 12, narrowed on
2026-10-11):
- coordination: a pair correlation and, if the real claims keep the global coordination level, a
  measure of a set of genes, with the set coupling probe of v0.4 change 9;
- mean expression of a marker present in only some cells or absent in one group, through the module
  probe (v1 ran it only on the useless N6c; v0.4 changes 1 and 2).

Metrics of noise and entropy, which would need probes of their own that v0.4 does not have, move to
the next version.

## 7. Real positive controls checked by independent code before registration

R1, the one real positive control expected to succeed, failed in two of its three claims: B2's
donor ids were samples, 6 of 50 pooling both sexes, the sexes never shared a development stage, and
the protocol counted ids per sex over all stages (D1; `b2_sex_structure.log`; Review 1 B8, B14).
Requirements 13–15.

- Each control's data structure is checked by code independent of the engine before the tag: each
  id one animal of one sex (by Xist and the Y genes), both sexes in one stratum, the strata the
  engine will form; the allowed set computed by that code, as the oracle computes the panel's.
- Each control runs end to end on its real data before the tag, with the run's loader, gene cut and
  library size (D1: a sparse matrix the engine could not read; a cut to 2,000 genes that changed
  CP10k's denominator).
- Each probe is checked by arithmetic against the data's composition before registration (for the
  module probe: the share of cells with the marker times its response there, against delta_min).
- **The source of the data is checked for every comparison** (the owner's second decisions of
  2026-10-18): in TMS FACS the 3-month mice come, by every sign in the metadata, from the original
  Tabula Muris and the older ones from TMS (`FACS.selection` empty at 3 months and filled after, the
  cells' names in two formats, GSE132042 beginning with GSE109774's samples, and the young male units
  with 1,298 more genes and 189,000 fewer counts per cell, exact p 0.016 and 0.008;
  `real_claims.md`, "The source in Tabula Muris Senis FACS"). Age is taken as confounded with the
  source in any comparison of 3-month mice with older ones: such a comparison is not a positive
  control, a null or a real claim of age in v2 unless the source is shown shared. The authors' text
  confirms it (the fourth run, 2026-10-19): the 3-month time point is the previously published Tabula
  Muris, and it was the reference of the later time points' annotation, so a comparison within a cell
  type compares labels given by two procedures.

## 8. Development data, conditions and reporting

- Real development backgrounds at the confirmatory run's settings, and every comparison with v1 on
  the same denominators (requirement 16).
- No outcome fixed by a test's discreteness: each design's smallest attainable p checked before
  registration (requirement 17; p07 and N4's 3 against 3).
- N2 tests the correction of depth: a level where the oracle detects its raw difference with power
  ≥ 0.9; the real effects need a key dose with Δ\* clearly above the SESOI; otherwise the run stops
  before the key (requirement 7).
- No condition on one background only; chemistry recorded and declared as a stratum where it is
  mixed (requirement 8). GATE 6 and the content estimand in the panel (requirement 10).
- Operating characteristics computed for a validator sound within the engine's correction family,
  so that what the criteria predict for an out-of-family artifact is explicit (requirement 5).
- The establishable cards reported per stratum, and the cards at risk next to every rate
  (requirement 19).
- **The errors a correction makes are counted** (requirement 20): on every artifact condition, a cell
  for the cards whose corrected effect is detected while the raw one is not, and for those where the
  correction removed a raw false effect. On v1's null cards at risk the gates removed 41 of the
  donor-level test's 74 false effects and added 21 (`v1_analyzers.log`, section 5): 18 with a raw
  difference the replicate-level test did not detect, which the correction enlarged (15, all 7
  false SUPPORTED of N3 at f = 0.4 on the high level among them), shrank (1) or reversed (2); and 3
  of N4 at 3 against 3, where the engine's t interval detected the raw difference and (b)'s exact
  test cannot go below p = 2/20 = 0.1 (`v1_correction_rule.log`, section 1). v0.4's rule (change 6)
  would have turned the 18 into INCONCLUSIVE, and 3 of v1's 275 true SUPPORTED.
- **Absences the correction makes** (the owner, 2026-10-12; v0.4 change 6): real effects with an
  artifact outside the correction's family (E + N3, E + N8), with a cell for the NO DETECTABLE
  EFFECT that only the corrected data show. v0.4 marks them in the report; if v2 finds errors among
  them, v0.5 extends the rule to absences. On v1 the correction made 85 of the 1,676 NO DETECTABLE
  EFFECT, none wrong by v1's rules (`v1_correction_rule.log`, section 3).
- **The verdict does not depend on the claim's wording** (the owner, 2026-10-12; v0.4 changes 1 and
  2): every marker and pair claim also in its mirrored wording, the order of the groups and the
  direction swapped ("Xist is higher in females", "Xist is lower in males"). The verdict and its
  cause must coincide on every such pair of cards, without exception.
- N8 with the power to decide whether v0.5 needs a correction for per-cell capture (v0.4 change 5).
