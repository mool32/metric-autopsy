# Validation v2: requirements (draft for the owner)

A draft, not a pre-registration. It makes operational the requirements collected in `PROJECT.md`
("Requirements for v2", items 1–19, each with its source) and the owner's list of 2026-10-11. The
engine v2 tests is v0.4 (`validation/plans/v0.4_spec.md`); each of its changes names its test here.
Numbers come from the logs named, each from a script in this repository. v1's datasets are
development data now (v1.md section 6), so v2's panel is drawn afresh, from backgrounds v1 did not
use.

## 1. The criteria count the cards at risk

v1's criteria of false SUPPORTED were counted over all of a condition's cards, and most of them
had a blind metric: in each key null 253–271 of 790 were at risk (a valid or ambiguous metric), so
S1's cell failed with probability 0.95 only at a false-SUPPORTED rate of 14.3–15.3% among them
(`v1_null_cells.log`). Requirement 1.

- Every criterion of false SUPPORTED is a cell of the cards at risk, per condition and level, sized
  so that it fails with probability ≥ 0.95 at a rate the protocol names in advance (for instance
  twice the nominal) and passes a sound validator with probability ≥ 0.95, at the design effect
  expected by donor and by gene pair (section 3).
- The blind and useless cards form their own cells (a false SUPPORTED there needs GATE 4 or 5 to
  pass the metric first), so that dilution cannot make a criterion easy.
- The pilot that measures the mix of truths before the key resizes the cells (v1's pilot showed
  the mix and did not: `oc_pilot.log`, a valid metric on 2,453 of 8,910 pair-metric cards).
- Every analysis a statement rests on is pre-registered, the test of each high-level null cell
  against its nominal with its correction over the cells included; everything else is marked post
  hoc (requirement 2). A prediction gets a cell that can confirm it, and names every criterion it
  binds (requirement 3).

## 2. A truth about validity independent of GATE 4's mechanism

v1's oracle called a metric valid when it responded to the injection GATE 4 itself uses, so S5
measured GATE 4's sampling error and could not find a probe built the wrong way; R1's module
probe was one (Review 1 B7; D1). Requirement 6.

Candidates, for the owner's choice (they can be combined):
- **A generative model other than GATE 4's**: the coupling planted by a shared latent program in a
  simulator fitted to the real background (marginals and a copula), with the truth the metric's
  population response to the program's strength; GATE 4 keeps its thinning injection.
- **External biology**: pairs and modules whose co-regulation is known outside the data (members of
  one protein complex, or co-expressed across many independent datasets), against pairs known not
  to be, with the truth from that knowledge, not from any injection.
- **Physical standards**: mixtures with known proportions (cell lines, spike-ins) for composition
  and content claims.

The truth's mechanism and its thresholds are fixed before the pilot, and the oracle that applies
them never imports the engine, as in v1.

## 3. Pairs beyond the most correlated, and enough of them

v1 drew the 8 most correlated pairs of each expression level (Review 1 B12). The design effect of
cards sharing a gene pair, which v1 did not compute, is large where the outcome depends on the pair:
3.91 for S1's N8 cell (6 of its 11 false SUPPORTED on one pair), 29.8 for S5's null cell, 5.93 for
S3's null stratum (`v1_pair_deff.log`); and v1's sensitivity rests on 3 pairs (118 of the 121 valid
real-effect cards with |Δ\*| ≥ SESOI SUPPORTED, all on 3 pairs). Requirements 4 and 11.

- Pairs drawn at random within each expression level, next to the most correlated and to weakly
  correlated ones, with the stratum of each pair recorded.
- Enough pairs that every cell's effective sample size, after the design effect by pair, reaches
  its size in section 1; the design effect is computed by donor, by pair and by both
  (`v1_pair_deff.py`'s estimator) for every cell.

## 4. The engine as a user runs it

In v1 the SESOI, GATE 4's dose and delta_min came from the pilot, that is, from the oracle
(requirement 9). v2 runs v0.4 with its own rules (v0.4 changes 3 and 4) and with the settings a
user would choose, never with the oracle's. Where the rule leaves a choice to the user, v2 includes
cards with defaults and cards with a typical user's choice, as separate conditions.

## 5. Metrics of the types the target claims use

v1 tested one metric, the log-normalized Pearson correlation of a gene pair, and three useless
metrics. Requirement 12. v2 includes, each with its own injected signal for GATE 4 and its own
truth (section 2):
- coordination and correlation: a pair correlation, mutual information, the coherence of a module;
- noise: the coefficient of variation or the Fano factor of a gene, and a residual variance after
  the mean–variance trend;
- entropy: the entropy of a cell's expression profile, and of the population's;
- mean expression of a marker present in only some cells, through the module probe (v1 ran it only
  on the useless N6c; v0.4 change 1).

## 6. Real positive controls checked by independent code before registration

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

## 7. The analyzers of the decisive experiment in the panel

What the gates add has to be measured, not argued (requirement 18). Every card of v2 is also run
through the three analyzers of the decisive experiment on v1's data (`v1_analyzers_plan.md`): (a) a
test over cells, (b) a test over donors without gates, (c) the engine with GATE 0, 2, 4 and 5 off.
v2 pre-registers, as secondary outcomes with their intervals, the differences between (b) and the
engine in false "effect" on the cards at risk and on the blind and useless ones, in sensitivity on
the valid cards with |Δ\*| ≥ SESOI, and in how often each says "no effect" where the metric is blind
to a real effect.

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
