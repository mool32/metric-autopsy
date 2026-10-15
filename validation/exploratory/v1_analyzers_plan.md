# What the gates add over a donor-level test: the plan (written before the run)

Exploratory. Validation v1's datasets are development data now (v1.md section 6); nothing here
changes v1's result or the engine. This file is committed before the analyzers run on the
panel, as the owner decided on 2026-10-11 ("size chosen and recorded before the run"). The code
(`v1_analyzers.py`) is written after this commit and tested on the first 20 datasets of the
key's order; any change to what is below is listed in its log as a deviation, with its reason.

## The question

What do v0.3's gates add over a careful test at the level of the donors?

## The sample: every card of v1

All 10,540 claim cards on all 9,750 datasets, no subsample. The datasets are rebuilt from the
key with the tag's frozen `panel.build` and the run's compact backgrounds (`blind.load_prepared`);
every rebuilt dataset must have the manifest's `data_sha256`. A dataset that does not is left out
of every analyzer and counted in the log. The claim cards (metric, gene pair, groups, direction,
SESOI) are the run's.

## The four analyzers, on the same data and cards

- **(a) A test over cells**, as most papers do. The claim's metric on all cells of each group, and
  a two-sided permutation test that relabels cells at random, keeping the group sizes (999
  relabelings; p = (1 + #{|Δ\*\*| ≥ |Δ|}) / 1000, Δ = metric(B) − metric(A)). Donors are ignored.
- **(b) A test over donors, without gates**: the best practice (the analogue of a pseudobulk
  test for a metric of a population of cells). The claim's metric on each donor's cells, one value
  per donor, and a two-sided permutation test over donors of the difference of the groups' means:
  exact over all relabelings where there are at most 20,000 (8 against 8: 12,870; N4's 3 against
  3: 20, so its smallest p is 0.1), else 9,999 random ones (N7's 12 against 12). Where every
  donor has cells in both groups (N5) the value is the donor's difference B − A, and the test flips
  the signs of the donors (all 256). No correction of depth, no check of the metric. Written with
  numpy only, without the engine.
- **(c) The engine v0.3 with GATE 0, 2, 4 and 5 off.** The engine's own functions, unchanged:
  GATE 1 (a diagnostic), its replicate design and graded rule (`effect.assess_design`,
  `effect._rep_values`, `effect._infer`: a permutation over replicates with at least 4 per group,
  a t interval with 3, no effect verdict with 2 or fewer or without a declared replicate unit),
  on the uncorrected data; its power and TOST against the SESOI on the observed scale; and its
  verdict rule (`report.decide_cause`), with the metric's validity taken as passed (the gates
  that check it are off) and no replication field. Its verdicts: SUPPORTED, NOT SUPPORTED
  (opposite direction), NO DETECTABLE EFFECT, INCONCLUSIVE, UNIDENTIFIABLE. The engine's seed is the
  card's (`run_panel._seed`).
- **(d) The full engine v0.3**: the published reports (`results/panel-v1`, results commit
  3cb3c211), read as v1's scoring reads them (`score.card_rows`).

(a) and (b) read the dataset, not the card's replicate column, so on N4 they give one answer for
both of a dataset's cards; (c) and (d) follow the card (N4's first card has no replicate unit).

## The outcomes, for each analyzer

- **"Effect" (in the claimed direction)**: (a), (b) p < 0.05 with Δ in the claimed direction;
  (c), (d) SUPPORTED. On null data a false "effect"; nominal 2.5%.
- **"Effect against the claim"**: (a), (b) p < 0.05 the other way; (c), (d) NOT SUPPORTED, opposite
  direction. Nominal 2.5% on null data; with "effect", 5%.
- What only (c) and (d) can say: NO DETECTABLE EFFECT, INCONCLUSIVE, UNIDENTIFIABLE; only (d):
  metric invalid, explained by depth, DEGENERATE METRIC, GATE 0's refusal.

Reported, each rate as k/n with its 95% Clopper-Pearson interval:

1. **False "effect" on null data** (N1–N5, N7, N8), per condition and variant, per expression
   level, and on the cards at risk (a valid or ambiguous metric by the oracle's truth) against the
   blind ones; also "effect" or "effect against the claim" together.
2. **Real effects** (E1 at its four doses, E2, E3): the share "effect", per condition, level and
   truth; the main cell is the cards with a valid metric and |Δ\*| ≥ SESOI.
3. **Blind and useless metrics**: every analyzer's outcomes on the cards with a blind metric (null
   data and real effects apart) and on N6a (random numbers), N6b (a constant) and N6c (a score of
   random genes).
4. **A check of (c)**: where (d)'s report holds the engine's raw replicate-level test (`raw_detected`,
   `p_raw`, with an exact permutation), (c)'s must equal it.

## Two more numbers for v1, from the published reports

- **v1's sensitivity**: the share SUPPORTED among the cards with a valid metric and |Δ\*| ≥ SESOI
  (and among all cards with a valid metric), with its interval.
- **The design effect by gene pair**: for every cell of v1's criteria (S1–S6, S3's strata and
  conditions) and for the main rates above, the design effect of cards that share a gene pair,
  with the estimator v1 used for shared donors (`score.overlap_interval`, a pair as the shared
  unit), next to v1's design effect by donor; and both together, Cov = β_donor × (shared donors) +
  β_pair × (same pair), fitted by least squares on the residuals' cross-products.

Nothing here is a criterion: there is no pass or fail, and no number is compared with a limit
fixed in advance.
