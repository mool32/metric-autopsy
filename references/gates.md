# The Gates — full definitions

*Progressive-disclosure reference for the `metric-autopsy` skill. Loaded on demand;
the thin protocol lives in `../SKILL.md`. This file is also the source text for the
methods section of the preprint — one text, two purposes.*

Before you claim that a computed metric reveals a biological signal, the metric has to
survive these gates **in order**. A failure at any gate means **STOP and fix before
proceeding** — a downstream gate cannot rescue an upstream failure. Each gate exists
because a real analysis died on it.

The gates are **metric-agnostic**. You supply a callable `metric(data) -> float` and the
names of your factorial `obs` columns; the gates treat the metric as a black box and probe
the *data* and the metric's *response to controlled perturbations of the data*. That is the
whole design: the harness knows single-cell QC, not your metric.

---

## The three errors these gates encode

Every gate traces to one of three failures we actually shipped weeks of work into.

**Error 1 — Mathematical dependency (ACP).** Claim: intracellular entropy and intercellular
entropy are anticorrelated during aging (ρ = −0.54). Reality: when cells become internally
uniform they also become mutually similar — a structural tendency of the two metrics sharing
one matrix, not biology. The signal vanished on 10x (ρ = −0.09) and *reversed* at low depth
(ρ = +0.38 at 500 genes). → **GATE 0**.

**Error 2 — Physical constant as biology (β-heart).** Claim: the ECG spectral exponent β
encodes cardiac health. Reality: β tracks conduction-system geometry — a biophysical constant
of the medium. The 12-lead β-vector reconstructs anatomy (AUC 0.98 LBBB vs RBBB) *because it
is the anatomy*. → **GATE 4**.

**Error 3 — Technical confound (MI coupling).** Claim: SMAD→ECM mutual information declines
with age while housekeeping coupling is preserved. Reality: male old cells detect 2.4× fewer
genes (3670 → 1540, p = 9.5e-56); female cells are stable. MI with a zero-bin is driven by
detection rate. The "pathway-specific coupling loss" was QC in the sex×age interaction —
invisible in the young-vs-old marginal. → **GATES 1, 2, 5, 6**.

---

## GATE 0 — Mathematical independence  · *auto*

**Question.** Can the metric move between conditions purely from a change in a summary
statistic (mean, variance, sparsity, dimensionality, n) with **no change** in the biological
quantity of interest?

**How the engine tests it.** It perturbs *your own data* — it does not synthesize a separate
null. First it builds a **bootstrap baseline**: resample cells with replacement and re-evaluate
the metric ~60× to characterise the metric's own sampling spread (its estimator noise with the
biology intact). Then it applies each technical nuisance ~20× to the same data and re-evaluates:
`extra_dropout` (zero out a fraction of detected entries), `depth_downsample` (binomial thinning
of counts; per-cell scaling for non-count input), `library_scale` (per-cell library factors),
and — for whole-matrix metrics — optional `variance_inflation` and `gene_subsample` (the genes
the metric is bound to are protected from subsampling). A nuisance is flagged as confounding
only if the perturbed mean is **both** meaningfully large relative to the metric's scale
(relative change > `tol`, default 0.25, with a floored denominator so a near-zero baseline can't
explode) **and** statistically separated from the baseline (z > `z_thresh`, default 4). The
two-part rule is the point: it separates a real shift in the metric's *expectation* (a confound)
from mere estimator *noise*, so a noisy-but-unbiased metric is not wrongly failed.

**Classic confounded forms.** MI with a zero-bin (sparsity); eigenvalue ratios (variance);
correlations across differing zero-inflation; any ratio whose numerator/denominator scale
differently with library size.

**Pass.** The metric's expectation is stable under every nuisance perturbation (no shift beyond
estimator noise) — **or** you ship an explicit correction for every statistic it responds to.

---

## GATE 1 — QC parity across ALL factorial combinations  · *auto · crown jewel*

**Question.** Do the groups you compare have equivalent technical quality — in *every*
stratification you will use, not just the primary axis?

**How the engine tests it.** For each combination of `age × sex × batch × tissue × cell_type`
(whichever `obs` columns you name), the gate flags any stratum where the two compared groups'
median `n_genes_by_counts` differs by more than **1.5×** (`thresh`) **or** whose `n_genes`
distributions barely overlap — where overlap is the shared fraction of the *narrower* group's
10th–90th-percentile range, flagged below `overlap_min` (default 0.2). `per_cell_qc` also
reports `total_counts` and, where identifiable from `var_names`, mito/ribo fractions for
inspection, but only n_genes ratio and overlap drive the pass/fail; there is no per-gene
detection-rate criterion.

**Why it is the crown jewel.** This is the check that kills the most work for the least effort.
We compared young vs old pooling sexes; the confound lived in sex×age — only male-old was
degraded. A 45-second QC check would have caught what three weeks of MI analysis did not.

**Pass.** All compared factorial pairs < 1.5× QC difference — **or** you restrict every metric
computation to QC-matched subsets (GATE 2).

---

## GATE 2 — n_genes matching preserves the signal  · *auto · crown jewel*

**Question.** When you equalize cell quality between groups, does the signal survive?

**How the engine tests it.** Find the overlapping n_genes range (10th–90th pct) between groups;
keep only cells inside it; confirm the matched subsets are balanced (median n_genes within
~1.1×); recompute the metric on the matched subset via your callable; compare to the unmatched
result. Guards both ends: an effect below a permutation-null floor has nothing to preserve, and
one that *amplifies* far beyond the unmatched value flags selection on the quality axis.

**Read-out.**
- Survives matching (> 50% of original effect) → likely biological.
- Disappears with matching → likely a QC artifact.
- Distributions do not overlap → groups are **incomparable. STOP.** (In MI-coupling, male
  young vs old had *3 cells* in the overlap.)

**Pass.** Signal survives n_genes matching with > 50% of the original effect size.

---

## GATE 3 — Visible in the raw data  · *auto (exports the plot data) + judgment*

**Question.** Can you see the effect in raw values, without the metric detecting it for you?

**How the engine tests it.** For a pairwise metric it exports the scatter of the raw inputs
(gene A vs gene B) split into the two compared groups, and auto-hints when the between-group
shape change is dropout-driven — when the zero-fraction gap between the groups exceeds 0.15.
Optionally it writes a two-panel A/B scatter PNG. You look. (Stratifying by sex and coloring by
n_genes yourself is the recommended manual follow-up; the control-pair scatters live in GATE 5.)

**What to see.** Elongation along the diagonal → real coupling; circular cloud → none; a shape
change between conditions → a coupling change — *unless* one condition simply has more points
piled on the zero axes, in which case the "shape change" is dropout. If the effect follows the
n_genes color gradient, it is a QC artifact. We computed MI for months before plotting; when we
finally looked, the coupling was not there.

**Pass.** The effect is visible in the raw scatter, split by group, and the between-group shape
change is not merely dropout — without needing the metric to surface it.

---

## GATE 4 — The metric measures what you think  · *judgment (the skill asks)*

**Question.** Does a change in the metric correspond **uniquely** to a change in the biology?

Enumerate every scenario that could move the metric and decide which your result is consistent
with:

| Scenario | Metric moves? | Biology moves? | Verdict |
|---|---|---|---|
| Real effect | yes | yes | true positive |
| Detection-rate shift | yes | no | false positive (QC) |
| Mean-expression shift | maybe | maybe | ambiguous |
| Variance change | maybe | no | math artifact |
| Sample-size change | maybe | no | statistical artifact |
| Batch effect | yes | no | technical confound |
| Composition shift | yes | maybe | Simpson's paradox |

If more than one row explains your result, you owe a test that separates them. This is also
where β-heart dies: a metric that perfectly tracks a physical/geometric property is *measuring
that property*, not a separate biological process.

**Pass.** You can rule out every non-biological row — or you state the ambiguity in the paper.

---

## GATE 5 — Controls behave in ALL strata  · *auto*

**Question.** Does a known-positive control show the effect and a known-negative control not —
in **every** factorial combination, not just pooled?

**How the engine tests it.** Runs your metric on a positive-control pair (housekeeping,
Actb↔Gapdh — should always couple) and a negative-control pair (random — should not), **per
stratum**. The null band is anchored to the *positive* control's magnitude (not the negative
control's own value, which is the quantity under test): the negative must stay below ~20% of the
positive signal, and the positive must clear that band. A control that "passes" only after
averaging over a confound is worthless: our HK control looked stable pooled, but declined in
males once stratified — confounded the same way as the test.

**Pass.** Positive control positive, negative control null, in all strata.

---

## GATE 6 — Cross-platform / cross-species replication  · *auto if a 2nd dataset is supplied*

**Question.** Does the effect replicate on independent data with different technical
characteristics — after QC matching?

**How the engine tests it.** Re-runs GATES 1–2 and the metric on a second `data` object you pass
(different platform, lab, or species). Platform-specific artifacts (SmartSeq2 dropout, 10x UMI
saturation) mimic biology but do not replicate: mouse TMS (SmartSeq2) showed coupling loss;
human skin (10x) showed it unmatched, but n_genes-matched it vanished even in the young.

**Pass.** Effect replicates on ≥ 1 independent dataset after QC matching. *(In v1 this gate is a
documented procedure demonstrated in `examples/`; supply a second dataset to run it.)*

---

## GATE 7 — Effect size is biologically meaningful  · *judgment (the skill asks)*

**Question.** Is the effect large enough to matter for cell function?

Calibrate against the metric value expected for a known, validated interaction; against test–
retest variability; against whether the change would move downstream expression enough to be
phenotypically relevant. **Red flag:** if you need 10,000+ cells to reach significance, it may be
real but biologically irrelevant — biology runs at single-cell scale.

**Pass.** Effect exceeds test–retest variability and is in range for a real regulatory interaction.

---

## Ordering and the mandatory sequence

```
Step 0  State the hypothesis in one sentence.
Step 1  QC across ALL factorial combinations (age×sex×batch×tissue×cell_type). Flag >1.5× now.
Step 2  Plot QC distributions — do they overlap? If not, those groups are incomparable. Stop.
Step 3  Plan the n_genes-matched analysis from the start.
Step 4  Define positive control, negative control, and what would disprove the claim.
Step 5  NOW compute the metric.
Step 6  Scatter the raw data, stratified by sex and QC quartile.
Step 7  Recompute on the QC-matched subset.
Step 8  Vary the main hyperparameter 2×.
Step 9  Replicate on a second platform/species if available.
```

The order is the point. We spent three weeks and 500+ lines on MI before checking QC; the QC
check took 45 seconds and invalidated most of the work. The cost of the gates is about an hour.
The cost of skipping them is weeks and a false conclusion.

---

## Retroactive audit (what the gates would have said)

| Metric | First failing gate | Verdict |
|---|---|---|
| ACP eigenvalue / entropy | GATE 0 (variance-structure function) | not biology |
| Cardiac β | GATE 4 (physical constant) | real, but physics not biology |
| MI SMAD→ECM (zero-bin) | GATE 0 (first applicable) | FAIL at the first applicable gate (GATE 0) and additionally 1, 2, 3, 6, with GATE 5 a partial fail; judgment gates 4/7 not scored — **not measuring biology** |

*Re-read this before starting any new metric analysis.*
