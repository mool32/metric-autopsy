# The Gates — full definitions

*Progressive-disclosure reference for the `metric-autopsy` skill. Loaded on demand;
the thin protocol lives in `../SKILL.md`. This file describes the engine in this checkout
(v0.3.0.dev0, under validation). The preprint's methods (`paper/manuscript.md`) still describe
v0.1.1 and will be rewritten after the confirmatory freeze.*

The gates are **metric-agnostic**. You supply a callable `metric(data) -> float` and the
names of your factorial `obs` columns; the gates treat the metric as a black box and probe
the *data* and the metric's *response to controlled perturbations of the data*. The harness
knows single-cell QC, not your metric.

> **Status: v0.x under validation.** The probes in `validation/probes/` are a development
> set: the fixes described here were developed against them, so passing them proves nothing
> about the validator's operating characteristics. Do not present a verdict as validated.

---

## How the verdict is decided

A metric that changes between conditions is not a finding. The verdict is not "the first
gate that fails" any more. It is decided by one rule (`report.decide`, shared by the Python
API, the CLI and the MCP server) from **four independent fields**, each answering a
different question. The verdict comes with its **cause** (`report.decide_cause`; `cause` in the
JSON report): one code per line of the rule below, e.g. `metric_invalid_gate4` (shown blind to the
injected signal), `metric_invalid_gate0` (a nuisance bias), `explained_by_depth`,
`opposite_direction`.

| Field | Question | Statuses |
|---|---|---|
| `metric_validity` | Does the metric respond to its construct and resist nuisance? | PASS · FAIL · UNTESTED · DEGENERATE (flags: LEVEL_SHIFT, DEPTH_BIAS_CORRECTED, BIAS_BELOW_TOLERANCE, BIAS_UNSIZED) |
| `design_adequacy` | Is the comparison identifiable, and can it detect the SESOI? | ADEQUATE · CORRECTED · INSUFFICIENT_REPLICATION · UNIDENTIFIABLE (flags: PARAMETRIC_ONLY, UNDERPOWERED, ATTENUATION, POWER_NOT_ASSESSED) |
| `effect` | Is there a difference at equal depth, across biological replicates? | DETECTED · NO_DETECTABLE_EFFECT · INCONCLUSIVE · NOT_ESTIMABLE · NOT_RUN |
| `replication` | Does it hold on independent data? | REPLICATED · NOT_REPLICATED · INCONCLUSIVE · NOT_RUN |

`decide` reads them in this order; the first matching line is the verdict:

1. metric DEGENERATE → **DEGENERATE METRIC**
2. metric FAIL (a nuisance bias beyond its SESOI tolerance, a failed control, a response to an injected signal shown below `delta_min`) → **NOT SUPPORTED — metric invalid**
3. design UNIDENTIFIABLE (no estimand, content estimand without spike-ins, partially crossed replicates, groups confounded with a stratifier) → **UNIDENTIFIABLE**
4. effect not evaluated → **INCONCLUSIVE**; the raw difference is explained by depth/capture → **NOT SUPPORTED**
5. design INSUFFICIENT_REPLICATION → **INCONCLUSIVE — insufficient replication**
6. effect NO_DETECTABLE_EFFECT → **NO DETECTABLE EFFECT**, only if the metric is PASS (absence of an effect cannot be claimed with a metric whose response is untested)
7. effect not DETECTED → **INCONCLUSIVE**
8. metric UNTESTED → **INCONCLUSIVE** (an effect seen by a metric whose response to signal was never demonstrated)
9. a nuisance bias that could not be sized (no SESOI) → **INCONCLUSIVE**
10. no pre-registered `direction` → **INCONCLUSIVE**; the effect is in the direction opposite to
    the pre-registered one → **NOT SUPPORTED**
11. replication NOT_REPLICATED → **NOT SUPPORTED**
12. judgment gates 4 and 7 pending (the default) → **INCONCLUSIVE**
13. → **SUPPORTED — replicated**, or **SUPPORTED (provisional until replicated)**

**Claims are directional.** The pre-registration states the claimed change of the metric from
`groups[0]` to `groups[1]`: `increase`, `decrease`, or `two-sided`. SUPPORTED needs the effect in
that direction (the test stays two-sided at alpha, so a null gives a false SUPPORTED at alpha/2);
a non-directional claim is allowed and the verdict marks it. An effect whose sign the correction
reverses is INCONCLUSIVE, as before.

The verdict names the assumption it rests on: *parametric only* (too few replicates for the
permutation test to reach alpha), *underpowered relative to the SESOI* and *non-directional
claim*.

**No rescue language.** A FAIL in an earlier field is never softened by a later one. With
`stop_on_first_fail` (the default) an invalid metric stops the analysis before the effect is
estimated; design adequacy, which does not depend on the metric, is still reported.

**Provenance.** Every report carries `data_sha256` (X, var_names and the obs columns used),
`prereg_sha256`, a `claim_id` (data + pre-registration + comparison), versions and the seed
(`Autopsy.to_json()`, strict JSON). Every run with a pre-registration — through the Python
API, the CLI or MCP — is appended to the run log (`metric_autopsy_runs.jsonl`, or
`$METRIC_AUTOPSY_LOG`; `log_path="off"` disables), and the report states which attempt at this
claim it is and what the earlier attempts said. Runs without a pre-registration (they cannot
reach an effect verdict) and the demo are not logged unless a path is given.

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
invisible in the young-vs-old marginal. → **GATES 1, 2, 5, 6**. *(These numbers come from the
by-hand analysis; the automated real-data run is under audit in `validation/flagship_audit/`.)*

---

## GATE 0 — Mathematical independence: bias vs attenuation  · *auto*

**Question.** Can a technical nuisance *create or inflate* a difference in the metric with
**no change** in the biology? Or does it only *shrink* the metric's signal?

**How the engine tests it.** On *your own data*:

1. a **bootstrap baseline** (cells resampled with replacement, 60×): the metric's sampling spread;
2. an **automatic null** (10×): for count input every gene, independently, takes the count of
   a random cell among the 20 next cells at least as deep, binomially thinned to the cell's own
   depth (the depth-matched draw of GATE 5's positive control). Gene-gene structure is
   destroyed; each gene's dependence on depth is kept cell by cell. *Signal* = baseline − null.
   (Until 2026-10-08 genes were permuted within depth deciles. The depth variation left inside
   each decile made a pair coupled only through cell size look coupled above the null, by more
   than 0.05 of a raw correlation of 0.7–0.9; with the new null the signal stays below 0.01.
   Non-count input keeps the decile shuffle; a perturbation that leaves the counts non-integer,
   `library_scale`, is compared with the decile null of the unperturbed data);
3. each **nuisance perturbation** (20×): `extra_dropout` (20% of detected entries zeroed),
   `depth_downsample` (binomial thinning to half depth; per-cell scaling for non-count
   input), `library_scale` (per-cell factors 0.5–2), and for whole-matrix metrics
   `gene_subsample` (80% of genes, the bound genes protected).

A shift counts only if it is beyond estimator noise (z > `z_thresh`, default 4). It is then
classified:

- **bias** — the perturbation moves the *null* (beyond noise) by more than `tol` (0.25) of the
  signal or more than the SESOI tolerance below, or reverses the signal, or inflates it by more
  than either. A nuisance that does this can create a difference with no biology. Whether it
  **blocks** depends on whether anything else handles it and whether it matters for the claim:
  - a bias in **depth** is reported, not failed, when the declared correction removes depth
    between the groups (composition: thinning to equal depth; content with spike-ins: thinning
    to equal capture) — the verdict then rests on the effect at equal depth;
  - any other bias (dropout, library scale; depth without such a correction) **FAILs** only if
    its size exceeds `bias_tolerance` × SESOI (pre-registered, default 0.5) **and** so does the
    lower bound of its 95% interval; otherwise it is reported;
  - without a SESOI a bias cannot be sized against the claim: it is reported as unsized, the
    effect is still estimated, and SUPPORTED is withheld (INCONCLUSIVE); the message says to
    declare a SESOI, against which the bias is then judged.
  At a realistic scale (probe p13, 800 cells per donor) log-normalized Pearson on a truly
  coupled pair has a resolved depth bias — the CP10k ratio correlation of its null grows as
  depth falls — that the earlier rule failed in 4 of 4 datasets although the composition
  correction removes it between groups.
- **attenuation** — the signal shrinks toward the null while the null stays put. This is
  reliability, not confounding (Spearman): uniform attenuation pulls an effect toward zero
  and cannot create one. It is **reported, not failed** — the fraction of signal lost goes to
  `design_adequacy` as a power check (below).
- **level_shift** — the metric has no material gene-gene structure (e.g. a mean score, or a
  signal above the null smaller than `tol` × the effect scale), so attenuation and bias
  cannot be told apart. A shift larger than `tol` × the effect scale (observed difference or
  SESOI) is flagged and reported; a between-group version is removed by the effect field's
  equalization.

A metric that returns the same value on every resample, null and perturbation is
**DEGENERATE**. Every rule is invariant to affine re-expressions `a·m + b` of the metric
(v0.1 divided the shift by |baseline|, so adding a constant passed a confounded metric).

**Classic confounded forms.** MI with a zero-bin (sparsity); eigenvalue ratios (variance);
correlations across differing zero-inflation; any ratio whose numerator/denominator scale
differently with library size. On QC-matched synthetic data `mi_3bin` loses ~63% of its signal
under extra dropout and ~25% under depth halving (attenuation); `norm_pearson` loses none
measurably.

**Pass.** No nuisance biases the metric. Attenuation and level shifts are reported with their
size. **WARN** — biases exist but none blocks (removed by the correction, within the tolerance,
or unsized); each is named with its size.

---

## GATE 1 — QC parity across ALL factorial combinations  · *auto · a design diagnostic*

**Question.** Do the groups you compare have equivalent technical quality — in *every*
stratification you will use, not just the primary axis?

**How the engine tests it.** For each combination of the `obs` columns you name (`age × sex ×
batch × tissue × cell_type`), it compares the two groups' median `n_genes_by_counts` and the
matchability overlap of their n_genes ranges. A stratum is flagged only when the breach is
**confident**: the (1 − α/K) bootstrap interval of the log median ratio lies beyond
±log(1.5), or the upper bootstrap bound of the overlap is below 0.2, where K is the number of
assessable strata (Bonferroni). Strata with fewer than `min_cells` (10) cells in either group
are reported, not assessed. (v0.1 flagged 82% of null datasets at 64 strata × 20 cells.)

**Why it is not a kill switch any more.** Biology moves QC: cycling cells carry ~2× more RNA,
and Xist is higher in female cells whatever the depth. A terminal GATE 1 killed both (probes
p02, p04). GATE 1 now returns **WARN** for a confident imbalance and hands it to the
estimand-dependent correction (GATE 2); **STOP** only when no stratum contains both groups.

**Why it is still the first thing to look at.** It is the check that finds the most for the
least effort. We compared young vs old pooling sexes; the confound lived in sex×age — only
male-old was degraded. A 45-second stratified QC check would have shown what three weeks of
MI analysis did not. Always stratify by every factor: pooling erases interaction confounds.

---

## GATE 2 — Estimand-dependent correction, and the effect  · *auto*

**Question.** When the technical difference between the groups is removed, in the way that is
valid for what the metric claims to measure, is there still a difference — across biological
replicates?

**The correction follows the pre-registered estimand.**

- **composition** (relative expression: correlations, module scores, MI of normalized
  counts): the deeper group of every stratum is **binomially thinned** by one common ratio,
  the ratio of the two groups' mean depths. Thinning keeps every cell and preserves expected
  composition; under a capture difference it is exact for every gene and keeps their joint
  distribution. (Matching each cell to the other group's depth *quantiles* does not: the
  keep-probability then depends on the cell's own total, which includes the genes the metric
  reads, and 18% of a pure MI artifact survived the correction in the dev data.)
- **content** (amount of RNA: total counts, genes detected, CytoTRACE-like scores): thinning
  to equal depth would delete the signal itself. The groups are thinned to equal **spike-in
  (ERCC) capture**. Without spike-ins, capture and content cannot be separated:
  **UNIDENTIFIABLE**.
- no estimand declared → **UNIDENTIFIABLE**.

*(v0.1 matched cells on n_genes. n_genes is downstream of biology — cell size, cycling, RNA
content — so matching on it is a bad control and deleted real effects; probe p03.
`gate2_ngenes_matching` is kept for backward compatibility, unused by `run_autopsy`.)*

**Replicate-level inference.** The metric is computed per biological replicate (`replicate_col`:
mouse, donor, plate), nested in groups and strata, or paired within replicates. The unit of
inference is the replicate, never the cell (v0.1 permuted cells and certified mouse-to-mouse
noise as an effect; probe p07). Graded rule, by replicates per group (or pairs):

- **>= 4**: an exact permutation over replicates (stratified; Monte Carlo with an MC error when
  enumeration is too large) is the test; the t interval estimates the effect. If the design
  cannot reach alpha (e.g. 5 pairs: min p = 0.0625), the t interval is the test, marked
  *parametric only*;
- **3**: the t interval on replicate values is the test, marked *parametric only*; the
  permutation p is reported with the note that alpha is unattainable by construction;
- **<= 2**, fewer than the pre-registered minimum, or no replicate unit:
  **INSUFFICIENT_REPLICATION** — estimates only, no effect verdict.

**Read-out.**
- **DETECTED** — the corrected effect is significant at the replicate level.
- **NO_DETECTABLE_EFFECT** — the (1 − 2α) interval lies inside ±(attenuated) SESOI (TOST).
  Needs a pre-registered SESOI; without one the result is INCONCLUSIVE.
- **explained by depth / capture** — the raw difference is detected across replicates and,
  after the correction, is not, with less than half retained → **NOT SUPPORTED**. A raw
  difference that is not itself detected across replicates cannot be "explained"; the effect
  is INCONCLUSIVE and the reason states how much of it the correction left.
- **reversed by the correction** — a detected raw difference becomes a detected difference of
  the *opposite* sign at equal depth → INCONCLUSIVE: the technical difference is larger than
  the effect, so the direction depends on how exactly the correction removed it. (In the dev
  data this turned three false SUPPORTED verdicts on a pure depth artifact into INCONCLUSIVE.)

**Power.** The SESOI is declared on the construct scale. GATE 0's attenuation under depth
halving gives a reliability-model estimate λ of how much of a construct-scale difference
survives at the analysed depth (noise ∝ 1/depth; the equalization's depth ratio is included).
The design is **UNDERPOWERED** when the minimum detectable effect exceeds λ × SESOI; the same
attenuated SESOI is used for the equivalence test.

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

## GATE 4 — The metric measures what you think  · *auto (injected signal) + judgment*

**Question.** Does a change in the metric correspond **uniquely** to a change in the biology?

**Automatic part: response to an injected signal.** A metric earns `metric_validity` PASS only
by *responding* to signal, not merely by resisting nuisance — a metric that returns random
numbers resists every nuisance. Supply `signal_test=` (`injected_signal.coupling(a, b, strength)`
or `injected_signal.module(genes)`): a known construct change is planted by binomial thinning
(seqgendiff-style, so the data stay valid counts with real technical noise) 200 times, each
**relative to a matched sham** — the same thinning without the signal (independent
keep-probabilities for a coupling; all genes thinned for a module) — and the mean response gets
its two-sided 95% t interval (decided 2026-10-08):

- **PASS** — the lower bound in the declared direction is above 0: a response is shown.
- **FAIL** — the upper bound is below `delta_min`, the smallest response that matters
  (pre-registered `delta_min`; default 0.5 × SESOI): the metric is shown blind. This holds also
  when the whole interval lies above 0: a response below `delta_min` is blindness by its definition.
- **UNTESTED** (the gate reports WARN) — neither: absence of evidence, not invalidity, as
  INCONCLUSIVE is to NO DETECTABLE EFFECT. Without a `delta_min` (no SESOI) the gate cannot FAIL.

The former rule (z ≥ 3 over 10 injections) failed a valid metric whose response is real but
weak (probe p14). Against the untouched data, instead of a sham, the thinning noise is confounded
with the signal: injecting coupling into an already strongly coupled pair *lowered* a valid
correlation and failed it. A positive control that beats its null (GATE 5) is the other way to
demonstrate response, except where GATE 4 ran on the analysed construct and is UNTESTED: a control
on another pair does not stand in for it. With neither, the metric is **UNTESTED** and no
positive verdict is reachable.

**Judgment part.** Enumerate every scenario that could move the metric and decide which your
result is consistent with:

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
that property*, not a separate biological process. The judgment is pending by default
(`judgment_pending`), which holds the verdict at INCONCLUSIVE.

---

## GATE 5 — Controls behave in ALL strata, against empirical nulls  · *auto*

**Question.** Does a known-positive control show the effect and a known-negative control not —
in **every** factorial combination, not just pooled?

**Which metric the controls test.** The controls run `pair_metric(data, gene_a, gene_b)` on the
control pairs, so they are evidence about the judged metric only if `pair_metric` bound to the
analysed gene pair *is* that metric. `run_autopsy` checks it by value on the data (both are black
boxes); on a mismatch, or without a gene pair, GATE 5 is SKIP and the controls count for nothing.
Before this check a metric that ignores gene b, given `norm_pearson`'s controls, was certified
and SUPPORTED in 3 of 3 dev datasets, and a random-number metric passed metric validity.

**How the engine tests it.** Per stratum (with at least `min_cells` cells), against empirical
nulls instead of a fixed band, Bonferroni across strata:

- the **negative control** pair is compared with 200 distinct unrelated pairs (never
  repeated) from the same expression neighbourhoods — the 20 genes closest in mean
  expression, or 5% of all genes; it must not stand out. The null centre is reported: a centre far
  from zero means the metric reports association between unrelated genes (e.g. closure from
  library-size normalization, which made the "robust" `norm_pearson` fail its own negative
  control on null data in v0.1; probe p06);
- the **positive control** pair is compared with 200 self-nulls of the *same* pair: gene b is
  replaced, cell by cell, by the count of a random neighbouring cell at least as deep, thinned
  binomially to the cell's own depth. Coupling is destroyed and the dependence on depth is kept
  exactly. It must stand out. (A shuffle *within depth bins* is not enough: the depth variation
  left inside each bin let a pair coupled only through cell size pass as a positive control in
  5 of 5 null datasets. An expression-matched pair null is unsuitable too: other truly coupled
  genes contaminate it.)

Depth bins (used by GATE 0's null and for non-count input) hold at least 10 cells: a one-cell
bin cannot be shuffled, which in an earlier draft made the positive control's null equal the
data.

p values are rank-based Monte Carlo p values, `(1 + #at least as extreme) / (1 + n)`, valid
at any n. When a control sits in the extreme tail of the first 200 draws and alpha/K is
below their resolution, the null is extended to 2K/alpha draws (at most 5000). No parametric
tail is assumed: a normal tail is wrong for skewed nulls such as MI's. With few genes the pool
of distinct unrelated pairs bounds the resolution; when alpha/K is below it, the negative
control cannot fail, and the gate's message says so.

**Was the positive control's silence informative?** In a stratum where the positive control
does not beat its null, the engine asks GATE 4's question of the control's genes (decided
2026-10-08, probe p15): it removes the pair's own coupling (gene b replaced by the depth-matched
draw), injects a coupling of known dose (`injected_signal.coupling`; pre-registered
`positive_control_dose`, default 2.0) 200 times against its sham, and gives the metric's own
response its two-sided interval at alpha/K. The metric is shown blind when the bound in the
declared direction is below `delta_min`; otherwise the silence says that the control is not
coupled here (when the metric responds to the injection), or nothing. The former rule failed a
silent control wherever a reference detector had power ≥ 0.8 for the dose, and so failed a valid
metric whose control was coupled, but weakly, in 16 of 20 p15 datasets.

**Read-out.**
- **FAIL** — the negative control stands out in some stratum (the metric reports association
  where there is none), **or** the positive control is silent in a stratum where the metric is
  shown blind to a coupling injected into its genes.
- **WARN** — the negative control is fine everywhere, and the positive control is silent
  somewhere without that proof: absence of evidence. If it fires nowhere, the metric's response
  stays undemonstrated (`metric_validity` UNTESTED, unless an injected signal shows it).
- **PASS** — positive control beats its null and negative control stays inside its null, in
  all strata.

In the dev data a metric that ignores gene b fails GATE 5 in 20 of 20 p15 datasets (with a
`delta_min`), and is UNTESTED in a 10-cell stratum where the interval is too wide. A positive control
that fires in at least one stratum (Bonferroni across strata) is evidence of response for
`metric_validity`. Passing `pos_min` / `neg_max` selects the legacy fixed band (not the
default), where a silent positive control FAILs without an injection. A control that "passes" only after averaging over
a confound is worthless: our HK control looked stable pooled, but declined in males once
stratified — confounded the same way as the test.

---

## GATE 6 — Replication on independent data  · *auto if a 2nd dataset is supplied*

**Question.** Does the effect replicate on independent data with different technical
characteristics?

**How the engine tests it.** Re-estimates the effect on the second dataset with the *same*
estimand, correction, replicate rule and strata, and re-runs the stratified QC diagnostic
there. **REPLICATED** = effect detected with the primary sign; **NOT_REPLICATED** = equivalent
to zero within the SESOI, or detected with the opposite sign (→ NOT SUPPORTED); anything else
INCONCLUSIVE. Platform-specific artifacts (SmartSeq2 dropout, 10x UMI saturation) mimic
biology but do not replicate. *(The automated GATE 6 has not yet been run on a real second
platform; the TMS → human skin result in the worked example is the by-hand analysis.)*

---

## GATE 7 — Effect size is biologically meaningful  · *judgment (the skill asks)*

**Question.** Is the effect large enough to matter for cell function?

Calibrate against the metric value expected for a known, validated interaction; against test–
retest variability; against whether the change would move downstream expression enough to be
phenotypically relevant. Declare the answer *before* the analysis as the SESOI: the engine
then uses it for equivalence and power. **Red flag:** if you need 10,000+ cells to reach
significance, it may be real but biologically irrelevant — biology runs at single-cell scale.

**Pass.** Effect exceeds test–retest variability and the SESOI, and is in range for a real
regulatory interaction.

---

## Ordering and the mandatory sequence

```
Step 0  State the hypothesis in one sentence.
Step 1  Pre-register: estimand (composition | content), replicate unit, minimum replicates,
        SESOI, alpha, controls, and the simplest non-biological explanation.
Step 2  QC across ALL factorial combinations (age×sex×batch×tissue×cell_type).
Step 3  Check the design: enough replicates per group? spike-ins for a content estimand?
Step 4  Define positive control, negative control (or an injected signal), and what would
        disprove the claim.
Step 5  NOW compute the metric — per replicate, at equal depth/capture.
Step 6  Scatter the raw data, stratified by sex and QC quartile.
Step 7  Vary the main hyperparameter 2×.
Step 8  Replicate on a second platform/species if available.
```

The order is the point. We spent three weeks and 500+ lines on MI before checking QC; the QC
check took 45 seconds and invalidated most of the work. The cost of the gates is about an hour.
The cost of skipping them is weeks and a false conclusion.

---

## Retroactive audit (v0.1 judgments)

| Metric | First failing gate (v0.1 logic) | Verdict |
|---|---|---|
| ACP eigenvalue / entropy | GATE 0 (variance-structure function) | not biology |
| Cardiac β | GATE 4 (physical constant) | real, but physics not biology |
| MI SMAD→ECM (zero-bin) | GATE 0 (first applicable), and 1, 2, 3, 6, with GATE 5 a partial fail | **not measuring biology** — under re-audit with the v0.3 logic (`validation/flagship_audit/`) |

Under v0.3, mi_3bin's dropout sensitivity is attenuation (reported), and the depth confound is
handled by the correction; on the bundled demo the raw MI difference (+0.026) is not detected
across mice and shrinks to −13% of itself at equal depth (INCONCLUSIVE, not supported).

*Re-read this before starting any new metric analysis.*
