---
name: metric-autopsy
description: Use when validating a computed metric on single-cell RNA-seq data (mutual information, correlation, coupling, entropy, eigenvalue ratios, or any bring-your-own metric) to check whether an apparent biological signal is actually a QC, technical, or mathematical artifact before making a claim. Trigger on "is this signal real", "validate this metric", "could this be a batch/dropout/library-size artifact", "why does my effect vanish on 10x", or before writing up any scRNA metric result.
license: MIT
---

# Metric autopsy

A metric that changes between conditions is not yet a finding. This skill **inverts the
default from "compute → believe" to "state your commitments → red-team → then believe."**
It exists because three real analyses (entropy anticorrelation, cardiac β, SMAD→ECM mutual
information) each survived months of work before a 45-second QC check killed them.

Your job when this skill runs is **not** to compute a number and report it. It is to run
the metric through a gauntlet of gates, each designed to catch one way a metric fakes
biology, and to report the verdict the engine decides from four fields — metric validity,
design adequacy, effect, replication — with the numbers behind each.

> **Status: v0.x under validation.** The validator has known failures in both directions
> (`validation/probes/README.md`). Never present a verdict from this version as validated.

## Input contract

You need, from the user (elicit anything missing — do not guess):

- **Data:** an AnnData `.h5ad` (or a `metric_autopsy.SimpleData`) whose `obs` has the
  factor columns you will compare and stratify by — at minimum the grouping column
  (e.g. `age`) and the confounder axes (`sex`, `batch`, `tissue`, `cell_type`).
- **Metric:** one of the reference metrics (`mi_3bin`, `pearson`, `codetected_spearman`,
  `norm_pearson`, `spectral_entropy`) **or** a user callable `metric(data) -> float`.
- **Genes / groups:** the gene pair (for pairwise metrics), the group column and the two
  levels to compare, and positive/negative control gene pairs.
- **Replicates:** the `obs` column of the biological replicate (mouse, donor, plate). Without
  it there is no effect verdict — cells are not independent replicates.

## Protocol — follow in order

### 1. Elicit the pre-registration *first* (before any computation)

Walk the user through `references/prereg_template.md`. Fill every field: the one-sentence
hypothesis, the biological process, the formula, **the simplest non-biological explanation**,
and what would disprove the claim. Then the commitments the engine enforces:

- **estimand** — `composition` (relative expression) or `content` (amount of RNA). It decides
  the correction: thinning to equal depth, or to equal spike-in capture. No estimand, or a
  content estimand without spike-ins, is UNIDENTIFIABLE;
- **direction** — the claimed change of the metric from the first group to the second
  (`increase`, `decrease`, or `two-sided` for a non-directional claim, which the verdict marks).
  SUPPORTED needs the effect in this direction; a detected effect the other way is NOT SUPPORTED;
- **replicate unit** and **minimum replicates** per group;
- **SESOI** — the smallest effect size of interest (needed to claim "no detectable effect",
  for the power check, and to size nuisance biases: a dropout or library-size bias blocks only
  above `bias_tolerance` × SESOI, default 0.5; without a SESOI SUPPORTED is withheld while a
  bias is unsized);
- alpha, power, controls (or an injected signal and its direction).

An empty field is not "TBD" — it is the reason the analysis will fail. Save the commitments as
JSON for `--prereg`; the engine hashes them into the report. The form supplies the answers to
the judgment gates (4 and 7).

Do not proceed to computation until the pre-reg is filled. This step is the whole point.

### 2. Run the gates

Invoke the engine — the thin wrapper is `scripts/run_gates.py`:

```bash
python scripts/run_gates.py --h5ad <data.h5ad> \
    --metric <name> --gene-a <A> --gene-b <B> \
    --group-col <col> --groups <g1> <g2> --within <factor...> \
    --replicate-col <mouse> --prereg <prereg.json> \
    --pos-pair <A> <B> --neg-pair <A> <B> [--inject-signal coupling] \
    [--data2 <replicate.h5ad>] --json autopsy.json
```

or drive `metric_autopsy.run_autopsy(...)` directly for a bring-your-own callable. Try
`python scripts/run_gates.py --demo --no-stop` to see the four fields on the bundled
synthetic confound.

The gates, and what each catches (full detail in `references/gates.md`):

| Gate | Catches | Feeds |
|---|---|---|
| **0 Mathematical independence** | a nuisance (dropout, depth, library size) that *biases* the metric — FAIL beyond the SESOI tolerance; a depth bias the declared correction removes is reported; attenuation is reported, not failed | metric validity; attenuation → power |
| **1 QC parity** | groups differ in technical quality — *in any factorial stratum*; a diagnostic (WARN), not a kill switch | design adequacy |
| **2 Estimand-dependent correction** | the raw difference is depth or capture: thinning to equal depth (composition) or spike-in capture (content), then replicate-level inference | effect |
| **3 Raw visibility** | effect isn't visible in the raw scatter; a "shape change" is really dropout | export + judgment |
| **4 Measures what you think** | the metric does not respond to an injected signal (auto); more than one non-biological scenario explains the result (judgment — you ask) | metric validity |
| **5 Controls** | negative control fires, or positive control silent where the design had the power to show it (FAIL); positive control silent where it had not (WARN); per stratum, against empirical nulls; the controls must run the metric under judgment | metric validity |
| **6 Replication** | the corrected, replicate-level effect does not hold on independent data | replication |
| **7 Effect size** | statistically real but biologically negligible — declare it as the SESOI | judgment — you ask |

An invalid metric (GATE 0 bias, a failed negative control, no response to an injected
signal) is **blocking**: by default the effect is not even estimated. Do not let a later gate
rescue it. GATE 1 is the highest-yield look at the data — always stratify by every factor.

### 3. Resolve the judgment gates (4, 7) with the user

These are not computable from data. Ask directly, using the pre-reg answers: *Which
non-biological scenarios in the GATE 4 table fit your result, and what test separates them?*
and *Is the effect larger than test–retest variability and in range for a real regulatory
interaction?* Only when the user rules the alternatives out, re-run with
`--resolve-judgment` (or `prereg["judgment_pending"] = False`). Otherwise the verdict stays
INCONCLUSIVE.

### 4. Emit the autopsy — numbers and a verdict, no rescue language

Print `Autopsy.to_markdown()`: the verdict, the four fields with their reasons, the
gate-by-gate table, and the provenance hashes. Keep the JSON report (`--json`) with the data
and pre-registration hashes. The verdict is one of SUPPORTED (provisional until replicated) ·
SUPPORTED — replicated · NOT SUPPORTED · NO DETECTABLE EFFECT · INCONCLUSIVE · UNIDENTIFIABLE ·
DEGENERATE METRIC, with its qualifiers (*parametric only*, *underpowered*). Follow the
research-flow rule: a failed metric is recorded as **killed**, plainly, with the number that
killed it — never softened. SUPPORTED is provisional ("real until replicated"), not a victory
lap. If the run log shows earlier attempts at the same claim, say so.

## Progressive disclosure

Keep this file thin. Load depth on demand:
- `references/gates.md` — full definition and rationale of every gate (also the preprint methods).
- `references/red_flags.md` — the 12-row red-flag table, each with a real example.
- `references/prereg_template.md` — the elicitation form for step 1.

## What "good" looks like

The reference case, on the bundled demo (`--demo --no-stop`): `mi_3bin` on SMAD→ECM, biology
identical everywhere, male-old capture degraded. GATE 0 finds no bias but heavy attenuation
(dropout −61%, depth halving −24%), which goes to the power check. GATE 1 flags the male
stratum (1.94× QC ratio, n_genes overlap 0.00) and the correction thins young males to the old
males' depth. Across 16 mice the raw MI difference (+0.026) is not detected (p = 0.52) and
shrinks to −13% of itself at equal depth: **INCONCLUSIVE — not supported**. Pooled, the male-only
confound is invisible; stratified, it is found and removed. **Stratify, or the artifact hides in
the interaction.**
