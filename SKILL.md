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
biology, and to stop at the first one it fails.

## Input contract

You need, from the user (elicit anything missing — do not guess):

- **Data:** an AnnData `.h5ad` (or a `metric_autopsy.SimpleData`) whose `obs` has the
  factor columns you will compare and stratify by — at minimum the grouping column
  (e.g. `age`) and the confounder axes (`sex`, `batch`, `tissue`, `cell_type`).
- **Metric:** one of the reference metrics (`mi_3bin`, `pearson`, `codetected_spearman`,
  `norm_pearson`, `spectral_entropy`) **or** a user callable `metric(data) -> float`.
- **Genes / groups:** the gene pair (for pairwise metrics), the group column and the two
  levels to compare, and positive/negative control gene pairs.

## Protocol — follow in order

### 1. Elicit the pre-registration *first* (before any computation)

Walk the user through `references/prereg_template.md`. Fill every field: the one-sentence
hypothesis, the biological process, the formula, **the simplest non-biological explanation**,
and what would disprove the claim. An empty field is not "TBD" — it is the reason the
analysis will fail. Record the filled form; it becomes the header of the autopsy report and
supplies the answers to the judgment gates (4 and 7).

Do not proceed to computation until the pre-reg is filled. This step is the whole point.

### 2. Run the auto gates in order, stopping at the first blocking failure

Invoke the engine — the thin wrapper is `scripts/run_gates.py`:

```bash
python scripts/run_gates.py --h5ad <data.h5ad> \
    --metric <name> --gene-a <A> --gene-b <B> \
    --group-col <col> --groups <g1> <g2> --within <factor...> \
    --pos-pair <A> <B> --neg-pair <A> <B> [--data2 <replicate.h5ad>]
```

or drive `metric_autopsy.run_autopsy(...)` directly for a bring-your-own callable. Try
`python scripts/run_gates.py --demo` to watch the reference `mi_3bin` die at GATE 0 (add
`--no-stop` to run every gate rather than halting at the first blocking failure).

The gates, and what each catches (full detail in `references/gates.md`):

| Gate | Catches | Auto? |
|---|---|---|
| **0 Mathematical independence** | metric moves under a nuisance stat (dropout, depth, variance) with no biology change | yes |
| **1 QC parity** | groups differ in technical quality — *in any factorial stratum*, not just the main axis | yes · **crown jewel** |
| **2 n_genes matching** | effect vanishes once cell quality is equalized; the "groups don't even overlap → STOP" verdict fires only when GATE 2 is run *within* a flagged stratum | yes · **crown jewel** |
| **3 Raw visibility** | effect isn't visible in the raw scatter; a "shape change" is really dropout | export + judgment |
| **4 Measures what you think** | more than one non-biological scenario explains the result | judgment — you ask |
| **5 Controls** | positive control silent or negative control fires, in any stratum | yes |
| **6 Replication** | effect doesn't survive an independent platform/species after QC matching | yes, if 2nd dataset |
| **7 Effect size** | statistically real but biologically negligible | judgment — you ask |

A `FAIL` or `STOP` is **blocking**: do not let a later gate rescue it. Report the death and
stop. GATE 1 is the highest-yield check — always run it first, stratified by every factor.

### 3. Resolve the judgment gates (4, 7) with the user

These are not computable from data. Ask directly, using the pre-reg answers: *Which
non-biological scenarios in the GATE 4 table fit your result, and what test separates them?*
and *Is the effect larger than test–retest variability and in range for a real regulatory
interaction?* If the user cannot rule the alternatives out, the verdict is INCONCLUSIVE, not
PASS.

### 4. Emit the autopsy — numbers and a verdict, no rescue language

Print the gate-by-gate table (`Autopsy.to_markdown()`) with the actual numbers, then a
single verdict decided by the first blocking gate. Follow the research-flow rule: a failed
metric is recorded as **killed**, plainly, with the number that killed it — never softened.
Cleared gates yield a *provisional* PASS ("real until replicated"), not a victory lap.

## Progressive disclosure

Keep this file thin. Load depth on demand:
- `references/gates.md` — full definition and rationale of every gate (also the preprint methods).
- `references/red_flags.md` — the 12-row red-flag table, each with a real example.
- `references/prereg_template.md` — the elicitation form for step 1.

## What "good" looks like

The reference failure: `mi_3bin` on SMAD→ECM dies at GATE 0 (expectation shifts ~60% under
dropout) and GATE 1 (male stratum 1.9× QC gap). GATE 1's per-stratum flag is the cue to act:
run pooled, GATE 2 finds nothing (the apparent effect is ~0 — pooling old = male+female dilutes
the male-only confound away), so the automated run reports GATE 2 PASS. Re-run GATE 2 *within*
the flagged male stratum and it STOPs — young/old n_genes don't even overlap there, so the
groups are incomparable. The confound is invisible pooled and lethal stratified: **stratify, or
the artifact hides in the interaction** — GATE 2 is only decisive once you restrict it to the
stratum GATE 1 flags.
