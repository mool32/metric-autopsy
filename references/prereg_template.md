# Pre-registration = the elicitation form

*This file does double duty. As a **reference** it is the pre-registration you fill in before
writing any code. As a **skill behavior** it is the form the agent walks you through *before* it
computes anything — the whole point of `metric-autopsy` is to invert the default from
"compute → believe" to "state your commitments → red-team → then believe".*

Fill every field. Empty fields are not "TBD"; they are the reason the analysis will fail.

```
Metric name:                    ___
Hypothesis (one sentence):      ___
Biological process it measures: ___
Mathematical formula:           ___
Simplest non-biological explanation: ___
What would disprove it:         ___
```

## Commitments the engine enforces

These fields change what the engine computes, so they must be fixed before any outcome is
seen. Pass them as JSON (`--prereg prereg.json`, or `run_autopsy(prereg=...)`); explicit CLI
flags override the file. The engine hashes the whole object (`prereg_sha256`) into the report
and the run log.

```
Estimand:            composition | content     ← required; no estimand = UNIDENTIFIABLE
  composition = relative expression (correlations, module scores, MI of normalized counts):
                the deeper group is thinned to equal sequencing depth.
  content     = amount of RNA (total counts, genes detected, CytoTRACE-like scores):
                thinning to equal depth would delete the signal, so groups are thinned to
                equal spike-in (ERCC) capture; without spike-ins the comparison is
                UNIDENTIFIABLE.
Replicate unit (obs column):  ___  (mouse, donor, plate — the unit of inference, never the cell)
Minimum replicates per group: ___  (default 3; the graded rule below applies on top)
SESOI:               ___  smallest effect size of interest, on the construct scale
                          (needed for "no detectable effect" and for the power check)
Alpha / power:       ___ / ___  (defaults 0.05 / 0.8)
Signal direction:    increase | decrease   (for an injected-signal test)
Positive control pair: ___   Negative control pair: ___
Spike-in prefix:     ___  (default "ERCC-")
```

Graded replicate rule (fixed in the engine, not chosen per analysis):
- **>= 4 per group** (or pairs): an exact permutation over replicates is the test; the t
  interval estimates the effect.
- **3 per group**: the t interval on replicate-level values is the test and the verdict is
  marked *parametric only*; the permutation p is reported with the note that alpha is
  unattainable by construction.
- **<= 2 per group**, or no replicate unit: no effect verdict
  (`design_adequacy = INSUFFICIENT_REPLICATION`).

```json
{"estimand": "composition", "min_replicates": 3, "sesoi": 0.1, "alpha": 0.05,
 "power": 0.8, "signal_direction": "increase",
 "hypothesis": "Smad3-Col1a1 coupling declines with age in fibroblasts",
 "simplest_non_biological_explanation": "old cells are sequenced shallower"}
```

---

## GATE 0 — Mathematical independence
```
Variables in the formula:                       ___
Which are confounded by QC (sparsity, library size, variance, n)? ___
Simulation result (metric on nuisance-perturbed synthetic data):  ___
```

## GATE 1 — QC parity (a diagnostic)
```
Factorial obs columns to check (must include age × sex × batch × tissue × cell_type as available): ___
Strata confidently beyond 1.5× (bootstrap CI, Bonferroni):  ___
Strata with < min_cells per group (not assessable):         ___
```

## GATE 2 — Estimand-dependent correction
```
Estimand (from the commitments above):           ___
Correction applied (depth / spike-in capture):   ___
Raw effect ___  vs corrected effect ___  (retained ___ %)
Explained by depth/capture? (raw detected across replicates, corrected not): ___
```
(The v0.1 n_genes matching is deprecated: n_genes is moved by biology — cell size, cycling,
RNA content — so matching on it can delete a real effect.)

## GATE 3 — Visible in raw data
```
Scatter shows the effect?                ___
Still visible after stratifying by sex?  ___
Still visible after coloring by n_genes? ___
```

## GATE 4 — Alternative explanations (judgment) and construct response (auto)
```
Injected signal (e.g. coupling of the gene pair) and the response it must produce: ___
Non-biological scenarios that fit this result:
  1. ___
  2. ___
  3. ___
Tests that separate them from the biological explanation: ___
```

## GATE 5 — Controls (against empirical nulls)
```
Positive control pair ___   beats its self-shuffled null in every stratum? ___
Negative control pair ___   inside its expression-matched pair null?       ___
Checked per factorial combination?       ___
```

## GATE 6 — Replication
```
Independent dataset (platform / species): ___
Same estimand, correction, replicate unit and strata?  ___
Replicate-level result (REPLICATED / NOT_REPLICATED / INCONCLUSIVE): ___
```

## GATE 7 — Effect size (judgment)
```
Observed effect:                          ___
Expected for a known validated interaction: ___
Larger than test–retest variability?      ___
```

---

*Record the hash of the filled form in `PROJECT.md` before you look at any outcome. A pre-reg you
edit after seeing results is not a pre-reg. The engine helps: every report carries
`prereg_sha256`, and the run log (`metric_autopsy_runs.jsonl` for the CLI and MCP server)
counts how many times the same claim — same data, pre-registration and comparison — has been
run, so "re-run until it passes" is visible.*
