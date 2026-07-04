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

---

## GATE 0 — Mathematical independence
```
Variables in the formula:                       ___
Which are confounded by QC (sparsity, library size, variance, n)? ___
Simulation result (metric on nuisance-perturbed synthetic data):  ___
```

## GATE 1 — QC parity
```
Factorial obs columns to check (must include age × sex × batch × tissue × cell_type as available): ___
Max QC ratio found across combinations:  ___
Unusable combinations (> 1.5×):          ___
```

## GATE 2 — n_genes matching
```
Overlap range (10th–90th pct):                          ___
Matched sample sizes:                                   ___
Effect with matching ___  vs without ___
Effect preserved (retained >= 50% of unmatched)?        ___
Effect not amplified (retained <= 300% / 3x)?           ___
Matched subset balanced (median n_genes ratio <= 1.1x)? ___
```

## GATE 3 — Visible in raw data
```
Scatter shows the effect?                ___
Still visible after stratifying by sex?  ___
Still visible after coloring by n_genes? ___
```

## GATE 4 — Alternative explanations (judgment)
```
Non-biological scenarios that fit this result:
  1. ___
  2. ___
  3. ___
Tests that separate them from the biological explanation: ___
```

## GATE 5 — Controls
```
Positive control pair ___   result ___
Negative control pair ___   result ___
Checked per factorial combination?       ___
```

## GATE 6 — Replication
```
Independent dataset (platform / species): ___
Result after QC matching:                 ___
```

## GATE 7 — Effect size (judgment)
```
Observed effect:                          ___
Expected for a known validated interaction: ___
Larger than test–retest variability?      ___
```

---

*Record the hash of the filled form in `PROJECT.md` before you look at any outcome. A pre-reg you
edit after seeing results is not a pre-reg.*
