# Diagnostics on the data of the pilot's second run (the plan, committed before the run)

The owner's decision of 2026-10-17 (item 1): before the code of v0.4, four checks run on the data of
the pilot's second run. **Status: exploratory, on development data.** These data stay outside the
confirmatory set of the real track (`validation/plans/real_claims.md`) for good. One attempt; the
outputs go to the branch `results/real-pilot-3`.

## What it is for

One question to every check: **does it separate the drop of MI on §4 from the sex differences?** A
check that flags §4's claim and lets Xist and the Y genes pass is a candidate mechanism for v0.4: the
translation of the SESOI into the metric's units is approved only together with a mechanism that
passes this (`validation/plans/v0.4_spec.md`, change 3), and the decision on GATE 5 waits for it
(change 13). Without such a mechanism a SESOI in the metric's units is tiny for an insensitive
metric: MI's response to GATE 4's coupling at strength 1 was +0.0005 against a drop of 0.020 in
males, and the claim on males alone (p = 0.016) would get SUPPORTED.

## The data and the engine

- **The second run's files, unchanged,** from `results/real-pilot-2` (commit bad785a):
  `data/case1.npz` (its parts `case1.npz.part-00` to `-03` concatenated; sha256
  `13e7a563d4ed8fa482f976ead77b492189e8c729eb3ab103a203ccfdf8f24a44`), `data/case2.npz` (sha256
  `7191ad731063a45110991a7fa7d3b8c15d43325bf7e9b6d95072165dc1655beb`) and `select/select.json`
  (sha256 `63032110bb981a7c1615d680848d116271b1a52c01bee1e9e0b1c585634493ac`; the cases' groups,
  strata and units). The hashes are checked before anything is computed. Nothing is read from the
  Census again.
- **The engine of the tag `v0.3.0-prereg`,** installed from a worktree with `src/` checked against
  its git tree `97cbab3c437c9ddb118a0690e6e0d24f8a16069f`, in v1's analysis environment with v1's
  numerical settings, as in the second run. Check (i) runs `run_autopsy`. Checks (ii) to (iv) use
  the engine's own pieces: the metrics (`metrics.mi_3bin`, `metrics.norm_pearson`, and R1's
  `anchors.log_cp10k_mean`), the replicate-level statistics (`stats.nested_t`,
  `stats.permutation_test_nested`) and the injectors (`injected_signal.coupling`,
  `injected_signal.module`). They compute a metric on the genes it reads and one more column, the
  rest of each cell's counts: the metrics read nothing else (the named genes and the row sum), so the
  values are the full matrix's; the script checks this on every claim and mouse, and compares the
  per-mouse values with the second run's (b).
- **STRING v12.0, mouse (taxon 10090),** for check (ii): `10090.protein.info.v12.0.txt.gz` and
  `10090.protein.links.v12.0.txt.gz` from `https://stringdb-downloads.org/download/`, downloaded by
  the workflow (this container cannot reach the site). Their sha256 are recorded.

## The claims

| id | data | metric | groups (first → second) | strata | mice | direction |
|---|---|---|---|---|---|---|
| 1a | case 1, all cells | `mi_3bin(Smad3, Col1a1)` | young → old | sex | females 5 → 3, males 5 → 4 | decrease |
| 1a-males | case 1, male cells | `mi_3bin(Smad3, Col1a1)` | young → old | none | 5 → 4 | decrease |
| 1b | case 1, all cells | `norm_pearson(Smad3, Col1a1)` | young → old | sex | as 1a | decrease |
| 2a | case 2 | `log_cp10k_mean(Xist)` | female → male | none | 3 → 4 | decrease |
| 2b | case 2 | `log_cp10k_mean(Ddx3y + Eif2s3y + Kdm5d + Uty)` | female → male | none | 3 → 4 | increase |

The mouse is the unit. **A claim's difference Δ** is that of the second run's analysis (b): the
metric on the cells of each mouse; Δ is the second group's minus the first's; over several strata
the engine's pooled-stratified estimate (`stats.nested_t`), in one stratum the difference of the
group means. A stratum's own Δ is the difference of its group means. "Detected" in a stratum means
the exact permutation of the group labels within it gives p < 0.05, as in the second run.

## The checks

### (i) The engine with GATE 4 and GATE 5 off: the raw and the corrected effect

`run_autopsy` with GATE 0, GATE 1 and the correction of depth (GATE 2) on, and neither an injected
signal nor controls: as the second run's (d), but `signal_test=None`, no control pairs, and
`stop_on_first_fail=False` so that the effect is estimated whatever GATE 0 says (the verdict is the
engine's either way). The composition estimand, the claim's direction, `judgment_pending=False`,
seed 0, run log off. The SESOI: for 1a and 1a-males the design's floor on the observed scale, as for
the second run's 1a-floor (SE × (t₀.₉₅,df + t₀.₉₇₅,df) of the claim's own uncorrected difference),
because 0.1 is wrong in scale for `mi_3bin`; for 1b 0.1; for 2a and 2b 0.5 (R1's). The SESOI
enters only the TOST behind `explained_by_depth`, not the outcome below.

Reported: the raw and the corrected effect, `p_raw`, `p_perm`, `raw_detected`, `retained`,
`explained_by_depth`, `reversed_by_correction`, the TOST, the depth ratios of the thinning, the
gates, the verdict and its cause.

**Outcome.** *Not applicable* where the engine does not detect the raw difference
(`raw_detected` false) or does not estimate the effect. *Flags* (the correction of depth accounts
for the difference) where the raw difference is detected and the corrected effect is not detected
in the raw's direction: `explained_by_depth`, `reversed_by_correction`, or not detected at equal
depth; which of these is reported. *Does not flag* where the corrected effect is detected with the
raw's sign.

### (ii) Specificity: the same donor-level difference for 100 pairs (genes) of similar detection

**The rule of selection, fixed here.** The detection of a gene is the share of the case file's
cells in which its count is above 0 (all cells of `case1.npz` for 1a, 1a-males and 1b; all cells of
`case2.npz` for 2a and 2b). The candidates are the file's genes except:

- the claim's genes; in case 1 also §4's controls (Actb, Gapdh, Tst, Lrrc42);
- genes with a known difference between the sexes, a fixed list (not a chromosome annotation): the
  X-inactivation centre (Xist, Tsix, Jpx, Ftx); Y genes (Ddx3y, Eif2s3y, Kdm5d, Uty, Zfy1, Zfy2, Sry,
  Usp9y, Ube1y1); their X homologues that escape inactivation (Ddx3x, Eif2s3x, Kdm5c, Kdm6a, Zfx,
  Usp9x, Uba1); and pseudoautosomal genes (Mid1, Erdr1, Sts);
- for the pairs of case 1, genes whose symbol is no preferred name of STRING v12.0 mouse (exact,
  case-sensitive): their interactions are unknown, not absent.

For 2a and 2b, a single gene and a sum of genes, "without known interactions" is read as "without a
known difference between the sexes": the list above. Each draw starts its own
`numpy.random.default_rng(20261017)`, over candidates listed in the file's gene order.

- **1a, 1a-males and 1b: 100 pairs, the same for the three claims.** For Smad3 and for Col1a1, the
  50 candidates nearest by |detection difference| (ties by the file's order). The eligible pairs:
  one gene from Smad3's 50 and one from Col1a1's, two different genes, without a STRING link of
  combined score ≥ 400 between any of their proteins, each unordered pair once. 100 are drawn
  without replacement (all, if fewer are eligible: reported).
- **2a: 100 genes**, drawn without replacement from the 200 candidates nearest to Xist by detection.
- **2b: 100 sets of four genes**, one from each Y gene's 50 nearest candidates, four different
  genes, each set once: drawn by repetition, one gene uniformly from each of the four lists, a set
  with a repeated gene or drawn before skipped. The metric sums the four as 2b sums the Y genes.

**The statistic** of every control is the claim's Δ, on the same cells, mice and strata, with the
claim's metric.

**Outcome.** The claim's position among the controls in its direction: p = (1 + the number of
controls whose Δ is at least as far in the claimed direction as the claim's) / 101. *Flags* (the
difference is not specific to the claim's genes) where p ≥ 0.05; *does not flag* where p < 0.05
(the claim is beyond at least 96 of the 100). Reported beside it: the two-sided position, the
controls' quantiles of Δ, their detection against the claim genes', and their first group's mean
against the claim's (the controls are matched on detection, not on the level of the metric).

### (iii) Plausibility: the response to the extremes of the construct against the observed difference

Per mouse, on its own cells, with the engine's injectors and the response as GATE 4 defines it (the
metric on injected data minus the metric on sham data, at the same thinning):

- **pair claims (1a, 1a-males, 1b):** *the maximal dose of coupling* is the largest group-mean
  response over coupling strengths 1, 3, 10 and 50 (at 50 the shared keep-probability is near 0 or
  near 1 in most cells; the script logs the share), because the response need not grow with the
  dose (`norm_pearson` reads co-detected cells, which an all-or-nothing factor leaves unscaled); 20
  injections and 20 shams per strength and mouse. *The full decoupling* is Col1a1's counts permuted
  across the mouse's cells, 20 permutations: the metric then minus the metric observed.
- **level claims (2a, 2b), by analogy:** the maximal dose is the largest group-mean response of the
  module probe over folds 2, 10 and 100, in every cell (20 each); the full loss is the genes' counts
  set to 0 (permuting one gene does not change a mean).

The draws come from generators seeded with 20261017, the measure and the mouse, so that 1a and
1a-males share their males' responses.

In every stratum, with h the group whose mean of the claim's metric is higher and l the lower:
(a) the decoupling (loss) lowers h's mean by at least |Δ| of the stratum; (b) the maximal dose
raises l's mean by at least |Δ|. **Outcome.** *Flags* (the difference is larger than the construct
can produce) where (a) or (b) fails in any stratum of the claim; *does not flag* otherwise. Both
parts and their ratios to |Δ| are reported.

### (iv) Agreement of measures: the sign and size of the change by two measures, by stratum

- **Pair claims:** `mi_3bin` and `norm_pearson` of the same pair (the other measure of 1a and
  1a-males is `norm_pearson`; of 1b, `mi_3bin`).
- **Level claims, by analogy:** the claim's metric and the share of the mouse's cells in which the
  genes (summed) are detected.

Per stratum: each measure's Δ, its exact permutation p, and its size standardized by the pooled
within-group SD of the per-mouse values. **Outcome.** *Flags* (the measures disagree) where, in any
stratum of the claim, the two Δ have opposite signs and at least one is detected; *does not flag*
otherwise. The sizes are reported; they do not enter the outcome.

## Outcomes and "separates"

Every check gives every claim one outcome: flags, does not flag, or not applicable. A check
**separates** §4 from the sex differences where it flags 1a and 1a-males and flags neither 2a nor
2b; **partly** where it flags one of 1a and 1a-males and neither 2a nor 2b; it **does not
separate** where it flags neither 1a nor 1a-males, or flags 2a or 2b. 1b's outcomes are reported
with the rest and do not enter this.

## Known before the run

From the second run's outputs (`validation/exploratory/real_pilot_2.log`), so not outcomes of this
run:

- **(i), 1a: not applicable.** The engine's raw test of 1a is (b)'s exact permutation over the same
  per-mouse values: p = 0.0536 in the second run, not detected. 1a-males: males alone gave
  p = 0.016 (exact, 126 relabelings), so the raw difference is expected detected; its correction
  was never computed (GATE 4 failed first).
- **(i), 1b: does not flag,** if the effect field reproduces the second run's (d), as expected (it
  draws from its own generator seeded 0, and its permutation is exact): corrected −0.3664 (first
  minus second) detected, raw −0.3687 detected, the same sign. 2a and 2b: never computed.
- **(iv), 1a, 1a-males and 1b: flags.** In males MI fell by 0.0197 (p = 0.016) and the Pearson rose
  by 0.1305 (p = 0.048): opposite signs, both detected. In females both rose (MI +0.0017, p = 0.68;
  Pearson +0.651, p = 0.089). 2a and 2b: not computed.
- **(iii), 2a and 2b: expected to flag, by construction.** The lower group has about no counts of
  the genes (per mouse, Xist 0 in 3_10_M and 0.0001 in 3_11_M; the Y genes 0 in 3_39_F and
  3_56_F). Thinning cannot create counts, so the response of l to any dose is about 0 and part (b)
  fails. Part (a) cannot fail for these claims (the loss takes h's mean to 0). For 1b, the young
  females' mean correlation is negative (−0.043): decoupling takes a correlation to about 0, so part
  (a) is expected to fail in females. For `mi_3bin`, part (a) fails only where l's mean is below
  h's decoupled MI.
- **(ii), 1b: expected to flag.** Its Δ (+0.369) is against the claimed decrease; it is beyond 96
  of the controls in the claimed direction only if 96 controls rise by more. Otherwise nothing of the
  controls has been computed.

## The prediction (the owner's, 2026-10-17, verbatim)

> если §4 — артефакт детекции, его отметят проверки (ii) и (iv), а Xist и Y-гены пройдут

("If §4 is an artifact of detection, checks (ii) and (iv) will flag it, and Xist and the Y genes
will pass.") Scored as: (ii) and (iv) each flag 1a and 1a-males and flag neither 2a nor 2b. The part
of (iv) on §4 is known before the run (above); the report says so.

## One attempt

`.github/workflows/real-pilot-3.yml`, started by the label `pilot-run-3` on the pilot's pull
request, commits whatever the checks produce, failures included, to the new branch
`results/real-pilot-3`; a second start finds the branch and stops. If the preparation fails before
any check has run (the second run's files, their hashes, the STRING download), it may be fixed and
started again, and the report says so. Nothing is tuned after the run: no rule above changes in
the light of its numbers.

## What is reported

`real_pilot_3.log` and `real_pilot_3.json`: for every check and claim the numbers and the outcome;
for every check whether it separates; the prediction's score; the controls of (ii) with their
detection and STRING status; the hashes of the inputs; the engine's reports of (i). One report to
the owner, by check and by claim: separates or not.

## Not in this

The code of v0.4, any change of the engine, a confirmatory claim, the first run's case 3.
