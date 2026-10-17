# A second diagnostics on the data of the pilot's second run (the plan, committed before the run)

The owner's decision of 2026-10-17, on the first diagnostics' report (`real_pilot_3.log`), item 2.
**Status: exploratory, on development data**, outside the confirmatory set of the real track for
good. One attempt; the outputs go to the branch `results/real-pilot-4`.

## What it is for

The changes 3, 4 and 13 of v0.4 wait for it. Two checks:

- **(v) the composition of cells on §4** (change 7 tried on §4 itself): does the drop of MI in
  males disappear within cell types?
- **(ii′) specificity again,** on the standardized difference and with the full pool of genes: are
  Xist and the Y genes extreme among matched controls, and §4 in males typical?

## The data and the engine

- **The second run's files,** unchanged and checked by their hashes as in the first diagnostics
  (`real_pilot_3_inputs.sh`: `results/real-pilot-2` at bad785a): `case1.npz`, `case2.npz` and
  `select.json`. (v) reads `case1.npz` alone: its cells carry the Census's `cell_type` and `tissue`.
- **For (ii′), the full set of genes.** The case files hold D1's cut of 2,000 genes. The counts of
  every gene that the Census marks as measured in every dataset of the case (21,025 genes for case
  1, as the second run's log counts them) are read again from the CELLxGENE Census, release
  2025-11-08, for the same cells (by their `soma_joinid` in the case files), 10,000 cells at a time,
  and checked to be raw counts (`real_pilot_select.read_blocks`). Only the claims' cells are read:
  the males of case 1 (for 1a-males) and the cells of case 2. Symbols that name more than one
  feature are left out of the pool.
- **STRING v12.0, mouse**, the files of the first diagnostics (sha256
  `ededc5cefddea16013c4d8f7d5a1080efea22d4b3f1f090c27019eb1d5209c08` and
  `803775e1d1a47864561e41944bf1f988094a836991e7cf13e4229233fc8cff3d`), downloaded again and
  checked against these hashes.
- **The engine of the tag `v0.3.0-prereg`** for the metrics (`metrics.mi_3bin`,
  `metrics.norm_pearson`, R1's `anchors.log_cp10k_mean`) and the statistics (`stats.nested_t`,
  `stats.permutation_test_nested`, `stats.clopper_pearson`), installed as in the earlier runs.
  `run_autopsy` does not run.

## (v) The composition of cells on §4

The cells of claim 1a (case 1: 76,797 cells, 17 units), within each sex. The units are those of the
second run, the two composite ids included.

- **Per unit:** the share of its cells of each cell type and of each tissue (the Census's
  `cell_type` and `tissue`).
- **The difference of composition, per sex:** for each type (and each tissue) the two groups' mean
  shares (each unit weighted equally) and their difference, old minus young, with its exact
  permutation p (descriptive, no correction for the number of types). As one number per sex, the
  total variation distance between the groups' mean compositions, TVD = ½ Σ |mean share old −
  mean share young| over the types, with its exact permutation p (every relabeling of the sex's
  units into the groups' sizes: 126 in males, 56 in females). The same for the tissues.
- **MI within cell types:** per sex, the types with at least 20 cells in at least 3 units of each
  group. In each such type, `mi_3bin(Smad3, Col1a1)` on the cells of the type in each unit that has
  at least 20 of them; Δ_type is the old units' mean minus the young units' mean, with the exact
  permutation p over the type's units. **Combined:** Δ_within, the mean of Δ_type over the
  qualifying types with equal weights. Its exact permutation p is over every relabeling of the
  sex's units, with each type's Δ recomputed from its own units under that relabeling. A type left
  with an empty group under a relabeling drops out of that relabeling's mean. The p is two-sided
  (|Δ| at least as large). `norm_pearson` of the pair is computed the same way and reported beside
  MI (v0.4 change 15), not scored.
- **A mouse's MI and its share of fibroblast-like cells.** Fibroblast-like: the cells whose Census
  `cell_type` label contains "fibroblast", "mesenchymal" or "stromal" (any case). This is fixed
  here, by the annotation and not by expression, since Col1a1 is a gene of the claim. The labels it
  takes are listed with their cells. A unit's MI is the second run's (b) value, over all its cells.
  Reported: Pearson r and Spearman ρ across the units of each sex, with exact permutation p (every
  ordering: 9! in males, 8! in females). Males are primary.
- **The outcome:** composition accounts for the drop in males where Δ_within is not detected
  (p ≥ 0.05) and is less than half the overall drop (Δ_within > ½ × Δ, Δ the males' overall
  difference of the second run, −0.0197): "not detected" alone could come from the fewer cells
  within a type.

**The prediction (the owner's, 2026-10-17, verbatim):**

> внутри типов падение MI у самцов исчезает

("Within cell types the males' drop of MI disappears.") Scored by the outcome above: it holds where
composition accounts for the drop in males.

## (ii′) Specificity again, on the donor-level t, with the full pool of genes

Claims 1a-males (`mi_3bin` of Smad3 and Col1a1, young → old, the males of case 1), 2a (Xist) and 2b
(the four Y genes summed), female → male in case 2.

- **The statistic** of the claim and of every control: t = Δ / SE. Δ is the donor-level difference
  (second group minus first) and SE its standard error, both from `stats.nested_t` (Welch, one
  stratum). A control whose SE is 0 has no t: it is counted and listed, not ranked.
- **The claims' values:** 1a-males per mouse as before (they must equal the second run's (b): the
  counts of Smad3 and Col1a1 are read again). For 2a and 2b the CP10k's denominator is now the sum
  of the cell's counts over the full set of genes, the same for the claim and its controls. The
  second run's values, whose denominator was D1's 2,000 genes, are reported beside them.
- **The pool:** every gene of the full set with at least one count in the claim's cells, except
  the claim's genes, §4's controls (Actb, Gapdh, Tst, Lrrc42; case 1), the first diagnostics' list
  of genes known to differ between the sexes, and, for the pairs of 1a-males, the genes whose
  symbol is no STRING v12.0 preferred name.
- **Matching on mean and detection:** for every gene, m is its mean count per cell and d the share
  of cells with a count, over the claim's cells (the males of case 1; the cells of case 2). The
  distance to a claim gene is √(((log₁₀ m − log₁₀ m₀) / s_m)² + ((d − d₀) / s_d)²), where s_m and
  s_d are the SDs of log₁₀ m and of d over the pool. Ties are broken by the Census's order of
  features.
- **The controls, K = 1,000 per claim.** Each draw starts its own
  `numpy.random.default_rng(20261018)`, over candidates listed in the Census's order.
  - 1a-males: for Smad3 and for Col1a1, the 100 nearest. The pairs take one gene from each list,
    two different genes, without a STRING link of combined score ≥ 400 between any of their
    proteins, each unordered pair once; 1,000 are drawn without replacement.
  - 2a: 1,000 drawn without replacement from the 2,000 nearest to Xist.
  - 2b: 1,000 sets of four, one from each Y gene's 100 nearest, four different genes, each set once.
    They are drawn by repetition (one gene uniformly from each list; a set with a repeated gene, or
    drawn before, is skipped). The metric sums the four as 2b sums the Y genes.
- **Positions.**
  - In the claimed direction: k₁ is the number of controls whose t is at least as far in that
    direction as the claim's.
  - By absolute value: k₂ is the number with |t| at least the claim's.
  - Each share k / n, n the controls with a t (K less those without one), has its 95%
    Clopper-Pearson interval. **Specific** where the interval lies below 0.05; **not specific**
    where it lies at or above 0.05; **undecided** where it crosses 0.05 (the owner's rule for
    GATE 5, 2026-10-17).
  - Reported beside them: p = (1 + k) / (1 + n), the controls at least as far (both ways, with their
    genes and t), and the matched genes' mean and detection against the claim genes'.

**The prediction (the owner's, 2026-10-17, verbatim):**

> Xist и Y-гены — крайние в обе стороны, §4 по самцам — типичный

("Xist and the Y genes extreme both ways, §4 in males typical.") Scored as: 2a and 2b specific
both in the claimed direction and by absolute value; 1a-males not specific in the claimed
direction.

## Known before the run

From the earlier runs' logs, so not outcomes of this run:

- **The units' tissues** (`real_pilot_2_select.log`). The old females' 728 cells are of the mammary
  gland only, so within types the females compare mammary types. The males' units span 16 to 20
  tissues, except the composite 3_10_M/3_11_M (63 cells, one tissue). The cell types per unit have
  not been tabulated before.
- **The claims' differences and errors** (`real_pilot_3.log`): 1a-males Δ −0.019747, SE 0.005327
  (df 4.442); 2a Δ −0.2215, SE 0.04034 and 2b Δ +0.67001, SE 0.02188, both with D1's denominator.
- **The first diagnostics' (ii)** on raw differences with D1's 2,000 genes. 8 of 100 matched pairs
  fell in males as far as Smad3–Col1a1, several of them pairs of fibroblast and matrix genes
  (Htra1–Ogn, Lum–Aebp1, Gas1–Sfrp1). Xist: 0 of 100. The Y genes: 2 of 100.

## One attempt

`.github/workflows/real-pilot-4.yml`, started by the label `pilot-run-4` on the pilot's pull
request, commits whatever the checks produce, failures included, to the new branch
`results/real-pilot-4`; a second start finds the branch and stops. If the preparation or the
extraction from the Census fails, no check starts and nothing is committed: it may be fixed and
started again, and the report says so. Nothing is tuned after the run.

## What is reported

`real_pilot_4.log` and `real_pilot_4.json`:
- for (v), the composition tables, the composition differences, MI and `norm_pearson` within types,
  the association, the outcome and the prediction's score;
- for (ii′), the positions, the outcomes, the controls at least as far, and the prediction's score;
- the hashes of the inputs.

One report to the owner.
