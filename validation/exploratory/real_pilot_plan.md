# The pilot on real data with v0.3.0's engine (the plan, committed before the run)

The owner's decisions of 2026-10-14 (item 5) and 2026-10-15: v0.3.0 is not released before a check
on real data; the pilot runs now, in a pull request of its own based on the branch of pull request
1, and the code of v0.4 follows its report. Every dataset the pilot reads is development data and
stays outside the confirmatory set of the real track (`validation/plans/real_claims.md`) for good.

## What it is for

Not a test of the engine. Three development cases on real data show what the real track's runs do
before v0.4 is coded: the conclusion of (a), (b) and (d) on each claim; where (d) changes (b)'s
conclusion, and by which gate; the running time; what breaks on real data; what the specification
of v0.4 and the rule of rightness need.

## The engine, the environment, the data

- **(d) is the engine of the tag `v0.3.0-prereg`** (its `src/` git tree
  `97cbab3c437c9ddb118a0690e6e0d24f8a16069f`, checked by the workflow, installed from a worktree of
  the tag). Every SUPPORTED is also read under the known issue of v0.3.0: with `raw_detected` false
  it reads INCONCLUSIVE. Both are reported.
- **Environments:** v1's, pinned, on ubuntu-24.04 with CPython 3.13.16: `requirements-select.txt`
  to read the Census, `requirements-panel.txt` for the analyses, with v1's numerical settings (one
  BLAS thread, the same CPU dispatch).
- **Data:** the CELLxGENE Census, release 2025-11-08 (v1's), read in GitHub Actions
  (`.github/workflows/real-pilot.yml`). The extraction does not import the engine.
- **One attempt.** The workflow commits its outputs to the branch `results/real-pilot`; a second
  start finds the branch and stops. Once the analyses have started, whatever they produced is
  committed, failures included. A fix of the extraction before any analysis has run is recorded.

## Loading, as deviation D1 loads a background

For each case: the cells of one cell type, at most 500 per donor and condition (drawn with the seed
20261015 where there are more), dense. The genes: the case's named genes first, then the most
expressed by mean count over the extracted cells, 2,000 in all, among the genes the Census marks
as measured in every dataset the case reads. The library size of CP10k is the sum over the 2,000
genes, as in D1. Donor ids that name two animals (as `3_38_F/3_39_F` in Tabula Muris Senis) are
excluded.

## Case 1: MI coupling on Tabula Muris Senis (FACS)

Preprint §4's claim as published: the mutual information of Smad3 and Col1a1 is lower in old mice
than in young ones, within sex, the mouse the unit. Young: 3 months. Old: the Census term
"20-month-old stage and over" (§4's "20m": males of 24 months; females of 21 months, from the
mammary gland only, by the audit, `validation/flagship_audit/REPORT.md`).

- **Data:** the TMS FACS datasets (Smart-seq2, primary cells) as the audit resolves them.
- **Cell type, by rule, before any metric is computed.** A (type, sex) is eligible where at least
  3 young and 3 old mice of the sex have at least 20 cells of the type. Among the eligible types,
  the one with the most cells in which both Smad3 and Col1a1 are detected, over the eligible
  sexes' mice of the two ages. The claim runs within the sexes in which that type is eligible.
- **Claims:** (1a) `mi_3bin(Smad3, Col1a1)`, a decrease from young to old; (1b)
  `norm_pearson(Smad3, Col1a1)`, a decrease: v1's validated type of metric on the same pair.
- **Settings, as a user runs the command line:** the composition estimand; a SESOI of 0.1 for
  both (for `mi_3bin` the value of the documented example, `AGENTS.md`; the same for
  `norm_pearson`); GATE 4 the coupling probe at its default strength, 1.0, and delta_min at its
  default, 0.5 × SESOI; GATE 5 the controls of §4, Actb–Gapdh positive and Tst–Lrrc42 negative,
  with the claim's metric as `pair_metric`.

## Case 2: Xist and the Y genes on Tabula Muris Senis (FACS)

Female against male mice within the 3-month stratum (4 against 4 by the audit, composite ids
excluded).

- **Structure, checked first by code that does not import the engine** (R1's lesson). Per mouse:
  its cells, the share with Xist above 0, the share with any Y gene above 0, its annotated sex, its
  tissues. A mouse is female where at least 50% of its cells have Xist and at most 5% a Y gene,
  male where the reverse holds. A mouse that is neither, or whose call contradicts its annotation,
  is excluded and reported; the case runs if at least 3 mice of each sex remain.
- **Cell type, by rule:** the type with at least 20 cells in every remaining mouse of the stratum
  and the most cells; if there is none, the type with at least 20 cells in the most mice (at least 3
  of each sex), then the most cells.
- **Claims, R1's:** (2a) the mean log1p CP10k of Xist, a decrease from female to male; (2b) the
  same of the four Y genes summed (Ddx3y, Eif2s3y, Kdm5d, Uty), an increase. R1's settings: SESOI
  0.5, delta_min 0.25, GATE 4's module probe 2-fold in every cell, which R1 showed defective; the
  pilot shows what it does here. No GATE 5 controls.

## Case 3: a negative control, the same donors in two technical conditions

- **The search, by rule, over the Census's primary cells** (human and mouse), before any metric. A
  qualifying set is one collection, one tissue and one cell type, normal tissue (disease
  "normal"), and two conditions that differ in assay or in suspension type (cell or nucleus), in
  each of which at least 4 donors have at least 50 cells. Excluded: Tabula Muris Senis (its mice
  measured by FACS and by droplet are a candidate truth for v2, `v2_requirements.md` section 3), and
  Kang et al. 2018 and Levy et al. 2020 (the confirmatory set), by collection name or DOI. The set
  chosen has the most such donors, then the most cells. If no set qualifies, the case does not
  run, and the report says so.
- **Claims, 20, all two-sided** (no difference exists by construction): (3a) `norm_pearson` on the
  8 most correlated pairs among the 100 most expressed genes (mean CP10k over the pooled cells;
  pairs without a shared gene, taken greedily); (3b) the same among the genes ranked 501 to 600;
  (3c) the mean log1p CP10k of each of the 4 most expressed genes.
- **Settings:** the composition estimand; a SESOI of 0.1 for the pairs and 0.5 for the means; GATE
  4 the coupling probe at strength 1.0 (pairs) or the module probe 2-fold in every cell (means); no
  GATE 5 controls; the paired design, each donor in both conditions.

## The three analyses of every claim

- **(a) over cells:** the claim's metric on all cells of each group, within strata where the claim
  has them (the mean of the strata's differences); the group labels permuted across the cells 999
  times (within strata); a two-sided p.
- **(b) over donors, without gates:** the metric per donor (per donor and condition where paired);
  an exact permutation over donors within strata (sign flips where paired), Monte Carlo with 9,999
  relabelings above 20,000; absence by a TOST at the SESOI, the 90% t interval of the donor-level
  difference (nested or paired) inside ±SESOI.
- **(d) the engine:** `run_autopsy` with the settings above, seed 0, run log off; its verdict,
  cause and four fields, `raw_detected`, and its reading under the known issue.
- **The conclusions compared:** an effect in the claimed direction (either, for a two-sided claim),
  one against it, an absence, none (a refusal, INCONCLUSIVE, explained by depth, metric invalid,
  UNIDENTIFIABLE). Where (d)'s conclusion differs from (b)'s, the report names the gate by the
  engine's cause.

## What is reported

One log per stage (`real_pilot_select.log`, `real_pilot.log`) and the engine's reports in JSON:
per case and claim, (a), (b), (d) and its reading; where (d) changes (b) and by which gate; the
time of each analysis and of each job; every exception with its traceback; the tables of the
selection and of the structure. Nothing is tuned: the settings are those above, in one attempt.

## Not in the pilot

v0.4, v2, the claims of level B and the controls of the confirmatory set. No rule of rightness is
applied; the pilot informs it.
