# The fourth run of v0.4 on the development data: §4 again as a recorded second attempt, and Xist, the Y genes and SEA-AD again (the plan, committed before the run)

The owner's second decisions of 2026-10-18, items 2, 8 and 10 (`validation/plans/v0.4_spec.md`, "Decided
(2026-10-18, second)"). **Status: exploratory, on development data**, outside the confirmatory set of
the real track for good. One attempt; the outputs go to the branch `results/real-pilot-6`. Nothing here
is tuned after the run: every claim, criterion and score is fixed below, and the claims are the third
run's, called as the third run called them.

## What it is for

- **§4 again, A3 and A4, as a recorded second attempt** (item 2), with the engine after the second
  decisions: change 12 (GATE 0 blocks only a bias it cannot show equal between the groups), the domain
  of a pair "in each group", the full dossier after the first FAIL, `norm_pearson`'s NaN and the
  technical covariates per replicate. The first attempts are the third run's (`real_pilot_5.log`:
  both NOT SUPPORTED by GATE 0).
- **Xist, the Y genes and SEA-AD again** (item 10), to check that the changes did not break them.
- **What the TMS authors write** about the origin of the 3-month FACS data (item 8), before the
  confound with the source is called a finding.
- **No control for coupling.** The search of item 5 found no control that could be registered: two
  candidates, mosaic loss of chromosome Y and promoter–reporter lines, each still to be checked on its
  data by code independent of the engine (`validation/plans/v2_requirements.md`, section 3); neither
  runs here. Marsh et al. 2022 is no longer a control for coupling and is not read.

## The data

The third run's, unchanged (`real_pilot_diagnostics_3_plan.md`, "The data"), read by the same scripts:
- **cases 1 and 2**, the second run's files (`results/real-pilot-2` at
  bad785a3c56e9de1ee24737953c1b4b63112aa7c): `case1.npz` (preprint §4's setup, 76,797 cells of 17
  units, sha256 `13e7a563…`), `case2.npz` (cardiac fibroblasts of 3 female and 4 male 3-month mice,
  2,088 cells, `7191ad73…`) and `select.json`;
- **case 3** (SEA-AD), the first run's files (`results/real-pilot` at
  6133dfdf4c05f083f012aaec4136e376715af05a): `case3.npz` (15,711 nuclei, `31c297fe…`) and its
  `select.json` with the 20 claims;
- each checked by its run's SHA256SUMS (`real_pilot_5_inputs.sh`, without GEO's matrices);
- **the full set of genes** for the cells of cases 1, 2 and 3 from the CELLxGENE Census, release
  2025-11-08, by the third run's code (`real_pilot_6_extract.py` calls `real_pilot_5_extract.run_census`
  unchanged), with its check that the case files' genes have the same counts; case 1's metadata as the
  third run fetched it. Marsh et al.'s matrices are not fetched.
- **The engine**: v0.4.0.dev0 at the head of the pull request of `claude/v0.4` when the run starts,
  in v1's analysis environment (`validation/prereg/requirements-panel.txt`); the run records the
  commit and the git tree of `src/`. The scripts: `real_pilot_5_inputs.sh`, `real_pilot_6_extract.py`
  (v1's selection environment), `real_pilot_6.py` (the analyses; it loads and calls each claim with
  `real_pilot_5.load_case` and `real_pilot_5.accept_run`, the third run's functions), and
  `real_pilot_6_text.py` (the authors' text).

## The claims

`run_autopsy` of `claude/v0.4` on each claim as the third run called it (`real_pilot_5.accept_run`): the
composition estimand; the claim's direction; the SESOI on the construct's scale with the claim's probe;
the judgment gates resolved; seed 0; the run log off; the check of specificity on, its pool the full set
of genes of the same cells; no second dataset (so SUPPORTED is provisional).

| id | claim | data | as in the third run |
|---|---|---|---|
| A1 | `log_cp10k_mean(Xist)`, female → male, decrease, SESOI 0.5 | case 2 | probe `module(["Xist"])`; composition `cell_type`; `known_to_differ` the 23 sex genes |
| A2 | `log_cp10k_mean` of Ddx3y, Eif2s3y, Kdm5d, Uty summed, female → male, increase, SESOI 0.5 | case 2 | probe `module` of the four; as A1 |
| A3 | `mi_3bin(Smad3, Col1a1)`, young → old, decrease, SESOI 0.1; **the second attempt** | the males of case 1 | probe `coupling`; controls Actb–Gapdh and Tst–Lrrc42; composition `cell_type`, `tissue` |
| A4 | the same, within sex; **the second attempt** | case 1 | as A3 |
| A8 | the 20 claims of case 3 (3a1–3a8, 3b1–3b8, 3c1–3c4), two-sided, SESOI 0.1 (pairs) and 0.5 (means) | case 3 | probes as the first run's; paired by donor; composition `cell_type` |

Every difference from the third run's verdicts comes from the engine: the claims, data, settings and
seed are the same.

## The criteria (the owner's, 2026-10-18), and how each is scored

A verdict is SUPPORTED where its label is SUPPORTED (provisional or replicated).
- **A1 and A2: SUPPORTED.**
- **A3 and A4: not SUPPORTED, with the dossier**: the label is not SUPPORTED, and the report holds the
  effect's estimate (its field has a status), the domain, and the technical covariates' rows.
- **A8: no SUPPORTED; INCONCLUSIVE allowed** (the owner's corrected criterion): met where all 20 claims
  ran (no exception) and none is SUPPORTED. Any other label — NO DETECTABLE EFFECT, EFFECT BELOW SESOI,
  INCONCLUSIVE, UNIDENTIFIABLE, NOT SUPPORTED — meets it.

## The owner's prediction, fixed before the run

> GATE 0 по-прежнему блокирует, потому что у молодых и старых самцов различаются число генов и счётов на клетку (nnz +1298, p 0,016; raw_sum p 0,008).

Scored on A3 (the males) and A4 (its male stratum), in two parts:
1. **GATE 0 blocks**: GATE 0 FAILs, the verdict's cause is `metric_invalid_gate0`, and change 12's test
   did not show the groups' technical difference within GATE 0's tolerance (`between_groups`
   not equivalent), on both claims;
2. **the young and old males differ in genes and counts per cell**: the dossier's technical
   covariates, all cells, the males (A3's one stratum, A4's male stratum): genes per cell higher in the
   young (young − old > 0) and counts per cell lower (young − old < 0), each with its permutation p
   below 0.05, on both claims.

**Holds** where both parts hold, **partly holds** where one does, **does not hold** otherwise. The
sign of the counts is the third run's (`raw_sum` −188,999 young − old); the engine's covariates are the
replicates' medians over every gene of the pool, the Census's `nnz` and `raw_sum` over its measured
genes, so the two need not be equal to the unit.

## Reported, not scored

- For A3 and A4, the dossier: the effect within the domain (corrected and raw), the whole mixture's,
  outside the domain and per type, with the claim genes' counts per type; the types where the pair
  reaches the threshold in one group only; GATE 0's responses and change 12's test (the technical
  difference, its 90% interval, the margin); s_m; the specificity line; the technical covariates on
  all cells and on the domain's.
- **A3's domain** against the list read off the third run's table (a1) by the rule "in each group"
  (`real_pilot_6_spec.A3_DOMAIN_EXPECTED`): fibroblast, fibroblast of cardiac tissue,
  granulocytopoietic cell, mesenchymal stem cell, mesenchymal stem cell of adipose tissue. The
  detections do not depend on random draws, so a different domain would mean a different rule or data.
- A8's verdicts that differ from the third run's (18 NO DETECTABLE EFFECT, 3c3 EFFECT BELOW SESOI, 3a8
  INCONCLUSIVE), and the units each claim left out as unmeasured (`norm_pearson` gives NaN below 20
  cells with both genes per unit; v1's donor values fell below 20 in 1.1%, 12.5% and 40.5% of the
  high, medium and low levels, `v1_fallbacks.log`).
- A1's and A2's technical covariates (female against male).

## What the TMS authors write

`real_pilot_6_text.py`, in a job of its own: the paper (Tabula Muris Consortium 2020, Nature
583:590–595, doi:10.1038/s41586-020-2496-1; its NIH author manuscript PMC8240505, without a Creative
Commons licence) from Europe PMC's full text and from NCBI's BioC, and the bioRxiv preprint
(doi:10.1101/661728) from its newest JATS text. Each source's URL, outcome, size and SHA256 are logged,
and every sentence matching one of `real_pilot_6_spec.TEXT_PATTERNS` (the 3-month stage, "previously
published", Tabula Muris without "Senis", batches, time points) is quoted with its section, at most 40
per source, each cut at 600 characters; the texts themselves are not written anywhere. **The rule:** the
authors are taken to state that the 3-month FACS data are the original Tabula Muris's where a quoted
sentence names the 3-month data (or time point) together with the Tabula Muris or with "previously
published"; to state otherwise where a quoted sentence gives another origin; and neither where no
sentence does, or no source could be fetched. The quotes are read after the run and reported as quoted.

## Known before the run

From the earlier runs' logs, so not outcomes of this one:
- **A3 and A4's first attempts** (`real_pilot_5.log`, `real_pilot_5_summary.log`): NOT SUPPORTED by
  GATE 0; extra dropout moved the null by −0.002584 (z 7.64) and −0.002458 (z 5.32) against the
  tolerances 0.0007749 and 0.0007367 (0.5 × s_m, s_m +0.00155 and +0.00147); the domain by the rule then
  ("in at least one group") 11 types in males; the effect within the domain not estimated.
- **The technical difference of the young and old males** (the third run's (b4)): the Census's `nnz`
  +1,298 (exact p 0.0159), `raw_sum` −188,999 (p 0.0079); in females neither (p 0.64 and 0.36).
- **A1 and A2** SUPPORTED (provisional), **A8** 18 NO DETECTABLE EFFECT, one EFFECT BELOW SESOI (3c3)
  and one INCONCLUSIVE (3a8).

## One attempt

`.github/workflows/real-pilot-6.yml`, started by the label `pilot-run-6` on the pull request of
`claude/v0.4`, commits whatever the analyses produce, failures included, to the new branch
`results/real-pilot-6`; a second start finds the branch and stops. If the preparation or the extraction
fails, no analysis starts and nothing is committed: it may be fixed and started again, and the report
says so. The analyses run in six parallel jobs (A1–A2; A3; A4; A8 in three parts), the authors' text in
a seventh; each job has one attempt. The third run's parts took 5 s to 21 min (`real_pilot_5.log`); the
dossier adds the effect's estimate after a FAIL, and change 12 an estimate of the metric on GATE 0's
null.

## What is reported

`real_pilot_6.log` and `real_pilot_6.json`: each claim's verdict, cause, gates and dossier; each
criterion's score and the prediction's; A3's domain against the expected list; A8's changes against the
third run; and `real_pilot_6_text.log` with the authors' quotes. The hashes of the inputs and the
engine's commit. One report to the owner.

## The run

Workflow run 39163054704 of `real-pilot-6.yml`, started by the label `pilot-run-6` at 20:23 UTC on
2026-10-18 on the pull request's head 57980d5 (the engine's `src/` git tree
d7f5fef250dd440c8110d09d2baa78f86917b59d), one attempt, every job green (the text 56 s, the extraction
639 s, the six parts 32 s to 20 min); its outputs on `results/real-pilot-6` at 1e2f824, every file as its
SHA256SUMS. Copied here: the log as `real_pilot_6.log`, the extraction's as `real_pilot_6_extract.log`,
the inputs' as `real_pilot_6_inputs.log`, the authors' quotes as `real_pilot_6_text.log`. The extraction
matched the case files' counts (62 of 62 genes checked in case 1, 55 of 55 in case 2, 50 of 50 in case 3).
**Accepted by the owner (2026-10-19): 5 of 5 criteria, the prediction scored as "partly holds"; no
further run on the development cases (A1–A8) before the tag of v0.4.**

| test | criterion | verdict | met |
|---|---|---|---|
| A1, Xist | SUPPORTED | SUPPORTED (provisional; parametric only) | yes |
| A2, the Y genes | SUPPORTED | SUPPORTED (provisional; parametric only) | yes |
| A3, §4 in males (attempt 2) | not SUPPORTED, with the dossier | INCONCLUSIVE (`effect_inconclusive`) | yes |
| A4, §4 within sex (attempt 2) | not SUPPORTED, with the dossier | INCONCLUSIVE (`effect_inconclusive`) | yes |
| A8, SEA-AD | no SUPPORTED; INCONCLUSIVE allowed | 18 NO DETECTABLE EFFECT, 1 EFFECT BELOW SESOI (3c3), 1 INCONCLUSIVE (3a8) | yes |

- **The owner's prediction: partly holds.** Its first part does not: GATE 0 passes on both claims,
  so change 12 was not reached. On the domain by the rule "in each group" (five types, as expected)
  extra dropout is an attenuation, −64% of the signal (−63% in A4), its null moving by −0.00117 (z 2.09,
  below 4) against a tolerance of 0.00136 (0.5 × s_m, s_m +0.00273); in the third run, on the
  eleven types of the rule "in at least one group", it moved the null by −0.00258 (z 7.6). The second
  part holds: the young male units have 1,298 more genes and 189,000 fewer counts per cell (exact p
  0.0159 and 0.0079 over 126 relabelings), on the domain's cells 1,166 more genes (p 0.0159).
- **A3 within the domain**: −0.0018 [−0.0179, +0.0142] (exact p 0.771 over 70 relabelings), the MDE
  0.0219; the whole mixture +0.0198 (p 0.0159), outside the domain +0.0172 (p 0.0286). What these
  show is written up in `real_pilot_6_summary.md`, computed from the run's JSON by
  `real_pilot_6_summary.py`.
- **A8**: every verdict as in the third run.
- **What the TMS authors write** (2 of 3 sources: Europe PMC answered HTTP 500; NCBI's BioC and the
  preprint's JATS text were read): by the rule above the authors state that the 3-month data are the
  original Tabula Muris's — "The previously published 3m time point, referred to as the Tabula Muris,
  ... was used as a basis to perform semi-automated cell type annotation of the additional time
  points", and "We used the Tabula Muris (3m time point) as a reference for the automated pipeline";
  no quoted sentence gives another origin.
