# Plan of the new manuscript (v0.3): for the owner's approval

Status: a plan only. No text is written until the owner approves it, and `paper/manuscript.md`
(the v0.1.1 preprint) stays untouched and is not patched. Every number below comes from a
script in this repository or from a workflow run, and its source is given as a file and a
commit.

Short names for the sources:

| Name | What | Where |
|---|---|---|
| **tag** | the frozen engine and protocol | tag `v0.3.0-prereg` = 90a7718 |
| **v1-records** | the run's records before the key | branch `prereg/panel-v1` at 7782254, which is the run tag `panel-v1-run` |
| **v1-results** | the run's results | branch `results/panel-v1`: results 3cb3c21, scores adc7d65, anchors 445c0dd |
| **by-level** | rates by expression level | `validation/exploratory/v1_by_level.log` (d241ea4), from scores.json (sha256 253de86d…) |
| **audit** | the flagship audit | branch `results/flagship-audit` at d2a5420 (data_sha256 d921bc22…, script_sha256 644e6ef5…) |
| **step1** | the audit's report | `validation/flagship_audit/REPORT.md` (c7d073d) |
| **mice** | section E without composite ids | `validation/flagship_audit/mouse_level.log` (7feadc2) |

## 1. The problem

"Compute, then believe": a metric that differs between conditions is reported as biology.

- **Content.** The three failures that started the project (entropy anticorrelation, cardiac
  β, SMAD–ECM mutual information), described qualitatively.
- **No numbers from ref 1.** It is an internal checklist with no script behind it: v0.1.1's
  §1 numbers, such as 2.4× and 3670 vs 1540, cannot be sourced. The TMS case is told with the
  audit's numbers in §5 instead.

## 2. v0.1 and how it failed its own probes

The released v0.1.1 errs in both directions.

| Number | Source |
|---|---|
| A random-number metric gets "PASS — cleared 3 auto gates" (p01) | `validation/probes/README.md`, `baseline_v0.1.1.log` (tag) |
| Xist (female > male, demo data) dies at GATE 1 (2.00× stratum ratio) although GATE 2 retains 100% (p04) | same |
| 3 vs 3 mice with no age effect: 22/40 runs report "effect survives matching"; the replicate-level exact test 0/40 (p07) | same |
| 18 failures in all, found by the probes | `validation/probes/README.md` (tag) |

The probes are a development set: they found the bugs, so they carry no confirmatory weight.

## 3. How v0.3 works

- The four-field verdict: metric validity, design adequacy, effect, replication. One rule
  decides it (`report.decide`).
- GATE 0 separates bias from attenuation. GATE 4 tests the metric's response to an injected
  signal by an interval rule. GATE 5 uses empirical nulls.
- GATE 2 corrects depth by thinning, by the estimand. The effect is inferred over biological
  replicates.
- The development result, labelled as having no confirmatory weight
  (`validation/probes/verdicts_v0.3.0.dev0.log`, tag):
  - verdict outside the allowed set: 1/190;
  - false SUPPORTED: 1/150;
  - definite verdicts: 130/150.

## 4. The blind validation v1

**Design** (v1.md, tag):
- 9,750 datasets and 10,540 claim cards.
- The conditions: N1–N8 and E1–E3.
- The criteria S1–S7b, each by the principle: P(pass | sound) ≥ 0.90, P(pass | doubled error)
  ≤ 0.05.
- P(all pass | sound): 0.912 / 0.913 / 0.912 in `validation/prereg/oc.log` (tag); 0.915
  (SE 0.0009) in `oc_pilot.log` (v1-records).
- The key: drand quicknet round 32929913 (`key.json`, v1-results).
- The backgrounds (`selection.md`, prereg/panel-v1 at a0f482b), with suspension type and assay
  from the Census (`validation/exploratory/census_metadata.py` 729a9fd, workflow run
  38033204255):
  - B1: human oligodendrocytes of the MSSM cohort; nuclei, 10x 3' v3.
  - B2: mouse islet beta cells; cells, 10x 3' v2 and v3.

**Result.** Every criterion passes (`scores.txt` / `scores.json`, v1-results adc7d65):

| Criterion | Count | Allowed |
|---|---|---|
| S1 | N1 5, N2 6, N5 8, N6c 0, N8 11 (of 790 each) | ≤ 29 |
| S2 | 5, 11, 3, 0, 5 (of 790); 2 of 1,285 | ≤ 30, 31, 36, 29, 30; ≤ 51 |
| S3 | 181/189, 986/1,075, 4,674/4,685; 16 of 16 conditions | ≥ 144, 778, 3,332 |
| S4 | 0/515, 0/364, 58/2,165, 0/1,147, 1/7,020, 0/4,261 | ≤ 18, 26, 102, 97, 320, 393 |
| S5 | 4/206, 92/2,252 | ≤ 32, 439 |
| S6 | 0/121, 5/7,020 | ≤ 15, 660 |
| S7a | 0 of 10,540 rule violations | 0 |
| S7b | 0 of 10,540 crashes | ≤ 18 |

**The N8 prediction, and why it failed.** The pre-registered prediction is reported as "the
prediction was not confirmed" (v1.md §6). The reason is dilution by blind pairs. These
numbers are exploratory, not a criterion (by-level):
- All 2,458 cards with a valid metric are at the high level. At the medium and low levels the
  metric is blind (55 medium-level cards are ambiguous), and no card there was SUPPORTED. So
  253 of N8's 790 cards were at risk.
- N8 at the high level:
  - false SUPPORTED 11/253 = 4.3% (95% CI 2.2–7.6%), at a nominal 2.5%;
  - SUPPORTED or against the direction 22/253 = 8.7%, at a nominal 5%.
- N3 f = 0.4 at the high level: 7/83 = 8.4% and 10/83 = 12.0%.
- N1, N5 and N7 at the high level: 5.9%, 3.8% and 4.0%. The pure nulls together: a false
  SUPPORTED on 18/812 = 2.2%.
- At N8's high-level rate on all 790 cards, about 34 false SUPPORTED would be expected against
  29, so S1 would pass with probability ≈ 0.20.
- The dev probe p16 on simulated data gave 11.3% (`p16_variable_capture_n8.log`, tag).

**Limitations.**
- The scope: one metric, two backgrounds, the binomial-thinning family and N8. Nothing else is
  validated, `mi_3bin` included.
- v0.3 does not correct per-cell capture (v1.md §6's sentence).
- The design effects reach 5.62 (scores.json `limitations`).
- GATE 0 refused 26 of N7's 276 high-level cards (9.4%; by-level).

**Checks after the run.**
- `score.py` on the published reports reproduces every criterion (PROJECT.md, e765f24).
- `blind.py verify` on the first 20 datasets, on another CPU model: 20 of 20 datasets
  identical, and 22 of 22 reports the same (`validation/exploratory/v1_verify_first20.log`,
  7feadc2).

**The anchors** (secondary; `anchors.json`, v1-results 445c0dd):
- R2a and R2b are INCONCLUSIVE and R2c is UNIDENTIFIABLE, each in its allowed set.
- R3 was dropped with B4.
- R1 did not run in v1 (a loader bug). Its one run after the results, outside v1
  (DEVIATIONS.md D1; workflow run 38033821978; `results/panel-v1-r1` at 36a04bf), gave:
  - Xist: NOT SUPPORTED, metric invalid by GATE 4 (+0.163 against δ_min 0.25);
  - Y genes: UNIDENTIFIABLE (B2's donor ids are partially crossed with sex);
  - the sham: NO DETECTABLE EFFECT.

  Two of the three claims are outside their allowed sets. The section reports this as a real
  case where v0.3 does not deliver a known difference, with both causes.

## 5. The corrected TMS case

The audit pinned to v0.1.1 (audit, step1). Its seven edits are accepted, with two more: (a) and
(b) below.

- **v0.1.1's verdict on `mi_3bin` reproduces:**
  - FAIL at GATE 0: dropout shift 43.0%, z 54.7;
  - GATE 1: the male ratio is 1.646×;
  - GATE 2 retains 43% pooled and 30% for males, and the sign flips;
  - the negative control's MI is 0.036, and GATE 5 fails.
- **That verdict is uninformative.** On the same data v0.1.1 also kills Xist and the Y-gene score
  at GATE 1 (2.01× in the 20-month stratum), although GATE 2 retains 100% and 104%. The sex
  markers match the label in 21 of 21 donor units.
- **No spike-ins:** 0 of 53,384 Census features, and 0 in each of the 17 source files. A
  technical drop in capture cannot be told apart from a biological drop in RNA content.
- **The design:**
  - the "old (20 months)" group is 4 males aged 24 months (31,551 cells) and 3 females aged
    21 months (728 cells, all from the mammary gland);
  - §4's age "correction" of ref 1 is wrong.
- **(a) Mice without the composite ids** (mice):
  - males 4 vs 4: median nnz and `mi_3bin` both separate young from old completely, exact
    p = 2/70 = 0.029 for each;
  - females 4 vs 3: p = 9/35 (nnz) and 14/35 (MI).
- **(b) The female arm.** The old females come from the mammary gland only, so any sex × age
  comparison across all tissues confounds age with tissue composition.
  - **Proposed:** drop the claims "female detection does not decline" and "a sex-by-age
    interaction".
  - The alternative is a comparison within the shared tissue (mammary gland). That would need a
    new scripted run on the TMS data in a workflow, and it would rest on 3 old mice of 222–267
    cells each.

## 6. Discussion

- The gates are necessary, not sufficient. A verdict is a statement about one dataset under
  the declared assumptions, not about biology.
- The lesson of the dilution: count false-SUPPORTED criteria by level, or only on cards whose
  metric is valid or ambiguous.
- v0.4's targets: per-cell capture (N8) and entry-level dropout (N3 f = 0.4).
- Real development backgrounds: p16 put N8's rate 2.6 times high.
- GATE 0's refusals on N7 need explaining (PROJECT.md, 3b3e7cb).
- Step 4, external verdicts, is next.

## Figures (each from a script, written after approval)

1. The design of v1: conditions, cards, the key, the criteria.
2. The criteria against their limits (scores.json).
3. The error rates by expression level: N8, N3 f = 0.4 and the nulls (by-level).
4. TMS: per-mouse `mi_3bin` against median nnz, males separated on both (audit `E.mice`; mice).
5. The positive control v0.1.1 kills: Xist and Y genes at GATE 1 (audit `F`).
