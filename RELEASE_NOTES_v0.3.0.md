# metric-autopsy v0.3.0 — validation v1 record

**What this release is.** The archive of validation v1: the engine frozen at the tag
`v0.3.0-prereg`, the protocol it was tested against, the blind run with its scores and anchors,
and the analyses made after the results. It carries no new engine. `src/` is byte-identical to
the tag (git tree `97cbab3c437c9ddb118a0690e6e0d24f8a16069f`), and the engine reports its version
as `0.3.0.dev0`.

**Status. v0.3.0-prereg passed the pre-registered blind validation v1, within its scope only.
The pass excludes only gross excesses of false SUPPORTED on the gene pairs where the metric
works. v0.3 confirmed no real biological difference.**

## The records

| What | Where |
|---|---|
| The frozen engine and protocol | tag [`v0.3.0-prereg`](https://github.com/mool32/metric-autopsy/tree/v0.3.0-prereg), commit `90a7718f77a46440ca11cd90150733a67f6b138d`; the protocol [`validation/prereg/v1.md`](https://github.com/mool32/metric-autopsy/blob/v0.3.0-prereg/validation/prereg/v1.md), approved by the owner on 2026-10-09 |
| The run | tag [`panel-v1-run`](https://github.com/mool32/metric-autopsy/tree/panel-v1-run), commit `7782254bc0c1edf1a1cd26f1ceba03054a0631c0` (branch `prereg/panel-v1`); [workflow run 37990113304](https://github.com/mool32/metric-autopsy/actions/runs/37990113304), one attempt, no shard re-run; the key is drand quicknet round 32929913 |
| Results, scores, anchors | branch [`results/panel-v1`](https://github.com/mool32/metric-autopsy/tree/results/panel-v1): the results at `3cb3c211a6b933254a89f62ffcf33acd3c4052ab`, the scores (`scores.json`, `scores.txt`) at `adc7d6589a4448b40262fc22d330745fb7a141c8`, the anchors at `445c0dd1e3f2c0909afbd7c71cbb60d2b23d1b42` |
| Anchor R1, once, after the results and outside v1 | deviation [D1](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/prereg/DEVIATIONS.md); [workflow run 38033821978](https://github.com/mool32/metric-autopsy/actions/runs/38033821978); branch [`results/panel-v1-r1`](https://github.com/mool32/metric-autopsy/tree/results/panel-v1-r1) at `36a04bfb85e5e1b3acb58b428fdeca4bec46d853`; run commit `b0ce92de1fb333b284bd93ed2e002c9749b6ce75`, the tag plus D1's files |
| The audit of v0.1.1's flagship example | [`validation/flagship_audit/REPORT.md`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/flagship_audit/REPORT.md); [workflow run 37990113298](https://github.com/mool32/metric-autopsy/actions/runs/37990113298); branch [`results/flagship-audit`](https://github.com/mool32/metric-autopsy/tree/results/flagship-audit) at `d2a54205cd7a0270472d8a594913622977fead6b` |
| Analyses after the results | [`validation/exploratory/`](https://github.com/mool32/metric-autopsy/tree/v0.3.0/validation/exploratory), each script with its log |
| The development set | [`validation/probes/`](https://github.com/mool32/metric-autopsy/tree/v0.3.0/validation/probes), the probes that found the bugs; no confirmatory weight |

## Scope

- **Data:** scRNA-seq counts only.
- **Metric:** one, the log-normalized Pearson correlation of a gene pair (`norm_pearson`), on the
  8 most correlated pairs of each expression level; and three useless metrics: random numbers
  (N6a), a constant (N6b) and a score of random genes (N6c).
- **Backgrounds**, from the CELLxGENE Census
  ([`census_metadata.log`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/exploratory/census_metadata.log)):
  - B1: human oligodendrocytes of the MSSM cohort; nuclei, 10x 3' v3. It carried every
    condition but N7: 9,750 of the 10,540 cards.
  - B2: mouse pancreatic islet beta cells; cells, 10x 3' v2 (39 donor ids) and v3 (11). Its 50
    donor ids are samples, 6 of them pools of both sexes; all single females are the NOD strain.
    It carried only the null N7 and the anchor R1.
- **Artifacts:** N2, N3 and N8: capture loss and entry-level dropout, of v1's binomial-thinning
  family, and per-cell capture, outside the engine's correction family. N2 did not test the
  correction of depth: the oracle's power to detect its raw difference was 0.21, 0.14 and 0.16 at
  the three levels, never the 0.9 that makes it informative
  ([`v1_pilot_table.log`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/exploratory/v1_pilot_table.log)).
- **Real effects:** no level had a key dose, so E1–E3 ran at the saturation dose (B1: 2.5, 2.75
  and 3.0). Of their 1,800 cards 206 had a valid metric, 189 of them establishable and 121 at or
  above the SESOI; all 1,285 with a blind metric were establishable too, with "metric invalid"
  as the correct outcome. At the high level Δ\* (0.258) is about the SESOI (0.255).
- **Not tested:** GATE 6 (no card had a second dataset, so every SUPPORTED was provisional) and
  the content estimand (only in the anchor R2).

Everything else is not validated: other metrics (`mi_3bin` included), other data and other
artifacts.

## The result

9,750 datasets, 10,540 claim cards. The criteria's numbers come from `scores.json` on
`results/panel-v1`.

| Criterion | What it counts | Result | Limit |
|---|---|---|---|
| S1 | false SUPPORTED on the key nulls, 790 datasets each | N1 5, N2 6, N5 8, N6c 0, N8 11 | ≤ 29 in each |
| S2 | false SUPPORTED in the other groups | N2's other steps 5, N3 11, N4 3, N6a 0, N7 5 (of 790 each); real effects 2 of 1,285 | ≤ 30, 31, 36, 29, 30; ≤ 51 |
| S3 | correct definite outcomes on establishable cards | real effects 181/189, nulls 986/1,075, invalid metrics 4,674/4,685; all 16 per-condition cells | ≥ 144, ≥ 778, ≥ 3,332 |
| S4 | NOT SUPPORTED against the direction / explained by depth where false; on nulls "against the direction" is allowed, but bounded as a false detection | real effects 0/515, 0/364; nulls 58/2,165 (allowed), 0/1,147; invalid metrics 1/7,020, 0/4,261 | ≤ 18, 26; 102, 97; 320, 393 |
| S5 | false "metric invalid" on a valid metric | real effects 4/206, nulls 92/2,252 | ≤ 32, 439 |
| S6 | false NO DETECTABLE EFFECT | real effects 0/121, invalid metrics 5/7,020 | ≤ 15, 660 |
| S7a | verdicts the engine's rules cannot give | 0 of 10,540 | 0 |
| S7b | crashes | 0 of 10,540 | ≤ 18 |

- **No useless metric was SUPPORTED:** 0 of 1,630 cards.
- **The prediction on N8 was not confirmed.** The protocol predicted that v0.3 fails S1 on N8
  (per-cell variable capture): 11 false SUPPORTED of 790, where S1 allows 29. Its S4 part was not
  confirmed either: S4's null cell held 58 of 2,165 (limit 102), 11 of them N8's.

## What the pass excludes

Computed after the results ([`v1_null_cells.log`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/exploratory/v1_null_cells.log)).
A false SUPPORTED needs a metric that responds, and in each key null of the pair metric (N1,
N2, N5, N8) only 253–271 of the 790 cards had a valid or ambiguous metric, all at the high
level. On those cards S1's cell fails with probability 0.95 only at a false-SUPPORTED rate of
14.3–15.3%; a validator wrong on 5% of them passes with probability 1.00, one wrong on 7.5%
with 0.98–0.99. S2's cells hold 242–331 cards at risk and fail with probability 0.95 at
12.1–16.9%. The pilot showed this mix before the key (`oc_pilot.log`: a valid metric on 2,453 of
the 8,910 expected pair-metric cards). At the rates the pilot measured, a sound validator, with
the oracle's true correction on N3 and N8, passes the criteria together with probability 1.000:
a property of the criteria, not a prediction for the engine.

## Known issue in v0.3.0

**A SUPPORTED whose raw replicate-level test did not detect the effect was made by the depth
correction: read it as INCONCLUSIVE.** The raw test, on the same data without the correction, is
in the JSON report as `fields.effect.detail.raw_detected` (with `raw_effect` and `p_raw`). The
correction thins the deeper group to equal depth, which is exact for a difference of capture;
where the groups differ otherwise, as by entry-level dropout, it can enlarge an undetected raw
difference until it is detected. On v1's cards (exploratory, after the results;
[`v1_correction_rule.log`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/exploratory/v1_correction_rule.log)): 18 of the 54 false
SUPPORTED on null data had `raw_detected` false, all on cards at risk, among them all 7 of N3 at
f = 0.4 on the high level. On 15 of the 18 the correction had enlarged the raw difference, on 1
shrunk it and on 2 (both N2) reversed its sign: v0.3 makes a sign reversed by the correction
INCONCLUSIVE only where the raw difference was detected. Of the 275 SUPPORTED on real effects, 3
had `raw_detected` false (none of the 118 valid cards with |Δ\*| ≥ SESOI). v0.4 is to make this
the rule: a correction can only remove a conclusion, never create one.

## Limitations

- **The SESOI and GATE 4's dose came from the pilot**, that is, from the oracle, on every card.
  A user sets both without a pilot.
- **The truth about validity is circular.** The oracle calls a metric valid when it responds to
  the injection that GATE 4 itself uses, so S5 measures GATE 4's sampling error and could not
  find a probe built the wrong way; R1's was one.
- **The metric worked only at the high expression level.** All 2,458 cards with a valid metric
  are there; at the medium and low levels the metric is blind (55 medium-level cards of N7 are
  ambiguous), and no card there was SUPPORTED (0 of 6,988)
  ([`v1_by_level.log`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/exploratory/v1_by_level.log); pre-registered secondary
  outcomes, v1.md sections 5 and 7).
- **v0.3 does not correct capture that varies from cell to cell, nor entry-level dropout.** Its
  depth correction brings whole groups to one depth. In the twelve high-level null cells, after
  Holm's correction over the twelve (decided after the results) no cell differs from the nominal
  2.5%; without it, N3 at f = 0.4's interval excludes it. The largest are N3 at f = 0.4, 7/83 =
  8.4% (95% CI 3.5–16.6%; Holm 0.058), and N8, 11/253 = 4.3% (2.2–7.6%; Holm 0.61). N8's is a
  hypothesis for v0.4 with a mechanism, not a measured defect. N3's was found after the results
  to be made by the depth correction (the known issue above): all 7 had raw differences the
  replicate-level test did not detect (p 0.06 to 0.34), which the correction enlarged 1.36 to 2.69
  times, and a test over donors without the correction gives none of the 83
  ([`v1_analyzers.log`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/exploratory/v1_analyzers.log), section 5).
- **Design effects** (datasets share donors) exceed 1.5 only for S3's null and invalid strata
  (4.74, 5.62), S5's null cell (5.28) and S3's E1 condition (1.54). In the false-SUPPORTED cells
  they are at most 1.055; at N5's, S1's cell passes a doubled error with probability 0.051, just
  above the principle's 0.05. Cards that share a gene pair, which v1 did not compute, give larger
  ones (exploratory, [`v1_pair_deff.log`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/exploratory/v1_pair_deff.log)): 3.91 for
  S1's N8 cell, 6 of whose 11 false SUPPORTED fell on one pair, and at that design effect the cell
  passes a doubled error with probability 0.211; 29.8 for S5's null cell; 5.93 for S3's null
  stratum.
- **v1's sensitivity rests on 3 gene pairs.** SUPPORTED on 118 of the 121 valid real-effect cards
  with |Δ\*| ≥ SESOI (97.5%, 95% CI 92.9–99.5%), and those 121 cards share 3 pairs.
- **GATE 0 refused 31 cards, all on N7** (26/276 = 9.4% at the high level, 5 at the medium).
  A refusal is allowed, never definite, but on these null data the refusals cost decisiveness
  (S3's N7 cell: 204 of 237). Decided after the results: 29 of the 31 fell on three of B2's gene
  pairs; eight tests found no association with B2's chemistry or depth (one at p = 0.049)
  ([`b2_depth_by_chemistry.log`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/exploratory/b2_depth_by_chemistry.log)).
- **The development data erred in both directions.** The development probe p16 (11.3% on
  simulated data) was 8.1 times N8's real rate over all levels, and at most 3.7 times on the
  datasets at risk; it ran with a SESOI of 0.1 and a dose of 2.0, where B1's high level had 0.255
  and 2.5, so simulation and settings are confounded.

## The failure of anchor R1

R1 (Xist and the Y genes, female against male mouse islets in B2) was the one real positive
control expected to succeed.

- **It did not run in v1.** `anchors.py` passed B2's sparse matrix to `SimpleData`, which accepts
  only dense arrays, so the engine was never reached.
- **It ran once, after the results and outside v1** (D1: the loader, and with it the cut to 2,000
  genes, which changed the denominator of CP10k). **Two of its three claims are outside their
  allowed sets.** The causes, read after the outcome
  ([`b2_sex_structure.log`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/exploratory/b2_sex_structure.log)):
  - **The Y genes (and the design of the Xist claim): UNIDENTIFIABLE, a correct reading of B2.**
    Its 50 donor ids are 34 single males, 10 single females (all NOD) and 6 pools of both sexes,
    and within a development stage the sexes meet only inside the pools, so sex is mixed with
    both stage and strain. The protocol erred: its allowed set counted the ids of each sex over
    all stages, and checked neither that an id is one animal of one sex nor that both sexes
    share a stratum.
  - **Xist: "metric invalid" by GATE 4, which is false: an engine defect** (next section).
  - **The sham, female against female: NO DETECTABLE EFFECT, which is correct.**
- v0.1.1 had not confirmed Xist either, because of GATE 1, on Tabula Muris Senis. In v0.3 GATE 1
  is a diagnostic and does not block by design, but v0.3 was not run on that data; on B2 another
  gate, GATE 4, failed the control, and B2's design could not have confirmed it either.

## Defects of GATE 4's module probe

- **It measures the response over all cells.** For a metric that is a mean over cells, that
  response is, by identity, the share of cells with the marker times the response in them, which
  is below ln 2 for a well-detected marker. R1's probe raised Xist 2-fold in every cell; with
  Xist's response of 0.631 in B2's Xist cells it needed at least 0.25 / 0.631 = 39.6% of cells
  with Xist, and B2 has 25.8%. Measured: +0.163 over all cells, below delta_min 0.25. With
  `anchors.py`'s constants the failure was predictable from B2's composition
  ([`gate4_module_analytic.log`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/exploratory/gate4_module_analytic.log)).
- **Its default cannot pass a well-detected mean-expression metric at a SESOI of 0.5:** 2-fold
  in 30% of cells gives at most 0.3 × ln 2 = 0.208 against delta_min 0.25, by arithmetic
  (confirmed by simulation in
  [`gate4_module_dilution.log`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/exploratory/gate4_module_dilution.log)).
- **On shallow data it can pass for the wrong reason.** At a depth of 2,000 counts it passes a
  marker of single counts, but of a response of 0.360 the sham's loss of the single counts makes
  0.178 and the 2-fold change 0.182, since the sham thins the marker and the injection does not.
- The module probe is reached through the Python API only; the CLI and the MCP server offer the
  coupling probe. v1's panel used the module probe only on the useless metric N6c.

## Checked after the run, on another CPU model

- `score.py` on the committed reports reproduces every criterion and number but one secondary
  value, `secondary.s1_at_realized_deff`: for a design effect even 2e-16 above 1 it switches to
  a beta-binomial, where scipy overflows (S1's N2, design effect 1 + 1.6e-15; re-scored on the
  other CPU, N8 too). At a design effect of 1 both are 0.982 / 0.046, N1's.
- `blind.py verify` re-ran the first 20 datasets
  ([`v1_verify_first20.log`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/exploratory/v1_verify_first20.log)): all 20 datasets
  identical, and all 22 reports the same (16 byte-identical, 6 within 1e-9, with the same
  verdicts and causes).

## What the gates add over a donor-level test (exploratory)

After the results, four analyzers read every card of v1 on v1's own datasets, which are
development data now (the plan was committed before the run:
[`v1_analyzers_plan.md`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/exploratory/v1_analyzers_plan.md); results:
[`v1_analyzers.log`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/exploratory/v1_analyzers.log); all 9,750 datasets rebuilt from the
key and identical to the run's): (a) a permutation test over cells, the common analysis; (b) a
permutation test over donors, without gates; (c) the engine with GATE 0, 2, 4 and 5 off, its own
code (it reproduces the engine's raw replicate-level test on 2,091 of 2,091 cards with an exact
permutation or a t interval); (d) the full engine, the published reports. "At risk": a valid or ambiguous metric.

| Share "effect" (SUPPORTED for (c), (d)) | (a) | (b) | (c) | (d) |
|---|---|---|---|---|
| pure nulls N1, N5, N7, metric at risk (867 cards) | 18.6% | 2.9% | 2.8% | 2.1% |
| pure nulls, blind metric (1,503) | 8.6% | 2.5% | 2.5% | 0% |
| artifacts N2, N3, N8, metric at risk (1,018) | 16.3% | 4.8% | 4.8% | 3.2% |
| artifacts, blind metric (2,142) | 17.6% | 8.2% | 8.2% | 0% |
| N4, 3 against 3 donors (1,580) | 12.4% | 0% | 0.9% | 0.2% |
| useless metrics N6a–c (1,630) | 18.2% | 2.7% | 1.8% | 0% |
| real effects, valid metric, \|Δ\*\| ≥ SESOI (121) | 100% | 99.2% | 99.2% | 97.5% |
| real effects, ambiguous metric (309) | 72.5% | 57.0% | 57.0% | 28.5% |
| real effects, blind metric (1,285) | 21.3% | 10.7% | 10.7% | 0.2% |

- **Where the metric works, a test over donors without gates is almost as good.** The large step
  is from cells to donors.
- **On all 7,110 null cards** the false "effect" is 4.0% for (b) and 0.8% for (d), and 212 of the
  232 fewer false effects are on metrics blind to the planted coupling, where the engine refuses
  to interpret ([`v1_correction_rule.log`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/exploratory/v1_correction_rule.log)).
- **The gain on depth artifacts over all levels is mostly GATE 4's refusal** at the medium and low
  levels, where every card had a blind metric (artifacts at the medium level 9.9% for (b) against
  0% for (d), at the low 6.5% against 0%). At risk the correction works both ways: N2 5.4% → 2.1%,
  N8 7.1% → 4.3% (on the same pairs as (b)'s), N3 1.2% → 4.5%.
- **The cost:** on an ambiguous metric the real effects are found half as often (57.0% against
  28.5%). On the real effects with a blind metric, (b) is not significant on 87.6%, (c) claims NO
  DETECTABLE EFFECT on 44.7% and (d) calls the metric invalid on 99.8%.
- **One mechanism builds the real effects, the truth about validity and GATE 4's probe.** That a
  refusal on a blind metric is right, and (c)'s "no detectable effect" there wrong, holds for
  metrics blind to that injection; whether they are blind to real biology only real data can
  tell.

## Next

v0.4 (a draft specification, [`v0.4_spec.md`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/plans/v0.4_spec.md)): GATE 4's
response measured with the claim's own replicate-level test, the signal injected in the group where
the probed genes have more counts per cell, whatever the claim's wording, with a module probe that
creates no response by losing counts; a rule for the SESOI and GATE 4's dose without an oracle; and
the rule of the known issue above, by which a correction can only remove a conclusion, never
create one. A correction for per-cell capture is postponed to v0.5; v0.4 only reports the spread
of per-cell depth. Then the test that decides whether the gates are needed at all, on published
claims with an independent verdict recorded before the runs ([`real_claims.md`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/validation/plans/real_claims.md)), and a new
pre-registration, validation v2, whose requirements are collected in
[`PROJECT.md`](https://github.com/mool32/metric-autopsy/blob/v0.3.0/PROJECT.md) ("Requirements for v2"). The draft manuscript on v0.3,
`paper/manuscript_v03.md`, is not for posting; it will be rewritten after v0.4, validation v2 and
the external verdicts.

## How to cite

Cite the version DOI Zenodo assigns to this release; the concept DOI
[10.5281/zenodo.21195679](https://doi.org/10.5281/zenodo.21195679) resolves to the latest
version. Full log: `CHANGELOG.md`.
