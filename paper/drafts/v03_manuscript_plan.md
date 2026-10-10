# Plan of the new manuscript (v0.3)

Status: the owner approved the plan on 2026-10-10 with edits, which this version carries out.
The text is written in `paper/manuscript_v03.md`: the abstract and §4 first, the rest after the
owner's reply. `paper/manuscript.md` (the v0.1.1 preprint) stays untouched and is not patched.
Every number below comes from a script in this repository or from a workflow run, and its
source is given as a file and a commit.

## Thesis

A check meant to tell biology from artifact in a single-cell metric must itself be validated
blind, on real data and with real positive controls: ours passed a pre-registered test against
planted truth within a narrow scope, yet missed its only real positive control, and the data it
was developed on had erred in both directions.

## Abstract skeleton

1. A metric that differs between single-cell conditions is often read as biology, although
   depth, dropout and the design can make the same difference.
2. metric-autopsy runs a pre-registered claim through gates aimed at named artifacts and decides
   one verdict from four fields: metric validity, design adequacy, an effect over biological
   replicates at equal depth, and replication.
3. Its first release, v0.1.1, erred in both directions on its own probes: it passed a random
   metric and blocked a known sex difference. We rebuilt the logic (v0.3) and froze it before a
   blind test.
4. In a pre-registered blind validation on two real CELLxGENE backgrounds with planted truth
   (9,750 datasets, 10,540 claim cards), v0.3 met every error criterion within its scope. The
   prediction that it would fail on per-cell variable capture was not confirmed; at high
   expression its false-SUPPORTED rate there was 4.3% at a nominal 2.5% (exploratory).
5. It did not deliver its only real positive control, the sex difference in Xist and the Y genes
   in mouse islets. The background let the sexes be compared only inside pooled samples, which
   the protocol had missed, and GATE 4 measured an injected Xist signal over all cells, three
   quarters of which lack Xist, and so called a valid metric invalid.
6. The development data had erred in both directions: they overstated one error rate 2.6-fold
   and hid the GATE 4 defect. So a validator has to be tested blind, on real backgrounds, with
   real positive controls and with the settings a user would choose. We report the scope, the
   failures and the targets of v0.4.

## Title options

1. Validating the validator: a blind, pre-registered test of artifact checks for single-cell
   metrics
2. metric-autopsy: artifact checks for single-cell metrics, their blind validation, and where
   they fail
3. Before believing a single-cell metric: a gate system and a pre-registered test of its own
   errors

## Sources

| Name | What | Where |
|---|---|---|
| **tag** | the frozen engine and protocol | tag `v0.3.0-prereg` = 90a7718 |
| **v1-records** | the run's records before the key | branch `prereg/panel-v1` at 7782254, which is the run tag `panel-v1-run` |
| **v1-results** | the run's results | branch `results/panel-v1`: results 3cb3c21, scores adc7d65, anchors 445c0dd |
| **by-level** | rates by expression level | `validation/exploratory/v1_by_level.log` (d241ea4), from scores.json (sha256 253de86d…) |
| **r1** | anchor R1, run once after the results | branch `results/panel-v1-r1` at 36a04bf (workflow run 38033821978; run commit b0ce92d); `validation/prereg/DEVIATIONS.md` D1 |
| **b2** | B2 by sex and GATE 4's Xist probe, after R1 | `validation/exploratory/b2_sex_structure.log` (12e2ef5, workflow run 38045140802) |
| **atlas** | B2's source: sample sex and how cells' sex was assigned | theislab/mouse_cross-condition_pancreatic_islet_atlas at 3af65e4 (Hrovatin et al. 2023) |
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
- The conditions: N1–N8 and E1–E3. N7 permutes the group labels over B2's donor ids: a null at
  the level of the id, not of the mouse.
- The criteria S1–S7b, each by the principle: P(pass | sound) ≥ 0.90, P(pass | doubled error)
  ≤ 0.05.
- P(all pass | sound): 0.912 / 0.913 / 0.912 in `validation/prereg/oc.log` (tag); 0.915
  (SE 0.0009) in `oc_pilot.log` (v1-records).
- The key: drand quicknet round 32929913 (`key.json`, v1-results).
- The backgrounds (`selection.md`, prereg/panel-v1 at a0f482b), with suspension type and assay
  from the Census (`validation/exploratory/census_metadata.py` 729a9fd, workflow run
  38033204255):
  - B1: human oligodendrocytes of the MSSM cohort; nuclei, 10x 3' v3.
  - B2: mouse islet beta cells; cells, 10x 3' v2 and v3; 50 donor ids (samples), 6 of them
    pools of both sexes; all single females are the NOD strain (b2).

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
- GATE 4's module probe measures the response over all cells, so for a marker present in only
  some cells it can call a valid metric invalid (§4.7). v1's panel used the module probe only on
  the useless metric N6c, so v1 did not test this path.

**Checks after the run.**
- `score.py` on the published reports reproduces every criterion (PROJECT.md, e765f24).
- `blind.py verify` on the first 20 datasets, on another CPU model: 20 of 20 datasets
  identical, and 22 of 22 reports the same (`validation/exploratory/v1_verify_first20.log`,
  7feadc2).

**The anchors R2 and R3** (secondary; `anchors.json`, v1-results 445c0dd):
- R2a and R2b are INCONCLUSIVE and R2c is UNIDENTIFIABLE, each in its allowed set. R2's known
  difference (the G2M score, G2M against G1) could not be confirmed by design: each phase is
  one capture batch.
- R3 was dropped with B4.
- R1 did not run in v1 (a loader bug); it is §4.7.

### 4.7 The only real positive control

R1 (DEVIATIONS.md D1; r1): Xist from female to male mice, a decrease, and the Y-gene score
(Ddx3y, Eif2s3y, Kdm5d, Uty), an increase, within the development stage, the donor id as the
mouse, SESOI 0.5 on the log1p-CP10k scale (δ_min 0.25); and a sham, female against female. It
ran once, after the results and outside v1.

**B2's structure** (b2; atlas):
- 50 donor ids, which are samples. By the share of their cells with Xist > 0:
  - 34 single males, Xist in 0.0–1.5% of their cells and a Y gene in 68.0–94.2% (Fltp_adult 4,
    STZ 7, VSG 8, spikein_drug 15);
  - 10 single females, 98.0–100.0% and 0.0–0.7% (NOD 3, NOD_elimination 7), all of the NOD
    strain;
  - 6 pools of both sexes, 53.0–80.0% and 33.8–45.8% (Fltp_P16 at 2 weeks, Fltp_2y at 20 months
    and over), each with 227–267 cells annotated female and 133–173 male of its 400.
- The atlas confirms both: its sample table lists the six as "mixed" and all twelve of its NOD
  samples as female (B2 holds ten). In the pools each cell's sex was assigned from a score of the
  Y-chromosome genes (the "data-driven" sex annotation; `2_annotate_Fltp_P16.py`,
  `2_annotate_Fltp_2y.py`).
- No development stage holds a single male and a single female. Within a stage the sexes meet
  only inside the pools (table from b2).

**The outcome** (r1): two of the three claims are outside their allowed sets.

| Claim | Verdict | Allowed |
|---|---|---|
| Xist | NOT SUPPORTED: metric invalid (GATE 4, +0.1627 against δ_min 0.25) | SUPPORTED |
| Y genes | UNIDENTIFIABLE: donor ids partially crossed with sex | SUPPORTED |
| Sham | NO DETECTABLE EFFECT | NO DETECTABLE EFFECT, INCONCLUSIVE |

**The causes differ.**
- **Y genes (and the design of the Xist claim): UNIDENTIFIABLE, a correct reading of B2.** The
  protocol erred: `r1_allowed` counted the ids of each sex over all stages (16 female, 40 male;
  a pool counts for both), and checked neither that an id is one animal of one sex nor that both
  sexes share a stratum. By the engine's own design rules the allowed set was UNIDENTIFIABLE,
  with the pools and without them (b2). Inside the pools a Y-gene comparison would be circular:
  the cells' sex there was assigned from the Y genes (atlas).
- **Xist: a false "metric invalid", an engine defect.** GATE 4's module probe raises Xist 2-fold
  in every cell and measures the metric over all cells. 74.2% of B2's cells have no Xist, where
  the probe changes nothing (b2):

  | Cells | Share with Xist | Response (95% CI) | GATE 4 |
  |---|---|---|---|
  | all (R1's Xist claim) | 25.8% | +0.1627 (+0.1626 to +0.1628) | FAIL |
  | with Xist > 0 (5,027 cells) | 100% | +0.6310 (+0.6304 to +0.6315) | PASS |
  | female (R1's sham; 4,976 cells) | 96.1% | +0.6149 (+0.6144 to +0.6154) | PASS |

  The all-cell response is close to the share of cells with Xist times the response there:
  0.258 × 0.631 = 0.163.

  v1 never tested this path. In the anchors the probe ran on four valid metrics and passed on
  three: R2a's G2M score, R1's Y-gene score and R1's sham, the same Xist metric on the female
  cells only.
- **The sham: correct.**

**v0.1.1 and v0.3 on their real positive controls.** v0.1.1 did not confirm Xist because of
GATE 1, on TMS (step1, section F). v0.3 did not either, because of GATE 4 and the design, on B2.
The manuscript does not say that v0.3 fixed v0.1.1's Xist error: it did not.

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
- **(b) The female arm.** The old females come from the mammary gland only, so a sex × age
  comparison across tissues confounds age with tissue composition. The claims "female detection
  does not decline" and "a sex-by-age interaction" are dropped (the owner's decision,
  2026-10-10). There is no new run within the mammary gland.

## 6. Related work

Short; four paragraphs. The references were checked on 2026-10-10 against the publisher's or
PubMed's record through web search, since Crossref and PubMed could not be read directly from
this environment. Items marked «проверить» need a look at the record before submission.

- **Pseudoreplication.** Cells of one animal are not independent replicates (Hurlbert 1984;
  Lazic 2010). In single-cell differential expression, tests over cells give false discoveries
  that tests over biological replicates avoid (Squair et al. 2021; Zimmerman et al. 2021;
  Crowell et al. 2020). v0.3 infers every effect over replicates. R1 adds the converse risk: a
  donor id need not be one animal.
- **Detection rate as a confounder.** The share of genes a cell detects varies for technical
  reasons and tracks the main axes of variation (Hicks et al. 2018); MAST models it as a
  covariate (Finak et al. 2015). Droplet counts fit sampling without zero inflation (Svensson
  2020), which supports binomial thinning both as the depth correction and as the generator of
  planted truth (Gerard 2020).
- **Spike-ins and normalization.** Scaling normalization assumes that most genes do not change.
  A global change in RNA content needs an external standard (Lovén et al. 2012), and spike-ins
  have their own limits in single cells (Risso et al. 2014; Lun et al. 2017; Vallejos et al.
  2017). Other normalizations (Lun, Bach & Marioni 2016; Hafemeister & Satija 2019) also leave
  this choice open. v0.3 makes it explicit through the estimand: composition is corrected by
  thinning to equal depth, content needs spike-ins, and without them it is UNIDENTIFIABLE (R2c).
- **Benchmarking versus validation** (`paper/drafts/validation_not_benchmarking.md`). Benchmarks
  compare methods across datasets with known truth (Weber et al. 2019; Soneson & Robinson 2018;
  Luecken et al. 2022). metric-autopsy asks whether one number on one dataset measures what the
  claim says: construct validity in its causal sense (Borsboom et al. 2004), argued rather than
  scored (Kane 2013), with negative controls (Lipsitch et al. 2010), equivalence against a SESOI
  (Lakens 2017) and pre-registration (Nosek et al. 2018). The validator itself is tested the
  same way: blind, against planted truth, with real positive controls.

References (verified as written unless marked):
- Borsboom D, Mellenbergh GJ, van Heerden J (2004). The concept of validity. *Psychol Rev*
  111(4):1061–1071. doi:10.1037/0033-295X.111.4.1061
- Buettner F, Natarajan KN, Casale FP, et al. (2015). Computational analysis of cell-to-cell
  heterogeneity in single-cell RNA-sequencing data reveals hidden subpopulations of cells.
  *Nat Biotechnol* 33:155–160. doi:10.1038/nbt.3102 (B3)
- CZI Cell Science Program, Abdulla S, Aevermann B, Assis P, et al. (2025). CZ CELLxGENE
  Discover: a single-cell data platform for scalable exploration, analysis and modeling of
  aggregated data. *Nucleic Acids Res* 53(D1):D886–D900. doi:10.1093/nar/gkae1142 — «проверить»
  the group byline (seen only in search snippets).
- Crowell HL, Soneson C, Germain P-L, et al. (2020). muscat detects subpopulation-specific state
  transitions from multi-sample multi-condition single-cell transcriptomics data. *Nat Commun*
  11:6077. doi:10.1038/s41467-020-19894-4
- Finak G, McDavid A, Yajima M, et al. (2015). MAST: a flexible statistical framework for
  assessing transcriptional changes and characterizing heterogeneity in single-cell RNA
  sequencing data. *Genome Biol* 16:278. doi:10.1186/s13059-015-0844-5
- Gerard D (2020). Data-based RNA-seq simulations by binomial thinning. *BMC Bioinformatics*
  21:206. doi:10.1186/s12859-020-3450-9
- Hafemeister C, Satija R (2019). Normalization and variance stabilization of single-cell RNA-seq
  data using regularized negative binomial regression. *Genome Biol* 20:296.
  doi:10.1186/s13059-019-1874-1
- Hicks SC, Townes FW, Teng M, Irizarry RA (2018). Missing data and technical variability in
  single-cell RNA-sequencing experiments. *Biostatistics* 19(4):562–578.
  doi:10.1093/biostatistics/kxx053 — «проверить» the end page (from citing records).
- Hrovatin K, Bastidas-Ponce A, Bakhti M, et al. (2023). Delineating mouse β-cell identity during
  lifetime and in diabetes with a single cell atlas. *Nat Metab* 5:1615–1637.
  doi:10.1038/s42255-023-00876-x (B2)
- Hurlbert SH (1984). Pseudoreplication and the design of ecological field experiments.
  *Ecol Monogr* 54(2):187–211. doi:10.2307/1942661
- Kane MT (2013). Validating the interpretations and uses of test scores. *J Educ Meas*
  50(1):1–73. doi:10.1111/jedm.12000
- Lakens D (2017). Equivalence tests: a practical primer for t tests, correlations, and
  meta-analyses. *Soc Psychol Personal Sci* 8(4):355–362. doi:10.1177/1948550617697177
- Lazic SE (2010). The problem of pseudoreplication in neuroscientific studies: is it affecting
  your analysis? *BMC Neurosci* 11:5. doi:10.1186/1471-2202-11-5
- Lipsitch M, Tchetgen Tchetgen E, Cohen T (2010). Negative controls: a tool for detecting
  confounding and bias in observational studies. *Epidemiology* 21(3):383–388.
  doi:10.1097/EDE.0b013e3181d61eeb
- Lovén J, Orlando DA, Sigova AA, et al. (2012). Revisiting global gene expression analysis.
  *Cell* 151(3):476–482. doi:10.1016/j.cell.2012.10.012
- Luecken MD, Büttner M, Chaichoompu K, et al. (2022). Benchmarking atlas-level data integration
  in single-cell genomics. *Nat Methods* 19:41–50. doi:10.1038/s41592-021-01336-8
- Lun ATL, Bach K, Marioni JC (2016). Pooling across cells to normalize single-cell RNA
  sequencing data with many zero counts. *Genome Biol* 17:75. doi:10.1186/s13059-016-0947-7
- Lun ATL, Calero-Nieto FJ, Haim-Vilmovsky L, Göttgens B, Marioni JC (2017). Assessing the
  reliability of spike-in normalization for analyses of single-cell RNA sequencing data.
  *Genome Res* 27(11):1795–1806. doi:10.1101/gr.222877.117
- Nosek BA, Ebersole CR, DeHaven AC, Mellor DT (2018). The preregistration revolution.
  *Proc Natl Acad Sci USA* 115(11):2600–2606. doi:10.1073/pnas.1708274114
- Risso D, Ngai J, Speed TP, Dudoit S (2014). Normalization of RNA-seq data using factor
  analysis of control genes or samples. *Nat Biotechnol* 32:896–902. doi:10.1038/nbt.2931
- Soneson C, Robinson MD (2018). Bias, robustness and scalability in single-cell differential
  expression analysis. *Nat Methods* 15:255–261. doi:10.1038/nmeth.4612
- Squair JW, Gautier M, Kathe C, et al. (2021). Confronting false discoveries in single-cell
  differential expression. *Nat Commun* 12:5692. doi:10.1038/s41467-021-25960-2
- Svensson V (2020). Droplet scRNA-seq is not zero-inflated. *Nat Biotechnol* 38:147–150.
  doi:10.1038/s41587-019-0379-5
- Tabula Muris Consortium (2020). A single-cell transcriptomic atlas characterizes ageing tissues
  in the mouse. *Nature* 583:590–595. doi:10.1038/s41586-020-2496-1 (TMS)
- Vallejos CA, Risso D, Scialdone A, Dudoit S, Marioni JC (2017). Normalizing single-cell RNA
  sequencing data: challenges and opportunities. *Nat Methods* 14:565–571.
  doi:10.1038/nmeth.4292
- Weber LM, Saelens W, Cannoodt R, et al. (2019). Essential guidelines for computational method
  benchmarking. *Genome Biol* 20:125. doi:10.1186/s13059-019-1738-8
- Zimmerman KD, Espeland MA, Langefeld CD (2021). A practical solution to pseudoreplication bias
  in single-cell studies. *Nat Commun* 12:738. doi:10.1038/s41467-021-21038-1

The other references of `validation_not_benchmarking.md` (Campbell & Fiske 1959, Cinelli et al.
2022, Cronbach & Meehl 1955, Gagnon-Bartsch & Speed 2012, ICH Q2(R2), Messick 1995, Platt 1964,
Schuemie et al. 2014, Shadish et al. 2002, Spearman 1904, Stevens 1946, Stuart 2010) are not
checked yet: «проверить» each before it is cited.

## 7. Discussion

- The gates are necessary, not sufficient. A verdict is a statement about one dataset under
  the declared assumptions, not about biology.
- The lesson of the dilution: count false-SUPPORTED criteria by level, or only on cards whose
  metric is valid or ambiguous.
- **The development data erred in both directions.** The simulated probe p16 put N8's rate 2.6
  times above the real high-level rate, and no development probe showed GATE 4's failure on a
  marker present in only some cells (§4.7).
- **v1 measured the engine with the oracle's settings.** The GATE 4 dose and the SESOI came from
  the pilot, that is, from the oracle. A user has neither, so v2 must test the engine with its
  defaults, or with the rule by which a user picks them.
- v0.4's targets: per-cell capture (N8), entry-level dropout (N3 f = 0.4), and GATE 4, which
  must measure the response where the signal is claimed (within a group, or over the
  replicates), not over all cells.
- Backgrounds for a sex control: each id one animal of one sex, checked by Xist and the Y genes,
  both sexes in one stratum, the allowed set computed by the engine's own design rules.
- GATE 0's refusals on N7 need explaining (PROJECT.md).
- Step 4, external verdicts, is next.

## Figures (each from a script, written after approval)

1. The design of v1: conditions, cards, the key, the criteria.
2. The criteria against their limits (scores.json).
3. The error rates by expression level: N8, N3 f = 0.4 and the nulls (by-level).
4. TMS: per-mouse `mi_3bin` against median nnz, males separated on both (audit `E.mice`; mice).
5. The real positive controls of both versions, side by side:
   - v0.1.1 on TMS: Xist and the Y-gene score killed at GATE 1 (2.01× in the 20-month stratum),
     though GATE 2 retains 100% and 104% (audit `F`);
   - v0.3 on B2: Xist called invalid at GATE 4 (+0.163 over all cells, +0.631 over the cells
     with Xist, against δ_min 0.25), and the Y genes UNIDENTIFIABLE by the design (b2; r1).
