# Validating the validator: a blind, pre-registered test of artifact checks for single-cell metrics

*Draft of 2026-10-10: the abstract and §4 only, for the owner. The title is the first of the
plan's three options. The other sections wait for the owner's reply. The plan, with every
number's source, is [`paper/drafts/v03_manuscript_plan.md`](drafts/v03_manuscript_plan.md).
`paper/manuscript.md` is the superseded v0.1.1 preprint.*

Theodor Spiro | ORCID 0009-0004-5382-9346

## Abstract

A metric that differs between single-cell conditions is often read as biology, although
sequencing depth, dropout and the study design can produce the same difference. metric-autopsy
runs a pre-registered claim through gates aimed at named artifacts and decides one verdict from
four fields: metric validity, design adequacy, an effect over biological replicates at equal
depth, and replication. Its first release (v0.1.1) erred in both directions on its own probes:
it passed a metric that returns random numbers and blocked a known sex difference. We rebuilt
the logic (v0.3) and froze it before a blind test. In a pre-registered blind validation on two
real CELLxGENE backgrounds with planted truth (9,750 datasets, 10,540 claim cards), v0.3 met
every error criterion within its scope: one metric, two backgrounds, the binomial-thinning
family of artifacts and per-cell variable capture. Our prediction that it would fail on per-cell
variable capture was not confirmed. In an exploratory reading by expression level, its rate of
false support there was 4.3% at high expression, against a nominal 2.5%. v0.3 did not deliver
its only real positive control, the sex difference in Xist and the Y genes in mouse islets.
The background allowed the sexes to be compared only inside pooled samples, which the protocol
had missed. GATE 4 measured an injected Xist signal over all cells, three quarters of which
lack Xist, and so called a valid metric invalid. The development data had erred in both
directions: they overstated one error rate 2.6-fold and hid the GATE 4 defect. A validator must
therefore be tested blind, on real backgrounds, with real positive controls and with the
settings a user would choose. We report the scope, the failures and the targets of the next
version.

## 4. A blind, pre-registered validation of v0.3

### 4.1 Design

The engine was frozen as the tag `v0.3.0-prereg` (commit 90a7718) together with the protocol
(`validation/prereg/v1.md`). One workflow run executed the protocol from the tag, in one attempt
(workflow run 37990113304). The records it fixed before drawing the key are the run tag
`panel-v1-run` (commit 7782254). The key that chose the gene pairs and the claim cards was the
drand quicknet beacon's round 32929913, published after those records.

**Backgrounds.** A rule fixed in advance chose two backgrounds from the CELLxGENE Census
(release 2025-11-08; CZI Cell Science Program et al. 2025):
- **B1:** human oligodendrocytes of the MSSM cohort; nuclei, 10x 3' v3;
- **B2:** mouse pancreatic islet beta cells (Hrovatin et al. 2023); cells, 10x 3' v2 and v3.
  Its 50 donor ids are samples: 6 of them are pools of both sexes, and all single females are
  of the NOD strain (§4.7).

The suspension type and assay come from the Census metadata (workflow run 38033204255).

**Planted truth.** Known truth was planted into the backgrounds' real counts by binomial
thinning (Gerard 2020). Each background holds 24 gene pairs, 8 at each of three expression
levels. The one metric is the log-normalized Pearson correlation of a pair. The conditions:
- **Nulls** (N1–N8): donor-label permutation (N1); N1 with capture loss (N2) or extra dropout
  (N3) in one group; pseudoreplication, 3 against 3 donors (N4); a within-donor sham (N5);
  useless metrics: a random number, a constant, a score of random genes (N6a–c); a label
  permutation over B2's donor ids (N7), which is a null at the level of the id, not of the
  mouse; and N1 with capture that varies from cell to cell (N8), outside the engine's
  correction family.
- **Real effects** (E1–E3): an injected coupling of the pair in one group (E1), alone or with a
  capture loss against it (E2) or with it (E3).

In all there are 9,750 datasets and 10,540 claim cards. For each card an oracle, which never
imports the engine, measured the truth about the metric on its pair: **valid**, **blind** or
**ambiguous** by its response to GATE 4's injection, at least 1.2 δ_min, at most 0.8 δ_min, or
in between; the metrics of N6 are **useless**.

**Criteria.** Eight error criteria (S1–S7b) were each set by one principle: a sound validator
passes with probability at least 0.90, and one with twice the error passes with probability at
most 0.05. A sound validator passes all of them together with probability 0.912 (`oc.log`).
- S1 bounds false SUPPORTED on each of the five key nulls;
- S2 bounds false SUPPORTED in the other null groups and on the real effects;
- S3 asks for correct definite outcomes where the truth can be established;
- S4 bounds false answers "opposite direction" and "explained by depth";
- S5 bounds false "metric invalid" on a valid metric;
- S6 bounds false NO DETECTABLE EFFECT;
- S7a allows no verdict that the engine's rules cannot give;
- S7b bounds crashes.

Three anchors on real data whose truth is public (R1–R3) are reported one by one, not pooled
into rates.

### 4.2 Result

v0.3 passed every criterion (`scores.txt`, branch `results/panel-v1`, commit adc7d65):

| Criterion | Count | Allowed |
|---|---|---|
| S1 | N1 5, N2 6, N5 8, N6c 0, N8 11 (of 790 each) | ≤ 29 each |
| S2 | N2's other steps 5, N3 11, N4 3, N6a 0, N7 5 (of 790 each); real effects 2 of 1,285 | ≤ 30, 31, 36, 29, 30; ≤ 51 |
| S3 | real effects 181/189, nulls 986/1,075, invalid metrics 4,674/4,685; 16 of 16 conditions | ≥ 144, 778, 3,332 |
| S4 | real effects 0/515, 0/364; nulls 58/2,165, 0/1,147; invalid metrics 1/7,020, 0/4,261 | ≤ 18, 26; 102, 97; 320, 393 |
| S5 | real effects 4/206, nulls 92/2,252 | ≤ 32, 439 |
| S6 | real effects 0/121, invalid metrics 5/7,020 | ≤ 15, 660 |
| S7a | 0 of 10,540 | 0 |
| S7b | 0 of 10,540 | ≤ 18 |

We had predicted that v0.3 would fail S1 on N8, since per-cell variable capture is outside its
correction. The prediction was not confirmed: 11 false SUPPORTED of 790, where S1 allows 29.

### 4.3 Error rates by expression level

This reading is exploratory, not a criterion, and was decided after the results
(`validation/exploratory/v1_by_level.log`, from `scores.json`).
- On these backgrounds the metric is valid only at high expression. All 2,458 cards with a
  valid metric are at the high level. At the medium and low levels the metric is blind (55
  medium-level cards are ambiguous), and no card there was SUPPORTED. So only 253 of N8's 790
  cards could carry a false SUPPORTED.
- N8 at the high level gave a false SUPPORTED on 11/253 = 4.3% (95% CI 2.2–7.6%), at a nominal
  2.5%. It gave SUPPORTED or an answer against the direction on 22/253 = 8.7%, at a nominal 5%.
- N3 at f = 0.4 at the high level gave 7/83 = 8.4% and 10/83 = 12.0%.
- For comparison, N1, N5 and N7 at the high level gave SUPPORTED or an answer against the
  direction on 5.9%, 3.8% and 4.0%. These pure nulls together gave a false SUPPORTED on
  18/812 = 2.2%.
- Had N8 kept its high-level rate on all 790 cards, about 34 false SUPPORTED would be expected
  against S1's 29, and S1 would pass with probability ≈ 0.20. The pass on N8 owes to the
  blind pairs that dilute the denominator.

### 4.4 Limitations

- **The scope** is one metric, two backgrounds, the binomial-thinning family and N8. Nothing
  else is validated, `mi_3bin` included.
- **v0.3 does not correct capture that varies from cell to cell.** Its depth correction brings
  whole groups to one depth.
- **The design effects reach 5.62,** because datasets share donors (`scores.json`,
  `limitations`). On other donors the error rates are known less precisely.
- **GATE 0 refused 26 of N7's 276 high-level cards (9.4%).** Whether random groups of donor ids
  differ in depth, or GATE 0 is oversensitive, is open.
- **GATE 4's module probe measures the response over all cells.** For a marker present in only
  some cells it can call a valid metric invalid (§4.7). v1's panel used the module probe only on
  the useless metric N6c, so v1 did not test this path.

### 4.5 Checks after the run

- `score.py`, run again on the published reports, reproduces every criterion.
- `blind.py verify` rebuilt the first 20 datasets on another CPU model. All 20 datasets are
  identical, and all 22 reports agree: 16 byte for byte, 6 within 1e-9, with the same verdicts
  and causes (`validation/exploratory/v1_verify_first20.log`).

### 4.6 The anchors R2 and R3

R2 compares mouse embryonic stem cells sorted by DNA content (B3; Buettner et al. 2015),
which carry ERCC spike-ins:
- **R2a**, the G2M score, G2M against G1: INCONCLUSIVE;
- **R2b**, total RNA with the spike-ins: INCONCLUSIVE;
- **R2c**, total RNA without the spike-ins: UNIDENTIFIABLE.

Each is in its allowed set. R2's known difference could not be confirmed by design: each
phase is one capture batch, so there are no replicates. R3 was dropped with its background,
whose deposited counts are not raw.

### 4.7 The only real positive control

**The anchor.** R1 asks for a known difference on B2: the sex markers, female against male
mice, within the development stage, with the donor id as the mouse. The claims:
- Xist decreases from female to male;
- the Y-gene score (Ddx3y, Eif2s3y, Kdm5d, Uty) increases;
- a sham compares female mice split at random.

The metric is the mean over cells of the marker's log1p CP10k, with a SESOI of 0.5 and so
δ_min = 0.25. The allowed
outcomes were SUPPORTED for both markers, and NO DETECTABLE EFFECT or INCONCLUSIVE for the
sham. R1 did not run in v1, because of a loading bug. It ran once after the results, outside
v1, with only the loader fixed (DEVIATIONS.md D1; workflow run 38033821978; branch
`results/panel-v1-r1`).

**B2's structure** (`validation/exploratory/b2_sex_structure.log`, workflow run 38045140802).
B2's 50 donor ids are samples. By the share of their cells with Xist > 0:
- **34 single males**, Xist in 0.0–1.5% of their cells and a Y gene in 68.0–94.2%: Fltp_adult 4,
  STZ 7, VSG 8, spikein_drug 15;
- **10 single females**, 98.0–100.0% and 0.0–0.7%: NOD 3, NOD_elimination 7, all of the NOD
  strain;
- **6 pools of both sexes**, 53.0–80.0% and 33.8–45.8%: Fltp_P16 at 2 weeks, Fltp_2y at 20
  months and over. Of each pool's 400 cells, 227–267 are annotated female and 133–173 male.

The source atlas agrees. Its sample table lists the six pools as "mixed" and all twelve of its
NOD samples as female (B2 holds ten of them). In the pools each cell's sex was assigned from a score of the Y-chromosome
genes (the atlas's "data-driven" sex annotation; Hrovatin et al. 2023, code at
theislab/mouse_cross-condition_pancreatic_islet_atlas, commit 3af65e4). No development stage
holds a single male and a single female:

| Development stage (Census) | Single males | Single females | Pools of both sexes |
|---|---|---|---|
| 2 weeks | — | — | 3 (Fltp_P16) |
| 5 weeks | — | 3 (NOD) | — |
| 8 weeks | — | 1 (NOD_elimination) | — |
| 2 months | 15 (spikein_drug) | — | — |
| 3 months | — | 6 (NOD_elimination) | — |
| 4 months | 12 (Fltp_adult, VSG) | — | — |
| 6 months | 7 (STZ) | — | — |
| 20 months and over | — | — | 3 (Fltp_2y) |

Within a stage the sexes meet only inside the pools.

**The outcome.** Two of R1's three claims are outside their allowed sets:

| Claim | Verdict | Allowed |
|---|---|---|
| Xist | NOT SUPPORTED: metric invalid (GATE 4: +0.1627, 95% CI +0.1626 to +0.1628, against δ_min 0.25) | SUPPORTED |
| Y genes | UNIDENTIFIABLE: the donor ids are partially crossed with sex | SUPPORTED |
| Sham, female against female | NO DETECTABLE EFFECT | NO DETECTABLE EFFECT, INCONCLUSIVE |

The two failures have different causes.

**The Y genes: UNIDENTIFIABLE, which reads B2 correctly; the protocol erred.** Within a stage
the sexes can be compared only inside the pools, so neither a nested nor a paired comparison of
mice is valid. The engine's design rules say so:
- with the pools, the ids are partially crossed with sex (`effect._replicate_design`);
- without them, no stage holds both sexes, and GATE 1 stops.

Inside the pools a comparison of the Y genes would be circular, since the cells' sex there was
assigned from the Y genes. The protocol's allowed set, SUPPORTED, was wrong. It counted the ids
of each sex over all stages (16 female and 40 male: a pool counts for both). It checked neither
that an id is one animal of one sex nor that both sexes share a stratum. The design of the Xist
claim is UNIDENTIFIABLE for the same reason, but metric validity decides first.

**Xist: a false "metric invalid", an engine defect.** GATE 4's module probe raises Xist 2-fold
in every cell (the other genes are thinned to half, against a sham that thins every gene) and
measures the metric over all cells. 74.2% of B2's cells have no Xist, and there the probe
changes nothing. The same probe, with the same 200 injections, gives on three sets of cells:

| Cells | Cells with Xist | Response (95% CI) | GATE 4 |
|---|---|---|---|
| all: R1's Xist claim | 25.8% | +0.1627 (+0.1626 to +0.1628) | FAIL |
| with Xist > 0 (5,027 cells) | 100% | +0.6310 (+0.6304 to +0.6315) | PASS |
| female: R1's sham claim (4,976 cells) | 96.1% | +0.6149 (+0.6144 to +0.6154) | PASS |

So the metric responds to Xist where Xist is, and the all-cell response is diluted by the cells
without it: it is close to the share of cells with Xist times the response there,
0.258 × 0.631 = 0.163. v1 never tested this path, because its panel used the module probe only
on the useless metric N6c. In the anchors the probe ran on four valid metrics and passed on
three: R2a's G2M score, R1's Y-gene score, and R1's sham, which is the same Xist metric on the
female cells only.

**The sham was right.** Female against female mice, it gave NO DETECTABLE EFFECT (TOST within
±0.5).

**Both versions on their real positive controls.** v0.1.1 did not confirm Xist because of
GATE 1. On Tabula Muris Senis it killed Xist and the Y-gene score at GATE 1 (2.01× in the
20-month stratum), although its own GATE 2 retained 100% and 104% of the effects
(`validation/flagship_audit/REPORT.md`, section F). v0.3 did not confirm Xist either, because
of GATE 4 and the design, on B2. v0.3 did not fix v0.1.1's error on this control. The targets
of v0.4 are a GATE 4 that measures the response where the signal is claimed (within a group, or
over the replicates), and backgrounds for a sex control whose ids are each one animal of one
sex, with both sexes in one stratum.

## References (of the abstract and §4)

- Buettner F, Natarajan KN, Casale FP, et al. (2015). Computational analysis of cell-to-cell
  heterogeneity in single-cell RNA-sequencing data reveals hidden subpopulations of cells.
  *Nat Biotechnol* 33:155–160. doi:10.1038/nbt.3102
- CZI Cell Science Program, Abdulla S, Aevermann B, Assis P, et al. (2025). CZ CELLxGENE
  Discover: a single-cell data platform for scalable exploration, analysis and modeling of
  aggregated data. *Nucleic Acids Res* 53(D1):D886–D900. doi:10.1093/nar/gkae1142 —
  «проверить» the group byline.
- Gerard D (2020). Data-based RNA-seq simulations by binomial thinning. *BMC Bioinformatics*
  21:206. doi:10.1186/s12859-020-3450-9
- Hrovatin K, Bastidas-Ponce A, Bakhti M, et al. (2023). Delineating mouse β-cell identity during
  lifetime and in diabetes with a single cell atlas. *Nat Metab* 5:1615–1637.
  doi:10.1038/s42255-023-00876-x
