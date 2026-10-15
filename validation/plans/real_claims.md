# Real claims: the main track (draft for the owner)

**If the tool shows no use on real tasks, it is not needed** (the owner's decision of 2026-10-11).
The simulations of v1 and v2 say what the engine does where the truth is planted; this track asks
whether it changes what an analyst would conclude from published data, and whether it is right
when it does. Nothing of the confirmatory set (level B and the controls) is run before the tag of
v0.4; its selection, data and verdicts are fixed before any run, so that v0.4 cannot be tuned to
these claims. A pilot on development cases runs now on v0.3.0's engine, the tag `v0.3.0-prereg`,
without a release (level A; the owner, 2026-10-15).

## What v1 leaves open

On v1's cards (exploratory, `validation/exploratory/v1_analyzers.log`) a test over donors without
gates did almost as well as the engine where the metric works (false "effect" on the pure nulls at
risk 2.9% against 2.1%; sensitivity on valid metrics 99.2% against 97.5%). The engine's unique use
was on metrics that see depth but not the planted coupling: it refuses to interpret them, where a
donor-level test stays silent or, with a TOST, claims an absence. How often real claims rest on
such metrics is not known; that is what this track measures. It asks two questions in order: how
often the engine changes the donor-level test's conclusion, and, only if that is often enough,
whether its gates pay off by v2's exchange, measured on real controls.

## Level A: development, not proof

The tool was built around these cases, so it must catch them, but catching them proves nothing.
They run first on v0.3.0's engine, the tag `v0.3.0-prereg`, in a pilot that runs now, without a
release (the owner, 2026-10-15), in a pull request of its own, with its plan
(`validation/exploratory/real_pilot_plan.md`) committed before the run; and again after the tag of
v0.4, reported apart from level B (the owner, 2026-10-14). The pilot adds one negative
control, the same donors (at least 4) in two technical conditions, if an open dataset has them.
All of the pilot's data are development data and stay outside the confirmatory set for good. Kang
et al. 2018 and Levy et al. 2020 are not in the pilot: they belong to the confirmatory set.

1. **The three documented failures**, where the data are available:
   - *MI coupling on Tabula Muris Senis (FACS).* Scripted already for v0.1.1 by the flagship audit
     (`validation/flagship_audit/audit_tms.py`, `REPORT.md`): in males the lower MI and the lower
     detection separate the same mice (exact p = 0.016 for each, 5 against 4 units).
   - *The entropy anticorrelation (ACP).* Its numbers come from the internal checklist and have no
     script; the data and the code have to be found first. It is a correlation between two metrics,
     so it fits the engine only recast as a claim about groups (for instance, an entropy that
     changes with age).
   - *The spectral exponent of the ECG (β-heart)* is not single-cell data; the engine does not apply.
2. **A positive control: Xist and the Y genes on Tabula Muris Senis**, female against male within an
   age stratum, the donor one mouse. In the FACS data, by the audit, the design without the two
   composite ids is 10 males (4, 2 and 4 at 3, 18 and 24 months) and 6 females (4 and 2 at 3 and 18
   months), plus 3 females at 21 months from the mammary gland only: the 3-month stratum holds 4
   against 4 mice, enough for the exact permutation (smallest p 2 / 70 = 0.029). The structure is
   checked by code independent of the engine before registration, R1's lesson: each id one mouse of
   one sex by its Xist and Y genes, both sexes in one stratum, the tissues of each mouse (v0.4
   change 7 makes tissue a stratum). The droplet data of the atlas, if they hold more mice per sex
   and age, may be the better background; the same checks apply.

Interferon-β stimulation (Kang et al. 2018), listed here before, is now a positive control of the
second question (below).

## Level B: confirmation

**At least 20 papers with published claims about shifts of coordination or coexpression** between
conditions (fewer cannot answer the first question; see the runs), selected by rules fixed before
any run and before any look at the engine's output on them.

Selection rules:
- the claim's metric is of a type v0.4 can check (`v0.4_spec.md`): a pair correlation, a measure of
  coordination of a gene set if v0.4 has the set coupling probe (change 9), or a marker's mean
  expression through the module probe;
- the data are public, with donor (or animal) ids, and at least 3 to 5 donors per group;
- an independent source of the verdict exists and is recorded before the run: a causal experiment,
  a published re-analysis (for claims of noise, of the kind v2 leaves to a later version,
  Ibañez-Solé et al. 2022, eLife, is one), or a replication that does not share the suspected
  artifact: another technology, or a design in which the artifact cannot arise. A repeat by the
  same method on data of the same type repeats the artifact too and is not a verdict (the owner,
  2026-10-12);
- the search that finds the claims (sources, terms, dates) is written down before it is run, and
  every claim it finds is listed with the reason it is kept or left out.

**The first candidate: Levy et al. 2020** (Nat Metab 2:1305–1315), the loss of gene-to-gene
coordination with age, by the global coordination level (GCL).
- *What GCL reads.* In the method papers of the same group (Sci Rep 2021,
  doi:10.1038/s41598-021-90353-w; Amit et al. 2022, Sci Rep, doi:10.1038/s41598-022-11507-y) the
  genes, the top 2,000 expressed in Amit et al. 2022, are split at random into two complementary
  halves, the distance correlation between the halves is computed over cells, and it is averaged
  over 20 to 50 random splits, with a jackknife over cells. This comes from the papers' abstracts
  and search summaries: this environment cannot open the journals, so the gene set Levy et al. used,
  the datasets (listed in their supplementary table), whether the datasets have donor ids, and how
  the ages were compared (over cells, or over animals) are the first things to check against the
  paper.
- *The hypothesis about the tool's use.* GCL's genes are not the few most expressed: in v1's
  oligodendrocyte background (nuclei, 10x 3' v3), 71.1% of the 1,500 most expressed genes are at
  the expression levels where v1's pair metric was blind; in its islet background, 38.9%
  (`validation/exploratory/gcl_gene_levels.log`). If GCL is blind there too, the engine would refuse
  to read a difference of GCL as a difference of coordination where a donor-level test would report
  it. That is testable, and it is the tool's unique use, not a property of donor-level testing.
- *The data.* By the press release (the owner, 2026-10-12) the paper holds about 20 experiments
  from 6 labs: mouse, human and Drosophila; brain, pancreas and haematopoietic stem cells. Each is
  checked against the selection rules (donor ids, at least 3 to 5 donors per group) before it is
  kept, and the paper's unit of comparison, cells or animals, against its Methods.
- *What the engine needs.* GCL as a black-box metric (`metric(data) -> float`) is enough for GATE 0,
  GATE 2 and the effect; GATE 4 needs a probe that couples a set of genes (v0.4 change 9, coded only
  after GCL's definition is checked against the paper and its data pass the selection rules; the
  rest of v0.4 does not wait for it); GATE 5's pair controls do not apply, so validity rests on
  GATE 4. The comparison must be per animal, with enough cells per animal for a distance
  correlation.

## Controls: the truth by construction

Selected by rules of their own and registered, each with its data, claims and structure check,
before it is run (the owner, 2026-10-12). They measure the exchange of the second question on real
data; a control dataset is one unit, its claims taken together.
- **Negative controls:** the same donors measured in two technical conditions (the chemistry, cells
  against nuclei, the protocol, the same libraries sequenced again): the artifacts are real, a
  biological difference there is none. Public data with donor ids and at least 3 to 5 donors; the
  claims formed by a rule fixed before the run, of the metric types of level B.
- **Positive controls:** experimental perturbations with a known response, an effect clearly above
  the SESOI only (the owner, 2026-10-14), such as interferon-β stimulation (Kang et al. 2018, Nat
  Biotechnol 36:89–94: the same donors' cells with and without stimulation, the donor as the pair;
  the claims the coexpression of interferon-response genes, the truth its known induction). The
  structure (donors demultiplexed by genotype, cells per donor, cell type and condition) is checked
  first.

## The runs, after the tag of v0.4

On every claim, three analyses: (a) the paper's own, as published (cells or pooled, its statistic);
(b) a test over donors (the donor-level analogue of the paper's metric, a permutation over donors, a
TOST at the SESOI for absence); (d) v0.4 with its own rules, never with settings fitted to the claim.
A paper is one unit: it counts as changed where (d) changes (b)'s conclusion on more than half of
its selected claims (claims of one paper are not independent; the owner, 2026-10-12).

**First question, frequency (no verdict needed).** The share of papers where (d) changes (b)'s
conclusion, a refusal to interpret included. **Threshold: 0.20, with at least 20 papers** (approved
by the owner, 2026-10-14). It lies midway between "almost never", one paper in ten, and "often",
three in ten, and with 20 papers the rule (at least 4 changed) passes a true frequency of 0.10 with
probability 0.133 and fails one of 0.30 with probability 0.107; with 30 papers 0.073 and 0.077;
with 10 papers 0.264 and 0.149 (`validation/exploratory/real_claims_power.log`, section 2). Below
the threshold the gates change almost nothing on real tasks: the product reduces to a test over
donors plus a report of the metric's sensitivity, and the track stops there.

**Second question, rightness (only above the threshold): v2's exchange on real controls** (the
owner, 2026-10-14). On the controls, where the truth is known by construction, the quantities of
v2's primary outcome (`v2_requirements.md`, section 1), the engine against (b), per claim:
- on negative controls, r: the false conclusions the gates remove, less those they add (an effect
  claimed in either direction is false). An absence is true there and does not enter r; the share
  of (d)'s refusals on them is reported apart;
- on positive controls, l: the true detections the gates lose, less those they find; and r₁: the
  false conclusions they remove there, less those they add. (b)'s "no effect" is a false
  conclusion on a positive control, and an effect in the wrong direction too, so the engine's
  refusal or detection where (b) claims an absence removes one;
- **π\*_real = r / (r + l/4 − r₁)**, with its interval over the control datasets.

π, the share of true effects among the selected papers' claims, is estimated from their
independent verdicts, with its interval. **The gates are justified on real data if π is below
π\*_real, both intervals taken into account.** The exact form of the rule, its operating
characteristics and the number of control datasets it needs are proposed after the pilot. The mix
of negative and positive controls is not fixed: r is taken per negative claim and l and r₁ per
positive one, so π\*_real does not depend on it.

A sign test on the controls was dropped (the owner, 2026-10-14): on v1's cards its outcome was set
by the share of negative controls (at 24 claims per dataset, negative controls gave the engine
1.1% of their deciding votes and positive ones 92.9%; `v1_correction_rule.log`, section 6), and it
weighs a false conclusion removed and a detection lost alike, against the exchange.

**Descriptive: a sign test on the papers with a deciding verdict** (a one-sided test at 0.05, a
paper one vote, the side on which most of its deciding disagreements fall, none on a tie). It does
not decide: with few votes it can say little (`real_claims_power.log`, section 1; the probability
that it passes when (d) is right in a share p of the disagreements):

| deciding votes | (d)'s votes needed | p = 0.80 | p = 0.85 | p = 0.90 |
|---|---|---|---|---|
| 5 | 5 | 0.328 | 0.444 | 0.590 |
| 7 | 7 | 0.210 | 0.321 | 0.478 |
| 10 | 9 | 0.376 | 0.544 | 0.736 |
| 13 | 10 | 0.747 | 0.882 | 0.966 |
| 15 | 12 | 0.648 | 0.823 | 0.944 |
| 20 | 15 | 0.804 | 0.933 | 0.989 |

**Separately:** the claims that fail already at (b). That is pseudoreplication, a known problem, and
no merit of the tool.

**Secondary:** π\*_real against v2's π\* in the same stratum (whether the simulation predicts the
real exchange), and π against v2's π\*. The decision rests on the frequency and on π against
π\*_real.

## Now, in parallel with v0.4

First the pilot on v0.3.0's engine (level A), now and without a release, in a pull request of its
own. Then selecting the claims and the controls, getting the data, checking their structure and
recording the independent verdicts. The data are read in GitHub Actions, as the Census was for v1 (this
environment's network cannot reach most data hosts). No public write-up of another group's finding
before the confirmatory runs. The plan needs the owner's approval before the selection starts.
