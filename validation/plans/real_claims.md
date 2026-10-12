# Real claims: the main track (draft for the owner)

**If the tool shows no use on real tasks, it is not needed** (the owner's decision of 2026-10-11).
The simulations of v1 and v2 say what the engine does where the truth is planted; this track asks
whether it changes what an analyst would conclude from published data, and whether it is right
when it does. Nothing here is run before the tag of v0.4; the selection, the data and the
verdicts are fixed before any run, so that v0.4 cannot be tuned to these claims.

## What v1 leaves open

On v1's cards (exploratory, `validation/exploratory/v1_analyzers.log`) a test over donors without
gates did almost as well as the engine where the metric works (false "effect" on the pure nulls at
risk 2.9% against 2.1%; sensitivity on valid metrics 99.2% against 97.5%). The engine's unique use
was on metrics that see depth but not the planted coupling: it refuses to interpret them, where a
donor-level test stays silent or, with a TOST, claims an absence. How often real claims rest on
such metrics is not known; that is what this track measures. Its primary outcome is where the
engine changes the donor-level test's conclusion, and on whose side an independent verdict falls.

## Level A: development, not proof

The tool was built around these cases, so it must catch them, but catching them proves nothing.
Each is run, like level B, only after the tag of v0.4, and reported apart from level B.

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
3. **A paired design: interferon-β stimulation** (Kang et al. 2018, Nat Biotechnol 36:89–94; the
   same donors' cells with and without stimulation): the coexpression of interferon-response genes
   in stimulated against control cells, the donor as the pair. The engine supports paired designs
   (v1's N5 used them). The truth is the known induction of the interferon response; the structure
   (donors demultiplexed by genotype, cells per donor, cell type and condition) is checked first.

## Level B: confirmation

**10 to 30 published claims about shifts of coordination or coexpression** between conditions,
selected by rules fixed before any run and before any look at the engine's output on them.

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

## The runs, after the tag of v0.4

On every claim, three analyses: (a) the paper's own, as published (cells or pooled, its statistic);
(b) a test over donors (the donor-level analogue of the paper's metric, a permutation over donors, a
TOST at the SESOI for absence); (d) v0.4 with its own rules, never with settings fitted to the claim.

- **Primary outcome:** the number of claims where (d) changes (b)'s conclusion, a refusal to
  interpret included, and in each of them the side of the independent verdict: (d)'s if it avoids a
  conclusion the verdict rejects or finds one the verdict confirms, (b)'s if it loses a conclusion
  the verdict confirms or makes one it rejects; a verdict that cannot decide is reported, not
  counted.
- **Claims of one paper are not independent** (the owner, 2026-10-12). In the sign test a paper is
  one vote: the side on which most of its deciding disagreements fall, none on a tie; so k and the
  test count papers, and the selection takes claims from as many papers as it can. (The owner's
  other option, a cap on the claims taken from one paper, is the alternative if the first leaves
  too few votes.)
- **Separately:** the claims that fail already at (b). That is pseudoreplication, a known problem,
  and no merit of the tool.
- **The decision rule, fixed before the run.** The gates justify themselves only if the papers whose
  disagreements meet a deciding verdict number at least k = 5 and a one-sided sign test on their
  votes at 0.05 favors (d). With 5 votes that needs all 5 on (d)'s side (0.5⁵ = 0.031); a simple
  majority, 3 against 2, cannot be told from chance (p = 16/32 = 0.5). Otherwise the product
  reduces to a test over donors plus a report of the metric's sensitivity.
- **k = 5** (the owner, 2026-10-12). Five is the smallest number of votes on which even a unanimous
  verdict for (d) can be told apart from chance by that test (four give 0.5⁴ = 0.0625); with fewer
  the track cannot decide and is extended before it is read.
- **The share of true effects among the selected claims**, estimated from their independent
  verdicts with its interval, is compared with v2's threshold π\* in the stratum the claims fall in
  (`v2_requirements.md`, section 1): below π\*, the gates pay off on such claims by v2's exchange;
  above it, a donor-level test does better.

## Now, in parallel with v0.4

Selecting the claims, getting the data, checking their structure and recording the independent
verdicts. The data are read in GitHub Actions, as the Census was for v1 (this environment's network
cannot reach most data hosts). No public write-up of another group's finding before the
confirmatory runs. The plan needs the owner's approval before the selection starts.
