# Level B: the protocol for finding the claims

**Status: approved by the owner on 2026-10-19 (the second decisions of that day), with the edits recorded
in section 10.** The draft (2360cc0, a04216d) went to the owner before any search. At this version nothing
of level B has been searched or run: the order is section 11's, and the engine is not run on level B
(`real_claims.md`: nothing of the confirmatory set runs before the tag of v0.4).

Level B is the confirmatory track of `real_claims.md`: at least 20 papers with published claims about
shifts of coordination or coexpression between conditions, selected by rules fixed before any run and
before any look at the engine's output on them. This protocol fixes how those papers are found, how each
is kept or left out, and where each kept claim's verdict comes from — independently of the engine and
before it runs.

## 1. What is searched for

A **claim** is a finding the authors state in the abstract or the results: between two conditions
(ages, disease and control, treated and untreated, genotypes), the coordination of genes changes in a
stated direction. **The claims kept are about pairs**: the correlation or coexpression of a named pair.
Claims about the coordination of a named set, or a global measure of coordination (Levy et al. 2020's
GCL), are **listed and counted at screening, not kept**: they need v0.4's change 9 (the set probe),
which is not coded, and the owner decides about it from those counts alone (section 9). Claims on a
metric of level (a marker's or a module's expression) are not searched for. A **paper** is the unit of
the first question (`real_claims.md`): it counts as changed where the engine changes the donor-level
test's conclusion on more than half of its selected claims.

## 2. Sources

- **PubMed** (MEDLINE and PubMed Central's records), through NCBI's E-utilities (`esearch` by date
  windows of at most 9,999 records, then `efetch` by PMID), the query of section 3, run once.
- **Europe PMC**, its REST search (`/search`, `resultType=core`, cursor pages), the same query in its
  syntax: it adds preprints (bioRxiv, medRxiv) and full text where the licence allows.
- **Forward citations** of method papers, from Europe PMC's citations endpoint, run once, each citing
  record's abstract then fetched from Europe PMC: Levy et al. 2020 (Nat Metab 2:1305–1315; GCL), its
  method papers as `real_claims.md` names them (Sci Rep 2021, doi:10.1038/s41598-021-90353-w; Amit et
  al. 2022, Sci Rep, doi:10.1038/s41598-022-11507-y), CS-CORE (Su et al. 2023, Nat Commun) and
  Skinnider, Squair and Foster 2019 (Nat Methods: measures of association for single-cell
  transcriptomics). CS-CORE and Skinnider et al. were named from memory: every method paper is found in
  Europe PMC by its own query and checked against its expected journal, year and title **before the
  search** (the lookup, section 11; `validation/level_b/method_papers.json`), and their Europe PMC
  records are written into that file before the search runs.
- **The lookup** (workflow run 39268500625, `results/level-b-lookup` at bdb6c45, 2026-10-19), each paper
  one journal article in Europe PMC matching its journal, year and title words or DOI, now its record in
  `method_papers.json`: Levy et al. 2020, MED:33139959 (doi:10.1038/s42255-020-00304-4, cited by 48);
  Sci Rep 2021, MED:34040065 ("A top-down measure of gene-to-gene coordination for analyzing cell-to-cell
  variability", cited by 3); Amit et al. 2022, MED:35534606 ("Global coordination level in single-cell
  transcriptomic data", cited by 5); CS-CORE, MED:37563115 (Nat Commun 2023, doi:10.1038/s41467-023-40503-7,
  "Cell-type-specific co-expression inference from single cell RNA-sequencing data", cited by 91);
  Skinnider et al. 2019, MED:30962620 (Nat Methods, doi:10.1038/s41592-019-0372-4, "Evaluating measures of
  association for single-cell transcriptomics", cited by 127).
- **Not used:** Google Scholar (not reproducible), search engines, recommendations, papers known to us
  from earlier work (Kang et al. 2018 is a positive control of the second question; the pilot's data are
  level A and stay outside), and any list compiled after a look at an engine's output.

The lookup and the search run in GitHub Actions (this environment cannot reach PubMed or Europe PMC), each
once, by a script committed before its run (`validation/level_b/search.py`, its dry run on fabricated
responses on every push); every raw response is committed as it came (gzip) to a results branch with its
SHA256 before any screening. A second run is allowed only if the first failed technically, and the report
says so.

## 3. Terms

PubMed (the conditions' block of the draft is removed: the screening selects by condition, section 5;
"differential correlation" is added):

```
("single-cell"[tiab] OR "single cell"[tiab] OR "single-nucleus"[tiab] OR "single nucleus"[tiab]
 OR scRNA-seq[tiab] OR snRNA-seq[tiab])
AND
(coexpression[tiab] OR co-expression[tiab] OR "gene-gene correlation*"[tiab] OR "gene coordination"[tiab]
 OR "transcriptional coordination"[tiab] OR "coordinated expression"[tiab] OR decoupling[tiab]
 OR "loss of coordination"[tiab] OR "network rewiring"[tiab] OR "differential coexpression"[tiab]
 OR "co-regulation"[tiab] OR "transcriptional coherence"[tiab] OR "differential correlation"[tiab])
AND English[lang] AND ("2017/01/01"[dp] : "2026/09/30"[dp])
```

Europe PMC: the same two blocks, each term searched in the title or the abstract (`TITLE:` or
`ABSTRACT:`), with `LANG:"eng"` and `FIRST_PDATE:[2017-01-01 TO 2026-09-30]`, every source included
(preprints too). Europe PMC does not truncate inside a phrase, so "gene-gene correlation*" is searched as
its singular and plural. The script's query strings are the run's record (`search.py`, `PUBMED_QUERY` and
`EUROPEPMC_QUERY`; a test checks that the PubMed one is this block). PubMed's translation of the query and
its warnings are written to the run's log. Terms of noise and variability ("transcriptional noise",
"variability") are left out: v2 leaves claims of noise to a later version (`real_claims.md`).

## 4. Dates

- **Publication window:** 1 January 2017 to 30 September 2026 (single-cell counts with donor ids at
  scale begin about 2017; the window closes before this protocol, so no paper can be added after a look).
- **The search's date:** the day the approved script runs, recorded with the database's own version
  stamp where it gives one.

## 5. Screening: inclusion

The records of the three sources are merged into one list (the same PMID, else the same DOI, else the
same title once lower-cased and stripped to letters and digits, is one record, with every source that
found it) and numbered in a fixed order (PMID, then DOI, then title). If they are more than 5,000, the
screening reads a random sample of 5,000 of them (section 9); every record stays in the list.

**Two independent passes**, both recorded (a PRISMA-style count at every step), each by its own session:
the second session sees none of the first's decisions. Each pass reads titles and abstracts first, then
the full text of the records it kept (fetched in Actions from Europe PMC where the licence allows; a
record without an accessible full text is marked so and goes to the owner). Every record the two passes
decide differently, and every record either pass marks unclear, goes to the owner; the agreement of the
two passes (Cohen's kappa, at each step) is in the report. A paper is kept where all of these hold:
1. **A claim as defined in section 1** (about a pair), with its direction, in the abstract or the results.
2. **A pair's metric v0.4 can check** (`v0.4_spec.md`): a correlation or coexpression of the pair
   computed from the cells' counts (Pearson or Spearman on normalized counts, mutual information,
   proportionality). Claims about sets or global measures are listed and counted (section 1).
3. **The data are public as counts** (UMI or full-length; raw or the authors' count matrices), on GEO,
   ArrayExpress, CELLxGENE, the Human Cell Atlas or a stable repository named in the paper.
4. **Donor (or animal) ids** in the cells' metadata, or recoverable from the samples, and **enough donors
   for an exact two-sided permutation test to reach p < 0.05** in the claim's comparison: the smallest p
   such a test can give, 2 over the number of relabelings, is below 0.05. For two groups of donors that is
   4 against 4 (2/70) or 3 against 5 (2/56), not 3 against 3 (2/20) or 3 against 4 (2/35); for the same
   donors in both conditions, 6 pairs (2/64), not 5 (2/32).
5. **The comparison is between groups of donors or conditions**, not between cell states within donors
   (trajectories, pseudotime, clusters).

## 6. Exclusion

- No claim with a direction; a claim about composition, expression levels alone, noise or variability;
  bulk or pseudobulk data only; spatial or imaging data only.
- A claim about a set's coordination or a global measure: listed and counted, not kept (section 1).
- Data not public, not as counts, or without donor ids; fewer donors than section 5.4 allows.
- The development data (Tabula Muris Senis as used in level A, SEA-AD's MTG, Marsh et al. 2022) and the
  positive and negative controls of the second question (Kang et al. 2018 and those registered later).
- Duplicates: the same claim on the same data in several papers counts once (the first published).
- Comparisons in TMS FACS of 3-month mice with older ones (age confounded with the source, established
  2026-10-19, `real_claims.md`), unless the source is shown shared.

Every record screened is listed with the pass at which it left and the rule that excluded it; a paper the
rules cannot decide is listed for the owner, never decided by an engine's output.

## 7. The independent source of the verdict, per claim

**A condition of the second question only** (rightness, `real_claims.md`): the first question, the
frequency, needs no verdict, and a kept claim without one stays in it. Every verdict found is recorded
before the engine runs on the claim, with its citation and a quote, as one of:
- **a causal experiment** in the same or another paper (a perturbation of the coupling's proposed cause
  whose effect on the pair's coordination is measured);
- **a published re-analysis** of the claim (a later paper that reproduces or refutes it on the same or
  other data);
- **a replication that does not share the suspected artifact**: another technology (single-molecule
  FISH, imaging-based spatial transcriptomics with cell-level counts, a protein measure of both genes in
  the same cells) or a design where the artifact cannot arise (the conditions measured in one batch, at
  equal depth and detection). A repeat by the same method on data of the same type repeats the artifact
  and is not a verdict (the owner, 2026-10-12).

The verdict is coded as *holds*, *does not hold* or *undecided*, with the quote that decides it. It is
found by a search of its own, written into this protocol before it runs: forward citations of the claim's
paper, and the claim's genes with the terms of section 3 and of perturbation ("knockdown", "knockout",
"overexpression", "smFISH", "in situ").

## 8. How many claims per paper

All claims of section 1 in the abstract, then in the results in order of appearance, up to 5 per paper
(the frequency rule counts a paper by the majority of its claims, so a paper with many claims must not
weigh more). The cap and the order are fixed before screening.

## 9. Stopping, samples and the number of papers

- **Records:** if the merged list holds more than 5,000 records, the screening reads 5,000 drawn without
  replacement by `numpy.random.default_rng(20261019).choice(n, 5000, replace=False)` over the records'
  numbers (section 5), in the numpy of v1's analysis environment (`requirements-panel.txt`); the seed is
  fixed here before the search.
- **Papers:** the target is at least 20 papers (`real_claims.md`, the frequency threshold 0.20 with at
  least 20) and at most 30. If more than 30 are kept, 30 are drawn by
  `numpy.random.default_rng(20261020).choice(n, 30, replace=False)` over the kept papers in the order of
  their records, and the others are listed.
- **Too few:** if fewer than 20 papers with pairs are kept, the report says so; if the papers with claims
  about sets (counted, section 1) would reach the target, the owner decides about change 9 from those
  counts alone. The search is not widened after the screening started.

## 10. The owner's decisions on the draft (2026-10-19, the second decisions)

1. The sources, the window and the terms as drafted, with the conditions' block of the query removed
   (the screening selects by condition), the term "differential correlation" added, and the forward
   citations of Skinnider, Squair and Foster 2019 (Nat Methods) added; CS-CORE and Skinnider et al.
   checked against their Europe PMC records before the search.
2. More than 5,000 records: a random sample of 5,000 is screened, its seed recorded here before the
   search (section 9); every record is in the list.
3. The minimum of donors is a rule, not a number: an exact two-sided permutation test can give p < 0.05
   (section 5.4).
4. Sets of genes: listed and counted at screening, not included; with fewer than 20 papers on pairs and
   enough on sets, the owner decides about change 9 from the counts alone. Levy et al. 2020 leaves the
   first queue (`real_claims.md`).
5. The limits: 5 claims per paper, 30 papers.
6. Screening: two independent passes, the second session not seeing the first's decisions; the
   disagreements and the unclear records to the owner; Cohen's kappa in the report.
7. Metrics of level: not searched for; `real_claims.md`'s rule of selection narrowed to coordination.
8. The independent verdict: a condition of the second question (rightness) only, not of the frequency;
   the verdicts found are recorded before any run.
9. The order of section 11; the search and the screening run alongside the engine's change 17'; the
   engine is not run on level B.

One consequence for the owner (not decided by this protocol): the engine needs at least 3 replicates
per group (`min_replicates`), so a design such as 2 against 8 donors (2/45 < 0.05) passes section 5.4 but
not the engine; such papers are listed apart for the owner.

## 11. The order

1. The dry run on fabricated responses, on every push (`.github/workflows/level-b-search.yml`).
2. The lookup, once (the label `level-b-lookup`): the method papers' Europe PMC records, checked against
   `method_papers.json`, to `results/level-b-lookup`; their ids are then written into
   `method_papers.json` and committed.
3. The search, once (the label `level-b-search`): every raw response, the merged list and the sample,
   with SHA256SUMS, to `results/level-b-search`.
4. The screening (section 5), its two passes, then the full texts of the records kept (fetched in Actions).
5. The list of claims to the owner, each with the reason it is kept or left out, the counts of
   section 1 and the verdicts found (section 7).
6. After the owner's approval of the list: data and structure checks; the runs, after the tag of v0.4.
