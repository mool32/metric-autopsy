# Level B: the protocol for finding the claims (draft for the owner)

**Status: a draft, not approved. Nothing of level B has been searched, downloaded or run** (the owner's
decisions of 2026-10-19: the draft goes to the owner before any search, and nothing of level B is
searched or run before it is approved). The terms, dates and rules below are proposals; every number in
them is a choice to approve or change, not a result.

Level B is the confirmatory track of `real_claims.md`: at least 20 papers with published claims about
shifts of coordination or coexpression between conditions, selected by rules fixed before any run and
before any look at the engine's output on them. This protocol fixes how those papers are found, how each
is kept or left out, and where each kept claim's verdict comes from — independently of the engine and
before it runs.

## 1. What is searched for

A **claim** is a finding the authors state in the abstract or the results: between two conditions
(ages, disease and control, treated and untreated, genotypes), the coordination of genes changes in a
stated direction — the correlation or coexpression of a named pair, the coordination of a named set, or
a global measure of coordination (Levy et al. 2020's GCL). A **paper** is the unit of the first question
(`real_claims.md`): it counts as changed where the engine changes the donor-level test's conclusion on
more than half of its selected claims.

## 2. Sources

- **PubMed** (MEDLINE and PubMed Central's records), through NCBI's E-utilities (`esearch` then `efetch`),
  the query below, run once.
- **Europe PMC**, its REST search (`/search`, `resultType=core`), the same query in its syntax: it adds
  preprints (bioRxiv, medRxiv) and full text where the licence allows.
- **Forward citations** of method papers, from Europe PMC's citation endpoint, run once: Levy et al.
  2020 (Nat Metab 2:1305–1315; GCL), its method papers as `real_claims.md` names them (Sci Rep 2021,
  doi:10.1038/s41598-021-90353-w; Amit et al. 2022, Sci Rep, doi:10.1038/s41598-022-11507-y), and
  one or two methods of coexpression in single cells whose users compare conditions over donors —
  proposed: CS-CORE (Su et al. 2023, Nature Communications), named from memory and to be checked against
  the paper before approval; the owner may add or replace.
- **Not used:** Google Scholar (not reproducible), search engines, recommendations, papers known to us
  from earlier work (Kang et al. 2018 is a positive control of the second question; the pilot's data are
  level A and stay outside), and any list compiled after a look at an engine's output.

The searches run in GitHub Actions (this environment cannot reach PubMed or Europe PMC), once, by a
script committed before the run; their raw results (IDs, titles, abstracts, dates) are committed to a
results branch with their SHA256 before any screening. A second run is allowed only if the first failed
technically, and the report says so.

## 3. Terms

PubMed (to be checked against PubMed's syntax in the script's dry run, with fabricated responses, before
the real run):

```
("single-cell"[tiab] OR "single cell"[tiab] OR "single-nucleus"[tiab] OR "single nucleus"[tiab]
 OR scRNA-seq[tiab] OR snRNA-seq[tiab])
AND
(coexpression[tiab] OR co-expression[tiab] OR "gene-gene correlation*"[tiab] OR "gene coordination"[tiab]
 OR "transcriptional coordination"[tiab] OR "coordinated expression"[tiab] OR decoupling[tiab]
 OR "loss of coordination"[tiab] OR "network rewiring"[tiab] OR "differential coexpression"[tiab]
 OR "co-regulation"[tiab] OR "transcriptional coherence"[tiab])
AND
(aging[tiab] OR ageing[tiab] OR age-related[tiab] OR disease[tiab] OR patients[tiab] OR stimulation[tiab]
 OR treatment[tiab] OR knockout[tiab] OR condition*[tiab])
AND English[lang] AND ("2017/01/01"[dp] : "2026/09/30"[dp])
```

Europe PMC: the same three blocks on `TITLE_ABS`, `LANG:"eng"`, `FIRST_PDATE:[2017-01-01 TO 2026-09-30]`,
with `SRC:PPR` (preprints) included. Terms of noise and variability ("transcriptional noise",
"variability") are left out: v2 leaves claims of noise to a later version (`real_claims.md`).

## 4. Dates

- **Publication window:** 1 January 2017 to 30 September 2026 (single-cell counts with donor ids at
  scale begin about 2017; the window closes before this protocol, so no paper can be added after a look).
- **The search's date:** the day the approved script runs, recorded with the database's own version
  stamp where it gives one.
- **Order:** the protocol approved → the search run → screening → the list of claims, each with the
  reason it is kept or left out, and each kept claim's independent verdict recorded → the owner's
  approval of the list → data and structure checks → the runs, after the tag of v0.4.

## 5. Screening: inclusion

Two passes, both recorded (a PRISMA-style count at every step): titles and abstracts, then full text. A
paper is kept where all of these hold:
1. **A claim as defined in section 1**, with its direction, in the abstract or the results.
2. **A metric v0.4 can check** (`v0.4_spec.md`): a pair's correlation or coexpression computed from the
   cells' counts (Pearson or Spearman on normalized counts, mutual information, proportionality); a
   measure of a set's coordination only if v0.4 has the set probe (change 9 is not coded: such claims are
   listed and deferred); GCL only through change 9 as well (`real_claims.md`, its first candidate).
   `real_claims.md`'s selection rule also admits "a marker's mean expression through the module probe";
   that is a metric of level, not of coordination, and section 1 leaves it out (section 10, item 6).
3. **The data are public as counts** (UMI or full-length; raw or the authors' count matrices), on GEO,
   ArrayExpress, CELLxGENE, the Human Cell Atlas or a stable repository named in the paper.
4. **Donor (or animal) ids** in the cells' metadata, or recoverable from the samples, and **at least 3
   donors per group** in the claim's comparison (`real_claims.md`: at least 3 to 5; the owner to fix 3 or
   5); a paired design (the same donors in both conditions) counts by its pairs.
5. **The comparison is between groups of donors or conditions**, not between cell states within donors
   (trajectories, pseudotime, clusters).
6. **An independent source of the verdict exists** (section 7) and can be recorded before the run.

## 6. Exclusion

- No claim with a direction; a claim about composition, expression levels alone, noise or variability;
  bulk or pseudobulk data only; spatial or imaging data only.
- Data not public, not as counts, or without donor ids; fewer donors than section 5.4.
- The development data (Tabula Muris Senis as used in level A, SEA-AD's MTG, Marsh et al. 2022) and the
  positive and negative controls of the second question (Kang et al. 2018 and those registered later).
- Duplicates: the same claim on the same data in several papers counts once (the first published).
- Comparisons in TMS FACS of 3-month mice with older ones (age confounded with the source, established
  2026-10-19, `real_claims.md`), unless the source is shown shared.

Every record found is listed with the pass at which it left and the rule that excluded it; a paper the
rules cannot decide is listed for the owner, never decided by an engine's output.

## 7. The independent source of the verdict, per claim

Recorded before the engine runs on the claim, with its citation and a quote, as one of:
- **a causal experiment** in the same or another paper (a perturbation of the coupling's proposed cause
  whose effect on the pair's or set's coordination is measured);
- **a published re-analysis** of the claim (a later paper that reproduces or refutes it on the same or
  other data);
- **a replication that does not share the suspected artifact**: another technology (single-molecule
  FISH, imaging-based spatial transcriptomics with cell-level counts, a protein measure of both genes in
  the same cells) or a design where the artifact cannot arise (the conditions measured in one batch, at
  equal depth and detection). A repeat by the same method on data of the same type repeats the artifact
  and is not a verdict (the owner, 2026-10-12).

The verdict is coded as *holds*, *does not hold* or *undecided*, with the quote that decides it. It is
found by a search of its own, written into the protocol's approved version before it runs: forward
citations of the claim's paper, and the claim's genes with the same terms of section 3 and of
perturbation ("knockdown", "knockout", "overexpression", "smFISH", "in situ"). A claim without an
independent verdict is left out of level B (section 5.6) and listed.

## 8. How many claims per paper

All claims of section 1 in the abstract, then in the results in order of appearance, up to 5 per paper
(a proposal: the frequency rule counts a paper by the majority of its claims, so a paper with many claims
must not weigh more). The cap and the order are fixed before screening.

## 9. Stopping and the number of papers

The target is at least 20 papers (`real_claims.md`, the frequency threshold 0.20 with at least 20). If
the search yields fewer, the report says so and the owner decides; the search is not widened after the
screening started. If it yields more than 30, 30 are drawn at random with a seed fixed before the draw
(as v1's key, `validation/prereg/v1.md`), and the others are listed.

## 10. What the owner decides before the search

1. The databases, the terms and the window of sections 2–4, and the method papers whose forward
   citations are searched (section 2: Levy et al. 2020 and its two method papers; CS-CORE or another).
2. The minimum donors per group: 3 or 5.
3. Whether claims about a set's coordination wait for change 9 or are left out.
4. The cap of claims per paper and the maximum number of papers.
5. Who screens: the agent by these rules, every undecided record to the owner, or a second screener.
6. Whether claims on a metric of level (a marker's mean expression, through the module probe), which
   `real_claims.md`'s selection rule admits, are searched for too: then section 1's definition and
   section 3's terms widen, or the rule in `real_claims.md` narrows to coordination.

How many records each term brings is not known before the run (no database is reachable from here;
the script's dry run uses fabricated responses only), so no term is judged by its yield before the
owner approves the list.
