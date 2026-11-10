# Level B, step S4: the independent verdicts (the protocol's section 7), the rules before their search

`../plans/level_b_protocol.md`, section 7: a verdict is needed for the second question only (rightness); a kept
claim without one stays in the first question. Every verdict found is recorded before the engine runs on the
claim, with its citation and a quote, as one of: **a causal experiment** (a perturbation of the coupling's
proposed cause whose effect on the pair's coordination is measured), **a published re-analysis** (a later
paper that reproduces or refutes the claim), **a replication that does not share the suspected artifact**
(another technology — single-molecule FISH, imaging-based spatial transcriptomics with cell-level counts, a
protein measure of both genes in the same cells — or a design where the artifact cannot arise). A repeat by the
same method on data of the same type is not a verdict. Coded *holds*, *does not hold* or *undecided*, with the
quote that decides it. Nothing below is fetched or read before this file is committed.

## 1. The search (a run in Actions, once, by the label `level-b-verdicts`; the branch `results/level-b-verdicts`)

For every kept paper (S3's list) and each of its kept claims, Europe PMC, every source (preprints included),
no window of dates, every raw response kept with its SHA256:
- **Forward citations of the claim's paper:** its citations endpoint (by its PMID, else its PMCID, else its
  preprint's id), every page, then each citing record's core record (title and abstract), as the search did
  for the method papers (`search.py`'s citations).
- **The claim's genes with the terms:** for a pair A–B, A and B each in the title or the abstract, and one of
  section 3's terms of coordination or of perturbation (knockdown, knockout, overexpression, smFISH, "in situ")
  in the title or the abstract. For a set: its name as a phrase where the paper names it so (a pathway, a
  program, a named module), with the same terms; a set without such a name is searched by its forward
  citations only (said so).
- The claim's own paper is left out of its records.

## 2. Which records are read

- **A fixed filter before any reading:** a record is read for a claim where its title or abstract names both
  genes of the pair (each as a word, case ignored), or, for a set, its name or at least two of its genes. The
  others are counted per claim, not read. (A re-analysis that names the genes only in its full text is missed:
  a limit of the rule, said with the counts.)
- **Two passes**, each by its own session, the second not seeing the first's codes, one session at a time:
  each reads the filtered records' titles and abstracts, and the full text (Europe PMC's open access, fetched in
  the same run for the records a filter keeps) of each record that may hold a verdict; it codes each record as
  a verdict of one of the three kinds (holds / does not hold / undecided, with the citation, the kind and the
  quote that decides it) or as no verdict (with why: a repeat by the same method, a mention without data, about
  other cells or conditions...).
- **A claim's verdict:** the passes' verdicts where they agree; a record they code differently goes to the owner
  with the list, both codes and quotes beside it. Several verdicts of one claim that disagree are all listed;
  none is chosen here.

## 3. What the list shows (S7)

Per claim: the records found (per source), filtered, read; its verdicts with their kind, citation and quote;
the records for the owner. The share of kept claims with a verdict, and of them the share that holds, are the
inputs of π (`real_claims.md`); they are not computed before the list is approved.
