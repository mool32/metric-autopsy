# Level B, step S1: the plan of its two runs (committed before they run)

The owner's decisions of 2026-11-10 (`JOURNAL.md`, J1) and amendment 1 (`../plans/level_b_protocol.md`,
section 12). Both runs are in Actions (this environment reaches neither Europe PMC nor the publishers), each
once, by a label of its own, every response kept as it came on a results branch with SHA256SUMS
(`.github/workflows/level-b-search.yml`); a second run only after a technical failure, recorded in the journal.

## 1. The second fetch: the label `level-b-fetch2`, the branch `results/level-b-fetch2`

The records and their actions are `fetch2_list.tsv` (`search.py fetch2`; a test checks the list is J1's):
- **Europe PMC again** for R00014 (PMC8392141), R00173 (PMC7322959) and R00394 (PMC8728743), which answered
  HTTP 500 to six attempts on 2026-10-19: the first fetch's request (`fullTextXML`), up to six attempts on
  HTTP 429 and 5xx.
- **One attempt at the publisher** for R00438, R00582, R01062, R01338, R02068 and R02134: one GET of
  `https://doi.org/<DOI>`, its redirects followed, with the search's User-Agent (no imitation of a browser,
  no login, no cookie kept); and, where that page names a PDF in a `citation_pdf_url` meta tag, one GET of
  that PDF. Neither is retried.
- **Where the counts could be public**, for R03285 and R03433: NCBI's links from the paper's PMID to GEO
  (`elink`, pubmed to gds), Europe PMC's supplementary files of its PMCID, and for R03285 NGDC's page of the
  BioProject its text names (PRJCA044651). A repository that cannot be reached is recorded, not a failed run.

**Judged after the run by these rules, fixed now; each judgement goes into the journal with its evidence (the
file, its size, what it holds) before any pass reads a text:**
- **A full text is obtained** where a response holds the paper's body: its sections (a Results or a Methods
  heading) with their text beyond the abstract. A paywall, a login or cookie page, an abstract alone, a
  refusal (HTTP 403 and the like), or a file from which no text can be read is not. Obtained: both passes read
  it in S3. Not: the record is excluded (J1).
- **Counts are public** where a response names a count matrix (counts per cell and gene, raw or the authors')
  that anyone can download without asking the authors: a GEO series with such files, an NGDC (OMIX or GSA)
  entry released with them, a supplementary file that is one. Raw reads alone are recorded as such: S5 reads
  count matrices only (the journal records that the protocol's rule 3 says "raw or the authors' count
  matrices" and the full-text instructions "raw reads or the authors' count matrices"). Public: the record is
  kept (both passes found its other rules held). Not: excluded (J1).

## 2. The preprints: the label `level-b-preprints`, the branch `results/level-b-preprints`

- **The query** (`search.py`, `PREPRINT_QUERY`): section 3's two blocks, each term in the title or the
  abstract, and its window, with `SRC:PPR` in place of `LANG:"eng"` (a test checks it is the search's Europe
  PMC query with that one change); cursor pages of 1,000 core records, every page as it came, the records
  parsed as the search parsed Europe PMC's (`preprints.jsonl.gz`).
- **The merge, after the run, here** (a script and its log, committed): the corrected list rebuilt by
  `search.py remerge` from `results/level-b-search` at 7d4938d (its SHA256 as `runs/remerge.SHA256SUMS`), and
  section 5's rule (the same PMID, else the same DOI, else the same title once lower-cased and stripped to
  letters and digits). A preprint already in the list is not screened again; the others are numbered from
  R03507 in the order of their Europe PMC ids and are S2's batches of preprints.

## 3. The screening instructions (12.1 and 12.2), before any pass of S2 and S3

`screening.md` and `screening_fulltext.md` get the set's rules and the language's code in sections dated
2026-11-10; their other rules stay as they were.

## What follows

S2's plan (`plan_s2.md`, its batches and passes) is written after these runs, before its first pass.
