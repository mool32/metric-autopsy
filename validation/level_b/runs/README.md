# Level B: the logs of the lookup and the search

Copies of the runs' logs, each identical to its results branch's file (SHA256 as in that branch's
SHA256SUMS); the raw responses stay on the branches.

| file | from | SHA256 |
|---|---|---|
| `lookup.log`, `lookup.json` | `results/level-b-lookup` at bdb6c45 (workflow run 39268500625) | `c857bdeea9bbfac84e4fd347d8a0d7b040756058a5820494b5a3ffb82f6dd157`, `6b7dc28900a4f8b0dc440f9f3597a691c43d805efada71810fa2779e23eb526d` |
| `search.log`, `search.json` | `results/level-b-search` at 7d4938d (workflow run 39269131987) | `0a4fa3c5e3b1dbdb1af1ca65bb97467920eb101544093e56c7c6a87c1c660a3b`, `2717be8325459cda955315658d016108678aa8d92f7d479ac7b8eaa0ea7267b3` |
| `remerge.log`, `remerge.json` | `search.py remerge` on that branch, here (below) | `remerge.SHA256SUMS` |

The search's merged list, `records.jsonl.gz` (SHA256 `18b051a371d09fc7f8886eb5ef7be37e296a47ed4b9160f207affefe5908dd40`), holds 3,394 records; its sample,
`sample.txt` (`64acb6384d1509ff3feb59c9a15638c860d12b22e96b88bf8f200ffaa3a6915c`), is every record (no more than 5,000). The first
screening read them (`../screening.md`, `../screening.py batches`).

## The corrected list (`remerge.log`, `remerge.json`, `remerge.SHA256SUMS`)

**That list is wrong, by a bug of `search.py`'s PubMed parser** (found on 2026-10-19, after the first
screening): it took the record's ids from every `ArticleId` under `PubmedData`, and PubMed's `ReferenceList`
lists the references' ids there too, so a record's DOI and PMCID were those of its last reference (of
PubMed's 3,359 records, 2,183 got another DOI and 2,596 another PMCID). Records whose last references shared
a DOI were merged as one paper: 64 records of the list held 176 papers (2 to 10 PMIDs each), and 112 papers
were never shown to a screener. The parser now reads `PubmedData/ArticleIdList` only (a test with
references' ids, and fabricated records whose references carry the next record's ids, guard it). The
numbers are `remerge_compare.log` (`../remerge_compare.py` on both lists).

The list is rebuilt **from the run's own responses, nothing fetched**: `search.py remerge` runs the
search's code again with each of the 51 recorded requests answered by the response recorded for it (matched
by method, URL, parameters and SHA256, in `requests.tsv`'s order). The corrected list holds **3,506
records** (6,762 from the sources, as before; 45 preprints only, as before). 47 of them hold two PMIDs
(one three), and in each one of the PubMed records is a preprint (its type, or a preprint server's DOI): in
46 the records share one title; the 47th is joined through its preprint's Europe PMC record, which has the
PubMed preprint's DOI and the journal version's later title. Its
`records.jsonl.gz` is reproducible byte for byte (two replays, the same SHA256
`48709ef7af34992616a929f498b0326f0aca58721cdf9ed5735af9fae4b523ca`; its content's, unzipped,
`e184565df0982a0fdc535738328c08fa209c9ad954acdca044c04d9583e2da15`); `sample.txt` (`c40f26152d75ca367c1fa6a108a4eb3857f6d7d94f130f15dbae36ad21fa5eeb`) is every
record. To reproduce: `python -I validation/level_b/search.py remerge --search <results/level-b-search at
7d4938d> --out <dir>`.

Every record of the first list keeps its text as the screeners read it (title, journal, year, types,
preprint, abstract): the 3,394 texts are all in the corrected list, unchanged, so the screening's decisions
on them carry over (`../screening.py carry`), and the 112 new records are screened by two new independent
passes (`../screening/README.md`). The full-text list is rebuilt from the corrected list (the first one's
PMCIDs were the references').

`search_checks.log` (`../search_checks.py` on `results/level-b-search`), two findings for the owner, not
technical failures of the run:

1. **Europe PMC returned no preprint.** Its query's pages hold 3,157 records, MED 3,118, AGR 28, PMC 11,
   PPR 0, every one with language "eng". Europe PMC's preprint records carry no language (73 of 73 among
   the forward citations' core records), so the protocol's `LANG:"eng"` left every preprint out, while the
   protocol meant "every source included (preprints too)". The 45 records that are preprints only came
   through the forward citations.
2. **PubMed dropped the phrase "loss of coordination"**: it is not in PubMed's phrase index
   (`quotedphrasesnotfound`), so PubMed's query ran without that term; Europe PMC's held it.
