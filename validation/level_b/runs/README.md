# Level B: the logs of the lookup and the search

Copies of the runs' logs, each identical to its results branch's file (SHA256 as in that branch's
SHA256SUMS); the raw responses stay on the branches.

| file | from | SHA256 |
|---|---|---|
| `lookup.log`, `lookup.json` | `results/level-b-lookup` at bdb6c45 (workflow run 39268500625) | `c857bdeea9bbfac84e4fd347d8a0d7b040756058a5820494b5a3ffb82f6dd157`, `6b7dc28900a4f8b0dc440f9f3597a691c43d805efada71810fa2779e23eb526d` |
| `search.log`, `search.json` | `results/level-b-search` at 7d4938d (workflow run 39269131987) | `0a4fa3c5e3b1dbdb1af1ca65bb97467920eb101544093e56c7c6a87c1c660a3b`, `2717be8325459cda955315658d016108678aa8d92f7d479ac7b8eaa0ea7267b3` |

The search's merged list, `records.jsonl.gz` (SHA256 `18b051a371d09fc7f8886eb5ef7be37e296a47ed4b9160f207affefe5908dd40`), holds 3,394 records; its sample,
`sample.txt` (`64acb6384d1509ff3feb59c9a15638c860d12b22e96b88bf8f200ffaa3a6915c`), is every record (no more than 5,000). The screening reads them
(`../screening.md`, `../screening.py batches`).

`search_checks.log` (`../search_checks.py` on `results/level-b-search`), two findings for the owner, not
technical failures of the run:

1. **Europe PMC returned no preprint.** Its query's pages hold 3,157 records, MED 3,118, AGR 28, PMC 11,
   PPR 0, every one with language "eng". Europe PMC's preprint records carry no language (73 of 73 among
   the forward citations' core records), so the protocol's `LANG:"eng"` left every preprint out, while the
   protocol meant "every source included (preprints too)". The 45 records that are preprints only came
   through the forward citations.
2. **PubMed dropped the phrase "loss of coordination"**: it is not in PubMed's phrase index
   (`quotedphrasesnotfound`), so PubMed's query ran without that term; Europe PMC's held it.
