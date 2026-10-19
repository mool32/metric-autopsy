# Level B: the search for claims

- Workflow run 39269131987 of .github/workflows/level-b-search.yml at the pull request's head bc7f3129a7746eb3aaa0761c1d3bb37cad190771; validation/level_b/search.py (blob c591bed785da96f2e30e15043bfbef35ab05c08c) by validation/plans/level_b_protocol.md there (approved by the owner on 2026-10-19).
- search.log and search.json: the queries, each source's counts, the merged list and the sample; requests.tsv: every request (its URL, parameters, status, size and SHA256).
- raw/: every response as it came (gzip): PubMed's esearch by date window and efetch by PMID, Europe PMC's search pages, each method paper's citations and its citing records' core records.
- records.jsonl.gz: the merged records, numbered (rid), each with every source that found it; sample.txt: the rids the screening reads (all of them, or 5,000 drawn with the protocol's seed).
- method_papers.json: the method papers as verified by the lookup (results/level-b-lookup).
