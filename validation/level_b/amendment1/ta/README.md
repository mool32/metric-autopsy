# Level B, amendment 1, S2: titles and abstracts of the sets' records and the new preprints — the record

The plan is `../../plan_s2.md` (81feb76, committed before the first pass); the instructions are `../../screening.md`
as amended on 2026-11-10; the tools are `../../screening.py`. Everything below happened on 2026-11-10, in this
order.

## The batches

665 records (the 210 sets' records of `../set_records.txt`, then the 455 new preprints of
`../preprints/preprints_new.jsonl.gz`), 4 batches of 167, 167, 167 and 164 (`batches_manifest.json`, with each
batch's SHA256 and the inputs').

## The passes

Each batch was screened twice, once in each pass, by a fresh agent session, one session at a time, with the
task text of `plan_s2.md` (section 4) word for word but for its pass, batch, size and paths. `pass1/` and
`pass2/` hold the decisions as the screeners wrote them; each file passed `screening.py check-one` when its
session ended.

(The sessions, their counts and the independence check are added when both passes are complete.)
