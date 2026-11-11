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

Pass 1's four sessions ran first, one at a time; then its files were moved out of the repository's working tree
(to a directory of the scratchpad that no instruction named, their SHA256 kept) while pass 2's four sessions ran,
one at a time, and came back unchanged (`sha256sum -c`). None was stopped by the usage limit, none was screened
again. Each session's counts, as its final message gave them:

| Pass | Batch | KEEP | UNCLEAR | EXCLUDE (codes) |
|---|---|---|---|---|
| 1 | 01 | 40 | 3 | 124 (X2 2, X4 112, X5 1, X7 6, X8 2, X10 1) |
| 1 | 02 | 22 | 0 | 145 (X2 43, X4 100, X7 2) |
| 1 | 03 | 8 | 0 | 159 (X1 1, X2 70, X4 86, X5 2) |
| 1 | 04 | 6 | 0 | 158 (X1 3, X2 19, X4 133, X5 2, X8 1) |
| 2 | 01 | 41 | 8 | 118 (X2 1, X4 108, X7 6, X8 2, X10 1) |
| 2 | 02 | 20 | 0 | 147 (X2 50, X4 96, X8 1) |
| 2 | 03 | 8 | 0 | 159 (X1 1, X2 49, X4 109) |
| 2 | 04 | 6 | 0 | 158 (X1 4, X2 57, X4 95, X5 1, X8 1) |

**Independence** (`independence.log`, `../../independence.py`): no session's transcript names the other pass's
directory, the directory that held pass 1 during pass 2, a batch other than its own, or the earlier decisions
(`validation/level_b/screening/`, `fulltext_screening`, `claims.md`, `claims.json`, `decision_sheet`,
`set_records.txt`, `JOURNAL.md`, `plan_s2.md`): 8 sessions, all independent.

## The comparison (`merge.json`, `all.tsv`, `to_owner.tsv`, `fulltext.txt`; `screening.py merge`)

- **Agreement over the 665 records:** Cohen's kappa over the three decisions **0.861** (raw agreement 0.970);
  keep-or-unclear against exclude 0.859; where both exclude, the same code in 0.869.
- **The two passes' decisions** (pass 1 / pass 2 → records): EXCLUDE / EXCLUDE 574, KEEP / KEEP 68, KEEP /
  EXCLUDE 8, EXCLUDE / KEEP 7, EXCLUDE / UNCLEAR 5, UNCLEAR / UNCLEAR 3.
- **Per pass:** pass 1 KEEP 76, UNCLEAR 3, EXCLUDE 586 (X1 4, X2 134, X4 431, X5 5, X7 8, X8 3, X10 1); pass 2 KEEP
  75, UNCLEAR 8, EXCLUDE 582 (X1 5, X2 157, X4 408, X5 1, X7 6, X8 4, X10 1). X6 in neither (amendment 1 retired
  it); X11 (not in English) in neither.
- **Sets and global measures:** records whose `kinds` name a set 83 and 91; a global measure 14 and 12; X7 8 and 6,
  both 6.
- **Decided differently or unclear in either: 23 records** (`to_owner.tsv`, both notes). By the deviation of
  `plan_s2.md` (J3) they are not sent to the owner now: with every record either pass keeps or marks unclear,
  **91 records** (`fulltext.txt`), they go to S3, where both passes read their full texts.
