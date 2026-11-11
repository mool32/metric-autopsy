# Level B, S4: the independent verdicts, both passes — the record

The plan is `../plan_s4.md` (36c11cd, committed before the search; its precisions 046275b, 8882583 and 4dd6653, each
before the search's run or its outcome); the instructions `../verdicts_coding.md` (8d62a76, before any record was
read); the tools `../verdicts.py` (the search), `../verdicts_codes.py` (the batches, the check, the comparison) and
`../independence.py`. The journal's entries J23, J30 and J31 say what happened, in that order, on 2026-11-11.

## The batches

`batches/` (b33f141, committed before either pass): 48 lines (claim, record), R00189's 34 in batches 01 and 02,
R01225's 10 and R02693's 4 in batch 03; 26 full texts as text (`V00001`–`V00026`, outside the repository, their
SHA256 in `batches/texts.sha256`). Every set's claim, and R00189:2 and 4, R01225:1, 4 and 5, R01962's two and
R02693:1, have no line: the search's filter kept no record for them (J30). R00189:1's 185 kept records exceed the
cap of 100: its 26 naming one of R00189's terms are read, the 159 others counted, not read (`verdicts_terms.json`).

## The sessions

Each batch was coded once in each pass by a fresh session with `plan_s4.md`'s task text word for word, one session at
a time: pass 1's batches 01–03, then pass 2's; each pass's files out of the working tree while the other coded (back
unchanged by `sha256sum -c`). Every file passed `verdicts_codes.py check-codes`; none was coded again, none stopped by
the usage limit. In pass 2's batch 02 one response of the session, after its ninth line was written, was stopped by
the harness's classifier (the harness's own message in the transcript); the session looked at its file and went on
in the same session (J31). Each session's counts, from its file (`merge/sessions.tsv`):

| Pass | Batches | Lines | VERDICT | NO_VERDICT |
|---|---|---|---|---|
| 1 | 01–03 | 48 | 0 | 48 |
| 2 | 01–03 | 48 | 0 | 48 |

`decisions.sha256`: the six files as committed.

## Independence

`independence.log`, two parts:

1. **Strict, as `plan_s4.md` section 2 states it:** 6 sessions to check. Each transcript names `claims.md`,
   `claims.json`, `fulltext_screening`, `JOURNAL.md` and the plans — all of them inside one thing: the harness's task
   reminder (attachments of type `task_reminder`, 1 to 6 per session), which put this session's own task list (20
   items, their subjects and descriptions) into each coder's context. Its one content, the same in the six sessions,
   is `task_reminder.jsonl`: it names paths, and holds no decision's content, no code of either pass, no batch and no
   text. S2's and S3's 28 sessions had no such attachment.
2. **The task list counted apart** (`--exclude-attachment task_reminder`, added after the passes for this, J31): all
   independent — no session names the other pass's directory, the directory that held pass 1's files, a batch or a
   text outside its own list, or an earlier decision; the plans and the protocol the instructions cite are named in
   no session's own call.

The criterion as the plan states it is not met: **a deviation, for the owner** (J31).

## The comparison

`merge/` (`verdicts_codes.py merge`; `merge/merge.log`): both passes code all 48 lines `NO_VERDICT`. The agreement of
the decisions and of the full codes is 1.000; Cohen's kappa is undefined (one category in both passes, the chance
agreement 1). No verdict is agreed (`merge/verdicts.tsv` holds its header only) and no line is coded differently
(`merge/to_owner.tsv`, its header only). The passes' notes: other cells or conditions, the genes' levels only, other
outcomes (perturbations read out on viral replication), the same method, a review, and one record both passes call
the claim's own paper (PPR:PPR182850, a ChemRxiv preprint of R00189's). Per claim, `merge/merge.json`.

No share of the claims with a verdict is computed here: π's inputs wait for the owner's approval of the list
(`plan_s4.md`, section 3).
