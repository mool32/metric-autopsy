# Level B, amendment 1, S3: the full texts, both passes — the record

The plan is `../../plan_s3.md` (0c72c4a, committed before the fetch and the passes; its section 4's rule of a kept
paper's claims added before the merge, 0f3c508); the instructions `../../screening_fulltext.md` with
`../../screening.md` as amended on 2026-11-10; the tools `../../screening.py`, `../../independence.py`,
`../../kept.py` and `../../ft_record.py`. The journal's entries J11, J13, J14, J18, J20 and J21 say what happened, in
that order, on 2026-11-10 and 2026-11-11.

## The records and their texts

- **(a)** S2's 91 records (`../ta/fulltext.txt`): Europe PMC gave 45 full texts, the publishers 10 more and, asked
  again with a pace, 1 (`results/level-b-fulltext2` at 15c3c61, `results/level-b-publisher2` at a450cca; J18, J20).
- **(b), (c)** the ten F2 records and D's two (R00599, R01382), their texts the first fetch's (`results/level-b-fulltext`
  at e373c68); R00220, also among (a), read once, with the first fetch's text (J13).
- **(d)** none: all nine of A's records were left without a full text by the second fetch (J7).
- **Without a full text, 35** (`batches/without_fulltext.txt`): 15 refused by their publishers (bot challenges, a
  redirect page, a client challenge), R01281 (no DOI), and 19 bioRxiv preprints refused twice for the requests' rate;
  they go to the owner with the list, excluded by default (J1's rule).

## The batches and the sessions

67 records in 10 batches, the same for both passes, in three groups as their texts came (`../ft_group1.txt`,
`../ft_group2.txt`, `../ft_group3.txt`; batches 01–02, 03–09, 10; each group's texts' SHA256 in `batches/`), each
group's batches committed before their passes (26c55b5, 4e52e95, 3384411). Each batch was read once in each pass by a
fresh session with `plan_s3.md`'s task text word for word, one session at a time: group 1 by pass 1 then pass 2,
then groups 2 and 3 by pass 1 then pass 2; each pass's files out of the working tree while the other read (back
unchanged by `sha256sum -c`). Every file passed `screening.py check-fulltext`; none was read again, none stopped by
the usage limit. Each session's counts, from its file (`sessions.tsv`, `ft_record.py`):

| Pass | Batches | KEEP | UNCLEAR | EXCLUDE | Claim lines |
|---|---|---|---|---|---|
| 1 | 01–10 (67 records) | 5 | 4 | 58 | 16 |
| 2 | 01–10 (67 records) | 5 | 6 | 56 | 14 |

`decisions.sha256`: the 20 files as committed.

**Independence** (`independence.log`, `independence.py --calls-only level_b_protocol --lists`): no session's
transcript names the other pass's directory, the directories that held a pass's files while the other read, a batch
or a full text outside its own list, or the earlier decisions (`screening/`, `fulltext_screening`, `amendment1/ta`,
`claims.md`, `claims.json`, `decision_sheet`, `JOURNAL.md`, the plans); the protocol, which the instructions cite, is
named in no session's own calls (J14). 20 sessions, all independent.

## The comparison (`merge/`; `screening.py merge-fulltext`, the titles from `merge_inputs/all.tsv`)

- **Read by both:** 67 records; Cohen's kappa **0.831** (raw agreement 0.955).
- **The passes' decisions** (pass 1 / pass 2 → records): EXCLUDE / EXCLUDE 56, KEEP / KEEP 4, UNCLEAR / UNCLEAR 4,
  EXCLUDE / KEEP 1, EXCLUDE / UNCLEAR 1, KEEP / UNCLEAR 1.
- **Per pass** (`merge/ft_merge.json`): pass 1 F1 30, F2 1, F3 22, F5 2, F7 1, F10 2; pass 2 F1 31, F2 1, F3 21, F5 1,
  F10 2.
- **Kept by both, 4:** R00653, R00765, R01795, R02175 — every one with claims about sets only.
- **To the owner, 42** (`merge/ft_to_owner.txt`, both passes' notes and claims): the 7 decided differently or
  unclear in either (R00444, R00508, R01186, R01215, R02858, R03206, R03277) and the 35 without a full text.

`merge_inputs/all.tsv` is `../ta/all.tsv` with the rows of the 2026-10-19 screening (`../../screening/all.tsv`) for
the 11 records of the first group S2 did not read (the F2 records but R00220, and D's two): it gives the merge each
record's title; its counts "kept at titles and abstracts" count its rows, of both screenings.

## The kept papers and their claims (`../kept_claims.json`, `kept.py`; `plan_s3.md`, section 4)

8 papers: the 4 kept here (sets) and the 4 of 2026-10-19 (R00189, R01225, R01962, R02693; pairs), 34 claims — 15
pairs and 19 sets; 23 listed by both passes, 11 by one only (to the owner with the list); no pair whose symbols
cannot be told. **Fewer than 20 papers with pairs** (4; the protocol's section 9): said so with the list. The set
claims' genes are not yet taken (`plan_s5.md`, section 5). `../kept_claims_first_rule.json` is the same assembly
by the rule as first written, before J21's correction of its reading of parentheses (34 against 42 claims).
