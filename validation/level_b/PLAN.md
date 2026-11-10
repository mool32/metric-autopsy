# Level B after amendment 1: the steps to the list (the owner, 2026-11-10)

The protocol is `../plans/level_b_protocol.md` (approved 2026-10-19, amended 2026-11-10, section 12); the
decisions and the deviations are `JOURNAL.md`. The owner's way: no stops until the list with the controls,
which is the one stop, for the owner's approval; every step with its plan committed before its run (a run:
anything that fetches from the network, or a pass of the screening); deviations into the journal.

| Step | What | Its plan, before its run |
|---|---|---|
| S1 | Two runs in Actions: **fetch2** (Europe PMC asked again for A's three records; one attempt at an open full text at the publisher for A's six; C's two repositories checked for public counts) and **preprints** (amendment 12.2: Europe PMC's preprints, the section 3 terms and window, without `LANG:"eng"`, merged with the list); the set's rules added to `screening.md` and `screening_fulltext.md` (12.1) | `plan_s1.md` |
| S2 | Titles and abstracts: two passes (sessions one at a time, independent as before) over the 210 records that were X6 or X7 in either pass, and over the new preprints; the merge | `plan_s2.md` |
| S3 | The full texts of the records S2 keeps (a run, as the first fetch), with A's fetched records and D's R00599 and R01382: two passes; the kept papers, pairs and sets apart | `plan_s3.md` |
| S4 | The independent verdicts (section 7): a search of its own per kept claim (forward citations of the claim's paper; the claim's genes with the terms of section 3 and of perturbation), then each verdict coded with its citation and quote | `plan_s4.md` |
| S5 | The data of every kept claim, in Actions: counts raw integers; donor ids, each one individual (sex by Xist and the Y genes); the groups and their donors; the claim's genes present; the source of the cells shared or balanced between the groups (12.3's check of origin) | `plan_s5.md` |
| S6 | The controls (`real_claims.md`, "Controls"): the rules of selection registered first, then the negative controls (the same donors in two technical conditions) and the positive ones (perturbations of known response, such as Kang et al. 2018), each with its structure checked as in S5 | `plan_s6.md` |
| S7 | **The list with the controls to the owner** (the stop): the kept claims, pairs and sets apart, with the counts of section 1, their verdicts, their data's checks, 12.3's counts so far, the controls and the journal | — |

**Progress (2026-11-10):** every step's plan is committed (`plan_s1.md` to `plan_s6.md`, before their runs); S1's
preprints ran (J3), its second fetch is running; S2's passes run one session at a time; S6's selection script and
workflow are ready for their one run (`controls.py`, `level-b-controls.yml`); S5's checks are coded
(`data_checks.py`), its loaders follow the list.

After the owner's approval of the list: the runs on the engine of the tag `v0.4.0.dev0` (commit 5a89f98), with
(b) and (d) on every claim (`real_claims.md`, "The runs"), by a script and a plan committed before them.

**The estimate given to the owner on 2026-11-10** (not a commitment): 7 to 10 working days to the list, and 1
to 2 days after its approval to the first run on 5a89f98. The risks: the formats of the data on GEO (R objects
that Python does not read; a loader per dataset), the account's usage limit on the passes' sessions, the
number of preprints the search finds, and about 50 minutes of CI for every push.
