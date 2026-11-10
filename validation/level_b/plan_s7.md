# Level B, step S7: the list for the owner — what it holds and how it is put together (before it is assembled)

The owner (2026-11-10): «Остановка одна — список с контролями на моё утверждение». The list is the one stop: the
runs on the engine of the tag `v0.4.0.dev0` (commit 5a89f98) start only after the owner approves it. What each of
its parts shows was fixed by its step's plan (`plan_s3.md` section 4, `plan_s4.md` section 3, `plan_s5.md` section
3, `plan_s6.md` section 4) and by the protocol (`../plans/level_b_protocol.md`, sections 8, 9, 12.1 and 12.3); this
file gathers them, fixes the list's form and names what goes to the owner. The list is assembled by a script from
the repository's files and the results branches, every number in it from that script, with its log; nothing in it
reads a claim's metric. Committed before any of S3's kept papers, S4's verdicts or S5's later checks are known.

## 1. Its parts, in this order

1. **The search and the screening, counted** (the protocol's section 1): the records per source; titles and
   abstracts per pass and Cohen's kappa (the screening of 2026-10-19 and amendment 1's); full texts per pass and
   kappa; the papers kept; the papers with pairs and those with sets.
2. **The kept papers:** those both passes keep at the full texts, with the four of 2026-10-19. More than 30: the 30
   of the protocol's draw (section 9, `numpy.random.default_rng(20261020)` over the kept papers in the order of
   their records), the others listed. Fewer than 20 papers with pairs: said so (section 9), with the papers with
   sets counted.
3. **Each paper's claims, pairs and sets apart** (at most 5 a paper, section 8; both passes' claims and quotes): the
   genes (a set: its name, its number of genes, where the paper lists them), the species and cells, the
   comparison, the direction, the metric, the data, the donors and rule 4. **Every set's line carries 12.1's two
   caveats:** v1 does not cover a set's path (`set_norm_pearson` with the `coupling_set` probe), and the set's
   definition was not checked against Levy et al. 2020. The designs with fewer donors than the graded rule needs
   are listed apart.
4. **Each claim's data** (S5): the deposit and its files (URL, size, SHA256); the six checks with their numbers;
   "the authors' cells not identifiable from the deposit" where so; the cells left out because no donor is named,
   where so; a claim whose data cannot be checked, with why.
5. **Each claim's verdicts** (S4): the records found per source, filtered and read; each verdict with its kind,
   citation, quote and code; the records the two passes code differently.
6. **The controls** (S6): per control its source and version, its donors and cells per condition, its checks, its
   claims (genes and direction) and what was left out and why; and `plan_s6.md`'s check that no control's
   collection is a kept paper's (a control so found replaced by the next candidate of the ranked table).
7. **12.3's two counts so far**, as claims and as papers, with the reasons: the claims that **fail the check of
   origin** (a source that one group's cells do not share with the other's, an id that names two individuals,
   counts that are not raw integers); and the claims whose **test over donors cannot run** (fewer donors than the
   graded rule needs, or the metric not measured on the data). Whether the test supports a claim is known only
   after the runs, and is not guessed here.
8. **For the owner** (the owner decides; nothing here): the records the passes decide differently or mark unclear
   at the full texts, and those without a full text (the default, J1's rule: excluded); the verdicts the passes
   code differently; the journal's open points (J7: the three records Europe PMC refuses, which NCBI may hold; J8:
   the disease 'unknown' of the second negative control, and the one pair of the positive control; J15: R00189's
   samples 189 and 190); and every deviation the journal records since the owner's last decisions.

## 2. Its form

`amendment1/list/list.md` (the owner's page) and `list.json` (the same, for the runs' script), naming the commit
and the results branches it was assembled from, with their SHA256, and the journal's entries by number. The runs'
own plan (the engine of 5a89f98 on every kept claim, (b) and (d): `real_claims.md`, "The runs") is written after
the owner approves the list, not with it.
