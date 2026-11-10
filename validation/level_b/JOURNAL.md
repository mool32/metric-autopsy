# Level B: the journal of decisions and deviations (after amendment 1)

Every decision of the owner on level B after the amendment of 2026-11-10 (`../plans/level_b_protocol.md`,
section 12), and every deviation from a plan committed before its run, in order. A deviation is recorded
when it is found, with what changed, why and what it touches; nothing is re-run to remove one.

## J1 — 2026-11-10: the owner's decisions on the sheet, and the way on

The owner, on `decision_sheet.md`: «A: R00014, R00173, R00394 — запроси Europe PMC ещё раз; по шести
остальным одна попытка открытого текста у издателя, иначе исключить. B: принимаю исключения. C: проверь
сейчас, доступны ли счётчики; нет — исключить. D: R00222 и R01395 исключить, R00599 и R01382 читают оба
прохода. Дальше иди без остановок: наборы (два прохода), независимые вердикты, структура данных, контроли.
Остановка одна — список с контролями на моё утверждение. Каждый шаг с планом до прогона, отклонения в
журнал. Дай оценку срока до первого прогона на 5a89f98.»

What follows from it:
- **A** (9 records without a full text): Europe PMC is asked again for R00014, R00173 and R00394 (each
  answered HTTP 500 to six attempts on 2026-10-19); R00438, R00582, R01062, R01338, R02068 and R02134 get
  one attempt at an open full text at their publishers. A record with a full text is read by both passes;
  one without is excluded.
- **B** (8 records read at the full text by one pass, which excluded each): excluded, as that pass decided.
- **C** (R03285, R03433, unclear in both passes on rule 3 alone): their repositories are checked once for
  public counts; a record with them is kept (both passes found the other rules held), one without is
  excluded.
- **D**: R00222 and R01395 excluded; R00599 and R01382 read at the full text by both passes.
- **The way on, without stops** (`PLAN.md`): the sets' two passes (and the preprints', amendment 1), the
  independent verdicts, the data's structure and origin, the controls; the one stop is the list with the
  controls, for the owner's approval. Every step has its plan committed before its run.
- The owner did not comment on section 12's text; the work proceeds under it as written.
- The sessions of the passes run one at a time (the owner, 2026-11-10: «Параллельных агентов не нужно»).

## J2 — 2026-11-10: rule 3's two wordings (found while writing `plan_s1.md`)

The protocol's rule 3 (section 5) reads "the data are public as counts (UMI or full-length; raw or the
authors' count matrices)"; `screening_fulltext.md`'s rule 3, which the passes of 2026-10-19 applied, reads
"raw reads or the authors' count matrices". Nothing is re-screened for it. S5 reads count matrices only (it
has no alignment of reads); a kept claim whose data are raw reads alone is listed so for the owner, and C's
check (`plan_s1.md`) records which of the two a repository holds.

## J3 — 2026-11-10: the preprints' run, and the full texts read by both passes (a deviation)

**The run** (S1, `plan_s1.md` section 2): workflow run 42852477112 on fc6ee3b, one attempt; `results/level-b-
preprints` at f616d1a, every file as its SHA256SUMS. Europe PMC gave 761 preprints for the query (one page),
all without a language. Merged with the corrected list here (`merge_preprints.py`;
`amendment1/preprints/preprints_merge.log`): 299 already in the list, 7 the same paper as an earlier preprint,
**455 new (R03507 to R03961)**; they are S2's, with the 210 sets' records (`plan_s2.md`).

**A deviation, fixed before S2's passes** (`plan_s2.md`, section 5): every record either pass keeps or marks
unclear at titles and abstracts is read at its full text by **both** passes. The protocol's section 5 has each
pass read the full texts of the records it kept and sends a record the passes decide differently to the owner;
under it, the first screening left 8 records read by one pass only and 4 read by none (the sheet's groups B and
D). The owner wants one stop, the list (J1), and decided D's two records this way («R00599 и R01382 читают оба
прохода»); a difference at titles and abstracts is settled by both passes reading the full text, and what stays
undecided after it goes to the owner with the list. It reads more full texts, never fewer.

## J4 — 2026-11-10: the positive control's files looked at before its loader was written

After `plan_s6.md` was committed (e857481) and before any selection, GSE96583's listing on GEO (the series' and
its two batch-2 samples' file names and sizes), the authors' table of cells for batch 2, its genes and the two
barcode lists were read here (in the scratchpad, as untrusted files) to write the loader, nothing else of the
data: the table has 29,065 cells (control 14,619, stimulated 14,446) with the columns `ind`, `stim`, `cluster`,
`cell`, `multiplets`; a stimulated cell whose barcode is also a control cell's carries "1" appended (313
cells), which maps every barcode of the two lists to the table. The counts were not fetched, and nothing in the
plan changed after it. The run fetches every file again in Actions, with its SHA256.

## J5 — 2026-11-10: the controls' script against its plan, before its run

`controls.py` implements `plan_s6.md` (e857481) with these precisions and one restriction, each fixed before
the run:
- **Primary tissue only** (`tissue_type == 'tissue'`, as the real pilot's case 3): an organoid or a culture
  of one "donor" measured in two conditions may be two differentiations or two passages, a difference of
  biology, which a negative control must not hold. The plan did not say it.
- **"tabula muris"** in a collection's name excludes Tabula Muris (2018) with Tabula Muris Senis: the 3-month
  mice of TMS are Tabula Muris's (`real_pilot_6_summary.md`), development data.
- **The directions' coin** is keyed by the seed 20261111 and the control's number, so that the three
  controls do not share one sequence; **the donors' sample** (only where more than 200 qualify) and the
  cells' sample use one generator with the seed 20261110, donors first (the plan said "the public seed, as
  v1" for the donors).
- **A donor's unit** in a condition is its cells of their commonest disease, sex and stage taken together
  (a donor's cells of another stage are another time point); the candidate's donors are those whose two
  units have at least 200 such cells and the same three values.
- **Reported, not rules:** each negative control's units' sex by the X-inactivation and Y genes (as Kang's),
  and the donors whose call differs between the two conditions.
- **Kang:** a lane's barcodes must all be in the authors' table and labelled with the lane's condition, and
  every name unique, else the structure check fails; the detection and the induction of the genes are
  measured on the kept cells (at most 400 per donor and condition), as the engine will see them.

## J6 — 2026-11-10: the four kept papers' deposits looked at, for S5's loaders (after `plan_s5.md`)

After `plan_s5.md` was committed (99be881) and before any of their counts was read, the GEO listings of the
four papers kept on 2026-10-19 were read here to plan S5's loaders (file names and sizes), and the barcode list of
GSE178265 (names only): GSE134809 (R00189) — per sample 10x MatrixMarket files, 31 samples, a 607 MB archive;
GSE189955 (R01225) — a count matrix and a table of cells as CSV; GSE266919 (R01962) — Seurat objects per cell
type (`Myeloid.rds.gz`, 244 MB), so R in Actions (`plan_s5.md`); GSE178265 (R02693) — one 4.7 GB MatrixMarket
file of 434,340 nuclei from 97 libraries whose names carry the donors' numbers, read in a stream. Nothing about a
claim's metric, groups or cells was read.
