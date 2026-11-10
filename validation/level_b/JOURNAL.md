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
