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
GSE189955 (R01225) — a count matrix and a table of cells as CSV (Smart-seq); GSE266919 (R01962) — Seurat objects per cell
type (`Myeloid.rds.gz`, 244 MB), so R in Actions (`plan_s5.md`); GSE178265 (R02693) — one 4.7 GB MatrixMarket
file of 434,340 nuclei from 97 libraries whose names carry the donors' numbers, read in a stream. Nothing about a
claim's metric was read. GSE189955's table of cells (`RNA_meta_data.csv`, 5,329 cells) was read too: 22 patients
(`orig.ident`), the sites (`location`: PT 3,148 cells, FT 741, and seven others) and the authors' cell types;
whether a patient gives cells to both of the claim's groups is S5's check (`donors_check`), made in its run.
To write that dataset's line of S5, its sites were then tabulated against the authors' cell types and patients
(metadata only): the primary tumours' epithelial cells are 1,940 (16 patients) and the fallopian tubes'
secretory and epithelial cells 520 (6 other patients), the paper's two numbers; S5's run makes the checks.
GSE134809's series matrix was read too (the samples' titles and characteristics: tissue, "status: Involved" or
"Uninvolved", "10x chemistry: V1" or "V2", a patient id for the later samples only).
GSE178265's series matrices (per platform) were read too: per library its title (with the donor's number),
the tissue (substantia nigra pars compacta or caudate nucleus), sex, age, post-mortem interval and disease
("Ctrl" or "Disease").

## J7 — 2026-11-10: the second fetch, judged by `plan_s1.md`'s rules

**The run:** workflow run 42852476682 on fc6ee3b, one attempt; `results/level-b-fetch2` at 660d486. Every file
matches SHA256SUMS (`raw/fetch2/R03285-supplementary.zip.gz`, 48,457,509 bytes, is on the branch in two parts
of the publish step's split, joined before the check). The texts below are `screening.py publisher-text`'s.

**A — no full text for any of the nine; all nine excluded (J1):**
- R00014 (PMC8392141), R00173 (PMC7322959), R00394 (PMC8728743): Europe PMC answered HTTP 500 to all six
  attempts again, as on 2026-10-19. (NCBI's own services hold PMC articles too; the owner's decision named
  Europe PMC, so none was asked.)
- R00438 (OUP) HTTP 403, R01062 (Wiley) HTTP 403, R01338 (Wiley's literatum) HTTP 403: each a 5.7 kB page
  whose text is "Just a moment..." (a bot challenge). R00582 (Nature) and R02134 (Springer) HTTP 200: a 3.0 kB
  page, "Client Challenge … A required part of this site couldn't load". R02068 (Elsevier) HTTP 200: a 2.7 kB
  page, "Redirecting". None names a `citation_pdf_url`; none holds the paper's body.

**C — no public counts for either; both excluded (J1):**
- R03285: NCBI's links from PMID 42675628 to GEO: none. NGDC's page of PRJCA044651: "The data you are
  retrieving is not released in BioProject: PRJCA044651". Europe PMC's supplementary files (28): tables of
  83,900 cells' metadata (`s006.xlsx`: barcode, sample, group, counts per cell), of differential genes and
  gene lists, a clinical table, figures and PDFs — no count matrix.
- R03433: NCBI's links from PMID 42865737 to GEO: none. Europe PMC's supplementary files (10): `Table_1.xlsx`
  (cell types' proportions, differential genes, mitochondrial pathways, gene pairs' correlations per cluster)
  and figures — no count matrix.

With B (8 excluded) and D (R00222 and R01395 excluded; R00599 and R01382 to both passes in S3), the sheet's 23
records are decided: 2 go on to S3, 21 are excluded.

## J8 — 2026-11-10: the controls' selection (S6), and two points for the owner with the list

**The run:** workflow run 42866887133 on c7a9ab6, one attempt; `results/level-b-controls` at 9e20beb, every file as
its SHA256SUMS; the extractions in the run's artifact `level-b-controls-data` (90 days), their SHA256 and joinids
on the branch. The Census's `stable` release was 2025-11-08 (v1's). Counted: 93,203,254 human and 18,384,205 mouse
primary cells of tissues.

**Negative controls:** 26 qualifying candidates, 14 of them SEA-AD's (excluded); 6 other collections, the first 3
tried and chosen (each with raw counts and v1's pool of 8 pairs per level, 26 claims: 24 pairs and 2 sets; no
donor's sex call differs between its two conditions):
- **N1** — "A spatially resolved single cell genomic atlas of the adult human breast": breast fibroblasts, 10x 3'
  v3 cells against nuclei, 8 donors, 6,400 cells extracted.
- **N2** — "Early human lung immune cell development and its role in epithelial cell fate": lung, the cell type
  the Census labels **"unknown"**, 10x 3' v2 against 10x 5', 7 donors, 5,600 cells.
- **N3** — "Mapping the developing human immune system across organs": liver macrophages, 10x 3' v2 against 10x
  5' v1, 7 donors, 5,600 cells.

**The positive control (Kang et al. 2018):** its structure passes (29,065 cells, 24,679 singlets; counts integers;
the 8 donors in both conditions with CD14+ monocytes; no label problem; no sex call in conflict, every donor
female or unclear). All twelve genes are induced in all 8 donors, but **only 3 are detected in at least 10% of the
control monocytes** (ISG15 26%, IFI6 23%, ISG20 23%; OAS1 9.8%, MX1 6.9%, the others 2% to 5%): the engine's
domain rule leaves **one pair (ISG20–IFI6)** and no set (fewer than 4 genes).

**For the owner with the list (not decided here):** (1) N2's cell type is the Census's "unknown", unannotated cells
whose make-up may differ between the two chemistries (a difference of composition, not of technique); the rule
chose it, and the ranked table holds the next qualifying collections (rank 4: "Molecular Signatures of Resilience
to Alzheimer's Disease…", glutamatergic neurons, 10x 3' v3 nuclei against Drop-seq, 6 donors). (2) The positive
control gives one claim: with one positive dataset and one claim, l and r₁ rest on one pair; more positive
controls, or another cell type of Kang's with more basal expression of these genes, are the owner's to decide
(`plan_s6.md`, section 2).

## J9 — 2026-11-10: S2 done; S3's records

S2 ran as `plan_s2.md` says (`amendment1/ta/README.md`): 8 sessions, one at a time, all independent by their
transcripts; Cohen's kappa 0.861 over the 665 records; kept by both 68; 23 decided differently or unclear in
either; **91 records to the full texts** (`amendment1/ta/fulltext.txt`). No deviation.

S3's records (`plan_s3.md`, section 1): (a) the 91; (b) the ten F2 records; (c) R00599 and R01382; (d) none of A's
(J7). The fetch's list (`amendment1/fulltext_list.tsv`, `screening.py fulltext-list` with the new preprints) holds
(a); (b) and (c) are read from the first fetch's texts (`results/level-b-fulltext` at e373c68), where they are.

## J10 — 2026-11-10: S5's first batch, its specs before the run

`data_specs/R01225.json` and `R02693.json`, and `run1.txt`, written before the run (`plan_s5.md`, section 4). To
write them, besides what J6 records, the names of GSE189955's cells were tabulated against the sites (their
protocol after the sort gate: the fallopian tubes' 520 cells mostly `Split_FACS`, the tumours' 1,940 mostly
`Split_mouth` or `Smart_mouth`), which is why the spec takes that protocol as the cells' source; and GSE178265's
97 library names (each with its donor's number) and its series' columns (`disease`, `sex`, `facs classification`,
`10x version`, `status`). No count was read.
- (J10, continued.) R01962's deposit (GSE266919) is described first (`plan_s5.md`, section 4): `rdescribe1.txt`
  names its myeloid cells' Seurat object, `r_export.R describe` writes its slots and its cells' metadata, by the
  label `level-b-rdescribe-1` (`level-b-data.yml`); its spec is written from that description.

## J11 — 2026-11-10: S3's fetch failed technically; its rule of an outage changed; a second run

**The first run** (workflow run 42878734705 on afd99bf, the label `level-b-fulltext2`) stopped before publishing
anything: Europe PMC gave 45 of the 91 full texts, and answered HTTP 500 to all six attempts at 35 records (most of
them preprints, by their PPR ids: R00242, R00272, R00289, R00740, R01963, R02036, R03132, R03468, R03474, R03475,
R03510, R03515, R03520, R03557, R03569, R03573, R03583, R03593, R03627, R03633, R03635, R03642, R03644, R03733,
R03739, R03742, R03760, R03794, R03857, R03866, R03874, R03915, R03916, R03952, R03953). The first fetch's rule
(fb76f44) stops the run when more than a quarter of the records answer with server errors, taking that for an
outage; here 45 full texts came, so it was not one: Europe PMC answers HTTP 500 to some full texts it does not hold,
as it did to R00014, R00173 and R00394 (J7).

**The fix** (`search.py`, `outage`): the run stops for an outage only where more than a quarter of the records
answer with server errors **and** fewer than a quarter get a full text. A record refused on every attempt is
recorded as without a full text and goes on to the one attempt at its publisher, as before. The job's time limits
are longer (the first run took 75 of its 80 minutes, most of them in the six attempts at the 35). Nothing else
changes: the same list, the same requests.

**The second run**, by the label `level-b-fulltext2` removed and added again (a new run, its first attempt), on the
commit that carries the fix; its results are S3's.

## J12 — 2026-11-10: S5's first batch, and R01962's object described again

**The batch** (workflow run 42878735292 on afd99bf, the label `level-b-data-1`; `results/level-b-data-1` at
8b3a7c8, every file as its SHA256SUMS): both datasets read and checked.
- **R01225** (GSE189955): counts non-negative integers; nested, 16 patients' tumour cells (1,940) against 5
  patients' fallopian-tube cells (520), rule 4 passes; every patient female or unclear by the X and Y genes (no
  conflict); the claims' genes present; **the source partly**: the tumours' cells mostly by the protocol
  `Split_mouth` (77%) or `Smart_mouth` (17%), the tubes' mostly `Split_FACS` (83%), only `Split_mouth` shared.
- **R02693** (GSE178265): counts non-negative integers (5.1 GB read in a stream); nested, 10 donors with disease
  against 9 controls (393,873 nuclei of the substantia nigra), rule 4 passes; every donor's sex by the X and Y
  genes the one its library gives; the claims' genes present; **the source partly**: the nuclei's sort (the
  deposit's `facs classification`) Negative 66% and Positive 34% in the disease, Negative 55%, Positive 42% and
  NeuN 3% in the controls.

**R01962's description failed technically** (workflow run 42878734814, the label `level-b-rdescribe-1`;
`results/level-b-rdescribe-1` at 834edbb holds its log): the object is a SingleCellExperiment, whose package was not
installed, so R could not read it. Fixed: `r-bioc-singlecellexperiment` from Ubuntu's archive; the job now says
when a description fails. Described again by the label `level-b-rdescribe-2` (`rdescribe2.txt`, the same object).
