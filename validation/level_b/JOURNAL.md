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

## J13 — 2026-11-10: S3's reading begins with the texts on hand (an order, not a rule)

While S3's fetch runs again (J11), both passes begin with the 12 records whose full texts the first fetch holds
(`results/level-b-fulltext` at e373c68, as text by `screening.py fulltext-text`): the ten F2 records and D's two
(`amendment1/ft_group1.txt`), in 2 batches (`amendment1/ft/batches/batch_01.list`, `batch_02.list`, 8 and 4
records; the texts' SHA256 in `texts_group1.sha256`). R00220, one of the ten, is also among S2's 91; it is read
here, with the first fetch's text, and not again. The batches of S2's other records follow, numbered from 03, when
their texts come. `plan_s3.md` has one list in the order of the records' numbers; here the order is in two groups.
Nothing else changes: the same instructions, task text, batch size, the same batches for both passes, one session
at a time, each pass's files out of the working tree while the other pass reads.

## J14 — 2026-11-10: S3's first group read by both passes; one precision of the independence check

**The four sessions** of J13's two batches ran one at a time (pass 1's two, then pass 2's two, pass 1's files out of
the working tree meanwhile and back unchanged by `sha256sum -c`); each file passed `screening.py check-fulltext`.
The decisions stay out of the repository until S3's last pass ends (`.git/info/exclude`), each pass's files out of
the working tree while the other pass reads, as in S2.

**The precision** (`independence.py --calls-only`, tested in `test_independence.py`): `plan_s3.md` forbids naming
the plans; the two instructions every session reads cite the protocol (`screening.md` line 4,
`screening_fulltext.md` line 3), so every transcript holds `level_b_protocol` twice. Counted over the whole
transcript, as the other forbidden paths are, it flagged all four sessions; a session now fails on the protocol
only where its own tool calls (the inputs of its requests) name it. Every other path is counted as before. On the
four sessions: the protocol named twice in each transcript and 0 times in its own calls; every other count 0; no
text named outside its own list; all four independent (the log goes into `amendment1/ft/independence.log` with
S3's other sessions).

## J15 — 2026-11-10: where R00189's pairs are published (for S5's spec; metadata only)

`plan_s5.md` (section 4) left R00189's spec waiting for where the pairing of GSE134809's ileal samples is published.
It was looked for here, before the spec, in the deposit's paper, Martin et al. 2019 (Cell 178:1493, PMID 31474370).
Europe PMC holds it as not open access (its full text answered HTTP 500, its supplementary files "not open
access"); PMC's page of the author manuscript was read once (it names Tables S1 and S2, and holds the methods);
PMC's download of the supplementary files asks the browser for a proof of work, which was not attempted; the
publisher's two files were fetched once each (`https://ars.els-cdn.com/content/image/1-s2.0-S0092867419308967-mmc1.pdf`,
SHA256 b20e5171ba5e36dbe08e680e84cd7c2b4bd4c48e460b08825d3b471204136792, and `…-mmc2.xlsx`, SHA256
030e5057afa1ebdfb821a197607ad3e462d598e626ec9483fa9f4f7d0c32dec5).
- **Table S1** gives the 14 patients' characteristics as a whole (no sample).
- **Table S2, sheet 1** gives each sample's number (`Sample_ID`, the number that ends its GEO title), tissue,
  status, patient (`Patient.ID`) and chemistry: 10 of the 11 pairs of ileal samples with their patient (pat. 5, 6, 7,
  8, 10, 11, 12, 14, 15, 16), and 9 PBMC samples. Samples 189 and 190 (GSM3972023, GSM3972024) are in no row; they
  appear in sheet 2 (each cluster's cells per sample) without a patient.
- Read: sheet 1; every sheet's name and first row; every sheet's first three columns, for 189 and 190; and the
  methods' text on how cells were called (barcodes with at least 800 counts from Cell Ranger's raw matrices, then
  mitochondrial, epithelial and hemoglobin filters). Nothing of the counts was read.

The spec (`data_specs/R00189.json`, batch 2 of `plan_s5.md`) follows: samples 189 and 190 are not paired by their
numbers (an inference, not the authors' statement); their cells are listed, not checked. `data_run.py`'s reader of
an Excel sheet (standard library; the run's environment has no Excel reader) gives sheet 1 cell for cell as
openpyxl does here. The spec was tried here on the real series matrix and Table S2 with fabricated matrices (no
count of the deposit): the 10 pairs found, samples 189 and 190 listed, both chemistries in both groups.

## J16 — 2026-11-10: R01962's object described again, in part; its spec (batch 2)

**The description** (workflow run 42897274410 on 26c55b5, the label `level-b-rdescribe-2`; `results/level-b-rdescribe-2`
at eccf6d8, every file as its SHA256SUMS): R read the object (`GSE266919_Myeloid.rds.gz`, 255,365,848 bytes, SHA256
54b63e7a…7d54) and wrote its cells' metadata whole (`meta_summary.tsv`: 56,180 cells; Barcode, Sample (78), Patient
(44), Group (Pre-treatment 32,995, Post-treatment 23,185), Tissue (breast, chest_wall, liver, lymph_node, lung),
Treatment (4), Efficacy (PR, SD, PD), Response (R, NR), n_genes, n_counts, Subset (20), SubCluster (5),
MajorCluster (Myeloid)), then stopped on its counts: the assay is a dense matrix, and `r_export.R` read a sparse
one's values (`m@x`), so `describe.json` (class, slots, assays, dimensions) was not written. **Not run again:** the
spec needs the metadata, which the description wrote; the dimensions and the counts' check come with the data run.
**Fixed** in `r_export.R`: a dense matrix is described by its first 1,000 cells (said so) and exported sparse 5,000
cells at a time; the export also writes a SingleCellExperiment's rowData and names its genes' first three in its
log. The data job now exports an R object before `data_run.py` reads it.

**The spec** (`data_specs/R01962.json`, batch 2 with R00189): the cells are the paper's own filter, read here in its
text (the TNBC data set's paragraph): breast biopsies, nab-paclitaxel with or without atezolizumab, responders,
patients with at least 30 myeloid cells before and after — 8 patients and 5,018 cells (2,215 before, 2,803 after)
in the paper. The deposit's columns give each (Tissue, Treatment, Response, Group, Patient); the minimum of cells is
`data_run.py`'s `min_cells_each_group`. No source is named: the deposit names no batch, platform or chemistry
beyond its samples.

**Tried here before the run:** R 4.3.3 with SingleCellExperiment 1.24.0 and Matrix 1.6.5, from Ubuntu's archive (the
job's versions), installed in this environment; on a fabricated object with a dense matrix of 11,420 cells, the
export in three blocks gives the dense values exactly, and `data_run.py` reads it, applies the cells' filter and the
minimum (a patient with 20 cells after is left out and listed), and names the genes by a rowData column where the
rownames are ids. Nothing of the deposit was read here beyond the description's metadata.

## J17 — 2026-11-10: a fault the dry run caught on 706f59e (pandas 2's missing values), fixed before any run used it

The data checks' dry job (workflow run 42902907844 on 706f59e) failed one test: in v1's environment (pandas 2.3.3)
`astype(str)` turns a missing value into the string "nan", so the cells of R00189's samples 189 and 190, which the
paper's table names no patient for, would have been checked as the cells of a donor "nan" (this environment's pandas
3.0.5, where the test had passed, keeps them missing). **Fixed** (`data_run.py`, `_labelled`): a missing value stays
missing in every column a spec names (the donor, the source, the sex label). **Batch 1 is not affected:** its reports
(`results/level-b-data-1` at 8b3a7c8) hold no value "nan" (R01225's sources include "unmatched", the regular
expression's misses, as intended). From now on level B's tests also run here in v1's environment (pandas 2.3.3,
numpy 2.5.3, as the jobs install it) before a push; batch 2's specs were tried again in it (R00189's on its real
metadata with fabricated matrices: 10 pairs, samples 189 and 190 listed; R01962's chain on a fabricated object).

## J18 — 2026-11-10: S3's fetch run again (J11) and judged; the second group's batches; bioRxiv's refusals asked again, paced

**The run** (workflow run 42897268539 on 26c55b5, the label `level-b-fulltext2` removed and added again;
`results/level-b-fulltext2` at 15c3c61, all 114 files as their SHA256SUMS): Europe PMC gave 45 of the 91 full texts
(the first run's 45) and answered HTTP 500 on every attempt for the first run's 35 records; J11's rule did not stop
it. Then one attempt at the publisher for the 45 records without a full text that have a DOI (R01281 has neither a
PMCID nor a DOI).

**Judged by `plan_s1.md`'s rule** (the texts by `screening.py publisher-text`):
- **A full text, 10:** the PDFs of R02036, R03510, R03573, R03635, R03742, R03866 (bioRxiv), R03132 (PLOS Biology),
  R03627, R03760 and R03916 (Research Square), each with its Results and Methods. Where the page holds the text too
  (R03132, R03627, R03916), the PDF is read, as for the others: the paper as published, its pages in order.
- **No full text, 15:** HTTP 403 with a 5.7 kB page, the bot challenge (R00242 and R00304, OUP; R00250, R00562 and
  R02185, Wiley; R00740, science.org); Elsevier's 2.7 kB redirect page (R00289, R00846, R01426, R01787, R02150,
  R02159, R03041); Nature's 3.0 kB "Client Challenge" (R00272, R01963). With R01281, these 16 go to the owner with
  the list, marked so, excluded by default (J1's rule).
- **Refused for the requests' rate, 20:** bioRxiv answered HTTP 429 (a 12 kB page) to 19 records' pages and to
  R03583's PDF (its page answered 200): R03468, R03474, R03475, R03515, R03520, R03557, R03569, R03583, R03593,
  R03633, R03642, R03644, R03733, R03739, R03794, R03857, R03874, R03915, R03952, R03953. **A deviation:** the
  one-attempt requests went out unpaced, and a 429 refuses their rate, not what they ask (bioRxiv's preprints are
  open). They are asked once more, by the label `level-b-publisher2` (`amendment1/publisher_again.tsv`), ten seconds
  between one-attempt requests (`search.py`'s `ONCE_PAUSE`, tested); R03583's page is asked again with its PDF. A
  record still without a full text then goes to the owner with the list, as the 16.

**The second group** (`amendment1/ft_group2.txt`, as J13's first): the 54 records with a full text now — Europe
PMC's 44 other than R00220 (read in the first group) and the publishers' 10 — in batches 03 to 09 (`ft/batches/`, 8
a batch, the last 6; the texts' SHA256 in `texts_group2.sha256`), the same for both passes. Pass 1 reads them first,
pass 2's files of the first group out of the working tree meanwhile. The records asked again follow as a third
group, its batches numbered from 10.

## J19 — 2026-11-10: S5's second batch (R00189, R01962)

**The run** (workflow run 42903541582 on b737f50, the label `level-b-data-2`; `results/level-b-data-2` at 1a8c6e1,
every file as its SHA256SUMS): both datasets read and checked.
- **R00189** (GSE134809; 607 MB of Cell Ranger's raw matrices and Martin et al.'s Table S2): counts non-negative
  integers; 90,894 ileal cells with at least 800 counts; **paired, 10 patients** in both groups (the table's pairs),
  rule 4 passes; the pair the table does not name (samples 189 and 190: 5,129 and 5,227 cells) listed, not checked;
  every patient one sex over both samples (7 male, 3 female; none split, none against a label — the deposit gives
  none); the claims' genes present; **the source balanced** (10x chemistry V1 and V2, each in both groups).
- **R01962** (GSE266919; the object exported by R: 21,035 genes by 56,180 cells, a dense matrix): counts
  non-negative integers; **the paper's stated filter gives 10 patients and 17,461 cells** (10,647 before, 6,814
  after; three patients left out for fewer than 30 cells in a group), **not the paper's 8 patients and 5,018 cells**
  (2,215 before, 2,803 after): the authors' cells are not identifiable exactly from the deposit, and the checks are
  made on the cells the stated filter gives (`plan_s5.md`, section 1). Paired, rule 4 passes; the Y genes are not
  among the deposit's genes (Xist in 35–56% of each patient's cells), so the sex check cannot tell one individual
  from two; the claims' genes present; **no source named** (the deposit names none beyond its samples). The filter
  is not changed to meet the paper's numbers.

## J20 — 2026-11-10: bioRxiv's refusals asked again, paced (J18): one full text; the third group

**The run** (workflow run 42913339645 on fe1f6f8, the label `level-b-publisher2`; `results/level-b-publisher2` at
a450cca, all 34 files as their SHA256SUMS): ten seconds between one-attempt requests, bioRxiv still answered HTTP
429 to 19 of the 20 records — 11 pages, and the PDFs of the 8 whose pages it answered (R03474, R03583, R03593,
R03644, R03733, R03739, R03857 and, with its PDF, R03874). **A full text, 1:** R03874's PDF (10,580 words, its Results
and Methods). bioRxiv refuses GitHub's runners whatever their pace, so **no third attempt**: the 19 (R03468, R03474,
R03475, R03515, R03520, R03557, R03569, R03583, R03593, R03633, R03642, R03644, R03733, R03739, R03794, R03857, R03915,
R03952, R03953) go to the owner with the list, marked "refused twice for the requests' rate", excluded by default as
the 16 of J18; another way to their texts (bioRxiv's own API, another network) is the owner's to decide, not taken
here.

**The third group** (`amendment1/ft_group3.txt`): R03874 alone, batch 10 (`ft/batches/batch_10.list`, its text's
SHA256 in `texts_group3.sha256`), the same for both passes, read after batch 09 in each pass.

## J21 — 2026-11-11: S3 complete — both passes, their comparison, the kept papers; one correction of kept.py

**The passes** (`amendment1/ft/README.md`): the 67 records with a full text in 10 batches, read by both passes, 20
sessions one at a time, every file passing its check, all 20 independent (`ft/independence.log`). Kappa 0.831
(agreement 0.955); kept by both 4 (R00653, R00765, R01795, R02175, claims about sets only); to the owner 42 (7
decided differently or unclear, 35 without a full text).

**The kept papers** (`amendment1/kept_claims.json`): 8 — those 4 and the 4 of 2026-10-19 — with 34 claims, 15 pairs
and 19 sets; 4 papers with pairs, fewer than the protocol's 20 (section 9: said so with the list).

**A correction of `kept.py`, after its first output and before any use of it.** Its rule (`plan_s3.md`, section 4,
0f3c508) matched the passes' claims by their genes as written; on the merge's claims it met forms it did not foresee:
a gene with its alias in parentheses ("HSPA5 (BIP)", "ERN1 (IRE1A)": one claim made two, and a "gene" no matrix
holds), a set's name followed by a description in parentheses ("Mitochondria NDCG (mitochondrial electron transport
chain)"), a comma inside those parentheses cutting the name, and an MSigDB set named inside a longer name
("HALLMARK_INTERFERON_GAMMA_RESPONSE weighted network"). Corrected (`kept.py`, tested): what parentheses hold is
dropped from a gene and from a set's name, a set's name ends at the first comma outside parentheses, and a set whose
name holds an MSigDB-style identifier is matched by it. The papers kept do not change, and no claim is dropped: 42
claims (17 pairs, 25 sets; 15 by both passes) became 34 (15 pairs, 19 sets; 23 by both). The first rule's output is
kept as `amendment1/kept_claims_first_rule.json`.

## J22 — 2026-11-11: the deposits of S3's kept papers looked at, for S5's batch 3 (names and listings only)

As J6 did for the first four, before any spec: GEO's file listings of GSE150728 (Wilk et al. 2020: per-sample R
objects, `cell.counts.matrices.rds.gz`, no cell types; a 390 MB archive), GSE153855 (per-sample count and RPKM tables
of Smart-seq2) and GSE214517 (the listing answered HTTP 404: looked at again with its spec); and
covid19cellatlas.org's page of patients' datasets, which links the Blish lab's annotated Seurat object
(`https://covid19.cog.sanger.ac.uk/submissions/release1/blish_awilk_covid_seurat.rds`, md5 52794c83…3d73d4c05 on the
page) — the object R01795 names. It is described first (`rdescribe3.txt`, label `level-b-rdescribe-3`; the job now
checks a file against the md5 its host names and unpacks only a gzip). Nothing of any count was read.
