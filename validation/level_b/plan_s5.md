# Level B, step S5: the data of every kept claim — their structure and origin (the rules, before any data is read)

`../plans/level_b_protocol.md` (rules 3 to 5 of section 5; section 12.3's check of origin) and
`../plans/real_claims.md` ("The source in Tabula Muris Senis FACS"). These checks are made by code independent
of the engine, before any run of it, and nothing in them reads a claim's metric. Each dataset's loader is
written after the list of kept claims is known (S3), committed with its own plan line (which files, from where)
before it runs; the reads run in Actions, once per dataset, every file with its SHA256 on a results branch
(`results/level-b-data`), by a label of their own. A check that cannot be made is said so, with why. Every
deviation goes into `JOURNAL.md`.

## 1. What is read

- **The authors' deposit only:** the count matrices the paper names (GEO, ArrayExpress, CELLxGENE, the Human
  Cell Atlas, Zenodo, NGDC or another stable repository), with the cells' metadata deposited with them. No
  alignment of reads (J2: a claim whose data are raw reads alone is listed so, and not checked).
- **Formats:** MatrixMarket (with its barcodes and genes), 10x HDF5, h5ad, loom, text tables, and R objects
  (`.rds`, `.rda`): an R object is read by R in Actions (`SeuratObject` or `SingleCellExperiment`, their counts
  slot, exported to MatrixMarket with the cells' metadata), its versions recorded. A format none of these reads
  is said so.
- **The claim's cells:** the cell type or cluster the paper names, by the authors' own labels where deposited;
  where they are not, the claim is listed as "the authors' cells not identifiable from the deposit", and its
  checks below are made on every cell that is (no labels made here).

## 2. The checks, per claim

1. **Counts:** the matrix holds non-negative integers (raw counts or the authors' counts); a normalized or
   scaled matrix alone fails rule 3 (`raw_counts`, as v1).
2. **Donors:** each cell's donor (or animal), from the metadata or from the samples (a sample that is one donor,
   by the deposit's own description); the donors per group of the claim's comparison, and whether they reach
   rule 4 (4 against 4, 3 against 5, or 6 pairs). Ids that name several donors (a pool, "/") are listed.
3. **Each id one individual:** where both sexes are among a dataset's donors, each donor's sex called from the
   X-inactivation gene and the Y genes (`controls.sex_call`: the other sex's genes in at most 5% of the cells,
   its own in at least 10 times that share and 5%) is one and agrees with the metadata where it gives one; a
   donor whose cells split between the two calls, or whose call contradicts its label, fails.
4. **The claim's genes** present in the matrix (by symbol, or by the deposit's ids mapped by its own table).
5. **The source of the cells (12.3):** the condition is not confounded with the data's source. A source is a
   study, a platform or chemistry (10x 3' v2 against v3, cells against nuclei), or a processing or sequencing
   batch, as the deposit's metadata or the paper's methods name it — not a sample or a library, which holds one
   donor in one condition by design. Per group: **fails** where every cell of one group comes from a source
   that no cell of the other group comes from (the lesson of Tabula Muris Senis), **balanced** where every
   source holds both groups, and **partly** otherwise (listed with the sources' shares); where the deposit and
   the methods name no source beyond the samples, "no source named" (not a pass). A comparison of the same
   donors in two conditions processed together passes by construction.
6. **The design as the paper ran it**, recorded beside the checks: cells or pooled cells, its normalization,
   its statistic (the analysis (a) of `real_claims.md`).

## 3. What the list shows (S7)

Per claim: its data (where, which files, their SHA256), the six checks with their numbers, and the two counts
of 12.3 so far (the claims that fail the check of origin; the claims whose test over donors cannot run, by rule
4 or an unmeasured metric — whether it supports a claim is known only after the runs).

## 4. How the runs go (added 2026-11-10, before any run of S5)

- **In batches:** the specs of a batch are named in `data_specs/run<n>.txt`, each spec committed before the
  batch's run; the label `level-b-data-<n>` runs every one of them once, each in a job of its own
  (`.github/workflows/level-b-data.yml`), and their reports go to `results/level-b-data-<n>`. The deposits' files
  are not kept; their URLs, sizes and SHA256 are.
- **A deposit held only in an R object** (a Seurat or SingleCellExperiment object) is first described by a run of
  its own (its slots, the names and values of its cells' metadata), and its spec is written from that
  description; no check is made in that run.
- **Batch 1** (the four papers kept on 2026-10-19, as far as their deposits allow without R): R01225 (GSE189955,
  a counts table and a table of cells) and R02693 (GSE178265, a 4.7 GB MatrixMarket file read in a stream, the
  libraries' characteristics from GEO's series matrices). R01962's deposit is Seurat objects (GSE266919): it is
  described first. R00189's deposit (GSE134809) names no patient for its ileal samples (`JOURNAL.md`, J6): its
  spec waits for where the pairing of the samples is published, and the list says so if it is nowhere public.
- **Batch 2** (added 2026-11-10, before its run): R00189 (GSE134809). The pairing is published in the deposit's
  paper, Martin et al. 2019, Table S2, sheet 1 (each sample's patient; `JOURNAL.md`, J15), the publisher's file read
  in the run with its SHA256: it names 10 of the 11 pairs; the cells of the pair it does not name (samples 189 and
  190) are listed, not checked, and not paired by their numbers. The deposit holds Cell Ranger's raw matrices (every
  barcode): a cell is a barcode with at least 800 counts, the deposit's authors' first filter. The source: each
  sample's 10x chemistry (GEO). For it, `data_run.py` gains the threshold (`min_umis`), the paper's tables
  (`tables`: a CSV, or an Excel sheet read with the standard library) and, in every report, the cells no donor
  names; tested on fabricated deposits.
  **R01962** (GSE266919) joins batch 2, its spec written from its object's description (`level-b-rdescribe-2`; J16):
  the cells are the paper's own filter (breast biopsies; nab-paclitaxel with or without atezolizumab; responders;
  patients with at least 30 myeloid cells before and after: 8 patients and 5,018 cells in the paper), the donors
  its patients, no source named (the deposit names none beyond its samples). The data job exports an R object by
  `r_export.R` first (R from Ubuntu's archive; a dense matrix made sparse 5,000 cells at a time), and `data_run.py`
  gains the paper's minimum of cells per donor and group (`min_cells_each_group`, the donors under it listed) and an
  export's genes named by a rowData column where its rownames are not symbols; tried here on a fabricated object
  with R 4.3.3, SingleCellExperiment 1.24.0 and Matrix 1.6.5 (the versions of Ubuntu's archive the job installs).
- **Batch 3** (added 2026-11-11, before its runs): the papers S3 kept (`JOURNAL.md`, J21, J22). R01795 and R00653
  analyse one dataset, Wilk et al. 2020's blood cells (GSE150728), through the authors' annotated Seurat object that
  covid19cellatlas.org hosts (the object R01795 names; GSE150728's per-sample objects hold no cell types): described
  first (`level-b-rdescribe-3`, its host's md5 checked), then each paper's spec. R02175's islets (GSE153855,
  GSE214517, E-MTAB-5061) and R00765's two datasets (Morabito et al. 2021 on Synapse; Jerber et al. 2021 on Zenodo)
  follow, each with its own look at its deposit first.
  **Precised on 2026-11-11, after those looks (J22–J24) and before any spec of R02175 or R00765:** R02175's beta cells
  are the paper's clusters of its merged datasets (Conos), which no deposit holds: each of its three deposits is
  checked as it stands, its cells said not to be the paper's. GSE153855 (dataset 1; `tab_samples`): every cell of
  its 11 donors, the deposit naming no cell type; GSE214517 (dataset 2, new donors; `table_series`): every cell, its
  table RPKM only; E-MTAB-5061 (dataset 2, the published donors; `rpkm_counts_table`): the cells its authors inferred
  to be beta cells (the SDRF's "inferred cell type"), the nearest the deposits come to the paper's. R00765's claim 1
  on GSE174367 (`tenx_h5`; the oligodendrocytes, its table's "ODC"; groups its diagnosis, donors its samples, source
  its batch, sex its label); claim 2's dataset (Jerber et al. 2021) has no public counts (J23): listed, not checked.
  **Precised on 2026-11-11, after the object's description (J27) and before its specs:** R01795's and R00653's cells
  are those of GSE150728's 13 samples (Wilk et al. 2020's deposit) in the object covid19cellatlas.org hosts, by its
  `orig.ident`; the donors its `record_id`, the groups its `Status` (COVID, Healthy), the sex its `Sex`, the source its
  `Array.type`; the cell types its `cell.type.coarse` — R01795 "CD4 T" (spec `R01795`) and "B" (`R01795_b`, its
  claims 2 and 6), R00653 "CD14 mono" and "CD16 mono" (the paper's "monocytes"). The export keeps a spec's cells
  only (`r_export.R`, `keep.json`).

## 5. The sets' genes (added 2026-11-10, before any set's list is fetched)

A kept set's claim names where the paper lists its genes (`SET: <name>, <n> genes, <where>`; `plan_s3.md`,
section 4). Its genes are needed by S4's filter (a record that names two of them) and by S5's check 4, so they
are taken before both, the same way for every set:
- **The files:** for each kept paper with a set's claim, Europe PMC's supplementary files (its archive of the
  open-access paper's supplement, by the PMCID), fetched in Actions once, by the label `level-b-setgenes` (the
  branch `results/level-b-setgenes`), every file with its SHA256. A paper Europe PMC holds no supplement of is
  listed so, with what its text says of where the list is.
- **Which list is the set's:** per set, a line of `set_specs.json` naming the file, its sheet and column, and the
  rows (a module's label in another column, say) that give its genes — written after the files' names and first
  rows are looked at (as J6 looked at the deposits), committed before the genes are extracted; the extraction is a
  script's, its output each set's genes and their number against the number the paper states.
- **A set the paper names from a public collection without listing its genes** (added 2026-11-11, before the run;
  R01795's MSigDB Hallmark sets, which the paper took from MSigDB on 2022-10-20): its genes from the collection's
  file of the release current at the paper's date (MSigDB v2022.1.Hs, `h.all.v2022.1.Hs.symbols.gmt`, published
  2022-09-02), named beforehand in `amendment1/set_extra.tsv` and fetched in the same run with its SHA256.
- **A set the paper lists in its text, not in a table** (added 2026-11-11, before any set's genes are extracted;
  R00765's claim 2: "a set of 13 genes (COX5A, …)"): its genes are the list the sentence holds — the sentence quoted
  in `set_specs.json`, found by the script in the paper's text as S3 read it (that text's SHA256 checked against
  `amendment1/ft/batches/texts_group*.sha256`), the list cut between its brackets at its commas; nothing typed by
  hand. A set the paper shows only in a figure is "not listed", with where it is shown.
- **Another version of a kept record** (added 2026-11-11, before its fetch; `JOURNAL.md`, J24): where the version S3
  read lists no set's genes and the record joins another version that Europe PMC holds open access, that version's
  supplement is fetched the same way, once (`amendment1/setgenes<n>.tsv`, the label `level-b-setgenes-<n>`), and a
  set's spec names its archive (`rid:pmcid`). Its list is the set's only where it is the module the claim names;
  else the set is not listed. R00653: S3 read the preprint (PMC9774209), whose supplement's Table S13 gives the
  sub-modules' enrichment, not their genes; the journal article (PMC10415381) holds "Supplementary Data 1-4"
  (`setgenes2.tsv`).
- A set whose list is not found or cannot be read is said so: its claim stays in the list without its genes, S4
  searches it by its name only, and S5 reports its genes as not listed.
