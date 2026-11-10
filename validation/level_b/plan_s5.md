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
