# A third diagnostics, and v0.4's acceptance tests, on the development data (the plan, committed before the run)

The owner's decisions of 2026-10-18, items 3, 4 and 6 (`validation/plans/v0.4_spec.md`, "The owner's
decisions of 2026-10-18"). **Status: exploratory, on development data**, outside the confirmatory set
of the real track for good. One attempt; the outputs go to the branch `results/real-pilot-5`. Nothing
here is tuned after the run: every threshold, list, contrast and score is fixed below.

## What it is for

- **(b) Technical counts on §4.** In the males' cell types that are not fibroblasts, do the young mice
  carry more counts of the fibroblasts' matrix genes ("foreign" counts), and does the drop of MI go
  with them? Does a batch or a source in the metadata separate the 3-month mice from the old?
- **(a) The rule of the domain**, fixed here from principle, tried on two cases: on §4 in males only
  the fibroblast-like types should remain; in Marsh et al. 2022 the induced coupling must not drop out
  for a low detection in the group without the induction.
- **The positive control of coupling, Marsh et al. 2022** (GSE152183), registered here with the
  owner's conditions (item 4).
- **The acceptance tests of v0.4** (item 6), with the engine of the branch `claude/v0.4`.

## The data

- **Cases 1 and 2**: the second run's files, `results/real-pilot-2` at
  bad785a3c56e9de1ee24737953c1b4b63112aa7c, checked by that run's SHA256SUMS
  (`real_pilot_3_inputs.sh`): `case1.npz` (preprint §4's setup, 76,797 cells of 17 units, sha256
  `13e7a563…`), `case2.npz` (cardiac fibroblasts of 3 female and 4 male 3-month mice, 2,088 cells,
  `7191ad73…`) and `select.json`.
- **Case 3** (SEA-AD, a negative control by construction): the first run's files, `results/real-pilot`
  at 6133dfdf4c05f083f012aaec4136e376715af05a, checked by that run's SHA256SUMS: `case3.npz` (15,711
  nuclei of L2/3-6 intratelencephalic neurons of 16 donors, each by 10x 3' v3 and by 10x multiome,
  sha256 `31c297fe…`) and `select.json` with its 20 claims (3a1–3a8, 3b1–3b8, 3c1–3c4).
- **The full set of genes** for the cells of cases 1, 2 and 3, read again from the CELLxGENE Census,
  release 2025-11-08, by their `soma_joinid`, 10,000 cells at a time: every gene the Census marks as
  measured in every dataset of the case (21,025 for cases 1 and 2, 35,483 for case 3, as the earlier
  runs counted them), checked to be raw counts (`real_pilot_select.read_blocks`). Symbols that name
  more than one feature are left out. The counts of the case files' genes must equal those of the
  case files; a mismatch stops the run before any analysis.
- **The metadata of case 1** for (b4): every obs column of the Census's schema for case 1's cells;
  the Census's row of the dataset 98e5ea9f-16d6-47ec-a529-686e76515e39 (Tabula Muris Senis, FACS);
  the obs of that dataset's source h5ad (`cellxgene_census.download_source_h5ad`, the release's
  version; its matrix is never read); and GEO's records (SOFT, family) of GSE132042 (Tabula Muris
  Senis) and GSE109774 (Tabula Muris, 2018).
- **Marsh et al. 2022**: the 12 sample matrices of GSE152183 from GEO (`GSE152183_RAW.tar`, 76,861,440
  bytes by GEO's file list), each file checked against the file list's names and sizes; the
  metadata of `results/marsh-metadata` at 757c2ab48af792f63b1992c8eb98b05a17d4272f.
- **The engine**: v0.4.0.dev0 at the head of the pull request of `claude/v0.4` when the run starts,
  installed in v1's analysis environment (`validation/prereg/requirements-panel.txt`); the run
  records the commit and the git tree of `src/`. The scripts: `real_pilot_5_inputs.sh`,
  `real_pilot_5_extract.py` (v1's selection environment) and `real_pilot_5.py` (the analyses).

## Marsh et al. 2022: the registration

Checked before this plan, on the metadata only (`results/marsh-metadata`, workflow run 39110775480):
GEO's records of GSE152184 and its subseries, Europe PMC's full text (CC BY) and the supplementary
archive. No count has been read.

- **The samples** (GSE152183's SOFT records and Supplementary Table 1): 12 male C57BL/6J mice at P90;
  whole brain; CD45+/CD11b+ cells sorted by FACS; 10x 3' v2; Cell Ranger 3.0.0 on mm10. Four arms of
  3 mice: enzymatic digestion without inhibitors (ENZ-NONE: GSM4605688, GSM4605692, GSM4605696;
  GEO's "Miltenyi_NO"), enzymatic with the inhibitor cocktail (ENZ-INHIB: GSM4605690, GSM4605694,
  GSM4605698; "Miltenyi_IN"), cold Dounce without (DNC-NONE: GSM4605687, GSM4605691, GSM4605695) and
  with (DNC-INHIB: GSM4605689, GSM4605693, GSM4605697).
- **The mice**: GEO's sample names reuse the numbers 01–06 across the two methods (Dounce_NO_sample_01
  and Miltenyi_IN_sample_01), but Supplementary Table 1 numbers the samples 1 to 12 and the text has
  "19,563 cells from n = 12 mice (n = 3 per group)": each GSM is one mouse, and the replicate unit is
  the GSM.
- **Batches** (the records' `batch`; Supplementary Table 1): batch 1 holds ENZ-NONE 2, ENZ-INHIB 1,
  DNC-NONE 2, DNC-INHIB 1; batch 2 the rest (1, 2, 1, 2). Recorded; the primary analysis is not
  stratified by batch, and each mouse's value is listed with its batch.
- **The matrices**: per sample, Cell Ranger's barcodes, features and matrix (12 to 15 KB of barcodes,
  so the called cells, not every droplet). Every barcode is kept, without a further filter; the
  barcodes per sample are reported beside Supplementary Table 1's post-QC numbers (991 to 2,209).
- **The gene set, from the paper**: Supplementary Table 4, "Micro/Myeloid Shared Act. Score", 25 genes:
  Rgs1, Hist2h2aa1, Hist1h4i, Nfkbiz, Klf2, Junb, Dusp1, Ccl3, Hspa1a, Hsp90aa1, Fos, Hspa1b, Jun, Jund,
  Nfkbid, Gem, Ccl4, Ier5, Txnip, Hist1h2bc, Zfp36, Hist1h1c, Egr1, Atf3, Rhob. **The claim's pair is
  Fos and Jun**: c-Fos and c-Jun, the AP-1 pair, the immediate early genes the paper names first for the
  ex vivo activated microglia ("several immediate early genes (IEGs) (for example, Fos, Jun)"), both in
  the set and in the all-CNS activation score (its rows 1 and 4). Chosen by this reason, before any count.
- **The data for the engine**: the cells of the two arms of a contrast; D1's cut over all 12 samples'
  cells (Fos and Jun, then the most expressed by mean count, 2,000 genes); `cell_type` "CD45+/CD11b+
  sorted cells" for every cell (GEO gives no annotation, so the domain has one type); the full matrix
  as the pool of the check of specificity.
- **The contrasts** (the owner's item 4), each for the coupling of Fos and Jun:
  - **M1, primary**: ENZ-INHIB → ENZ-NONE, `increase` (the same exposure, the induction blocked in one
    arm);
  - **M2, secondary**: DNC-NONE → ENZ-NONE, `increase`;
  - **M3, null**: DNC-NONE → DNC-INHIB, `two-sided`.
  The metric: `norm_pearson` (primary) and `mi_3bin` (secondary, reported). The SESOI: 0.1 on the
  construct's scale (a change of 0.1 in the correlation of the pair's log expected expression), the
  probe `injected_signal.coupling("Fos", "Jun")`, dosed by the engine. No GATE 5 controls.
- **3 mice per arm: the test is parametric only.** The exact permutation over 3 against 3 mice has
  C(6, 3) = 20 relabelings, its smallest two-sided p 0.1, so it cannot reach α = 0.05; the engine's
  replicate-level t decides and the verdict carries *parametric only*. Recorded here, before the run.
- **The caveat, in the report as it is**: the genes' mean and detection move together with their
  coupling. The control checks that the engine does not refuse a real change of coupling, not that it
  tells a change of coupling from a change of detection.

## (b) Technical counts on §4

The males of case 1: 56,754 cells, 9 units (young 3_8_M, 3_9_M, 3_10_M, 3_11_M and the composite
3_10_M/3_11_M; old 24_58_M to 24_61_M), as the second run. The types are the second diagnostics'
(check (v)): the Census's `cell_type` with at least 20 cells in at least 3 units of each group
(35 types in males); **fibroblast-like** where the label contains "fibroblast", "mesenchymal" or
"stromal" (any case), the rule of the second diagnostics; **non-fibroblast** the others. A unit
counts for a type with at least 20 of its cells. Every difference below is young minus old.

**The genes, fixed here from the literature:**
- the matrix genes that Muhl et al. 2020 (Nat Commun 11:3953, doi:10.1038/s41467-020-17740-1; the full
  text, CC BY, fetched by the metadata workflow) name in their short-list of markers common to the
  fibroblasts of four organs: "many ECM genes, such as Col1a1, Col1a2, Col5a1, Loxl1, Lum, Fbln1, and
  Fbln2, as well as the cell surface receptors Cd34 and Pdgfra" (the two receptors are left out: not
  matrix genes, and Cd34 marks endothelial and haematopoietic cells too);
- the owner's examples not among them: Col3a1, Dcn, Sparc.

The list: **Col1a1, Col1a2, Col3a1, Col5a1, Dcn, Fbln1, Fbln2, Loxl1, Lum, Sparc.**

- **(b1) Detection and level, per gene.** Per (unit, type): the share of the cells with a count of the
  gene, and its mean log1p CP10k with the cell's total over the full set of genes. Per type, Δ_type
  is the young units' mean minus the old units' (each unit weighted equally); **combined**, the mean
  of Δ_type over the non-fibroblast types with equal weights. Its exact p is over all 126 relabelings
  of the 9 units into 5 and 4, each type's Δ recomputed from its own units (a type left with an empty
  group drops out of that relabeling's mean), two-sided (|Δ| at least as large). The same over the
  fibroblast-like types, reported.
- **(b2) The index of foreign counts, per cell:** F = log1p(10⁴ × the cell's counts of the index genes
  / the cell's total over the full set of genes). **The index genes are the list less Col1a1 and
  Sparc**: Col1a1 is a gene of the claim, and an index that holds it is tied to the claim's metric by
  construction; Sparc is expressed outside the fibroblasts by name in the literature (Marsh et al. 2022
  list it in the microglial identity score, their Supplementary Table 4). F⁺, with Col1a1, is
  reported. A unit's F in a type is the mean over its cells of the type; combined as (b1).
- **(b3) The association with MI.** For each (unit, type) of the non-fibroblast types: MI =
  `mi_3bin(Smad3, Col1a1)` on the unit's cells of the type (the second diagnostics' within-type MI),
  and F.
  - **(b3-i), scored: within type and age.** r, Pearson's correlation of MI and F after the mean of
    each (type, age) cell is subtracted from both; p from 9,999 permutations of F among the units of
    each (type, age) cell (`numpy.random.default_rng(20261018)`), (1 + #|r*| ≥ |r|) / 10,000. It asks
    whether, among mice of one age within one type, more foreign counts go with a higher MI. An
    association centred by type alone would follow from (b2)'s difference and the known drop of MI
    by themselves, so it is reported (b3-0) and not scored.
  - **(b3-ii), reported: across types.** r between the types' Δ of MI and Δ of F; p from 9,999
    permutations of the types' Δ of F.
- **(b4) The metadata.** For each obs column of the Census, and of the source h5ad, that is neither the
  age, nor the id of a unit, nor the cell type or tissue: whether it separates the 3-month units from
  the old in each sex (no value shared by a young and an old unit; numeric columns, such as the
  Census's `raw_sum` and `nnz`, as per-unit medians, their young − old difference with the exact p over
  the relabelings). A column with one value per unit and one unit per value is an id of the unit and is
  left out; a column whose values never cross the ages but each belong to one unit is "nested within
  units", not a batch. The units' sex and age are case 1's. GEO: whether the records of GSE109774 name
  the units of case 1, each id matched as a whole token with "_", "-" or "." between its parts, also
  with its sex before its number ("3-M-8"); and the lines of GSE132042's summary and overall design
  that mention the 3-month stage or Tabula Muris, quoted.

**The prediction (the owner's, 2026-10-18, verbatim):**

> у молодых в нефибробластных типах больше счётов генов матрикса, и разница MI идёт вместе с ними

("In the non-fibroblast types the young have more counts of the matrix genes, and the difference of
MI goes together with them.") Scored: its first part holds where (b2)'s combined Δ of F is positive
with p < 0.05; its second where (b3-i)'s r is positive with p < 0.05. The prediction holds where both
do, partly where one does.

**The owner's statement to check (verbatim):**

> Насколько я знаю, трёхмесячные взяты из исходного Tabula Muris — проверь.

("As far as I know, the 3-month mice are taken from the original Tabula Muris — check.") Scored by
GEO: it holds where GSE109774's records name every non-composite 3-month unit of case 1 (3_8_M, 3_9_M,
3_10_M, 3_11_M, 3_38_F, 3_39_F, 3_56_F, 3_57_F) and no old unit; it does not hold where they name
units of case 1 but miss a 3-month one or name an old one; undecided where they name none. (b4)'s
columns say whether the data carry a batch or a source that separates the ages.

## (a) The rule of the domain

**The rule** (`metric_autopsy.domain`, fixed before this run from principle, not from §4's data): a
cell type is in the domain of a claim about a pair where, **in at least one of the two groups**, both
genes are detected in at least **10%** of the type's cells; detection in a group is the mean over its
replicates, each weighted equally, of the share of the replicate's cells of the type with a count of
the gene. 10% is the lower bound of validation v1's expression levels (low: detected in 10–50% of
cells): no metric was tested below it. *In at least one group*, because a claim of a change may be the
coupling's very appearance, so that requiring both groups would remove the type where the change
itself moves a gene across the threshold, and pooling the cells would make the rule depend on the
groups' numbers of cells. For a set of genes, the set (a count of any of its genes) is detected so;
a claim without genes takes every type. A type also needs at least 20 cells in at least 3 replicates
of each group. The rule applies within each stratum. The detections per type on §4 have never been
computed, in this or an earlier run.

- **(a1) §4 in males**: the domain of `mi_3bin(Smad3, Col1a1)`, young → old, on the males of case 1,
  with every type's detections per group.
- **(a2) Marsh**: the domains of M1, M2 and M3, with the detections of Fos and Jun per arm, and, for
  the record and not used, whether the type would be in M1's domain under "in both groups" and under
  "over the pooled cells".

**The predictions (the owner's, 2026-10-18, verbatim):**

> В §4 по самцам должен остаться только фибробластоподобный домен.

Scored: holds where the domain of (a1) is not empty and each of its types is fibroblast-like by (b)'s
label rule.

> В Marsh индуцированная связь не должна выпасть из-за низкой детекции в неиндуцированной группе.

Scored: holds where M1's one type is in its domain.

## The acceptance tests of v0.4

`run_autopsy` of `claude/v0.4` on each claim, as a user of the API runs it: the composition estimand;
the claim's direction; the SESOI on the construct's scale with the claim's probe, dosed by the engine;
the judgment gates resolved (`judgment_pending` False, as in v1's panel); seed 0; the run log off;
`stop_on_first_fail` by default; the check of specificity on, its pool the full set of genes of the
same cells; no second dataset (so SUPPORTED is provisional). The claims:

| id | claim | data | settings beside the common ones |
|---|---|---|---|
| A1 | `log_cp10k_mean(Xist)`, female → male, decrease, SESOI 0.5 | case 2 | probe `module(["Xist"])`; composition `cell_type`; `known_to_differ` the first diagnostics' 23 sex genes |
| A2 | `log_cp10k_mean` of Ddx3y, Eif2s3y, Kdm5d, Uty summed, female → male, increase, SESOI 0.5 | case 2 | probe `module` of the four; as A1 |
| A3 | `mi_3bin(Smad3, Col1a1)`, young → old, decrease, SESOI 0.1 | the males of case 1 | probe `coupling`; controls Actb–Gapdh and Tst–Lrrc42 (as the pilot); composition `cell_type`, `tissue` |
| A4 | the same | case 1, within sex | as A3 |
| A5 | `norm_pearson(Fos, Jun)`, M1 | Marsh | as registered |
| A6 | `norm_pearson(Fos, Jun)`, M2 | Marsh | as registered |
| A7 | `norm_pearson(Fos, Jun)`, M3 | Marsh | as registered |
| A5m–A7m | `mi_3bin(Fos, Jun)`, M1–M3 | Marsh | reported, not scored |
| A8 | the 20 claims of case 3, two-sided, SESOI 0.1 (pairs) and 0.5 (means) | case 3 | probes as the first run's (`coupling`; `module` of the gene); paired by donor; composition `cell_type` |

**The criteria (the owner's, 2026-10-18, verbatim)**, each scored as below; a verdict is SUPPORTED
where its label is SUPPORTED (provisional or replicated):

> Xist и Y-гены — SUPPORTED;

A1 and A2 SUPPORTED.

> §4 по самцам — не SUPPORTED: внутри домена эффекта нет, а отчёт показывает, что падение в смеси идёт из типов вне домена;

A3: not SUPPORTED; the domain is not empty and the effect within it is not detected (the effect
field's status is not DETECTED); and the report's mixture and its types outside the domain both show
the decrease, each detected in the claimed direction by the engine's test.

> §4 по обеим стратам — не SUPPORTED;

A4 not SUPPORTED.

> Marsh: основной контраст — SUPPORTED, нулевой — не SUPPORTED;

A5 SUPPORTED, A7 not SUPPORTED. A6 and the `mi_3bin` claims are reported.

> SEA-AD (технический отрицательный контроль) — ни одного SUPPORTED, только EFFECT BELOW SESOI или отсутствие эффекта.

A8, two readings: **strict**, each of the 20 verdicts is EFFECT BELOW SESOI or NO DETECTABLE EFFECT;
**lenient**, none is SUPPORTED. Both are reported.

## Known before the run

From the earlier runs' logs, so not outcomes of this one:
- **§4 in males** (`real_pilot_4.log`): MI old minus young −0.01975; within the 35 types −0.01356
  (exact p 0.01587); the four fibroblast-like types +0.008123 on average (2 of 4 lower), the other 31
  −0.01636 (29 of 31 lower); the composition of types differs (TVD 0.3411, p 0.03175); 20 of 1,000
  matched control pairs at least as far on the donor-level t.
- **The sex claims** (`real_pilot_2.log`): the test over mice found −0.2215 (Xist) and +0.6700 (the Y
  genes), p 0.02857, the smallest attainable with 3 against 4 mice.
- **SEA-AD** (`real_pilot.log`): v0.3's engine gave SUPPORTED on 3a2 (`norm_pearson` of SYN3 and
  PDZD2, +0.0239) and 3c3 (RBFOX1, +0.1213), whose own TOST called them equivalent to zero within
  the SESOI (the evidence of change 14), and NO DETECTABLE EFFECT or INCONCLUSIVE on the other 18.
- **Marsh**: only the paper's own results (the ex vivo activated cluster, nearly all from ENZ-NONE);
  no count of GSE152183 has been read.

## One attempt

`.github/workflows/real-pilot-5.yml`, started by the label `pilot-run-5` on the pull request of
`claude/v0.4`, commits whatever the analyses produce, failures included, to the new branch
`results/real-pilot-5`; a second start finds the branch and stops. If the preparation or an
extraction (the Census, the source h5ad, GEO) fails, no analysis starts and nothing is committed: it
may be fixed and started again, and the report says so. The analyses run in nine parallel jobs: (b);
(a); A1–A2; A3; A4; A5–A7 with A5m–A7m; A8 in three parts (3a, 3b, 3c). Each job has one attempt. The
time and the memory of the heaviest job, A4, were measured before the run on fabricated counts of
case 1's size (`real_pilot_5_timing.py` and its log; nothing real).

## What is reported

`real_pilot_5.log` and `real_pilot_5.json`: for (b), the tables per gene and type, the index, the
associations, the metadata's columns and GEO's matches, and the prediction's and the statement's
scores; for (a), the domains with their detections and the two predictions' scores; for the
acceptance tests, each claim's verdict, cause, fields, the domain, the mixture and the outside
estimates, s_m and the dose, the line of specificity, and each criterion's score; the hashes of the
inputs and the engine's commit. One report to the owner.
