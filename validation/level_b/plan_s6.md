# Level B, step S6: the controls — the rules of their selection, registered before any is selected

`real_claims.md`, "Controls: the truth by construction": the controls measure the second question's exchange on
real data (π\*_real); a control dataset is one unit, its claims taken together. They are run only if the first
question passes its threshold, but they are selected now, so that the owner approves them with the list (the
owner, 2026-11-10: «Остановка одна — список с контролями на моё утверждение»). Nothing below is selected, fetched
or run before this file is committed. The selection runs in Actions, once, by a label of its own; its log, its
tables and the SHA256 of everything it extracts go to a results branch. Every deviation goes into `JOURNAL.md`.
The engine is not run here: these rules choose data and claims, and nothing in them reads a claim's metric.

## 1. Negative controls: the same donors in two technical conditions

**The source.** The CELLxGENE Census, its `stable` release at the run (recorded), human and mouse, primary data
(`is_primary_data`), read as v1's backgrounds were selected (`validation/prereg/select_backgrounds.py`, imported
unchanged, in its pinned environment `requirements-select.txt`).

**A technical condition** is a cell's assay and suspension type (`assay`, `suspension_type`): 10x 3' v2 against
v3, cells against nuclei, a droplet method against a plate one. **A candidate** is a collection (the Census's
table of datasets gives each dataset's collection), an organism, a tissue, a cell type and two technical
conditions A and B, A the first of the two in the order of their labels. A donor is a `donor_id` within the
collection: the same string in two datasets of one collection is one donor. A candidate qualifies where:
- **at least 6 donors have at least 200 cells in each of A and B** — 200 is v1's cells per donor
  (`panel.CELLS_PER_DONOR`, which v1's rule of pairs below needs); 6 pairs is the smallest number at which the
  engine's exact paired test can reach p < 0.05 (2/2^6 = 0.031; at 5 pairs its smallest p is 0.0625);
- **the conditions differ in the technique only:** for each such donor, the commonest disease, sex and
  development stage of its cells are the same in A and in B (the tissue and the cell type are the candidate's);
- **raw counts:** non-negative integers on the extraction (v1's `raw_counts`);
- **v1's rule of pairs gives a pool** on the extraction (below);
- **not data the engine was developed or validated on, and not the positive control:** not the collections of
  v1's backgrounds B1 and B2 (`results/panel-v1`), of the real pilot (Tabula Muris Senis; SEA-AD; Marsh et al.
  2022) or of Kang et al. 2018, matched by the collection's DOI (the DOIs are listed in the script); and not a
  collection whose DOI is that of a paper kept in level B (checked at S7, when the list is final: a control so
  found is replaced by the next candidate that qualifies, from the ranked table the run writes whole).

**The choice.** The candidates ranked by their donors that qualify (more first), then by the median of their
cells per donor and condition (more first), then by the collection's id, the tissue, the cell type and the
conditions' labels in order; one candidate per collection (its first). The raw counts and the pool are checked
on the extraction, in this order, over at most the first 10 collections; **the first three that pass are the
negative controls.** Fewer than three: the list says so, and nothing is put in their place by judgement.

**The extraction**, per candidate tried: every qualifying donor (at most 200, a sample with the public seed, as
v1), at most 400 of its cells in each condition (v1's cap; a sample of its sorted joinids with the public seed
20261110), all genes, the raw counts; with each cell its donor, assay, suspension type, tissue, cell type,
disease, sex, development stage and dataset. Written as v1's npz, with its SHA256 and the cells' joinids, so
that the same Census release gives the same cells again.

**The claims of a negative control** (by rule, on the extraction, before any metric of a claim is computed):
- **Pairs: v1's rule** (`panel.plan_background`, its code unchanged) with each donor in each condition as v1's
  donor: genes neither mitochondrial nor ribosomal, detected in at least 5% of the cells of every donor in each
  condition, in v1's three levels of expression; per level, the k most correlated disjoint pairs by the pooled
  within-unit correlation of log-normalized expression with the log library size partialled out (200 cells per
  unit, v1's public seed), k = 8 or the largest k ≥ 4 that every level has. Every pair of the pool is a claim
  (at most 24). v1's plan gives each its GATE 5 controls (a negative pair; one positive pair per control),
  which the runs use if their plan gives every claim GATE 5's controls.
- **Sets:** per level, the genes of that level's pairs form one set (2k genes), for the high and the medium
  levels only — a low level's genes are detected in 10% to 50% of the cells, so most of the set's pairs would
  have fewer than the 20 cells with both genes that `norm_pearson`, and so `set_norm_pearson`, needs. At most
  2 sets.
- **Directions:** each claim says the coexpression is higher in A, or higher in B, by a fair coin with the
  seed 20261111, in the order of the pool's pairs, then the sets (high, medium).
- **The truth:** no difference of biology between A and B in one donor; an effect claimed in either direction
  is false (`real_claims.md`). Cells and nuclei hold different shares of a gene's transcripts: a difference of
  coexpression between them is the technique's, which is what a negative control is for.

## 2. Positive control: interferon-β (Kang et al. 2018)

**Named, as v1 named B3 and B4** (`real_claims.md`): Kang et al. 2018 (Nat Biotechnol 36:89–94), GEO GSE96583,
its second batch: peripheral blood mononuclear cells of 8 lupus patients pooled, half the pool stimulated with
interferon-β for 6 hours, the donors demultiplexed by their genotypes; the donor is the pair. Fetched once in
Actions from GEO (the batch's two count matrices, its genes, its barcodes, the authors' table of cells), every
file with its SHA256.

**Its structure is checked first** (S5's checks): the authors' singlets only; the donors in both conditions and
their cells per condition in CD14+ monocytes (the authors' label); the counts non-negative integers; each
donor's sex from XIST and the Y genes the same in both conditions. Fewer than 6 donors in both conditions: the
control fails, and the list says so. At most 400 cells per donor and condition (the public seed 20261110), as
the negative controls.

**Its claims** (in CD14+ monocytes):
- **The genes, fixed now from the literature and not from the data:** twelve canonical interferon-stimulated
  genes — ISG15, IFI6, IFIT1, IFIT2, IFIT3, MX1, MX2, OAS1, OAS2, RSAD2, IFI44L, ISG20.
- **Eligible:** a gene in the data, detected in at least 10% of the cells in each condition (the engine's
  domain rule for a pair: the mean over the donors of the share of their cells with a count), and induced in
  every donor (its mean log-normalized expression higher in the stimulated cells than in the donor's control
  cells): the induction is the truth these claims rest on, so it is checked, gene by gene, on the level of
  expression, which is not the claims' metric.
- **Pairs:** a random matching of the eligible genes with the seed 20261112 (disjoint pairs; with an odd number
  one gene is left unmatched): at most 6 claims. **A set:** the eligible genes, where there are at least 4.
- **Direction:** the coexpression is higher in the stimulated cells. **The truth:** the genes' induction by
  interferon-β, known and checked above.

**For the owner, with the list** (not decided here):
- The owner's rule of 2026-10-14 (a positive control's effect clearly above the SESOI) cannot be checked on
  these claims before the runs without measuring them; this plan meets it through the perturbation (a strong
  induction, checked on the genes' level), not through their coexpression.
- One positive control dataset gives no interval over positive datasets. Candidates for more, none verified
  here (neither their data nor these citations): the 1M-scBloodNL data (Oelen et al. 2022, Nat Commun: blood
  cells of many donors stimulated with pathogens) and the influenza-infected blood cells of Randolph et al. 2021
  (Science).

## 3. The run (Actions, once, by the label `level-b-controls`; the branch `results/level-b-controls`)

The script `controls.py` and its workflow are committed before the run, with a dry run on fabricated Census
tables and fabricated GEO files on every push. The run writes: the Census release; every qualifying candidate
with the facts that decide it, ranked (`candidates.tsv`); each candidate tried and why it was or was not chosen;
for each control its cells' joinids, its extraction's SHA256 and its claims (`claims.json`); for Kang, its
files' SHA256, its structure check and its claims. The extractions themselves go to a workflow artifact (90 days,
as v1's backgrounds); the joinids and the release reproduce them.

## 4. What the list shows (S7)

For each control: its source and version, its donors and cells per condition, its checks, its claims with
their genes and directions, and the reason any candidate, gene or claim was left out.
