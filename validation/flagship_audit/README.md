# Flagship-data audit: preprint §4 (Tabula Muris Senis FACS)

Validation-plan step 1. The audit re-derives §4 **as computed**, before any change to the
verdict logic, and checks the data design under it: ages, mice per group, ERCC, and the
mouse as the unit of analysis. It does not edit the manuscript.

**Status (2026-10-07): script ready, not yet run on real data.** The cloud environment's
network policy blocks every data host this audit needs. `--dry-run-synthetic` runs the
full code path offline on fabricated TMS-shaped data, to check the code only.

## Run it

```bash
# 1. pin the engine §4 was computed with
git worktree add ../metric-autopsy-v0.1.1 v0.1.1
pip install -e ../metric-autopsy-v0.1.1 cellxgene-census

# 2. run from this checkout (the script lives here; the engine comes from the worktree)
python validation/flagship_audit/audit_tms.py --out validation/flagship_audit/out
python validation/flagship_audit/audit_tms.py --source-h5ad --out validation/flagship_audit/out  # + ERCC in source files (GB)

# offline check of the code path (SYNTHETIC; says nothing about TMS)
python validation/flagship_audit/audit_tms.py --dry-run-synthetic --out /tmp/audit_dry
```

Hosts it needs:
- `census.cellxgene.cziscience.com` (Census release directory)
- `cellxgene-census-public-us-west-2.s3.us-west-2.amazonaws.com` (Census data, and the source `.h5ad` for `--source-h5ad`)

Output: `out/audit_report.md` (tables) and `out/audit_results.json` (everything). The
script writes numbers, not conclusions. The human-written step-1 report decides which
§4 claims still hold.

## What it checks

| | Question | Why it matters for §4 |
|---|---|---|
| A | raw `development_stage` strings × sex, with cells and unique mice | §4's "old (20 months)" group and its replicate count |
| B | how the age regex maps those strings, and the age encoded in each mouse id | does "20m" merge 21- and 24-month mice? |
| C | ERCC spike-ins in the Census var and in the source `.h5ad` | the only way to tell a technical capture drop in male-old cells from a biological RNA-content drop |
| D | §4 recomputed with v0.1.1 at the cell level, side by side with the stated numbers | is §4 reproducible at all? (no script for it was ever committed) |
| E | the mouse as the unit: per-mouse QC and mi_3bin, exact permutation within sex, tissue composition | is the "sex × age interaction" distinguishable from individual mice? |
| F | Xist and Y genes, female vs male within age; per-mouse sex-marker check | the same-data positive control: v0.1.1's prediction is death at GATE 1 |

## Already established offline

1. **The primary description conflicts with §4.** Zhang et al. 2021 (*eLife* 10:e62293)
   describe the TMS FACS data as 16 mice (10 male, 6 female). Males are at 3/18/24 months;
   females only at 3/18 months. §4 reports a female "20-month" group of 728 cells and
   "corrects" the ages to "3/18/20 months, not 3/24".
2. **The regex is faithful to its input.** On test strings, `download_data.py:62-65` maps
   "24-month-old stage" → `24m` and "20 month-old stage and over" → `20m`. A `20m` token
   therefore means the Census strings say "20 month". The working hypothesis (check A/B):
   CELLxGENE files 21- and 24-month mice under one "20 month-old stage and over"-type
   term. If so, §4's "old" group mixes 24-month males with a 21-month female group.
3. **Gene names.** `download_data.py` never sets `var_names` to gene symbols. Census
   AnnData objects default to a RangeIndex (`tiledbsoma._util._df_set_index`), so the
   README command `--gene-a Smad3` would raise `KeyError` on the downloaded file. The §4
   run must have included an undocumented step. The audit records
   `census_var_names_were_symbols`.
4. **No script behind §4.** Commit `df5fd76`, which added §4, changed only `paper/`.
   This audit is the first scripted reproduction.
5. **Caveat for section F.** On this gene-subset object, v0.1.1's GATE 0 is not evaluable
   for library-normalized metrics: its perturbations change X but not the obs library
   size. F therefore reports GATES 1–2, the first blocking gate among them, and the
   per-mouse values.
