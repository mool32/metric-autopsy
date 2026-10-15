# Step-1 report: which claims of preprint §4 hold

The flagship audit (validation-plan step 1, [README](README.md)) ran once, on the push of the
tag `v0.3.0-prereg` (commit 90a7718), in workflow run
[37990113298](https://github.com/mool32/metric-autopsy/actions/runs/37990113298). Its outputs
are on the branch `results/flagship-audit` at commit d2a5420: `out/audit_report.md`,
`out/audit_results.json` and the `--source-h5ad` pair (`out/*_source_h5ad.*`). The engine is
v0.1.1, the one §4 was computed with.

This report reads those files and decides which statements of §4 still hold. It does not edit
the manuscript: under the validation plan the paper stays untouched until step 3. The edits
listed at the end are for the author to decide.

Every number below comes from those outputs: `data_sha256`
d921bc225bc963f970a6b5b0fe4d9f7605b252a2c5f662257f9a70e994100b06, `script_sha256`
644e6ef5155f47ec62b2d65ccaa1af9f22336729dd9dc549ae51fec44e2845f7, Census `stable` (build
2025-11-08, schema 2.4.0). "Line" means a line of `paper/manuscript.md` at the tag.

**Provenance note.** Both outputs record `git_dirty: true`, but the code that ran is the tag's.
`script_sha256` and `download_module_sha256` (e1617f3b…) equal the sha256 of
`validation/flagship_audit/audit_tms.py` and `examples/mi_coupling_tms/download_data.py` at the
tag. The flag comes from untracked files in the checkout:
- The workflow passes `--cache audit_cache`, so the cache's `census_meta_stable.json` sits in the
  checkout before the provenance is taken. `.gitignore` covers only `*.h5ad`, not this file.
- The second run also finds `out/` in the checkout.

## In short

- **The engine's verdict reproduces.**
  - v0.1.1 kills `mi_3bin` on `Smad3`–`Col1a1` at GATE 0, and GATES 1, 2 and 5 fail as §4 says.
  - §4's other cell-level numbers reproduce at the stated precision, with three exceptions:
    the GATE 0 shift (43% against §4's 42%), its z (54.7 against 63.2) and the negative-control
    MI (0.036 against 0.038).
- **The design §4 describes does not hold.**
  - Its "old (20 months)" group is 4 males aged 24 months and 3 females aged 21 months. All 728
    cells of the 3 females come from the mammary gland.
  - §4's "correction" of the ages in ref 1 is wrong.
  - Two claims are not established at the level of mice: "female detection does not decline"
    and "a sex-by-age interaction".
- **"Detection artifact" is not established.** Neither the Census nor any of the 17 source files
  has ERCC spike-ins. So a technical drop in capture in old males cannot be told apart from a
  biological drop in RNA content.
- **On the same data, v0.1.1 also kills a true difference.**
  - The positive controls are Xist and the Y-gene score, female against male. Every mouse's
    markers match its sex label.
  - Both die at GATE 1. So a GATE 1 failure on this data is no evidence against biology.

## §4 claim by claim

| # | Claim (line) | Audit | Holds? |
|---|---|---|---|
| 1 | 110,824 Smart-seq2 cells, all primary (154) | `n_cells` 110,824 | yes |
| 2 | young 3 months against old 20 months (154; table header, 158) | The `20m` token is a single Census term, "20-month-old stage and over". It holds 4 males aged 24 months by their ids (31,551 cells) and 3 females aged 21 months (728 cells). | no |
| 3 | the ages are "3 / 18 / 20 months, not 3 / 24", correcting ref 1 (163) | males 3 / 18 / 24; females 3 / 18 / 21 | no: for males, ref 1's 24 months is right |
| 4 | median `nnz`: male 2799 / 1882 / 1701, female 2495 / 2291 / 3413 (160–161) | 2799 / 1882 / 1701 and 2495 / 2291 / 3413.5. The 18-month values are from section F's GATE 1 table. | yes |
| 5 | male young/old ratio 1.65x; GATE 1 flags males and passes females (163) | male 1.646, flagged; female 1.368 | yes |
| 6 | female detection does not decline; the 20-month female group is 728 cells "and noisy" (163) | The 728 cells are 3 mice of 222–267 cells each, all from the mammary gland; young females span 22 tissues. Per mouse, female young against old on `median_nnz`: p = 0.64 (5 against 3 units). | not established: the "20-month" contrast sets one tissue against many, and the 18-month point has 2 mice |
| 7 | the confound is a sex-by-age interaction (163) | Males: every old unit detects fewer genes than every young unit (`median_nnz` 1574–1975 against 2327–3970; exact p = 0.016, 5 against 4 units). Females: no tissue-matched contrast (row 6). | male age effect on detection: yes; the interaction: not established |
| 8 | GATE 0: a 42% shift under dropout, z = 63.2; FAIL at GATE 0 (167) | 43.0%, z = 54.7; FAIL at GATE 0 | verdict: yes; shift: close; z: no |
| 9 | GATE 2: the pooled effect retains 43%; the male one retains 30% and flips sign (167) | pooled 0.432 (0.006787 → −0.002933); male 0.300 (0.01185 → −0.003556) | yes; the pooled effect flips sign too |
| 10 | negative control `Tst`–`Lrrc42` has MI 0.038; GATE 5 fails (167) | MI 0.0362; GATE 5 fails in 2 of 2 strata | GATE 5: yes; MI: 0.036, not 0.038 |
| 11 | zero fraction 0.92 → 0.97 (169) | 0.917 → 0.971 | yes |
| 12 | "fresh automated-gate results" (154) | No script for §4 was ever committed (README, offline finding 4). The Census `var_names` are not gene symbols (`census_var_names_were_symbols` false), so the documented command cannot run as written. | the numbers now have a script: this audit |
| 13 | "pooled across all 23 tissues" (163) | The audit lists tissues only for the 3- and 20-month groups (22 tissues). | not checked |
| 14 | "confounded by a sex-by-age detection artifact rather than measuring biology" (13, the abstract); "not measuring biology" (173) | No spike-ins: `census_var_ercc` is 0 of 53,384 features, and 0 in each of the 17 source files. See also rows 6–7. | GATE 0's verdict: yes; "artifact": no; "sex-by-age": no |
| 15 | §3: "the real fact that female cells are stable while only male-old is degraded" (136) | rows 6–7 | not established |

## The mouse as the unit (section E)

§4 is computed over cells. The audit also tests per mouse (exact permutation within sex, young
against old):

| sex | units, young against old | `mi_3bin` mean | p | `median_nnz` mean | p |
|---|---|---|---|---|---|
| male | 5 against 4 | 0.0286 against 0.0088 | 0.016 | 3002 against 1704 | 0.016 |
| female | 5 against 3 | 0.0119 against 0.0136 | 0.68 | 3134 against 3392 | 0.64 |

In males, the lower MI and the lower detection separate exactly the same mice: every old unit
is below every young unit on both. So at the level of mice, the MI drop cannot be told apart
from the detection drop.

**Composite ids.** Two of the young "units" are composite ids: `3_10_M/3_11_M` (63 cells) and
`3_38_F/3_39_F` (75 cells). Each names two mice that are also counted on their own.
- The audit counts them as units. Each young group is therefore 4 mice, not 5, and the p-values
  above are those for 5 units.
- Without them, the design is 10 males (4 / 2 / 4 at 3 / 18 / 24 months) and 6 females
  (4 / 2 at 3 / 18 months). This matches the primary description's 16 mice (README, offline
  finding 1). On top of them come the 3 females at 21 months, mammary gland only.

## The same-data positive control (section F)

Xist and the Y-gene score, female against male within age, are a known difference: the sex
markers match the sex label in 21 of 21 donor units.
- v0.1.1 kills both at GATE 1. In the 20-month stratum the ratio is 2.01x, with 728 female
  cells against 31,551 male cells.
- GATE 2 keeps both effects under matching: Xist 0.4199 → 0.4216 (100%), Y score 104%.
- This is the death the README predicted for v0.1.1.

So a GATE 1 failure on this data does not tell an artifact from biology, and §4's GATE 1
failure (row 5) is no evidence against biology on its own. GATE 0, the gate that kills
`mi_3bin`, cannot be evaluated for the sex markers on the audit's gene-subset object.

## What the audit does not settle

- Whether old males detect fewer genes for technical or biological reasons: there are no
  spike-ins.
- Tissue-matched contrasts, such as female mammary gland young against old: not computed.
- The 18-month mice at the level of mice: section E compares 3 against 20 months only.
- The v0.3 engine's verdict on this data: the audit is pinned to v0.1.1 by design.

## Edits for the author to decide (nothing here is edited)

1. **Ages.** Write "old: 24-month males; 21-month females, mammary gland only" in place of
   "old (20 months)" (154, 158, figure caption 165). Withdraw the correction of ref 1's ages
   (163).
2. **Female arm.** State that it is 3 mice from one tissue. Either drop "female detection does
   not decline" (163, 165, 136) or limit it to the cell-level 18-month comparison, which has
   2 mice.
3. **"Artifact".** Write "detection difference" in place of "detection artifact" (13; also
   132, §3's account of ref 1): without spike-ins, the difference cannot be called technical.
4. **Numbers.** GATE 0's z is 54.7, not 63.2, and the negative-control MI is 0.036, not 0.038
   (167).
5. **Mouse level.** Add that the male MI drop and the male detection drop separate the same
   mice: exact p = 0.016 for each, 5 against 4 units, including one pooled id.
6. **Positive control.** Add that v0.1.1 also kills Xist and the Y-gene score at GATE 1 on the
   same data.
7. **Script.** Point §4's numbers to `validation/flagship_audit/audit_tms.py` and its
   `data_sha256`.
