# Flagship-data audit (§4, TMS FACS) — engine v0.1.1


## Provenance

```json
{
  "census_version": "stable",
  "census_summary": {
    "census_schema_version": "2.4.0",
    "census_build_date": "2025-11-08",
    "dataset_schema_version": "7.0.0",
    "total_cell_count": "217768036",
    "unique_cell_count": "125463259"
  },
  "tms_dataset_ids": [
    "524179b0-b406-4723-9c46-293ffa77ca81",
    "de4e7a0c-91b2-44e4-b382-87da74c9efb6",
    "e80d4e1c-672f-496a-8f32-37eab34f727d",
    "ec6c52b8-3368-4f72-a416-1ade0dab97bf",
    "4546e757-34d0-4d17-be06-538318925fcd",
    "170ce19f-7a2f-4926-a1cc-adcad99e7474",
    "3f4fe86f-aced-4d10-b174-ee35b9f46b9d",
    "bf12f9c6-4211-4c91-9c71-22019f29f516",
    "bc7466d7-ff13-4ff2-9c3d-7a1d208bd492",
    "05e6f6e3-0473-4b85-9f94-bcc5f1b5e04b",
    "e3b8c485-7811-407e-99ed-c7d574be9d7c",
    "c9096ac4-ea44-4cf9-82f4-af05cb83eb24",
    "1d29fd10-c8b3-4611-b0ac-3c578125adbf",
    "1efd4700-87dd-4b45-8762-11ba3fea7a65",
    "c08f8441-4a10-4748-872a-e70c0bcccdba",
    "db55b719-6102-493a-9251-404bc501d0de",
    "98e5ea9f-16d6-47ec-a529-686e76515e39"
  ],
  "tms_dataset_titles": [
    "Kidney - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "Large intestine - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "Spleen - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "Limb muscle - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "Liver - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "Thymus - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "Trachea - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "Bladder lumen - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "Mammary gland - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "Lung - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "Pancreas - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "Tongue - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "Skin of body - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "Heart - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "Brain myeloid cells - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "Bone marrow - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2",
    "All - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2"
  ],
  "census_var_n": 53384,
  "census_var_ercc": 0,
  "census_var_names_were_symbols": false,
  "engine_version": "0.1.1",
  "n_cells": 110824,
  "provenance": {
    "timestamp_utc": "2026-10-09T20:55:48+00:00",
    "git_commit": "90a7718f77a46440ca11cd90150733a67f6b138d",
    "git_dirty": true,
    "script_sha256": "644e6ef5155f47ec62b2d65ccaa1af9f22336729dd9dc549ae51fec44e2845f7",
    "download_module_sha256": "e1617f3bd14edab5649b29321a28d6c0d3946b21559234c826741a22f3112e66",
    "data_sha256": "d921bc225bc963f970a6b5b0fe4d9f7605b252a2c5f662257f9a70e994100b06",
    "data_hash_obs_columns": [
      "donor_id",
      "sex",
      "development_stage",
      "tissue",
      "assay",
      "age",
      "n_genes_by_counts",
      "total_counts"
    ],
    "cached_h5ad_sha256": "893c692a2a851787eaf0942983fdfdbb6208a34d07dcee44c06f763cc38f392e",
    "params": {
      "census_version": "stable",
      "young": "3m",
      "old": "20m",
      "source_h5ad": false,
      "dry_run_synthetic": false,
      "allow_engine_version": false
    },
    "params_sha256": "bfeb72acffb9f8930eb710b46fd788b25963f59a0ac20b2bf49e4fb6fd00695b",
    "versions": {
      "python": "3.12.15",
      "numpy": "2.5.3",
      "pandas": "2.3.3",
      "anndata": "0.13.4",
      "scipy": "1.18.1",
      "cellxgene_census": "1.18.0",
      "metric_autopsy": "0.1.1"
    }
  }
}
```

## A. development_stage × sex (cells, unique mice)

| development_stage | sex | n_cells | n_mice |
|---|---|---|---|
| 18-month-old stage | female | 19163 | 2 |
| 18-month-old stage | male | 14864 | 2 |
| 20-month-old stage and over | female | 728 | 3 |
| 20-month-old stage and over | male | 31551 | 4 |
| 3-month-old stage | female | 19315 | 5 |
| 3-month-old stage | male | 25203 | 5 |

## B. Age regex (download_data.py) and mouse ids

Mapping raw → token:

```json
{
  "18-month-old stage": "18m",
  "20-month-old stage and over": "20m",
  "3-month-old stage": "3m"
}
```

Collisions (one token, several raw strings): `{}`

Groups as the §4 analysis sees them (sex × token):

| sex | age_token | n_mice | n_cells | donor_ages | mice |
|---|---|---|---|---|---|
| female | 18m | 2 | 19163 | [18] | ['18_46_F', '18_47_F'] |
| female | 20m | 3 | 728 | [21] | ['21_48_F', '21_54_F', '21_55_F'] |
| female | 3m | 5 | 19315 | [3] | ['3_38_F', '3_38_F/3_39_F', '3_39_F', '3_56_F', '3_57_F'] |
| male | 18m | 2 | 14864 | [18] | ['18_45_M', '18_53_M'] |
| male | 20m | 4 | 31551 | [24] | ['24_58_M', '24_59_M', '24_60_M', '24_61_M'] |
| male | 3m | 5 | 25203 | [3] | ['3_10_M', '3_10_M/3_11_M', '3_11_M', '3_8_M', '3_9_M'] |

Mice whose id-encoded age disagrees with their token:

| donor_id | sex | age_token | donor_age_months | stages | n_cells |
|---|---|---|---|---|---|
| 21_48_F | female | 20m | 21 | ['20-month-old stage and over'] | 222 |
| 21_54_F | female | 20m | 21 | ['20-month-old stage and over'] | 267 |
| 21_55_F | female | 20m | 21 | ['20-month-old stage and over'] | 239 |
| 24_58_M | male | 20m | 24 | ['20-month-old stage and over'] | 8031 |
| 24_59_M | male | 20m | 24 | ['20-month-old stage and over'] | 7897 |
| 24_60_M | male | 20m | 24 | ['20-month-old stage and over'] | 8353 |
| 24_61_M | male | 20m | 24 | ['20-month-old stage and over'] | 7270 |

## C. ERCC spike-ins

```json
{
  "census_var_ercc": 0
}
```

## D. §4 recomputed with the pinned engine (cell level; groups ('3m', '20m'))

# Metric autopsy — mi_3bin

## Gates

| Gate | Name | Status | Finding |
|---|---|---|---|
| 0 | Mathematical independence | **FAIL** | metric's expectation shifts 43% (z=54.7) under 'extra_dropout' — confounded by that nuisance |
| 1 | QC parity | **FAIL** | 1/2 strata fail QC parity: 1 exceed 1.5x QC ratio (worst 1.65x at {'sex': 'male'}) — unusable without n_genes matching (GATE 2) |
| 2 | n_genes matching | **FAIL** | effect collapses under matching: 0.006787 -> -0.002933 (43% retained, need >= 50%) — likely a QC artifact |
| 3 | Raw visibility | **JUDGMENT** | zero-fractions comparable (0.92 vs 0.97) — inspect the scatter for real coupling |
| 5 | Controls | **FAIL** | 2/2 strata fail controls (e.g. {'sex': 'female'}: pos=0.0626, neg=0.0443; need \|pos\|>0.0169, \|neg\|<=0.0169) |

## Verdict

**FAIL — died at GATE 0 (Mathematical independence)**


| quantity | §4 says | recomputed |
|---|---|---|
| GATE 0 dropout rel / z | 0.42 / 63.2 | 0.43016952708443573 / 54.72077558287852 |
| GATE 2 pooled retained | 0.43 | 0.4321142620538533 (effect collapses under matching: 0.006787 -> -0.002933 (43% retained, need >= 50%) — likely a QC artifact) |
| GATE 2 male retained | 0.3 | 0.29993413885965653 (effect collapses under matching: 0.01185 -> -0.003556 (30% retained, need >= 50%) — likely a QC artifact) |
| negative control MI (Tst-Lrrc42) | 0.038 | 0.03625 |
| zero-fraction young/old | (0.92, 0.97) | (0.917, 0.971) |

Median nnz by sex × age (cells): {"female:20m": 3413.5, "female:3m": 2495.0, "male:20m": 1701.0, "male:3m": 2799.0}

GATE 1 table:

| stratum | n_a | n_b | median_n_genes_a | median_n_genes_b | n_genes_ratio | overlap | ratio_breach | overlap_breach | flagged |
|---|---|---|---|---|---|---|---|---|---|
| {'sex': 'female'} | 19315 | 728 | 2495 | 3414 | 1.368 | 1 | False | False | False |
| {'sex': 'male'} | 25203 | 31551 | 2799 | 1701 | 1.646 | 0.7607 | True | False | True |

## E. Mouse as the unit

| donor_id | sex | age | donor_age_months | n_cells | n_tissues | median_nnz | mi_3bin |
|---|---|---|---|---|---|---|---|
| 24_60_M | male | 20m | 24 | 8353 | 20 | 1670 | 0.006265 |
| 24_59_M | male | 20m | 24 | 7897 | 20 | 1574 | 0.01112 |
| 24_61_M | male | 20m | 24 | 7270 | 18 | 1975 | 0.01075 |
| 24_58_M | male | 20m | 24 | 8031 | 20 | 1598 | 0.00724 |
| 3_39_F | female | 3m | 3 | 7985 | 21 | 2626 | 0.007005 |
| 3_38_F | female | 3m | 3 | 7814 | 21 | 2100 | 0.004896 |
| 3_56_F | female | 3m | 3 | 3041 | 10 | 3025 | 0.01636 |
| 3_9_M | male | 3m | 3 | 7008 | 19 | 2830 | 0.02207 |
| 3_10_M | male | 3m | 3 | 7244 | 18 | 2868 | 0.01656 |
| 3_8_M | male | 3m | 3 | 7071 | 19 | 3018 | 0.02604 |
| 3_11_M | male | 3m | 3 | 3817 | 16 | 2327 | 0.03137 |
| 3_10_M/3_11_M | male | 3m | 3 | 63 | 1 | 3970 | 0.04692 |
| 3_57_F | female | 3m | 3 | 400 | 2 | 3679 | 0.009046 |
| 3_38_F/3_39_F | female | 3m | 3 | 75 | 1 | 4239 | 0.022 |
| 21_55_F | female | 20m | 21 | 239 | 1 | 3303 | 0.00963 |
| 21_48_F | female | 20m | 21 | 222 | 1 | 3368 | 0.02127 |
| 21_54_F | female | 20m | 21 | 267 | 1 | 3507 | 0.009892 |

Exact permutation over mice within sex (young vs old):

| test | n_young | n_old | mean_young | mean_old | p_exact | min_attainable_p |
|---|---|---|---|---|---|---|
| female:median_nnz | 5 | 3 | 3134 | 3392 | 0.6429 | 0.01786 |
| female:mi_3bin | 5 | 3 | 0.01186 | 0.0136 | 0.6786 | 0.01786 |
| male:median_nnz | 5 | 4 | 3002 | 1704 | 0.01587 | 0.007937 |
| male:mi_3bin | 5 | 4 | 0.02859 | 0.008844 | 0.01587 | 0.007937 |

Cells per tissue × (sex, age):

| tissue | female:20m | female:3m | male:20m | male:3m |
|---|---|---|---|---|
| aorta | 0 | 87 | 224 | 279 |
| bladder lumen | 0 | 328 | 434 | 978 |
| bone marrow | 0 | 1649 | 4664 | 3420 |
| brain | 0 | 2960 | 5180 | 5014 |
| brown adipose tissue | 0 | 278 | 848 | 435 |
| diaphragm | 0 | 435 | 705 | 468 |
| gonadal fat pad | 0 | 562 | 1067 | 902 |
| heart | 0 | 2499 | 3185 | 1934 |
| kidney | 0 | 189 | 663 | 313 |
| large intestine | 0 | 1979 | 1955 | 2008 |
| limb muscle | 0 | 524 | 1232 | 578 |
| liver | 0 | 161 | 948 | 570 |
| lung | 0 | 716 | 1877 | 665 |
| mammary gland | 728 | 1912 | 0 | 0 |
| mesenteric fat pad | 0 | 539 | 773 | 648 |
| pancreas | 0 | 733 | 963 | 855 |
| skin of body | 0 | 696 | 1122 | 1650 |
| spleen | 0 | 624 | 1110 | 1078 |
| subcutaneous adipose tissue | 0 | 788 | 1002 | 933 |
| thymus | 0 | 668 | 1242 | 691 |
| tongue | 0 | 403 | 1358 | 1015 |
| trachea | 0 | 585 | 999 | 769 |

## F. Same-data positive controls (female vs male, within age)

### Xist

| Gate | Status | Finding |
|---|---|---|
| 0 | not evaluable | gene-subset object: perturbations do not update the obs library size |
| 1 | FAIL | 1/3 strata fail QC parity: 1 exceed 1.5x QC ratio (worst 2.01x at {'age': '20m'}) — unusable without n_genes matching (GATE 2) |
| 2 | PASS | effect survives matching: 0.4199 -> 0.4216 (100% retained, subset balanced 1.05x) |

**v0.1.1 verdict (GATES 1-2): FAIL — died at GATE 1 (QC parity)**

### Y_score

| Gate | Status | Finding |
|---|---|---|
| 0 | not evaluable | gene-subset object: perturbations do not update the obs library size |
| 1 | FAIL | 1/3 strata fail QC parity: 1 exceed 1.5x QC ratio (worst 2.01x at {'age': '20m'}) — unusable without n_genes matching (GATE 2) |
| 2 | PASS | effect survives matching: -0.1194 -> -0.1247 (104% retained, subset balanced 1.05x) |

**v0.1.1 verdict (GATES 1-2): FAIL — died at GATE 1 (QC parity)**

Per-mouse sex markers:

| donor_id | sex | age | n_cells | xist | y_score | label_consistent |
|---|---|---|---|---|---|---|
| 18_53_M | male | 18m | 7014 | 0.009416 | 0.1036 | True |
| 18_45_M | male | 18m | 7850 | 0.01133 | 0.09643 | True |
| 18_47_F | female | 18m | 9562 | 0.4582 | 0.0005749 | True |
| 18_46_F | female | 18m | 9601 | 0.3943 | 0.0006623 | True |
| 24_60_M | male | 20m | 8353 | 0.00191 | 0.08897 | True |
| 24_59_M | male | 20m | 7897 | 0.001841 | 0.08939 | True |
| 24_61_M | male | 20m | 7270 | 0.002613 | 0.09567 | True |
| 24_58_M | male | 20m | 8031 | 0.003273 | 0.1016 | True |
| 3_39_F | female | 3m | 7985 | 0.5012 | 9.248e-05 | True |
| 3_38_F | female | 3m | 7814 | 0.3206 | 0.0002225 | True |
| 3_56_F | female | 3m | 3041 | 0.4926 | 0.0001424 | True |
| 3_9_M | male | 3m | 7008 | 0.003521 | 0.172 | True |
| 3_10_M | male | 3m | 7244 | 0.002119 | 0.1636 | True |
| 3_8_M | male | 3m | 7071 | 0.003646 | 0.1768 | True |
| 3_11_M | male | 3m | 3817 | 0.001389 | 0.1254 | True |
| 3_10_M/3_11_M | male | 3m | 63 | 0.0001569 | 0.2321 | True |
| 3_57_F | female | 3m | 400 | 0.3763 | 3.41e-05 | True |
| 3_38_F/3_39_F | female | 3m | 75 | 0.2617 | 0.001398 | True |
| 21_55_F | female | 20m | 239 | 0.4292 | 0 | True |
| 21_48_F | female | 20m | 222 | 0.3724 | 0 | True |
| 21_54_F | female | 20m | 267 | 0.3829 | 0 | True |
