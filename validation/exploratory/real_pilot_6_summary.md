# Preprint §4 on the fourth run of v0.4 (exploratory)

Written after the run by `real_pilot_6_summary.py` from the run's JSON (`results/real-pilot-6` at 1e2f824; workflow run 39163054704 on 57980d5, the engine v0.4.0.dev0; the plan `real_pilot_6_plan.md`). The claim (preprint §4, the setup of the v0.1.1 audit): `mi_3bin(Smad3, Col1a1)` falls from the 3-month to the 24-month mice of Tabula Muris Senis FACS. The differences below are young − old, the engine's convention (groups[0] − groups[1]): the claim's direction is a positive difference. A3 is the claim in males, A4 within sex. Only the mixture without the composite unit is computed here (the engine 0.4.0.dev0's `nested_t` and `permutation_test_nested`, as the run's mixture line); every other number is the run's. Exploratory: not in the preprint, not in the manuscript.

The files read, as the run's SHA256SUMS: `run/parts/part-case1m/reports/A3.json` 1d846ab09f51…; `run/parts/part-case1/reports/A4.json` 7004d1837d45…; `text/real_pilot_6_text.json` 309e7da9f157…; `extract/source_obs_case1.csv.gz` d29b931a29d9…; `README.md` 2e6e8f978e23….

- **A3:** INCONCLUSIVE [cause `effect_inconclusive`]: within the domain -0.001844 [-0.0179, +0.0142], exact p 0.771.
- **A4:** INCONCLUSIVE [cause `effect_inconclusive`]: within the domain -0.0003159 [-0.0168, +0.0162], exact p 1.

## 1. The whole mixture

At equal depth +0.0198 [+0.00561, +0.034] (raw +0.01975), exact p 0.0159 over 126 relabelings, 5 young against 4 old units: the smallest two-sided p these units can reach (2/126). The units, at equal depth:

| unit | group | at equal depth | raw | note |
|---|---|---|---|---|
| 3_10_M/3_11_M | young | +0.04692 | +0.04692 | composite: its id names 3_10_M and 3_11_M, each also a unit of its own; the largest value |
| 3_11_M | young | +0.03137 | +0.03137 |  |
| 3_8_M | young | +0.02604 | +0.02604 |  |
| 3_9_M | young | +0.02207 | +0.02207 |  |
| 3_10_M | young | +0.01656 | +0.01656 |  |
| 24_59_M | old | +0.01159 | +0.01112 |  |
| 24_61_M | old | +0.01098 | +0.01075 |  |
| 24_58_M | old | +0.006958 | +0.00724 |  |
| 24_60_M | old | +0.005636 | +0.006265 |  |

The same code on the nine units gives +0.0198 [+0.00561, +0.034], exact p 0.0159 — the run's estimate.

The composite unit `3_10_M/3_11_M` (the engine's flag `COMPOSITE_UNITS`) holds 63 cells, all in heart (63), on plate(s) MAA100037 (63); each cell's name gives one mouse: 3_10_M 31, 3_11_M 32. It enters only the mixture: no type has 20 of its cells, so the domain and the types outside it are estimated over 4 against 4 units.

**Without the composite unit** (4 young against 4 old units): +0.01522 [+0.00585, +0.0246], exact p 0.0286 over 70 relabelings — again the smallest p these units can reach (0.0286): every young unit lies above every old one. The mixture's difference is smaller by 0.00458.

## 2. Within the domain

The domain by the rule "in each group" (types where, in each group, both genes are detected in at least 10% of the type's cells; estimable with at least 20 cells in at least 3 replicates of each group): fibroblast, fibroblast of cardiac tissue, granulocytopoietic cell, mesenchymal stem cell, mesenchymal stem cell of adipose tissue. Combined, the types weighted equally: -0.001844 [-0.0179, +0.0142] (raw -0.004141), exact p 0.771 over 70 relabelings, 4 young against 4 old units: INCONCLUSIVE. The MDE at 80% power is 0.0219, larger than the mixture's effect (+0.0198, without the composite unit +0.0152): within the domain this design could not have detected a drop the size of the mixture's. The SESOI in the metric's units is s_m +0.00273, below the design's floor 0.0287.

| type | units | effect at equal depth | p | raw | Smad3 young | Smad3 old | Col1a1 young | Col1a1 old |
|---|---|---|---|---|---|---|---|---|
| fibroblast | 7 | -0.03323 | 0.0889 (parametric) | -0.03036 | 52.3% | 34.0% | 99.5% | 72.6% |
| fibroblast of cardiac tissue | 8 | -0.006174 | 0.686 (permutation) | -0.009285 | 43.5% | 25.8% | 97.0% | 71.6% |
| granulocytopoietic cell | 6 | +0.01915 | 0.282 (parametric) | +0.01179 | 19.3% | 13.0% | 10.9% | 10.4% |
| mesenchymal stem cell | 7 | +0.009297 | 0.176 (parametric) | +0.007052 | 35.4% | 19.5% | 98.7% | 53.1% |
| mesenchymal stem cell of adipose tissue | 8 | +0.001732 | 0.629 (permutation) | +9.85e-05 | 51.8% | 34.2% | 99.6% | 91.4% |

In both fibroblast types the difference is negative, against the claim: fibroblast -0.03323 (p 0.0889); fibroblast of cardiac tissue -0.006174 (p 0.686). Positive in 3 of 5 types of the domain.

## 3. Outside the domain

Combined over the 30 types outside the domain, the types weighted equally: +0.01717 [+0.0141, +0.0202] (raw +0.01651), exact p 0.0286 over 70 relabelings, 4 young against 4 old units (the smallest p they can reach). Per type, the claim genes' detection (the share of the type's cells with a count) in each group:

| type | units | effect at equal depth | p | raw | Smad3 young | Smad3 old | Col1a1 young | Col1a1 old |
|---|---|---|---|---|---|---|---|---|
| B cell | 8 | +0.001764 | 0.2 (permutation) | +0.001333 | 10.0% | 4.7% | 5.4% | 2.1% |
| CD4-positive, alpha-beta T cell | 8 | +0.03764 | 0.0857 (permutation) | +0.03427 | 19.4% | 13.6% | 8.5% | 1.2% |
| CD8-positive, alpha-beta T cell | 8 | -2.026e-05 | 1 (permutation) | -0.0002567 | 23.9% | 16.1% | 2.6% | 1.0% |
| DN4 thymocyte | 7 | +0.008751 | 0.22 (parametric) | +0.008699 | 17.0% | 18.3% | 4.1% | 0.2% |
| T cell | 8 | +0.0129 | 0.257 (permutation) | +0.0101 | 25.3% | 20.5% | 9.4% | 2.1% |
| basal cell of epidermis | 8 | +0.003893 | 0.0286 (permutation) | +0.00348 | 34.7% | 21.1% | 8.2% | 3.2% |
| bladder urothelial cell | 6 | +0.01936 | 0.172 (parametric) | +0.01936 | 77.1% | 42.4% | 3.4% | 0.0% |
| bronchial smooth muscle cell | 7 | +0.001209 | 0.427 (parametric) | +0.001209 | 1.7% | 0.3% | 8.9% | 0.7% |
| cardiac valve cell | 7 | +0.09738 | 0.0223 (parametric) | +0.09668 | 8.6% | 0.9% | 50.7% | 28.2% |
| endothelial cell | 8 | +0.004131 | 0.429 (permutation) | +0.004041 | 13.0% | 4.1% | 9.0% | 2.9% |
| endothelial cell of coronary artery | 8 | +0.008743 | 0.457 (permutation) | +0.009685 | 4.1% | 1.7% | 7.1% | 2.7% |
| enterocyte of epithelium of large intestine | 8 | +0.01063 | 0.114 (permutation) | +0.01186 | 36.1% | 16.4% | 11.9% | 1.0% |
| epidermal cell | 7 | +0.003867 | 0.159 (parametric) | +0.003867 | 16.4% | 9.8% | 6.0% | 1.3% |
| epithelial cell of large intestine | 8 | +0.009542 | 0.0286 (permutation) | +0.008971 | 11.5% | 9.5% | 12.3% | 2.5% |
| granulocyte | 7 | +0.01933 | 0.195 (parametric) | +0.02075 | 22.7% | 9.7% | 10.5% | 2.7% |
| hematopoietic stem cell | 7 | +0.002979 | 0.0326 (parametric) | +0.002582 | 55.3% | 25.5% | 7.0% | 0.5% |
| intestinal crypt stem cell | 7 | +0.04041 | 0.236 (parametric) | +0.04041 | 50.0% | 23.7% | 9.2% | 0.0% |
| keratinocyte | 7 | +0.01195 | 0.0503 (parametric) | +0.01184 | 28.9% | 17.8% | 11.0% | 4.3% |
| keratinocyte stem cell | 8 | +0.00958 | 0.0286 (permutation) | +0.008537 | 11.1% | 6.5% | 7.4% | 3.2% |
| large intestine goblet cell | 8 | +0.02956 | 0.0286 (permutation) | +0.0272 | 31.5% | 23.2% | 8.3% | 2.6% |
| macrophage | 8 | +0.01453 | 0.343 (permutation) | +0.01243 | 17.4% | 10.4% | 7.5% | 1.7% |
| microglial cell | 7 | +0.005 | 0.0201 (parametric) | +0.005171 | 42.9% | 19.3% | 6.3% | 0.9% |
| monocyte | 8 | +0.01876 | 0.0571 (permutation) | +0.02027 | 22.5% | 11.7% | 9.1% | 7.2% |
| myeloid cell | 8 | +0.01445 | 0.0571 (permutation) | +0.01479 | 22.3% | 13.9% | 9.4% | 5.8% |
| naive B cell | 7 | +0.01597 | 0.123 (parametric) | +0.01551 | 7.3% | 4.4% | 8.9% | 2.4% |
| natural killer cell | 6 | +0.006731 | 0.497 (parametric) | -0.0004627 | 49.3% | 17.3% | 8.5% | 3.0% |
| precursor B cell | 6 | +0.0295 | 0.0396 (parametric) | +0.0295 | 15.1% | 4.9% | 10.8% | 0.0% |
| promonocyte | 7 | +0.03362 | 0.0451 (parametric) | +0.03096 | 24.4% | 10.1% | 7.6% | 1.7% |
| skeletal muscle satellite cell | 8 | +0.02464 | 0.0286 (permutation) | +0.02445 | 34.5% | 10.8% | 12.1% | 1.1% |
| thymocyte | 6 | +0.01839 | 0.0979 (parametric) | +0.01814 | 17.2% | 20.2% | 6.6% | 0.2% |

- The sign: positive (the claim's direction) in 29 of 30 types at equal depth (28 of 30 raw); nominally p < 0.05 in 10.
- Col1a1 is detected more often in the young in 30 of 30 types, Smad3 in 28 of 30; the index of the claim genes' counts per cell is higher in the young in 26 of 30.
- Outside the domain Col1a1, a gene of fibroblasts, is detected in more than 5% of the young cells but under 2% of the old in 12 types, among them CD4-positive, alpha-beta T cell (8.5% against 1.2%), bronchial smooth muscle cell (8.9% against 0.7%), enterocyte of epithelium of large intestine (11.9% against 1.0%), epidermal cell (6.0% against 1.3%), hematopoietic stem cell (7.0% against 0.5%), intestinal crypt stem cell (9.2% against 0.0%).

## 4. The technical covariates per unit (the dossier's line)

The median cell of each unit, the young units' mean against the old's, with the exact permutation p (`Autopsy.technical`; genes and counts over the pool's 21,025 genes).

| claim | cells | covariate | stratum | young | old | difference | p (relabelings) | units, young / old |
|---|---|---|---|---|---|---|---|---|
| A3 | all cells | counts per cell | all | 427,512.6 | 616,511.1 | -188,998.5 | 0.00794 (126) | 5 / 4 |
| A3 | all cells | genes per cell | all | 3,002.5 | 1,704.2 | +1,298.2 | 0.0159 (126) | 5 / 4 |
| A3 | the domain's cells | counts per cell | all | 703,755.5 | 1,052,672.2 | -348,916.8 | 0.119 (126) | 5 / 4 |
| A3 | the domain's cells | genes per cell | all | 3,503.7 | 2,337.5 | +1,166.2 | 0.0159 (126) | 5 / 4 |
| A4 | all cells | counts per cell | male | 427,512.6 | 616,511.1 | -188,998.5 | 0.00794 (126) | 5 / 4 |
| A4 | all cells | counts per cell | female | 610,674.5 | 535,471.2 | +75,203.3 | 0.357 (56) | 5 / 3 |
| A4 | all cells | counts per cell | combined | 519,093.5 | 581,779.7 | -68,092.6 | 0.241 (7056) | 10 / 7 |
| A4 | all cells | genes per cell | male | 3,002.5 | 1,704.2 | +1,298.2 | 0.0159 (126) | 5 / 4 |
| A4 | all cells | genes per cell | female | 3,133.8 | 3,392.5 | -258.7 | 0.643 (56) | 5 / 3 |
| A4 | all cells | genes per cell | combined | 3,068.2 | 2,427.8 | +585.7 | 0.114 (7056) | 10 / 7 |
| A4 | the domain's cells | counts per cell | male | 703,755.5 | 1,052,672.2 | -348,916.8 | 0.119 (126) | 5 / 4 |
| A4 | the domain's cells | counts per cell | combined | 703,755.5 | 1,052,672.2 | -348,916.8 | 0.119 (126) | 5 / 4 |
| A4 | the domain's cells | genes per cell | male | 3,503.7 | 2,337.5 | +1,166.2 | 0.0159 (126) | 5 / 4 |
| A4 | the domain's cells | genes per cell | combined | 3,503.7 | 2,337.5 | +1,166.2 | 0.0159 (126) | 5 / 4 |

With p < 0.05: A3 all cells, counts per cell, all (-188,999); A3 all cells, genes per cell, all (+1,298); A3 the domain's cells, genes per cell, all (+1,166); A4 all cells, counts per cell, male (-188,999); A4 all cells, genes per cell, male (+1,298); A4 the domain's cells, genes per cell, male (+1,166); A4 the domain's cells, genes per cell, combined (+1,166); no female row. GATE 1's ratio of median genes passed the males at 1.32.

## 5. Attenuation, the probe and the correction on §4

- GATE 0: extra dropout removes 63.7% of the signal above the null (an attenuation, not a bias: the null moves by -0.00117, z 2.09, against the tolerance 0.00136 = 0.5 × s_m). GATE 0 passes an attenuation on as a loss of power, but where the groups differ in detection (section 4) it can create a difference at equal coupling (the differential sensitivity, decided 2026-10-19).
- GATE 4 injected the coupling into the young only, the group with more probed counts per cell (1605 against 652.9): s_m +0.00273 [+0.00119, +0.00427] is the young's response; the old's was not measured.
- The correction thinned the old cells to the young's depth; the gap in genes per cell (young − old, relative) grew from 24.4% to 29.1% on all cells, and from 16.9% to 21.4% on the domain's: equal depth makes the groups' detection less equal.

## 6. Within sex (A4)

The domain by stratum: males fibroblast, fibroblast of cardiac tissue, granulocytopoietic cell, mesenchymal stem cell, mesenchymal stem cell of adipose tissue; females none — their old units have no cell of these types (fibroblast 233 against 0; fibroblast of cardiac tissue 1354 against 0; granulocytopoietic cell 31 against 0; mesenchymal stem cell 262 against 0; mesenchymal stem cell of adipose tissue 884 against 0 cells), so the domain's estimate is the males'. Within the domain -0.0003159 [-0.0168, +0.0162] (exact p 1); the mixture +0.009945 (p 0.0522, exact permutation over replicates (n=7056)); outside the domain +0.01073 (p 0.00204, exact permutation over replicates (n=2450)). GATE 4's response, injected into the young, was not established (WARN).

## 7. What the TMS authors write

- tms-epmc (the Nature paper, Europe PMC's full text (an NIH author manuscript)): FAIL HTTPError: HTTP Error 500: .
- tms-bioc (the Nature paper, NCBI's BioC (PMC open access and author manuscripts)): ok 200; 108,731 bytes, sha256 1f19e3bc7c09….
- tms-biorxiv (the bioRxiv preprint, its JATS full text (the newest version, by bioRxiv's API)): ok 200; version 3, licence cc_by_nc_nd; 153,512 bytes, sha256 86f47ce066c1….

By the plan's rule (a sentence that names the 3-month data together with the Tabula Muris or with "previously published") the authors state that the 3-month data are the original Tabula Muris's: 6 such sentences, in 2 sources. Quoted, with the batch correction:

> The previously published 3m time point, referred to as the Tabula Muris, represents ~20% of the cells in the entire dataset and was used as a basis to perform semi-automated cell type annotation of the additional time points (Figure 1b, Extended Data Figure 4b).
>
> — the Nature paper (NCBI's BioC), INTRO

> To investigate whether cell annotations were consistent across the entire organism, we used bbknn to correct for method-associated batch effects (Supplementary Table 3).
>
> — the Nature paper (NCBI's BioC), INTRO

> When performing batch correction to remove the technical artifacts introduced by the technologies, we replaced the neighborhood graph computation with bbknn.
>
> — the Nature paper (NCBI's BioC), METHODS

> We used the Tabula Muris (3m time point) as a reference for the automated pipeline and the annotations were manually curated by tissue experts. c,d, UMAP plot of all cells, colored by organ and overlaid with the Louvain cluster numbers (c) and age (d); n = 356,213 individual cells.
>
> — the Nature paper (NCBI's BioC), FIG

The preprint (bioRxiv, version 3) has 4 such sentences, each matched to the paper's by its similarity (difflib's ratio, 1 = the same text; the figure and reference numbers differ): Overview of Tabula Muris Senis. (0.98); Overview of Tabula Muris Senis. (0.98); Overview of Tabula Muris Senis. (0.92); Data pre-processing (0.99).

bbknn corrects the neighbourhood graph between the two technologies (FACS and droplet); it does not change the counts, which the engine reads, and no quoted sentence names a correction between the 3-month data and the later time points. The later time points' cells were annotated with the 3-month cells as the reference: a comparison within a cell type compares labels given by two procedures.

## Conclusion

Most likely a difference in co-detection between two generations of data (Tabula Muris 2018 against Tabula Muris Senis), not a coupling; the engine does not establish this — the decisive steps are the biology (in which cells Col1a1 is expressed) and the authors' text. (The owner's wording, 2026-10-19: «скорее всего, различие в со-детекции между двумя поколениями данных, а не связь; движок этого не устанавливает — решающие шаги: биология (в каких клетках экспрессируется Col1a1) и текст авторов».)
