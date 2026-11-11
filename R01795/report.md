# Level B, S5: R01795's data and their checks (data_run.py; plan_s5.md)

Files: blish_awilk_covid_seurat.rds 5,985,045,048 bytes sha256 e81600242a5a3070740c158907704daf98c952eb06b5e0e6542a7c8da5ffa889

1. Counts: non-negative integers True.
2. Donors: nested, per group {'COVID': 6, 'Healthy': 6}; ids naming several donors none; rule 4 passes (smallest p 0.0021645021645021645).
3. Each id one individual: both sexes among the donors True; failing none.
4. The claims' genes: 1 missing CCN1, PLPP3; 2 missing WARS1, MARCHF1; 3 missing WARS1, MARCHF1; 4 missing WARS1, TENT5A; 5 missing RBSN; 7 missing WARS1, MARCHF1.
5. The source: balanced (sources ['PBMC']).
6. The paper's design (from the spec): the paper's: GSNCA's test of a Hallmark set's co-expression (CS-CORE correlations) between the CD4+ T cells of 7 hospitalized COVID-19 patients and 6 healthy donors (Wilk et al. 2020).

Cells: 5773 ({'COVID': 4535, 'Healthy': 1238}).
