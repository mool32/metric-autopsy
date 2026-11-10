# Level B, S5: R01962's data and their checks (data_run.py; plan_s5.md)

Files: GSE266919_Myeloid.rds.gz 255,365,848 bytes sha256 54b63e7ae0a1d797d650e844be70780f9bd9c1b8d5e6c09bac2f35437ed57d54

1. Counts: non-negative integers True.
2. Donors: paired, per group {'Post-treatment': 10, 'Pre-treatment': 10}; ids naming several donors none; rule 4 passes (smallest p 0.001953125); donors with fewer than 30 cells in a group, left out as the paper did: {'P031': {'Post-treatment': 0, 'Pre-treatment': 132}, 'P039': {'Post-treatment': 1090, 'Pre-treatment': 15}, 'P053': {'Post-treatment': 17, 'Pre-treatment': 747}}.
3. Each id one individual: both sexes among the donors False; failing none.
4. The claims' genes: 1 present; 2 present.
5. The source: no source named.
6. The paper's design (from the spec): the paper's: scCOSMiX's copula correlation with the patients as random effects, the myeloid cells of 8 responders' breast biopsies before and after nab-paclitaxel with or without atezolizumab, 5,018 cells (2,215 before, 2,803 after).

Cells: 17461 ({'Pre-treatment': 10647, 'Post-treatment': 6814}).
