# Level B, S5: R01225's data and their checks (data_run.py; plan_s5.md)

Files: GSE189955_RNA_Count_matrix.csv.gz 34,493,403 bytes sha256 ff71d5589547cde643a250a8e03ddaea3783e4427d05f1489aa62aa64f5ec1e8; GSE189955_RNA_meta_data.csv.gz 57,585 bytes sha256 dc74681915241668c20535ca6d20a3d22b74eede2c9e0fa31e9304e42fc7de3a

1. Counts: non-negative integers True.
2. Donors: nested, per group {'tumour': 16, 'normal': 5}; ids naming several donors none; rule 4 passes (smallest p 9.828492800629023e-05).
3. Each id one individual: both sexes among the donors False; failing none.
4. The claims' genes: 1 present; 2 present; 3 present; 4 present; 5 present.
5. The source: partly (sources ['Smart_flow', 'Smart_mouth', 'SplitCool', 'SplitIVMsnmC', 'Split_FACS', 'Split_mouth', 'unmatched']).
6. The paper's design (from the spec): the paper's: a correlation r per condition over the pooled cells (DiffCorr), 1,940 primary tumour cells against 520 epithelial cells of fallopian tubes.

Cells: 2460 ({'tumour': 1940, 'normal': 520}).
