# Level B, S5: R00765's data and their checks (data_run.py; plan_s5.md)

Files: GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5 273,975,534 bytes sha256 6ba98a1af8772af08c8cdd5e5e63eebb0890daed48cac2e1b8961dcc67069b77; GSE174367_snRNA-seq_cell_meta.csv.gz 435,170 bytes sha256 ab1a029deb43196e2bb1fea5907d838885750cfff3850c600c07588c7c7cdb2b

1. Counts: non-negative integers True.
2. Donors: nested, per group {'AD': 11, 'Control': 7}; ids naming several donors none; rule 4 passes (smallest p 6.284565108094519e-05).
3. Each id one individual: both sexes among the donors True; failing none.
4. The claims' genes: .
5. The source: balanced (sources ['1.0', '2.0', '3.0']).
6. The paper's design (from the spec): the paper's: Dozer's donor-level co-expression networks of the oligodendrocytes of Morabito et al. 2021 (snRNA-seq, prefrontal cortex), module density compared between late-stage AD and control donors by a t-test.

Cells: 37052 ({'AD': 22059, 'Control': 14993}).
