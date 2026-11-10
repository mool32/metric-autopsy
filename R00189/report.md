# Level B, S5: R00189's data and their checks (data_run.py; plan_s5.md)

Files: GSE134809_RAW.tar 607,467,520 bytes sha256 875d8232fcd6d489037b4d3bddf81c588fdfcedde3190bf6eb4522d0a040d84a; GSE134809_series_matrix.txt.gz 5,137 bytes sha256 b980b1d829ca680afddc924c72a3917d21ddbd935619bbdc8669ba9ac7eff198; 1-s2.0-S0092867419308967-mmc2.xlsx 7,540,057 bytes sha256 030e5057afa1ebdfb821a197607ad3e462d598e626ec9483fa9f4f7d0c32dec5

1. Counts: non-negative integers True.
2. Donors: paired, per group {'Involved': 10, 'Uninvolved': 10}; ids naming several donors none; rule 4 passes (smallest p 0.001953125); cells no donor names, left out: {'Ileal Involved 190 / Involved': 5227, 'Ileal Uninvolved 189 / Uninvolved': 5129}.
3. Each id one individual: both sexes among the donors True; failing none.
4. The claims' genes: 1 present; 2 present; 3 present; 4 present; 5 present.
5. The source: balanced (sources ['V1', 'V2']).
6. The paper's design (from the spec): the paper's: a Pearson correlation of the two genes over the absorptive enterocytes' cells pooled per condition (Fig 6G, 6I), inflamed against uninflamed ileum of 8 Crohn disease patients with both samples.

Cells: 90894 ({'Uninvolved': 46461, 'Involved': 44433}).
