# Level B, S5: R02693's data and their checks (data_run.py; plan_s5.md)

Files: GSE178265_Homo_matrix.mtx.gz 5,089,750,754 bytes sha256 d39c5da8b7746e7c9e83db0f29a66e5d60a4e8131274f200e1ef4e9d9302aa9f; GSE178265_Homo_features.tsv.gz 295,300 bytes sha256 fa2b4e9cc742ef81747ccd9c348b5dc55095525f790e329168ff4f1daff8cbf9; GSE178265_Homo_bcd.tsv.gz 2,732,185 bytes sha256 5e8ccf574fb1afcedd9d258b7432b1c2f4e88c8747406f901d70f0121623751c; GSE178265-GPL20084_series_matrix.txt.gz 2,888 bytes sha256 62e26b61c082ab2fd9d7856971d659851b02439b666e579f89e42bc29cd20f31; GSE178265-GPL21120_series_matrix.txt.gz 3,110 bytes sha256 346ae8c998c7a1929fc535d4666c995e2e299dd787317bb6cc19a61e86e2e217; GSE178265-GPL24676_series_matrix.txt.gz 6,717 bytes sha256 9f0f0d75e50b796709e5c10e961ee4167b75b2122c0fb0f56d776bc0f763c74a; GSE178265-GPL30282_series_matrix.txt.gz 2,893 bytes sha256 b9e9f9cece2c1e35b19f030affcf25772c607a1a9e68b9e43eea99eab3f1d087

1. Counts: non-negative integers True.
2. Donors: nested, per group {'Disease': 10, 'Ctrl': 9}; ids naming several donors none; rule 4 passes (smallest p 2.165017644893806e-05).
3. Each id one individual: both sexes among the donors True; failing none.
4. The claims' genes: 1 present; 2 present; 3 present.
5. The source: partly (sources ['Negative', 'NeuN', 'Positive']).
6. The paper's design (from the spec): the paper's: Pearson correlation of log-normalized expression over the pooled nuclei of the substantia nigra, PD against neurologically healthy controls (11 against 10 in the paper).

Cells: 393873 ({'Disease': 202810, 'Ctrl': 191063}).
