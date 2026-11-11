# Level B, S5: the data checks, batch 4

- Workflow run 42942044094 of .github/workflows/level-b-data.yml at the pull request's head f7559076995718e21b1221a16b8d1cef802d217f; validation/level_b/data_run.py (blob 98966b7002eb23d41f1337262fa9a340d756e861) on each spec of run4.txt; the plan validation/level_b/plan_s5.md (blob 91f46896c7f5195cb9551452eebc5af581a3400e), committed before the run. The data job: success.
- Per dataset: spec.json, report.json and report.md (the six checks), files.tsv (each file's URL, size and SHA256), data_run.log; for an R object also fetch/ (its download), export.log (r_export.R export) and r_versions.txt; a dataset without report.md failed (its log, if any, says where).
