# Level B, S5: the data checks, batch 2

- Workflow run 42903541582 of .github/workflows/level-b-data.yml at the pull request's head b737f50e2c155388a81da7fcfc0226fce7f56f2f; validation/level_b/data_run.py (blob 05dfae3f7c84c72d3edc89fbd38402fc8253393e) on each spec of run2.txt; the plan validation/level_b/plan_s5.md (blob e24ea20b72d23de3963ef8572cc6532edfc36a51), committed before the run. The data job: success.
- Per dataset: spec.json, report.json and report.md (the six checks), files.tsv (each file's URL, size and SHA256), data_run.log; for an R object also fetch/ (its download), export.log (r_export.R export) and r_versions.txt; a dataset without report.md failed (its log, if any, says where).
