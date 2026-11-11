# Level B, S5: the data checks, batch 3

- Workflow run 42939298430 of .github/workflows/level-b-data.yml at the pull request's head b8552b36cb22914e61d3b303cd2d8c817b059551; validation/level_b/data_run.py (blob 612ecaa2f6d1692ed1768972d1b29388e922b431) on each spec of run3.txt; the plan validation/level_b/plan_s5.md (blob 91f46896c7f5195cb9551452eebc5af581a3400e), committed before the run. The data job: success.
- Per dataset: spec.json, report.json and report.md (the six checks), files.tsv (each file's URL, size and SHA256), data_run.log; for an R object also fetch/ (its download), export.log (r_export.R export) and r_versions.txt; a dataset without report.md failed (its log, if any, says where).
