# The third diagnostics and v0.4's acceptance tests on the development data (exploratory)

- Workflow run 39120635021 (attempt 1) of `.github/workflows/real-pilot-5.yml` (blob f2c63093cbcaf1359e2c929042118948819ac00c), at the pull request's head f5e1df0e94f7e86590fe3c2e2ebe978d6e522144; the plan: `validation/exploratory/real_pilot_diagnostics_3_plan.md` there.
- The engine: v0.4.0.dev0, this head's src/ (git tree fe2596bd0469cca1e3c4f9a8269d394c2b789367), installed in v1's analysis environment.
- The jobs: prepare success; extract success; run (the nine parts) success.
- The data: the second run's files on results/real-pilot-2 at bad785a3c56e9de1ee24737953c1b4b63112aa7c and the first run's case 3 on results/real-pilot at 6133dfdf4c05f083f012aaec4136e376715af05a (inputs/), not copied here; the full set of genes for their cells and case 1's metadata from the CELLxGENE Census, release 2025-11-08 (CC BY 4.0), the counts not copied here; Marsh et al. 2022's matrices from GEO (GSE152183), not copied here. extract/: the extraction's log and JSON, case 1's Census obs and its dataset's source obs.
- run/: real_pilot_5.log and real_pilot_5.json (every criterion and prediction scored as the plan fixes it), real_pilot_5_summary.log (the tallies), and each part's log, JSON and the engine's reports (parts/).
