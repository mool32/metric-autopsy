# The pilot on real data with v0.3.0's engine (exploratory, development data)

- Workflow run 38892540972 (attempt 1) of `.github/workflows/real-pilot.yml` (blob 4cd1c13ef7950fc04d2be8a71cb496336130f95c), at the pull request's head 08b04fdd415eb480689e8f1d445fae169689115b; the plan: `validation/exploratory/real_pilot_plan.md` there.
- The engine: the tag v0.3.0-prereg, src/ git tree 97cbab3c437c9ddb118a0690e6e0d24f8a16069f (sha256 of its files' sha256 list 17e9a04f6fdd82e98606da89ad98cdce9a7f22fc35c1b31235efcb0637b132e8), installed from a worktree.
- The jobs: select success; run (the five parts) success.
- select/: the extraction's log, select.json (every choice and the facts that decided it), the environment. data/: the cells of each case (panel.save_npz; files over 45 MB split into .part-NN, to be concatenated). run/: the merged log and JSON (real_pilot.log, real_pilot.json), the engine's reports (reports/) and each part's own outputs (parts/).
- Data: the CELLxGENE Census, release 2025-11-08 (CC BY 4.0); the datasets' ids and DOIs are in select/select.json.
