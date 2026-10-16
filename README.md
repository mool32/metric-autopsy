# The pilot on real data with v0.3.0's engine, second run: cases 1 and 2 (exploratory, development data)

- Workflow run 38932685766 (attempt 1) of `.github/workflows/real-pilot-2.yml` (blob e7ddaee69a1709aa8ff3a39bcd3f3254676ef9ff), at the pull request's head 6646514dd77f5c8529ab0ad938803356728229bf; the plan: `validation/exploratory/real_pilot_plan.md` there, its section "Deviation: a second run of cases 1 and 2". The first run is on the branch results/real-pilot.
- The engine: the tag v0.3.0-prereg, src/ git tree 97cbab3c437c9ddb118a0690e6e0d24f8a16069f (sha256 of its files' sha256 list 17e9a04f6fdd82e98606da89ad98cdce9a7f22fc35c1b31235efcb0637b132e8), installed from a worktree.
- The jobs: select success; run (the four parts) success.
- select/: the extraction's log, select.json (every choice and the facts that decided it), the environment. data/: the cells of each case (panel.save_npz; files over 45 MB split into .part-NN, to be concatenated). run/: the merged log and JSON (real_pilot.log, real_pilot.json), the engine's reports (reports/) and each part's own outputs (parts/).
- Data: the CELLxGENE Census, release 2025-11-08 (CC BY 4.0); the datasets' ids and DOIs are in select/select.json.
