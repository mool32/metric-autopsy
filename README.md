# The pilot's diagnostics on the second run's data, with v0.3.0's engine (exploratory, development data)

- Workflow run 39029096682 (attempt 1) of `.github/workflows/real-pilot-3.yml` (blob cab2c71918f9b3647fb0252216c55b464f52245c), at the pull request's head 8780ce7a6b3ca0cb391dec3f9ba62becee3caf31; the plan: `validation/exploratory/real_pilot_diagnostics_plan.md` there.
- The engine: the tag v0.3.0-prereg, src/ git tree 97cbab3c437c9ddb118a0690e6e0d24f8a16069f (sha256 of its files' sha256 list 17e9a04f6fdd82e98606da89ad98cdce9a7f22fc35c1b31235efcb0637b132e8), installed from a worktree.
- The jobs: prepare success; run (case 1 and case 2) success.
- The data: the second run's files on results/real-pilot-2 at bad785a3c56e9de1ee24737953c1b4b63112aa7c (inputs/inputs.log: each checked against that run's SHA256SUMS, copied to inputs/run2-SHA256SUMS); not copied here. STRING v12.0 (CC BY 4.0), mouse: https://stringdb-downloads.org/download/protein.info.v12.0/10090.protein.info.v12.0.txt.gz and .../protein.links.v12.0/10090.protein.links.v12.0.txt.gz, their sha256 in inputs/inputs.log; not copied here.
- run/: the merged log and JSON (real_pilot_3.log, real_pilot_3.json: the outcomes by check and claim, whether each check separates, the owner's prediction scored), the engine's reports of check (i) (reports/) and each part's own outputs (parts/).
