# metric-autopsy — status card

> The single source of truth for where this project is. Update it at every stage transition.
> Lifecycle and rules: ../_meta/RESEARCH_FLOW.md

**Stage:** 0 Seed · 1 Pre-reg · 2 Execute · 3 Verdict · 4 Write-up · **5 Publish** · 6 Archive  ← current
**One-liner:** A gate system — shipped as a Claude Code skill, a pip package, *and* an MCP server — that red-teams a computed single-cell metric to tell biological signal apart from QC/technical/mathematical artifacts.
**Started:** 2026-07-04   **Last update:** 2026-07-04

## Links
- GitHub: https://github.com/mool32/metric-autopsy (public) · CI green (3.9–3.12)
- Preprint: none yet — manuscript ready in `paper/manuscript.md`; bioRxiv posting pending
- Zenodo DOI: **10.5281/zenodo.21195679** (concept) · v0.1.0 `…680` — badge live
- Portfolio entry: pending — `mool32.github.io/_data/publications.yml` + `papers.bib`

## Origin
Direct descendant of `perceptual_modules/paper1/oscilatory/docs/metric_validation_checklist.md`
— the "work over the errors" written after ACP, β-heart, and MI-coupling each failed. This repo
turns that checklist into runnable behavior.

## Seed (stage 0)
- **Question:** Can the "compute → believe" default be replaced by a reusable "state commitments
  → red-team → then believe" workflow that catches QC/technical/mathematical artifacts *before* a
  claim is made — packaged so it reaches both the Claude-skill audience and the pip audience?
- **Why it matters:** The three failures cost weeks each and would have produced false papers.
  The checks cost ~1 hour. A tool that enforces them, with dual distribution, is high-leverage.
- **What would falsify it:** If the gates cannot separate a known-confounded metric (mi_3bin on a
  planted QC confound) from a known-clean one on the same data — i.e. if they neither catch the
  artifact nor pass the real signal. *(Currently: 32 tests show the separation holds, including
  adversarial-audit regressions and SimpleData↔AnnData agreement.)*

## Design decisions (locked for v1)
- **Metric as a plugin.** Gates take `metric(data) -> float` as a black box; they know scRNA-seq
  QC, not your metric. Keeps focus on single-cell while accepting any metric.
- **Triple purpose, one engine.** Skill (`SKILL.md`) + pip package (`metric-autopsy`) + worked
  example (`examples/`). `src/metric_autopsy/` is the single source of truth.
- **Scripted in v1:** GATES 0, 1, 2, 3, 5 (1 & 2 are the crown jewels). GATE 6 runs when a second
  dataset is supplied; it is fully demonstrated only in the worked example. GATES 4 & 7 are
  judgment — the skill elicits them, they are not scripted.

## Self-test (this is a tool, not a hypothesis test)
- **Result:** engine **35/35 tests pass** (synthetic gate tests + 24-finding adversarial-audit
  regressions + SimpleData↔AnnData compatibility). The falsification criterion holds: the gates
  separate a known-confounded metric (`mi_3bin`) from a known-clean one (`norm_pearson`) on the
  same planted data. Engine + docs both passed an adversarial multi-agent review.

## Data
- See DATASETS.md. All public (Tabula Muris Senis, human skin CELLxGENE). Nothing irreplaceable.

## Open threads / next
- [x] Worked-example notebook `examples/mi_coupling_tms/notebook.ipynb` — executed on synthetic
  `demo_data`, with committed outputs and figures. (Optional: re-run on real TMS after data pull.)
- [x] `references/gates.md` reworked into the preprint methods (`paper/manuscript.md`).
- [x] Public repo + `gh` metadata (description, homepage→DOI, topics incl. `tool`/`single-cell`).
- [x] Release v0.1.0 → Zenodo DOI (concept + version) → DOI badge; CI green across Python 3.9–3.12.
- [ ] Post preprint to bioRxiv (manuscript ready); then fill preprint DOI in README/CITATION/manuscript.
- [ ] Portfolio: add to `mool32.github.io/_data/publications.yml` + `papers.bib`.
- [ ] v1.1: turnkey GATE 6 second-platform replication, more example datasets.
