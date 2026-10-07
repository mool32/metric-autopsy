# Dev-set journal: changes to frozen expectations

The expected verdicts in `test_probes.py` were frozen before any fix, and a fix must meet them
without editing them. This journal records every exception: what changed, why, who decided it,
and the evidence. An entry is written *before* the new expectation is tested.

The dev set has no confirmatory weight either way (see `README.md`). Confirmatory validation
runs on a frozen tag with a new panel; its own deviations go to `validation/prereg/`.

---

## D1 — p11 split by design (2026-10-07; decided by the project owner)

**Frozen expectation.** p11 (demo male stratum: identical biology, old-male capture 30%,
4 vs 4 mice of 100 cells): NOT SUPPORTED, explained by depth.

**Why it changes.** The expectation asked for more than the design can establish. With 4 vs 4
mice the exact permutation over replicates has 70 assignments, so a raw difference reaches
p < 0.05 only at complete separation of the mice (p = 0.029; the next attainable values are
0.057 and 0.086). The raw MI difference is positive in every split and the corrected one is near
zero, but "explained by depth" requires the raw difference to be detected across mice, and it is
in only 3 of 10 splits of the same cells (`verdicts_v0.3.0.dev0.log`) and 5 of 40 independent
datasets (`p11b_design.log`). The test passed for the split it was written with. A verdict that
depends on what is establishable under the design is the tool's own principle; it now applies to
its tests.

**New expectations.**
- **p11a — 4 vs 4 mice:** allowed verdicts {NOT SUPPORTED, INCONCLUSIVE}; SUPPORTED is the only
  error. Tested over 10 splits.
- **p11b — N vs N mice:** NOT SUPPORTED with the diagnosis "explained by depth". N is the
  smallest number of mice per group (from 4, 6, 8, 10, 12, 16) whose lower 95% Clopper-Pearson
  bound of P(explained by depth) over 40 independent datasets is ≥ 0.80. The rule was committed
  before it was run (`p11b_design.py`, commit `cda7dd9`); the result is in `p11b_design.log`.

**Not changed.** The data generator and the truth (identical biology, capture 30%).

**Amendment D1a (2026-10-07).** No candidate up to 16 met the rule: P(explained by depth) was
5/40 at N = 4, 25/40 at 6 and 8, 26/40 at 10, 31/40 at 12 and 37/40 at 16 (lower 95% bound
0.796, just below 0.80). SUPPORTED was 0/40 at every N. The per-mouse MI from 100 cells is
noisy (between-mouse SD 0.04-0.08 against a raw difference of 0.035-0.075), so power grows
slowly with N. Rather than relax the criterion, the candidate set was extended to
{20, 24, 32} with the same criterion; the amendment was committed before those runs.

## D2 — the p01 rescaling and p10 lock-in tests (2026-10-07; decided by the project owner)

Revised when the verdict scheme changed, so that they test the v0.3 meaning instead of
conserving the v0.1.1 one: the rescaling test compares classifications and |shift| in scale
units (v0.1.1 compared statuses only); p10 keeps mi_3bin's dropout sensitivity but expects it to
be classified as attenuation and reported, not failed. Recorded here for completeness; the change
was made in `d3f05fb`.
