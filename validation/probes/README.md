# Validator probes: exploratory / development set

> **Status: exploratory dev set. These probes have no confirmatory weight.**
> They found the bugs listed below, and the fixes are developed against them, so
> passing them after a fix proves nothing about the validator's operating
> characteristics. Confirmatory validation needs a **frozen code tag** and a **new
> panel** with no overlap with this directory, assessed blind (plan step 3).

Each probe gives the gates a case whose truth is known by construction and checks
whether the verdict matches it. `baseline_v0.1.1.log` is the output against the
v0.1.1 engine, i.e. the "before" numbers. `after_v0.3.0.dev0.log` is the same scripts
against this branch. The scripts call the gates with the v0.1.1 signatures (no replicate
unit, no estimand), so their gate-level numbers compare directly, while their verdict lines
mostly read UNIDENTIFIABLE or INSUFFICIENT_REPLICATION by design. `verdicts_v03.py` re-runs
every truth through the v0.3 API, as the tests do; its output is
`verdicts_v0.3.0.dev0.log`.

```bash
./run_all.sh                 # ~15 min on this branch; prints to stdout
WITH_MEMORY=1 ./run_all.sh   # also runs p12 (peaks at ~6 GB RAM)
python verdicts_v03.py       # ~8 min; the v0.3 verdict on every truth
```

`test_probes.py` turns every probe into a regression test that asserts the
**correct expected verdict**. That verdict comes from the truth, never from what
the engine currently prints. Tests marked `xfail(strict=True)` documented the known
failures of v0.1.1. When a fix made one pass, strict mode reported XPASS as an error
until the marker was removed, so every flip shows up in a diff. Expected verdicts are
frozen when the test is written; a fix must meet them without editing them. All 18
markers have been removed in v0.3.0.dev0 (each replaced by a `# fixed in` comment).

`probe_common.py` finds the checkout by walking up to `src/metric_autopsy` and
reuses the 40-gene generators from `tests/test_gates.py`. `probe_sim.py` is a
1,500-gene negative-binomial simulator with 10x-like depth and a cell-cycle
module, in which cycling cells can carry more RNA. The 40-gene generator plants
exactly the confound that GATE 1 measures, so it cannot test specificity.

## What the probes found (v0.1.1)

| Probe | Case (truth) | v0.1.1 verdict / result | Mechanism |
|---|---|---|---|
| p01 | metric returns random numbers (useless) | **PASS — cleared 3 auto gates** | GATE 0 tests invariance only; "no effect" counts as PASS; the Python API leaves judgment unset |
| p01 | `mi_3bin + 0.5` (same confounded metric) | GATE 0 PASS (rel 62% → 10%) | shift is divided by \|baseline\| (`gates.py:132`), so the verdict depends on the metric's location |
| p01 | constant metric (degenerate) | FAIL at GATE 2, "likely a QC artifact" | `0 < 0` is false in the null-floor check; wrong diagnosis |
| p02 | sorted G1 vs G2M, normalized G2M score (real, huge) | PASS at G2M RNA content ×1.0–1.9; **FAIL at GATE 1** at ×2.4 while GATE 2 retains 99% | GATE 1 is data-level and terminal |
| p03 | proliferation 35% → 5% cycling; cycling cells carry more RNA (real) | ×1.9: PASS, but matching removes 32–41% of the effect; ×2.4: **FAIL at GATE 2**, "likely a QC artifact" (27–30% retained) | n_genes is downstream of the biology, so matching on it is a bad control |
| p04 | Xist female > male on the demo data (certain) | **FAIL at GATE 1** (old stratum 2.00×); GATE 2 retains 100% | §4's real TMS table gives 3413/1701 = 2.01× in the old stratum |
| p05 | null data, no confound, K strata | GATE 1 false-FAIL 28% (16 strata × 20 cells) and 82% (64 × 20); GATE 5 (norm_pearson) false-FAIL 22% pooled at 400 cells, 70% at 4 strata, 98% at 16 | no min_cells, no sampling uncertainty, no multiplicity control, one negative pair |
| p06 | independent genes under CP10k normalization | mean norm_pearson r = +0.084 (raw co-detected Spearman −0.003) | compositional closure: spurious correlation of ratios |
| p07 | 3 vs 3 mice, no age effect, mouse-level variation | 22/40 runs: "effect survives matching"; replicate-level exact test: 0/40 (min p = 0.10) | the permutation unit is the cell, not the mouse |
| p08 | truly coupled pair, standard log-normalized Pearson | GATE 0 FAIL at every depth (0.3–10 UMI/cell) | attenuation (reliability) is scored as confounding; extra_dropout zeroes large UMI counts |
| p09 | moderate true effect (coupling 1.5 vs 1.2) | 21/40 "survives", 19/40 "the groups simply do not differ" | null floor = 95th pct of 20 permutations (0.047–0.123 across seeds) |
| p10 | mi_3bin under smaller and larger perturbations | FAIL everywhere (27% shift at only 5% extra dropout) | **in the tool's favour**: the headline sensitivity of mi_3bin is robust |
| p11 | depth equalization by binomial thinning vs n_genes cell selection | proliferation 100% retained (selection 34%, FAIL); Xist 99% (selection STOP); mi_3bin artifact 5% (correctly removed) | thinning removes depth without selecting cells on a biology-affected variable; valid only for a *composition* estimand |
| p12 | 5,000 × 20,000 sparse AnnData (40 MB nnz) | GATE 0 peak RSS 5.6 GB for 14 metric evaluations (the default is 120) | the whole matrix is densified and copied per bootstrap/perturbation; TMS FACS ≈ 20 GB per copy |

Caveat on p11: thinning to equal depth is correct only when the pre-registered
estimand is about **composition** (relative expression). For a **content** estimand
(total RNA, number of detected genes, e.g. CytoTRACE), thinning destroys the
signal itself. The correction must follow from the estimand declared in the
pre-registration: thinning for composition; an external standard (ERCC) or an
UNIDENTIFIABLE verdict for content.

## After the rework (v0.3.0.dev0, development result)

Source: `verdicts_v0.3.0.dev0.log` (`verdicts_v03.py`). The fixes were developed against
these probes, so none of this is confirmatory.

| Probe | Truth | v0.1.1 | v0.3.0.dev0 |
|---|---|---|---|
| p01 | random-number metric | PASS | metric UNTESTED; never SUPPORTED |
| p01 | constant metric | FAIL at GATE 2, "QC artifact" | DEGENERATE METRIC |
| p01 | `mi_3bin + C`, `k · mi_3bin` | GATE 0 status depends on C | identical classification and shifts |
| p02 | sorted G2M vs G1, 2.4× RNA | FAIL at GATE 1 | SUPPORTED (provisional), QC gap corrected |
| p03 | proliferation 35% → 5% | FAIL at GATE 2 (27–30% retained) | SUPPORTED (provisional), 101% retained after thinning |
| p04 | Xist female > male, demo data | FAIL at GATE 1 | SUPPORTED (provisional), 99% retained, CORRECTED |
| p05 | null strata, 64 × 20 cells | GATE 1 flags 82% | 0/30 |
| p05 | null controls, 4 × 400 cells | GATE 5 fails 70% | FAIL 1/20 (only the negative control can fail) |
| p05 | null controls, 16 strata, `mi_3bin` | GATE 5 fails 100 / 100 / 83 / 0% at 10 / 30 / 100 / 400 cells | FAIL 0/20 at 100 and at 400 cells (`after_v0.3.0.dev0.log`: 3 / 0 / 0 / 2%) |
| p06 | unrelated genes under CP10k closure | negative control fails | PASS; the null centre (+0.09 female, +0.29 male) is reported |
| p07 | no age effect, 3 vs 3 mice | 22/40 "effect survives" | 0/20 detected (parametric only); no verdict without a replicate unit |
| p07 | no age effect, 6 vs 6 mice | — | 1/30 detected (exact permutation over mice) |
| p08 | truly coupled pair, log-normalized Pearson | GATE 0 FAIL | attenuation reported (dropout −50%, depth −20%), PASS |
| p09 | coupling 1.5 vs 1.2 | 19/40 "groups do not differ" | DETECTED 3/3; on null data NO DETECTABLE EFFECT only with a SESOI; stable over tool seeds |
| p10 | mi_3bin at 5% extra dropout | FAIL (lock-in) | sensitivity kept, now classified as attenuation |
| p11 | pure depth artifact, 4 vs 4 mice | — (no thinning) | NOT SUPPORTED, explained by depth, in 3/10 mouse splits; INCONCLUSIVE in 7/10; never SUPPORTED |
| p12 | 5,000 × 20,000 sparse matrix | silent densification | `DenseMemoryWarning` before densifying |

Found while checking the robustness of these results, and fixed:

- **GATE 5 in small strata.** Before the fix the reworked GATE 5 still failed 100% of null
  datasets at 10 cells per stratum and 97% at 4 × 30: a one-cell depth bin cannot be shuffled,
  so the positive control's null equalled the data, and a silent positive control in an
  underpowered stratum counted as FAIL. Depth bins now hold ≥ 10 cells, and only the negative
  control can fail the gate; a silent positive control is WARN (absence of evidence).
  After: FAIL 0/20 at 1 × 10 and 0/20 at 4 × 30 cells.
- **GATE 5 with a skewed metric and many strata.** At 16 strata alpha/K = 0.003 is below the
  resolution of 200 null draws. The engine then assumed a normal tail, which is wrong for MI's
  right-skewed null: `mi_3bin` failed 17–33% of null datasets (`run_all.sh` at commit
  `96d4e53`, p05). Extending the null instead exposed a second flaw: the unrelated pairs were
  drawn with replacement from a few hundred candidates in this 40-gene panel, so an extended
  null repeated pairs and overstated its resolution. Pairs are now distinct, the null is
  extended only as far as distinct pairs allow, and no parametric tail is assumed (after: the
  16-strata rows of `verdicts_v0.3.0.dev0.log` and the p05 rows of `after_v0.3.0.dev0.log`).
- **GATE 0 on level metrics.** The gene shuffle removes cell-size covariance and moves a mean
  log total by ~0.5%. GATE 0 took that as the metric's "signal" and judged depth shifts in its
  units (30–137×), failing a content metric even on null data. A null signal smaller than
  `tol` × the effect scale no longer makes a metric a structure metric.

Known fragility: **p11.** With 4 vs 4 mice the exact permutation has 70 assignments, so a
raw difference reaches p < 0.05 only at complete separation of the mice (p = 0.029); the next
attainable values are 0.057 and 0.086. The raw MI difference is positive in all 10 splits and
the corrected one is near zero (−8% to +17% retained), but "explained by depth" needs the
raw difference to be detected across mice, so the frozen expectation (NOT SUPPORTED) holds in
3 of 10 splits; the test passes for the split it was written with. A per-mouse test of the
depth component (raw − corrected) was tried and rejected: the untouched group's component is
zero by construction, which makes that permutation test anti-conservative (complete
separation under the null with probability 1/8).
