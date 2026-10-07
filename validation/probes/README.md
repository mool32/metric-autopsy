# Validator probes: exploratory / development set

> **Status: exploratory dev set. These probes have no confirmatory weight.**
> They found the bugs listed below, and the fixes are developed against them, so
> passing them after a fix proves nothing about the validator's operating
> characteristics. Confirmatory validation needs a **frozen code tag** and a **new
> panel** with no overlap with this directory, assessed blind (plan step 3).

Each probe gives the gates a case whose truth is known by construction and checks
whether the verdict matches it. `baseline_v0.1.1.log` is the output against the
v0.1.1 engine, i.e. the "before" numbers.

```bash
./run_all.sh                 # ~3 min; prints to stdout
WITH_MEMORY=1 ./run_all.sh   # also runs p12 (peaks at ~6 GB RAM)
```

`test_probes.py` turns every probe into a regression test that asserts the
**correct expected verdict**. That verdict comes from the truth, never from what
the engine currently prints. Tests marked `xfail(strict=True)` document a known
failure. When a fix makes one pass, strict mode reports XPASS as an error until
the marker is removed, so every flip shows up in a diff. Expected verdicts are
frozen when the test is written; a fix must meet them without editing them.

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
