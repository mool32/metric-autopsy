"""The "after" numbers: every probe's truth run through the v0.3 API (four-field verdict).

The pNN_*.py scripts call the gates with the v0.1.1 signatures (no replicate unit, no
estimand), so their verdict lines cannot show the v0.3 verdict on each truth. This script
re-runs each truth exactly as `test_probes.py` does and prints the fields, plus the null
rates and the robustness checks behind them.

The last section has two axes per case, each a rate with its 95% Clopper-Pearson interval over
independent datasets:

* **errors** — verdicts outside the set that is correct *given the design* (e.g. INCONCLUSIVE
  is allowed for a pure depth artifact with 4 vs 4 mice, SUPPORTED never is), and false
  SUPPORTED on nulls, artifacts and useless metrics;
* **decisiveness** — on cases whose truth is establishable under their design, the share of
  definite verdicts (anything but INCONCLUSIVE; an UNTESTED metric makes any effect claim
  INCONCLUSIVE, so "not UNTESTED" is included) and of correct definite verdicts.

A validator that never errs because it never commits is useless; both axes are reported.
Deterministic; ~25 min.

    python validation/probes/verdicts_v03.py > validation/probes/verdicts_v0.3.0.dev0.log

Development set: these probes found the bugs and the fixes were developed against them, so
the numbers below have no confirmatory weight.
"""
from __future__ import annotations

import platform
import subprocess
import sys
from functools import partial
from pathlib import Path

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import test_probes as tp  # noqa: E402  (builders and the frozen calls)
from test_gates import add_mice, make_clean  # noqa: E402
from test_v03_verdict import (  # noqa: E402
    blind_pair_metric, capture_confound_mice, log_total_endogenous, mean_lognorm_gene5, with_ercc,
)
from metric_autopsy.core import SimpleData  # noqa: E402
from metric_autopsy.stats import fmt_rate  # noqa: E402
from metric_autopsy import (  # noqa: E402
    GateStatus, __version__, gate0_independence, gate1_qc_parity, gate5_controls, injected_signal,
    metrics, run_autopsy,
)
from metric_autopsy.cli import demo_data  # noqa: E402

AGE = dict(group_col="age", groups=("young", "old"))


def fields(a) -> str:
    return " | ".join(f"{k}={v.status}" + (f"[{','.join(v.flags)}]" if v.flags else "")
                      for k, v in a.fields().items())


def show(label, a):
    print(f"  {label}")
    print(f"    {fields(a)}")
    print(f"    VERDICT: {a.verdict}")


def header():
    try:
        rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, capture_output=True,
                             text=True).stdout.strip() or "?"
    except Exception:
        rev = "?"
    print(f"# metric-autopsy {__version__} — v0.3 verdicts on the probe truths — git {rev}, "
          f"python {platform.python_version()}, numpy {np.__version__}")


def main():
    header()

    print("=== p01 useless, constant and offset metrics")
    rng = np.random.default_rng(0)
    show("random-number metric (make_clean)",
         run_autopsy(lambda data: float(rng.normal()), tp.make_clean(), within=["sex"], **AGE))
    show("constant metric", run_autopsy(lambda data: 0.0, tp.make_clean(), **AGE))
    d = tp.make_clean()
    base = gate0_independence(tp.MI, d)
    for c in (0.5, 10.0):
        g = gate0_independence(lambda data, c=c: tp.MI(data) + c, d)
        same = all(abs(g.detail["responses"][k]["shift_std"] - r["shift_std"]) < 1e-9
                   for k, r in base.detail["responses"].items())
        print(f"  GATE 0 on mi_3bin + {c:g}: {g.status.value} (same as mi_3bin: {g.status == base.status}; "
              f"identical shifts: {same})")

    print("=== p02 sorted G2M vs G1 (2.4x RNA in G2M), plates as replicates")
    show("mean G2M score + injected module", run_autopsy(
        tp.probe_sim.mean_g2m_score, tp._sorted_cell_cycle(), group_col="sorted_phase",
        groups=("G2M", "G1"), replicate_col="plate",
        signal_test=injected_signal.module(tp.G2M_GENES, fold=2.0, frac=0.3), prereg=tp.TRUE_DECREASE))

    print("=== p03 proliferation 35% -> 5% cycling, 3 mice per sex x age")
    a = run_autopsy(tp.probe_sim.mean_g2m_score, tp._proliferation(), within=["sex"],
                    replicate_col="mouse",
                    signal_test=injected_signal.module(tp.G2M_GENES, fold=2.0, frac=0.3),
                    prereg=tp.TRUE_DECREASE, **AGE)
    show(f"mean G2M score, retained {a.effect.detail['retained']:.0%} after depth thinning", a)

    print("=== p04 Xist female > male on the demo data (male-old capture 30%)")
    a = run_autopsy(tp._mean_lognorm_xist, tp._xist_demo(), group_col="sex", groups=("female", "male"),
                    within=["age"], replicate_col="mouse",
                    signal_test=injected_signal.module(["Xist"], fold=2.0, frac=0.3), prereg=tp.TRUE_DECREASE)
    show(f"mean log Xist, retained {a.effect.detail['retained']:.0%}", a)

    print("=== p05 null strata: flag / FAIL rates")
    flagged = sum(gate1_qc_parity(tp._null_qc_strata(64, 20, s), "age", ("young", "old"),
                                  within=["stratum"]).status != GateStatus.PASS for s in range(30))
    print(f"  GATE 1, 64 strata x 20 cells: {flagged}/30 null datasets not PASS (v0.1.1: 82%)")
    for name, k, n, n_sets in (("norm_pearson", 4, 400, 20), ("norm_pearson", 4, 30, 20),
                               ("norm_pearson", 1, 10, 20), ("norm_pearson", 16, 100, 20),
                               ("mi_3bin", 16, 100, 20), ("mi_3bin", 16, 400, 20)):
        res = [gate5_controls(getattr(metrics, name), tp._null_controls(k, n, s), ("Actb", "Gapdh"),
                              ("Gene0", "Gene1"), within=["stratum"]) for s in range(n_sets)]
        fails = sum(r.status == GateStatus.FAIL for r in res)
        warns = sum(r.status == GateStatus.WARN for r in res)
        print(f"  GATE 5 {name}, {k} x {n} cells: FAIL {fails}/{n_sets}, WARN {warns}/{n_sets} "
              f"(positive control not demonstrated in some stratum)")

    print("=== p06 negative control under CP10k closure (demo data)")
    r = gate5_controls(metrics.norm_pearson, demo_data(), ("Actb", "Gapdh"), ("Gene0", "Gene1"), within=["sex"])
    centres = ", ".join(f"{row['stratum']['sex']} {row['null_center']:+.3f}" for row in r.detail["rows"])
    print(f"  GATE 5: {r.status.value}; null centre of unrelated pairs: {centres}")

    print("=== p07 pseudoreplication: no age effect, mouse-level variation")
    det3 = sum(run_autopsy(tp.NPR, tp._mice(3, mouse_sd=0.35, seed=s), replicate_col="mouse",
                           prereg=tp.COMPOSITION, **AGE).effect.status == "DETECTED" for s in range(20))
    cells = {run_autopsy(tp.NPR, tp._mice(3, mouse_sd=0.35, seed=s), prereg=tp.COMPOSITION,
                         **AGE).design_adequacy.status for s in range(20)}
    det6 = sum(run_autopsy(tp.NPR, tp._mice(6, mouse_sd=0.35, seed=s), replicate_col="mouse",
                           prereg=tp.COMPOSITION, **AGE).effect.status == "DETECTED" for s in range(30))
    print(f"  3 vs 3 mice: DETECTED {det3}/20 (parametric only); without replicate_col: {sorted(cells)} "
          f"(v0.1.1: 22/40 'effect survives')")
    print(f"  6 vs 6 mice: DETECTED {det6}/30 (exact permutation over mice)")

    print("=== p08 attenuation of a truly coupled pair")
    dcp, ga, gb = tp._coupled_pair()
    g0 = gate0_independence(partial(tp._lognorm_pearson, gene_a=ga, gene_b=gb), dcp, protect_genes=(ga, gb))
    resp = g0.detail["responses"]
    print(f"  GATE 0: {g0.status.value}; " + ", ".join(
        f"{k} {v['classification']}" + (f" −{v['signal_loss']:.0%}" if "signal_loss" in v else "")
        for k, v in resp.items()))

    print("=== p09 equivalence needs a SESOI; stability across tool seeds")
    real = [run_autopsy(tp.NPR, tp._mice(6, c_old=1.2, seed=s), replicate_col="mouse",
                        prereg=tp.COMPOSITION, **AGE).effect.status for s in range(3)]
    null = tp._mice(6, seed=0)
    no_s = run_autopsy(tp.NPR, null, replicate_col="mouse", prereg=tp.COMPOSITION, **AGE).effect.status
    with_s = run_autopsy(tp.NPR, null, replicate_col="mouse", prereg={**tp.COMPOSITION, "sesoi": 0.15},
                         **AGE).effect.status
    seeds = sorted({run_autopsy(tp.NPR, tp._mice(6, c_old=1.2, seed=0), replicate_col="mouse",
                                prereg=tp.COMPOSITION, seed=s, **AGE).effect.status for s in range(4)})
    print(f"  coupling 1.5 vs 1.2, 3 datasets: {real}")
    print(f"  null data: no SESOI -> {no_s}; SESOI 0.15 -> {with_s}")
    print(f"  one dataset, tool seeds 0-3: {seeds}")

    print("=== p10 mi_3bin sensitivity is reported as attenuation (demo data, default perturbations)")
    g0 = gate0_independence(tp.MI, demo_data())
    print(f"  GATE 0: {g0.status.value}; attenuation " + ", ".join(
        f"{k} −{v:.0%}" for k, v in g0.detail["attenuation"].items()))

    print("=== p11 pure depth artifact (demo male stratum, 4 vs 4 mice): robustness over mouse splits")
    d = demo_data()
    males = d[np.asarray(d.obs["sex"]) == "male"]
    verdicts = []
    for split in range(10):
        a = run_autopsy(tp.MI, tp._with_mice(males, 4, ["age"], seed=split), replicate_col="mouse",
                        prereg=tp.COMPOSITION, **AGE)
        e = a.effect.detail
        verdicts.append(a.verdict.split(" — ")[0])
        print(f"  split {split}: raw {e['raw_effect']:+.4f} (p={e['p_raw']:.3f}) -> at equal depth "
              f"{e['effect']:+.4f} ({e['retained']:+.0%} retained): {a.verdict.split(' — ')[0]}")
    print(f"  NOT SUPPORTED in {verdicts.count('NOT SUPPORTED')}/10 splits, otherwise INCONCLUSIVE; "
          "never SUPPORTED: " + str("SUPPORTED" not in verdicts and
                                    not any(v.startswith("SUPPORTED") for v in verdicts)))

    print("=== demo (CLI defaults): mi_3bin and norm_pearson")
    for name in ("mi_3bin", "norm_pearson"):
        fn = getattr(metrics, name)
        show(name, run_autopsy(partial(fn, gene_a="Smad3", gene_b="Col1a1"), demo_data(), within=["sex"],
                               gene_pair=("Smad3", "Col1a1"), replicate_col="mouse", pair_metric=fn,
                               pos_pair=("Actb", "Gapdh"), neg_pair=("Gene0", "Gene1"),
                               prereg=tp.COMPOSITION, stop_on_first_fail=False, **AGE))


# --------------------------------------------------------------------------- #
# errors and decisiveness, case by case
# --------------------------------------------------------------------------- #
P11B_MICE = 20  # p11b_design.log (journal D1a/D1b)


@dataclass
class Case:
    name: str
    make: Callable  # seed -> (metric, data, run_autopsy kwargs)
    allowed: tuple  # verdict prefixes that are correct given the design
    definite: Optional[tuple] = None  # expected definite verdicts if the truth is establishable
    null: bool = False  # a null, an artifact or a useless metric: SUPPORTED is a false positive
    n: int = 10


def _matches(verdict: str, prefixes) -> bool:
    return any(verdict.startswith(p) for p in prefixes)


def _xist_demo_seeded(seed: int, mice: int = 4):
    d = demo_data(seed=seed)
    rng = np.random.default_rng(seed + 7)
    female = np.asarray(d.obs["sex"]) == "female"
    old_male = (~female) & (np.asarray(d.obs["age"]) == "old")
    xist = rng.poisson(np.where(female, 25.0, 0.0) * np.where(old_male, 0.30, 1.0)).astype(float)
    dx = SimpleData(np.column_stack([d.X, xist]), d.obs, list(d.var_names) + ["Xist"])
    return tp._with_mice(dx, mice, ["sex", "age"], seed=seed)


def _demo_males(seed: int, n_mice: int):
    d = demo_data(seed=seed, n=100 * n_mice, mice_per_block=n_mice)
    return d[np.asarray(d.obs["sex"]) == "male"]


def _random_metric(seed):
    rng = np.random.default_rng(seed)
    return lambda data: float(rng.normal())


SIGNAL = injected_signal.coupling("Smad3", "Col1a1")
G2M = injected_signal.module(tp.G2M_GENES, fold=2.0, frac=0.3)
# JOURNAL.md, D3: real effects declare their true direction; nulls and artifacts declare
# 'two-sided', so a detected effect in either direction could become SUPPORTED
ANY = {**tp.RESOLVED, "direction": "two-sided"}


def cases():
    out = [
        Case("p01 random-number metric, full design (mice, estimand, injected signal)",
             lambda s: (_random_metric(s), add_mice(make_clean(seed=s)),
                        dict(within=["sex"], replicate_col="mouse", signal_test=SIGNAL, prereg=ANY, **AGE)),
             allowed=("NOT SUPPORTED", "INCONCLUSIVE"), definite=("NOT SUPPORTED",), null=True),
        Case("p01 constant metric",
             lambda s: (lambda data: 0.0, add_mice(make_clean(seed=s)),
                        dict(replicate_col="mouse", prereg=ANY, **AGE)),
             allowed=("DEGENERATE METRIC",), definite=("DEGENERATE METRIC",), null=True),
        Case("blind pair metric (ignores gene b), its own controls, 800 cells per stratum",
             lambda s: (partial(blind_pair_metric, gene_a="Smad3", gene_b="Col1a1"), add_mice(make_clean(seed=s)),
                        dict(within=["sex"], replicate_col="mouse", gene_pair=("Smad3", "Col1a1"),
                             pair_metric=blind_pair_metric, pos_pair=("Actb", "Gapdh"),
                             neg_pair=("Gene0", "Gene1"), prereg=ANY, **AGE)),
             allowed=("NOT SUPPORTED", "INCONCLUSIVE"), definite=("NOT SUPPORTED",), null=True),
        Case("p02 sorted G2M vs G1, 4 plates each (real, huge)",
             lambda s: (tp.probe_sim.mean_g2m_score, tp._sorted_cell_cycle(seed=s + 1),
                        dict(group_col="sorted_phase", groups=("G2M", "G1"), replicate_col="plate",
                             signal_test=G2M, prereg=tp.TRUE_DECREASE)),
             allowed=("SUPPORTED", "INCONCLUSIVE"), definite=("SUPPORTED",)),
        Case("p03 proliferation 35% -> 5%, 3 mice per sex x age (real)",
             lambda s: (tp.probe_sim.mean_g2m_score, tp._proliferation(seed=s + 2),
                        dict(within=["sex"], replicate_col="mouse", signal_test=G2M, prereg=tp.TRUE_DECREASE, **AGE)),
             allowed=("SUPPORTED", "INCONCLUSIVE"), definite=("SUPPORTED",)),
        Case("p04 Xist female > male, demo data, 4 mice per block (real, QC gap)",
             lambda s: (tp._mean_lognorm_xist, _xist_demo_seeded(s),
                        dict(group_col="sex", groups=("female", "male"), within=["age"], replicate_col="mouse",
                             signal_test=injected_signal.module(["Xist"], fold=2.0, frac=0.3), prereg=tp.TRUE_DECREASE)),
             allowed=("SUPPORTED", "INCONCLUSIVE"), definite=("SUPPORTED",)),
        Case("p07 no effect, 3 vs 3 mice with mouse variance (not establishable)",
             lambda s: (tp.NPR, tp._mice(3, mouse_sd=0.35, seed=s), dict(replicate_col="mouse", prereg=tp.COMPOSITION, **AGE)),
             allowed=("INCONCLUSIVE", "NO DETECTABLE EFFECT", "NOT SUPPORTED"), null=True, n=20),
        Case("p09 no effect, 6 vs 6, SESOI 0.15, injected signal (establishable absence)",
             lambda s: (tp.NPR, tp._mice(6, seed=s),
                        dict(replicate_col="mouse", signal_test=SIGNAL,
                             prereg={**ANY, "sesoi": 0.15}, **AGE)),
             allowed=("NO DETECTABLE EFFECT", "INCONCLUSIVE", "NOT SUPPORTED"), definite=("NO DETECTABLE EFFECT",),
             null=True, n=20),
        Case("p09 moderate effect (coupling 1.5 vs 1.2), 6 vs 6, injected signal (real)",
             lambda s: (tp.NPR, tp._mice(6, c_old=1.2, seed=s),
                        dict(replicate_col="mouse", signal_test=SIGNAL, prereg=tp.TRUE_DECREASE, **AGE)),
             allowed=("SUPPORTED", "INCONCLUSIVE"), definite=("SUPPORTED",)),
        Case("p11a pure depth artifact, 4 vs 4 mice (not establishable)",
             lambda s: (tp.MI, _demo_males(s, 4), dict(replicate_col="mouse", prereg=tp.COMPOSITION, **AGE)),
             allowed=("NOT SUPPORTED", "INCONCLUSIVE"), null=True, n=20),
        Case("level metric on a pure depth artifact, 6 vs 6, injected module",
             lambda s: (mean_lognorm_gene5, capture_confound_mice(seed=s),
                        dict(replicate_col="mouse", signal_test=injected_signal.module(["Gene5"], fold=2.0, frac=0.3),
                             prereg=ANY, **AGE)),
             allowed=("NOT SUPPORTED", "INCONCLUSIVE"), definite=("NOT SUPPORTED",), null=True, n=20),
        Case("content estimand, capture halved, ERCC present (artifact)",
             lambda s: (log_total_endogenous, with_ercc(capture=(1.0, 0.5), seed=s),
                        dict(replicate_col="mouse", prereg={"estimand": "content"}, **AGE)),
             allowed=("NOT SUPPORTED", "INCONCLUSIVE"), definite=("NOT SUPPORTED",), null=True),
        Case("content estimand, capture halved, no spike-ins",
             lambda s: (log_total_endogenous, with_ercc(capture=(1.0, 0.5), spikes=False, seed=s),
                        dict(replicate_col="mouse", prereg={"estimand": "content"}, **AGE)),
             allowed=("UNIDENTIFIABLE",), definite=("UNIDENTIFIABLE",), null=True),
    ]
    if P11B_MICE:
        out.insert(10, Case(f"p11b pure depth artifact, {P11B_MICE} vs {P11B_MICE} mice (establishable)",
                            lambda s: (tp.MI, _demo_males(s, P11B_MICE),
                                       dict(replicate_col="mouse", prereg=tp.COMPOSITION, **AGE)),
                            allowed=("NOT SUPPORTED", "INCONCLUSIVE"), definite=("NOT SUPPORTED",), null=True, n=20))
    return out


def decisiveness():
    print("=== errors and decisiveness (95% Clopper-Pearson; independent datasets per case)")
    tot = dict(n=0, out=0, null_n=0, false_sup=0, est_n=0, definite=0, correct=0)
    for c in cases():
        verdicts = []
        for s in range(c.n):
            metric, data, kw = c.make(s)
            verdicts.append(run_autopsy(metric, data, log_path="off", **kw).verdict)
        out = sum(not _matches(v, c.allowed) for v in verdicts)
        line = f"  {c.name}\n    outside the allowed set: {fmt_rate(out, c.n)}"
        tot["n"] += c.n
        tot["out"] += out
        if c.null:
            fs = sum(v.startswith("SUPPORTED") for v in verdicts)
            line += f"; false SUPPORTED: {fmt_rate(fs, c.n)}"
            tot["null_n"] += c.n
            tot["false_sup"] += fs
        if c.definite:
            dfn = sum(not v.startswith("INCONCLUSIVE") for v in verdicts)
            cor = sum(_matches(v, c.definite) for v in verdicts)
            line += f"; definite: {fmt_rate(dfn, c.n)}; correct definite: {fmt_rate(cor, c.n)}"
            tot["est_n"] += c.n
            tot["definite"] += dfn
            tot["correct"] += cor
        else:
            line += "; not establishable under this design (decisiveness not scored)"
        kinds = sorted({v.split(" — ")[0].split(" [")[0] for v in verdicts})
        print(line + f"\n    verdicts seen: {kinds}", flush=True)
    print("  POOLED (dev set):")
    print(f"    outside the allowed set: {fmt_rate(tot['out'], tot['n'])}")
    print(f"    false SUPPORTED on nulls, artifacts and useless metrics: {fmt_rate(tot['false_sup'], tot['null_n'])}")
    print(f"    decisiveness on establishable cases: {fmt_rate(tot['definite'], tot['est_n'])}; "
          f"correct definite: {fmt_rate(tot['correct'], tot['est_n'])}")


if __name__ == "__main__":
    main()
    decisiveness()
