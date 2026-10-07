"""The "after" numbers: every probe's truth run through the v0.3 API (four-field verdict).

The pNN_*.py scripts call the gates with the v0.1.1 signatures (no replicate unit, no
estimand), so their verdict lines cannot show the v0.3 verdict on each truth. This script
re-runs each truth exactly as `test_probes.py` does and prints the fields, plus the null
rates and the robustness checks behind them. Deterministic; ~8 min.

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

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import test_probes as tp  # noqa: E402  (builders and the frozen calls)
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
        signal_test=injected_signal.module(tp.G2M_GENES, fold=2.0, frac=0.3), prereg=tp.RESOLVED))

    print("=== p03 proliferation 35% -> 5% cycling, 3 mice per sex x age")
    a = run_autopsy(tp.probe_sim.mean_g2m_score, tp._proliferation(), within=["sex"],
                    replicate_col="mouse",
                    signal_test=injected_signal.module(tp.G2M_GENES, fold=2.0, frac=0.3),
                    prereg=tp.RESOLVED, **AGE)
    show(f"mean G2M score, retained {a.effect.detail['retained']:.0%} after depth thinning", a)

    print("=== p04 Xist female > male on the demo data (male-old capture 30%)")
    a = run_autopsy(tp._mean_lognorm_xist, tp._xist_demo(), group_col="sex", groups=("female", "male"),
                    within=["age"], replicate_col="mouse",
                    signal_test=injected_signal.module(["Xist"], fold=2.0, frac=0.3), prereg=tp.RESOLVED)
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


if __name__ == "__main__":
    main()
