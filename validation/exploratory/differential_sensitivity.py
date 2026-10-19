#!/usr/bin/env python3
"""The differential sensitivity of a metric of structure on synthetic data (exploratory; the owner's decisions
of 2026-10-19, change 17 of the v0.4 spec): the numbers behind its acceptance tests and its known issues.

1. The acceptance scenarios of ``tests/test_v04_differential.py`` (mice nested in age, the old group deeper but
   less complex, as on preprint §4), for ``mi_3bin`` and ``norm_pearson``: the verdict with the check and
   without it, the groups' responses k, their ratio, the metric minus its null M, the difference at an equal
   coupling A with its interval, and the effect less the bound.
2. Where each group's response is measured: GATE 4's injection into the data (as coded) against an injection
   into GATE 0's null of the claim's genes (a common coupling, none), at the SESOI's change of the construct:
   on the data a metric that is not linear in the coupling responds more where the coupling is larger; on the
   null the response of ``mi_3bin`` is of second order.
3. Why the check needs a difference in detection: on ``tests/test_gates.py``'s clean data (the coupling
   changes the genes' variance with it, the technology is equal) the groups' responses differ by their
   biology, and the check would remove a true effect.
4. How large a real decrease must be to pass with a gap of detection as on §4, per metric.
5. How often the test of detection (genes per cell, the claim's genes left out, p < 0.05) opens at an equal
   detection, over 40 datasets.

    python validation/exploratory/differential_sensitivity.py > validation/exploratory/differential_sensitivity.log
"""
from __future__ import annotations

import platform
import subprocess
import sys
import time
from functools import partial
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO / "tests"))

import metric_autopsy  # noqa: E402
from metric_autopsy import SimpleData, injected_signal, metrics, run_autopsy  # noqa: E402
from metric_autopsy import gates as G  # noqa: E402
from metric_autopsy import report as R  # noqa: E402
from metric_autopsy.combined import Estimator, one_group_mean  # noqa: E402
from metric_autopsy.composition import technical_covariates  # noqa: E402
from metric_autopsy.injected_signal import at_construct_change  # noqa: E402
from test_v04_differential import METRICS, PREREG, SCENARIOS, deeper_less_complex  # noqa: E402

N_REP = 50


def f(x, n=4):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return str(x)
    return f"{x:+.{n}g}" if np.isfinite(x) else "—"


def ci(c, n=3):
    return f"[{f(c[0], n)}, {f(c[1], n)}]" if c else "—"


def autopsy(d, metric, check=True, prereg=PREREG, gene_pair=("A", "B")):
    saved = R.STRUCTURE_PROBES
    R.STRUCTURE_PROBES = saved if check else ()
    try:
        return run_autopsy(metric, d, group_col="age", groups=("young", "old"), replicate_col="mouse",
                           gene_pair=gene_pair, signal_test=injected_signal.coupling(*gene_pair), prereg=dict(prereg),
                           log_path="off", specificity=False)
    finally:
        R.STRUCTURE_PROBES = saved


def ds_line(a) -> str:
    d = a.effect.detail.get("differential_sensitivity") or {}
    if not d:
        return "no differential sensitivity"
    if "responses" not in d:
        return f"not sized: {d.get('error')}"
    gap = d["detection"]
    rows = "; ".join(f"{r['stratum']} {f(r['difference'], 3)} (p {r['p']:.3g})" for r in gap.get("rows", [])
                     if r.get("difference") is not None)
    adj = d.get("adjusted") or {}
    return (f"genes per cell {rows} -> differs {gap.get('differs')}; k young {f(d['responses']['young']['response'])}, "
            f"old {f(d['responses']['old']['response'])}; ratio k_{d['poor']}/k_{d['rich']} {f(d['ratio'], 3)} "
            f"{ci(d['ratio_ci'])}; M {f(d['signal'])}; A {f(d['artifact'])} {ci(d['artifact_ci'])}; tolerance "
            f"{f(d.get('tolerance'), 3)}; applies {d.get('applies')}"
            + (f"; the effect less {f(d['offset'], 3)}: {f(adj.get('effect'))} {ci(adj.get('ci'))}, p "
               f"{f(adj.get('p'), 3)}, {adj.get('status')}, passes {adj.get('passes')}" if adj else ""))


def summary(d) -> str:
    """Per group: the mean depth, genes per cell and the pair's detection (the share of cells with a count)."""
    X, grp = d.X, np.asarray(d.obs["age"])
    out = []
    for g in ("young", "old"):
        x = X[grp == g]
        out.append(f"{g} depth {x.sum(1).mean():.0f}, genes per cell {(x > 0).sum(1).mean():.1f}, detection A "
                   f"{(x[:, 0] > 0).mean():.1%} B {(x[:, 1] > 0).mean():.1%}")
    return "; ".join(out)


def section1():
    print("## 1. The acceptance scenarios (tests/test_v04_differential.py)\n")
    for metric, (fn, mu) in METRICS.items():
        for name, kw in SCENARIOS.items():
            t = time.time()
            d = deeper_less_complex(mu_claim=mu, **kw)
            print(f"- data: {metric}'s expression level {mu}, {name}: {summary(d)}")
            a, b = autopsy(d, fn), autopsy(d, fn, check=False)
            e = a.effect.detail
            print(f"- {metric}, {name} {kw}: with the check {a.verdict.split(' [')[0][:90]} [{a.cause}]; without "
                  f"{b.verdict.split(' [')[0][:60]} [{b.cause}]; the effect {f(e.get('effect'))} {ci(e.get('ci'))}, "
                  f"s_m {f(e.get('effective_sesoi'), 3)} ({time.time() - t:.0f} s)")
            print(f"  - {ds_line(a)}")
    print()


def responses(metric, d, change, on, rng, n_rep=N_REP):
    """Each group's combined response to `change` of the construct injected into its units' data or null."""
    units, _ = G._replicate_units(d, "age", ("young", "old"), "mouse", [])
    w = Estimator(units).weights
    out = {}
    for k in (0, 1):
        var = G._bio_variances(d.X, [0, 1], units, k)
        probe, _ = at_construct_change(injected_signal.coupling("A", "B"), change, var_a=var[0], var_b=var[1])
        r = np.full(units.n_units, np.nan)
        used = np.zeros(units.n_units, bool)
        for i in range(units.n_units):
            if G._unit_group(units, i) != k:
                continue
            idx = G._unit_cells(units, i, k)
            used[i] = True
            sub = SimpleData(np.array(d.X[idx], copy=True), d.obs.iloc[idx].copy(), list(d.var_names))
            base = sub if on == "data" else G.claim_null(sub, ["A", "B"], rng)
            v, _ = G._replicate_delta(metric, base, probe, n_rep, rng)
            r[i] = np.nan if v is None else v
        mu, se, _, _ = one_group_mean(r[used][None, :], units.subset(used), k, w)
        out[k] = (float(mu[0]), float(se[0]))
    return out


def section2():
    print("## 2. Where the responses are measured: the data (as coded) or GATE 0's null of the claim's genes\n")
    print(f"The change of the construct 0.1 (the SESOI), {N_REP} injections per unit; the ratio's SE by the delta "
          "method from the two groups' standard errors.\n")
    print("| metric | scenario | on | k young (SE) | k old (SE) | ratio old/young (SE) |")
    print("|---|---|---|---|---|---|")
    for metric, (fn, mu) in METRICS.items():
        for name in ("gap, equal coupling", "equal detection, decrease"):
            d = deeper_less_complex(mu_claim=mu, **SCENARIOS[name])
            for on in ("data", "null"):
                r = responses(fn, d, 0.1, on, np.random.default_rng(1))
                (ky, sy), (ko, so) = r[0], r[1]
                ratio = ko / ky
                se_r = abs(ratio) * np.sqrt((so / ko) ** 2 + (sy / ky) ** 2) if ko and ky else float("nan")
                print(f"| {metric} | {name} | {on} | {f(ky)} ({f(sy, 2)}) | {f(ko)} ({f(so, 2)}) | {f(ratio, 3)} ({f(se_r, 2)}) |")
    print()


def section3():
    from test_gates import add_mice, make_clean
    print("## 3. The clean data of tests/test_v03_verdict.py: a real difference, an equal technology\n")
    npr = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
    d = add_mice(make_clean())
    pre = {"estimand": "composition", "direction": "decrease", "judgment_pending": False}
    kw = dict(group_col="age", groups=("young", "old"), within=["sex"], replicate_col="mouse",
              gene_pair=("Smad3", "Col1a1"), signal_test=injected_signal.coupling("Smad3", "Col1a1"), prereg=pre,
              log_path="off", specificity=False)
    a = run_autopsy(npr, d, **kw)
    print(f"- with the check (the test of detection decides): {a.verdict.split(' [')[0][:100]} [{a.cause}]")
    print(f"  - {ds_line(a)}")
    saved = R.detection_gap
    R.detection_gap = lambda tech, cells: dict(saved(tech, cells), differs=True)
    try:
        b = run_autopsy(npr, d, **kw)
    finally:
        R.detection_gap = saved
    print(f"- were the test of detection open: {b.verdict.split(' [')[0][:160]} [{b.cause}]")
    print("  (the coupling of these data changes with the genes' variance, 2.6 against 0.3 loadings of one factor, so "
          "the metric's reliability differs between the groups by their biology, at an equal technology)\n")


def section4():
    print("## 4. A real decrease with a gap of detection as on §4: what passes\n")
    for metric, (fn, mu) in METRICS.items():
        for rho in ((0.6, 0.3), (0.7, 0.2), (0.8, 0.1)):
            d = deeper_less_complex(mu_claim=mu, rho=rho)
            a = autopsy(d, fn)
            print(f"- {metric}, coupling {rho[0]} against {rho[1]} (a decrease of {rho[0] - rho[1]:.1f}, "
                  f"{(rho[0] - rho[1]) / 0.1:.0f} x SESOI): {a.verdict.split(' [')[0][:80]} [{a.cause}]")
            print(f"  - {ds_line(a)}")
    print()


def section5(n=40):
    print("## 5. The test of detection at an equal detection: how often it opens\n")
    opened = 0
    for s in range(n):
        d = deeper_less_complex(seed=1000 + s, **SCENARIOS["equal detection, equal coupling"])
        tech = technical_covariates(d, "age", ("young", "old"), "mouse", alpha=0.05,
                                    rng=np.random.default_rng(s), exclude=["A", "B"])
        opened += R.detection_gap(tech, "all cells")["differs"]
    print(f"- {opened} of {n} datasets (8 against 8 mice, equal technology): the test's own rate, alpha 0.05\n")


def main() -> int:
    rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "--", str(REPO / "src")], cwd=HERE, capture_output=True,
                           text=True).stdout.strip()
    print(f"# The differential sensitivity on synthetic data — git {rev}{' (with uncommitted changes to src)' if dirty else ''}, "
          f"metric-autopsy {metric_autopsy.__version__}, python {platform.python_version()}, numpy {np.__version__}\n")
    t = time.time()
    section1()
    section2()
    section3()
    section4()
    section5()
    print(f"{time.time() - t:.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
