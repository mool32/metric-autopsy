#!/usr/bin/env python3
"""The acceptance of change 17' on synthetic data (the owner's second decisions of 2026-10-19; v0.4 spec,
section 17'): the detection of a pair's genes equalized after equal depth, with GATE 4's s_m, the effect and
the line of norm_pearson on the equalized data.

The criteria below were fixed in this script, committed before its run, and are printed at the head of its
log before any result. The generator is `tests/detection_data.py` (that of `differential_sensitivity.py`),
20 seeds per scenario, `mi_3bin` and `norm_pearson`; "v0.4 without change 17" is the engine of 1337ead (the
commit before change 17), run from a worktree on the same data; "the new engine" is this checkout's src.
Changes 18 and 19, which the owner kept, change only the dossier (the sign of its rows, a line on the
correction), so 1337ead's verdicts are those of v0.4 without change 17.

    git worktree add <dir> 1337ead
    python validation/exploratory/detection_equalization.py --base-src <dir>/src --workers 4 \\
        > validation/exploratory/detection_equalization.log
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
SEEDS = list(range(20))
METRIC_NAMES = ("mi_3bin", "norm_pearson")
GAP = dict(dom=(0.1, 0.75), depth=(3000.0, 4200.0))  # the old group deeper but less complex, as on §4
EQUAL = dict(dom=(0.1, 0.1), depth=(3000.0, 3000.0))
COUPLINGS = {3: (0.6, 0.3), 5: (0.7, 0.2), 7: (0.8, 0.1)}  # decreases of 3, 5 and 7 x SESOI (0.1)
# (vi): one group's cells vary more in capture (sd of log capture 0.9 against 0.3), the mean depth equal; a
# wider spread lowers the group's detection and couples its genes through the cells' shared capture. In
# vi-a the young group's (the artifact in the claim's direction: young - old > 0, a "decrease"), in vi-b the
# old group's.
WIDE = 3000.0 * math.exp((0.3 ** 2 - 0.9 ** 2) / 2)
CAPTURE_YOUNG = dict(dom=(0.1, 0.1), capture_sd=(0.9, 0.3), depth=(WIDE, 3000.0))
CAPTURE_OLD = dict(dom=(0.1, 0.1), capture_sd=(0.3, 0.9), depth=(3000.0, WIDE))
SCENARIOS = {
    "i: equal coupling, gap as on §4": dict(rho=(0.5, 0.5), **GAP),
    "eq0: equal detection, equal coupling": dict(rho=(0.5, 0.5), **EQUAL),
    **{f"eq{k}: equal detection, decrease of {k} x SESOI": dict(rho=r, **EQUAL) for k, r in COUPLINGS.items()},
    **{f"gap{k}: gap as on §4, decrease of {k} x SESOI": dict(rho=r, **GAP) for k, r in COUPLINGS.items()},
    "vi-a: equal coupling, the young group's capture varying more": dict(rho=(0.5, 0.5), **CAPTURE_YOUNG),
    "vi-b: equal coupling, the old group's capture varying more": dict(rho=(0.5, 0.5), **CAPTURE_OLD),
}
EQUAL_DETECTION = [s for s in SCENARIOS if s.startswith("eq")]
CAPTURE_SCENARIOS = [s for s in SCENARIOS if s.startswith("vi")]
BASE_SCENARIOS = ["i: equal coupling, gap as on §4", *EQUAL_DETECTION, *CAPTURE_SCENARIOS]

CRITERIA = """The criteria (fixed before the run; the owner's second decisions of 2026-10-19, item 3), for each
of mi_3bin and norm_pearson; s_m is GATE 4's response to the SESOI (detail["sesoi_metric"]), the effect the
verdict's (groups[0] - groups[1], young - old: a decrease is positive):
(i)   an equal coupling with a gap of detection as on §4 (scenario i), the new engine: at most 2 of 20
      SUPPORTED (all false), and |mean of the effects| <= 0.5 x mean of s_m over the 20 seeds.
(ii)  an equal detection (scenarios eq0, eq3, eq5 and eq7: an equal coupling and decreases of 3, 5 and 7 x
      SESOI): the verdict (its class and its cause) of the new engine equals that of v0.4 without change 17
      (the engine of 1337ead) in every seed: 80 of 80.
(iii) a decrease of 3 x SESOI at an equal detection (eq3): the share of SUPPORTED of the new engine is at
      least that of the engine of 1337ead less 1/20.
(iv)  decreases of 3, 5 and 7 x SESOI with a gap as on §4 (gap3, gap5, gap7), the new engine: the mean over
      the seeds of effect / s_m is at least 0.7 x the same mean at an equal detection with the same
      coupling (eq3, eq5, eq7, the new engine; 1337ead's mean is reported beside it); the shares of
      SUPPORTED are reported. Seeds where s_m is not positive (GATE 4 failed) are left out of the means
      and counted.
(v)   the clean data of tests/test_v03_verdict.py (test_gates.make_clean with add_mice, norm_pearson of
      Smad3 and Col1a1 with the coupling probe), the new engine: SUPPORTED.
(vi)  a gap made by cell-to-cell variation of capture (scenarios vi-a and vi-b, both engines): reported only,
      as a known limitation.
The acceptance passes where (i)-(v) pass for both metrics. If it fails: the report, and stop (no further
design)."""


# --------------------------------------------------------------------------- #
# the worker: one engine, its own sys.path
# --------------------------------------------------------------------------- #
def _verdict_class(v: str) -> str:
    for sep in (" —", " (", " ["):
        v = v.split(sep)[0]
    return v.strip()


def _one(task):
    import numpy as np
    import metric_autopsy as m
    from metric_autopsy import injected_signal, metrics, run_autopsy
    from functools import partial
    t0 = time.time()
    if task["scenario"] == "v":
        from test_gates import add_mice, make_clean
        d = add_mice(make_clean())
        fn = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
        a = run_autopsy(fn, d, group_col="age", groups=("young", "old"), within=["sex"], replicate_col="mouse",
                        gene_pair=("Smad3", "Col1a1"), signal_test=injected_signal.coupling("Smad3", "Col1a1"),
                        prereg={"estimand": "composition", "direction": "decrease", "judgment_pending": False},
                        log_path="off", specificity=False)
    else:
        from detection_data import METRICS, PREREG, deeper_less_complex
        fn, mu = METRICS[task["metric"]]
        d = deeper_less_complex(mu_claim=mu, seed=task["seed"], **SCENARIOS[task["scenario"]])
        a = run_autopsy(fn, d, group_col="age", groups=("young", "old"), replicate_col="mouse", gene_pair=("A", "B"),
                        signal_test=injected_signal.coupling("A", "B"), prereg=dict(PREREG), log_path="off",
                        specificity=False)
    e = a.effect.detail if a.effect is not None else {}
    g4 = next((r for r in a.results if r.gate == 4), None)
    p = e.get("p_perm") if e.get("tier") == "permutation" else e.get("p_t")
    ed = e.get("equal_depth_effect") or {}

    def fl(x):
        return float(x) if x is not None and np.isfinite(x) else None
    return dict(task, verdict=a.verdict, cls=_verdict_class(a.verdict), cause=a.cause, effect=fl(e.get("effect")),
                raw_effect=fl(e.get("raw_effect")), p=fl(p), s_m=fl(g4.detail.get("sesoi_metric")) if g4 else None,
                g4=g4.status.value if g4 else None, g4_data=g4.detail.get("data") if g4 else None,
                equal_depth_effect=fl(ed.get("effect")),
                detection=[dict(gene=r["gene"], raw=r.get("raw"), after=r["after"], keep_p=r["keep_p"])
                           for r in (e.get("detection") or {}).get("records", [])],
                engine=m.__file__, seconds=round(time.time() - t0, 1))


def worker(src: str, tasks: list, workers: int) -> list:
    sys.path[:0] = [src, str(REPO / "tests")]
    import metric_autopsy
    assert str(Path(metric_autopsy.__file__).resolve()).startswith(str(Path(src).resolve())), metric_autopsy.__file__
    if workers > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(workers) as pool:
            return pool.map(_one, tasks, chunksize=1)
    return [_one(t) for t in tasks]


# --------------------------------------------------------------------------- #
# the scores
# --------------------------------------------------------------------------- #
def _mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else float("nan")


def score(new: list, base: list) -> tuple[list, dict]:
    lines, passed = [], {}
    by = {}
    for r in new:
        by.setdefault(("new", r["scenario"], r["metric"]), []).append(r)
    for r in base:
        by.setdefault(("base", r["scenario"], r["metric"]), []).append(r)

    def rows(engine, scen, metric):
        return sorted(by.get((engine, scen, metric), []), key=lambda r: r["seed"])

    def sup(rs):
        return sum(r["cls"] == "SUPPORTED" for r in rs)
    for metric in METRIC_NAMES:
        lines.append(f"### {metric}")
        # (i)
        rs = rows("new", "i: equal coupling, gap as on §4", metric)
        n_sup, me, ms = sup(rs), _mean([r["effect"] for r in rs]), _mean([r["s_m"] for r in rs])
        ok = n_sup <= 2 and abs(me) <= 0.5 * ms
        passed[(metric, "i")] = ok
        lines.append(f"- (i) {'PASS' if ok else 'FAIL'}: {n_sup} of {len(rs)} SUPPORTED (at most 2); mean effect "
                     f"{me:+.4g}, mean s_m {ms:.4g}, |mean effect| / mean s_m = {abs(me) / ms:.3f} (at most 0.5)")
        # (ii)
        agree, total, diffs = 0, 0, []
        for scen in EQUAL_DETECTION:
            nb = {r["seed"]: r for r in rows("base", scen, metric)}
            for r in rows("new", scen, metric):
                b = nb.get(r["seed"])
                total += 1
                if b is not None and (b["cls"], b["cause"]) == (r["cls"], r["cause"]):
                    agree += 1
                else:
                    diffs.append(f"{scen.split(':')[0]} seed {r['seed']}: {b and (b['cls'], b['cause'])} -> "
                                 f"{(r['cls'], r['cause'])}")
        ok = agree == total == 20 * len(EQUAL_DETECTION)
        passed[(metric, "ii")] = ok
        lines.append(f"- (ii) {'PASS' if ok else 'FAIL'}: {agree} of {total} verdicts (class and cause) as the "
                     "engine of 1337ead's" + (f"; differ: {'; '.join(diffs)}" if diffs else ""))
        # (iii)
        rn = rows("new", "eq3: equal detection, decrease of 3 x SESOI", metric)
        rb = rows("base", "eq3: equal detection, decrease of 3 x SESOI", metric)
        ok = sup(rn) / 20 >= sup(rb) / 20 - 1 / 20
        passed[(metric, "iii")] = ok
        lines.append(f"- (iii) {'PASS' if ok else 'FAIL'}: SUPPORTED {sup(rn)} of {len(rn)} (the new engine) against "
                     f"{sup(rb)} of {len(rb)} (1337ead's), at least {sup(rb) - 1}")
        # (iv)
        oks = []
        for k in COUPLINGS:
            gap = rows("new", f"gap{k}: gap as on §4, decrease of {k} x SESOI", metric)
            eq = rows("new", f"eq{k}: equal detection, decrease of {k} x SESOI", metric)
            eqb = rows("base", f"eq{k}: equal detection, decrease of {k} x SESOI", metric)
            rg, re_, rb_ = ([r["effect"] / r["s_m"] for r in rs if r["s_m"] and r["s_m"] > 0 and r["effect"] is not None]
                            for rs in (gap, eq, eqb))
            mg, mq = _mean(rg), _mean(re_)
            ok = mg >= 0.7 * mq
            oks.append(ok)
            lines.append(f"- (iv) {k} x SESOI {'PASS' if ok else 'FAIL'}: mean effect / s_m {mg:.3f} with the gap "
                         f"({len(rg)} seeds with s_m > 0) against {mq:.3f} at an equal detection ({len(re_)}), ratio "
                         f"{mg / mq:.3f} (at least 0.7); SUPPORTED {sup(gap)} of {len(gap)} with the gap, {sup(eq)} of "
                         f"{len(eq)} at an equal detection; 1337ead's at an equal detection: mean effect / s_m "
                         f"{_mean(rb_):.3f} ({len(rb_)}), SUPPORTED {sup(eqb)} of {len(eqb)}")
        passed[(metric, "iv")] = all(oks)
        # (vi)
        for scen in CAPTURE_SCENARIOS:
            for engine, label in (("new", "the new engine"), ("base", "1337ead's")):
                rs = rows(engine, scen, metric)
                if rs:
                    ratios = [r["effect"] / r["s_m"] for r in rs if r["s_m"] and r["s_m"] > 0 and r["effect"] is not None]
                    lines.append(f"- ({scen.split(':')[0]}) {label} (reported only): SUPPORTED {sup(rs)} of {len(rs)}; "
                                 f"mean effect {_mean([r['effect'] for r in rs]):+.4g}, mean effect / s_m "
                                 f"{_mean(ratios):.3f}")
        lines.append("")
    v = [r for r in new if r["scenario"] == "v"]
    ok = bool(v) and v[0]["cls"] == "SUPPORTED"
    passed[("norm_pearson", "v")] = ok
    vb = [r for r in base if r["scenario"] == "v"]
    lines.append(f"- (v) {'PASS' if ok else 'FAIL'}: the clean data of tests/test_v03_verdict.py: the new engine "
                 f"{v[0]['verdict'][:110] if v else '—'} [{v[0]['cause'] if v else '—'}]; 1337ead's "
                 f"{vb[0]['verdict'][:60] if vb else '—'}")
    return lines, passed


def table(rs: list, engine: str) -> list:
    out = [f"| scenario | metric | {engine}: SUPPORTED | verdicts | mean effect | mean s_m | mean effect at equal depth "
           "alone | detection after (A young/old) |", "|---|---|---|---|---|---|---|---|"]
    groups = {}
    for r in rs:
        groups.setdefault((r["scenario"], r["metric"]), []).append(r)
    for (scen, metric), g in groups.items():
        counts = {}
        for r in g:
            counts[f"{r['cls']} [{r['cause']}]"] = counts.get(f"{r['cls']} [{r['cause']}]", 0) + 1
        det = [x for r in g for x in r["detection"] if x["gene"] in ("A", "Smad3")]
        da = (f"{_mean([x['after']['young'] for x in det]):.3f}/{_mean([x['after']['old'] for x in det]):.3f}"
              if det else "—")
        out.append(f"| {scen} | {metric} | {sum(r['cls'] == 'SUPPORTED' for r in g)} of {len(g)} | "
                   + "; ".join(f"{k} {v}" for k, v in sorted(counts.items())) + " | "
                   f"{_mean([r['effect'] for r in g]):+.4g} | {_mean([r['s_m'] for r in g]):.4g} | "
                   f"{_mean([r['equal_depth_effect'] for r in g]):+.4g} | {da} |")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-src", required=False, help="the src of a worktree at 1337ead")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seeds", type=int, default=len(SEEDS), help=argparse.SUPPRESS)  # a smoke test of the machinery
    ap.add_argument("--worker", nargs=3, metavar=("SRC", "TASKS", "OUT"), help=argparse.SUPPRESS)
    a = ap.parse_args()
    if a.worker:
        src, tasks, out = a.worker
        res = worker(src, json.loads(Path(tasks).read_text()), a.workers)
        Path(out).write_text(json.dumps(res))
        return 0
    rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "--", "src"], cwd=REPO, capture_output=True, text=True).stdout.strip()
    base_rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=Path(a.base_src).parent, capture_output=True,
                              text=True).stdout.strip()
    import numpy as np
    print(f"# The acceptance of change 17' on synthetic data — the new engine at git {rev}"
          f"{' (with uncommitted changes to src)' if dirty else ''}, the engine without change 17 at {base_rev}; "
          f"python {platform.python_version()}, numpy {np.__version__}; {a.workers} workers\n")
    print(CRITERIA)
    print()
    sys.stdout.flush()
    tasks_new = [dict(scenario=s, metric=m, seed=k) for s in SCENARIOS for m in METRIC_NAMES for k in SEEDS[:a.seeds]]
    tasks_new.append(dict(scenario="v", metric="norm_pearson", seed=0))
    tasks_base = [t for t in tasks_new if t["scenario"] in BASE_SCENARIOS or t["scenario"] == "v"]
    tmp = Path(os.environ.get("TMPDIR", "/tmp")) / f"detection_equalization_{os.getpid()}"
    tmp.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONHASHSEED="0")
    results = {}
    t0 = time.time()
    for name, src, tasks in (("new", str(REPO / "src"), tasks_new), ("base", a.base_src, tasks_base)):
        (tmp / f"{name}.tasks.json").write_text(json.dumps(tasks))
        subprocess.run([sys.executable, "-I", str(Path(__file__).resolve()), "--workers", str(a.workers), "--worker", src,
                        str(tmp / f"{name}.tasks.json"), str(tmp / f"{name}.json")], check=True, env=env)
        results[name] = json.loads((tmp / f"{name}.json").read_text())
        print(f"{name}: {len(results[name])} runs, engine {results[name][0]['engine']}, {time.time() - t0:.0f} s")
    print()
    print("## The runs")
    print()
    for line in table(results["new"], "the new engine"):
        print(line)
    print()
    for line in table(results["base"], "1337ead's"):
        print(line)
    print()
    same = sum(1 for r in results["new"] for b in results["base"]
               if (r["scenario"], r["metric"], r["seed"]) == (b["scenario"], b["metric"], b["seed"])
               and r["equal_depth_effect"] is not None and b["effect"] is not None
               and abs(r["equal_depth_effect"] - b["effect"]) <= 1e-12 * max(1.0, abs(b["effect"])))
    both = sum(1 for r in results["new"] if r["equal_depth_effect"] is not None and r["scenario"] in BASE_SCENARIOS)
    print(f"The new engine's effect at equal depth alone equals 1337ead's effect in {same} of {both} runs on the same "
          "data (the thinning to equal depth draws as before).")
    print()
    print("## The criteria")
    print()
    lines, passed = score(results["new"], results["base"])
    for line in lines:
        print(line)
    print()
    ok = all(passed.values())
    print(f"The acceptance {'PASSES' if ok else 'FAILS'}: " + ", ".join(f"{m} ({c}) {'pass' if v else 'FAIL'}"
                                                                        for (m, c), v in passed.items()))
    print(f"\n{time.time() - t0:.0f} s")
    (HERE / "detection_equalization.json").write_text(json.dumps(dict(new=results["new"], base=results["base"],
                                                                      passed={f"{m} {c}": v for (m, c), v in passed.items()}),
                                                                 indent=0) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
