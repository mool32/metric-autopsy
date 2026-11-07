"""JOURNAL.md D8, the owner's item 7 of 2026-10-18: what p14's low-level expectation encodes.

p14 (``p14_gate4_by_expression_level.log``, git 446cd57) measured, per pool pair, the population response
of log-normalized Pearson to an injected coupling at the level's saturation dose (24 datasets x 10
injections), and called a pair blind below 0.8 x delta_min = 0.04 (v0.3's rule, delta_min = 0.5 x SESOI).
Its low-level expectation, GATE 4 FAILs with probability >= 0.8, rests on that threshold. v0.4 change 3
(decided 2026-10-17 and 2026-10-18) drops delta_min: GATE 4 FAILs only where the response's upper bound is
at or below 0. This script runs the engine installed, v0.4's GATE 4, on p14's 40 low-level null datasets
(the same datasets, seeds and dose as p14's step 3) and sets each outcome beside its pair's population
response from p14's log, classed by v0.4's own definition with the log's standard error: **responds**
where r - 1.96 SE > 0, **blind** where r + 1.96 SE <= 0, **undecided** otherwise. A PASS on a blind pair
certifies a response the population does not have (an error by the probe's truth); a FAIL on a pair that
responds calls a responsive metric blind.

    python d8_p14_low_level.py > d8_p14_low_level.log   # about 5 min on 2 cores
"""
import platform
import re
import subprocess
import sys
from functools import partial
from multiprocessing import Pool
from pathlib import Path

from probe_common import *  # noqa: F401,F403  (repo path, numpy, pandas, metrics, gates)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "prereg"))
import panel as PANEL  # noqa: E402
import simulate  # noqa: E402
import metric_autopsy  # noqa: E402
from metric_autopsy import SimpleData, gate4_signal_response, injected_signal, metrics  # noqa: E402

LEVEL, DOSE, GATE_DATASETS, BASE = "low", 2.25, 40, 9000  # p14's step 3 at the low level
P14_LOG = HERE / "p14_gate4_by_expression_level.log"
_BG = {}


def _setup():
    if "bgs" not in _BG:
        _BG["bgs"] = simulate.dry_backgrounds()
        _BG["pilot"] = simulate.dry_pilot(_BG["bgs"])
    return _BG["bgs"], _BG["pilot"]


def population_responses() -> dict:
    """{pool index: (response, SE)} of the low level, from p14's log (its section 2)."""
    out = {}
    pat = re.compile(r"pool\s+(\d+) low\s+: ([+-]\d\.\d+) \(SE (\d\.\d+)")
    for line in P14_LOG.read_text().splitlines():
        m = pat.search(line)
        if m:
            out[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
    return out


def gate(job):
    pair, i = job
    bgs, pilot = _setup()
    e = dict(id=f"G{i}", condition="N1", variant="null", index=i, side="A", pair=pair, seed=BASE + i)
    X, obs, genes, cards = PANEL.build(e, bgs, pilot)
    ga, gb = cards[0]["gene_pair"]
    r = gate4_signal_response(partial(metrics.norm_pearson, gene_a=ga, gene_b=gb), SimpleData(X, obs, genes),
                              injected_signal.coupling(ga, gb, strength=DOSE), delta_min=0.05, seed=i)
    d = r.detail
    return pair, i, d.get("outcome", r.status.value), d.get("mean_response"), d.get("signed_lower"), d.get("signed_upper")


def main():
    rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "--", str(REPO / "src")], cwd=HERE, capture_output=True,
                           text=True).stdout.strip()
    rev += " (with uncommitted changes to src)" if dirty else ""
    print(f"# JOURNAL D8: p14's low level under the engine installed — git {rev}, metric-autopsy "
          f"{metric_autopsy.__version__}, python {platform.python_version()}, numpy {np.__version__}")
    truth = population_responses()
    cls = {k: ("responds" if r - 1.96 * se > 0 else "blind" if r + 1.96 * se <= 0 else "undecided")
           for k, (r, se) in truth.items()}
    bgs, _ = _setup()
    low = [pe["index"] for pe in bgs["B1"].plan["pool"] if pe["level"] == LEVEL]
    if sorted(low) != sorted(truth):
        raise SystemExit(f"the pool's low pairs {low} are not the log's {sorted(truth)}")
    print(f"# p14's population responses (its log, section 2) classed by v0.4's definition (r -/+ 1.96 SE against 0):")
    for k in low:
        r, se = truth[k]
        print(f"#   pool {k}: {r:+.4f} (SE {se:.4f}, 95% [{r - 1.96 * se:+.4f}, {r + 1.96 * se:+.4f}]) -> {cls[k]}")
    jobs = [(low[i % len(low)], i) for i in range(GATE_DATASETS)]
    with Pool(2) as workers:
        out = workers.map(gate, jobs)
    print(f"1. GATE 4 (v0.4) on p14's {GATE_DATASETS} low-level null datasets at dose {DOSE:g}, per dataset")
    for pair, i, oc, m, lo, hi in out:
        print(f"  dataset {i:2d} pool {pair}: {oc:8s} response {m:+.4f} [{lo:+.4f}, {hi:+.4f}]  ({cls[pair]})")
    print("2. by the pair's class (v0.4's definition of blind)")
    for c in ("responds", "undecided", "blind"):
        rows = [r for r in out if cls[r[0]] == c]
        counts = {o: sum(r[2] == o for r in rows) for o in ("PASS", "FAIL", "UNTESTED")}
        print(f"  {c:9s}: {len(rows)} datasets — " + ", ".join(f"{o} {n}" for o, n in counts.items()))
    counts = {o: sum(r[2] == o for r in out) for o in ("PASS", "FAIL", "UNTESTED")}
    print(f"3. all {len(out)}: " + ", ".join(f"{o} {n}" for o, n in counts.items()) + "; p14's expectation (v0.3's "
          f"threshold 0.8 x delta_min): FAIL in at least 80%")


if __name__ == "__main__":
    main()
