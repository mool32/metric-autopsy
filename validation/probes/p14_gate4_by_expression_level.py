"""GATE 4 on the valid metric by the expression level of the pair (found on v0.3.0.dev0).

The confirmatory panel (validation/prereg/v1.md) draws its analysed pairs from three expression
levels. On simulated backgrounds of the panel's planned size (validation/prereg/simulate.py:
negative-binomial counts with cell-size variation, coupled pairs at every level), GATE 4 injects
a coupling of the pair (``injected_signal.coupling``) into null datasets of 2 x 8 donors x 200
cells. The first version of this probe (git 317080e) used the former rule, z >= 3 over 10
injections, at strengths 1.0 and 2.0: at the medium level a real but weak response failed in
8/12 and 2/12 datasets, at the low level the metric does not respond at all and failed 12/12
(log-normalized Pearson is computed on co-detected cells, where low counts are mostly 1 and
their log-normalized values follow the cell's total, so a coupling that acts through
co-detection is invisible to it).

Decided by the project owner on 2026-10-08 (JOURNAL.md, D5): GATE 4 PASSes when the lower 95%
bound of the response is above 0, FAILs only when the upper bound is below delta_min (default
0.5 x SESOI), and is UNTESTED otherwise; the injection is the largest dose of the grid without
saturation. This probe measures the rule as a regression:

1. per level, the response of the metric to the injection (injected minus sham, on whole
   datasets, as GATE 4 measures it) over the dose grid 0.25-4.0 of the pilot (a population
   simulation), and the level's saturation dose: the smallest grid dose whose response reaches
   95% of the grid maximum (the largest grid dose where the maximum is not above 0);
2. per pool pair, its population response at that dose against delta_min = 0.05 (0.5 x the
   stand-in SESOI 0.1 of ``simulate.dry_pilot``): valid above 1.2 x delta_min, blind below
   0.8 x delta_min, ambiguous between;
3. GATE 4 (the engine, its default number of injections) on fresh null datasets at the level's
   saturation dose: PASS / FAIL / UNTESTED.

Required (the owner's regression): the high level PASSes; at the medium level FAIL is no more
frequent than alpha = 0.05 on the valid pairs; at the low level FAIL has probability >= 0.8.
``test_probes.py`` holds a smaller copy of step 3 at the doses of this log.

    python p14_gate4_by_expression_level.py > p14_gate4_by_expression_level.log   # ~15 min, 4 cores
"""
import platform
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
from metric_autopsy import SimpleData, gate4_signal_response, injected_signal, metrics  # noqa: E402
from metric_autopsy import gates as G  # noqa: E402

DOSE_GRID = tuple(round(0.25 * k, 2) for k in range(1, 17))  # the pilot's grid (oracle.DOSE_GRID)
SESOI = 0.1
DELTA_MIN = G.DELTA_MIN_FRACTION * SESOI
CURVE_DATASETS, CURVE_REPS = 6, 4      # per pair and dose, for the saturation dose
TRUTH_DATASETS, TRUTH_REPS = 24, 10    # per pair, at the saturation dose
GATE_DATASETS = 40                     # per level, GATE 4 itself
_BG = {}


def _setup():
    if "bgs" not in _BG:
        _BG["bgs"] = simulate.dry_backgrounds()
        _BG["pilot"] = simulate.dry_pilot(_BG["bgs"])
    return _BG["bgs"], _BG["pilot"]


def dataset(pair, i, base):
    bgs, pilot = _setup()
    e = dict(id=f"G{i}", condition="N1", variant="null", index=i, side="A", pair=pair, seed=base + i)
    X, obs, genes, cards = PANEL.build(e, bgs, pilot)
    return SimpleData(X, obs, genes), cards[0]["gene_pair"]


def response(job):
    """Mean of injected-minus-sham log-normalized Pearson on one dataset (`reps` injections)."""
    pair, i, base, dose, reps = job
    sd, (ga, gb) = dataset(pair, i, base)
    inj = injected_signal.coupling(ga, gb, strength=dose)
    rng = np.random.default_rng([base, i, int(dose * 100)])
    m = partial(metrics.norm_pearson, gene_a=ga, gene_b=gb)
    return job, float(np.mean([m(inj(sd, rng)) - m(inj.sham(sd, rng)) for _ in range(reps)]))


def gate(job):
    pair, i, dose = job
    sd, (ga, gb) = dataset(pair, i, 9000)
    r = gate4_signal_response(partial(metrics.norm_pearson, gene_a=ga, gene_b=gb), sd,
                              injected_signal.coupling(ga, gb, strength=dose), delta_min=DELTA_MIN, seed=i)
    return job, r.detail.get("outcome", r.status.value), r.detail.get("mean_response"), r.detail.get("signed_upper")


def saturation_dose(curve: dict) -> float:
    top = max(curve.values())
    if top <= 0:
        return DOSE_GRID[-1]
    return next(d for d in DOSE_GRID if curve[d] >= 0.95 * top)


def main():
    rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, capture_output=True, text=True).stdout.strip()
    print(f"# metric-autopsy probe p14 — git {rev}, python {platform.python_version()}, numpy {np.__version__}")
    bgs, _ = _setup()
    pool = bgs["B1"].plan["pool"]
    for pe in pool:
        print(f"# pool {pe['index']:2d} {pe['level']:6s} {pe['pair']} mean {pe['mean'][0]:.2f}/{pe['mean'][1]:.2f} "
              f"detected {pe['detection'][0]:.2f}/{pe['detection'][1]:.2f}")
    print(f"# delta_min {DELTA_MIN:g} (0.5 x SESOI {SESOI:g}); GATE 4 with {G.GATE4_N_REP} injections, 95% CI")
    by_level = {lv: [pe["index"] for pe in pool if pe["level"] == lv] for lv in PANEL.LEVELS}
    level_of = {pe["index"]: pe["level"] for pe in pool}
    with Pool(4) as workers:
        jobs = [(k, i, 1000, d, CURVE_REPS) for lv in PANEL.LEVELS for k in by_level[lv]
                for d in DOSE_GRID for i in range(CURVE_DATASETS)]
        curve = {}
        for (k, i, _, d, _), v in workers.map(response, jobs):
            curve.setdefault(level_of[k], {}).setdefault(d, []).append(v)
        dose = {}
        print("1. response by dose (level mean over its pairs and datasets) and the saturation dose")
        for lv in PANEL.LEVELS:
            c = {d: float(np.mean(v)) for d, v in curve[lv].items()}
            dose[lv] = saturation_dose(c)
            print(f"  {lv:6s}: " + " ".join(f"{d:g}:{c[d]:+.3f}" for d in DOSE_GRID) + f"  -> dose {dose[lv]:g}")
        jobs = [(k, i, 5000, dose[level_of[k]], TRUTH_REPS) for lv in PANEL.LEVELS
                for k in by_level[lv] for i in range(TRUTH_DATASETS)]
        per_pair = {}
        for (k, *_), v in workers.map(response, jobs):
            per_pair.setdefault(k, []).append(v)
        truth = {}
        print("2. population response per pair at the level's dose, against delta_min")
        for k in sorted(per_pair):
            v = np.asarray(per_pair[k])
            r = float(v.mean())
            truth[k] = ("valid" if r >= 1.2 * DELTA_MIN else "blind" if r <= 0.8 * DELTA_MIN else "ambiguous")
            print(f"  pool {k:2d} {level_of[k]:6s}: {r:+.4f} (SE {v.std(ddof=1) / np.sqrt(len(v)):.4f}, "
                  f"between-dataset SD {v.std(ddof=1):.4f}) -> {truth[k]}")
        jobs = [(by_level[lv][i % len(by_level[lv])], i, dose[lv]) for lv in PANEL.LEVELS
                for i in range(GATE_DATASETS)]
        out = workers.map(gate, jobs)
    print(f"3. GATE 4 on {GATE_DATASETS} null datasets per level at its dose")
    for lv in PANEL.LEVELS:
        rows = [(job, oc, m, u) for job, oc, m, u in out if level_of[job[0]] == lv]
        counts = {o: sum(r[1] == o for r in rows) for o in ("PASS", "FAIL", "UNTESTED")}
        valid = [r for r in rows if truth[r[0][0]] == "valid"]
        print(f"  {lv:6s} (dose {dose[lv]:g}): " + ", ".join(f"{o} {n}/{len(rows)}" for o, n in counts.items())
              + f"; on the valid pairs FAIL {sum(r[1] == 'FAIL' for r in valid)}/{len(valid)}"
              + f"; response median {np.median([r[2] for r in rows]):+.4f}, upper bound median "
              f"{np.median([r[3] for r in rows]):+.4f}")


if __name__ == "__main__":
    main()
