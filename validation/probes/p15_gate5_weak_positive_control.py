"""GATE 5's silent-positive-control rule on weak but valid controls (found on v0.3.0.dev0).

The question of the project owner (2026-10-08): GATE 5 FAILs a metric whose positive control is
silent where a reference detector had power >= 0.8 to show an injected coupling of dose 2.0 (the
power rule of 2026-10-07). Does that rule FAIL a valid metric more often than alpha when the
control is coupled, but weakly? If so, the silent control is to be judged as GATE 4 is: FAIL only
when the metric's own response to a coupling injected into the decoupled control pair is shown
to be below delta_min (upper 95% bound < delta_min), else untested.

Data of the p14 type (validation/prereg/simulate.py: negative-binomial counts with cell-size
variation, coupled pairs at three expression levels), null datasets of 2 x 8 donors x 200 cells,
the negative control the level's own. Positive controls:

* planted: a planted coupled pair of the level (a shared per-cell factor, exp(0.5 z) on both);
* weak, dose d: an uncoupled high-level pair (a negative-control pair of the pool) with a
  coupling of dose d injected into the dataset (``injected_signal.coupling``): real, but weak.
  The injection thins both genes to about half, which keeps high-level genes high (at the medium
  level it would push them into the low level, where the metric is blind: a first version of
  this probe did that and is not kept).

Metrics: log-normalized Pearson, and a blind metric: the raw-count Pearson of gene a with a fixed
third gene, which never reads gene b. The truth of every case is the metric's own response to the
coupling (dose 2.0) injected into the decoupled control pair, measured on every dataset (the
response rule's estimate, averaged over the 20 datasets): valid above 1.2 x delta_min, blind below
0.8 x delta_min. delta_min = 0.05 (0.5 x the stand-in SESOI 0.1 of ``simulate.dry_pilot``),
alpha = 0.05. A FAIL is false where the metric is valid.

Three rules on every dataset: the power rule as the engine had it until D6 (a frozen copy below,
so that the comparison reproduces at any commit), the engine's GATE 5 at this commit (with
delta_min), and the probe's own copy of the response rule, which also measures the truth.

    python p15_gate5_weak_positive_control.py > p15_gate5_weak_positive_control.log   # ~15 min
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
from metric_autopsy import SimpleData, gate5_controls, injected_signal, metrics  # noqa: E402
from metric_autopsy import gates as G  # noqa: E402
from metric_autopsy.stats import extend_null  # noqa: E402

DATASETS = 20
DELTA_MIN = 0.05
N_REP = G.GATE4_N_REP
CASES = [("planted", "high", None), ("planted", "medium", None), ("planted", "low", None),
         ("weak", "high", 0.25), ("weak", "high", 0.35), ("weak", "high", 0.5)]
_BG = {}


# --------------------------------------------------------------------------- #
# the power rule of 2026-10-07, as removed from gates.py by D6 (a frozen copy)
# --------------------------------------------------------------------------- #
POWER_MIN, POWER_REPS = 0.8, 40


def _ranks(x):
    _, inv, counts = np.unique(x, return_inverse=True, return_counts=True)
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x))
    r[order] = np.arange(len(x))
    return (np.bincount(inv, weights=r) / counts)[inv]


def _ref_coupling(a, b):
    """The removed rule's reference detector: Spearman correlation of raw counts."""
    ra, rb = _ranks(np.asarray(a, float)), _ranks(np.asarray(b, float))
    ra, rb = ra - ra.mean(), rb - rb.mean()
    den = float(np.sqrt((ra * ra).sum() * (rb * rb).sum()))
    return float((ra * rb).sum() / den) if den > 0 else 0.0


def power_rule_power(sd, pair, dose, alpha_s, rng, n_rep=POWER_REPS, n_null=200, threshold=POWER_MIN):
    """The removed ``_positive_control_power``: the stratum's power to establish a coupling of `dose`
    injected into the decoupled control pair, with the reference detector against the
    depth-matched self-null at alpha/K (stopping once ``power >= threshold`` is decided)."""
    X = np.asarray(sd.X, float)
    names = [str(g) for g in sd.var_names]
    ia, ib = names.index(pair[0]), names.index(pair[1])
    a0, b0, rest = X[:, ia], X[:, ib], G._pair_depth(X, ia, ib)
    obs = pd.DataFrame(index=range(len(a0)))
    inject = injected_signal.coupling("a", "b", strength=dose)
    need = int(np.ceil(threshold * n_rep))
    hits = done = 0
    for _ in range(n_rep):
        base = SimpleData(np.column_stack([a0, G._depth_matched_draw(b0, rest, rng), rest]), obs, ["a", "b", "rest"])
        inj = np.asarray(inject(base, rng).X)
        a, b = inj[:, 0], inj[:, 1]

        def draw(k, a=a, b=b):
            return np.array([_ref_coupling(a, G._depth_matched_draw(b, rest, rng)) for _ in range(int(k))])
        _, pval, _ = extend_null(_ref_coupling(a, b), draw, n_null, alpha_s)
        hits += bool(np.isfinite(pval) and pval < alpha_s)
        done += 1
        if hits >= need or hits + (n_rep - done) < need:
            break
    return hits / done


def blind_to_b(data, *, gene_a, gene_b, third):
    """Pearson of gene a's raw counts with a fixed third gene's: never reads gene b."""
    return metrics.pearson(data, gene_a=gene_a, gene_b=third)


def response_rule(pair_metric, sd, pair, dose, rng):
    """The GATE 4 rule on the control pair: gene b decoupled (a depth-matched draw), a coupling of
    `dose` injected against its sham, the metric's own response (``gates.response_interval``,
    ``gates.judge_response``): FAIL if the upper 95% bound is below delta_min."""
    X = np.asarray(sd.X, float).copy()
    names = [str(g) for g in sd.var_names]
    ia, ib = names.index(pair[0]), names.index(pair[1])
    b0, rest = X[:, ib].copy(), G._pair_depth(X, ia, ib)
    inject = injected_signal.coupling(pair[0], pair[1], strength=dose)
    deltas = []
    for _ in range(N_REP):
        X[:, ib] = G._depth_matched_draw(b0, rest, rng)
        base = SimpleData(X, sd.obs, names)
        deltas.append(pair_metric(inject(base, rng), gene_a=pair[0], gene_b=pair[1])
                      - pair_metric(inject.sham(base, rng), gene_a=pair[0], gene_b=pair[1]))
    iv = G.response_interval(deltas)
    return G.judge_response(iv, DELTA_MIN), iv["mean_response"], iv["signed_upper"]


def one(job):
    kind, lvl, dose, mname, i = job
    if "bgs" not in _BG:
        _BG["bgs"] = simulate.dry_backgrounds()
        _BG["pilot"] = simulate.dry_pilot(_BG["bgs"])
    bgs, pilot = _BG["bgs"], _BG["pilot"]
    pool = bgs["B1"].plan["pool"]
    level_pool = [pe for pe in pool if pe["level"] == lvl]
    e = dict(id=f"G{i}", condition="N1", variant="null", index=i, side="A",
             pair=next(pe["index"] for pe in pool if pe["level"] == "low"), seed=2000 + i)
    X, obs, genes, _ = PANEL.build(e, bgs, pilot)
    pe = level_pool[i % len(level_pool)]
    neg = tuple(pe["neg_pair"])
    rng = np.random.default_rng(10_000 + i)
    sd = SimpleData(X, obs, genes)
    if kind == "planted":
        pos = tuple(pe["pair"])
    else:  # a real but weak coupling: the uncoupled negative-control pair of another pool entry
        other = level_pool[(i + 1) % len(level_pool)]
        pos = tuple(other["neg_pair"])
        sd = injected_signal.coupling(*pos, strength=dose)(sd, rng)
    pm = metrics.norm_pearson
    if mname == "blind":
        third = next(g for q in pool for g in q["pair"] if g not in pos + neg)
        pm = partial(blind_to_b, third=third)
    r = gate5_controls(pm, sd, pos, neg, seed=i, delta_min=DELTA_MIN)
    row = r.detail["rows"][0]
    engine = ("FAIL (silent, shown blind)" if row.get("pos_blind") else
              "FAIL (negative control)" if not row["neg_ok"] else
              "fires" if row["pos_fires"] else "silent, WARN")
    # the removed power rule, on the same controls (a firing control never fails)
    power = (None if row["pos_fires"] else
             power_rule_power(sd, pos, G.POSITIVE_CONTROL_DOSE, 0.05, np.random.default_rng(30_000 + i)))
    power_rule = ("FAIL (silent, power >= 0.8)" if power is not None and power >= POWER_MIN else
                  "FAIL (negative control)" if not row["neg_ok"] else
                  "fires" if row["pos_fires"] else "silent, WARN")
    # the response rule runs on every dataset: its estimate is the truth; its FAIL counts only where
    # the control is silent (a firing control never fails)
    resp = response_rule(pm, sd, pos, G.POSITIVE_CONTROL_DOSE, np.random.default_rng(20_000 + i))
    return job, power_rule, engine, row["neg_ok"], row["pos_fires"], resp


def main():
    rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, capture_output=True, text=True).stdout.strip()
    print(f"# metric-autopsy probe p15 — git {rev}, python {platform.python_version()}, numpy {np.__version__}")
    print(f"# {DATASETS} null datasets per row; delta_min {DELTA_MIN}; response rule {N_REP} injections of "
          f"dose {G.POSITIVE_CONTROL_DOSE:g} into the decoupled control pair")
    jobs = [(k, lvl, d, m, i) for k, lvl, d in CASES for m in ("norm_pearson", "blind")
            if not (m == "blind" and (k, lvl) != ("planted", "high")) for i in range(DATASETS)]
    with Pool(4) as pool:
        out = pool.map(one, jobs)
    rows = {}
    for job, power_rule, engine, neg_ok, fires, resp in out:
        rows.setdefault(job[:4], []).append(dict(power=power_rule, engine=engine, neg_ok=neg_ok, fires=fires,
                                                 resp=resp))
    print("# per row: the truth (the metric's mean response); the removed power rule; the engine's GATE 5 at "
          "this commit; the probe's copy of the response rule on the silent controls; false (or correct) FAILs "
          "of a silent control by each rule")
    for (kind, lvl, dose, m), rr in rows.items():
        what = f"{kind}{'' if dose is None else f' dose {dose:g}'}"
        r = float(np.mean([x["resp"][1] for x in rr]))
        t = "valid" if r >= 1.2 * DELTA_MIN else "blind" if r <= 0.8 * DELTA_MIN else "ambiguous"

        def tally(key):
            c = {}
            for x in rr:
                c[x[key]] = c.get(x[key], 0) + 1
            return ", ".join(f"{k} {v}/{len(rr)}" for k, v in sorted(c.items()))
        oc = {}
        for x in rr:
            if not x["fires"]:
                oc[x["resp"][0]] = oc.get(x["resp"][0], 0) + 1
        power_fail = sum(x["power"].startswith("FAIL (silent") for x in rr)
        engine_fail = sum(x["engine"].startswith("FAIL (silent") for x in rr)
        resp_fail = oc.get("FAIL", 0)
        what_fail = "false FAIL" if t == "valid" else "FAIL (correct where blind)" if t == "blind" else "FAIL"
        print(f"{m:12s} PC {what:13s} {lvl:6s}: response {r:+.4f} -> metric {t} | power rule: {tally('power')} | "
              f"engine GATE 5: {tally('engine')} | response rule on the silent ones: "
              + (", ".join(f"{k} {v}/{len(rr)}" for k, v in sorted(oc.items())) or "none silent")
              + f" | {what_fail} of a silent control: power rule {power_fail}/{len(rr)}, engine {engine_fail}/{len(rr)}, "
              f"response rule {resp_fail}/{len(rr)} | negative control outside its null "
              f"{sum(not x['neg_ok'] for x in rr)}/{len(rr)}")


if __name__ == "__main__":
    main()
