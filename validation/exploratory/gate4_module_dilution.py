"""GATE 4's module probe on a mean-expression metric: how its response scales with the probe's
`frac` and with the share of cells that carry the marker. Exploratory, after R1's result, on
synthetic counts (and, with --compact, on B2's own Xist cells); the engine is not changed.

The metric is anchor R1's: the mean over cells of log1p CP10k of one gene (anchors.log_cp10k_mean).
GATE 4 is the engine's own (gates.gate4_signal_response, 200 injections against the sham, seed 0)
with delta_min 0.25, which is 0.5 x SESOI 0.5, R1's values. The probe is
injected_signal.module([marker], fold 2, frac): frac 0.3 is the module's default, frac 1.0 is
R1's. The module raises the marker 2-fold in a random `frac` of cells, against a sham that thins
every gene in another random `frac` of cells, so in expectation the response over all cells is
frac x (share of cells with the marker) x (the response in those cells).

* The synthetic grid: depth per cell about 20,000 counts (B2's cells: a median of 23,336 over the
  2,000 genes kept) and 2,000; the marker's CP10k in the cells that carry it about 1 (low), 6 (as
  B2's Xist: a median of 5.9) or 50 (high); the share of cells with the marker 25%, 50% and 100%;
  frac 0.3 and 1.0. A cell with the marker has at least one count, so the share is exact.
* The smallest share at which R1's probe (frac 1.0) does not fail: delta_min over the response at
  a share of 100% (the response scales with the share), checked by runs just below and above it.
* With --compact (compact/B2.npz of results/panel-v1), the same on B2's cells: 1,000 cells drawn
  with a given share of cells with Xist > 0.

    python validation/exploratory/gate4_module_dilution.py [--compact compact/B2.npz]
"""
from __future__ import annotations

import argparse
import hashlib
from functools import partial
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation" / "prereg"))
import anchors  # noqa: E402  (sets the numerical environment before numpy loads)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from metric_autopsy import SimpleData, __version__, gates, injected_signal  # noqa: E402

DELTA_MIN = 0.5 * anchors.SESOI  # 0.25, R1's
SEED = 0  # run_autopsy's GATE 4 seed
DATA_SEED = 20261010
N_CELLS, N_GENES = 1000, 300
DEPTHS = (20000, 2000)
LEVELS = (("low", 1.0), ("as B2's Xist", 6.0), ("high", 50.0))  # the marker's CP10k where present
SHARES = (0.25, 0.50, 1.00)
FRACS = (0.3, 1.0)


def synthetic(depth: float, cp10k: float, share: float, seed: int) -> SimpleData:
    """Counts of N_CELLS cells: N_GENES background genes with a skewed profile at a per-cell depth
    around `depth`, and a marker in a `share` of cells at about `cp10k` (at least one count)."""
    rng = np.random.default_rng(seed)
    profile = rng.gamma(0.6, 1.0, N_GENES)
    profile /= profile.sum()
    d = depth * np.exp(rng.normal(0.0, 0.3, N_CELLS))
    X = rng.poisson(np.outer(d, profile)).astype(float)
    marker = np.zeros(N_CELLS)
    pos = rng.permutation(N_CELLS)[: int(round(share * N_CELLS))]
    lam = d[pos] * cp10k / 1e4
    marker[pos] = 1 + rng.poisson(np.maximum(lam - 1.0, 0.0))
    genes = ["Marker", *[f"G{j}" for j in range(N_GENES)]]
    return SimpleData(np.column_stack([marker, X]), pd.DataFrame(index=range(N_CELLS)), genes)


def gate4(data: SimpleData, gene: str, frac: float) -> dict:
    cols = anchors._cols(data.var_names, [gene])
    metric = partial(anchors.log_cp10k_mean, cols=cols)
    g = gates.gate4_signal_response(metric, data, injected_signal.module([gene], fold=2.0, frac=frac),
                                    direction="increase", delta_min=DELTA_MIN, alpha=0.05, seed=SEED)
    d = g.detail
    return dict(mean=d["mean_response"], lo=d["ci"][0], hi=d["ci"][1], status=g.status.value)


def fmt(r: dict) -> str:
    return f"{r['mean']:+.4f} ({r['lo']:+.4f} to {r['hi']:+.4f}) {r['status']}"


def threshold_lines(label: str, full: float, check) -> list[str]:
    """The smallest share at frac 1.0 that does not fail, from the response at 100%, and runs at
    the share 5 points below and above it (`check(share) -> result`)."""
    if full <= DELTA_MIN:
        return [f"  {label}: the response at a share of 100% is {full:+.4f}, below delta_min {DELTA_MIN:g}: "
                "R1's probe fails at every share"]
    star = DELTA_MIN / full
    out = [f"  {label}: {DELTA_MIN:g} / {full:.4f} = {100 * star:.1f}% of cells"]
    for s in (star - 0.05, star + 0.05):
        if 0.0 < s <= 1.0:
            out.append(f"    at {100 * s:.1f}%: {fmt(check(s))}")
    return out


def compact_xist(path: Path, n: int = 1000) -> list[str]:
    z = np.load(path, allow_pickle=False)
    genes = [str(g) for g in z["genes"]]
    X = z["X"]
    j = genes.index("Xist")
    pos, neg = np.where(X[:, j] > 0)[0], np.where(X[:, j] == 0)[0]

    def draw(share: float) -> SimpleData:
        rng = np.random.default_rng(DATA_SEED)
        k = int(round(share * n))
        rows = np.sort(np.r_[rng.choice(pos, size=k, replace=False), rng.choice(neg, size=n - k, replace=False)])
        return SimpleData(X[rows], pd.DataFrame(index=range(n)), genes)

    full = {f: gate4(draw(1.0), "Xist", f) for f in FRACS}
    out = [f"## B2's Xist cells (compact B2, {path.name}, sha256 {hashlib.sha256(path.read_bytes()).hexdigest()})",
           f"{n} cells drawn (seed {DATA_SEED}) with a given share of cells with Xist > 0; the panel's 2,000 genes",
           f"  share 100%, frac 0.3: {fmt(full[0.3])}",
           f"  share 100%, frac 1.0: {fmt(full[1.0])}",
           "the smallest share of cells with Xist at which R1's probe (frac 1.0) does not fail:"]
    out += threshold_lines("Xist", full[1.0]["mean"], lambda s: gate4(draw(s), "Xist", 1.0))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--compact", type=Path, help="compact/B2.npz of results/panel-v1 (B2's own Xist cells)")
    args = ap.parse_args(argv)
    out = ["# GATE 4's module probe on a mean-expression metric: frac and the share of cells with the marker",
           "# Exploratory, after R1's result, on synthetic counts; the engine is not changed.",
           f"# metric_autopsy {__version__}; script sha256 {hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}; "
           f"injected_signal.module's defaults: fold 2.0, frac 0.3",
           f"# metric: mean log1p CP10k of the marker (anchors.log_cp10k_mean); GATE 4: {gates.GATE4_N_REP} injections, "
           f"seed {SEED}, delta_min {DELTA_MIN:g} (0.5 x SESOI {anchors.SESOI:g}); data seed {DATA_SEED}; "
           f"{N_CELLS} cells x ({N_GENES} genes + the marker)", ""]
    full = {}
    for depth in DEPTHS:
        out += [f"## Depth about {depth:,} counts per cell", "",
                f"  {'marker CP10k':22s} {'share':>6s}   {'frac 0.3: response (95% CI) GATE 4':42s} "
                f"{'frac 1.0: response (95% CI) GATE 4':42s} ratio"]
        for name, cp in LEVELS:
            for share in SHARES:
                data = synthetic(depth, cp, share, DATA_SEED)
                r = {f: gate4(data, "Marker", f) for f in FRACS}
                full[(depth, name, share)] = r
                ratio = r[0.3]["mean"] / r[1.0]["mean"] if r[1.0]["mean"] else float("nan")
                out.append(f"  {f'{name} ({cp:g})':22s} {100 * share:5.0f}%   {fmt(r[0.3]):42s} {fmt(r[1.0]):42s} "
                           f"{ratio:.3f}")
        out.append("")
        out.append("  the response against the share (frac 1.0), as a fraction of the response at 100%:")
        for name, cp in LEVELS:
            at = {s: full[(depth, name, s)][1.0]["mean"] for s in SHARES}
            out.append(f"    {name}: " + ", ".join(f"{100 * s:.0f}% -> {at[s] / at[1.0]:.3f}" for s in SHARES))
        out += ["", "  the smallest share of cells with the marker at which R1's probe (frac 1.0) does not fail:"]
        for name, cp in LEVELS:
            out += threshold_lines(f"{name} ({cp:g})", full[(depth, name, 1.0)][1.0]["mean"],
                                   lambda s, depth=depth, cp=cp: gate4(synthetic(depth, cp, s, DATA_SEED), "Marker", 1.0))
        best = max(full[(depth, name, 1.0)][0.3]["hi"] for name, _ in LEVELS)
        out += ["", f"  frac 0.3 at a share of 100%: the largest upper bound is {best:+.4f}, against delta_min "
                f"{DELTA_MIN:g}", ""]
    if args.compact:
        out += compact_xist(args.compact)
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
