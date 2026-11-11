#!/usr/bin/env python3
"""v0.4 change 10 at a realistic size: a claim about a pair run end to end on a sparse matrix, its time and peak memory.

Synthetic counts (no real data): mice nested in age, a coupled pair (Smad3, Col1a1), a coupled control pair
(Actb, Gapdh) and fillers, most of low expression (about `--density` of the cells with a count), one in twenty
expressed as the pairs are (matched controls for the check of specificity). The claim:
`norm_pearson` of the pair, GATE 5's controls, the coupling probe dosed to a SESOI of 0.1 on the construct's
scale, the check of specificity (1,000 matched controls) and the technical covariates. The matrix is passed as a
scipy.sparse CSR matrix; the log gives the matrix, its stored entries, what a dense float64 copy of it would take,
the run's wall time and the process's peak resident memory (``resource.getrusage``, the whole process: Python, the
data and the engine). Exploratory: one run per size, nothing tuned.

    python validation/exploratory/sparse_scale.py --cells 10000 --genes 10000 > validation/exploratory/sparse_scale.log
"""
from __future__ import annotations

import argparse
import platform
import resource
import subprocess
import sys
import time
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[1] / "src"
sys.path.insert(0, str(SRC))  # this checkout's engine

import metric_autopsy as _engine  # noqa: E402
from metric_autopsy import SimpleData, injected_signal, metrics, run_autopsy  # noqa: E402

if Path(_engine.__file__).resolve().parent.parent != SRC.resolve():
    sys.exit(f"the engine imported is {_engine.__file__}, not this checkout's ({SRC})")


def peak_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0  # kB on Linux


def make(cells: int, genes: int, mice: int, density: float, seed: int):
    """CSR counts built one mouse at a time (never dense whole)."""
    import scipy.sparse as sp
    rng = np.random.default_rng(seed)
    n_fill = genes - 4
    # a filler's mean count per cell: Poisson zero share 1 - density on average, means varied; one filler in
    # twenty expressed as the pairs are, so that the check of specificity has matched controls
    mu = -np.log(1 - density) * np.exp(rng.normal(0, 0.7, n_fill)) / np.exp(0.7 ** 2 / 2)
    high = rng.random(n_fill) < 0.05
    mu[high] = np.exp(rng.normal(1.2, 0.5, int(high.sum())))
    per = cells // (2 * mice)
    blocks, rows = [], []
    for age in ("young", "old"):
        for m in range(mice):
            lat, hk = rng.normal(size=(per, 1)), rng.normal(size=(per, 1))
            pair = rng.poisson(np.exp(1.2 + 0.9 * lat + rng.normal(0, 0.3, (per, 2))))
            ctrl = rng.poisson(np.exp(1.6 + 0.9 * hk + rng.normal(0, 0.2, (per, 2))))
            size = np.exp(rng.normal(0, 0.3, per))
            fill = sp.csr_matrix(rng.poisson(mu[None, :] * size[:, None]).astype(np.float64))
            blocks.append(sp.hstack([sp.csr_matrix(np.hstack([pair, ctrl]).astype(np.float64)), fill], format="csr"))
            rows += [dict(age=age, mouse=f"{age}{m}")] * per
    X = sp.vstack(blocks, format="csr")
    names = ["Smad3", "Col1a1", "Actb", "Gapdh"] + [f"G{i}" for i in range(n_fill)]
    return SimpleData(X, pd.DataFrame(rows), names)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=int, default=10000)
    ap.add_argument("--genes", type=int, default=10000)
    ap.add_argument("--mice", type=int, default=5, help="mice per group")
    ap.add_argument("--density", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "--", str(SRC)], cwd=HERE, capture_output=True,
                           text=True).stdout.strip()
    print(f"# metric-autopsy {_engine.__version__} — change 10 at scale — git {rev}{' (src modified)' if dirty else ''}, "
          f"python {platform.python_version()}, numpy {np.__version__}")
    t0 = time.time()
    d = make(a.cells, a.genes, a.mice, a.density, a.seed)
    nnz = int(d.X.nnz)
    print(f"matrix: {d.n_obs} cells x {d.n_vars} genes, {a.mice} mice per group; stored entries {nnz} "
          f"({nnz / (d.n_obs * d.n_vars):.1%}); CSR {(nnz * 12 + (d.n_obs + 1) * 4) / 1e6:.0f} MB; "
          f"a dense float64 copy would take {d.n_obs * d.n_vars * 8 / 1e6:.0f} MB")
    print(f"built in {time.time() - t0:.1f} s; peak resident memory so far {peak_mb():.0f} MB")
    t1 = time.time()
    res = run_autopsy(partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1"), d, group_col="age",
                      groups=("young", "old"), gene_pair=("Smad3", "Col1a1"), replicate_col="mouse",
                      pair_metric=metrics.norm_pearson, pos_pair=("Actb", "Gapdh"), neg_pair=("G0", "G1"),
                      signal_test=injected_signal.coupling("Smad3", "Col1a1"),
                      prereg=dict(estimand="composition", direction="decrease", sesoi=0.1), log_path="off")
    wall = time.time() - t1
    print(f"run_autopsy: {wall:.0f} s wall; peak resident memory of the process {peak_mb():.0f} MB")
    print("gates: " + "; ".join(f"GATE {r.gate} {r.status.value}" for r in res.results))
    spec = res.specificity or {}
    d_ = spec.get("in_direction")
    print("specificity: " + (f"{d_['k']} of {d_['n']} controls at least as far in the claim's direction, {d_['outcome']}"
                             if d_ else str(spec.get("error", spec.get("reason", "not run")))))
    print(f"VERDICT: {res.verdict}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
