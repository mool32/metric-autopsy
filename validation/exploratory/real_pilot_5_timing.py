#!/usr/bin/env python3
"""The time and the peak memory of the third diagnostics' heaviest part, the acceptance test A4 (§4 within
sex, case 1's size), on fabricated counts (nothing real): the 17 units of case 1 with their numbers of
cells, 40 cell types, D1's 2,000 genes dense and a pool of 21,025 genes sparse, run through
real_pilot_5.accept_run exactly as the part `case1` runs it.

    python validation/exploratory/real_pilot_5_timing.py [--scale 1.0] [--claim A4]

Prints the sizes, the seconds and the peak resident memory (ru_maxrss), for the workflow's limits (a job
of at most 350 minutes on a runner of 16 GB).
"""
from __future__ import annotations

import argparse
import resource
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import real_pilot_5 as R  # noqa: E402  (pins the numerics before numpy loads)
import real_pilot_5_spec as S  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import sparse  # noqa: E402

UNITS = {  # case 1's units and their cells (results/real-pilot-2, select.json)
    "3_38_F": ("female", "young", 7814), "3_38_F/3_39_F": ("female", "young", 75), "3_39_F": ("female", "young", 7985),
    "3_56_F": ("female", "young", 3041), "3_57_F": ("female", "young", 400), "21_48_F": ("female", "old", 222),
    "21_54_F": ("female", "old", 267), "21_55_F": ("female", "old", 239), "3_10_M": ("male", "young", 7244),
    "3_10_M/3_11_M": ("male", "young", 63), "3_11_M": ("male", "young", 3817), "3_8_M": ("male", "young", 7071),
    "3_9_M": ("male", "young", 7008), "24_58_M": ("male", "old", 8031), "24_59_M": ("male", "old", 7897),
    "24_60_M": ("male", "old", 8353), "24_61_M": ("male", "old", 7270)}


def fabricate(scale: float, seed: int = 0):
    rng = np.random.default_rng(seed)
    rows = []
    for u, (sex, age, n) in UNITS.items():
        rows += [dict(mouse=u, sex=sex, age=age)] * max(1, int(round(n * scale)))
    obs = pd.DataFrame(rows)
    n = len(obs)
    types = [f"fibroblast {k}" for k in range(4)] + [f"type {k}" for k in range(36)]
    obs["cell_type"] = rng.choice(types, size=n, p=np.r_[np.full(4, 0.025), np.full(36, 0.9 / 36)])
    obs["tissue"] = rng.choice([f"tissue {k}" for k in range(20)], size=n)
    obs["soma_joinid"] = np.arange(n)
    named = ["Smad3", "Col1a1", "Actb", "Gapdh", "Tst", "Lrrc42"]
    genes = named + [f"G{j:04d}" for j in range(2000 - len(named))]
    mu = np.exp(rng.normal(-1.0, 1.5, len(genes)))
    mu[:6] = [0.6, 0.3, 60.0, 40.0, 1.0, 1.0]
    lib = np.exp(rng.normal(0, 0.4, n))
    fib = obs["cell_type"].str.startswith("fibroblast").to_numpy()
    X = np.empty((n, len(genes)), dtype=np.int32)
    for start in range(0, n, 5000):  # in blocks: the Poisson draws of all cells at once would double the memory
        sl = slice(start, min(n, start + 5000))
        lam = np.outer(lib[sl], mu)
        lam[:, 1] *= np.where(fib[sl], 30.0, 0.3)
        X[sl] = rng.poisson(lam)
    extra, per_cell = 21025 - len(genes), 2700
    indptr = np.arange(0, (n + 1) * per_cell, per_cell, dtype=np.int64)
    idx = np.empty(n * per_cell, dtype=np.int32)
    for start in range(0, n, 2000):  # sorted, distinct columns per cell
        stop = min(n, start + 2000)
        block = np.sort(rng.random((stop - start, extra)).argpartition(per_cell, axis=1)[:, :per_cell], axis=1)
        idx[start * per_cell:stop * per_cell] = (block + len(genes)).ravel()
    data = rng.geometric(0.5, size=n * per_cell).astype(np.int32)
    tail = sparse.csr_matrix((data, idx, indptr), shape=(n, len(genes) + extra))
    head = sparse.csr_matrix(X)
    pool = sparse.hstack([head, tail[:, len(genes):]], format="csr").astype(np.int32)
    del tail, head, data, idx
    names = genes + [f"P{j:05d}" for j in range(extra)]
    return X, obs, genes, pool, names


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--scale", type=float, default=1.0)
    p.add_argument("--claim", default="A4", choices=["A3", "A4"])
    args = p.parse_args(argv)
    t0 = time.time()
    X, obs, genes, pool, names = fabricate(args.scale)
    print(f"# SYNTHETIC timing of {args.claim} (real_pilot_5_timing.py; nothing real)")
    print(f"fabricated {X.shape[0]} cells x {X.shape[1]} genes dense, the pool {pool.shape[1]} genes with "
          f"{pool.nnz:,} counts ({time.time() - t0:.0f} s); peak so far {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.2f} GB",
          flush=True)
    claim = next(c for c in S.ACCEPT if c["id"] == args.claim)
    o = obs.astype(str)
    keep = np.ones(len(o), bool) if not claim.get("males") else (o["sex"] == "male").to_numpy()
    data = R.SimpleData(np.asarray(X[keep], dtype=np.float64), o[keep].reset_index(drop=True), genes)
    del X
    pl = R.Pool(pool[np.where(keep)[0]] if claim.get("males") else pool, names, data.obs)
    del pool
    t1 = time.time()
    with tempfile.TemporaryDirectory() as tmp:
        rec = R.accept_run(claim, data, pl, Path(tmp), lambda *m: print(*m, flush=True))
    print(f"{args.claim}: {time.time() - t1:.0f} s in the engine; {rec['verdict'][:160]}")
    print(f"peak resident memory {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.2f} GB; "
          f"total {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
