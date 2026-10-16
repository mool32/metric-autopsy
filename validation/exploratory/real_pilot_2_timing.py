#!/usr/bin/env python3
"""The time and memory of one claim of the pilot's second run, case 1, on fabricated counts of its
size (exploratory; real_pilot_plan.md, "Deviation: a second run of cases 1 and 2"). Nothing here is
real data or a result: it sized the workflow's jobs before the run.

The case: the flagship audit's 17 donor units with their cells (audit_results.json, E.mice: 76,797
cells), 2,000 genes, Poisson counts around a gamma-distributed mean, Smad3 and Col1a1 sharing a
per-cell factor, old males at 0.6 of the depth. The claims and settings are the second run's case 1
(1a, 1a-floor, 1b). `make` writes the case and a select.json; `run` runs real_pilot.py on it as its
__main__ and prints the wall time and the peak resident memory (getrusage). The measurement in
real_pilot_2_timing.log ran the whole of case 1 in one process and was stopped after claim 1a;
memory was also sampled every 5 s with ps.

    python validation/exploratory/real_pilot_2_timing.py make OUT/data audit_results.json
    python validation/exploratory/real_pilot_2_timing.py run validation/exploratory/real_pilot.py \
        --data OUT/data --out OUT/out --part case1a
"""
from __future__ import annotations

import sys


def make(out_dir: str, audit_path: str) -> int:
    import hashlib
    import json
    from pathlib import Path

    import numpy as np
    import pandas as pd
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "prereg"))
    import panel as P

    out = Path(out_dir)
    audit = json.loads(Path(audit_path).read_text())["E"]["mice"]
    rng = np.random.default_rng(1)
    g = 2000
    base = np.exp(rng.normal(-9.5, 1.6, g))
    base /= base.sum()
    blocks, rows = [], []
    for m in audit:
        n = int(m["n_cells"])
        age = "young" if m["age"] == "3m" else "old"
        lib = np.exp(rng.normal(np.log(3e5), 0.6, n)) * (0.6 if (age == "old" and m["sex"] == "male") else 1.0)
        z = rng.normal(0, 1, n)
        lam = np.outer(lib, base) * rng.gamma(2.0, 0.5, (n, g)) * 0.2
        lam[:, 0] *= np.exp(0.6 * z)
        lam[:, 1] *= np.exp(0.6 * z)
        blocks.append(rng.poisson(lam).astype(np.int32))
        rows += [dict(mouse=m["donor_id"], age=age, sex=m["sex"])] * n
    X = np.vstack(blocks)
    genes = ["Smad3", "Col1a1", "Actb", "Gapdh", "Tst", "Lrrc42"] + [f"G{j}" for j in range(g - 6)]
    print("zeros", float((X == 0).mean()), "cells", X.shape)
    out.mkdir(parents=True, exist_ok=True)
    P.save_npz(out / "case1.npz", X, pd.DataFrame(rows), genes)
    sha = hashlib.sha256((out / "case1.npz").read_bytes()).hexdigest()
    base_claim = dict(genes=["Smad3", "Col1a1"], direction="decrease", signal=dict(kind="coupling", strength=1.0),
                      controls=dict(pos=["Actb", "Gapdh"], neg=["Tst", "Lrrc42"]), skip="")
    claims = [dict(base_claim, id="1a", metric="mi_3bin", sesoi=0.1),
              dict(base_claim, id="1a-floor", metric="mi_3bin", sesoi="floor", sesoi_scale="observed", seed_id="1a"),
              dict(base_claim, id="1b", metric="norm_pearson", sesoi=0.1)]
    spec = dict(synthetic=True, run="second", release="SYNTHETIC", cases=dict(case1=dict(
        status="run", data=dict(file="case1.npz", sha256=sha), group_col="age", groups=["young", "old"],
        within=["sex"], replicate="mouse", claims=claims, audit_mice=audit)))
    (out / "select.json").write_text(json.dumps(spec, indent=1))
    print("npz MB", (out / "case1.npz").stat().st_size / 1e6)
    return 0


def run(script: str, args: list) -> int:
    """`script` as __main__ with `args`; the wall time and peak RSS to stderr (stdlib only, so that
    the script's own numerical settings come before numpy)."""
    import resource
    import runpy
    import time
    t0 = time.time()
    sys.argv = [script, *args]
    code = 0
    try:
        runpy.run_path(script, run_name="__main__")
    except SystemExit as e:
        code = e.code
    print(f"WALL {time.time() - t0:.1f} s; MAXRSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.2f} GB; "
          f"exit {code}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    if sys.argv[1] == "make":
        raise SystemExit(make(sys.argv[2], sys.argv[3]))
    raise SystemExit(run(sys.argv[2], sys.argv[3:]))
