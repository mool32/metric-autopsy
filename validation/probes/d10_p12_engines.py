#!/usr/bin/env python3
"""D10 (JOURNAL.md): the new p12 test's first check against the engines of two checkouts.

    python -I validation/probes/d10_p12_engines.py SRC [SRC ...] > validation/probes/d10_p12_engines.log

Each SRC is a checkout's ``src`` (a worktree of another commit, or this checkout's); each runs in its own
``python -I`` process with that ``src`` first, and stops if another engine was imported. GATE 1 on p12's data
(a 2,000 x 500 sparse AnnData, a dense budget of 0.001 GB): (1) with a ``DenseMemoryWarning`` made an error;
(2) with it allowed, the matrix wrapped in a CSR class whose ``toarray`` refuses a dense copy of more than a
block of 256 rows by more than 256 columns.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def one(src: str) -> None:
    sys.path.insert(0, src)
    import os
    import warnings

    import pandas as pd
    import scipy.sparse as sparse
    import anndata

    import metric_autopsy as ma
    from metric_autopsy import core
    from metric_autopsy.gates import gate1_qc_parity

    if Path(ma.__file__).resolve().parent.parent != Path(src).resolve():
        sys.exit(f"the engine imported is {ma.__file__}, not {src}'s")
    os.environ["METRIC_AUTOPSY_DENSE_WARN_GB"] = "0.001"
    X = sparse.random(2000, 500, density=0.02, format="csr", random_state=0)
    obs = pd.DataFrame({"age": ["young", "old"] * 1000}, index=[f"c{i}" for i in range(2000)])
    var = pd.DataFrame(index=["mt-Co1"] + [f"g{i}" for i in range(499)])

    def whole(shape):
        return shape[0] > 256 and shape[1] > 256

    class Guarded(sparse.csr_matrix):
        def toarray(self, *a, **k):
            if whole(self.shape):
                raise AssertionError(f"the whole matrix densified: {self.shape}")
            return super().toarray(*a, **k)

    ad = anndata.AnnData(X=Guarded(X), obs=obs, var=var)
    for label, strict in (("(1) DenseMemoryWarning an error", True), ("(2) the warning allowed", False)):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error" if strict else "ignore", core.DenseMemoryWarning)
                r = gate1_qc_parity(ad, "age", ("young", "old"))
            print(f"  {label}: GATE 1 {r.status.value}, " + ("no warning and no dense copy of the matrix" if strict
                  else "no dense copy of the matrix"))
        except Exception as e:  # what the check records
            print(f"  {label}: {type(e).__name__}: {str(e).splitlines()[0][:120]}")


def main() -> int:
    if len(sys.argv) > 2 and sys.argv[1] == "--one":
        one(sys.argv[2])
        return 0
    import platform

    print("# D10: p12's first check (GATE 1 on a 2,000 x 500 sparse AnnData, budget 0.001 GB) against each engine; "
          f"python {platform.python_version()}")
    for src in sys.argv[1:]:
        rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=src, capture_output=True, text=True)
        dirty = subprocess.run(["git", "status", "--porcelain", "--", "."], cwd=src, capture_output=True, text=True)
        print(f"engine at git {rev.stdout.strip()}{' (src modified)' if dirty.stdout.strip() else ''}:")
        out = subprocess.run([sys.executable, "-I", str(Path(__file__).resolve()), "--one", src],
                             capture_output=True, text=True)
        print(out.stdout.rstrip())
        if out.returncode:
            print(f"  exit {out.returncode}: {out.stderr.strip().splitlines()[-1] if out.stderr.strip() else ''}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
