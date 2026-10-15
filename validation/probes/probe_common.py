"""Shared helpers for the validator probes. Finds the repo (for tests/test_gates.py generators)."""
import os, sys
from pathlib import Path
REPO = Path(os.environ.get("METRIC_AUTOPSY_REPO") or next(
    (q for q in Path(__file__).resolve().parents if (q / "src" / "metric_autopsy").exists()),
    Path("/home/user/metric-autopsy")))
sys.path.insert(0, str(REPO / "src")); sys.path.insert(0, str(REPO / "tests"))
from functools import partial
import numpy as np, pandas as pd
from metric_autopsy import (SimpleData, GateStatus, metrics, run_autopsy,
    gate0_independence, gate1_qc_parity, gate2_ngenes_matching, gate5_controls)
from test_gates import make_clean, make_confounded, GENES

def show(autopsy):
    for r in sorted(autopsy.results, key=lambda r: r.gate):
        print(f"  G{r.gate} {r.status.value:8s} {r.message[:150]}")
    print("  VERDICT:", autopsy.verdict)
