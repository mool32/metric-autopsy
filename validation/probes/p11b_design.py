"""p11b design: how many mice make a pure depth artifact *establishable* at the replicate level?

p11 (demo male stratum: identical biology, old-male capture 30%, 100 cells per mouse) reached
its frozen expectation, NOT SUPPORTED — explained by depth, in only 3 of 10 splits of the cells
into 4 vs 4 mice. With 4 vs 4 the exact permutation has 70 assignments and reaches p < 0.05 only
at complete separation, so the raw difference is often not detected and the correct verdict is
INCONCLUSIVE. The probe was therefore split by design (decided 2026-10-07, journal entry D1):

* p11a — 4 vs 4 mice: the truth is not establishable; allowed {NOT SUPPORTED, INCONCLUSIVE};
  only SUPPORTED is an error.
* p11b — N vs N mice: NOT SUPPORTED with the diagnosis "explained by depth" is expected.

Pre-specified rule for N, fixed before this script was run: the smallest N in
{4, 6, 8, 10, 12, 16} for which the lower 95% Clopper-Pearson bound of
P(NOT SUPPORTED and explained_by_depth) over 40 independent datasets (data seeds 0-39, each
with its own split into mice) is >= 0.80. Same generator as the demo (`demo_data`), 100 cells
per mouse, male stratum only, mi_3bin on Smad3-Col1a1, composition estimand.

    python validation/probes/p11b_design.py > validation/probes/p11b_design.log    # ~10 min
"""
from __future__ import annotations

import platform
import subprocess
import sys
from functools import partial
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from metric_autopsy import __version__, metrics, run_autopsy  # noqa: E402
from metric_autopsy.cli import demo_data  # noqa: E402
from metric_autopsy.stats import clopper_pearson, fmt_rate  # noqa: E402

MI = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")
CANDIDATES = (4, 6, 8, 10, 12, 16)
N_DATASETS = 40
CELLS_PER_MOUSE = 100
LOWER_BOUND_MIN = 0.80


def male_stratum(n_mice: int, seed: int):
    d = demo_data(seed=seed, n=CELLS_PER_MOUSE * n_mice, mice_per_block=n_mice)
    return d[np.asarray(d.obs["sex"]) == "male"]


def run(n_mice: int, seed: int):
    a = run_autopsy(MI, male_stratum(n_mice, seed), group_col="age", groups=("young", "old"),
                    replicate_col="mouse", prereg={"estimand": "composition"}, log_path="off")
    e = a.effect.detail
    return dict(explained=bool(e.get("explained_by_depth")), raw_detected=bool(e.get("raw_detected")),
                verdict=a.verdict.split(" — ")[0], supported=a.verdict.startswith("SUPPORTED"))


def main():
    rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, capture_output=True,
                         text=True).stdout.strip() or "?"
    print(f"# p11b design — metric-autopsy {__version__}, git {rev}, python {platform.python_version()}, "
          f"numpy {np.__version__}")
    print(f"# rule: smallest N with lower 95% CP bound of P(explained by depth) >= {LOWER_BOUND_MIN} "
          f"over {N_DATASETS} datasets; {CELLS_PER_MOUSE} cells per mouse")
    chosen = None
    for n in CANDIDATES:
        rows = [run(n, s) for s in range(N_DATASETS)]
        k = sum(r["explained"] for r in rows)
        kd = sum(r["raw_detected"] for r in rows)
        ks = sum(r["supported"] for r in rows)
        lo, _ = clopper_pearson(k, N_DATASETS)
        print(f"N={n:2d} vs {n:2d}: explained by depth {fmt_rate(k, N_DATASETS)}; raw detected "
              f"{fmt_rate(kd, N_DATASETS)}; SUPPORTED {fmt_rate(ks, N_DATASETS)}", flush=True)
        if chosen is None and lo >= LOWER_BOUND_MIN:
            chosen = n
    print(f"CHOSEN N = {chosen}")


if __name__ == "__main__":
    main()
