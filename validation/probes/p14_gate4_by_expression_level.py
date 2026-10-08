"""GATE 4 on the valid metric by the expression level of the pair (exploratory, found on v0.3.0.dev0).

The confirmatory panel (validation/prereg/v1.md) draws its analysed pairs from three expression
levels. On simulated backgrounds of the panel's planned size (validation/prereg/simulate.py:
negative-binomial counts with cell-size variation, coupled pairs at every level), GATE 4 injects
a coupling of the pair (`injected_signal.coupling`, the claim cards' default strength 1.0, and
2.0) into null datasets of 2 x 8 donors x 200 cells and requires log-normalized Pearson to move
by z >= 3. A metric that responds weakly fails as if it ignored its construct: GATE 4 has no
power rule (GATE 5's silent positive control has one, decided 2026-10-07). At the low level the
metric does not respond at all: log-normalized Pearson is computed on co-detected cells, where
low counts are mostly 1 and their log-normalized values follow the cell's total, so a coupling
that acts through co-detection is invisible to it.

No expected verdict is frozen: what GATE 4 should do with a weak response, and which verdicts
are correct for a real effect the metric cannot see, are put to the project owner.

    python p14_gate4_by_expression_level.py > p14_gate4_by_expression_level.log   # ~5 min
"""
import platform
import subprocess
import sys
from functools import partial
from pathlib import Path

from probe_common import *  # noqa: F401,F403  (repo path, numpy, pandas, metrics, gates)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "prereg"))
import panel as PANEL  # noqa: E402
import simulate  # noqa: E402
from metric_autopsy import SimpleData, gate4_signal_response, injected_signal, metrics  # noqa: E402

DATASETS = 12


def main():
    rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, capture_output=True, text=True).stdout.strip()
    print(f"# metric-autopsy probe p14 — git {rev}, python {platform.python_version()}, numpy {np.__version__}")
    bgs = simulate.dry_backgrounds()
    pilot = simulate.dry_pilot(bgs)
    for pe in bgs["B1"].plan["pool"]:
        print(f"# pool {pe['index']:2d} {pe['level']:6s} {pe['pair']} mean {pe['mean'][0]:.2f}/{pe['mean'][1]:.2f} "
              f"detected {pe['detection'][0]:.2f}/{pe['detection'][1]:.2f}")
    print(f"GATE 4 on log-normalized Pearson, null datasets 2 x 8 donors x 200 cells, {DATASETS} per row")
    for strength in (1.0, 2.0):
        for li, lvl in enumerate(PANEL.LEVELS):
            fails, zs, resp = 0, [], []
            for i in range(DATASETS):
                e = dict(id=f"G{i}", condition="N1", variant="null", index=i, side="A",
                         pair=li * PANEL.PAIRS_PER_LEVEL + i % PANEL.PAIRS_PER_LEVEL, seed=1000 + i)
                X, obs, genes, cards = PANEL.build(e, bgs, pilot)
                ga, gb = cards[0]["gene_pair"]
                r = gate4_signal_response(partial(metrics.norm_pearson, gene_a=ga, gene_b=gb),
                                          SimpleData(X, obs, genes),
                                          injected_signal.coupling(ga, gb, strength=strength), seed=i)
                fails += r.status.value == "FAIL"
                zs.append(r.detail.get("z", np.nan))
                resp.append(r.detail.get("mean_response", np.nan))
            print(f"  strength {strength:.1f}, {lvl:6s}: FAIL {fails}/{DATASETS}; z median {np.nanmedian(zs):5.1f} "
                  f"[{np.nanmin(zs):5.1f}, {np.nanmax(zs):5.1f}]; response median {np.nanmedian(resp):+.4f}")


if __name__ == "__main__":
    main()
