"""Timing pilot for the compute plan (validation/prereg/v1.md, section 4).

Runs the frozen engine through run_panel.py on a sample of claim cards built by panel.py from a
simulated background of B1's planned size (2 x 8 donors, 200 cells each, 2,000 genes), with the
condition mix of the panel, at 1 worker and at W workers, and extrapolates to the whole panel.
Re-run on the real B1 once it is downloaded (section 8, step 2).

    python validation/prereg/timing.py --workers 4 --cards 24 > validation/prereg/timing.log
"""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel as P  # noqa: E402
import run_panel as R  # noqa: E402
from test_prereg import simulated_background  # noqa: E402


def sample(entries, k: int, rng) -> list:
    """k entries with the panel's mix of conditions (probability proportional to cards)."""
    w = np.array([P.conditions()[e["condition"]].cards for e in entries], float)
    pick = rng.choice(len(entries), size=k, replace=False, p=w / w.sum())
    return [entries[i] for i in sorted(pick)]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--cards", type=int, default=24)
    p.add_argument("--genes", type=int, default=2500)
    args = p.parse_args(argv)
    rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, capture_output=True,
                         text=True).stdout.strip() or "?"
    print(f"# timing pilot — git {rev}, python {platform.python_version()}, numpy {np.__version__}")
    bgs = {"B1": simulated_background("B1", donors=24, cells=220, n_genes=args.genes),
           "B2": simulated_background("B2", donors=12, cells=220, n_genes=args.genes, seed=1)}
    pilot = dict(sesoi=0.1, key_dose=2.0)
    entries = sample(P.assign(1), args.cards, np.random.default_rng(0))
    print(f"# sample of {len(entries)} datasets: {dict(Counter(e['condition'] for e in entries))}")
    rows = {}
    with tempfile.TemporaryDirectory() as tmp:
        P.write_panel(1, bgs, pilot, Path(tmp) / "panel", only=[e["id"] for e in entries])
        shape = np.load(next((Path(tmp) / "panel").glob("*.npz")))["X"].shape
        print(f"# dataset shape (cells x genes): {shape}")
        for w in sorted({1, args.workers}):
            s = R.run(Path(tmp) / "panel", Path(tmp) / f"reports{w}", workers=w)
            rows[w] = s
            per = [json.loads(f.read_text()) for f in (Path(tmp) / f"reports{w}").glob("D*.json")]
            by = {}
            for rep in per:
                cond = next(e["condition"] for e in entries if rep["id"].startswith(e["id"]))
                by.setdefault(cond, []).append(rep["elapsed_seconds"])
            print(f"workers={w}: {s['run']} cards in {s['wall_seconds']:.0f} s wall; per card mean "
                  f"{s['mean_seconds']:.1f} s, median {s['median_seconds']:.1f} s; errors {s['errors']}")
            print("  per condition (mean s): " + ", ".join(f"{c} {np.mean(v):.1f}" for c, v in sorted(by.items())))
    n_cards = sum(P.conditions()[e["condition"]].cards for e in P.assign(1))
    w = args.workers
    per_card_wall = rows[w]["wall_seconds"] / rows[w]["run"]
    print(f"# whole panel: {n_cards} cards; at {w} workers {n_cards * per_card_wall / 3600:.1f} h wall "
          f"({per_card_wall:.1f} s per card); single-core {n_cards * rows[1]['mean_seconds'] / 3600:.1f} CPU-h")


if __name__ == "__main__":
    main()
