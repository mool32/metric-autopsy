"""Timing pilot for the compute plan (validation/prereg/v1.md, section 4).

Runs the frozen engine through run_panel.py on a sample of the panel's datasets, built on the
fly by panel.py, with the condition mix of the panel, at 1 worker and at W workers, and
extrapolates to the whole panel. Default: simulated backgrounds of the planned sizes
(simulate.py). In the pilot (section 8, step 3) it is re-run on the real backgrounds with the
pilot's SESOI and key dose, at the workers the run will use.

    python validation/prereg/timing.py --workers 4 --cards 24 > validation/prereg/timing.log
    python validation/prereg/timing.py --workers W --backgrounds backgrounds.json --pilot pilot.json
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")  # as run_panel.py: one BLAS thread per worker, set before numpy

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel as P  # noqa: E402
import run_panel as R  # noqa: E402
import simulate  # noqa: E402

TIMING_KEY = "7" * 32  # public: the timing sample is not the panel


def sample(entries, k: int, rng) -> list:
    """k entries with the panel's mix of conditions (probability proportional to cards)."""
    w = np.array([P.conditions()[e["condition"]].cards for e in entries], float)
    pick = rng.choice(len(entries), size=k, replace=False, p=w / w.sum())
    return [entries[i] for i in sorted(pick)]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--cards", type=int, default=24)
    p.add_argument("--genes", type=int, default=2500, help="genes of the simulated backgrounds")
    p.add_argument("--backgrounds", help="backgrounds JSON as panel.py takes it (default: simulated)")
    p.add_argument("--pilot", help="pilot.json from oracle.py (default: a stand-in, SESOI 0.1, dose 2.0)")
    p.add_argument("--data-dir", help="the downloaded files of backgrounds.json (default: next to it)")
    args = p.parse_args(argv)
    rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, capture_output=True,
                         text=True).stdout.strip() or "?"
    print(f"# timing pilot — git {rev}, python {platform.python_version()}, numpy {np.__version__}")
    print(f"# machine: {json.dumps(R.machine())}")
    if args.backgrounds:
        bgs = P.load_backgrounds(args.backgrounds, args.data_dir)
        print(f"# backgrounds: {args.backgrounds} (sha256 {P.sha256(Path(args.backgrounds))})")
    else:
        bgs = simulate.dry_backgrounds(args.genes)
        print("# backgrounds: simulated (validation/prereg/simulate.py, dry_backgrounds)")
    pilot = json.loads(Path(args.pilot).read_text()) if args.pilot else simulate.dry_pilot(bgs)
    entries = sample(P.assign(TIMING_KEY, pilot.get("dropped", ())), args.cards, np.random.default_rng(0))
    print(f"# sample of {len(entries)} datasets: {dict(Counter(e['condition'] for e in entries))}; "
          f"levels {dict(Counter(P.level_of_pair(e['pair']) for e in entries))}")
    print(f"# dataset shape (cells x genes): {P.build(entries[0], bgs, pilot)[0].shape}")
    rows = {}
    with tempfile.TemporaryDirectory() as tmp:
        for w in sorted({1, args.workers}):
            s = R.run(entries, bgs, pilot, Path(tmp) / f"out{w}", workers=w)
            rows[w] = s
            manifest = json.loads((Path(tmp) / f"out{w}" / "manifest.json").read_text())
            by = {}
            cond_of = {e["id"]: e["condition"] for e in entries}
            for row in manifest["datasets"]:
                for c in row["cards"]:
                    by.setdefault(cond_of[row["id"]], []).append(c["seconds"])
            print(f"workers={w}: {s['run']} cards in {s['wall_seconds']:.0f} s wall; per card mean "
                  f"{s['mean_seconds']:.1f} s, median {s['median_seconds']:.1f} s; errors {s['errors']}")
            print("  per condition (mean s): " + ", ".join(f"{c} {np.mean(v):.1f}" for c, v in sorted(by.items())))
    n_cards = P.n_cards(pilot.get("dropped", ()))
    w = args.workers
    per_card_wall = rows[w]["wall_seconds"] / rows[w]["run"]
    print(f"# whole panel: {n_cards} cards; at {w} workers {n_cards * per_card_wall / 3600:.1f} h wall "
          f"({per_card_wall:.1f} s per card); single-core {n_cards * rows[1]['mean_seconds'] / 3600:.1f} CPU-h")


if __name__ == "__main__":
    main()
