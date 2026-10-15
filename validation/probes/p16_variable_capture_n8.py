"""N8, per-cell variable capture, through the whole engine (the second review predicted that S1
fails there; the third asked for a script of the repository behind the numbers v1.md quotes).

N8 (validation/prereg/v1.md 3.2) thins every cell of one side by its own Beta(2, 2) capture: an
artifact outside the engine's correction family, which equalizes depth by one common ratio per
group. On the panel's simulated backgrounds (validation/prereg/simulate.py, dry_backgrounds: B1 of
24 donors, 2,000 kept genes; datasets of 2 x 8 donors x 200 cells), with the panel's own datasets
and claim cards (panel.build; pairs in turn over the whole pool; simulate.dry_pilot's stand-in
SESOI 0.1 and dose 2.0), the engine runs every card as the blind run does (run_panel.run). Per
condition the probe counts the outcomes and the false SUPPORTED (on null data every SUPPORTED is
false), against S1's allowance at the panel's 790 datasets of a key condition (oc.error_rule). N1,
the same design without an artifact, is the reference. The oracle's analysis of the effect with
the true correction (oracle.effect_outcome: it undoes the planted capture) on the same datasets
shows what a correction of the right family would give.

An indication only: the backgrounds are simulated, and the real run's SESOI and doses come from
its own pilot on the real B1.

    python p16_variable_capture_n8.py > p16_variable_capture_n8.log   # ~1 h at 3 workers
"""
import json
import platform
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "prereg"))
import run_panel as R  # noqa: E402  (sets the numerical environment before numpy loads)
import numpy as np  # noqa: E402
import oc  # noqa: E402
import oracle as O  # noqa: E402
import panel as P  # noqa: E402
import simulate  # noqa: E402

DATASETS = int(sys.argv[1]) if len(sys.argv) > 1 else 160  # a smaller number for a smoke test
WORKERS = 3
CONDITIONS = (("N1", "null"), ("N8", "beta(2,2)"))


def entries(name: str, variant: str, n: int, pool: int) -> list[dict]:
    """n datasets of the condition, the pairs in turn over the pool, sides and seeds from a public seed."""
    rng = np.random.default_rng([16, len(name), sum(map(ord, name))])
    return [dict(id=f"{name}-{i:04d}", condition=name, variant=variant, index=i, side=("A", "B")[int(rng.integers(2))],
                 seed=int(rng.integers(2 ** 62)), pair=i % pool) for i in range(n)]


def oracle_effect(e: dict, bgs: dict, pilot: dict) -> str:
    X, obs, genes, cards = P.build(e, bgs, pilot)
    pe = bgs["B1"].plan["pool"][e["pair"]]
    ia, ib = genes.index(pe["pair"][0]), genes.index(pe["pair"][1])
    return O.effect_outcome(e, X, obs, cards[-1], ia, ib, np.random.default_rng([e["seed"], 1]))


def main():
    rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, capture_output=True, text=True).stdout.strip()
    print(f"# metric-autopsy probe p16 — git {rev}, python {platform.python_version()}, numpy {np.__version__}")
    bgs = simulate.dry_backgrounds()
    pilot = simulate.dry_pilot(bgs)
    pool = len(bgs["B1"].plan["pool"])
    k_star, _, _ = oc.error_rule(P.KEY_N, oc.E_SUPPORTED)
    print(f"# {DATASETS} datasets per condition on simulated B1 (pool of {pool} pairs, in turn); S1 allows at most "
          f"{k_star} false SUPPORTED of {P.KEY_N} datasets ({k_star / P.KEY_N:.1%}) per key condition")
    for name, variant in CONDITIONS:
        es = entries(name, variant, DATASETS, pool)
        with tempfile.TemporaryDirectory() as tmp:
            R.run(es, bgs, pilot, Path(tmp), workers=WORKERS)
            got = Counter()
            for e in es:
                rep = json.loads((Path(tmp) / "reports" / f"{e['id']}.json").read_text())
                got[P.outcome(rep.get("verdict"), rep.get("cause")) if "error" not in rep else P.ERROR] += 1
        by_level = Counter()
        for e in es:
            by_level[bgs["B1"].plan["pool"][e["pair"]]["level"]] += 1
        oracle = Counter(oracle_effect(e, bgs, pilot) for e in es)
        sup, ora = got.get(P.SUPPORTED, 0), oracle.get(P.SUPPORTED, 0)
        print(f"{name} {variant}: engine {dict(sorted(got.items()))}")
        print(f"   false SUPPORTED by the engine {sup}/{DATASETS} = {sup / DATASETS:.3f} "
              f"({'above' if sup / DATASETS > k_star / P.KEY_N else 'within'} S1's {k_star / P.KEY_N:.3f}); "
              f"the oracle's effect analysis with the true correction: SUPPORTED {ora}/{DATASETS} = {ora / DATASETS:.3f}, "
              f"{dict(sorted(oracle.items()))}; pairs by level {dict(sorted(by_level.items()))}", flush=True)


if __name__ == "__main__":
    main()
