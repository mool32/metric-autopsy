"""Timing pilot for the compute plan (validation/prereg/v1.md, section 4), and GATE 0's refusals.

Runs the frozen engine through run_panel.py on a sample of the panel's datasets, built on the
fly by panel.py, stratified by condition (one dataset of every condition, the rest in proportion
to the conditions' cards), at 1 worker and at W workers, and extrapolates to the whole panel
condition by condition: every card at its condition's mean time. Default: simulated backgrounds
of the planned sizes (simulate.py). In the pilot (section 8, step 3) it is re-run on the real
backgrounds with the pilot's SESOI and key dose, at the workers the run will use, and decides
the drops of section 4. With --refusals N it then runs the engine on N establishable cards of
every S3 stratum (``refusal_sample``) and writes GATE 0's refusal share per stratum into the
pilot (the fourth review: the sound validator's model left GATE 0 out, so a sound engine could
fail S3; the oracle has no copy of GATE 0, so its share is the one model input measured with the
engine itself, before the key, on datasets the key does not draw).

    python validation/prereg/timing.py --workers 4 --cards 24 > validation/prereg/timing.log
    python validation/prereg/timing.py --workers W --backgrounds backgrounds.json --pilot pilot.json \
        --write-drops --refusals 100
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import frozen  # noqa: E402  (standard library only)

frozen.pin_numerics()  # as run_panel.py: one BLAS thread per worker and one numerical path, before numpy

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel as P  # noqa: E402
import run_panel as R  # noqa: E402
import simulate  # noqa: E402

TIMING_KEY = "7" * 64  # public: the timing sample is not the panel
REFUSAL_KEY = "8" * 64  # public: nor is the refusal sample
SHARDS, WORKERS = 20, 4  # the blind run: 20 shard jobs of 4 workers (validation.yml)
BUDGET_HOURS = 3.0       # a shard's expected time may not exceed this (a job may run 6 h)
MEMORY_SHARE = 0.75      # the run's workers at a worker's measured peak may use at most this share of the runner's memory


def runner_memory_mb() -> float | None:
    """The machine's physical memory in MB (None where the system does not say)."""
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 2 ** 20
    except (ValueError, OSError, AttributeError):
        return None


def check_memory(peak_mb: float, total_mb: float | None, workers: int = WORKERS, enforce: bool = True) -> str:
    """The memory rule of v1.md section 4 (the third review: the drops weighed time only, and a
    worker killed for memory after the key would leave no results): `workers` workers at the
    sample's peak memory of one worker (every condition is in the sample, N7 at its largest
    design) must fit in MEMORY_SHARE of the runner's memory; else, where the rule is enforced (the
    pilot step), the pilot stops before the key exists. Returns the line for the log."""
    need = peak_mb * workers
    if total_mb is None:
        if enforce:
            raise SystemExit("the runner's memory is unknown: the memory rule cannot be checked")
        return f"# memory: a worker's peak {peak_mb:.0f} MB x {workers} workers = {need:.0f} MB; the machine's memory is unknown"
    line = (f"# memory: a worker's peak {peak_mb:.0f} MB x {workers} workers = {need:.0f} MB; the runner has "
            f"{total_mb:.0f} MB (at most {MEMORY_SHARE:.0%} may be used)")
    if enforce and need > MEMORY_SHARE * total_mb:
        raise SystemExit(line[2:] + ": too much, the run stops before the key")
    return line


def shard_hours(seconds_per_card: dict, dropped=()) -> float:
    """A shard's expected wall time: every card of the design (after the drops) at its condition's
    mean time per card measured at WORKERS workers, spread over SHARDS shards of WORKERS workers
    (the key's order mixes the conditions over the shards). A condition without a measurement
    takes the mean over the measured ones."""
    fallback = float(np.mean(list(seconds_per_card.values())))
    total = sum(n * c.cards * seconds_per_card.get(c.name, fallback) for c in P.CONDITIONS
                for v, n in c.variants if not P._dropped(c.name, v, dropped))
    return total / WORKERS / SHARDS / 3600


def decide_drops(seconds_per_card: dict, budget_hours: float = BUDGET_HOURS, dropped=()) -> list:
    """Drops in the pre-registered order (panel.DROP_ORDER) until a shard's expected time
    (`shard_hours`) fits the budget (v1.md, section 4), after the drops already made (`dropped`: a
    background without a candidate, with its cases). Returns the order's drops."""
    before, out = list(dropped), []
    for name, variant in P.DROP_ORDER:
        if shard_hours(seconds_per_card, before + out) <= budget_hours:
            break
        out.append(f"{name}:{variant}")
    if shard_hours(seconds_per_card, before + out) > budget_hours:
        raise SystemExit("the key conditions and the real effects at the key dose do not fit the budget")
    return out


def sample(entries, k: int, rng) -> list:
    """k entries stratified by condition: one of every condition among the entries, the rest drawn
    with probability proportional to the condition's cards, so that every condition is timed."""
    conds = P.conditions()
    by = {}
    for i, e in enumerate(entries):
        by.setdefault(e["condition"], []).append(i)
    if k < len(by):
        raise SystemExit(f"--cards {k}: the sample takes one dataset of each of the {len(by)} conditions")
    pick = [int(rng.choice(by[c])) for c in sorted(by)]
    rest = [i for i in range(len(entries)) if i not in set(pick)]
    w = np.array([conds[entries[i]["condition"]].cards for i in rest], float)
    pick += [int(i) for i in rng.choice(rest, size=k - len(pick), replace=False, p=w / w.sum())]
    return [entries[i] for i in sorted(pick)]


def refusal_sample(pilot: dict, n: int, rng) -> list:
    """n datasets of every S3 stratum (fewer where it has fewer): drawn uniformly from the datasets
    of a public key (REFUSAL_KEY, after the drops) whose card is establishable in that stratum
    (``oc.card_model``), so in proportion to the conditions' share of the stratum's cards. Not the
    constant metric: the engine decides DEGENERATE METRIC before GATE 0."""
    import oc
    conds = P.conditions()
    by = {st: [] for st in oc.STRATA}
    for e in P.assign(REFUSAL_KEY, pilot.get("dropped", ()), pilot.get("pool_size")):
        c = conds[e["condition"]]
        row = oc.card_model(c, e["variant"], int(e["pair"]), pilot)
        if row["establishable"] and row["stratum"] in by and c.metric != "constant":
            by[row["stratum"]].append(e)
    out = []
    for st in oc.STRATA:
        if by[st]:
            pick = rng.choice(len(by[st]), size=min(n, len(by[st])), replace=False)
            out += [dict(by[st][int(i)], stratum=st) for i in sorted(pick)]
    return out


def measure_refusals(entries: list, bgs: dict, pilot: dict, out_dir: Path, workers: int) -> dict:
    """GATE 0's refusals on `refusal_sample`'s datasets, run by the frozen engine as in the blind
    run: per S3 stratum the number of cards, of refusals and their share, and every outcome's
    count (a crash included)."""
    R.run([{k: v for k, v in e.items() if k != "stratum"} for e in entries], bgs, pilot, out_dir, workers=workers)
    out = {}
    for e in entries:
        rec = out.setdefault(e["stratum"], dict(n=0, refused=0, outcomes={}))
        for cid in P.card_ids(e):
            path = out_dir / "reports" / f"{cid}.json"
            rep = json.loads(path.read_text()) if path.exists() else {}
            o = P.CRASH if "error" in rep or not rep else P.outcome(rep.get("verdict"), rep.get("cause"))
            rec["n"] += 1
            rec["refused"] += o == P.REFUSAL
            rec["outcomes"][o] = rec["outcomes"].get(o, 0) + 1
    for rec in out.values():
        rec["share"] = rec["refused"] / rec["n"]
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--cards", type=int, default=24)
    p.add_argument("--genes", type=int, default=2500, help="genes of the simulated backgrounds")
    p.add_argument("--backgrounds", help="backgrounds JSON as panel.py takes it (default: simulated)")
    p.add_argument("--pilot", help="pilot.json from oracle.py (default: a stand-in, SESOI 0.1, saturation dose 2.0)")
    p.add_argument("--data-dir", help="the downloaded files of backgrounds.json (default: next to it)")
    p.add_argument("--write-drops", action="store_true",
                   help="write the drops the projection needs into --pilot (the pilot step, before the key)")
    p.add_argument("--refusals", type=int, default=0,
                   help="then measure GATE 0's refusals on this many establishable cards per S3 stratum, after "
                        "the drops; with --write-drops written into --pilot")
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
    entries = sample(P.assign(TIMING_KEY, pilot.get("dropped", ()), pilot.get("pool_size")), args.cards,
                     np.random.default_rng(0))
    level = {(b, pe["index"]): pe["level"] for b, pool in pilot["pool"].items() for pe in pool}
    print(f"# sample of {len(entries)} datasets: {dict(Counter(e['condition'] for e in entries))}; levels "
          f"{dict(Counter(level[(P.conditions()[e['condition']].background, e['pair'])] for e in entries))}")
    print(f"# dataset shape (cells x genes): {P.build(entries[0], bgs, pilot)[0].shape}")
    rows, sec, peaks = {}, {}, {}
    with tempfile.TemporaryDirectory() as tmp:
        for w in sorted({1, args.workers}):
            s = R.run(entries, bgs, pilot, Path(tmp) / f"out{w}", workers=w)
            rows[w] = s
            manifest = json.loads((Path(tmp) / f"out{w}" / "manifest.json").read_text())
            by = {}
            cond_of = {e["id"]: e["condition"] for e in entries}
            runtime = {c["id"]: c for c in json.loads((Path(tmp) / f"out{w}" / "runtime.json").read_text())["cards"]}
            for row in manifest["datasets"]:
                for c in row["cards"]:
                    by.setdefault(cond_of[row["id"]], []).append(runtime[c["id"]]["seconds"])
            sec[w] = {c: float(np.mean(v)) for c, v in by.items()}
            peak = peaks[w] = max(c.get("peak_rss_mb") or 0 for c in runtime.values())
            print(f"workers={w}: {s['run']} cards in {s['wall_seconds']:.0f} s wall; per card mean "
                  f"{s['mean_seconds']:.1f} s, median {s['median_seconds']:.1f} s; errors {s['errors']}; "
                  f"peak memory of a worker {peak:.0f} MB")
            print("  per condition (mean s): " + ", ".join(f"{c} {v:.1f}" for c, v in sorted(sec[w].items())))
    dropped = pilot.get("dropped", ())
    n_cards = P.n_cards(dropped)
    w = args.workers
    cpu = shard_hours(sec[w], dropped) * WORKERS * SHARDS
    print(f"# whole panel: {n_cards} cards, each at its condition's mean time: {cpu:.1f} CPU-h at {w} workers "
          f"(single-core {shard_hours(sec[1], dropped) * WORKERS * SHARDS:.1f} CPU-h)")
    print(f"# blind run: {SHARDS} shards of {WORKERS} workers, a shard's expected time "
          f"{shard_hours(sec[w], dropped):.2f} h (budget {BUDGET_HOURS} h per shard; a job may run 6 h)")
    if args.write_drops and w != WORKERS:
        raise SystemExit(f"the drops are decided at the run's {WORKERS} workers")
    print(check_memory(peaks[w], runner_memory_mb(), WORKERS, enforce=args.write_drops))
    if args.write_drops:
        backgrounds = [d for d in pilot.get("dropped", ()) if d.endswith(":*")]  # the pilot's: kept
        dropped = backgrounds + decide_drops(sec[w], dropped=backgrounds)
        pilot["dropped"] = P.check_dropped(dropped)
        Path(args.pilot).write_text(json.dumps(pilot, indent=1))
        print(f"# drops written to {args.pilot}: {dropped or 'none'}")
    if args.refusals > 0:
        entries = refusal_sample(pilot, args.refusals, np.random.default_rng(1))
        with tempfile.TemporaryDirectory() as tmp:
            rec = measure_refusals(entries, bgs, pilot, Path(tmp) / "refusals", args.workers)
        for st, r in rec.items():
            print(f"# GATE 0 refusals, S3 stratum {st}: {r['refused']} of {r['n']} establishable cards "
                  f"({r['share']:.3f}); outcomes {dict(sorted(r['outcomes'].items()))}")
        if args.write_drops:
            pilot["gate0_refusals"] = rec
            Path(args.pilot).write_text(json.dumps(pilot, indent=1))
            print(f"# GATE 0's refusal shares written to {args.pilot}")


if __name__ == "__main__":
    main()
