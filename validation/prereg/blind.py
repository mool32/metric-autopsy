"""The blind run of the confirmatory panel (validation/prereg/v1.md, section 3.3 and 8).

Run by .github/workflows/validation.yml after the pilot: the workflow pushes the run tag
``panel-v1-run``, whose message names a drand round at least an hour later; when the round is
out, ``beacon.py`` fetches and verifies it, and its randomness is the key. Everything here is a
deterministic function of the frozen code, the backgrounds, pilot.json and the key: no dataset
is stored, each is built on the fly, and the manifest records its canonical sha256, so anyone
can rebuild and check every dataset and report from the public key.

    python blind.py guard   --frozen-tag v0.3.0-prereg            # the code is the frozen tag's
    python blind.py prepare --backgrounds backgrounds.json --data-dir DATA --out compact   # verify, plan
    python blind.py run --compact compact --pilot pilot.json --key-record key.json \\
        --run-tag panel-v1-run --shard 0 --shards 20 --workers 4 --out out
    python blind.py collect --shards-dir shards --pilot pilot.json --key-record key.json --out results
    python blind.py all --backgrounds backgrounds.json --data-dir DATA --pilot pilot.json \\
        --key-record key.json --workers 8 --out results       # the same steps on one machine
    python blind.py verify --results RESULTS --compact compact --pilot pilot.json --datasets 20 \\
        --out rerun                 # anyone: re-run datasets of a finished run, compare the reports

``--dry-run`` replaces the backgrounds by simulated ones (simulate.py) and, without a key
record, the key by a public test key; the workflow runs it on pull requests.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import frozen  # noqa: E402  (standard library only)

frozen.pin_numerics()  # one BLAS thread per worker and one numerical path, set before numpy
import panel as P  # noqa: E402
import run_panel as R  # noqa: E402
from frozen import guard  # noqa: E402  (frozen.py, which a job also runs on its own)

ROUND_LINE = re.compile(r"drand (\w+) round:\s*(\d+)")


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=HERE, capture_output=True, text=True, check=True).stdout


# --------------------------------------------------------------------------- #
# the key: the randomness of the drand round named in the run tag
# --------------------------------------------------------------------------- #
def round_from_tag(tag: str) -> tuple[str, int]:
    """The chain and round the run tag's message names ('drand quicknet round: <R>')."""
    m = ROUND_LINE.search(_git("tag", "-l", "--format=%(contents)", tag))
    if not m:
        raise SystemExit(f"no 'drand <chain> round: <R>' line in the message of tag {tag}")
    return m.group(1), int(m.group(2))


def check_key_record(rec: dict, expect_round: tuple[str, int] | None = None) -> str:
    """The key of a key record written by ``beacon.py wait``: the randomness, which must be the
    sha256 of the round's signature (beacon.py verified the signature itself) and belong to the
    round the run tag named."""
    key = P.check_key(rec["key"])
    if key != rec["randomness"] or key != hashlib.sha256(bytes.fromhex(rec["signature"])).hexdigest():
        raise SystemExit("the key record's randomness is not the sha256 of its signature")
    if expect_round is not None and (rec["chain"], int(rec["round"])) != expect_round:
        raise SystemExit(f"the key record is {rec['chain']} round {rec['round']}; the run tag names "
                         f"{expect_round[0]} round {expect_round[1]}")
    return key


def read_key(dry_run: bool, key_record: Path | None, run_tag: str | None) -> tuple[str, dict]:
    if key_record is None:
        if not dry_run:
            raise SystemExit("the blind run takes its key from the beacon's key record (--key-record)")
        import simulate
        return simulate.DRY_RUN_KEY, dict(dry_run=True, key=simulate.DRY_RUN_KEY)
    rec = json.loads(Path(key_record).read_text())
    expect = round_from_tag(run_tag) if run_tag else None
    if expect is None and not dry_run:
        raise SystemExit("the blind run checks the key record against the run tag's round (--run-tag)")
    # the manifest keeps the round and its value, not when or from where it was fetched (key.json has
    # those): the manifest is a function of the code, the backgrounds, the pilot and the round
    return check_key_record(rec, expect), {k: v for k, v in rec.items() if k not in ("answers", "fetched_utc")}


# --------------------------------------------------------------------------- #
# prepare: verify, plan and compact the backgrounds
# --------------------------------------------------------------------------- #
def prepare(backgrounds: Path | None, out: Path, data_dir: Path | None, dry_run: bool = False) -> dict:
    if dry_run:
        import simulate
        bgs = simulate.dry_backgrounds()
    else:
        bgs = P.load_backgrounds(backgrounds, data_dir)  # every file is checked against its sha256
    return write_compact(bgs, out)


def write_compact(bgs: dict, out: Path) -> dict:
    """The planned backgrounds as compact files with their sha256 and content sha256, and the
    digest of the whole (compact/SHA256), which every shard checks."""
    out.mkdir(parents=True, exist_ok=True)
    info = {}
    for name, bg in bgs.items():
        P.save_compact(bg, out / f"{name}.npz")
        info[name] = dict(file=f"{name}.npz", sha256=P.sha256(out / f"{name}.npz"),
                          content_sha256=P.background_sha256(bg), cells=int(bg.X.shape[0]),
                          genes=len(bg.genes), donors=len(bg.donors), pool=len(bg.plan["pool"]))
    (out / "compact.json").write_text(json.dumps(info, indent=1, sort_keys=True))
    digest = hashlib.sha256(json.dumps(info, sort_keys=True).encode()).hexdigest()
    (out / "SHA256").write_text(digest)
    return dict(backgrounds=info, sha256=digest)


def load_prepared(compact: Path, expect: str | None) -> tuple[dict, dict]:
    info = json.loads((compact / "compact.json").read_text())
    digest = hashlib.sha256(json.dumps(info, sort_keys=True).encode()).hexdigest()
    if expect and digest != expect:
        raise SystemExit(f"compact backgrounds {digest} differ from the prepare job's {expect}")
    bgs = {}
    for name, rec in info.items():
        if P.sha256(compact / rec["file"]) != rec["sha256"]:
            raise SystemExit(f"{rec['file']}: sha256 differs from the prepare job's")
        bgs[name] = P.load_compact(compact / rec["file"])
        if P.background_sha256(bgs[name]) != rec["content_sha256"]:
            raise SystemExit(f"{name}: content differs from the prepare job's")
    return bgs, info


# --------------------------------------------------------------------------- #
# run one shard, collect the shards
# --------------------------------------------------------------------------- #
def shard_entries(entries: list[dict], shard: int, shards: int) -> list[dict]:
    """Shard k of S takes the entries k, k + S, ... of the key's (shuffled) order."""
    if not 0 <= shard < shards:
        raise SystemExit(f"shard {shard} outside 0..{shards - 1}")
    return entries[shard::shards]


def pilot_sha256(pilot: dict) -> str:
    return hashlib.sha256(json.dumps(pilot, sort_keys=True).encode()).hexdigest()


def run_shard(compact: Path, pilot: dict, key: str, beacon: dict, shard: int, shards: int, out: Path,
              workers: int, expect: str | None = None, limit: int | None = None) -> dict:
    bgs, info = load_prepared(compact, expect)
    for name, sha in (pilot.get("backgrounds") or {}).items():  # the pilot ran on these backgrounds
        if info.get(name, {}).get("content_sha256") != sha:
            raise SystemExit(f"{name}: the prepared background differs from the one the pilot used")
    entries = P.assign(key, pilot.get("dropped", ()), pilot.get("pool_size"))
    mine = shard_entries(entries, shard, shards)[:limit]
    summary = R.run(mine, bgs, pilot, out, workers)
    meta = dict(shard=shard, shards=shards, key=key, beacon=beacon, entries=len(mine),
                backgrounds=info, pilot_sha256=pilot_sha256(pilot))
    (out / "shard.json").write_text(json.dumps(meta, indent=1, sort_keys=True))
    return dict(meta, summary=summary)


def collect(shards_dir: Path, out: Path, expected_datasets: int | None = None, reruns: list | None = None) -> dict:
    """Merge the shards: every shard ran with the same key, backgrounds and pilot; no dataset or
    card twice. Writes the reports, the manifest (deterministic), the merged runtime records and
    run log, the record of re-run shards, and SHA256SUMS of every file."""
    metas = [json.loads(p.read_text()) for p in sorted(shards_dir.glob("*/shard.json"))]
    if not metas:
        raise SystemExit(f"no shards under {shards_dir}")
    for field in ("key", "pilot_sha256", "shards", "backgrounds", "beacon"):
        if len({json.dumps(m[field], sort_keys=True) for m in metas}) != 1:
            raise SystemExit(f"shards disagree on {field}")
    got = sorted(m["shard"] for m in metas)
    if got != list(range(metas[0]["shards"])):
        raise SystemExit(f"shards present {got}, expected 0..{metas[0]['shards'] - 1}")
    (out / "reports").mkdir(parents=True, exist_ok=True)
    datasets, runtime, runlog, seen = [], [], [], set()
    for m in metas:
        d = shards_dir / f"shard-{m['shard']}"
        for row in json.loads((d / "manifest.json").read_text())["datasets"]:
            if row["id"] in seen:
                raise SystemExit(f"dataset {row['id']} appears twice")
            seen.add(row["id"])
            datasets.append(row)
            for c in row["cards"]:
                src = d / "reports" / f"{c['id']}.json"
                if hashlib.sha256(src.read_bytes()).hexdigest() != c["report_sha256"]:
                    raise SystemExit(f"report {c['id']} differs from its shard manifest")
                shutil.copyfile(src, out / "reports" / src.name)
        rt = json.loads((d / "runtime.json").read_text())
        runtime.append(dict(shard=m["shard"], **rt["summary"]))
        runlog += [json.loads(line) for line in (d / "runlog.jsonl").read_text().splitlines() if line]
    if expected_datasets is not None and len(datasets) != expected_datasets:
        raise SystemExit(f"{len(datasets)} datasets collected, the design has {expected_datasets}")
    datasets.sort(key=lambda r: r["id"])
    m0 = metas[0]
    manifest = dict(key=m0["key"], beacon=m0["beacon"], pilot_sha256=m0["pilot_sha256"],
                    backgrounds=m0["backgrounds"], shards=m0["shards"], datasets=datasets)
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True))
    (out / "runtime.json").write_text(json.dumps(dict(shards=runtime, reruns=reruns or []), indent=1, default=str))
    (out / "runlog.jsonl").write_text("".join(json.dumps(r, sort_keys=True) + "\n" for r in runlog))
    cards = [c for r in datasets for c in r["cards"]]
    summary = dict(datasets=len(datasets), cards=len(cards), errors=sum(c.get("error", False) for c in cards),
                   key=m0["key"], beacon_round=m0["beacon"].get("round"), reruns=len(reruns or []),
                   slowest_shard_seconds=max(r["wall_seconds"] for r in runtime))
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    sums = sorted(f"{P.sha256(p)}  {p.relative_to(out).as_posix()}" for p in out.rglob("*")
                  if p.is_file() and p.name != "SHA256SUMS")
    (out / "SHA256SUMS").write_text("\n".join(sums) + "\n")
    (out / "summary.md").write_text(
        f"### Panel v1 blind run\n\n- datasets {summary['datasets']}, claim cards {summary['cards']}, "
        f"engine errors {summary['errors']}\n- key: drand {m0['beacon'].get('chain', '?')} round "
        f"{summary['beacon_round']}, randomness `{summary['key']}`\n- shards re-run after an "
        f"infrastructure failure: {summary['reruns']}\n- slowest shard {summary['slowest_shard_seconds']:.0f} s\n"
        f"- sha256 of SHA256SUMS `{P.sha256(out / 'SHA256SUMS')}`\n")
    return summary


REL_TOL = 1e-9  # numbers of a re-run report agree with the published ones within this (relative)


def same_report(a, b, rel: float = REL_TOL) -> bool:
    """Whether two reports say the same: equal everywhere, numbers within `rel` of each other
    (relative; absolutely within rel where both are below 1). Byte-identity holds on the same kind
    of machine; between CPU models the last digits of floating-point numbers can differ (the second
    review: relative differences up to 8.4e-13 between GitHub's runners and an AVX-512 Xeon, with
    the same verdicts and causes)."""
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same_report(a[k], b[k], rel) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(same_report(x, y, rel) for x, y in zip(a, b))
    if isinstance(a, bool) or isinstance(b, bool):
        return type(a) is type(b) and a == b  # True is not 1 in a report
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return a == b or abs(a - b) <= rel * max(1.0, abs(a), abs(b))
    return a == b


def verify(results: Path, compact: Path, pilot: dict, out: Path, workers: int, datasets: int | None = None) -> dict:
    """Re-run datasets of a finished run and compare them with the run's: every dataset's sha256
    with the manifest's, and every report with the published one — byte-identical, or saying the
    same (`same_report`: the verdict, the cause and every field equal, numbers within REL_TOL),
    else different. Anyone can check that the published reports are what the frozen code gives
    for the key. The datasets are the first `datasets` of the key's order (all if None); the key,
    the pilot and the prepared backgrounds must be the run's, and `out` must be empty (the runner
    reuses a report it finds there)."""
    manifest = json.loads((results / "manifest.json").read_text())
    key = P.check_key(manifest["key"])
    record = results / "key.json"
    if record.exists() and check_key_record(json.loads(record.read_text())) != key:
        raise SystemExit("key.json and the manifest name different keys")
    if manifest["pilot_sha256"] != pilot_sha256(pilot):
        raise SystemExit("the pilot differs from the run's")
    if Path(out).exists() and any(Path(out).iterdir()):
        raise SystemExit(f"{out} is not empty: verify re-runs into an empty directory")
    bgs, info = load_prepared(compact, None)
    if json.dumps(info, sort_keys=True) != json.dumps(manifest["backgrounds"], sort_keys=True):
        raise SystemExit("the prepared backgrounds differ from the run's")
    want = {r["id"]: r for r in manifest["datasets"]}
    entries = [e for e in P.assign(key, pilot.get("dropped", ()), pilot.get("pool_size")) if e["id"] in want]
    entries = entries[:datasets] if datasets is not None else entries
    R.run(entries, bgs, pilot, out, workers)
    rerun = {r["id"]: r for r in json.loads((out / "manifest.json").read_text())["datasets"]}
    rows, data = [], []
    for e in entries:
        data.append(dict(dataset=e["id"], identical=rerun[e["id"]]["data_sha256"] == want[e["id"]]["data_sha256"]))
        for c in want[e["id"]]["cards"]:
            mine = (out / "reports" / f"{c['id']}.json").read_bytes()
            got = hashlib.sha256(mine).hexdigest()
            published = results / "reports" / f"{c['id']}.json"  # bound to the manifest by its sha256
            bound = published.exists() and hashlib.sha256(published.read_bytes()).hexdigest() == c["report_sha256"]
            same = got == c["report_sha256"] or (bound and same_report(json.loads(published.read_text()),
                                                                       json.loads(mine)))
            rows.append(dict(card=c["id"], published=c["report_sha256"], rerun=got,
                             identical=got == c["report_sha256"], same=bool(same)))
    return dict(key=key, datasets=len(entries), datasets_identical=sum(d["identical"] for d in data),
                cards=len(rows), identical=sum(r["identical"] for r in rows), same=sum(r["same"] for r in rows),
                different=sum(not r["same"] for r in rows), rows=rows, data=data)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("guard")
    g.add_argument("--frozen-tag", required=True)
    pr = sub.add_parser("prepare")
    pr.add_argument("--backgrounds", type=Path)
    pr.add_argument("--data-dir", type=Path)
    pr.add_argument("--out", type=Path, required=True)
    pr.add_argument("--dry-run", action="store_true")
    for name in ("run", "all"):
        r = sub.add_parser(name)
        r.add_argument("--pilot", type=Path)
        r.add_argument("--key-record", type=Path)
        r.add_argument("--run-tag")
        r.add_argument("--workers", type=int, default=max(1, os.cpu_count() or 1))
        r.add_argument("--out", type=Path, required=True)
        r.add_argument("--dry-run", action="store_true")
        r.add_argument("--limit", type=int, help="dry run: at most this many datasets per shard")
        if name == "run":
            r.add_argument("--compact", type=Path, required=True)
            r.add_argument("--expect-compact")
            r.add_argument("--shard", type=int, required=True)
            r.add_argument("--shards", type=int, required=True)
        else:
            r.add_argument("--backgrounds", type=Path)
            r.add_argument("--data-dir", type=Path)
    v = sub.add_parser("verify", help="re-run datasets of a finished run and compare the datasets and reports")
    v.add_argument("--results", type=Path, required=True, help="the run's results (manifest.json, key.json)")
    v.add_argument("--compact", type=Path, required=True, help="the backgrounds prepared by `prepare`")
    v.add_argument("--pilot", type=Path)
    v.add_argument("--dry-run", action="store_true")
    v.add_argument("--datasets", type=int, help="the first N datasets of the key's order (default: all)")
    v.add_argument("--workers", type=int, default=max(1, os.cpu_count() or 1))
    v.add_argument("--out", type=Path, required=True)
    c = sub.add_parser("collect")
    c.add_argument("--shards-dir", type=Path, required=True)
    c.add_argument("--out", type=Path, required=True)
    c.add_argument("--pilot", type=Path, help="checks the number of datasets against the design")
    c.add_argument("--reruns", type=Path, help="JSON list of the shards re-run after an infrastructure failure")
    args = p.parse_args(argv)

    if args.cmd == "guard":
        print(json.dumps(guard(args.frozen_tag), indent=1))
        return
    if args.cmd == "prepare":
        print(json.dumps(prepare(args.backgrounds, args.out, args.data_dir, args.dry_run), indent=1))
        return
    if args.cmd == "collect":
        pilot = json.loads(args.pilot.read_text()) if args.pilot else None
        n = P.n_datasets(pilot.get("dropped", ())) if pilot and not pilot.get("dry_run") else None
        reruns = json.loads(args.reruns.read_text()) if args.reruns and args.reruns.exists() else []
        print(json.dumps(collect(args.shards_dir, args.out, n, reruns), indent=1))
        return
    if args.cmd == "verify":
        if args.pilot is None:
            if not args.dry_run:
                raise SystemExit("verify needs the run's pilot.json (--pilot), or --dry-run for a dry run")
            import simulate
            pilot = simulate.dry_pilot(simulate.dry_backgrounds())
        else:
            pilot = json.loads(args.pilot.read_text())
        res = verify(args.results, args.compact, pilot, args.out, args.workers, args.datasets)
        for r in res["rows"]:
            word = "identical" if r["identical"] else f"the same within {REL_TOL:g}" if r["same"] else "DIFFERENT"
            print(f"{r['card']}  published {r['published']}  re-run {r['rerun']}  {word}")
        print(f"{res['datasets_identical']} of {res['datasets']} datasets identical; {res['identical']} of "
              f"{res['cards']} reports identical, {res['same']} the same (verdicts, causes, numbers within "
              f"{REL_TOL:g}), {res['different']} different (key {res['key']})")
        if res["different"] or res["datasets_identical"] != res["datasets"]:
            raise SystemExit(1)
        return
    if args.limit is not None and not args.dry_run:
        raise SystemExit("--limit is for the dry run only: the blind run takes every dataset")
    key, beacon = read_key(args.dry_run, args.key_record, args.run_tag)
    if args.pilot is None:
        if not args.dry_run:
            raise SystemExit("the blind run needs the oracle's pilot.json (--pilot)")
        import simulate
        pilot = simulate.dry_pilot(simulate.dry_backgrounds())
    else:
        pilot = json.loads(args.pilot.read_text())
        if pilot.get("dry_run") and not args.dry_run:
            raise SystemExit("the blind run needs the oracle's pilot.json, not a dry-run stand-in")
    if args.cmd == "run":
        meta = run_shard(args.compact, pilot, key, beacon, args.shard, args.shards, args.out, args.workers,
                         args.expect_compact, args.limit)
        print(json.dumps(meta["summary"], indent=1))
        return
    # all: prepare, one shard with every dataset, collect - the same code on one machine
    work = args.out.parent / (args.out.name + ".work")
    res = prepare(args.backgrounds, work / "compact", args.data_dir, args.dry_run)
    run_shard(work / "compact", pilot, key, beacon, 0, 1, work / "shards" / "shard-0", args.workers,
              res["sha256"], args.limit)
    n = None if (args.limit or pilot.get("dry_run")) else P.n_datasets(pilot.get("dropped", ()))
    print(json.dumps(collect(work / "shards", args.out, n), indent=1))


if __name__ == "__main__":
    main()
