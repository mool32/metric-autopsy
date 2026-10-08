"""The blind run of the confirmatory panel (validation/prereg/v1.md, section 8).

Run by .github/workflows/panel.yml when the project owner pushes the tag ``panel-v1-run``, or,
as the fallback, by the owner on one machine (``all``). The key is read from the environment
variable KEY_SEED (a repository secret in Actions) and checked against the commitment in the
frozen tag's message; it is never written anywhere. No dataset is stored: each is built on the
fly, and the manifest records its canonical sha256, so anyone can rebuild and check it once the
key is revealed.

    python blind.py guard   --frozen-tag v0.3.0-prereg            # the code is the frozen tag's
    python blind.py prepare --backgrounds backgrounds.json --out compact   # download, verify, plan
    KEY_SEED=... python blind.py run --compact compact --pilot pilot.json \\
        --commitment-tag v0.3.0-prereg --shard 0 --shards 20 --workers 4 --out out
    python blind.py collect --shards-dir shards --out results
    KEY_SEED=... python blind.py all --backgrounds backgrounds.json --pilot pilot.json \\
        --commitment-tag v0.3.0-prereg --workers 8 --out results       # the one-machine fallback

``--dry-run`` replaces the backgrounds by simulated ones and the key by a public test key
(simulate.py); the workflow runs it on pull requests, which never see the secret.
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
import time
import urllib.request
from pathlib import Path

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")  # as run_panel.py: one BLAS thread per worker, set before numpy

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel as P  # noqa: E402
import run_panel as R  # noqa: E402

FROZEN_PATHS = ("src", "pyproject.toml", "validation/prereg/*.py", "validation/prereg/requirements-panel.txt",
                ".github/workflows/panel.yml")
COMMITMENT = re.compile(r"key sha256:\s*([0-9a-f]{64})")


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=HERE, capture_output=True, text=True, check=True).stdout


# --------------------------------------------------------------------------- #
# guard: the run's code is the frozen tag's
# --------------------------------------------------------------------------- #
def guard(frozen_tag: str) -> dict:
    """The checked-out commit descends from the frozen tag and its engine, panel code, pinned
    requirements and workflow are byte-identical to the tag's; only data and documents
    (backgrounds.json, pilot.json, logs) may differ."""
    head = _git("rev-parse", "HEAD").strip()
    tag_commit = _git("rev-list", "-n", "1", frozen_tag).strip()
    if subprocess.run(["git", "merge-base", "--is-ancestor", tag_commit, head], cwd=HERE).returncode != 0:
        raise SystemExit(f"guard: HEAD {head[:12]} does not descend from {frozen_tag}")
    changed = [f for f in _git("diff", "--name-only", tag_commit, head, "--", *FROZEN_PATHS).split() if f]
    if changed:
        raise SystemExit(f"guard: frozen files differ from {frozen_tag}: {changed}")
    return dict(head=head, frozen_tag=frozen_tag, frozen_commit=tag_commit)


def commitment_from_tag(tag: str) -> str:
    m = COMMITMENT.search(_git("tag", "-l", "--format=%(contents)", tag))
    if not m:
        raise SystemExit(f"no 'key sha256: <64 hex>' line in the message of tag {tag}")
    return m.group(1)


def read_key(dry_run: bool, commitment_tag: str | None) -> str:
    if dry_run:
        import simulate
        return simulate.DRY_RUN_KEY
    key = P.check_key(os.environ.get("KEY_SEED", "").strip())
    if commitment_tag is None:
        raise SystemExit("the real run checks the key against the frozen tag's commitment (--commitment-tag)")
    if P.key_commitment(key) != commitment_from_tag(commitment_tag):
        raise SystemExit("KEY_SEED does not match the commitment in the frozen tag's message")
    return key


# --------------------------------------------------------------------------- #
# prepare: download, verify, plan and compact the backgrounds
# --------------------------------------------------------------------------- #
def fetch(spec: dict, data_dir: Path) -> Path:
    """The background's file, downloaded from spec['url'] if absent, checked against spec['sha256']."""
    data_dir.mkdir(parents=True, exist_ok=True)
    path = data_dir / spec["file"]
    if not path.exists():
        tmp = path.with_suffix(path.suffix + ".part")
        for attempt in range(4):
            try:
                with urllib.request.urlopen(spec["url"], timeout=120) as r, open(tmp, "wb") as fh:
                    shutil.copyfileobj(r, fh, length=1 << 22)
                break
            except OSError:
                if attempt == 3:
                    raise
                time.sleep(2 ** (attempt + 1))
        tmp.rename(path)
    got = P.sha256(path)
    if got != spec["sha256"]:
        raise SystemExit(f"{spec['file']}: sha256 {got} differs from backgrounds.json ({spec['sha256']})")
    return path


def prepare(backgrounds: Path | None, out: Path, data_dir: Path, dry_run: bool = False) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    if dry_run:
        import simulate
        bgs = simulate.dry_backgrounds()
    else:
        for s in json.loads(backgrounds.read_text()).values():
            fetch(s, data_dir)
        bgs = P.load_backgrounds(backgrounds, data_dir)
    info = {}
    for name, bg in bgs.items():
        P.save_compact(bg, out / f"{name}.npz")
        info[name] = dict(file=f"{name}.npz", sha256=P.sha256(out / f"{name}.npz"),
                          content_sha256=P.background_sha256(bg), cells=int(bg.X.shape[0]),
                          genes=len(bg.genes), donors=len(bg.donors))
    (out / "compact.json").write_text(json.dumps(info, indent=1))
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


def run_shard(compact: Path, pilot: dict, key: str, shard: int, shards: int, out: Path, workers: int,
              expect: str | None = None, limit: int | None = None, provenance: dict | None = None) -> dict:
    bgs, info = load_prepared(compact, expect)
    for name, sha in (pilot.get("backgrounds") or {}).items():  # the pilot ran on these backgrounds
        if info.get(name, {}).get("content_sha256") != sha:
            raise SystemExit(f"{name}: the prepared background differs from the one the pilot used")
    entries = shard_entries(P.assign(key, pilot.get("dropped", ())), shard, shards)[:limit]
    summary = R.run(entries, bgs, pilot, out, workers)
    meta = dict(shard=shard, shards=shards, key_sha256=P.key_commitment(key), entries=len(entries),
                backgrounds=info, pilot_sha256=hashlib.sha256(json.dumps(pilot, sort_keys=True).encode()).hexdigest(),
                provenance=provenance or {}, summary=summary)
    (out / "shard.json").write_text(json.dumps(meta, indent=1))
    return meta


def collect(shards_dir: Path, out: Path, expected_datasets: int | None = None) -> dict:
    """Merge the shards: every shard ran with the same key, backgrounds and pilot; no dataset or
    card twice; the reports, the merged run log, the manifest and SHA256SUMS of every file."""
    metas = [json.loads(p.read_text()) for p in sorted(shards_dir.glob("*/shard.json"))]
    if not metas:
        raise SystemExit(f"no shards under {shards_dir}")
    for field in ("key_sha256", "pilot_sha256", "shards"):
        if len({json.dumps(m[field], sort_keys=True) for m in metas}) != 1:
            raise SystemExit(f"shards disagree on {field}")
    if len({json.dumps(m["backgrounds"], sort_keys=True) for m in metas}) != 1:
        raise SystemExit("shards disagree on the backgrounds")
    got = sorted(m["shard"] for m in metas)
    if got != list(range(metas[0]["shards"])):
        raise SystemExit(f"shards present {got}, expected 0..{metas[0]['shards'] - 1}")
    (out / "reports").mkdir(parents=True, exist_ok=True)
    datasets, runlog, seen = [], [], set()
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
        runlog += [json.loads(line) for line in (d / "runlog.jsonl").read_text().splitlines() if line]
    if expected_datasets is not None and len(datasets) != expected_datasets:
        raise SystemExit(f"{len(datasets)} datasets collected, the design has {expected_datasets}")
    datasets.sort(key=lambda r: r["id"])
    runlog.sort(key=lambda r: r.get("timestamp_utc", ""))
    m0 = metas[0]
    manifest = dict(key_sha256=m0["key_sha256"], pilot_sha256=m0["pilot_sha256"], backgrounds=m0["backgrounds"],
                    shards=m0["shards"], provenance=m0["provenance"], datasets=datasets)
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1))
    (out / "runlog.jsonl").write_text("".join(json.dumps(r) + "\n" for r in runlog))
    cards = [c for r in datasets for c in r["cards"]]
    secs = [c["seconds"] for c in cards if "seconds" in c]
    summary = dict(datasets=len(datasets), cards=len(cards), errors=sum(c.get("error", False) for c in cards),
                   run_log_records=len(runlog), key_sha256=m0["key_sha256"],
                   mean_seconds_per_card=(sum(secs) / len(secs)) if secs else None,
                   shard_wall_seconds=[m["summary"]["wall_seconds"] for m in metas],
                   machine=m0["summary"]["machine"], workers=m0["summary"]["workers"])
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    sums = sorted(f"{P.sha256(p)}  {p.relative_to(out).as_posix()}" for p in out.rglob("*")
                  if p.is_file() and p.name != "SHA256SUMS")
    (out / "SHA256SUMS").write_text("\n".join(sums) + "\n")
    (out / "summary.md").write_text(
        f"### Panel v1 blind run\n\n- datasets {summary['datasets']}, claim cards {summary['cards']}, "
        f"engine errors {summary['errors']}\n- key sha256 `{summary['key_sha256']}`\n"
        f"- mean {summary['mean_seconds_per_card'] or float('nan'):.1f} s per card, {summary['workers']} workers "
        f"per shard on {summary['machine']['cpus']} cores; slowest shard {max(summary['shard_wall_seconds']):.0f} s\n"
        f"- sha256 of SHA256SUMS `{P.sha256(out / 'SHA256SUMS')}`\n")
    return summary


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("guard")
    g.add_argument("--frozen-tag", required=True)
    pr = sub.add_parser("prepare")
    pr.add_argument("--backgrounds", type=Path)
    pr.add_argument("--data-dir", type=Path, default=Path("data"))
    pr.add_argument("--out", type=Path, required=True)
    pr.add_argument("--dry-run", action="store_true")
    for name in ("run", "all"):
        r = sub.add_parser(name)
        r.add_argument("--pilot", type=Path)
        r.add_argument("--commitment-tag")
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
            r.add_argument("--data-dir", type=Path, default=Path("data"))
    c = sub.add_parser("collect")
    c.add_argument("--shards-dir", type=Path, required=True)
    c.add_argument("--out", type=Path, required=True)
    c.add_argument("--pilot", type=Path, help="checks the number of datasets against the design")
    args = p.parse_args(argv)

    if args.cmd == "guard":
        print(json.dumps(guard(args.frozen_tag), indent=1))
        return
    if args.cmd == "prepare":
        res = prepare(args.backgrounds, args.out, args.data_dir, args.dry_run)
        print(json.dumps(res, indent=1))
        return
    if args.cmd == "collect":
        pilot = json.loads(args.pilot.read_text()) if args.pilot else None
        n = P.n_datasets(pilot.get("dropped", ())) if pilot else None
        print(json.dumps(collect(args.shards_dir, args.out, n), indent=1))
        return
    if args.limit is not None and not args.dry_run:
        raise SystemExit("--limit is for the dry run only: the blind run takes every dataset")
    key = read_key(args.dry_run, args.commitment_tag)
    if args.dry_run:
        import simulate
        pilot = simulate.dry_pilot(simulate.dry_backgrounds()) if args.pilot is None else json.loads(args.pilot.read_text())
    else:
        pilot = json.loads(args.pilot.read_text())
        if pilot.get("dry_run"):
            raise SystemExit("the blind run needs the oracle's pilot.json, not a dry-run stand-in")
    prov = dict(git_head=_git("rev-parse", "HEAD").strip(), workflow_run=os.environ.get("GITHUB_RUN_ID"),
                dry_run=bool(args.dry_run))
    if args.cmd == "run":
        meta = run_shard(args.compact, pilot, key, args.shard, args.shards, args.out, args.workers,
                         args.expect_compact, args.limit, prov)
        print(json.dumps(meta["summary"], indent=1))
        return
    # all: the one-machine fallback - prepare, one shard with every dataset, collect
    work = args.out.parent / (args.out.name + ".work")
    res = prepare(args.backgrounds, work / "compact", args.data_dir, args.dry_run)
    run_shard(work / "compact", pilot, key, 0, 1, work / "shards" / "shard-0", args.workers, res["sha256"],
              args.limit, prov)
    print(json.dumps(collect(work / "shards", args.out, None if args.limit else P.n_datasets(pilot.get("dropped", ()))),
                     indent=1))


if __name__ == "__main__":
    main()
