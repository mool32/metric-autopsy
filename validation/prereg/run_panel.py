"""Run the frozen engine once on every claim card of the panel, each dataset built on the fly
(validation/prereg/v1.md, section 8). Used by blind.py (the blind run) and timing.py.

Every dataset is built from the key's entry by ``panel.build`` inside a worker process and never
written: the manifest records its canonical sha256 (``panel.dataset_sha256``), its donors, and
per claim card the card's and the report's sha256. A card whose report exists is not run again
(one attempt per card); an engine exception is written as {"error": ...} and scored as an
error.

The run is a deterministic function of the key, the backgrounds and the pilot: the engine's seed
comes from the card (``_seed``), one BLAS thread per worker, and a report holds only what the
engine computed, so the same card gives byte-identical reports on any machine, with any number
of workers (``test_the_same_cards_give_byte_identical_reports``). What varies from run to run -
the time, the machine, the worker's peak memory after each card, the engine's run log - goes to
``runtime.json`` and ``runlog.jsonl``.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import time
from functools import partial
from pathlib import Path

# One BLAS thread per worker process and one numerical path (frozen.NUMERIC_ENV). OpenBLAS and
# MKL size their thread pools and pick their kernels when numpy loads them, so the variables are
# set before numpy is imported (timing.py and blind.py do the same); a value already set wins.
# machine() records them and the thread count numpy's BLAS actually uses.
import frozen  # noqa: E402  (standard library only)

frozen.pin_numerics()

import numpy as np  # noqa: E402

import panel as P  # noqa: E402

MODULE_FOLD, MODULE_FRAC = 2.0, 0.3


def _seed(card_id: str) -> int:
    """The engine's seed for a card: the first 32 bits of the sha256 of its ID."""
    return int(hashlib.sha256(card_id.encode()).hexdigest()[:8], 16)


# Fields of a report that depend on when and where it ran, not on the card: they go to
# runtime.json, so that a report is a function of its card alone.
RUNTIME_FIELDS = ("timestamp_utc", "environment", "log")


def deterministic_report(rep: dict) -> dict:
    """The engine's report without the fields that record when and where it ran."""
    rep = dict(rep)
    rep["provenance"] = {k: v for k, v in (rep.get("provenance") or {}).items() if k not in RUNTIME_FIELDS}
    return rep


def make_run_args(card: dict):
    """The metric and the run_autopsy arguments a claim card specifies."""
    from metric_autopsy import injected_signal, metrics
    from metric_autopsy.core import as_dense
    ga, gb = card["gene_pair"]
    kw = dict(group_col=card["group_col"], groups=tuple(card["groups"]), replicate_col=card["replicate_col"],
              prereg=card["prereg"], seed=_seed(card["id"]))
    st = card["signal_test"]
    kw["signal_test"] = (injected_signal.coupling(*st["genes"], strength=float(st["strength"]))
                         if st["kind"] == "coupling"
                         else injected_signal.module(st["genes"], fold=MODULE_FOLD, frac=MODULE_FRAC))
    if card["metric"] == "norm_pearson":
        metric = partial(metrics.norm_pearson, gene_a=ga, gene_b=gb)
        kw.update(gene_pair=(ga, gb), pair_metric=metrics.norm_pearson,
                  pos_pair=tuple(card["pos_pair"]), neg_pair=tuple(card["neg_pair"]))
    elif card["metric"] == "random":
        rng = np.random.default_rng(_seed(card["id"]))
        metric = (lambda data: float(rng.normal()))
    elif card["metric"] == "constant":
        metric = (lambda data: 0.0)
    elif card["metric"] == "score":
        genes = list(card["score_genes"])

        def metric(data):
            X = as_dense(data.X)
            names = list(data.var_names)
            cols = [names.index(g) for g in genes]
            tot = X.sum(axis=1, keepdims=True)
            return float(np.log1p(X[:, cols] / np.where(tot > 0, tot, 1.0) * 1e4).mean())
    else:
        raise ValueError(card["metric"])
    return metric, kw


def _cpu_model() -> str | None:
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or None


def machine() -> dict:
    """Where the run ran, for the record (v1.md, section 3.3; the second review found the versions
    and the CPU missing): the CPU model and cores, the runner's image (GitHub's ImageOS and
    ImageVersion), the numerical environment as set (`numeric_env`) and as loaded (the third
    review): the BLAS libraries with their threads in use and the kernels they chose
    (threadpoolctl, when installed; None otherwise) and numpy's SIMD baseline and the dispatch
    targets it uses (`numpy_simd`); and the versions of Python, numpy, scipy, pandas and the
    engine. Between CPU models the last digits of a report's numbers can differ; its verdict and
    cause do not (blind.verify)."""
    try:
        from threadpoolctl import threadpool_info
        blas = [{k: i.get(k) for k in ("internal_api", "version", "architecture", "num_threads")}
                for i in threadpool_info() if i.get("user_api") == "blas"]
        used = sorted({b["num_threads"] for b in blas})
    except ImportError:
        blas = used = None
    versions = {"python": platform.python_version(), "platform": platform.platform()}
    for mod in ("numpy", "scipy", "pandas", "metric_autopsy"):
        try:
            versions[mod] = __import__(mod).__version__
        except ImportError:
            versions[mod] = None
    try:
        from numpy._core import _multiarray_umath as mu
        simd = dict(baseline=list(mu.__cpu_baseline__),
                    dispatch=[t for t in mu.__cpu_dispatch__ if mu.__cpu_features__.get(t)])
    except (ImportError, AttributeError):
        simd = None
    image = {k: os.environ[k] for k in ("ImageOS", "ImageVersion") if os.environ.get(k)}
    return dict(cpus=os.cpu_count(), cpu_model=_cpu_model(), runner_image=image or None,
                blas_threads_in_use=used, blas=blas, numpy_simd=simd,
                numeric_env={v: os.environ.get(v) for v in frozen.NUMERIC_ENV}, versions=versions)


_STATE: dict = {}  # backgrounds, pilot and output directory, inherited by forked workers


def peak_rss_mb() -> float | None:
    """The process's peak resident memory so far, in MB (None where the platform has no getrusage)."""
    try:
        import resource
    except ImportError:
        return None
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024  # KB on Linux


def run_entry(entry: dict) -> tuple[dict, list[dict]]:
    """Build one dataset, run the engine on each of its claim cards; return its manifest row
    (deterministic) and the cards' runtime records."""
    from metric_autopsy import SimpleData, run_autopsy
    bgs, pilot, out_dir = _STATE["bgs"], _STATE["pilot"], _STATE["out"]
    X, obs, genes, cards = P.build(entry, bgs, pilot)
    row = dict(id=entry["id"], data_sha256=P.dataset_sha256(X, obs, genes),
               donors=sorted(set(map(str, obs["donor"]))), cards=[])
    runtime = []
    for card in cards:
        out = out_dir / "reports" / f"{card['id']}.json"
        rec = dict(id=card["id"], card_sha256=P.card_sha256(card))
        if out.exists():
            rec.update(report_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),
                       error="error" in json.loads(out.read_bytes()))  # the report's own field, as below
            row["cards"].append(rec)
            runtime.append(dict(id=card["id"], skipped=True))
            continue
        t0 = time.time()
        log = out_dir / "runlog" / f"{card['id']}.jsonl"
        try:
            metric, kw = make_run_args(card)
            a = run_autopsy(metric, SimpleData(X.copy(), obs.copy(), list(genes)), log_path=str(log), **kw)
            full = a.to_dict()
            rep = dict(deterministic_report(full), id=card["id"])
            run_rec = {k: full.get("provenance", {}).get(k) for k in RUNTIME_FIELDS}
        except Exception as exc:  # scored as an error
            rep, run_rec = dict(id=card["id"], error=repr(exc)), {}
        text = json.dumps(rep, allow_nan=False, default=str, sort_keys=True)
        out.write_text(text)
        rec.update(report_sha256=hashlib.sha256(text.encode()).hexdigest(), error="error" in rep)
        row["cards"].append(rec)
        runtime.append(dict(id=card["id"], seconds=time.time() - t0, pid=os.getpid(),
                            peak_rss_mb=peak_rss_mb(), **run_rec))
    return row, runtime


def fork_map(fn, jobs: list, workers: int, chunksize: int = 1) -> list:
    """`fn` over `jobs` in `workers` forked processes, in order. A worker the system kills (out of
    memory) breaks the pool at once (BrokenProcessPool), where multiprocessing.Pool waits for it
    until the job's timeout (the third review)."""
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(workers, mp_context=mp.get_context("fork")) as ex:
        return list(ex.map(fn, jobs, chunksize=chunksize))


def run(entries: list[dict], bgs: dict, pilot: dict, out_dir: Path, workers: int = 1) -> dict:
    """Run every entry (in parallel with `workers` forked processes); write manifest.json (the
    datasets and the sha256 of every card and report, deterministic), runtime.json and
    runlog.jsonl (when, where and how long; the engine's run log) to out_dir and return the
    summary."""
    out_dir = Path(out_dir)
    (out_dir / "reports").mkdir(parents=True, exist_ok=True)
    (out_dir / "runlog").mkdir(parents=True, exist_ok=True)
    _STATE.update(bgs=bgs, pilot=pilot, out=out_dir)
    t0 = time.time()
    if workers <= 1:
        results = [run_entry(e) for e in entries]
    else:
        results = fork_map(run_entry, entries, workers)
    wall = time.time() - t0
    rows = sorted((r for r, _ in results), key=lambda r: r["id"])
    runtime = sorted((c for _, rt in results for c in rt), key=lambda c: c["id"])
    merged = []
    for f in sorted((out_dir / "runlog").glob("*.jsonl")):
        merged += [json.loads(line) for line in f.read_text().splitlines() if line]
    twice = sorted({r["claim_id"] for r in merged if r.get("claim_id")
                    and sum(x.get("claim_id") == r["claim_id"] for x in merged) > 1})
    with open(out_dir / "runlog.jsonl", "a") as fh:
        fh.write("".join(json.dumps(r, sort_keys=True) + "\n" for r in merged))
    shutil.rmtree(out_dir / "runlog")
    (out_dir / "manifest.json").write_text(json.dumps(dict(datasets=rows), indent=1, sort_keys=True))
    secs = [c["seconds"] for c in runtime if "seconds" in c]
    summary = dict(datasets=len(rows), cards=sum(len(r["cards"]) for r in rows), run=len(secs),
                   skipped=sum(c.get("skipped", False) for c in runtime),
                   errors=sum(c.get("error", False) for r in rows for c in r["cards"]), workers=workers,
                   machine=machine(), wall_seconds=wall,
                   mean_seconds=float(np.mean(secs)) if secs else None,
                   median_seconds=float(np.median(secs)) if secs else None,
                   peak_rss_mb=max((c.get("peak_rss_mb") or 0 for c in runtime), default=None),
                   logged_twice=twice)
    (out_dir / "runtime.json").write_text(json.dumps(dict(summary=summary, cards=runtime), indent=1, default=str))
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=1))
    return summary
