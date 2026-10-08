"""Run the frozen engine once on every claim card of the panel, each dataset built on the fly
(validation/prereg/v1.md, section 8). Used by blind.py (the blind run) and timing.py.

Every dataset is built from the key's entry by ``panel.build`` inside a worker process and never
written: the manifest records its canonical sha256 (``panel.dataset_sha256``), its donors, and
per claim card the card's and the report's sha256. A card whose report exists is not run again
(one attempt per card); an engine exception is written as {"error": ...} and scored as a verdict
outside the allowed set. Workers write their own run logs (``runlog.<pid>.jsonl``), merged into
``runlog.jsonl`` at the end; a claim logged twice is reported.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from functools import partial
from pathlib import Path

# One BLAS thread per worker process. OpenBLAS and MKL size their thread pools when numpy loads
# them, so the variables are set before numpy is imported (timing.py and blind.py do the same);
# a value already set wins. machine() records the thread count numpy's BLAS actually uses.
BLAS_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
for _var in BLAS_VARS:
    os.environ.setdefault(_var, "1")

import numpy as np  # noqa: E402

import panel as P  # noqa: E402

MODULE_FOLD, MODULE_FRAC = 2.0, 0.3


def _seed(card_id: str) -> int:
    return int(hashlib.sha256(card_id.encode()).hexdigest()[:8], 16)


def make_run_args(card: dict):
    """The metric and the run_autopsy arguments a claim card specifies."""
    from metric_autopsy import injected_signal, metrics
    from metric_autopsy.core import as_dense
    ga, gb = card["gene_pair"]
    kw = dict(group_col=card["group_col"], groups=tuple(card["groups"]), replicate_col=card["replicate_col"],
              prereg=card["prereg"], seed=_seed(card["id"]))
    st = card["signal_test"]
    kw["signal_test"] = (injected_signal.coupling(*st["genes"]) if st["kind"] == "coupling"
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


def machine() -> dict:
    """Cores and BLAS threads of the run, for the record (v1.md, section 4). The threads in use
    are read with threadpoolctl when it is installed (None otherwise)."""
    try:
        from threadpoolctl import threadpool_info
        used = sorted({i["num_threads"] for i in threadpool_info() if i.get("user_api") == "blas"})
    except ImportError:
        used = None
    return dict(cpus=os.cpu_count(), blas_threads_in_use=used,
                blas_env={v: os.environ.get(v) for v in BLAS_VARS})


_STATE: dict = {}  # backgrounds, pilot and output directory, inherited by forked workers


def run_entry(entry: dict) -> dict:
    """Build one dataset, run the engine on each of its claim cards, return its manifest row."""
    from metric_autopsy import SimpleData, run_autopsy
    bgs, pilot, out_dir = _STATE["bgs"], _STATE["pilot"], _STATE["out"]
    X, obs, genes, cards = P.build(entry, bgs, pilot)
    row = dict(id=entry["id"], data_sha256=P.dataset_sha256(X, obs, genes),
               donors=sorted(set(map(str, obs["donor"]))), cards=[])
    for card in cards:
        out = out_dir / "reports" / f"{card['id']}.json"
        rec = dict(id=card["id"], card_sha256=P.card_sha256(card))
        if out.exists():
            rec.update(skipped=True, report_sha256=hashlib.sha256(out.read_bytes()).hexdigest())
            row["cards"].append(rec)
            continue
        t0 = time.time()
        try:
            metric, kw = make_run_args(card)
            a = run_autopsy(metric, SimpleData(X.copy(), obs.copy(), list(genes)),
                            log_path=str(out_dir / f"runlog.{os.getpid()}.jsonl"), **kw)
            rep = a.to_dict()
            rep.update(id=card["id"], elapsed_seconds=time.time() - t0)
        except Exception as exc:  # scored as outside the allowed set
            rep = dict(id=card["id"], error=repr(exc), elapsed_seconds=time.time() - t0)
        text = json.dumps(rep, allow_nan=False, default=str)
        out.write_text(text)
        rec.update(report_sha256=hashlib.sha256(text.encode()).hexdigest(), seconds=rep["elapsed_seconds"],
                   error="error" in rep)
        row["cards"].append(rec)
    return row


def run(entries: list[dict], bgs: dict, pilot: dict, out_dir: Path, workers: int = 1) -> dict:
    """Run every entry (in parallel with `workers` forked processes); write manifest.json,
    runlog.jsonl and summary.json to out_dir and return the summary."""
    out_dir = Path(out_dir)
    (out_dir / "reports").mkdir(parents=True, exist_ok=True)
    _STATE.update(bgs=bgs, pilot=pilot, out=out_dir)
    t0 = time.time()
    if workers <= 1:
        rows = [run_entry(e) for e in entries]
    else:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(workers) as pool:
            rows = list(pool.imap_unordered(run_entry, entries, chunksize=1))
    wall = time.time() - t0
    rows.sort(key=lambda r: r["id"])
    merged, seen, twice = [], set(), []
    for f in sorted(out_dir.glob("runlog.*.jsonl")):
        for line in f.read_text().splitlines():
            rec = json.loads(line)
            if rec.get("claim_id") in seen:
                twice.append(rec.get("claim_id"))
            seen.add(rec.get("claim_id"))
            merged.append(rec)
        f.unlink()
    merged.sort(key=lambda r: r.get("timestamp_utc", ""))
    with open(out_dir / "runlog.jsonl", "a") as fh:
        fh.write("".join(json.dumps(r) + "\n" for r in merged))
    (out_dir / "manifest.json").write_text(json.dumps(dict(datasets=rows), indent=1))
    cards = [c for r in rows for c in r["cards"]]
    secs = [c["seconds"] for c in cards if "seconds" in c]
    summary = dict(datasets=len(rows), cards=len(cards), run=len(secs),
                   skipped=sum(c.get("skipped", False) for c in cards),
                   errors=sum(c.get("error", False) for c in cards), workers=workers,
                   machine=machine(), wall_seconds=wall,
                   mean_seconds=float(np.mean(secs)) if secs else None,
                   median_seconds=float(np.median(secs)) if secs else None, logged_twice=twice)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=1))
    return summary
