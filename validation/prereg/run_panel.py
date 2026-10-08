"""Run the frozen engine once on every claim card of the panel (v1.md, section 8, step 4).

    python validation/prereg/run_panel.py --panel DIR --out REPORTS --workers 8

Each card is run through the Python API with its pre-registration; the JSON report goes to
REPORTS/<id>.json with the run time. A card whose report exists is not run again (one attempt
per card); an engine exception is written as {"error": ...} and scored as a verdict outside the
allowed set. Workers write their own run logs (REPORTS/runlog.<worker>.jsonl), merged into
REPORTS/runlog.jsonl at the end; a claim logged twice is reported. Commit REPORTS with the
sha256 of every file before the key is revealed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from functools import partial
from pathlib import Path

import numpy as np

MODULE_FOLD, MODULE_FRAC = 2.0, 0.3


def _seed(card_id: str) -> int:
    return int(hashlib.sha256(card_id.encode()).hexdigest()[:8], 16)


def load(card: dict, panel_dir: Path):
    from metric_autopsy import SimpleData
    import pandas as pd
    z = np.load(panel_dir / card["data"], allow_pickle=False)
    obs = pd.DataFrame({k[4:]: z[k] for k in z.files if k.startswith("obs_")})
    return SimpleData(z["X"].astype(np.float64), obs, [str(g) for g in z["genes"]])


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


def run_card(job) -> dict:
    card_path, panel_dir, out_dir, worker_log = job
    from metric_autopsy import run_autopsy
    card = json.loads(Path(card_path).read_text())
    out = Path(out_dir) / f"{card['id']}.json"
    if out.exists():
        return dict(id=card["id"], skipped=True)
    t0 = time.time()
    try:
        data = load(card, Path(panel_dir))
        metric, kw = make_run_args(card)
        a = run_autopsy(metric, data, log_path=worker_log, **kw)
        rep = a.to_dict()
        rep.update(id=card["id"], elapsed_seconds=time.time() - t0)
    except Exception as exc:  # scored as outside the allowed set
        rep = dict(id=card["id"], error=repr(exc), elapsed_seconds=time.time() - t0)
    out.write_text(json.dumps(rep, allow_nan=False, default=str))
    return dict(id=card["id"], seconds=rep["elapsed_seconds"], error="error" in rep)


def _init_worker():
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"


def run(panel_dir: Path, out_dir: Path, workers: int = 1, limit: int | None = None) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((panel_dir / "manifest.json").read_text())
    cards = [panel_dir / d["card"] for d in manifest["datasets"]][:limit]
    jobs = [(str(c), str(panel_dir), str(out_dir), str(out_dir / f"runlog.{k % max(workers, 1)}.jsonl"))
            for k, c in enumerate(cards)]
    t0 = time.time()
    if workers <= 1:
        results = [run_card(j) for j in jobs]
    else:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(workers, initializer=_init_worker) as pool:
            results = list(pool.imap_unordered(run_card, jobs, chunksize=1))
    wall = time.time() - t0
    merged, seen, twice = [], set(), []
    for f in sorted(out_dir.glob("runlog.*.jsonl")):
        for line in f.read_text().splitlines():
            rec = json.loads(line)
            if rec.get("claim_id") in seen:
                twice.append(rec.get("claim_id"))
            seen.add(rec.get("claim_id"))
            merged.append(rec)
    merged.sort(key=lambda r: r.get("timestamp_utc", ""))
    (out_dir / "runlog.jsonl").write_text("".join(json.dumps(r) + "\n" for r in merged))
    secs = [r["seconds"] for r in results if "seconds" in r]
    summary = dict(cards=len(jobs), run=len(secs), skipped=sum(r.get("skipped", False) for r in results),
                   errors=sum(r.get("error", False) for r in results), workers=workers,
                   wall_seconds=wall, mean_seconds=float(np.mean(secs)) if secs else None,
                   median_seconds=float(np.median(secs)) if secs else None, logged_twice=twice)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--panel", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--limit", type=int)
    args = p.parse_args(argv)
    print(json.dumps(run(Path(args.panel), Path(args.out), args.workers, args.limit), indent=1))


if __name__ == "__main__":
    main()
