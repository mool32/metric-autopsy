#!/usr/bin/env python3
"""Command-line entry point: run the gates and print a Markdown autopsy.

Exposed as the ``metric-autopsy`` console script (pip) and driven by
``scripts/run_gates.py`` (source checkout / the skill).

Examples
--------
Demo on the bundled synthetic confound (no data needed):
    metric-autopsy --demo

Real analysis on your AnnData:
    metric-autopsy --h5ad data.h5ad \
        --metric norm_pearson --gene-a Smad3 --gene-b Col1a1 \
        --group-col age --groups young old --within sex --replicate-col mouse \
        --estimand composition --sesoi 0.1 \
        --pos-pair Actb Gapdh --neg-pair Gene1 Gene2 --json autopsy.json

``--metric`` names a function in ``metric_autopsy.metrics``. Bring-your-own metrics are
supported from Python via ``metric_autopsy.run_autopsy``; the CLI covers the reference set.
A pre-registration can be given as JSON (``--prereg``); explicit flags override it. Every
run with a pre-registration is appended to a JSON-lines run log (``--log``, default
``$METRIC_AUTOPSY_LOG`` or ``metric_autopsy_runs.jsonl``; ``--no-log`` disables) — the same
rule as the Python API. The demo is not a claim and logs only when ``--log`` is given.
"""
from __future__ import annotations

import argparse
import json
import sys
from functools import partial

import numpy as np
import pandas as pd

from .core import SimpleData
from . import injected_signal, metrics
from . import provenance as _prov
from .report import run_autopsy


def _load_h5ad(path: str):
    try:
        import anndata
    except ImportError:
        sys.exit("anndata is required to read .h5ad — `pip install anndata`, or use --demo")
    return anndata.read_h5ad(path)


def demo_data(seed: int = 0, n: int = 400, mice_per_block: int = 4) -> SimpleData:
    """Synthetic confound: biology identical everywhere, only male-old QC-degraded.

    The reference case for the whole project. Male-old cells have 30% capture, so
    `mi_3bin` is lower there, though coupling is identical. With depth equalization the
    difference disappears ("explained by depth"). Each (age, sex) block is split into
    `mice_per_block` mice (`obs["mouse"]`, a replicate column; it does not change X).
    """
    genes = ["Smad3", "Col1a1", "Actb", "Gapdh"] + [f"Gene{i}" for i in range(36)]
    rng = np.random.default_rng(seed)

    def block(coupling, eff, age, sex):
        lat = rng.normal(0, 1, n)
        smad = np.exp(1.2 + coupling * 0.6 * lat + rng.normal(0, 0.3, n))
        col = np.exp(1.2 + coupling * 0.6 * lat + rng.normal(0, 0.3, n))
        hk = rng.normal(0, 1, n)
        actb = np.exp(1.6 + 0.9 * hk + rng.normal(0, 0.2, n))
        gapdh = np.exp(1.6 + 0.9 * hk + rng.normal(0, 0.2, n))
        filler = np.exp(0.4 + rng.normal(0, 0.5, (n, 36)))
        lam = np.column_stack([smad, col, actb, gapdh, filler]) * eff
        return rng.poisson(lam).astype(float), pd.DataFrame({"age": [age] * n, "sex": [sex] * n})

    blocks = [
        block(1.5, 1.0, "young", "male"), block(1.5, 1.0, "young", "female"),
        block(1.5, 0.30, "old", "male"), block(1.5, 1.0, "old", "female"),
    ]
    X = np.vstack([b[0] for b in blocks])
    obs = pd.concat([b[1] for b in blocks], ignore_index=True)
    if mice_per_block:
        # separate RNG so X is identical to the v0.1 demo
        mrng = np.random.default_rng(seed + 1000)
        mouse = np.empty(len(obs), dtype=object)
        for (age, sex), idx in obs.groupby(["age", "sex"], sort=False).indices.items():
            mouse[idx] = [f"{age}-{sex}-m{j}" for j in mrng.permutation(np.arange(len(idx)) % mice_per_block)]
        obs["mouse"] = mouse
    return SimpleData(X, obs, genes)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="metric-autopsy",
                                description="Run the metric-autopsy gates and print a report.")
    p.add_argument("--h5ad", help="path to an AnnData .h5ad")
    p.add_argument("--demo", action="store_true", help="use the bundled synthetic confound")
    p.add_argument("--metric", default="mi_3bin", help="function name in metric_autopsy.metrics")
    p.add_argument("--gene-a"); p.add_argument("--gene-b")
    p.add_argument("--group-col", default="age")
    p.add_argument("--groups", nargs=2, default=("young", "old"))
    p.add_argument("--within", nargs="*", default=["sex"])
    p.add_argument("--replicate-col", help="obs column of the biological replicate (mouse, donor)")
    p.add_argument("--estimand", choices=["composition", "content"],
                   help="what the metric claims to measure: relative expression or RNA content")
    p.add_argument("--direction", choices=["increase", "decrease", "two-sided"],
                   help="the claimed change of the metric from the first group to the second "
                        "(required for SUPPORTED; 'two-sided' marks a non-directional claim)")
    p.add_argument("--sesoi", type=float, help="smallest effect size of interest (construct scale)")
    p.add_argument("--bias-tolerance", type=float,
                   help="GATE 0: a nuisance bias blocks only above this many SESOIs (default 0.5)")
    p.add_argument("--delta-min", type=float,
                   help="GATE 4/5: the smallest response to an injected signal that matters, on the "
                        "metric's scale (default 0.5 x SESOI); a metric is invalid only if shown below it")
    p.add_argument("--min-replicates", type=int, help="minimum replicates per group (>= 3)")
    p.add_argument("--prereg", help="pre-registration JSON (explicit flags override it)")
    p.add_argument("--inject-signal", choices=["coupling"],
                   help="test the metric's response to an injected coupling of the gene pair")
    p.add_argument("--pos-pair", nargs=2)
    p.add_argument("--neg-pair", nargs=2)
    p.add_argument("--data2", help="second .h5ad for GATE 6 replication")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--json", help="write the full JSON report (hashes, parameters) to this path")
    p.add_argument("--log", help=f"append to this run log (default ${{METRIC_AUTOPSY_LOG}} or {_prov.DEFAULT_LOG}; "
                                 "runs with a pre-registration are always logged)")
    p.add_argument("--no-log", action="store_true", help="do not write a run log")
    p.add_argument("--no-stop", action="store_true",
                   help="estimate the effect even when the metric is invalid")
    p.add_argument("--resolve-judgment", action="store_true",
                   help="mark judgment gates (4, 7) resolved, allowing a provisional SUPPORTED")
    return p


def build_prereg(args) -> dict:
    """JSON pre-registration (if any) overridden by explicit flags."""
    prereg = {}
    if args.prereg:
        with open(args.prereg) as fh:
            prereg = json.load(fh)
    for key, val in (("estimand", args.estimand), ("direction", args.direction), ("sesoi", args.sesoi),
                     ("bias_tolerance", args.bias_tolerance), ("delta_min", args.delta_min),
                     ("min_replicates", args.min_replicates)):
        if val is not None:
            prereg[key] = val
    if args.resolve_judgment:
        prereg["judgment_pending"] = False
    return prereg


def main(argv=None):
    args = build_parser().parse_args(argv)

    is_pairwise = args.metric != "spectral_entropy"
    if args.demo or not args.h5ad:
        data = demo_data()
        if is_pairwise:
            args.gene_a = args.gene_a or "Smad3"
            args.gene_b = args.gene_b or "Col1a1"
            args.pos_pair = args.pos_pair or ["Actb", "Gapdh"]
            args.neg_pair = args.neg_pair or ["Gene0", "Gene1"]
        args.replicate_col = args.replicate_col or "mouse"
        args.estimand = args.estimand or "composition"
    else:
        data = _load_h5ad(args.h5ad)

    metric_fn = getattr(metrics, args.metric, None)
    if metric_fn is None:
        sys.exit(f"unknown metric {args.metric!r}; choose from "
                 f"{[m for m in dir(metrics) if not m.startswith('_')]}")

    gene_pair = (args.gene_a, args.gene_b) if (is_pairwise and args.gene_a and args.gene_b) else None
    metric = partial(metric_fn, gene_a=args.gene_a, gene_b=args.gene_b) if gene_pair else metric_fn

    data2 = _load_h5ad(args.data2) if args.data2 else None
    pairwise_controls = bool(is_pairwise and args.pos_pair and args.neg_pair)
    signal = (injected_signal.coupling(*gene_pair) if (args.inject_signal == "coupling" and gene_pair)
              else None)
    is_demo = bool(args.demo or not args.h5ad)
    # None lets run_autopsy apply the shared rule (logged iff pre-registered)
    log_path = "off" if args.no_log else (args.log or ("off" if is_demo else None))

    autopsy = run_autopsy(
        metric, data,
        group_col=args.group_col, groups=tuple(args.groups),
        gene_pair=gene_pair,
        within=args.within,
        pair_metric=metric_fn if pairwise_controls else None,
        pos_pair=tuple(args.pos_pair) if pairwise_controls else None,
        neg_pair=tuple(args.neg_pair) if pairwise_controls else None,
        data2=data2,
        prereg=build_prereg(args),
        stop_on_first_fail=not args.no_stop,
        replicate_col=args.replicate_col,
        signal_test=signal,
        seed=args.seed,
        log_path=log_path,
    )
    autopsy.metric_name = args.metric
    print(autopsy.to_markdown())
    if args.json:
        autopsy.save_json(args.json)


if __name__ == "__main__":
    main()
