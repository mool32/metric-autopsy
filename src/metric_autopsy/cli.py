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
        --metric mi_3bin --gene-a Smad3 --gene-b Col1a1 \
        --group-col age --groups young old --within sex \
        --pos-pair Actb Gapdh --neg-pair Gene1 Gene2

``--metric`` names a function in ``metric_autopsy.metrics``. Bring-your-own metrics are
supported from Python via ``metric_autopsy.run_autopsy``; the CLI covers the reference set.
"""
from __future__ import annotations

import argparse
import sys
from functools import partial

import numpy as np
import pandas as pd

from .core import SimpleData
from . import metrics
from .report import run_autopsy


def _load_h5ad(path: str):
    try:
        import anndata
    except ImportError:
        sys.exit("anndata is required to read .h5ad — `pip install anndata`, or use --demo")
    return anndata.read_h5ad(path)


def demo_data(seed: int = 0, n: int = 400) -> SimpleData:
    """Synthetic confound: biology identical everywhere, only male-old QC-degraded.

    The reference case for the whole project — `mi_3bin` on this data dies at GATE 0 and
    GATE 1. Pooled, GATE 2 finds no effect (the male-only confound is diluted); re-run GATE 2
    within the male stratum and it STOPs (young/old n_genes do not overlap).
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
    p.add_argument("--pos-pair", nargs=2)
    p.add_argument("--neg-pair", nargs=2)
    p.add_argument("--data2", help="second .h5ad for GATE 6 replication")
    p.add_argument("--no-stop", action="store_true",
                   help="run all gates instead of stopping at the first blocking failure")
    p.add_argument("--resolve-judgment", action="store_true",
                   help="mark judgment gates (4, 7) resolved, allowing a provisional PASS")
    return p


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
    else:
        data = _load_h5ad(args.h5ad)

    metric_fn = getattr(metrics, args.metric, None)
    if metric_fn is None:
        sys.exit(f"unknown metric {args.metric!r}; choose from "
                 f"{[m for m in dir(metrics) if not m.startswith('_')]}")

    if is_pairwise and args.gene_a and args.gene_b:
        metric = partial(metric_fn, gene_a=args.gene_a, gene_b=args.gene_b)
    else:
        metric = metric_fn

    data2 = _load_h5ad(args.data2) if args.data2 else None
    pairwise_controls = bool(is_pairwise and args.pos_pair and args.neg_pair)

    autopsy = run_autopsy(
        metric, data,
        group_col=args.group_col, groups=tuple(args.groups),
        gene_pair=(args.gene_a, args.gene_b) if (is_pairwise and args.gene_a and args.gene_b) else None,
        within=args.within,
        pair_metric=metric_fn if pairwise_controls else None,
        pos_pair=tuple(args.pos_pair) if pairwise_controls else None,
        neg_pair=tuple(args.neg_pair) if pairwise_controls else None,
        data2=data2,
        prereg={"judgment_pending": not args.resolve_judgment},
        stop_on_first_fail=not args.no_stop,
    )
    autopsy.metric_name = args.metric
    print(autopsy.to_markdown())


if __name__ == "__main__":
    main()
