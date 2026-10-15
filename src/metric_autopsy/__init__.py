"""metric-autopsy — red-team a computed metric before you believe it.

Invert the default from "compute -> believe" to "state your commitments -> red-team ->
then believe". The engine is metric-agnostic: you pass a callable ``metric(data) -> float``
and the names of your factorial ``obs`` columns; the gates probe the data and the metric's
response to controlled perturbations, never the metric's internals.

The verdict has four independent fields — metric validity, design adequacy, effect and
replication — decided by ``report.decide`` (one rule for the API, the CLI and MCP).

Quick start
-----------
>>> from functools import partial
>>> from metric_autopsy import run_autopsy, metrics, injected_signal
>>> m = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
>>> autopsy = run_autopsy(
...     m, adata, group_col="age", groups=("young", "old"), within=["sex"],
...     replicate_col="mouse", gene_pair=("Smad3", "Col1a1"),
...     pair_metric=metrics.norm_pearson, pos_pair=("Actb", "Gapdh"), neg_pair=("Gene1", "Gene2"),
...     prereg={"estimand": "composition", "direction": "decrease", "sesoi": 0.1},
... )
>>> print(autopsy.to_markdown()); autopsy.save_json("autopsy.json")

`adata` may be an `anndata.AnnData` or a `metric_autopsy.SimpleData`.
"""
from .core import Assessment, DenseMemoryWarning, GateResult, GateStatus, SimpleData
from .gates import (
    gate0_independence,
    gate1_qc_parity,
    gate2_ngenes_matching,
    gate3_raw_visibility,
    gate4_signal_response,
    gate5_controls,
    gate6_replication,
)
from .effect import estimate_effect
from .report import Autopsy, decide, run_autopsy
from . import injected_signal, metrics, qc

__version__ = "0.3.0.dev0"

__all__ = [
    "SimpleData",
    "GateResult",
    "GateStatus",
    "Assessment",
    "DenseMemoryWarning",
    "gate0_independence",
    "gate1_qc_parity",
    "gate2_ngenes_matching",
    "gate3_raw_visibility",
    "gate4_signal_response",
    "gate5_controls",
    "gate6_replication",
    "estimate_effect",
    "Autopsy",
    "decide",
    "run_autopsy",
    "injected_signal",
    "metrics",
    "qc",
    "__version__",
]
