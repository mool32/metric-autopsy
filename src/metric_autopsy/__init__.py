"""metric-autopsy — red-team a computed metric before you believe it.

Invert the default from "compute -> believe" to "state your commitments -> red-team ->
then believe". The engine is metric-agnostic: you pass a callable ``metric(data) -> float``
and the names of your factorial ``obs`` columns; the gates probe the data and the metric's
response to controlled perturbations, never the metric's internals.

Quick start
-----------
>>> from functools import partial
>>> from metric_autopsy import run_autopsy, metrics
>>> m = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")
>>> autopsy = run_autopsy(
...     m, adata, group_col="age", groups=("young", "old"),
...     gene_pair=("Smad3", "Col1a1"), within=["sex"],
...     pair_metric=metrics.mi_3bin, pos_pair=("Actb", "Gapdh"), neg_pair=("Gene1", "Gene2"),
... )
>>> print(autopsy.to_markdown())

`adata` may be an `anndata.AnnData` or a `metric_autopsy.SimpleData`.
"""
from .core import GateResult, GateStatus, SimpleData
from .gates import (
    gate0_independence,
    gate1_qc_parity,
    gate2_ngenes_matching,
    gate3_raw_visibility,
    gate5_controls,
    gate6_replication,
)
from .report import Autopsy, run_autopsy
from . import metrics, qc

__version__ = "0.1.1"

__all__ = [
    "SimpleData",
    "GateResult",
    "GateStatus",
    "gate0_independence",
    "gate1_qc_parity",
    "gate2_ngenes_matching",
    "gate3_raw_visibility",
    "gate5_controls",
    "gate6_replication",
    "Autopsy",
    "run_autopsy",
    "metrics",
    "qc",
    "__version__",
]
