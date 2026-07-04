"""Core data container and result types shared across the gates.

The whole engine is *metric-agnostic*: it never inspects your metric, only probes the
data and the metric's response to controlled perturbations of the data. To do that it
needs a minimal, duck-typed view of single-cell data:

    - ``.X``          cells x genes matrix (numpy ndarray; sparse is densified on read)
    - ``.obs``        pandas DataFrame of per-cell metadata, indexed by STRING cell labels
                      (do not assume a positional integer index — use ``.iloc`` for that)
    - ``.var_names``  sequence of gene names, one per column of X; MUST be unique
    - ``data[mask]``  row-subsetting by boolean mask or integer index -> same type

`anndata.AnnData` already satisfies this contract, so real analyses pass an AnnData
directly. `SimpleData` is a dependency-light stand-in (numpy + pandas only) for tests,
examples, and users who have not installed anndata. Its obs carries a *string* index to
match AnnData, so a bring-your-own ``metric(data)`` that mis-uses ``obs.loc`` fails the
same way on both backends instead of passing on SimpleData and breaking on AnnData.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Sequence

import numpy as np
import pandas as pd


class GateStatus(str, Enum):
    """Outcome of a single gate."""

    PASS = "PASS"          # gate cleared
    FAIL = "FAIL"          # gate failed — the claimed signal is (partly) an artifact
    STOP = "STOP"          # analysis is impossible as posed (e.g. groups incomparable)
    JUDGMENT = "JUDGMENT"  # not auto-decidable; the analyst must rule
    SKIP = "SKIP"          # not run (e.g. no second dataset, or too few cells to test)


@dataclass
class GateResult:
    """The verdict of one gate, with the numbers that produced it."""

    gate: int
    name: str
    status: GateStatus
    message: str
    detail: dict[str, Any] = field(default_factory=dict)

    @property
    def blocking(self) -> bool:
        """A blocking result halts the sequence: you may not proceed past it."""
        return self.status in (GateStatus.FAIL, GateStatus.STOP)

    def __str__(self) -> str:
        return f"GATE {self.gate} [{self.status.value}] {self.name}: {self.message}"


class SimpleData:
    """Dependency-light AnnData stand-in: X (cells x genes) + obs + var_names.

    Supports the exact contract the gates rely on, including ``data[mask]`` row
    subsetting. Not a full AnnData — no layers, no var frame, no I/O — deliberately.
    var_names must be unique (name-based gene lookup would otherwise be ambiguous).
    """

    def __init__(self, X: np.ndarray, obs: pd.DataFrame, var_names: Sequence[str]):
        X = np.asarray(X, dtype=float)
        if X.ndim != 2:
            raise ValueError("X must be 2-D (cells x genes)")
        if X.shape[0] != len(obs):
            raise ValueError(f"X has {X.shape[0]} cells but obs has {len(obs)} rows")
        var_names = list(var_names)
        if X.shape[1] != len(var_names):
            raise ValueError(f"X has {X.shape[1]} genes but {len(var_names)} var_names")
        if len(set(var_names)) != len(var_names):
            dupes = sorted({g for g in var_names if var_names.count(g) > 1})
            raise ValueError(f"var_names must be unique; duplicates: {dupes[:10]}")
        self.X = X
        # Mirror AnnData: obs indexed by string cell labels, not a positional RangeIndex.
        obs = obs.copy()
        obs.index = obs.index.astype(str)
        self.obs = obs
        self.var_names = var_names
        self._gene_index = {g: i for i, g in enumerate(self.var_names)}

    # --- shape ---------------------------------------------------------------
    @property
    def n_obs(self) -> int:
        return self.X.shape[0]

    @property
    def n_vars(self) -> int:
        return self.X.shape[1]

    def __len__(self) -> int:
        return self.n_obs

    # --- access --------------------------------------------------------------
    def gene_vector(self, gene: str) -> np.ndarray:
        """Expression of one gene across all cells."""
        if gene not in self._gene_index:
            raise KeyError(f"gene {gene!r} not in var_names")
        return self.X[:, self._gene_index[gene]]

    def __getitem__(self, mask) -> "SimpleData":
        """Row-subset by boolean mask or integer index array -> new SimpleData."""
        idx = np.asarray(mask)
        if idx.dtype == bool:
            if len(idx) != self.n_obs:
                raise ValueError("boolean mask length must equal n_obs")
        return SimpleData(self.X[idx], self.obs.iloc[idx], self.var_names)

    def __repr__(self) -> str:
        return f"SimpleData(n_obs={self.n_obs}, n_vars={self.n_vars})"


def as_dense(x) -> np.ndarray:
    """Return a dense ndarray whether x came from a numpy array or a sparse matrix."""
    if hasattr(x, "toarray"):  # scipy sparse (AnnData.X is often sparse)
        return np.asarray(x.toarray(), dtype=float)
    return np.asarray(x, dtype=float)


def unique_col_index(var_names, gene: str) -> int:
    """Column index of `gene`, raising if the name is missing or duplicated.

    AnnData permits non-unique var_names (it only warns); silently taking the first match
    would compute the metric on an arbitrary column. Fail loudly instead.
    """
    var_names = list(var_names)
    matches = [i for i, g in enumerate(var_names) if g == gene]
    if not matches:
        raise KeyError(f"gene {gene!r} not in var_names")
    if len(matches) > 1:
        raise ValueError(
            f"gene {gene!r} appears in {len(matches)} columns (indices {matches}); "
            "var_names must be unique for name-based lookup"
        )
    return matches[0]


def gene_column(data, gene: str) -> np.ndarray:
    """Expression of one gene across cells, working for SimpleData *and* AnnData."""
    if hasattr(data, "gene_vector"):
        return data.gene_vector(gene)
    # AnnData path — duplicate-safe.
    j = unique_col_index(list(data.var_names), gene)
    return as_dense(data.X[:, j]).ravel()
