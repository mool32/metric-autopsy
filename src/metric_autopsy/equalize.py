"""Estimand-dependent removal of a technical depth/capture difference between two groups.

Matching cells on n_genes (the v0.1 GATE 2) conditions on a variable that biology itself
moves (cell size, cycling, RNA content), so it can delete real effects. These functions
keep every cell and remove only the technical difference, by binomial thinning of the
deeper group's counts:

* ``composition`` estimands (relative expression): thin to equal **sequencing depth**
  (per-cell total counts) by one common ratio within each stratum. Valid because thinning
  preserves expected composition.
* ``content`` estimands (total RNA, number of genes detected): thinning to equal depth
  would delete the signal itself. The only way to separate capture from content is an
  external standard. With spike-ins (e.g. ERCC) the groups are thinned to equal
  **spike-in capture**; without one the comparison is unidentifiable.
"""
from __future__ import annotations

import itertools
from typing import Sequence

import numpy as np
import pandas as pd

from .core import SimpleData, as_dense


def looks_like_counts(X) -> bool:
    """Non-negative integers (checked on a sample of entries for large matrices)."""
    X = as_dense(X) if not isinstance(X, np.ndarray) else X
    flat = X.ravel()
    if flat.size > 200_000:
        flat = flat[np.random.default_rng(0).integers(0, flat.size, 200_000)]
    return bool(np.all(flat >= 0) and np.allclose(flat, np.round(flat)))


def stratum_masks(obs: pd.DataFrame, within: Sequence[str]):
    """(label dict, boolean mask) for every combination of the `within` factors present."""
    within = list(within)
    if not within:
        yield {}, np.ones(len(obs), dtype=bool)
        return
    levels = [sorted(pd.unique(obs[f].dropna())) for f in within]
    for combo in itertools.product(*levels):
        mask = np.ones(len(obs), dtype=bool)
        for f, v in zip(within, combo):
            mask &= np.asarray(obs[f]) == v
        if mask.any():
            yield dict(zip(within, combo)), mask


def spikein_columns(var_names, prefix: str = "ERCC-") -> list[int]:
    return [j for j, g in enumerate(var_names) if str(g).startswith(prefix)]


def thin_to_match(data: SimpleData, group_col: str, groups: tuple, *, totals: np.ndarray,
                  within: Sequence[str] = (), rng=None) -> tuple[SimpleData, dict]:
    """Thin the higher-`totals` group of every stratum by one common ratio, the ratio of the two
    groups' mean `totals`.

    `totals` is per-cell sequencing depth (composition) or per-cell spike-in counts (content).
    Every cell of the deeper group keeps each molecule with probability
    ``p = mean(totals_shallow) / mean(totals_deep)``. When the groups differ by a capture
    factor this is exact: ``Binomial(Poisson(λ), p) = Poisson(pλ)`` for every gene, so the
    thinned group has the other group's depth with the joint distribution of the genes kept.

    Earlier drafts matched each cell to the other group's depth *quantiles*. That makes the
    keep-probability of a cell depend on its own total, which includes the genes the metric
    reads, and distorts their joint distribution: on a pure depth artifact 18% of the raw MI
    difference survived (7 standard errors) and was "detected" in 5 of 20 datasets with 24 mice
    per group (common ratio: 5%, 1 of 20). Shape differences between the groups' depth
    distributions are reported (``ks_after``), not removed.
    """
    rng = rng if rng is not None else np.random.default_rng(0)
    X = np.array(data.X, dtype=float, copy=True)
    totals = np.asarray(totals, dtype=float)
    obs = data.obs.copy()
    labels = np.asarray(obs[group_col])
    keep_p = np.ones(len(obs))
    info = []
    for stratum, smask in stratum_masks(obs, within):
        a = np.where(smask & (labels == groups[0]))[0]
        b = np.where(smask & (labels == groups[1]))[0]
        if len(a) == 0 or len(b) == 0:
            continue
        hi, lo = (a, b) if totals[a].mean() >= totals[b].mean() else (b, a)
        before = float(totals[hi].mean())
        p = float(np.clip(totals[lo].mean() / before, 0.0, 1.0)) if before > 0 else 1.0
        keep_p[hi] = p
        X[hi] = rng.binomial(X[hi].astype(np.int64), p).astype(float)
        after_t = totals[hi] * p
        info.append(dict(stratum=stratum or "pooled", thinned_group=str(labels[hi[0]]),
                         n_thinned=int(len(hi)), mean_before=before, mean_after=float(after_t.mean()),
                         mean_other=float(totals[lo].mean()), depth_ratio=p,
                         ks_after=_ks_distance(after_t, totals[lo])))
    if "total_counts" in obs:
        obs["total_counts"] = np.asarray(obs["total_counts"], dtype=float) * keep_p
    depth_ratio = min((r["depth_ratio"] for r in info), default=1.0)
    return SimpleData(X, obs, data.var_names), dict(strata=info, depth_ratio=depth_ratio)


def _ks_distance(x: np.ndarray, y: np.ndarray) -> float:
    """Two-sample Kolmogorov-Smirnov distance (a shape diagnostic, not a test)."""
    x, y = np.sort(np.asarray(x, float)), np.sort(np.asarray(y, float))
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    grid = np.concatenate([x, y])
    fx = np.searchsorted(x, grid, side="right") / len(x)
    fy = np.searchsorted(y, grid, side="right") / len(y)
    return float(np.max(np.abs(fx - fy)))
