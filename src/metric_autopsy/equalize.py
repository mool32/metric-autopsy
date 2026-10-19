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

After equal depth, a claim about a pair (a metric of structure) has its genes' detection equalized
too (``equalize_detection``, v0.4 change 17', the owner's decision of 2026-10-19): equal depth leaves
the groups' detection unequal where their libraries differ in complexity (on preprint §4 it widened the
gap in genes per cell from 24% to 29%), and a metric of structure reads detection.
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
        ma, mb = totals[a].mean(), totals[b].mean()
        # the deeper group; on an exact tie the first by label, whatever the wording's order
        first = ma > mb or (ma == mb and str(groups[0]) <= str(groups[1]))
        hi, lo = (a, b) if first else (b, a)
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


def detection_keep_probability(x, target: float) -> float:
    """The keep probability p at which binomial thinning of the counts `x` (one gene, one group's cells of
    one cell set) leaves an expected share `target` of the cells with a count:
    mean(1 - (1 - p)^x) = target, solved by bisection (the share rises with p, from 0 at p = 0 to the cells'
    share now at p = 1). 1 where the share already is at or below the target, 0 where the target is 0."""
    x = np.asarray(x, dtype=np.int64)
    n, pos = len(x), x[x > 0]
    if n == 0 or len(pos) == 0 or target >= len(pos) / n:
        return 1.0
    if target <= 0:
        return 0.0
    vals, counts = np.unique(pos, return_counts=True)
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if float((counts * (1.0 - (1.0 - mid) ** vals)).sum()) / n < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def equalize_detection(data: SimpleData, genes: Sequence[str], group_col: str, groups: tuple, *,
                       within: Sequence[str] = (), rng=None) -> list:
    """Equalize the claim genes' detection between the groups, after equal depth (v0.4 change 17', the
    owner's decision of 2026-10-19), in place on `data` (a copy the thinning made; its X dense).

    In every cell set of the `within` columns (the strata, and the cell type where the effect is estimated
    per type) and for each gene of `genes`, the group whose cells detect the gene in the larger share has
    the gene's counts thinned binomially, with one keep probability for the set
    (``detection_keep_probability``): its expected share of cells with a count becomes the other group's.
    No threshold: where the shares are equal the probability is 1, and near it the thinning is slight. The
    other genes, and the cells' depth otherwise, are left as they are. Returns one record per set and gene:
    both groups' shares before and after, the thinned group and its keep probability."""
    rng = rng if rng is not None else np.random.default_rng(0)
    X = data.X
    if not isinstance(X, np.ndarray):
        raise TypeError("equalize_detection works on the dense copy that thin_to_match makes")
    names = [str(v) for v in data.var_names]
    cols = [(str(g), names.index(str(g))) for g in genes if str(g) in names]
    labels = np.asarray(data.obs[group_col])
    g0, g1 = (str(g) for g in groups)
    out = []
    for stratum, smask in stratum_masks(data.obs, within):
        a = np.where(smask & (labels == groups[0]))[0]
        b = np.where(smask & (labels == groups[1]))[0]
        if len(a) == 0 or len(b) == 0:
            continue
        for gene, j in cols:
            da, db = float((X[a, j] > 0).mean()), float((X[b, j] > 0).mean())
            hi, target, thinned = (a, db, g0) if da > db else (b, da, g1)
            p = detection_keep_probability(X[hi, j], target) if da != db else 1.0
            if p < 1.0:
                X[hi, j] = rng.binomial(X[hi, j].astype(np.int64), p).astype(float)
            out.append(dict(cells=stratum or "pooled", gene=gene, before={g0: da, g1: db},
                            after={g0: float((X[a, j] > 0).mean()), g1: float((X[b, j] > 0).mean())},
                            thinned_group=thinned if p < 1.0 else None, keep_p=float(p),
                            n_cells={g0: int(len(a)), g1: int(len(b))}))
    return out


def detection_shares(data, genes: Sequence[str], group_col: str, groups: tuple, *, within: Sequence[str] = ()) -> list:
    """Per cell set of the `within` columns and gene of `genes`: each group's share of cells with a count
    (the dossier's uncorrected detection, beside ``equalize_detection``'s)."""
    names = [str(v) for v in data.var_names]
    cols = [(str(g), names.index(str(g))) for g in genes if str(g) in names]
    labels = np.asarray(data.obs[group_col])
    g0, g1 = (str(g) for g in groups)
    columns = {gene: np.asarray(as_dense(data.X[:, [j]]), dtype=float).ravel() for gene, j in cols}
    out = []
    for stratum, smask in stratum_masks(data.obs, within):
        a, b = smask & (labels == groups[0]), smask & (labels == groups[1])
        if not a.any() or not b.any():
            continue
        for gene, _ in cols:
            x = columns[gene]
            out.append(dict(cells=stratum or "pooled", gene=gene,
                            shares={g0: float((x[a] > 0).mean()), g1: float((x[b] > 0).mean())}))
    return out
