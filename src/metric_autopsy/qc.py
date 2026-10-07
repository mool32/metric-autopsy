"""Per-cell QC and factorial QC comparison — the metric-agnostic core of GATES 1 & 2.

These functions never touch your metric. They ask one question: do the groups you intend
to compare have equivalent technical quality, in every stratification you will use? That
question sank the MI-coupling analysis (male-old cells detected 2.4x fewer genes), and it
is the cheapest, highest-yield check in the whole harness.
"""
from __future__ import annotations

import itertools
from typing import Sequence

import numpy as np
import pandas as pd

from .core import as_dense, warn_if_dense_too_large
from .stats import norm_sf

MITO_PREFIXES = ("mt-", "MT-", "Mt-")
RIBO_PREFIXES = ("Rps", "Rpl", "RPS", "RPL")


def per_cell_qc(data) -> pd.DataFrame:
    """Compute per-cell QC metrics, preferring obs columns if already present.

    Returns a DataFrame (RangeIndex) aligned to cells with columns: n_genes, total_counts,
    and (if identifiable from var_names) pct_mito, pct_ribo.

    Footgun to know about: if ``obs`` already carries ``n_genes_by_counts`` / ``total_counts``
    (e.g. from scanpy's ``calculate_qc_metrics``) they are trusted verbatim, so GATES 1 & 2
    reflect whatever produced them. If ``X`` was normalized or subset *after* those columns
    were written they may be stale — drop the columns to force a recompute from ``X`` if unsure.
    """
    obs = data.obs
    n = len(obs)
    out = pd.DataFrame(index=range(n))

    X = None

    def _dense_X():
        warn_if_dense_too_large(getattr(data.X, "shape", (n, len(data.var_names))), "per-cell QC")
        return as_dense(data.X)

    if "n_genes_by_counts" in obs:
        out["n_genes"] = np.asarray(obs["n_genes_by_counts"])
    else:
        X = _dense_X() if X is None else X
        out["n_genes"] = (X > 0).sum(axis=1)

    if "total_counts" in obs:
        out["total_counts"] = np.asarray(obs["total_counts"])
    else:
        X = _dense_X() if X is None else X
        out["total_counts"] = X.sum(axis=1)

    var_names = list(data.var_names)
    # Sum over ALL columns matching a prefix (duplicate-safe, no name lookup).
    mito_cols = [i for i, g in enumerate(var_names) if str(g).startswith(MITO_PREFIXES)]
    ribo_cols = [i for i, g in enumerate(var_names) if str(g).startswith(RIBO_PREFIXES)]
    if mito_cols or ribo_cols:
        X = _dense_X() if X is None else X
        tot = X.sum(axis=1)
        tot_safe = np.where(tot == 0, 1.0, tot)
        if mito_cols:
            out["pct_mito"] = X[:, mito_cols].sum(axis=1) / tot_safe
        if ribo_cols:
            out["pct_ribo"] = X[:, ribo_cols].sum(axis=1) / tot_safe
    return out


def qc_by_strata(data, factors: Sequence[str]) -> pd.DataFrame:
    """Median QC per factorial combination of `factors` (obs columns)."""
    factors = list(factors)
    missing = [f for f in factors if f not in data.obs.columns]
    if missing:
        raise KeyError(f"factors not in obs: {missing}")
    qc = per_cell_qc(data)
    tbl = qc.copy()
    for f in factors:
        tbl[f] = np.asarray(data.obs[f])
    agg = (
        tbl.groupby(factors, dropna=False)
        .agg(
            n_cells=("n_genes", "size"),
            median_n_genes=("n_genes", "median"),
            median_total_counts=("total_counts", "median"),
        )
        .reset_index()
    )
    return agg


def _dist_overlap(a: np.ndarray, b: np.ndarray) -> float:
    """Matchability of two n_genes distributions: 0 = disjoint ranges, 1 = fully shared.

    The fraction of the *narrower* group's 10-90 percentile range that the two groups
    share. This answers GATE 1's real question — is there a common n_genes band with cells
    from both groups (so they can be n_genes-matched)? A tight distribution nested inside a
    broad one with the same location scores ~1.0 (fully matchable); genuinely separated
    distributions (the male young/old case) score 0.0. Residual within-band shape
    differences are left to GATE 2's balance check.
    """
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return 0.0
    a_lo, a_hi = np.percentile(a, 10), np.percentile(a, 90)
    b_lo, b_hi = np.percentile(b, 10), np.percentile(b, 90)
    lo, hi = max(a_lo, b_lo), min(a_hi, b_hi)
    if hi <= lo:
        # Degenerate: if the pooled range is also a single point the groups coincide (1.0);
        # otherwise the 10-90 ranges are genuinely disjoint (0.0).
        full_lo, full_hi = min(a_lo, b_lo), max(a_hi, b_hi)
        return 1.0 if full_hi <= full_lo else 0.0
    smaller = min(a_hi - a_lo, b_hi - b_lo)
    if smaller <= 0:
        return 1.0  # a zero-width (peaked) group sitting inside the other's range
    return float(min(1.0, (hi - lo) / smaller))


def qc_parity(
    data,
    group_col: str,
    groups: tuple,
    within: Sequence[str] = (),
    thresh: float = 1.5,
    overlap_min: float = 0.2,
    alpha: float = 0.05,
    min_cells: int = 10,
    n_boot: int = 200,
    seed: int = 0,
) -> dict:
    """GATE 1 core: compare two `groups` on QC within every combination of `within`.

    For each stratum, compute the ratio of median n_genes between the groups and the
    matchability overlap of their n_genes ranges. A stratum is flagged only when the breach
    is **confident**: the (1 - alpha/K) interval of the log median ratio lies beyond
    ±log(`thresh`), or the upper bootstrap bound of the overlap is below `overlap_min`.
    K is the number of assessable strata (Bonferroni). Strata where either group has fewer
    than `min_cells` cells are reported but not assessed. In v0.1 every small stratum could
    flag by chance, so stratifying by more factors guaranteed a failure on null data.
    Each flag records which criterion tripped.
    """
    within = list(within)
    obs = data.obs
    qc = per_cell_qc(data)
    ng = qc["n_genes"].values.astype(float)
    g_a, g_b = groups
    rng = np.random.default_rng(seed)
    rows = []

    if within:
        levels = [sorted(pd.unique(obs[f].dropna())) for f in within]
        combos = list(itertools.product(*levels))
    else:
        combos = [()]

    for combo in combos:
        stratum_mask = np.ones(len(obs), dtype=bool)
        for f, v in zip(within, combo):
            stratum_mask &= (np.asarray(obs[f]) == v)
        ga = stratum_mask & (np.asarray(obs[group_col]) == g_a)
        gb = stratum_mask & (np.asarray(obs[group_col]) == g_b)
        na, nb = int(ga.sum()), int(gb.sum())
        if na == 0 or nb == 0:
            continue
        a, b = ng[ga], ng[gb]
        med_a, med_b = float(np.median(a)), float(np.median(b))
        ratio = max(med_a, med_b) / max(min(med_a, med_b), 1e-9)
        overlap = _dist_overlap(a, b)
        assessable = na >= min_cells and nb >= min_cells
        row = dict(
            stratum=dict(zip(within, combo)) if within else "pooled",
            n_a=na, n_b=nb, median_n_genes_a=med_a, median_n_genes_b=med_b,
            n_genes_ratio=ratio, overlap=overlap, assessable=bool(assessable),
            ratio_breach=False, overlap_breach=False, flagged=False,
        )
        if assessable:
            ia = rng.integers(0, na, size=(n_boot, na))
            ib = rng.integers(0, nb, size=(n_boot, nb))
            ba, bb = a[ia], b[ib]
            log_r = np.log(np.maximum(np.median(ba, axis=1), 1e-9)) - np.log(np.maximum(np.median(bb, axis=1), 1e-9))
            ov = np.array([_dist_overlap(x, y) for x, y in zip(ba, bb)])
            row["_boot"] = (log_r, ov)
        rows.append(row)

    assessable_rows = [r for r in rows if r["assessable"]]
    k = max(len(assessable_rows), 1)
    alpha_s = alpha / k
    z = _norm_ppf_upper(alpha_s / 2.0)
    for r in assessable_rows:
        log_r, ov = r.pop("_boot")
        est = float(np.log(max(r["median_n_genes_a"], 1e-9)) - np.log(max(r["median_n_genes_b"], 1e-9)))
        se = float(np.std(log_r, ddof=1)) if len(log_r) > 1 else 0.0
        lo, hi = est - z * se, est + z * se
        r["log_ratio_ci"] = (lo, hi)
        r["ratio_ci"] = (float(np.exp(lo)), float(np.exp(hi)))
        r["ratio_breach"] = bool(lo > np.log(thresh) or hi < -np.log(thresh))
        ov_upper = float(np.quantile(ov, min(1.0, 1.0 - alpha_s / 2.0)))
        r["overlap_upper"] = ov_upper
        r["overlap_breach"] = bool(ov_upper < overlap_min)
        r["flagged"] = bool(r["ratio_breach"] or r["overlap_breach"])

    worst = max((r["n_genes_ratio"] for r in rows), default=float("nan"))
    flagged = [r for r in rows if r["flagged"]]
    return dict(
        table=rows,
        n_compared=len(rows),
        n_assessable=len(assessable_rows),
        worst_ratio=worst,
        flagged=flagged,
        thresh=thresh,
        overlap_min=overlap_min,
        alpha=alpha,
        alpha_per_stratum=alpha_s,
        min_cells=min_cells,
    )


def _norm_ppf_upper(p: float) -> float:
    """z such that P(Z > z) = p (bisection on the normal survival function)."""
    lo, hi = 0.0, 40.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if norm_sf(mid) > p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def match_by_ngenes(
    data,
    group_col: str,
    groups: tuple,
) -> dict:
    """GATE 2 core: restrict both groups to their overlapping n_genes range.

    Returns the retained-cell mask, the overlap range, per-group retained counts, and a KS
    p-value on the matched n_genes distributions (higher = better balanced). If the groups
    do not overlap, `matchable` is False and the analysis should STOP.
    """
    obs = data.obs
    qc = per_cell_qc(data)
    ng = qc["n_genes"].values
    g_a, g_b = groups
    a = np.asarray(obs[group_col]) == g_a
    b = np.asarray(obs[group_col]) == g_b
    if a.sum() == 0 or b.sum() == 0:
        return dict(matchable=False, reason="a group is empty", mask=None)

    lo = max(np.percentile(ng[a], 10), np.percentile(ng[b], 10))
    hi = min(np.percentile(ng[a], 90), np.percentile(ng[b], 90))
    if hi < lo:  # strictly disjoint ranges; hi == lo (peaked, coincident) is matchable
        return dict(matchable=False, reason="n_genes ranges do not overlap",
                    overlap_range=(float(lo), float(hi)), mask=None)

    in_range = (ng >= lo) & (ng <= hi)
    mask = in_range & (a | b)
    ks_p = _ks_2samp(ng[mask & a], ng[mask & b])
    balanced = bool(np.isnan(ks_p) or ks_p > 0.05)
    return dict(
        matchable=True,
        overlap_range=(float(lo), float(hi)),
        mask=mask,
        n_a=int((mask & a).sum()),
        n_b=int((mask & b).sum()),
        ks_p=ks_p,
        balanced=balanced,
    )


def _ks_2samp(a: np.ndarray, b: np.ndarray) -> float:
    """Two-sample KS p-value; uses scipy if available, else an asymptotic fallback."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    try:
        from scipy.stats import ks_2samp  # type: ignore

        return float(ks_2samp(a, b).pvalue)
    except Exception:
        allv = np.sort(np.concatenate([a, b]))
        cdf_a = np.searchsorted(np.sort(a), allv, side="right") / len(a)
        cdf_b = np.searchsorted(np.sort(b), allv, side="right") / len(b)
        d = float(np.max(np.abs(cdf_a - cdf_b)))
        if d == 0.0:
            return 1.0
        en = np.sqrt(len(a) * len(b) / (len(a) + len(b)))
        lam = (en + 0.12 + 0.11 / en) * d
        # Kolmogorov survival series; odd term count so truncation ends on a positive term.
        s = sum((-1) ** (k - 1) * np.exp(-2 * (k * lam) ** 2) for k in range(1, 102))
        return float(max(0.0, min(1.0, 2 * s)))
