"""The gates as metric-agnostic functions.

Each ``gateN_*`` takes your metric (or the raw data) as a black box and returns a
``GateResult``. ``report.run_autopsy`` combines them into four independent verdict fields
(metric validity, design adequacy, effect, replication); see ``report.decide``.

Auto gates:  0 (nuisance invariance: bias vs attenuation), 1 (QC parity — a design
             diagnostic, not a blocking gate), 3 (raw visibility), 4 (response to an injected
             signal, when supplied), 5 (controls against an empirical null),
             6 (replication, if a 2nd dataset is supplied).
GATE 2 is now the estimand-dependent correction inside ``effect.estimate_effect``;
``gate2_ngenes_matching`` is kept only for backward compatibility (it conditions on a
variable that biology moves, so ``run_autopsy`` no longer uses it).
Judgment gates 4 (beyond the injected-signal check) and 7 remain the analyst's.
See references/gates.md for the full rationale.
"""
from __future__ import annotations

from typing import Callable, Sequence

import numpy as np
import pandas as pd

from .core import GateResult, GateStatus, SimpleData, as_dense, warn_if_dense_too_large
from . import qc as _qc
from .stats import extend_null

Metric = Callable[[object], float]
# PairMetric: fn(data, *, gene_a, gene_b) -> float  (used by GATES 3 & 5)


def _as_simple(data) -> SimpleData:
    """Materialize any duck-typed data object as a mutable SimpleData copy."""
    if isinstance(data, SimpleData):
        return SimpleData(data.X.copy(), data.obs.copy(), data.var_names)
    warn_if_dense_too_large(getattr(data.X, "shape", (len(data.obs), len(data.var_names))),
                            "the gates (dense working copy)")
    return SimpleData(as_dense(data.X).copy(), data.obs.copy(), list(data.var_names))


def _safe_call(fn):
    """Evaluate a metric, returning None (not raising) on a perturbation it can't handle.

    The gates treat the metric as an opaque black box and must not crash if a user's callable
    throws on a perturbed input (e.g. a LinAlgError on a degenerate matrix, an IndexError after
    gene subsampling). Any exception is therefore swallowed to None; the gate reports SKIP when
    too few evaluations succeed rather than propagating a failure from inside the metric.
    """
    try:
        v = fn()
    except Exception:
        return None
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


# --------------------------------------------------------------------------- #
# GATE 0 — mathematical independence: bias vs attenuation
# --------------------------------------------------------------------------- #
NULL_DEPTH_BINS = 10


def _looks_like_counts(X: np.ndarray) -> bool:
    return bool(np.all(X >= 0) and np.allclose(X, np.round(X)))


def _perturb(data: SimpleData, kind: str, rng: np.random.Generator,
             protect: frozenset = frozenset()) -> tuple[SimpleData, dict]:
    """Return (perturbed data, note-dict) for one nuisance perturbation."""
    X = data.X.copy()
    note = {}
    if kind == "extra_dropout":
        nz = np.argwhere(X > 0)
        k = int(0.2 * len(nz))
        if k:
            pick = nz[rng.choice(len(nz), size=k, replace=False)]
            X[pick[:, 0], pick[:, 1]] = 0.0
    elif kind == "depth_downsample":
        if _looks_like_counts(X) and np.all(X <= np.iinfo(np.int64).max):
            X = rng.binomial(X.astype(np.int64), 0.5).astype(float)
        else:
            # Non-count input: per-CELL depth reduction (not a global scalar), so a
            # per-cell-depth-confounded metric is still probed rather than trivially invariant.
            per_cell = rng.uniform(0.25, 0.75, size=X.shape[0])[:, None]
            X = X * per_cell
            note["note"] = "approximate: input is not raw counts (per-cell scaling used)"
    elif kind == "library_scale":
        factors = rng.uniform(0.5, 2.0, size=X.shape[0])[:, None]
        X = X * factors
    elif kind == "gene_subsample":
        keep = rng.random(X.shape[1]) < 0.8
        for j, gname in enumerate(data.var_names):
            if gname in protect:  # never drop the genes the metric is bound to
                keep[j] = True
        if keep.sum() >= 2:
            return SimpleData(X[:, keep], data.obs,
                              [g for g, k in zip(data.var_names, keep) if k]), note
    else:
        raise ValueError(kind)
    return SimpleData(X, data.obs, data.var_names), note


MIN_CELLS_PER_DEPTH_BIN = 10


def _depth_bins(sd: SimpleData, n_bins: int = NULL_DEPTH_BINS) -> np.ndarray:
    """Depth-quantile bins with at least MIN_CELLS_PER_DEPTH_BIN cells each (fewer bins for
    small data: a bin of one cell cannot be shuffled, which made the null equal the data)."""
    tot = np.asarray(_qc.per_cell_qc(sd)["total_counts"], dtype=float)
    n_bins = max(1, min(n_bins, len(tot) // MIN_CELLS_PER_DEPTH_BIN))
    ranks = np.argsort(np.argsort(tot, kind="mergesort"), kind="mergesort")
    return (ranks * n_bins // max(len(tot), 1)).astype(int)


def _shuffle_null(sd: SimpleData, bins: np.ndarray, rng: np.random.Generator) -> SimpleData:
    """Permute every gene independently across cells *within depth bins*.

    Keeps each gene's distribution and the depth structure; destroys gene-gene co-variation.
    The metric's value on this data is its 'no biology, same technology' reference.
    """
    X = sd.X.copy()
    for b in np.unique(bins):
        idx = np.where(bins == b)[0]
        if len(idx) < 2:
            continue
        order = np.argsort(rng.random((len(idx), X.shape[1])), axis=0)
        X[idx] = np.take_along_axis(X[idx], order, axis=0)
    return SimpleData(X, sd.obs, sd.var_names)


def _mean_sd(vals):
    vals = np.asarray(vals, float)
    return float(vals.mean()), (float(vals.std(ddof=1)) if len(vals) > 1 else 0.0)


def gate0_independence(
    metric: Metric,
    data,
    tol: float = 0.25,
    abs_tol: float = 1e-3,
    z_thresh: float = 4.0,
    n_baseline: int = 60,
    n_perturb: int = 20,
    protect_genes: Sequence[str] = (),
    include_matrix_perturbations: bool = False,
    seed: int = 0,
    effect_scale: float | None = None,
    n_null: int = 10,
) -> GateResult:
    """Does a nuisance *create or inflate* signal (bias), or only *shrink* it (attenuation)?

    Three references are measured on the user's own data: a bootstrap baseline (the metric's
    sampling spread), an automatic null (every gene permuted within depth bins: no gene-gene
    biology, same technology), and each nuisance perturbation. Every shift is classified:

    * **bias** — the nuisance moves the null, reverses the signal, or inflates it. It can
      create a false effect, so GATE 0 FAILs;
    * **attenuation** — the signal shrinks toward the null while the null stays put. This
      is reliability, not confounding: uniform attenuation pulls toward zero and cannot
      create an effect. It is reported (``detail["attenuation"]``, fraction of signal lost)
      and passed to ``design_adequacy`` as a power check against the SESOI;
    * **level_shift** — no gene-gene structure to reference (e.g. a mean score), so
      attenuation and bias cannot be told apart. A shift larger than `tol` x the effect
      scale is flagged and reported, not failed. A between-group version is removed by the
      effect field's equalization.

    Structure metrics are sized in units of their signal above the null. Level metrics — and
    metrics whose signal above the null is smaller than `tol` x the effect scale — are sized
    in units of the effect (``effect_scale``: the observed difference or the SESOI, as passed
    by ``run_autopsy``), else of the null SD. Every rule is invariant to affine
    re-expressions ``a·m + b`` of the metric (v0.1 divided by |baseline|, so adding a
    constant passed a confounded metric). A metric that never varies is DEGENERATE.
    ``abs_tol`` is accepted for backward compatibility and ignored.
    """
    sd = _as_simple(data)
    rng = np.random.default_rng(seed)
    protect = frozenset(protect_genes)
    n = sd.n_obs

    base_vals = []
    for _ in range(n_baseline):
        idx = rng.integers(0, n, size=n)
        v = _safe_call(lambda: metric(sd[idx]))
        if v is not None:
            base_vals.append(v)
    if len(base_vals) < 3:
        return GateResult(
            0, "Mathematical independence", GateStatus.SKIP,
            "metric could not be evaluated on bootstrap resamples", dict(),
        )
    base_mean, base_sd = _mean_sd(base_vals)

    bins = _depth_bins(sd)
    null_vals = [v for v in (_safe_call(lambda: metric(_shuffle_null(sd, bins, rng)))
                             for _ in range(n_null)) if v is not None]
    if len(null_vals) >= 2:
        null_center, null_sd = _mean_sd(null_vals)
        signal = base_mean - null_center
        signal_se = float(np.sqrt(base_sd ** 2 / len(base_vals) + null_sd ** 2 / len(null_vals)))
        has_signal = abs(signal) > max(3.0 * signal_se, 1e-12 * max(1.0, abs(base_mean)))
    else:
        null_center = null_sd = signal = None
        has_signal = False

    # Structure metrics are judged against their own signal above the null (scale-free and
    # independent of the comparison). Level metrics (no structure in the null) are sized
    # against the effect scale (observed difference or SESOI), else the null SD. A signal
    # that is detectable but immaterial next to the effect under test (e.g. a mean log
    # total moved by the shuffle's loss of cell-size covariance) does not make a metric a
    # structure metric: sizing nuisance shifts in its units would fail any level metric.
    has_effect_scale = effect_scale is not None and np.isfinite(effect_scale) and effect_scale > 0
    if has_signal and has_effect_scale and abs(signal) < tol * float(effect_scale):
        has_signal = False
        immaterial_signal = True
    else:
        immaterial_signal = False
    if has_signal:
        scale, scale_source = abs(signal), "signal above null"
    elif has_effect_scale:
        scale, scale_source = float(effect_scale), "effect"
    elif base_sd > 0:
        scale, scale_source = base_sd, "null SD (size not judged)"
    else:
        scale, scale_source = None, "none"

    kinds = ["extra_dropout", "depth_downsample", "library_scale"]
    if include_matrix_perturbations:
        # gene_subsample is a whole-matrix perturbation (changes dimensionality): meaningful
        # only for metrics that read the whole matrix, so run_autopsy auto-enables it for
        # metrics with no bound gene pair.
        kinds += ["gene_subsample"]

    all_vals = list(base_vals) + list(null_vals)
    responses, attenuation, level_shifts = {}, {}, {}
    for kind in kinds:
        vals, note = [], {}
        for _ in range(n_perturb):
            pdata, note = _perturb(sd, kind, rng, protect)
            v = _safe_call(lambda: metric(pdata))
            if v is not None:
                vals.append(v)
        if len(vals) < 3:
            responses[kind] = dict(skipped=True, rel_change=0.0, z=0.0, shift_std=0.0,
                                   classification="none", confounded=False, **note)
            continue
        all_vals += vals
        pert_mean, pert_sd = _mean_sd(vals)
        shift = pert_mean - base_mean
        se = float(np.sqrt(base_sd ** 2 / len(base_vals) + pert_sd ** 2 / len(vals)))
        z = abs(shift) / se if se > 0 else (float("inf") if shift != 0 else 0.0)
        shift_std = shift / scale if scale else float("nan")
        r = dict(mean_perturbed=pert_mean, shift=float(shift), z=float(z),
                 shift_std=float(shift_std), rel_change=float(abs(shift_std)) if scale else 0.0,
                 **note)

        if z <= z_thresh:
            cls = "none"
        elif has_signal:
            nvals = []
            for _ in range(min(5, n_null)):
                pdata, _ = _perturb(sd, kind, rng, protect)
                pbins = bins if pdata.n_obs == len(bins) else _depth_bins(pdata)
                v = _safe_call(lambda: metric(_shuffle_null(pdata, pbins, rng)))
                if v is not None:
                    nvals.append(v)
            if len(nvals) >= 2:
                null_pert, null_pert_sd = _mean_sd(nvals)
                null_shift = null_pert - null_center
                null_se = float(np.sqrt(null_sd ** 2 / len(null_vals) + null_pert_sd ** 2 / len(nvals)))
                null_z = abs(null_shift) / null_se if null_se > 0 else (float("inf") if null_shift else 0.0)
                signal_pert = pert_mean - null_pert
                r.update(null_perturbed=null_pert, null_shift=float(null_shift),
                         null_shift_std=float(null_shift / scale), null_shift_z=float(null_z),
                         signal=float(signal), signal_perturbed=float(signal_pert))
                null_moves = abs(null_shift) / scale > tol and null_z > z_thresh
                reverses = np.sign(signal_pert) != np.sign(signal) and abs(signal_pert) > tol * abs(signal)
                change = (abs(signal_pert) - abs(signal)) / abs(signal)
                r["signal_change"] = float(change)
                if null_moves or reverses or change > tol:
                    cls = "bias"
                    r["bias_kind"] = ("null moves" if null_moves else
                                      "signal reverses" if reverses else "signal inflated")
                elif change < 0:
                    cls = "attenuation"
                    r["signal_loss"] = float(-change)
                    attenuation[kind] = float(-change)
                else:
                    cls = "none"
            else:
                cls = "unclassified"
        elif scale_source == "effect":
            # No gene-gene structure to tell attenuation from bias (e.g. a mean score pulled
            # toward its no-signal level by dropout). Reported as a flag, not failed: a
            # between-group version is removed by the effect field's equalization.
            cls = "level_shift" if abs(shift_std) > tol else "none"
            if cls == "level_shift":
                level_shifts[kind] = float(shift_std)
        else:
            cls = "unscaled"
        r["classification"] = cls
        r["confounded"] = cls == "bias"
        responses[kind] = r

    detail = dict(baseline=base_mean, baseline_sd=base_sd, null_center=null_center,
                  null_sd=null_sd, signal=signal, has_signal=bool(has_signal),
                  immaterial_signal=immaterial_signal, scale=scale,
                  scale_source=scale_source, effect_scale=effect_scale, responses=responses,
                  attenuation=attenuation, level_shifts=level_shifts, tol=tol, z_thresh=z_thresh)

    spread = float(np.ptp(all_vals)) if all_vals else 0.0
    if spread <= 1e-12 * max(1.0, abs(base_mean)):
        detail["degenerate"] = True
        return GateResult(
            0, "Mathematical independence", GateStatus.DEGENERATE,
            f"the metric returned the same value ({base_mean:.6g}) on every resample, null and "
            "perturbation — it carries no information about the data", detail)

    flagged = {k: r for k, r in responses.items() if r.get("confounded")}
    if flagged:
        worst = max(flagged, key=lambda k: abs(flagged[k].get("null_shift_std", flagged[k]["shift_std"])))
        w = flagged[worst]
        what = w.get("bias_kind", "level shift")
        size = abs(w.get("null_shift_std", w["shift_std"]))
        return GateResult(
            0, "Mathematical independence", GateStatus.FAIL,
            f"'{worst}' biases the metric ({what}): {size:.0%} of the {scale_source} scale "
            f"(z={w['z']:.1f}) — the nuisance can create a difference with no biology",
            detail,
        )
    parts = []
    if attenuation:
        parts.append("attenuation " + ", ".join(f"{k} −{v:.0%}" for k, v in attenuation.items())
                     + " of the signal (passed to design adequacy as a power check)")
    if level_shifts:
        parts.append("level shift " + ", ".join(f"{k} {v:+.0%}" for k, v in level_shifts.items())
                     + " of the effect scale (no structure to classify it; the effect field equalizes "
                     "depth between groups)")
    msg = ("no bias; " + "; ".join(parts) + " — reported, not failed") if parts else \
        "no bias or attenuation beyond estimator noise under any nuisance perturbation"
    return GateResult(0, "Mathematical independence", GateStatus.PASS, msg, detail)


# --------------------------------------------------------------------------- #
# GATE 1 — QC parity across factorial strata (a design diagnostic)
# --------------------------------------------------------------------------- #
def gate1_qc_parity(
    data,
    group_col: str,
    groups: tuple,
    within: Sequence[str] = (),
    thresh: float = 1.5,
    alpha: float = 0.05,
    min_cells: int = 10,
    seed: int = 0,
) -> GateResult:
    """Do the compared groups have equivalent QC in *every* stratum of `within`?

    A design diagnostic, not a blocking gate. A confident imbalance (see ``qc.qc_parity``:
    bootstrap intervals, Bonferroni across strata, ``min_cells``) returns WARN. The effect
    field then removes it in an estimand-appropriate way, or declares the comparison
    unidentifiable. STOP means no stratum contains both groups.
    """
    res = _qc.qc_parity(data, group_col, groups, within=within, thresh=thresh,
                        alpha=alpha, min_cells=min_cells, seed=seed)
    if res["n_compared"] == 0:
        return GateResult(
            1, "QC parity", GateStatus.STOP,
            f"no stratum of {list(within)} contains both groups {groups} — groups are "
            "confounded with the stratifier and cannot be compared",
            res,
        )
    if res["n_assessable"] == 0:
        return GateResult(
            1, "QC parity", GateStatus.SKIP,
            f"no stratum has >= {min_cells} cells in both groups — QC parity not assessable",
            res,
        )
    flagged = res["flagged"]
    if flagged:
        ratio_breaches = [r for r in flagged if r["ratio_breach"]]
        overlap_breaches = [r for r in flagged if r["overlap_breach"]]
        parts = []
        if ratio_breaches:
            wr = max(ratio_breaches, key=lambda r: r["n_genes_ratio"])
            parts.append(f"{len(ratio_breaches)} confidently beyond {thresh}x QC ratio "
                         f"(worst {wr['n_genes_ratio']:.2f}x at {wr['stratum']})")
        if overlap_breaches:
            wo = min(overlap_breaches, key=lambda r: r["overlap"])
            parts.append(f"{len(overlap_breaches)} with n_genes overlap confidently < {res['overlap_min']} "
                         f"(worst {wo['overlap']:.2f} at {wo['stratum']})")
        return GateResult(
            1, "QC parity", GateStatus.WARN,
            f"{len(flagged)}/{res['n_assessable']} assessable strata are QC-imbalanced: "
            + "; ".join(parts) + " — handled by the estimand-dependent correction (effect field)",
            res,
        )
    return GateResult(
        1, "QC parity", GateStatus.PASS,
        f"no confident QC imbalance in {res['n_assessable']} assessable strata "
        f"(worst point ratio {res['worst_ratio']:.2f}x; alpha {alpha} Bonferroni)",
        res,
    )


# --------------------------------------------------------------------------- #
# GATE 2 (deprecated) — n_genes matching
# --------------------------------------------------------------------------- #
def _group_effect(metric: Metric, data, group_col: str, groups: tuple) -> float:
    """Effect = metric(group A) - metric(group B)."""
    a = np.asarray(data.obs[group_col]) == groups[0]
    b = np.asarray(data.obs[group_col]) == groups[1]
    return float(metric(data[a]) - metric(data[b]))


def _permutation_effect_floor(metric, data, group_col, groups, n_perm, rng) -> float:
    """95th-percentile |effect| under shuffled group labels — a metric-scale-free floor
    below which an 'effect' is indistinguishable from label noise."""
    obs = data.obs
    a = np.asarray(obs[group_col]) == groups[0]
    b = np.asarray(obs[group_col]) == groups[1]
    union = np.where(a | b)[0]
    labels = np.where(a[union], 0, 1)
    vals = []
    for _ in range(n_perm):
        perm = rng.permutation(labels)
        ga = union[perm == 0]
        gb = union[perm == 1]
        if len(ga) == 0 or len(gb) == 0:
            continue
        v = _safe_call(lambda: metric(data[ga]) - metric(data[gb]))
        if v is not None:
            vals.append(abs(v))
    if not vals:
        return 0.0
    return float(np.percentile(vals, 95))


def _matched_ratio(data, group_col, groups, mask) -> float:
    """Median n_genes ratio between the two groups within the matched subset."""
    sub = data[mask]
    ng = _qc.per_cell_qc(sub)["n_genes"].values
    a = np.asarray(sub.obs[group_col]) == groups[0]
    b = np.asarray(sub.obs[group_col]) == groups[1]
    if a.sum() == 0 or b.sum() == 0:
        return float("inf")
    ma, mb = np.median(ng[a]), np.median(ng[b])
    return max(ma, mb) / max(min(ma, mb), 1e-9)


def gate2_ngenes_matching(
    metric: Metric,
    data,
    group_col: str,
    groups: tuple,
    keep_frac: float = 0.5,
    max_frac: float = 3.0,
    min_effect: float | None = None,
    balance_ratio: float = 1.1,
    n_perm: int = 20,
    seed: int = 0,
) -> GateResult:
    """DEPRECATED (v0.1 GATE 2): does the effect survive restricting cells to a shared n_genes range?

    run_autopsy no longer calls this. n_genes is downstream of biology (cell size,
    cycling, RNA content), so selecting cells on it is a bad control that can delete real
    effects. Use the estimand-dependent correction in effect.estimate_effect instead.
    Kept for backward compatibility.

    Does the between-group effect survive equalizing cell quality (n_genes)?

    Guards: (a) if the unmatched effect is below a permutation-null floor there is nothing
    to preserve; (b) the matched effect must keep its sign and stay within
    [keep_frac, max_frac] of the unmatched one — collapse *and* runaway amplification fail;
    (c) the matched subset must actually be balanced (median n_genes within `balance_ratio`).
    """
    rng = np.random.default_rng(seed)
    m = _qc.match_by_ngenes(data, group_col, groups)
    if not m["matchable"]:
        return GateResult(
            2, "n_genes matching", GateStatus.STOP,
            f"groups are incomparable: {m['reason']} — do not proceed",
            m,
        )
    unmatched = _group_effect(metric, data, group_col, groups)
    matched = _group_effect(metric, data[m["mask"]], group_col, groups)
    retained = abs(matched) / (abs(unmatched) + 1e-12)
    matched_ratio = _matched_ratio(data, group_col, groups, m["mask"])
    balanced = matched_ratio <= balance_ratio

    if min_effect is None:
        min_effect = _permutation_effect_floor(metric, data, group_col, groups, n_perm, rng)

    detail = dict(unmatched_effect=unmatched, matched_effect=matched,
                  overlap_range=m["overlap_range"], n_a=m["n_a"], n_b=m["n_b"],
                  ks_p=m["ks_p"], matched_ratio=matched_ratio, min_effect=min_effect,
                  retained_frac=float(retained))

    if abs(unmatched) < min_effect:
        return GateResult(
            2, "n_genes matching", GateStatus.PASS,
            f"no meaningful pre-matching effect (|{unmatched:.4g}| < null floor {min_effect:.4g}) "
            "— nothing to preserve; the groups simply do not differ on this metric",
            detail,
        )
    if np.sign(matched) != np.sign(unmatched) or retained < keep_frac:
        return GateResult(
            2, "n_genes matching", GateStatus.FAIL,
            f"effect collapses under matching: {unmatched:.4g} -> {matched:.4g} "
            f"({retained:.0%} retained, need >= {keep_frac:.0%}) — likely a QC artifact",
            detail,
        )
    if retained > max_frac:
        return GateResult(
            2, "n_genes matching", GateStatus.FAIL,
            f"effect AMPLIFIES under matching: {unmatched:.4g} -> {matched:.4g} "
            f"({retained:.0%} > {max_frac:.0%}) — possible selection on the quality axis",
            detail,
        )
    if not balanced:
        return GateResult(
            2, "n_genes matching", GateStatus.FAIL,
            f"matched subset still imbalanced (median n_genes {matched_ratio:.2f}x > "
            f"{balance_ratio}x) — residual QC confound not removed",
            detail,
        )
    return GateResult(
        2, "n_genes matching", GateStatus.PASS,
        f"effect survives matching: {unmatched:.4g} -> {matched:.4g} "
        f"({retained:.0%} retained, subset balanced {matched_ratio:.2f}x)",
        detail,
    )


# --------------------------------------------------------------------------- #
# GATE 3 — visible in the raw data
# --------------------------------------------------------------------------- #
def gate3_raw_visibility(
    data,
    gene_a: str,
    gene_b: str,
    group_col: str,
    groups: tuple,
    save_to: str | None = None,
) -> GateResult:
    """Export raw scatter data per group and flag if a shape change is dropout-driven.

    This gate is JUDGMENT: it hands you the numbers (and optionally a PNG) to eyeball. It
    auto-hints when the two groups differ mostly in how many points sit on the zero axes.
    """
    from .core import gene_column

    a_mask = np.asarray(data.obs[group_col]) == groups[0]
    b_mask = np.asarray(data.obs[group_col]) == groups[1]
    xa, ya = gene_column(data[a_mask], gene_a), gene_column(data[a_mask], gene_b)
    xb, yb = gene_column(data[b_mask], gene_a), gene_column(data[b_mask], gene_b)

    def _zero_frac(x, y):
        return float(np.mean((x <= 0) | (y <= 0)))

    zf_a, zf_b = _zero_frac(xa, ya), _zero_frac(xb, yb)
    dropout_driven = abs(zf_a - zf_b) > 0.15
    detail = dict(
        group_a=dict(n=int(a_mask.sum()), zero_frac=zf_a),
        group_b=dict(n=int(b_mask.sum()), zero_frac=zf_b),
        dropout_gap=abs(zf_a - zf_b),
        scatter={"a": (xa.tolist(), ya.tolist()), "b": (xb.tolist(), yb.tolist())}
        if save_to is None else "written",
    )

    if save_to is not None:
        try:  # optional plotting
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharex=True, sharey=True)
            for ax, (x, y, lab) in zip(axes, [(xa, ya, groups[0]), (xb, yb, groups[1])]):
                ax.scatter(x, y, s=6, alpha=0.3)
                ax.set_title(f"{lab} (zero-frac {(_zero_frac(x, y)):.2f})")
                ax.set_xlabel(gene_a); ax.set_ylabel(gene_b)
            fig.tight_layout(); fig.savefig(save_to, dpi=120); plt.close(fig)
            detail["figure"] = save_to
        except Exception as e:  # pragma: no cover
            detail["figure_error"] = str(e)

    msg = ("shape change between groups is likely dropout-driven "
           f"(zero-frac {zf_a:.2f} vs {zf_b:.2f}) — inspect the scatter"
           if dropout_driven else
           f"zero-fractions comparable ({zf_a:.2f} vs {zf_b:.2f}) — inspect the scatter for real coupling")
    return GateResult(3, "Raw visibility", GateStatus.JUDGMENT, msg, detail)




# --------------------------------------------------------------------------- #
# GATE 4 (automatable part) — does the metric respond to an injected signal?
# --------------------------------------------------------------------------- #
def gate4_signal_response(
    metric: Metric,
    data,
    inject: Callable,
    direction: str = "increase",
    n_rep: int = 10,
    z_min: float = 3.0,
    seed: int = 0,
) -> GateResult:
    """Inject a known construct change (``injected_signal``) and require the metric to move,
    reliably and in the declared direction. A metric that ignores its construct (a random
    number, a constant) cannot pass this.
    """
    from .injected_signal import describe

    sd = _as_simple(data)
    rng = np.random.default_rng(seed)
    base = _safe_call(lambda: metric(sd))
    if base is None:
        return GateResult(4, "Construct response (injected signal)", GateStatus.SKIP,
                          "metric could not be evaluated on the data", {})
    deltas = []
    for _ in range(n_rep):
        try:
            injected = inject(sd, rng)
        except Exception as e:  # e.g. non-count input
            return GateResult(4, "Construct response (injected signal)", GateStatus.SKIP,
                              f"signal could not be injected: {e}", {})
        v = _safe_call(lambda: metric(injected))
        if v is not None:
            deltas.append(v - base)
    if len(deltas) < 3:
        return GateResult(4, "Construct response (injected signal)", GateStatus.SKIP,
                          "metric could not be evaluated on injected data", {})
    deltas = np.asarray(deltas)
    mean, sd_ = float(deltas.mean()), float(deltas.std(ddof=1))
    se = sd_ / np.sqrt(len(deltas))
    z = mean / se if se > 0 else (float("inf") * np.sign(mean) if mean else 0.0)
    signed_z = z if direction == "increase" else -z
    detail = dict(injection=describe(inject), direction=direction, base=base,
                  mean_response=mean, sd_response=sd_, z=float(z), n_rep=len(deltas), z_min=z_min)
    if signed_z >= z_min:
        return GateResult(4, "Construct response (injected signal)", GateStatus.PASS,
                          f"responds to {describe(inject)}: {mean:+.4g} (z={z:.1f}, expected {direction})",
                          detail)
    return GateResult(4, "Construct response (injected signal)", GateStatus.FAIL,
                      f"does not respond to {describe(inject)} as declared ({direction}): "
                      f"{mean:+.4g} (z={z:.1f}, need z >= {z_min} in that direction)", detail)


# --------------------------------------------------------------------------- #
# GATE 5 — controls against an empirical null, in every stratum
# --------------------------------------------------------------------------- #
MATCHED_NEIGHBOURS_MIN = 20


def _neighbourhood(means: np.ndarray, names: list, target: str, exclude: set,
                   k_min: int = MATCHED_NEIGHBOURS_MIN, frac: float = 0.05) -> list:
    """The genes closest to `target` in mean expression (at least `k_min`, or `frac` of all)."""
    pos = names.index(target)
    k = max(k_min, int(np.ceil(frac * len(names))))
    dist = np.abs(means - means[pos])
    order = [j for j in np.argsort(dist, kind="mergesort") if names[j] not in exclude and j != pos]
    return [names[j] for j in order[:k]]


def _matched_pool(sub, pair, exclude, rng) -> list:
    """Distinct unrelated gene pairs from the expression neighbourhoods of `pair`, shuffled.

    Pairs are unordered and never repeated: the null's resolution is limited by the number of
    distinct pairs, and drawing the same pair twice would overstate it (a small panel offers
    only a few hundred)."""
    X = as_dense(sub.X)
    names = [str(g) for g in sub.var_names]
    means = X.mean(axis=0)
    na = _neighbourhood(means, names, str(pair[0]), exclude)
    nb = _neighbourhood(means, names, str(pair[1]), exclude)
    pool = sorted({tuple(sorted((a, b))) for a in na for b in nb if a != b})
    return [pool[i] for i in rng.permutation(len(pool))]


def _matched_null_draw(pair_metric, sub, pool: list):
    """`draw(k)` for ``stats.extend_null``: the metric on the next `k` pairs of the pool
    (without replacement; fewer when the pool runs out)."""
    state = dict(next=0)

    def draw(k):
        i = state["next"]
        state["next"] = i + int(k)
        vals = []
        for a, b in pool[i:i + int(k)]:
            v = _safe_call(lambda: pair_metric(sub, gene_a=a, gene_b=b))
            if v is not None:
                vals.append(v)
        return np.asarray(vals)
    return draw


DEPTH_SWAP_WINDOW = 20


def _depth_matched_draw(b: np.ndarray, depth: np.ndarray, rng, window: int = DEPTH_SWAP_WINDOW) -> np.ndarray:
    """Gene b's counts with its link to every other gene destroyed and its depth dependence kept.

    Each cell receives the count of a random cell among the next `window` cells at least as
    deep, binomially thinned to its own depth (``Binomial(b_j, depth_i / depth_j)``). Thinning
    a count to a lower depth is exact for sampling, so the draw follows b's depth dependence
    cell by cell, not just per bin. A shuffle *within depth bins* leaves the depth variation
    inside each bin, and a pair coupled only through depth then passed as a positive control
    (p = 0.005 in 5 of 5 null datasets with a log-depth SD of 0.8). The deepest cell has no
    deeper neighbour and keeps its own count (conservative: one cell keeps its coupling).
    """
    n = len(b)
    depth = np.asarray(depth, float)
    order = np.argsort(depth, kind="mergesort")
    rank = np.empty(n, dtype=int)
    rank[order] = np.arange(n)
    width_up = np.minimum(window, n - 1 - rank)
    u = rng.random(n)
    partner = np.where(width_up > 0, rank + 1 + np.floor(u * np.maximum(width_up, 1)).astype(int), rank)
    j = order[np.clip(partner, 0, n - 1)]
    ratio = np.where(depth[j] > 0, depth / np.where(depth[j] > 0, depth[j], 1.0), 1.0)
    return rng.binomial(np.round(np.asarray(b, float)[j]).astype(np.int64),
                        np.clip(ratio, 0.0, 1.0)).astype(float)


def _pair_depth(X: np.ndarray, ia: int, ib: int) -> np.ndarray:
    """Per-cell depth for the pair's null: total counts of all *other* genes."""
    return X.sum(axis=1) - X[:, ia] - X[:, ib]


def _shuffled_null(pair_metric, sub, pair, n_null, rng):
    """Self null: the same two genes, with gene b replaced by a depth-matched draw from
    neighbouring cells (``_depth_matched_draw``).

    Destroys the pair's biological coupling while keeping both genes' distributions and their
    dependence on depth. This is the reference for a positive control: does the metric register
    the known coupling above its technical baseline? An expression-matched pair null is
    unsuitable here, because other truly coupled genes contaminate it. Non-count input falls
    back to a shuffle within depth bins.
    """
    sd = _as_simple(sub)
    names = list(map(str, sd.var_names))
    ia, jb = names.index(str(pair[0])), names.index(str(pair[1]))
    original = sd.X[:, jb].copy()
    counts = _looks_like_counts(sd.X[:, [ia, jb]])
    depth = _pair_depth(sd.X, ia, jb)
    groups = None
    if not counts:
        bins = _depth_bins(sd)
        groups = [np.where(bins == k)[0] for k in np.unique(bins)]
    vals = []
    for _ in range(n_null):
        sd.X[:, jb] = (_depth_matched_draw(original, depth, rng) if counts
                       else _shuffle_within(original, groups, rng))
        v = _safe_call(lambda: pair_metric(sd, gene_a=pair[0], gene_b=pair[1]))
        if v is not None:
            vals.append(v)
    sd.X[:, jb] = original
    return np.asarray(vals)


POSITIVE_CONTROL_DOSE = 2.0  # induces a Spearman correlation of ~0.34 between genes with ~7 counts per cell
POSITIVE_CONTROL_POWER = 0.8


def _shuffle_within(col: np.ndarray, groups: list, rng) -> np.ndarray:
    out = col.copy()
    for idx in groups:
        if len(idx) > 1:
            out[idx] = col[rng.permutation(idx)]
    return out


def _ranks(x: np.ndarray) -> np.ndarray:
    _, inv, counts = np.unique(x, return_inverse=True, return_counts=True)
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x))
    r[order] = np.arange(len(x))
    return (np.bincount(inv, weights=r) / counts)[inv]


def _ref_coupling(a: np.ndarray, b: np.ndarray, rest=None) -> float:
    """Reference coupling detector for the power check: Spearman correlation of raw counts.

    Rank-based and free of normalization artifacts (closure); the depth-matched self-null
    carries the depth dependence. Among log-CP10k Pearson, raw Pearson, raw Spearman and
    log1p Pearson it had the highest power for an injected coupling in the dev data, together
    with log1p Pearson, and it is the least sensitive to outliers.
    """
    ra, rb = _ranks(np.asarray(a, float)), _ranks(np.asarray(b, float))
    ra, rb = ra - ra.mean(), rb - rb.mean()
    den = float(np.sqrt((ra * ra).sum() * (rb * rb).sum()))
    return float((ra * rb).sum() / den) if den > 0 else 0.0


def _positive_control_power(sub, pair, dose: float, alpha_s: float, rng, n_rep: int = 40,
                            n_null: int = 200, threshold: float = POSITIVE_CONTROL_POWER) -> float:
    """Power of this stratum to establish a coupling of known `dose` between the control genes.

    The pair's own coupling is destroyed (gene b replaced by a depth-matched draw), a coupling
    of `dose` is injected (``injected_signal.coupling``, the same injection as GATE 4), and a
    *reference* detector is tested against the same depth-matched self-null at the same
    alpha/K. The power belongs to the design, not to the metric under test: measured with the
    metric itself, a blind metric would always look underpowered and stay UNTESTED forever.
    Replicates stop as soon as ``power >= threshold`` is decided, so the returned estimate is
    exact for that decision but not unbiased.
    """
    from .injected_signal import coupling

    X = as_dense(sub.X)
    names = [str(g) for g in sub.var_names]
    ia, ib = names.index(str(pair[0])), names.index(str(pair[1]))
    if not _looks_like_counts(X[:, [ia, ib]]):
        return float("nan")  # the injection needs counts
    a0, b0 = X[:, ia].astype(float), X[:, ib].astype(float)
    rest = _pair_depth(X, ia, ib)
    obs = pd.DataFrame(index=range(len(a0)))
    inject = coupling("a", "b", strength=dose)
    need = int(np.ceil(threshold * n_rep))
    hits = done = 0
    for _ in range(n_rep):
        base = SimpleData(np.column_stack([a0, _depth_matched_draw(b0, rest, rng), rest]), obs,
                          ["a", "b", "rest"])
        inj = inject(base, rng).X
        a, b = inj[:, 0], inj[:, 1]
        x = _ref_coupling(a, b, rest)

        def draw(k, a=a, b=b):
            return np.array([_ref_coupling(a, _depth_matched_draw(b, rest, rng), rest)
                             for _ in range(int(k))])
        _, p, _ = extend_null(x, draw, n_null, alpha_s)
        hits += bool(np.isfinite(p) and p < alpha_s)
        done += 1
        # stop once the remaining replicates cannot change "power >= threshold"
        if hits >= need or hits + (n_rep - done) < need:
            break
    return hits / done


def gate5_controls(
    pair_metric: Callable,
    data,
    pos_pair: tuple,
    neg_pair: tuple,
    within: Sequence[str] = (),
    pos_min: float | None = None,
    neg_max: float | None = None,
    min_cells: int = 10,
    n_null: int = 200,
    alpha: float = 0.05,
    exclude: Sequence[str] = (),
    seed: int = 0,
    pos_dose: float = POSITIVE_CONTROL_DOSE,
    pos_power_min: float = POSITIVE_CONTROL_POWER,
    n_power: int = 40,
) -> GateResult:
    """Positive control must fire and negative control must stay null — in every stratum —
    judged against **empirical nulls** instead of a fixed band.

    v0.1 used one negative pair and a band of 20% of the positive control, with no sampling
    noise. Library-size normalization alone makes unrelated genes correlate (closure), so the
    "robust" reference metric failed its own negative control on null data.

    * Negative control: compared with `n_null` distinct unrelated pairs from the same
      expression neighbourhoods (the 20 genes closest in mean, or 5% of genes). It must not
      stand out (two-sided p >= alpha/K, K = assessable strata, Bonferroni). The null centre
      is reported: a centre far from zero means the metric reports association between
      unrelated genes. With few genes the pool of distinct pairs bounds the resolution; when
      alpha/K is below it, the negative control cannot fail and the message says so.
    * Positive control: compared with `n_null` self-nulls of the *same* pair, gene b replaced
      by a depth-matched draw from neighbouring cells (thinned to each cell's depth). Coupling
      is destroyed, the dependence on depth kept exactly. It must stand out (p < alpha/K).

    p values are rank-based Monte Carlo p values. When a control lies in the extreme tail of
    the first `n_null` draws and alpha/K is below their resolution, the null is extended to
    ``2K/alpha`` draws (at most 5000; ``stats.extend_null``) instead of assuming a normal tail.

    Status:

    * FAIL when the negative control stands out in any stratum (the metric reports
      association where there is none), or when the positive control does not stand out in a
      stratum that had the power to show it: the stratum establishes an injected coupling of
      dose `pos_dose` with power >= `pos_power_min` (``_positive_control_power``; defaults 2.0
      and 0.8, pre-registered as ``positive_control_dose`` / ``positive_control_power``). The
      metric is then insensitive to the coupling (or the control is not coupled there).
    * WARN when the positive control is silent only in strata without that power: absence of
      evidence. If it fires nowhere, ``detail["pos_demonstrated"]`` is False and the metric's
      response stays undemonstrated (UNTESTED) unless an injected signal shows it.
    * PASS when both controls behave in every stratum.

    Passing explicit `pos_min` / `neg_max` selects the legacy fixed band (not the default),
    where a silent positive control FAILs without a power check.
    """
    import itertools

    obs = data.obs
    within = list(within)
    if within:
        levels = [sorted(pd.unique(obs[f].dropna())) for f in within]
        combos = list(itertools.product(*levels))
    else:
        combos = [()]
    legacy = pos_min is not None or neg_max is not None
    rng = np.random.default_rng(seed)
    excl = set(map(str, pos_pair)) | set(map(str, neg_pair)) | set(map(str, exclude))

    def _val(sub, pair):
        return float(pair_metric(sub, gene_a=pair[0], gene_b=pair[1]))

    strata = []
    for combo in combos:
        mask = np.ones(len(obs), dtype=bool)
        for f, v in zip(within, combo):
            mask &= (np.asarray(obs[f]) == v)
        if mask.sum() >= min_cells:
            strata.append((dict(zip(within, combo)) if within else "pooled", mask))
    if not strata:
        return GateResult(
            5, "Controls", GateStatus.SKIP,
            f"no stratum had >= {min_cells} cells (checked {len(combos)}) — controls not evaluable",
            dict(rows=[]),
        )
    alpha_s = alpha / len(strata)

    rows, failures = [], []
    for label, mask in strata:
        sub = data[mask]
        pos, neg = _val(sub, pos_pair), _val(sub, neg_pair)
        row = dict(stratum=label, n_cells=int(mask.sum()), pos=pos, neg=neg)
        if legacy:
            pm = pos_min if pos_min is not None else 0.0
            nm = neg_max if neg_max is not None else float("inf")
            row.update(pos_fires=abs(pos) > pm, neg_ok=abs(neg) <= nm, method="fixed band")
        else:
            pool = _matched_pool(sub, neg_pair, excl, rng)
            null_neg, p_neg, m_neg = extend_null(neg, _matched_null_draw(pair_metric, sub, pool),
                                                 n_null, alpha_s)
            null_pos, p_pos, m_pos = extend_null(
                pos, lambda k, sub=sub: _shuffled_null(pair_metric, sub, pos_pair, k, rng),
                n_null, alpha_s)
            row.update(
                null_center=float(np.median(null_neg)) if len(null_neg) else float("nan"),
                null_sd=float(np.std(null_neg, ddof=1)) if len(null_neg) > 1 else float("nan"),
                pos_null_center=float(np.median(null_pos)) if len(null_pos) else float("nan"),
                n_null=(int(len(null_neg)), int(len(null_pos))), n_pairs_available=len(pool),
                p_neg=p_neg, p_pos=p_pos, method=f"{m_neg}/{m_pos}",
                pos_fires=bool(np.isfinite(p_pos) and p_pos < alpha_s),
                neg_ok=bool(not np.isfinite(p_neg) or p_neg >= alpha_s),
            )
        if not legacy and not row["pos_fires"]:
            power = _positive_control_power(sub, pos_pair, pos_dose, alpha_s, rng, n_rep=n_power,
                                            n_null=n_null, threshold=pos_power_min)
            row.update(pos_power=power, pos_power_dose=pos_dose, pos_power_reps=n_power,
                       pos_insensitive=bool(np.isfinite(power) and power >= pos_power_min))
        row["ok"] = bool(row["pos_fires"] and row["neg_ok"])
        rows.append(row)
        if not row["ok"]:
            failures.append(row)

    n_pos = sum(r["pos_fires"] for r in rows)
    detail = dict(rows=rows, alpha=alpha, alpha_per_stratum=alpha_s, n_null=n_null,
                  legacy_band=legacy, pos_min=pos_min, neg_max=neg_max,
                  pos_demonstrated=bool(n_pos > 0), n_pos_fires=int(n_pos),
                  pos_dose=pos_dose, pos_power_min=pos_power_min)
    centres = [r.get("null_center") for r in rows if r.get("null_center") is not None]
    centre_note = (f"; null centre {np.nanmin(centres):.3g}..{np.nanmax(centres):.3g}"
                   if centres and not legacy else "")
    weak_neg = [r for r in rows if "unattainable" in str(r.get("method", "")).split("/")[0]]
    detail["neg_test_unattainable"] = [r["stratum"] for r in weak_neg]
    if weak_neg:
        centre_note += (f"; in {len(weak_neg)}/{len(rows)} strata the negative-control test "
                        f"cannot reach alpha {alpha_s:.3g} (only {min(r['n_null'][0] for r in weak_neg)} "
                        "distinct unrelated pairs)")
    neg_bad = [r for r in rows if not r["neg_ok"]]
    pos_miss = [r for r in rows if not r["pos_fires"]]
    insensitive = [r for r in pos_miss if r.get("pos_insensitive")]
    if neg_bad or insensitive or (legacy and pos_miss):
        f0 = (neg_bad or insensitive or pos_miss)[0]
        why = []
        if not f0["neg_ok"]:
            why.append(f"negative control {f0['neg']:.3g} outside the null"
                       + (f" (p={f0['p_neg']:.3g}, null centre {f0['null_center']:.3g})" if not legacy else ""))
        if f0.get("pos_insensitive"):
            why.append(f"positive control {f0['pos']:.3g} inside its null (p={f0['p_pos']:.3g}) although "
                       f"the stratum establishes an injected coupling of dose {pos_dose:g} with power "
                       f"{f0['pos_power']:.2f} >= {pos_power_min:g} — the metric is insensitive to the "
                       "coupling (or the control is not coupled here)")
        if legacy and not f0["pos_fires"]:
            why.append(f"positive control {f0['pos']:.3g} below the band")
        n_bad = (len(failures) if legacy else
                 sum(1 for r in rows if not r["neg_ok"] or r.get("pos_insensitive")))
        return GateResult(
            5, "Controls", GateStatus.FAIL,
            f"{n_bad}/{len(rows)} strata fail controls (e.g. {f0['stratum']}: "
            + "; ".join(why) + ")",
            detail,
        )
    if pos_miss:
        where = ", ".join(f"{r['stratum']} (n={r['n_cells']}, p={r['p_pos']:.3g}, power "
                          f"{r.get('pos_power', float('nan')):.2f})" for r in pos_miss[:3])
        lead = (f"positive control beats its null in {n_pos}/{len(rows)} strata; not demonstrated in "
                if n_pos else "positive control not demonstrated in any stratum — ")
        return GateResult(
            5, "Controls", GateStatus.WARN,
            lead + where + ("…" if len(pos_miss) > 3 else "")
            + f" (power < {pos_power_min:g} for an injected coupling of dose {pos_dose:g} at alpha "
            f"{alpha_s:.3g}: absence of evidence, not a failure); negative control inside the null "
            f"in all strata{centre_note}",
            detail,
        )
    return GateResult(
        5, "Controls", GateStatus.PASS,
        f"positive fires and negative stays inside the empirical null in all {len(rows)} strata"
        + centre_note,
        detail,
    )


# --------------------------------------------------------------------------- #
# GATE 6 — replication of the effect on an independent dataset
# --------------------------------------------------------------------------- #
def gate6_replication(
    metric: Metric,
    data2,
    group_col: str,
    groups: tuple,
    within: Sequence[str] = (),
    thresh: float = 1.5,
    *,
    replicate_col: str | None = None,
    estimand: str | None = None,
    sesoi: float | None = None,
    primary_effect: float | None = None,
    alpha: float = 0.05,
    min_replicates: int = 3,
    n_perm: int = 1000,
    attenuation_half: float | None = None,
    seed: int = 0,
) -> GateResult:
    """Re-estimate the effect on an independent dataset with the same estimand, correction,
    replicate rule and strata. REPLICATED (PASS) = effect detected with the primary sign;
    NOT_REPLICATED (FAIL) = equivalent to zero within the SESOI, or detected with the opposite
    sign; anything else is INCONCLUSIVE (WARN). QC parity on the replication set is reported.
    """
    from .effect import estimate_effect

    if data2 is None:
        return GateResult(
            6, "Replication", GateStatus.SKIP,
            "no second dataset supplied — supply an independent AnnData to run this gate",
            dict(replication="NOT_RUN"),
        )
    g1 = gate1_qc_parity(data2, group_col, groups, within=within, thresh=thresh, seed=seed)
    eff, design = estimate_effect(
        metric, data2, group_col=group_col, groups=groups, within=within,
        replicate_col=replicate_col, estimand=estimand, sesoi=sesoi, alpha=alpha,
        min_replicates=min_replicates, n_perm=n_perm,
        qc_imbalanced=g1.status == GateStatus.WARN, attenuation_half=attenuation_half, seed=seed,
    )
    est = eff.detail.get("effect")
    same_sign = (primary_effect is None or est is None or not np.isfinite(est)
                 or np.sign(est) == np.sign(primary_effect))
    if g1.status == GateStatus.STOP or not design.get("identifiable", False):
        rep, status = "INCONCLUSIVE", GateStatus.WARN
        note = design.get("reason") or g1.message
    elif eff.status == "DETECTED" and same_sign:
        rep, status, note = "REPLICATED", GateStatus.PASS, eff.reason
    elif eff.status == "NO_DETECTABLE_EFFECT" or (eff.status == "DETECTED" and not same_sign):
        rep, status, note = "NOT_REPLICATED", GateStatus.FAIL, eff.reason
    else:
        rep, status, note = "INCONCLUSIVE", GateStatus.WARN, eff.reason
    return GateResult(
        6, "Replication", status,
        f"independent dataset: {rep} — {note} (QC parity {g1.status.value})",
        dict(replication=rep, effect=eff.detail, effect_status=eff.status, effect_flags=eff.flags,
             design=design, qc_parity=g1.detail, qc_status=g1.status.value, within=list(within)),
    )
