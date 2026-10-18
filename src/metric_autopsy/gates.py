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
from .stats import extend_null, t_ppf

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
BIAS_TOLERANCE = 0.5  # in SESOI units: a nuisance bias blocks only above this (and its 95% lower bound)


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
    sesoi: float | None = None,
    corrected: Sequence[str] = (),
    bias_tolerance: float | None = None,
) -> GateResult:
    """Does a nuisance *create or inflate* signal (bias), or only *shrink* it (attenuation)?

    Three references are measured on the user's own data: a bootstrap baseline (the metric's
    sampling spread), an automatic null (no gene-gene biology, same technology: for counts every
    gene takes, independently, the count of a random neighbouring cell at least as deep, thinned
    to the cell's depth, as in GATE 5; otherwise genes are permuted within depth bins), and each
    nuisance perturbation. Every shift is classified:

    * **bias** — the nuisance moves the null, reverses the signal, or inflates it. It can
      create a false effect. It blocks (GATE 0 FAILs) only where nothing else handles it and it
      is material for the claim: a bias in a nuisance that the declared correction removes
      between groups (`corrected`, e.g. depth under thinning) is reported, not failed; any
      other bias blocks only if its size exceeds `bias_tolerance` x `sesoi` and so does the
      lower bound of its 95% interval. Without a SESOI a bias cannot be sized for the claim:
      it is reported as unsized (and ``report.decide`` then withholds SUPPORTED);
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
    mode = _null_mode(sd)
    null_vals = [v for v in (_safe_call(lambda: metric(_gene_null(sd, rng, bins, mode)))
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

    if bias_tolerance is None:
        bias_tolerance = BIAS_TOLERANCE
    tol_abs = (float(bias_tolerance) * abs(float(sesoi))
               if sesoi is not None and np.isfinite(sesoi) and sesoi != 0 else None)
    corrected = set(corrected)
    # reference nulls by construction: a perturbation that leaves the counts non-integer
    # (library_scale) gets the depth-bin null, and is compared with the depth-bin null of the
    # unperturbed data, never with the depth-matched one
    refs = {mode: (null_center, null_sd, len(null_vals))} if null_center is not None else {}

    def _ref(m):
        if m not in refs:
            vals = [v for v in (_safe_call(lambda: metric(_gene_null(sd, rng, bins, m)))
                                for _ in range(n_null)) if v is not None]
            refs[m] = (*_mean_sd(vals), len(vals)) if len(vals) >= 2 else (None, None, 0)
        return refs[m]

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
            nvals, pmode = [], None
            for _ in range(min(5, n_null)):
                pdata, _ = _perturb(sd, kind, rng, protect)
                pmode = pmode or _null_mode(pdata)
                pbins = bins if pdata.n_obs == len(bins) else _depth_bins(pdata)
                v = _safe_call(lambda: metric(_gene_null(pdata, rng, pbins, pmode)))
                if v is not None:
                    nvals.append(v)
            ref_center, ref_sd, ref_n = _ref(pmode) if pmode else (None, None, 0)
            if len(nvals) >= 2 and ref_center is not None:
                ref_signal = base_mean - ref_center
                null_pert, null_pert_sd = _mean_sd(nvals)
                null_shift = null_pert - ref_center
                null_se = float(np.sqrt(ref_sd ** 2 / ref_n + null_pert_sd ** 2 / len(nvals)))
                null_z = abs(null_shift) / null_se if null_se > 0 else (float("inf") if null_shift else 0.0)
                signal_pert = pert_mean - null_pert
                r.update(null_perturbed=null_pert, null_shift=float(null_shift),
                         null_shift_std=float(null_shift / scale), null_shift_z=float(null_z),
                         signal=float(ref_signal), signal_perturbed=float(signal_pert),
                         null_construction=pmode)
                material = (lambda x, rel: abs(x) > rel or (tol_abs is not None and abs(x) > tol_abs))
                null_moves = null_z > z_thresh and material(null_shift, tol * scale)
                reverses = (np.sign(signal_pert) != np.sign(ref_signal)
                            and abs(signal_pert) > tol * abs(ref_signal))
                gain = abs(signal_pert) - abs(ref_signal)
                change = gain / abs(ref_signal) if ref_signal else float("inf")
                r["signal_change"] = float(change)
                if null_moves or reverses or material(max(gain, 0.0), tol * abs(ref_signal)):
                    cls = "bias"
                    r["bias_kind"] = ("null moves" if null_moves else
                                      "signal reverses" if reverses else "signal inflated")
                    size, size_se = ((abs(null_shift), null_se) if null_moves else (abs(shift), se))
                    r.update(bias_size=float(size), bias_se=float(size_se),
                             bias_ci_low=float(size - 1.96 * size_se))
                    if kind in corrected:
                        r["handling"] = "corrected"
                    elif tol_abs is None:
                        r["handling"] = "unsized"
                    else:
                        r["bias_tolerance"] = tol_abs
                        r["handling"] = ("blocking" if size > tol_abs and size - 1.96 * size_se > tol_abs
                                         else "below tolerance")
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
        r["confounded"] = cls == "bias" and r.get("handling") == "blocking"
        responses[kind] = r

    detail = dict(baseline=base_mean, baseline_sd=base_sd, null_center=null_center,
                  null_sd=null_sd, signal=signal, has_signal=bool(has_signal),
                  immaterial_signal=immaterial_signal, scale=scale,
                  scale_source=scale_source, effect_scale=effect_scale, responses=responses,
                  attenuation=attenuation, level_shifts=level_shifts, tol=tol, z_thresh=z_thresh,
                  null_construction=mode, sesoi=sesoi, bias_tolerance=float(bias_tolerance),
                  bias_tolerance_abs=tol_abs, corrected=sorted(corrected),
                  bias_handling={k: r["handling"] for k, r in responses.items() if r.get("handling")})

    spread = float(np.ptp(all_vals)) if all_vals else 0.0
    if spread <= 1e-12 * max(1.0, abs(base_mean)):
        detail["degenerate"] = True
        return GateResult(
            0, "Mathematical independence", GateStatus.DEGENERATE,
            f"the metric returned the same value ({base_mean:.6g}) on every resample, null and "
            "perturbation — it carries no information about the data", detail)

    flagged = {k: r for k, r in responses.items() if r.get("confounded")}
    if flagged:
        worst = max(flagged, key=lambda k: flagged[k]["bias_size"] / flagged[k]["bias_tolerance"])
        w = flagged[worst]
        return GateResult(
            0, "Mathematical independence", GateStatus.FAIL,
            f"'{worst}' biases the metric ({w['bias_kind']}): {w['bias_size']:.4g} "
            f"(95% lower bound {w['bias_ci_low']:.4g}) against a tolerance of {w['bias_tolerance']:.4g} "
            f"({float(bias_tolerance):g} x SESOI; z={w['z']:.1f}) — the nuisance can create a "
            "difference with no biology",
            detail,
        )
    handled = {k: r for k, r in responses.items() if r.get("classification") == "bias"}
    notes = []
    for k, r in handled.items():
        if r["handling"] == "corrected":
            notes.append(f"'{k}' biases the metric ({r['bias_kind']}, {r['bias_size']:.4g}); the declared "
                         "correction removes it between groups")
        elif r["handling"] == "below tolerance":
            notes.append(f"'{k}' biases the metric ({r['bias_kind']}, {r['bias_size']:.4g}, 95% lower "
                         f"bound {r['bias_ci_low']:.4g}) within the tolerance {r['bias_tolerance']:.4g}")
        else:
            notes.append(f"'{k}' biases the metric ({r['bias_kind']}, {r['bias_size']:.4g}); without a "
                         "SESOI it cannot be sized against the claim: declare a SESOI (the smallest "
                         "effect that matters, on the metric's scale) in the pre-registration, and the "
                         "bias is judged against bias_tolerance x SESOI")
    parts = []
    if attenuation:
        parts.append("attenuation " + ", ".join(f"{k} −{v:.0%}" for k, v in attenuation.items())
                     + " of the signal (passed to design adequacy as a power check)")
    if level_shifts:
        parts.append("level shift " + ", ".join(f"{k} {v:+.0%}" for k, v in level_shifts.items())
                     + " of the effect scale (no structure to classify it; the effect field equalizes "
                     "depth between groups)")
    if notes:
        msg = "; ".join(notes + parts) + " — reported, not failed"
        return GateResult(0, "Mathematical independence", GateStatus.WARN, msg, detail)
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
GATE4_N_REP = 200    # injections (each against its own reference) of the all-cell response
GATE4_N_REP_REPLICATE = 50  # injections into each unit's cells for the claim's statistic (v0.4)
GATE4_N_REP_CURVE = 20  # injections into each unit's cells at the other doses of the dose-response curve
GATE4_CURVE = (0.5, 2.0)  # the curve's construct changes, as multiples of the SESOI (beside the SESOI itself)
DELTA_MIN_FRACTION = 0.5  # v0.3's default delta_min = 0.5 x SESOI; no longer used (v0.4 change 3)


def response_interval(deltas, direction: str = "increase", alpha: float = 0.05) -> dict:
    """Mean of the paired responses (injected minus sham) with its two-sided (1 - alpha) t
    interval, and both bounds signed so that the declared direction is positive."""
    if direction not in ("increase", "decrease"):
        raise ValueError(f"the signal direction is 'increase' or 'decrease', not {direction!r}")
    d = np.asarray(deltas, float)
    n = len(d)
    mean = float(d.mean())
    sd_ = float(d.std(ddof=1)) if n > 1 else float("nan")
    se = sd_ / np.sqrt(n)
    half = t_ppf(1.0 - alpha / 2.0, n - 1) * se if se > 0 else 0.0
    lo, hi = mean - half, mean + half
    sign = 1.0 if direction == "increase" else -1.0
    s_lo, s_hi = sorted((sign * lo, sign * hi))
    return dict(mean_response=mean, sd_response=sd_, n_rep=n, ci=(lo, hi), ci_level=1.0 - alpha,
                signed_lower=s_lo, signed_upper=s_hi, z=float(mean / se) if se > 0 else
                (float("inf") * np.sign(mean) if mean else 0.0))


def _interval(mean: float, se: float, df: float, n: int, direction: str, alpha: float) -> dict:
    """``response_interval``'s fields for a mean with a given standard error and degrees of freedom
    (the interval over replicates of GATE 4's claim-level response)."""
    if direction not in ("increase", "decrease"):
        raise ValueError(f"the signal direction is 'increase' or 'decrease', not {direction!r}")
    if not (np.isfinite(se) and df > 0):
        return dict(mean_response=mean, se=se, df=df, n_rep=n, ci=(float("nan"), float("nan")),
                    ci_level=1.0 - alpha, signed_lower=float("nan"), signed_upper=float("nan"), z=float("nan"))
    half = t_ppf(1.0 - alpha / 2.0, df) * se if se > 0 else 0.0
    lo, hi = mean - half, mean + half
    sign = 1.0 if direction == "increase" else -1.0
    s_lo, s_hi = sorted((sign * lo, sign * hi))
    return dict(mean_response=mean, se=se, df=df, n_rep=n, ci=(lo, hi), ci_level=1.0 - alpha,
                signed_lower=s_lo, signed_upper=s_hi, z=float(mean / se) if se > 0 else
                (float("inf") * np.sign(mean) if mean else 0.0))


def judge_response(iv: dict, delta_min: float | None = None) -> str:
    """GATE 4 decides only whether the metric is blind (v0.4 change 3, decided 2026-10-17): PASS when
    a response is shown, the interval's lower bound in the declared direction above 0; FAIL when the
    metric is shown not to respond, the whole interval at or below 0 in that direction; otherwise
    UNTESTED, absence of evidence. `delta_min` is not used any more: what a change that matters is
    in the metric's units is the response itself (s_m), which the effect field uses as the SESOI."""
    if iv["signed_lower"] > 0:
        return "PASS"
    if iv["signed_upper"] <= 0:
        return "FAIL"
    return "UNTESTED"


def _check_delta_min(delta_min):
    if delta_min is None:
        return None
    delta_min = float(delta_min)
    if not (np.isfinite(delta_min) and delta_min > 0):
        raise ValueError(f"delta_min must be a positive number on the metric's scale, not {delta_min!r}")
    return delta_min


def _column_local(inject) -> bool:
    """An injector that rewrites only its genes' columns (``inject.genes``, ``inject.columns``),
    against a sham that does the same (``coupling``) or against the untouched data (``module``,
    which has no sham)."""
    sham = getattr(inject, "sham", None)
    return (getattr(inject, "genes", None) is not None and callable(getattr(inject, "columns", None))
            and (sham is None or callable(getattr(sham, "columns", None))))


def injected_deltas(metric, data, inject, n_rep: int, rng, *, decouple: str | None = None) -> list:
    """`n_rep` paired responses: the metric on injected data minus the metric on the injector's
    sham (or on the untouched data if the injector has no sham). Raises ValueError if the signal
    cannot be injected (e.g. non-count input).

    A column-local injector (``inject.genes`` and ``inject.columns``: ``injected_signal.coupling``
    and ``module``) rewrites its genes in one private copy of the data and restores them after
    every evaluation, with the same random draws as injecting into copies: the metric sees exactly
    the data a copy would hold, without copying the matrix twice per injection.
    `decouple` names the gene replaced, before every injection, by a depth-matched draw from
    neighbouring cells (``_depth_matched_draw``): GATE 5 injects into the positive control's
    genes with their own coupling removed. It needs a column-local injector with a sham.
    """
    sham = getattr(inject, "sham", None)
    if not _column_local(inject):
        if decouple is not None:
            raise ValueError("decoupling needs a column-local injector")
        sd = _as_simple(data)
        base = None if sham is not None else _safe_call(lambda: metric(sd))
        out = []
        for _ in range(n_rep):
            injected = inject(sd, rng)
            reference = sham(sd, rng) if sham is not None else None
            v = _safe_call(lambda: metric(injected))
            ref = base if reference is None else _safe_call(lambda: metric(reference))
            if v is not None and ref is not None and np.isfinite(v - ref):
                out.append(v - ref)
        return out
    if decouple is not None and sham is None:
        raise ValueError("decoupling needs an injector with a sham")
    work = SimpleData(np.array(as_dense(data.X), copy=True), data.obs.copy(), list(data.var_names))
    return _column_deltas(metric, work, inject, n_rep, rng, decouple=decouple)


def _column_deltas(metric, work: SimpleData, inject, n_rep: int, rng, *, decouple: str | None = None) -> list:
    """The column-local responses on `work`, a private copy that is restored before returning."""
    from .core import unique_col_index

    sham = getattr(inject, "sham", None)
    X, obs = work.X, work.obs
    cols = [unique_col_index(work.var_names, g) for g in inject.genes]
    original = X[:, cols].copy()
    if not _looks_like_counts(original):
        raise ValueError("the injection needs raw counts (non-negative integers)")
    tot = np.asarray(obs["total_counts"], dtype=float).copy() if "total_counts" in obs else None
    if decouple is not None:
        kb = list(inject.genes).index(decouple)
        depth = X.sum(axis=1) - original.sum(axis=1)
    base = _safe_call(lambda: metric(work)) if sham is None else None
    out = []
    try:
        for _ in range(n_rep):
            before = original.copy()
            if decouple is not None:
                before[:, kb] = _depth_matched_draw(original[:, kb], depth, rng)
            vals = []
            for fn in ((inject.columns,) if sham is None else (inject.columns, sham.columns)):
                after = fn(before, rng)
                X[:, cols] = after
                if tot is not None:
                    obs["total_counts"] = tot - (original - after).sum(axis=1)
                vals.append(_safe_call(lambda: metric(work)))
            v, ref = (vals[0], base) if sham is None else vals
            if v is not None and ref is not None and np.isfinite(v - ref):
                out.append(v - ref)
    finally:
        X[:, cols] = original
        if tot is not None:
            obs["total_counts"] = tot
    return out


def _replicate_delta(metric, sub: SimpleData, inject, n_rep: int, rng) -> tuple[float | None, int]:
    """The mean response of one unit's cells (`sub`, a private copy) over `n_rep` injections."""
    if _column_local(inject):
        deltas = _column_deltas(metric, sub, inject, n_rep, rng)
    else:
        deltas = injected_deltas(metric, sub, inject, n_rep, rng)
    return (float(np.mean(deltas)) if deltas else None), len(deltas)


def _replicate_units(sd: SimpleData, group_col, groups, replicate_col, within):
    """Without cell types, a unit is a replicate's cells (a paired replicate: both its groups' cells)."""
    from .combined import Units
    from .effect import _replicate_design

    design = _replicate_design(sd.obs, group_col, groups, replicate_col, list(within))
    if design["kind"] == "mixed":
        return None, design
    rep = np.asarray(sd.obs[replicate_col]).astype(str)
    grp = np.asarray(sd.obs[group_col])
    reps = list(design["reps"])
    strata = list(dict.fromkeys(design["strata"][r] for r in reps))
    s_index = {s: i for i, s in enumerate(strata)}
    first = np.array([1 if groups[0] in design["membership"][r] else 0 for r in reps]) if design["kind"] == "nested" \
        else np.zeros(len(reps), int)
    cells = []
    for r in reps:
        m = rep == r
        cells.append((np.where(m & (grp == groups[0]))[0], np.where(m & (grp == groups[1]))[0])
                     if design["kind"] == "paired" else np.where(m)[0])
    units = Units(design["kind"], reps, first, np.array([s_index[design["strata"][r]] for r in reps]), strata,
                  ["all cells"], np.arange(len(reps)), np.zeros(len(reps), int), cells)
    return units, design


def _unit_cells(units, i: int, k: int) -> np.ndarray:
    """Unit i's cells of group index k (a nested unit holds one group's cells)."""
    c = units.cells[i]
    return c[k] if units.kind == "paired" else c


def _unit_group(units, i: int) -> int:
    """0 for a unit of groups[0], 1 for groups[1] (nested)."""
    return 0 if units.rep_first[units.unit_rep[i]] == 1 else 1


def _bio_variances(X: np.ndarray, cols, units, k: int) -> list:
    """Each probed gene's biological variance of log expression, the mean over group k's units."""
    tot = X.sum(axis=1)
    out = []
    for j in cols:
        vals = []
        for i in range(units.n_units):
            if units.kind == "nested" and _unit_group(units, i) != k:
                continue
            idx = _unit_cells(units, i, k)
            if len(idx) >= 3:
                from .injected_signal import biological_log_variance
                v = biological_log_variance(X[idx, j], tot[idx])
                if np.isfinite(v):
                    vals.append(v)
        out.append(float(np.mean(vals)) if vals else float("nan"))
    return out


def _unit_responses(metric, sd: SimpleData, X: np.ndarray, units, k: int, inject, n_rep: int, rng):
    """Group k's units: the mean response of each over `n_rep` injections into its cells (NaN where
    fewer than half gave one), the share of their counts the injection removes, and the units used."""
    resp = np.full(units.n_units, np.nan)
    used = np.zeros(units.n_units, bool)
    depth = []
    for i in range(units.n_units):
        if units.kind == "nested" and _unit_group(units, i) != k:
            continue
        idx = _unit_cells(units, i, k)
        if len(idx) == 0:
            continue
        used[i] = True
        sub = SimpleData(np.array(X[idx], copy=True), sd.obs.iloc[idx].copy(), list(sd.var_names))
        v, n_ok = _replicate_delta(metric, sub, inject, n_rep, rng)
        if v is not None and np.isfinite(v) and n_ok >= max(2, n_rep // 2):
            resp[i] = v
        if callable(getattr(inject, "depth_loss", None)):
            depth.append(inject.depth_loss(sub))
    return resp, (float(np.mean(depth)) if depth else None), used


def _claim_response(metric, sd: SimpleData, inject, group_col, groups, replicate_col, within,
                    direction: str, alpha: float, n_rep: int, rng, *, domain=None,
                    construct_change: float | None = None, curve=GATE4_CURVE,
                    n_rep_curve: int = GATE4_N_REP_CURVE) -> dict:
    """GATE 4 on the claim's own statistic (v0.4 changes 1 and 3).

    The units are the claim's: with cell types, a replicate's cells of one type of the domain
    (``domain.build_units``), else a replicate's cells. The signal goes into the group whose units have
    more of the probed genes' counts per cell (each unit weighted equally; on an exact tie into each
    group, and the weaker response is judged). With a SESOI on the construct's scale the probe is
    re-dosed to a construct change of that size (``injected_signal.at_construct_change``; a coupling's
    strength from the injected group's biological variances). The response is the change of the
    injected group's combined value, types equally weighted within a stratum and the strata by the
    effect test's weights, as the effect field combines them; its interval is over the replicates."""
    from .combined import Estimator, one_group_mean
    from .core import unique_col_index
    from .injected_signal import at_construct_change, describe

    genes = getattr(inject, "genes", None)
    if not genes:
        return dict(error="the injector names no genes (inject.genes), so the group to inject into "
                          "cannot be chosen by their counts")
    if domain is not None:
        from .domain import build_units
        units = build_units(domain, sd, group_col=group_col, groups=groups, replicate_col=replicate_col)
        if units is None:
            return dict(error="the claim's domain holds no type with enough cells and replicates")
        kind = units.kind
    else:
        units, design = _replicate_units(sd, group_col, groups, replicate_col, within)
        if units is None:
            return dict(error=design["reason"])
        kind = units.kind
    X = as_dense(sd.X)
    cols = [unique_col_index(sd.var_names, g) for g in genes]
    probe = np.asarray(X[:, cols], dtype=float).sum(axis=1)
    level = {}
    for k, gname in enumerate(groups):
        per = []
        for i in range(units.n_units):
            if kind == "nested" and _unit_group(units, i) != k:
                continue
            idx = _unit_cells(units, i, k)
            if len(idx):
                per.append(float(probe[idx].mean()))
        level[str(gname)] = float(np.mean(per)) if per else float("nan")
    a, b = level[str(groups[0])], level[str(groups[1])]
    # on an exact tie both groups, in the order of their labels (a mirrored wording, the same draws)
    targets = [0] if a > b else [1] if b > a else sorted((0, 1), key=lambda k: str(groups[k]))
    weights = Estimator(units).weights if kind == "nested" else np.ones(1)
    runs = []
    for k in targets:
        gname = groups[k]
        probe_k, dose = inject, None
        if construct_change is not None:
            var = _bio_variances(X, cols, units, k) if getattr(inject, "kind", None) == "coupling" else [None, None]
            probe_k, dose = at_construct_change(inject, construct_change, var_a=var[0], var_b=var[1] if len(var) > 1 else None)
            if probe_k is None:
                return dict(error=f"no rule translates a SESOI on the construct's scale into {describe(inject)}'s dose")
        resp, depth, used = _unit_responses(metric, sd, X, units, k, probe_k, n_rep, rng)
        sub = units.subset(used)
        mean, se, df, n_reps = one_group_mean(resp[used][None, :], sub, k, weights)
        mean, se, df = float(mean[0]), float(se[0]), float(df[0])
        iv = _interval(mean, se, df, n_reps, direction, alpha)
        sign = 1.0 if k == 0 else -1.0  # the effect test's difference is groups[0] - groups[1]
        lo, hi = iv["ci"]
        label = (lambda i: str(units.reps[units.unit_rep[i]]) if len(units.types) == 1 and domain is None
                 else f"{units.reps[units.unit_rep[i]]} | {units.types[units.unit_type[i]]}")
        run = dict(group=str(gname), group_index=k, interval=iv, statistic_change=sign * mean,
                   statistic_ci=tuple(sorted((sign * lo, sign * hi))) if np.isfinite(lo) else (lo, hi),
                   depth_loss=depth, replicate_responses={label(i): float(resp[i]) for i in np.where(used)[0]
                                                          if np.isfinite(resp[i])},
                   probe=describe(probe_k), dose=dose, n_units=int(used.sum()))
        if construct_change is not None and curve:
            points = []
            for mult in curve:
                change = mult * construct_change
                if getattr(inject, "kind", None) == "coupling" and not 0 < change < 1:
                    continue
                pk, dk = at_construct_change(inject, change, var_a=var[0], var_b=var[1] if len(var) > 1 else None)
                r_c, _, used_c = _unit_responses(metric, sd, X, units, k, pk, n_rep_curve, rng)
                m_c, se_c, df_c, n_c = one_group_mean(r_c[used_c][None, :], units.subset(used_c), k, weights)
                iv_c = _interval(float(m_c[0]), float(se_c[0]), float(df_c[0]), n_c, direction, alpha)
                points.append(dict(construct_change=change, dose=dk, response=iv_c["mean_response"], ci=iv_c["ci"],
                                   n_rep=n_rep_curve))
            points.append(dict(construct_change=construct_change, dose=dose, response=mean, ci=(lo, hi), n_rep=n_rep))
            run["dose_response"] = sorted(points, key=lambda p: p["construct_change"])
        runs.append(run)
    run = min(runs, key=lambda x: (x["interval"]["signed_lower"] if np.isfinite(x["interval"]["signed_lower"])
                                   else -np.inf, x["interval"]["mean_response"] if direction == "increase"
                                   else -x["interval"]["mean_response"]))
    return dict(run, levels=level, tie=len(targets) == 2, design=kind,
                units=dict(n=units.n_units, types=list(units.types), domain=domain is not None),
                runs=[{k2: v2 for k2, v2 in x.items() if k2 != "replicate_responses"} for x in runs]
                if len(runs) > 1 else None)


def gate4_signal_response(
    metric: Metric,
    data,
    inject: Callable,
    direction: str | None = None,
    delta_min: float | None = None,
    n_rep: int = GATE4_N_REP,
    alpha: float = 0.05,
    seed: int = 0,
    *,
    group_col: str | None = None,
    groups: tuple | None = None,
    replicate_col: str | None = None,
    within: Sequence[str] = (),
    n_rep_replicate: int = GATE4_N_REP_REPLICATE,
    domain: dict | None = None,
    construct_change: float | None = None,
    curve=GATE4_CURVE,
    n_rep_curve: int = GATE4_N_REP_CURVE,
) -> GateResult:
    """Inject a known construct change (``injected_signal``) and judge whether the metric responds.

    With the claim's groups and replicate unit (``run_autopsy`` passes them), the response is
    measured with the claim's own statistic (v0.4 change 1): on the claim's units (with cell types, a
    replicate's cells of one type of the claim's ``domain``: GATE 4 runs on the domain's cells, not on
    the whole mixture), the signal goes into the group where the probed genes have more counts per
    cell, whatever the claim's wording, and the response is the change of that group's combined value,
    its (1 - `alpha`) interval taken over the replicates (`n_rep_replicate` injections within each
    unit). With `construct_change`, a SESOI on the construct's scale, the probe is re-dosed to a
    change of the construct of that size, and the response s_m is the SESOI in the metric's units
    (``detail["sesoi_metric"]``, v0.4 change 3); the doses of `curve` x the SESOI give the
    dose-response curve beside it. The response is judged against the direction the injection itself
    must move the metric (`direction`, by default the injector's ``inject.direction``), not the
    claim's. The all-cell response (v0.3's GATE 4: the metric over all the claim's cells, `n_rep`
    injections) stays in the detail as a diagnostic, and is the judged response when no replicate unit
    is given.

    The reference is the injector's sham (``coupling``: the same thinning without the signal) or, for
    an injector without one (``module``, which thins the module's own genes), the untouched data.

    * PASS — a response is shown: the lower bound in the declared direction is above 0.
    * FAIL — the metric is shown not to respond: the whole interval is at or below 0.
    * WARN (UNTESTED) — neither: absence of evidence, not invalidity.

    GATE 4 decides only whether the metric is blind (decided 2026-10-17). `delta_min`, v0.3's smallest
    response that matters, is accepted and not used: the response itself is the SESOI's size in the
    metric's units.
    """
    from .injected_signal import describe

    name = "Construct response (injected signal)"
    delta_min = _check_delta_min(delta_min)
    direction = direction or getattr(inject, "direction", None) or "increase"
    sd = _as_simple(data)
    claim_mode = group_col is not None and groups is not None and replicate_col is not None
    if claim_mode:
        gcol = np.asarray(sd.obs[group_col])
        keep = (gcol == groups[0]) | (gcol == groups[1])
        if domain is not None and not domain.get("error"):
            keep &= _domain_mask(sd.obs, domain, within)
        else:
            domain = None
        if not keep.any():
            return GateResult(4, name, GateStatus.SKIP, "no cell of the claim's domain to inject into", {})
        sd = sd[keep]
    rng = np.random.default_rng(seed)
    base = _safe_call(lambda: metric(sd))
    if base is None:
        return GateResult(4, name, GateStatus.SKIP, "metric could not be evaluated on the data", {})
    sham = getattr(inject, "sham", None)
    reference = "sham" if sham is not None else "untouched data"
    claim = None
    if claim_mode:
        try:
            claim = _claim_response(metric, sd, inject, group_col, groups, replicate_col, within,
                                    direction, alpha, n_rep_replicate, rng, domain=domain,
                                    construct_change=construct_change, curve=curve, n_rep_curve=n_rep_curve)
        except Exception as e:  # reported, never silent
            claim = dict(error=f"{type(e).__name__}: {e}")
    judged = inject
    if claim is not None and claim.get("error") is None and claim.get("dose") is not None:
        from .injected_signal import at_construct_change
        judged = at_construct_change(inject, construct_change,
                                     var_a=claim["dose"].get("var_a"), var_b=claim["dose"].get("var_b"))[0] or inject
    try:
        deltas = injected_deltas(metric, sd, judged, n_rep, rng)
    except Exception as e:  # e.g. non-count input
        return GateResult(4, name, GateStatus.SKIP, f"signal could not be injected: {e}", {})
    all_cells = response_interval(deltas, direction, alpha) if len(deltas) >= 3 else None
    detail = dict(injection=describe(judged), declared_injection=describe(inject), direction=direction, base=base,
                  reference=reference, alpha=alpha, all_cells=all_cells, construct_change=construct_change,
                  delta_min=None, delta_min_declared=delta_min)
    lvl = f"{100 * (1 - alpha):g}%"
    if claim is not None and claim.get("error") is None:
        iv = claim["interval"]
        detail.update(mode="claim statistic", injected_group=claim["group"], probed_counts_per_cell=claim["levels"],
                      tie=claim["tie"], replicate_responses=claim["replicate_responses"],
                      statistic_change=claim["statistic_change"], statistic_ci=claim["statistic_ci"],
                      depth_loss=claim["depth_loss"], n_rep_replicate=n_rep_replicate, runs=claim["runs"],
                      dose=claim["dose"], dose_response=claim.get("dose_response"), units=claim["units"], **iv)
        lo, hi = iv["ci"]
        others = {g: v for g, v in claim["levels"].items() if g != claim["group"]}
        where = ("on a tie of the probed counts the weaker of the two groups' responses is judged" if claim["tie"]
                 else f"{claim['levels'][claim['group']]:.4g} probed counts per cell against "
                      + ", ".join(f"{v:.4g} in {g}" for g, v in others.items()))
        on = (f"the domain's {len(claim['units']['types'])} type(s)" if claim["units"]["domain"] else "the replicates")
        resp = (f"{claim['probe']} injected into {claim['group']} ({where}), its values move by {iv['mean_response']:+.4g} "
                f"against the {reference} ({lvl} CI [{lo:+.4g}, {hi:+.4g}] over {iv['n_rep']} replicates, on "
                f"{on}, {n_rep_replicate} injections per unit, expected {direction}); the replicate-level difference "
                f"{groups[0]} - {groups[1]} moves by {claim['statistic_change']:+.4g}")
        if claim.get("depth_loss") is not None:
            resp += f"; the injection removes {claim['depth_loss']:.2%} of their counts"
        if claim.get("dose") is not None:
            sign = 1.0 if direction == "increase" else -1.0
            s_m = sign * iv["mean_response"]
            detail["sesoi_metric"] = float(s_m) if np.isfinite(s_m) and s_m > 0 else None
            detail["sesoi_metric_ci"] = tuple(sorted((sign * lo, sign * hi)))
            resp += (f"; the probe is a change of {construct_change:g} on the construct's scale "
                     f"({claim['dose']['scale']}), so the SESOI in the metric's units is s_m = {s_m:+.4g}")
        if not np.isfinite(iv["signed_lower"]):
            detail["outcome"] = "UNTESTED"
            return GateResult(4, name, GateStatus.WARN,
                              f"response to {describe(judged)} not established: {resp}; an interval over "
                              "the replicates needs at least two replicates with a response in one stratum "
                              "— untested, not invalid", detail)
    else:
        if all_cells is None:
            return GateResult(4, name, GateStatus.SKIP, "metric could not be evaluated on injected data", detail)
        iv = all_cells
        detail.update(mode="all cells", **iv)
        if claim is not None:
            detail["claim_statistic_error"] = claim["error"]
        lo, hi = iv["ci"]
        resp = (f"{iv['mean_response']:+.4g} against the {reference} over all cells "
                f"({lvl} CI [{lo:+.4g}, {hi:+.4g}], {iv['n_rep']} injections, expected {direction})")
        if claim is not None:
            resp += f"; the claim's statistic could not be used: {claim['error']}"
    outcome = judge_response(iv)
    detail["outcome"] = outcome
    if outcome == "FAIL":
        return GateResult(4, name, GateStatus.FAIL,
                          f"blind to {describe(judged)}: {resp}; the whole interval is at or below 0 in the "
                          "injection's direction", detail)
    if outcome == "PASS":
        return GateResult(4, name, GateStatus.PASS, f"responds to {describe(judged)}: {resp}", detail)
    return GateResult(4, name, GateStatus.WARN,
                      f"response to {describe(judged)} not established: {resp}; the interval includes 0 — "
                      "untested, not invalid", detail)


def _domain_mask(obs, domain: dict, within) -> np.ndarray:
    """The cells of the domain's (stratum, type) cells."""
    lab = np.asarray(obs[domain["column"]]).astype(str)
    if within:
        keys = obs[list(within)].astype(str).agg("\x1f".join, axis=1).to_numpy()
    else:
        keys = np.full(len(obs), "", dtype=object)
    mask = np.zeros(len(obs), bool)
    for st, types in domain["domain"].items():
        if types:
            mask |= (keys == st) & np.isin(lab, types)
    return mask


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


def _depth_matched_matrix(X: np.ndarray, depth: np.ndarray, rng, window: int = DEPTH_SWAP_WINDOW,
                          chunk: int = 256) -> np.ndarray:
    """`_depth_matched_draw` for every column of a count matrix, independently per gene.

    Cell i takes gene g's count from a random cell among the next `window` cells at least as
    deep (a fresh draw for each gene, so no cell's genes travel together), thinned to cell i's
    depth. Every gene keeps its dependence on depth; every gene-gene link is destroyed.
    """
    n, n_genes = X.shape
    depth = np.asarray(depth, float)
    order = np.argsort(depth, kind="mergesort")
    rank = np.empty(n, dtype=np.int64)
    rank[order] = np.arange(n)
    width_up = np.minimum(window, n - 1 - rank)
    up = (width_up > 0)[:, None]
    width = np.maximum(width_up, 1)[:, None]
    safe = np.where(depth > 0, depth, 1.0)
    out = np.empty((n, n_genes), dtype=float)
    for c0 in range(0, n_genes, chunk):
        c1 = min(n_genes, c0 + chunk)
        u = rng.random((n, c1 - c0))
        partner = np.where(up, rank[:, None] + 1 + np.floor(u * width).astype(np.int64), rank[:, None])
        j = order[np.clip(partner, 0, n - 1)]
        counts = np.round(X[j, np.arange(c0, c1)[None, :]]).astype(np.int64)
        ratio = np.where(depth[j] > 0, depth[:, None] / safe[j], 1.0)
        out[:, c0:c1] = rng.binomial(counts, np.clip(ratio, 0.0, 1.0))
    return out


def _null_mode(sd: SimpleData) -> str:
    return "depth-matched" if _looks_like_counts(sd.X) else "depth bins"


def _gene_null(sd: SimpleData, rng, bins: np.ndarray | None = None, mode: str | None = None) -> SimpleData:
    """GATE 0's 'no biology, same technology' reference.

    Counts: every gene replaced by a depth-matched draw (`_depth_matched_matrix`, depth = the
    cell's total). A shuffle within depth bins leaves the depth variation inside each bin, so
    genes coupled only through cell size looked coupled above the null (the same flaw let a
    depth-only pair pass as GATE 5's positive control). Non-count input: genes permuted within
    depth bins.
    """
    if (mode or _null_mode(sd)) == "depth-matched":
        X = np.asarray(sd.X, float)
        return SimpleData(_depth_matched_matrix(X, X.sum(axis=1), rng), sd.obs, sd.var_names)
    return _shuffle_null(sd, bins if bins is not None else _depth_bins(sd), rng)


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


def _shuffle_within(col: np.ndarray, groups: list, rng) -> np.ndarray:
    out = col.copy()
    for idx in groups:
        if len(idx) > 1:
            out[idx] = col[rng.permutation(idx)]
    return out


def _control_response(pair_metric, sub, pair, dose: float, n_rep: int, rng) -> list:
    """The metric's own response to a coupling of `dose` injected into the control pair with its
    own coupling removed: gene b replaced by a depth-matched draw before every injection
    (``injected_deltas(decouple=...)``), injected against the sham, as GATE 4 measures it."""
    from .injected_signal import coupling

    metric = (lambda d: pair_metric(d, gene_a=pair[0], gene_b=pair[1]))
    return injected_deltas(metric, sub, coupling(str(pair[0]), str(pair[1]), strength=dose), n_rep, rng,
                           decouple=str(pair[1]))


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
    delta_min: float | None = None,
    direction: str = "increase",
    n_rep: int = GATE4_N_REP,
) -> GateResult:
    """The positive control must fire in every stratum, judged against an **empirical null**; the
    negative control is reported (v0.4 change 13: the check of specificity in the report takes its
    place, ``specificity``).

    * Negative control: compared with `n_null` distinct unrelated pairs from the same expression
      neighbourhoods (the 20 genes closest in mean, or 5% of genes), two-sided p against alpha/K
      (K = assessable strata, Bonferroni), and the null centre: a centre far from zero means the
      metric reports association between unrelated genes. Reported, not a FAIL since v0.4: an
      unrelated pair is exchangeable with the null's pairs, so one pair fails a valid metric in about
      alpha of the claims by chance alone (the pilot's claim 1b, p = 0.0498, decided 2026-10-17).
      The legacy fixed band (explicit `pos_min` / `neg_max`) still fails it.
    * Positive control: compared with `n_null` self-nulls of the *same* pair, gene b replaced
      by a depth-matched draw from neighbouring cells (thinned to each cell's depth). Coupling
      is destroyed, the dependence on depth kept exactly. It must stand out (p < alpha/K).
    * A silent positive control is judged as GATE 4 judges a response (decided 2026-10-08,
      validation/probes/JOURNAL.md D6, probe p15; the rule of 2026-10-17 since v0.4): a coupling of
      dose `pos_dose` (pre-registered ``positive_control_dose``, default 2.0) is injected `n_rep`
      times into the control pair with its own coupling removed, and the metric's own response gets
      its two-sided (1 - alpha/K) interval (``response_interval``). The metric is shown blind to the
      coupling when the whole interval lies at or below 0 in the declared `direction`
      (``judge_response``); otherwise the silence says that the control is not coupled here, or
      nothing. `delta_min` is accepted and not used.

    p values are rank-based Monte Carlo p values. When a control lies in the extreme tail of
    the first `n_null` draws and alpha/K is below their resolution, the null is extended to
    ``2K/alpha`` draws (at most 5000; ``stats.extend_null``) instead of assuming a normal tail.

    Status:

    * FAIL when the positive control is silent in a stratum where the metric is shown blind to an
      injected coupling of its genes (and, with the legacy band only, when the negative control is
      outside it).
    * WARN when the positive control is silent somewhere without that proof: absence of evidence.
      If it fires nowhere, ``detail["pos_demonstrated"]`` is False and the metric's response stays
      undemonstrated (UNTESTED) unless an injected signal shows it.
    * PASS when the positive control fires in every stratum.

    The former rule (2026-10-07) failed a silent control wherever a reference detector had
    power >= 0.8 for the dose: it failed a valid metric whose control was coupled, but weakly
    (p15). Passing explicit `pos_min` / `neg_max` selects the legacy fixed band (not the
    default), where a silent positive control FAILs without an injection.
    """
    import itertools

    obs = data.obs
    within = list(within)
    delta_min = _check_delta_min(delta_min)
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
            try:
                deltas = _control_response(pair_metric, _as_simple(sub), pos_pair, pos_dose, n_rep, rng)
            except ValueError as e:  # the injection needs counts
                deltas, row["pos_response_note"] = [], str(e)
            if len(deltas) >= 3:
                iv = response_interval(deltas, direction, alpha_s)
                outcome = judge_response(iv, delta_min)
                row.update(pos_response=iv["mean_response"], pos_response_ci=iv["ci"],
                           pos_response_reps=iv["n_rep"], pos_response_outcome=outcome)
            else:
                outcome = "UNTESTED"
                row.update(pos_response_outcome=outcome)
            row["pos_blind"] = outcome == "FAIL"
        row["ok"] = bool(row["pos_fires"] and row["neg_ok"])
        rows.append(row)
        if not row["ok"]:
            failures.append(row)

    n_pos = sum(r["pos_fires"] for r in rows)
    detail = dict(rows=rows, alpha=alpha, alpha_per_stratum=alpha_s, n_null=n_null,
                  legacy_band=legacy, pos_min=pos_min, neg_max=neg_max,
                  pos_demonstrated=bool(n_pos > 0), n_pos_fires=int(n_pos),
                  pos_dose=pos_dose, delta_min=delta_min, direction=direction, n_rep=n_rep)
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
    blind = [r for r in pos_miss if r.get("pos_blind")]
    detail["negative_control_outside_null"] = [r["stratum"] for r in neg_bad]
    neg_note = ""
    if neg_bad and not legacy:
        f = neg_bad[0]
        neg_note = (f"; the negative control {f['neg']:.3g} lies outside its null in {len(neg_bad)}/{len(rows)} strata "
                    f"(e.g. {f['stratum']}: p={f['p_neg']:.3g}, null centre {f['null_center']:.3g}) — reported, not a "
                    "failure: the check of specificity in the report judges the claim instead")
    if blind or (legacy and (pos_miss or neg_bad)):
        f0 = (blind or neg_bad or pos_miss)[0]
        why = []
        if legacy and not f0["neg_ok"]:
            why.append(f"negative control {f0['neg']:.3g} outside the band")
        if f0.get("pos_blind"):
            lo, hi = f0["pos_response_ci"]
            why.append(f"positive control {f0['pos']:.3g} inside its null (p={f0['p_pos']:.3g}) and the "
                       f"metric is blind to a coupling of dose {pos_dose:g} injected into its genes: "
                       f"response {f0['pos_response']:+.4g} ({100 * (1 - alpha_s):g}% CI [{lo:+.4g}, "
                       f"{hi:+.4g}]), the whole interval at or below 0 in the declared direction")
        if legacy and not f0["pos_fires"]:
            why.append(f"positive control {f0['pos']:.3g} below the band")
        n_bad = (len(failures) if legacy else len(blind))
        return GateResult(
            5, "Controls", GateStatus.FAIL,
            f"{n_bad}/{len(rows)} strata fail controls (e.g. {f0['stratum']}: "
            + "; ".join(why) + ")" + neg_note,
            detail,
        )
    if pos_miss:
        def _where(r):
            resp = ""
            if "pos_response" in r:
                lo, hi = r["pos_response_ci"]
                what = ("the metric responds to it, so the control is not coupled here"
                        if r["pos_response_outcome"] == "PASS" else "neither a response nor blindness shown")
                resp = (f"; a coupling of dose {pos_dose:g} injected into its genes: {r['pos_response']:+.3g} "
                        f"[{lo:+.3g}, {hi:+.3g}], {what}")
            return f"{r['stratum']} (n={r['n_cells']}, p={r['p_pos']:.3g}{resp})"
        where = ", ".join(_where(r) for r in pos_miss[:3])
        lead = (f"positive control beats its null in {n_pos}/{len(rows)} strata; not demonstrated in "
                if n_pos else "positive control not demonstrated in any stratum — ")
        return GateResult(
            5, "Controls", GateStatus.WARN,
            lead + where + ("…" if len(pos_miss) > 3 else "")
            + " (the metric is not shown blind to an injected coupling of the control's genes: "
            f"absence of evidence, not a failure){centre_note}{neg_note}",
            detail,
        )
    return GateResult(
        5, "Controls", GateStatus.PASS,
        f"the positive control beats its empirical null in all {len(rows)} strata" + centre_note + neg_note,
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
    sesoi_scale: str = "construct",
) -> GateResult:
    """Re-estimate the effect on an independent dataset with the same estimand, correction,
    replicate rule, strata and SESOI scale. REPLICATED (PASS) = effect detected with the primary sign;
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
        qc_imbalanced=g1.status == GateStatus.WARN, attenuation_half=attenuation_half,
        sesoi_scale=sesoi_scale, seed=seed,
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
