"""The gates as metric-agnostic functions.

Each ``gateN_*`` takes your metric (or the raw data) as a black box and returns a
``GateResult``. Run them in order; a result whose ``.blocking`` is True halts the sequence.

Auto gates:  0 (independence), 1 (QC parity), 2 (n_genes matching), 3 (raw visibility),
             5 (controls), 6 (replication, if a 2nd dataset is supplied).
Judgment gates 4 and 7 are not computable from data alone; the skill elicits them and
`report.autopsy` records them. See references/gates.md for the full rationale.
"""
from __future__ import annotations

from typing import Callable, Sequence

import numpy as np
import pandas as pd

from .core import GateResult, GateStatus, SimpleData, as_dense
from . import qc as _qc

Metric = Callable[[object], float]
# PairMetric: fn(data, *, gene_a, gene_b) -> float  (used by GATES 3 & 5)


def _as_simple(data) -> SimpleData:
    """Materialize any duck-typed data object as a mutable SimpleData copy."""
    if isinstance(data, SimpleData):
        return SimpleData(data.X.copy(), data.obs.copy(), data.var_names)
    return SimpleData(as_dense(data.X).copy(), data.obs.copy(), list(data.var_names))


def _safe_call(fn):
    """Evaluate a metric, returning None (not raising) on a perturbation it can't handle."""
    try:
        v = fn()
    except (KeyError, ValueError, TypeError, ZeroDivisionError):
        return None
    return v if (v is not None and np.isfinite(v)) else None


# --------------------------------------------------------------------------- #
# GATE 0 — mathematical independence
# --------------------------------------------------------------------------- #
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
    elif kind == "variance_inflation":
        g = rng.choice(X.shape[1], size=max(1, X.shape[1] // 10), replace=False)
        mu = X[:, g].mean(axis=0, keepdims=True)
        X[:, g] = mu + (X[:, g] - mu) * 3.0
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
) -> GateResult:
    """Does the metric's EXPECTATION move when a nuisance statistic changes but biology does not?

    Separates confounding (a real shift of the metric's expectation) from estimator noise:
    a bootstrap baseline gives the metric's own sampling spread, and each nuisance is
    flagged only if the perturbed mean is both (a) meaningfully large relative to the metric's
    scale (rel > `tol`, with the denominator floored at `abs_tol` so a near-zero baseline does
    not explode) and (b) statistically separated from the baseline (z > `z_thresh`). This
    stops noisy-but-unbiased metrics and near-null pairs from wrongly failing.
    """
    sd = _as_simple(data)
    rng = np.random.default_rng(seed)
    protect = frozenset(protect_genes)
    n = sd.n_obs

    # Bootstrap baseline: resample cells with replacement (biology preserved) to get the
    # metric's own sampling mean/spread.
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
    base_vals = np.array(base_vals)
    base_mean = float(base_vals.mean())
    base_sd = float(base_vals.std(ddof=1)) if len(base_vals) > 1 else 0.0
    scale = max(abs(base_mean), base_sd, abs_tol)

    kinds = ["extra_dropout", "depth_downsample", "library_scale"]
    if include_matrix_perturbations:
        kinds += ["variance_inflation", "gene_subsample"]

    responses = {}
    for kind in kinds:
        vals, note = [], {}
        for _ in range(n_perturb):
            pdata, note = _perturb(sd, kind, rng, protect)
            v = _safe_call(lambda: metric(pdata))
            if v is not None:
                vals.append(v)
        if len(vals) < 3:
            responses[kind] = dict(skipped=True, rel_change=0.0, z=0.0,
                                   confounded=False, **note)
            continue
        vals = np.array(vals)
        pert_mean = float(vals.mean())
        pert_sd = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
        shift = abs(pert_mean - base_mean)
        rel = shift / scale
        se = np.sqrt(base_sd ** 2 / len(base_vals) + pert_sd ** 2 / len(vals))
        z = shift / se if se > 0 else (float("inf") if shift > 0 else 0.0)
        confounded = bool(rel > tol and z > z_thresh)
        responses[kind] = dict(mean_perturbed=pert_mean, rel_change=float(rel),
                               z=float(z), shift=float(shift), confounded=confounded, **note)

    flagged = {k: r for k, r in responses.items() if r.get("confounded")}
    detail = dict(baseline=base_mean, baseline_sd=base_sd, scale=scale,
                  responses=responses, tol=tol, z_thresh=z_thresh)
    if flagged:
        worst_kind = max(flagged, key=lambda k: flagged[k]["rel_change"])
        w = flagged[worst_kind]
        return GateResult(
            0, "Mathematical independence", GateStatus.FAIL,
            f"metric's expectation shifts {w['rel_change']:.0%} (z={w['z']:.1f}) under "
            f"'{worst_kind}' — confounded by that nuisance",
            detail,
        )
    return GateResult(
        0, "Mathematical independence", GateStatus.PASS,
        "expectation stable under all nuisance perturbations (no shift beyond estimator noise)",
        detail,
    )


# --------------------------------------------------------------------------- #
# GATE 1 — QC parity across factorial strata
# --------------------------------------------------------------------------- #
def gate1_qc_parity(
    data,
    group_col: str,
    groups: tuple,
    within: Sequence[str] = (),
    thresh: float = 1.5,
) -> GateResult:
    """Do the compared groups have equivalent QC in *every* stratum of `within`?"""
    res = _qc.qc_parity(data, group_col, groups, within=within, thresh=thresh)
    if res["n_compared"] == 0:
        return GateResult(
            1, "QC parity", GateStatus.STOP,
            f"no stratum of {list(within)} contains both groups {groups} — groups are "
            "confounded with the stratifier and cannot be compared; do not proceed",
            res,
        )
    flagged = res["flagged"]
    if flagged:
        ratio_breaches = [r for r in flagged if r["ratio_breach"]]
        overlap_breaches = [r for r in flagged if r["overlap_breach"]]
        parts = []
        if ratio_breaches:
            wr = max(ratio_breaches, key=lambda r: r["n_genes_ratio"])
            parts.append(f"{len(ratio_breaches)} exceed {thresh}x QC ratio "
                         f"(worst {wr['n_genes_ratio']:.2f}x at {wr['stratum']})")
        if overlap_breaches:
            wo = min(overlap_breaches, key=lambda r: r["overlap"])
            parts.append(f"{len(overlap_breaches)} have n_genes overlap < {res['overlap_min']} "
                         f"(worst {wo['overlap']:.2f} at {wo['stratum']})")
        return GateResult(
            1, "QC parity", GateStatus.FAIL,
            f"{len(flagged)}/{res['n_compared']} strata fail QC parity: " + "; ".join(parts)
            + " — unusable without n_genes matching (GATE 2)",
            res,
        )
    return GateResult(
        1, "QC parity", GateStatus.PASS,
        f"all {res['n_compared']} strata within {thresh}x QC (worst {res['worst_ratio']:.2f}x)",
        res,
    )


# --------------------------------------------------------------------------- #
# GATE 2 — n_genes matching preserves the signal
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
    """Does the between-group effect survive equalizing cell quality (n_genes)?

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
# GATE 5 — controls behave in all strata
# --------------------------------------------------------------------------- #
def gate5_controls(
    pair_metric: Callable,
    data,
    pos_pair: tuple,
    neg_pair: tuple,
    within: Sequence[str] = (),
    pos_min: float | None = None,
    neg_max: float | None = None,
    min_cells: int = 10,
) -> GateResult:
    """Positive control must fire and negative control must stay null — in every stratum.

    `pair_metric(data, *, gene_a, gene_b) -> float`. The null band is anchored to the
    POSITIVE control's magnitude (not the negative control's own value, which is the very
    quantity under test): by default the negative must be < 20% of the positive signal, and
    the positive must clear that band. Magnitudes are used so signed metrics with strong
    negative-going controls still pass.
    """
    import itertools

    obs = data.obs
    within = list(within)
    if within:
        levels = [sorted(pd.unique(obs[f].dropna())) for f in within]
        combos = list(itertools.product(*levels))
    else:
        combos = [()]

    def _val(sub, pair):
        return float(pair_metric(sub, gene_a=pair[0], gene_b=pair[1]))

    pos_pooled = abs(_val(data, pos_pair))
    if neg_max is None:
        neg_max = 0.2 * pos_pooled + 1e-6           # independent of the neg control itself
    if pos_min is None:
        pos_min = neg_max                           # positive must clear the null band

    rows, failures = [], []
    for combo in combos:
        mask = np.ones(len(obs), dtype=bool)
        for f, v in zip(within, combo):
            mask &= (np.asarray(obs[f]) == v)
        if mask.sum() < min_cells:
            continue
        sub = data[mask]
        pos, neg = _val(sub, pos_pair), _val(sub, neg_pair)
        ok = abs(pos) > pos_min and abs(neg) <= neg_max
        row = dict(stratum=dict(zip(within, combo)) if within else "pooled",
                   pos=pos, neg=neg, ok=ok)
        rows.append(row)
        if not ok:
            failures.append(row)

    detail = dict(rows=rows, pos_min=pos_min, neg_max=neg_max, pos_pooled=pos_pooled)
    if not rows:
        return GateResult(
            5, "Controls", GateStatus.SKIP,
            f"no stratum had >= {min_cells} cells (checked {len(combos)}) — controls not evaluable",
            detail,
        )
    if failures:
        f0 = failures[0]
        return GateResult(
            5, "Controls", GateStatus.FAIL,
            f"{len(failures)}/{len(rows)} strata fail controls "
            f"(e.g. {f0['stratum']}: pos={f0['pos']:.3g}, neg={f0['neg']:.3g}; "
            f"need |pos|>{pos_min:.3g}, |neg|<={neg_max:.3g})",
            detail,
        )
    return GateResult(
        5, "Controls", GateStatus.PASS,
        f"positive fires and negative null in all {len(rows)} strata",
        detail,
    )


# --------------------------------------------------------------------------- #
# GATE 6 — cross-platform / cross-species replication
# --------------------------------------------------------------------------- #
def gate6_replication(
    metric: Metric,
    data2,
    group_col: str,
    groups: tuple,
    thresh: float = 1.5,
) -> GateResult:
    """Re-run QC parity + n_genes matching on an independent dataset."""
    if data2 is None:
        return GateResult(
            6, "Replication", GateStatus.SKIP,
            "no second dataset supplied — supply an independent AnnData to run this gate",
            {},
        )
    g1 = gate1_qc_parity(data2, group_col, groups, thresh=thresh)
    g2 = gate2_ngenes_matching(metric, data2, group_col, groups)
    replicated = g2.status == GateStatus.PASS
    status = GateStatus.PASS if replicated else GateStatus.FAIL
    return GateResult(
        6, "Replication", status,
        f"independent dataset: QC parity {g1.status.value}, matched effect {g2.status.value}",
        dict(qc_parity=g1.detail, matching=g2.detail),
    )
