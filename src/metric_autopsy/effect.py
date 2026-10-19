"""The effect field: a between-group difference, corrected for the declared estimand,
inferred at the level of biological replicates.

Steps (``estimate_effect``):

1. **Estimand-dependent correction** (``equalize``). For a *composition* estimand the
   deeper group is thinned to equal depth. For a *content* estimand it is thinned to equal
   spike-in capture, or the comparison is UNIDENTIFIABLE when there is no external
   standard. With no declared estimand the comparison is UNIDENTIFIABLE too.
2. **Replicates.** The metric is computed per replicate (mouse, donor, plate…), nested in
   groups (and strata) or paired within replicates. The unit of inference is the
   replicate, never the cell.
3. **Graded rule** (minimum replicates per group, or pairs):
   * >= 4: an exact permutation over replicates (stratified; Monte Carlo with >= 1000 draws
     and an MC error when enumeration is too large) is the test; the t interval estimates
     the effect;
   * 3: the t interval on replicate-level values is the test, and the result is flagged
     PARAMETRIC_ONLY; the permutation p is reported with its minimum attainable value,
     because alpha cannot be reached by construction;
   * <= 2 (or fewer than the pre-registered minimum, or no replicate unit declared):
     no effect verdict; design adequacy = INSUFFICIENT_REPLICATION.
4. **Decision.** DETECTED, or NO_DETECTABLE_EFFECT (TOST: the (1-2 alpha) interval lies
   inside ±SESOI), or INCONCLUSIVE. ``explained_by_depth`` marks a raw difference that is
   detected at the replicate level and disappears after the correction (corrected effect
   not detected, less than half retained, and - with a SESOI - shown smaller than the SESOI
   by the TOST, journal D7). ``reversed_by_correction`` marks a detected raw
   difference whose sign the correction reverses, also detected: INCONCLUSIVE, because the
   direction then depends on how exactly a technical difference larger than the effect was
   removed. A raw difference that is not itself detected
   cannot be "explained": the effect is then INCONCLUSIVE, and the reason reports how much
   of it the correction left. (Testing raw minus corrected per replicate instead would be
   anti-conservative: the untouched group's component is zero by construction.)
   A correction can only remove a conclusion, never create one (v0.4 change 6): DETECTED needs
   the effect detected, in one direction, by the raw replicate-level test (the same data without
   the correction) and by the corrected one; ``detected_only_after_correction`` marks the rest,
   INCONCLUSIVE. The TOST runs on the uncorrected data too, and an absence that only the corrected
   data show is flagged ABSENCE_ONLY_AFTER_CORRECTION (reported, not blocking).
5. **Power.** The SESOI is declared on the construct scale. GATE 0's attenuation under
   depth halving gives a reliability-model estimate of how much of a construct-scale
   difference survives at the analysed depth (lambda). The design is UNDERPOWERED when the
   minimum detectable effect exceeds lambda x SESOI. The same attenuated SESOI is used for TOST.
   A SESOI pre-registered on the metric's observed scale (``sesoi_scale="observed"``) is used as
   it is (lambda = 1).
"""
from __future__ import annotations

from typing import Sequence

import numpy as np

from .core import Assessment, SimpleData
from .equalize import detection_shares, equalize_detection, looks_like_counts, spikein_columns, thin_to_match
from .gates import _as_simple, _safe_call
from . import qc as _qc
from .stats import mde, nested_t, paired_t, permutation_test_nested, permutation_test_paired, t_ppf

THINNING = ("depth_thinning", "capture_thinning")
FLOOR_POWER = 0.95  # change 4: the TOST's power at a true difference of 0, fixed in advance (v1's rule)


def design_floor(se: float, df: float, alpha: float, power: float = FLOOR_POWER) -> float:
    """The design's floor (v0.4 change 4): the smallest SESOI at which this design confirms absence by
    the TOST with probability `power` when the true difference is 0, in the metric's units:
    se x (t(1 - alpha, df) + t((1 + power) / 2, df))."""
    if not (np.isfinite(se) and se > 0 and df > 0):
        return float("nan")
    return float(se * (t_ppf(1.0 - alpha, df) + t_ppf((1.0 + power) / 2.0, df)))


def equalizing_correction(data, estimand: str | None, spikein_prefix: str = "ERCC-") -> str | None:
    """The thinning the estimand implies, if any: depth thinning for a composition estimand and
    capture thinning for a content estimand with spike-ins, both on raw counts. Either removes
    a between-group depth (capture) difference, so GATE 0 does not block a depth bias it
    removes."""
    if not looks_like_counts(data.X):
        return None
    if estimand == "composition":
        return "depth_thinning"
    if estimand == "content" and spikein_columns(data.var_names, spikein_prefix):
        return "capture_thinning"
    return None


SESOI_SCALES = ("construct", "observed")


def attenuation_lambda(a_half: float | None, depth_ratio: float = 1.0,
                       sesoi_scale: str = "construct") -> tuple[float, dict]:
    """Fraction of a construct-scale difference retained at the analysed depth.

    Reliability model: measurement noise variance ∝ 1/depth, so reliability at depth d is
    1 / (1 + k/d). If halving depth removes a fraction ``a`` of the signal, then
    k/d = a / (1 - 2a). Equalization further scales the thinned group's depth by
    ``depth_ratio``. Returns 1 when no attenuation was measured (nothing to correct), and when
    the SESOI is pre-registered on the observed scale (``sesoi_scale="observed"``: it already
    is a difference of the metric as measured at this depth).
    """
    if sesoi_scale == "observed":
        return 1.0, dict(model="observed scale (pre-registered sesoi_scale): no attenuation applied")
    if a_half is None or not np.isfinite(a_half) or a_half <= 0:
        return 1.0, dict(model="none (no attenuation measured)")
    if a_half >= 0.5:
        return 0.0, dict(model="reliability", a_half=a_half,
                         note="halving depth removed >= 50% of the signal: noise-dominated")
    x = a_half / (1.0 - 2.0 * a_half)
    s = depth_ratio if (depth_ratio and depth_ratio > 0) else 1.0
    lam = 1.0 / (1.0 + x / s)
    return float(lam), dict(model="reliability: noise ∝ 1/depth", a_half=a_half,
                            k_over_depth=x, depth_ratio=s)


def _cell_effect(metric, data, group_col, groups):
    g = np.asarray(data.obs[group_col])
    a, b = g == groups[0], g == groups[1]
    if not a.any() or not b.any():
        return None
    va = _safe_call(lambda: metric(data[a]))
    vb = _safe_call(lambda: metric(data[b]))
    return None if va is None or vb is None else float(va - vb)


def _replicate_design(obs, group_col, groups, replicate_col, within):
    """Classify replicates as nested in groups or paired within replicates; check strata."""
    rep = np.asarray(obs[replicate_col]).astype(str)
    grp = np.asarray(obs[group_col])
    reps = list(dict.fromkeys(rep))
    membership = {r: set(grp[rep == r]) & set(groups) for r in reps}
    if all(len(m) == 1 for m in membership.values()):
        kind = "nested"
    elif all(len(m) == 2 for m in membership.values()):
        kind = "paired"
    else:
        return dict(kind="mixed", reps=reps, reason=(
            "replicates are partially crossed with groups (some contain both groups, some one) "
            "— neither a nested nor a paired analysis is valid"))
    strata = {}
    if within:
        keys = obs[list(within)].astype(str).agg("\x1f".join, axis=1).to_numpy()
        for r in reps:
            ks = set(keys[rep == r])
            if len(ks) != 1:
                return dict(kind="mixed", reps=reps, reason=(
                    f"replicate {r!r} spans several strata of {list(within)}"))
            strata[r] = next(iter(ks))
    else:
        strata = {r: "" for r in reps}
    return dict(kind=kind, reps=reps, membership=membership, strata=strata)


_UNIT_SEPARATORS = "/+,;&|"


def composite_units(reps) -> list:
    """Units whose id names several animals (parts joined by one of ``/ + , ; & |``) and the other
    units they overlap with (decided 2026-10-17, v0.4 change 16: on preprint section 4's setup two
    units held two mice each, 3_10_M/3_11_M and 3_38_F/3_39_F, whose parts were units of their own).
    A check of the ids' text, not of the cells: it warns, it does not exclude."""
    import re
    split = re.compile("\\s*[" + re.escape(_UNIT_SEPARATORS) + "]\\s*")
    ids = [str(r) for r in reps]
    parts = {r: {x for x in split.split(r) if x} or {r} for r in ids}
    out = []
    for r in ids:
        if len(parts[r]) >= 2:
            out.append(dict(unit=r, names=sorted(parts[r]),
                            overlaps=sorted(u for u in ids if u != r and parts[u] & parts[r])))
    return out


def _tier(n_eff: int, min_replicates: int) -> str:
    threshold = max(3, int(min_replicates or 3))
    if n_eff <= 2 or n_eff < threshold:
        return "insufficient"
    return "parametric" if n_eff == 3 else "permutation"


def _rep_values(metric, data, group_col, groups, replicate_col, design):
    rep = np.asarray(data.obs[replicate_col]).astype(str)
    grp = np.asarray(data.obs[group_col])
    out = {}
    for r in design["reps"]:
        m = rep == r
        if design["kind"] == "nested":
            out[r] = _safe_call(lambda: metric(data[m]))
        else:
            va = _safe_call(lambda: metric(data[m & (grp == groups[0])]))
            vb = _safe_call(lambda: metric(data[m & (grp == groups[1])]))
            out[r] = None if va is None or vb is None else va - vb
    return out


def _n_eff(reps, design: dict, groups) -> int:
    """The replicates per group (the smaller group) or the pairs among `reps`."""
    if design["kind"] == "paired":
        return len(reps)
    first = sum(1 for r in reps if groups[0] in design["membership"][r])
    return min(first, len(reps) - first)


def _measured(raw_vals: dict, corr_vals: dict, design: dict, groups) -> tuple[dict, dict, list, int, int]:
    """The replicates where the metric gives a value on the uncorrected and on the corrected data. The
    others are left out of both tests and listed (v0.4, the owner's decision of 2026-10-18: a value that
    was not measured is not a value of 0): (raw, corrected, the unmeasured, replicates per group before
    and after)."""
    def fin(v):
        return v is not None and np.isfinite(v)
    reps = list(design["reps"])
    ok = {r: fin(raw_vals.get(r)) and fin(corr_vals.get(r)) for r in reps}
    unmeasured = [dict(replicate=r, raw=raw_vals.get(r) if fin(raw_vals.get(r)) else None,
                       corrected=corr_vals.get(r) if fin(corr_vals.get(r)) else None) for r in reps if not ok[r]]
    return ({r: raw_vals[r] if ok[r] else None for r in reps}, {r: corr_vals[r] if ok[r] else None for r in reps},
            unmeasured, _n_eff(reps, design, groups), _n_eff([r for r in reps if ok[r]], design, groups))


def _unmeasured_reason(n_left: int, n_total: int, n_after: int, min_replicates: int, what: str) -> str:
    return (f"the metric gives no value (not measured) on {n_left} of the {n_total} {what}, which leaves "
            f"{n_after} replicate(s) per group, fewer than {min_replicates}: the comparison cannot be made on "
            "these data (the unmeasured are listed)")


def _infer(values: dict, design: dict, groups, tier: str, alpha: float, n_perm: int, rng,
           max_exact: int = 20000, min_replicates: int = 3) -> dict:
    """Run the graded test on replicate values ({replicate: value}); NaNs dropped."""
    reps = [r for r in design["reps"] if values.get(r) is not None and np.isfinite(values[r])]
    notes = []
    if design["kind"] == "paired":
        d = np.array([values[r] for r in reps], float)
        te = paired_t(d)
        perm = permutation_test_paired(d, n_perm=n_perm, rng=rng, max_exact=max_exact)
        n_eff = len(d)
        n_counts = dict(pairs=n_eff)
    else:
        g = np.array([1 if groups[0] in design["membership"][r] else 0 for r in reps])
        labels = {s: i for i, s in enumerate(dict.fromkeys(design["strata"][r] for r in reps))}
        s = np.array([labels[design["strata"][r]] for r in reps])
        y = np.array([values[r] for r in reps], float)
        te = nested_t(y, g, s)
        perm = permutation_test_nested(y, g, s, n_perm=n_perm, rng=rng, max_exact=max_exact)
        n_counts = {str(groups[0]): int(g.sum()), str(groups[1]): int((1 - g).sum())}
        n_eff = int(min(g.sum(), (1 - g).sum()))
    return _tiered(te, perm, n_eff, n_counts, tier, alpha, min_replicates, notes)


def _tiered(te, perm, n_eff, n_counts, tier, alpha, min_replicates, notes) -> dict:
    """The graded rule on a t estimate and a relabeling test: which of them is the test."""
    if tier != "insufficient":
        tier = _tier(n_eff, min_replicates)  # replicates whose metric failed are dropped
        if tier == "insufficient":
            notes.append(f"only {n_eff} replicate(s) per group had a finite metric value")
    if tier == "permutation" and not (perm.min_attainable_p < alpha):
        notes.append(f"permutation p = {perm.p:.3g} reported only: with {n_eff} replicates per "
                     f"group under this design alpha={alpha} is unattainable (min attainable p = "
                     f"{perm.min_attainable_p:.3g}), so the t interval is the test (parametric only)")
        tier = "parametric"
    elif tier == "parametric" and perm.min_attainable_p >= alpha:
        notes.append(f"permutation p = {perm.p:.3g} reported only: alpha={alpha} is unattainable "
                     f"by construction (min attainable p = {perm.min_attainable_p:.3g})")
    if tier == "permutation":
        detected = bool(perm.p < alpha)
        test = f"{perm.method} permutation over replicates (n={perm.n})"
    elif tier == "parametric":
        detected = bool(np.isfinite(te.p) and te.p < alpha)
        test = f"t interval on replicate values ({te.method}; parametric only)"
    else:
        detected, test = False, "none (insufficient replication)"
    return dict(est=te.est, se=te.se, df=te.df, t_method=te.method, p_t=te.p,
                ci=te.ci(1 - alpha), ci_2alpha=te.ci(1 - 2 * alpha), perm=perm, tier=tier,
                detected=detected, test=test, notes=notes, n_counts=n_counts, n_eff=n_eff,
                te=te)


def _units_inference(units, values, groups, tier, alpha, n_perm, rng, max_exact=20000, min_replicates=3) -> dict:
    """The graded test on the domain's combined estimate (``combined``): the same rule as on
    replicate values, with the estimate's relabeling test and its cluster-robust t interval."""
    from .combined import Estimator, t_estimate
    est = Estimator(units)
    V = np.asarray(values, float)[None, :]
    D, se, df = est.estimate(V)
    te = t_estimate(D[0], se[0], df[0], units.kind)
    perm = est.permutation(V, n_perm=n_perm, rng=rng, max_exact=max_exact)[0]
    c = est.counts()
    if units.kind == "paired":
        n_counts, n_eff = dict(pairs=c["pairs"]), c["pairs"]
    else:
        n_counts = {str(groups[0]): c["first"], str(groups[1]): c["second"]}
        n_eff = min(c["first"], c["second"])
    return _tiered(te, perm, n_eff, n_counts, tier, alpha, min_replicates, [])


def _restrict(data, group_col, groups):
    sd_all = _as_simple(data)
    gcol = np.asarray(sd_all.obs[group_col])
    return sd_all[(gcol == groups[0]) | (gcol == groups[1])]


def assess_design(
    data,
    *,
    group_col: str,
    groups: tuple,
    within: Sequence[str] = (),
    replicate_col: str | None = None,
    estimand: str | None = None,
    min_replicates: int = 3,
    spikein_prefix: str = "ERCC-",
    qc_imbalanced: bool = False,
    _sd: SimpleData | None = None,
) -> dict:
    """Design adequacy that does not depend on the metric: is the estimand's correction
    possible on these data, and how many biological replicates does each group have?

    Returns a dict with ``identifiable`` / ``reason``, ``correction`` (+ ``totals`` for
    thinning), ``replication`` and the graded ``tier`` (permutation | parametric |
    insufficient) with ``tier_reason``.
    """
    within = list(within)
    sd = _sd if _sd is not None else _restrict(data, group_col, groups)
    design = dict(identifiable=True, reason="", correction="none", estimand=estimand,
                  replicate_col=replicate_col, notes=[], totals=None)

    # ---- estimand-dependent correction -------------------------------------------------
    if estimand == "composition":
        if looks_like_counts(sd.X):
            design["totals"] = np.asarray(_qc.per_cell_qc(sd)["total_counts"], dtype=float)
            design["correction"] = "depth_thinning"
        elif qc_imbalanced:
            design.update(identifiable=False, reason=(
                "non-count input: a depth difference between the groups cannot be removed by "
                "thinning (supply raw counts)"))
        else:
            design["notes"].append("non-count input; QC balanced, no correction applied")
    elif estimand == "content":
        cols = spikein_columns(sd.var_names, spikein_prefix)
        if not cols:
            design.update(identifiable=False, reason=(
                f"a content estimand (total RNA, genes detected) needs an external standard; no "
                f"spike-ins with prefix {spikein_prefix!r} were found, so capture cannot be "
                "separated from content"))
        elif not looks_like_counts(sd.X):
            design.update(identifiable=False, reason="spike-in capture correction needs raw counts")
        else:
            totals = sd.X[:, cols].sum(axis=1)
            frac_with = float(np.mean(totals > 0))
            design["spikein_columns"] = len(cols)
            design["frac_cells_with_spikeins"] = frac_with
            if frac_with < 0.9:
                design.update(identifiable=False, reason=(
                    f"spike-ins detected in only {frac_with:.0%} of cells — capture is not "
                    "estimable per cell"))
            else:
                design.update(correction="capture_thinning", totals=totals)
    elif estimand is None:
        design.update(identifiable=False, reason=(
            "no estimand pre-registered: declare 'composition' (relative expression) or "
            "'content' (total RNA / genes detected) — the correction depends on it"))
    else:
        design.update(identifiable=False, reason=f"unknown estimand {estimand!r}")

    # ---- replicates ---------------------------------------------------------------------
    design["rep_design"] = None
    if replicate_col is None:
        design.update(tier="insufficient", replication=dict(kind="undeclared"), tier_reason=(
            "no replicate unit declared (replicate_col): a cell-level comparison treats cells "
            "as independent, so no effect verdict is given (estimates only)"))
        return design
    rep_design = _replicate_design(sd.obs, group_col, groups, replicate_col, within)
    if rep_design["kind"] == "mixed":
        design.update(identifiable=False, reason=rep_design["reason"], tier="insufficient",
                      replication=dict(kind="mixed"))
        return design
    if rep_design["kind"] == "paired":
        n_eff = len(rep_design["reps"])
        counts = dict(pairs=n_eff)
    else:
        na = sum(1 for r in rep_design["reps"] if groups[0] in rep_design["membership"][r])
        nb = len(rep_design["reps"]) - na
        n_eff, counts = min(na, nb), {str(groups[0]): na, str(groups[1]): nb}
    tier = _tier(n_eff, min_replicates)
    design.update(rep_design=rep_design, tier=tier,
                  replication=dict(kind=rep_design["kind"], counts=counts, n_eff=n_eff))
    composite = composite_units(rep_design["reps"])
    if composite:
        design["composite_units"] = composite
        design["notes"].append(
            "units whose ids name several animals: " + "; ".join(
                f"{c['unit']}" + (f" (overlaps {', '.join(c['overlaps'])})" if c["overlaps"] else "")
                for c in composite) + " — each counts as one replicate; check that it is one")
    if tier == "insufficient":
        design["tier_reason"] = (f"{n_eff} replicate(s) per group (need >= "
                                 f"{max(3, int(min_replicates or 3))}): no effect verdict")
    return design


def capture_check(raw: SimpleData, equalized: SimpleData, group_col: str, groups, replicate_col=None) -> dict:
    """The data against the capture model, for the report (not a gate; v0.4 change 6, decided).

    Under the model behind the thinning the groups differ by capture alone, so at equal depth their
    cells detect as many genes. Per group: the mean of the replicates' (else the cells') mean depth
    and genes detected per cell, before and after the thinning, and the groups' relative difference.
    A difference in genes detected that the thinning leaves is either biology or an artifact the
    capture model does not describe (entry-level dropout, say); the check cannot tell which."""
    out = {}
    for name, d in (("raw", raw), ("equalized", equalized)):
        X = d.X if isinstance(d.X, np.ndarray) else np.asarray(d.X)
        depth = X.sum(axis=1)
        genes = (X > 0).sum(axis=1)
        grp = np.asarray(d.obs[group_col])
        rep = np.asarray(d.obs[replicate_col]).astype(str) if replicate_col else None
        row = {}
        for g in groups:
            m = grp == g
            if rep is None:
                row[str(g)] = dict(depth=float(depth[m].mean()) if m.any() else float("nan"),
                                   genes=float(genes[m].mean()) if m.any() else float("nan"), unit="cells")
                continue
            units = [r for r in dict.fromkeys(rep[m])]
            row[str(g)] = dict(depth=float(np.mean([depth[m & (rep == r)].mean() for r in units])) if units else float("nan"),
                               genes=float(np.mean([genes[m & (rep == r)].mean() for r in units])) if units else float("nan"),
                               unit="replicates")
        a, b = row[str(groups[0])], row[str(groups[1])]
        row["relative_difference"] = {k: (a[k] / b[k] - 1.0 if b[k] else float("nan")) for k in ("depth", "genes")}
        out[name] = row
    return out


def public_design(design: dict) -> dict:
    """The design dict without internal arrays/objects (for reports)."""
    return {k: v for k, v in design.items() if k not in ("totals", "rep_design")}


def _totals(data, correction: str, spikein_prefix: str) -> np.ndarray:
    if correction == "capture_thinning":
        cols = spikein_columns(data.var_names, spikein_prefix)
        return np.asarray(data.X[:, cols].sum(axis=1), dtype=float)
    return np.asarray(_qc.per_cell_qc(data)["total_counts"], dtype=float)


def _stratum_keys(obs, within) -> np.ndarray:
    if within:
        return obs[list(within)].astype(str).agg("\x1f".join, axis=1).to_numpy()
    return np.full(len(obs), "", dtype=object)


def _summary(inf: dict | None, alpha: float) -> dict | None:
    """An inference's numbers for the report."""
    if inf is None:
        return None
    p = inf["perm"].p if inf["tier"] == "permutation" else inf["p_t"]
    return dict(effect=inf["est"], ci=inf["ci"], ci_level=1 - alpha, se=inf["se"], df=inf["df"], p=p,
                test=inf["test"], tier=inf["tier"], detected=inf["detected"], n_replicates=inf["n_counts"])


def _mean_over_sets(sets) -> np.ndarray:
    """Per unit, the mean of its finite values over the equalized draws (NaN where none is finite)."""
    stack = np.array(sets)
    ok = np.isfinite(stack)
    return np.where(ok.any(axis=0), np.where(ok, stack, 0.0).sum(axis=0) / np.maximum(ok.sum(axis=0), 1), np.nan)


def _mean_over_dicts(per_set: list, reps) -> dict:
    """Per replicate, the mean of its values over the equalized draws (None where it has none)."""
    out = {}
    for r in reps:
        xs = [p[r] for p in per_set if p.get(r) is not None]
        out[r] = float(np.mean(xs)) if xs else None
    return out


def _domain_effect(metric, sd: SimpleData, dom: dict, *, group_col, groups, within, replicate_col, correction,
                   tier, alpha, n_perm, n_equalize, rng, max_exact, min_replicates, spikein_prefix,
                   claim_genes=(), det_rng=None, dossier_rng=None, reference=None) -> dict:
    """The claim within its domain (types equally weighted), raw and corrected, and beside it the types
    outside the domain and every type on its own. The correction thins within each stratum and type; for a
    pair (`claim_genes`, change 17') each draw at equal depth then has the claim genes' detection equalized
    within each stratum and type (``equalize_detection``, from `det_rng`), and the same estimate at equal
    depth alone is kept for the dossier. `reference` (the engine's ``norm_pearson``) is computed on the same
    data, raw and corrected, for the line beside the metric."""
    from .domain import build_units, unit_values

    lab = np.asarray(sd.obs[dom["column"]]).astype(str)
    keys = _stratum_keys(sd.obs, within)
    qualifying = {k: sorted(set(dom["domain"].get(k, [])) | set(dom["outside"].get(k, []))) for k in dom["domain"]}
    mask = np.zeros(sd.n_obs, bool)
    for k, ts in qualifying.items():
        if ts:
            mask |= (keys == k) & np.isin(lab, ts)
    if not mask.any():
        return dict(empty=True)
    typed = sd[mask]
    units = build_units(dom, typed, group_col=group_col, groups=groups, replicate_col=replicate_col, which="all")
    if units is None:
        return dict(empty=True)
    rep_key = {r: dom["_design"]["strata"][r] for r in units.reps}
    in_dom = np.array([units.types[t] in dom["domain"].get(rep_key[units.reps[r]], [])
                       for r, t in zip(units.unit_rep, units.unit_type)], bool)
    if not in_dom.any():
        return dict(empty=True)
    raw = unit_values(metric, typed, units)
    ref_raw = unit_values(reference, typed, units) if reference is not None else None
    eq_info, cc, det_info, corr_depth, ref_corr = None, None, None, None, None
    cells_by = list(within) + [dom["column"]]
    if correction in THINNING:
        totals = _totals(typed, correction, spikein_prefix)
        sets, sets_depth, ref_sets = [], [], []
        for i in range(max(1, n_equalize)):
            eq, info = thin_to_match(typed, group_col, groups, totals=totals, within=cells_by, rng=rng)
            eq_info = eq_info or info
            if i == 0:
                cc = capture_check(typed, eq, group_col, groups, replicate_col)
            if claim_genes:
                sets_depth.append(unit_values(metric, eq, units))
                det = equalize_detection(eq, claim_genes, group_col, groups, within=cells_by, rng=det_rng)
                det_info = det_info or det
            sets.append(unit_values(metric, eq, units))
            if reference is not None:
                ref_sets.append(unit_values(reference, eq, units))
            del eq
        corr = _mean_over_sets(sets)
        corr_depth = _mean_over_sets(sets_depth) if sets_depth else None
        ref_corr = _mean_over_sets(ref_sets) if ref_sets else None
    else:
        corr = raw
        ref_corr = ref_raw
    # the units where the metric gives no value, uncorrected or at equal depth, are left out and listed;
    # an Estimator counts a replicate by its units, so they must go before it is built
    measured = np.isfinite(raw) & np.isfinite(corr)
    unit_table = [dict(replicate=str(units.reps[r]), type=units.types[t], in_domain=bool(d),
                       raw=float(v) if np.isfinite(v) else None, corrected=float(c) if np.isfinite(c) else None,
                       measured=bool(m))
                  for r, t, d, v, c, m in zip(units.unit_rep, units.unit_type, in_dom, raw, corr, measured)]
    unmeasured = [u for u in unit_table if not u["measured"]]
    n_before = _units_n_eff(units.subset(in_dom))
    keep = measured
    units, raw, corr, in_dom = units.subset(keep), raw[keep], corr[keep], in_dom[keep]
    corr_depth = corr_depth[keep] if corr_depth is not None else None
    ref_raw = ref_raw[keep] if ref_raw is not None else None
    ref_corr = ref_corr[keep] if ref_corr is not None else None
    n_after = _units_n_eff(units.subset(in_dom)) if in_dom.any() else 0
    detection = (dict(genes=list(claim_genes), cells="each stratum and cell type of the domain",
                      records=_detection_records(typed, claim_genes, group_col, groups, cells_by, det_info,
                                                 [t for ts in dom["domain"].values() for t in ts], dom["column"]))
                 if det_info is not None else None)
    if unmeasured and (not in_dom.any() or n_after < min_replicates <= n_before):
        return dict(empty=False, unmeasured=unmeasured, unmeasured_insufficient=True, n_before=n_before,
                    n_after=n_after, n_units=int(len(unit_table)), n_cells=int(typed.n_obs), units=unit_table,
                    types=sorted({u["type"] for u in unit_table if u["in_domain"]}), eq_info=eq_info,
                    capture_check=cc, outside=None, per_type=[], detection=detection)
    dom_units = units.subset(in_dom)
    inf = _units_inference(dom_units, corr[in_dom], groups, tier, alpha, n_perm, rng, max_exact, min_replicates)
    inf_raw = _units_inference(dom_units, raw[in_dom], groups, tier, alpha, n_perm, rng, max_exact, min_replicates)
    inf_depth = (_units_inference(dom_units, corr_depth[in_dom], groups, tier, alpha, n_perm, dossier_rng, max_exact,
                                  min_replicates) if corr_depth is not None else None)
    outside = None
    if (~in_dom).any():
        out_units = units.subset(~in_dom)
        outside = dict(types=sorted({units.types[t] for t in units.unit_type[~in_dom]}),
                       corrected=_summary(_units_inference(out_units, corr[~in_dom], groups, tier, alpha, n_perm, rng,
                                                           max_exact, min_replicates), alpha),
                       raw=_summary(_units_inference(out_units, raw[~in_dom], groups, tier, alpha, n_perm, rng,
                                                     max_exact, min_replicates), alpha))
    per_type = []
    det = {r["type"]: r for s in dom["strata"] for r in s["types"]}
    for ti, t in enumerate(units.types):
        sel = units.unit_type == ti
        if not sel.any():
            continue
        sub = units.subset(sel)
        rec = dict(type=t, in_domain=bool(in_dom[sel].any()), units=int(sel.sum()),
                   detection=det.get(t, {}).get("detection"), index=det.get(t, {}).get("index"),
                   corrected=_summary(_units_inference(sub, corr[sel], groups, tier, alpha, n_perm, rng, max_exact,
                                                       min_replicates), alpha),
                   raw=_summary(_units_inference(sub, raw[sel], groups, tier, alpha, n_perm, rng, max_exact,
                                                 min_replicates), alpha))
        per_type.append(rec)
    ref = (dict(kind="units", units=dom_units, metric=corr[in_dom], reference=ref_corr[in_dom],
                raw_metric=raw[in_dom], raw_reference=ref_raw[in_dom], corrected=correction in THINNING)
           if reference is not None else None)
    return dict(empty=False, inf=inf, inf_raw=inf_raw, inf_depth=inf_depth, outside=outside, per_type=per_type,
                eq_info=eq_info, capture_check=cc, n_units=int(in_dom.sum()), n_cells=int(typed.n_obs),
                units=unit_table, types=sorted({units.types[t] for t in units.unit_type[in_dom]}),
                unmeasured=unmeasured, unmeasured_insufficient=False, n_before=n_before, n_after=n_after,
                detection=detection, _reference=ref)


def _detection_records(raw_data, genes, group_col, groups, cells_by, det_info, keep_types=None, type_col=None) -> list:
    """The dossier's detection of the claim's genes per cell set: each group's share of cells with a count,
    uncorrected, at equal depth and after the equalization (the first draw's), with the thinned group and
    its keep probability; with `keep_types`, only the cell sets of those types."""
    raw = {(r["cells"] if isinstance(r["cells"], str) else tuple(sorted(r["cells"].items())), r["gene"]): r["shares"]
           for r in detection_shares(raw_data, genes, group_col, groups, within=cells_by)}
    out = []
    for r in det_info or []:
        cells = r["cells"]
        if keep_types is not None and type_col is not None and isinstance(cells, dict) and cells.get(type_col) not in keep_types:
            continue
        key = (cells if isinstance(cells, str) else tuple(sorted(cells.items())), r["gene"])
        out.append(dict(cells=cells, gene=r["gene"], raw=raw.get(key), equal_depth=r["before"], after=r["after"],
                        thinned_group=r["thinned_group"], keep_p=r["keep_p"], n_cells=r["n_cells"]))
    return out


def equalized_for_probe(data, *, group_col: str, groups: tuple, within: Sequence[str] = (), type_col: str | None = None,
                        estimand: str | None = None, spikein_prefix: str = "ERCC-", claim_genes: Sequence[str] = (),
                        seed: int = 0) -> tuple:
    """The data GATE 4 measures its response on (the owner's decision of 2026-10-19: the SESOI in the metric's
    units is measured on the data the effect is judged on, at equal depth and, for a pair, at equal detection).
    One draw of the effect's correction: the comparison's cells thinned to equal depth (or spike-in capture)
    within each stratum, and each cell type with `type_col` (the domain's), then for a pair the claim genes'
    detection equalized in the same cell sets (``equalize.equalize_detection``), from streams keyed by
    `seed`. Without a thinning correction (no estimand that implies one, or no counts) the data as they
    are. Returns (data, a record of what was done)."""
    sd = _restrict(data, group_col, groups)
    correction = equalizing_correction(sd, estimand, spikein_prefix)
    if correction not in THINNING:
        return data, dict(correction=None, data="as they are (no thinning correction)")
    cells_by = list(within) + ([type_col] if type_col else [])
    eq, info = thin_to_match(sd, group_col, groups, totals=_totals(sd, correction, spikein_prefix), within=cells_by,
                             rng=np.random.default_rng([int(seed), 4]))
    genes = [str(g) for g in (claim_genes or ())]
    det = (equalize_detection(eq, genes, group_col, groups, within=cells_by, rng=np.random.default_rng([int(seed), 4, 17]))
           if genes else None)
    what = "equal depth" if correction == "depth_thinning" else "equal spike-in capture"
    return eq, dict(correction=correction, cells_by=cells_by, depth=info, detection=det,
                    data=what + (", the claim genes' detection equalized" if genes else ""))


def _units_n_eff(units) -> int:
    """The contributing replicates per group (the smaller group) or the pairs of a set of units."""
    from .combined import Estimator
    if units.n_units == 0:
        return 0
    c = Estimator(units).counts()
    return c["pairs"] if units.kind == "paired" else min(c["first"], c["second"])


def estimate_effect(
    metric,
    data,
    *,
    group_col: str,
    groups: tuple,
    within: Sequence[str] = (),
    replicate_col: str | None = None,
    estimand: str | None = None,
    sesoi: float | None = None,
    alpha: float = 0.05,
    power: float = 0.8,
    min_replicates: int = 3,
    n_perm: int = 1000,
    n_equalize: int = 5,
    spikein_prefix: str = "ERCC-",
    qc_imbalanced: bool = False,
    attenuation_half: float | None = None,
    sesoi_scale: str = "construct",
    seed: int = 0,
    max_exact: int = 20000,
    domain: dict | None = None,
    sesoi_label: str | None = None,
    claim_genes: Sequence[str] = (),
    reference=None,
) -> tuple[Assessment, dict]:
    """Estimate the between-group effect; return (effect assessment, design info).

    The SESOI: on the construct's scale by default (equivalence within ±λ x SESOI, λ the measured
    attenuation), or on the metric's scale with ``sesoi_scale="observed"`` (``run_autopsy`` passes
    GATE 4's response s_m so, v0.4 change 3, with `sesoi_label` naming it); without one the design's
    floor (``design_floor``, v0.4 change 4) is the margin of absence, and the report says that absence
    is confirmed only for effects of at least the floor. With a `domain` (``domain.claim_domain``) the
    effect is the claim's combined estimate within the domain's cell types, types weighted equally, and
    the whole mixture's estimate goes into the detail beside it (v0.4, decided 2026-10-18).

    `claim_genes` (a claim about a pair, a metric of structure: v0.4 change 17', the owner's decision of
    2026-10-19): after equal depth, each draw has the claim genes' detection equalized
    (``equalize.equalize_detection``) in every cell set where the effect is estimated (each stratum; with a
    domain, each stratum and type), from a random stream of its own keyed by `seed`, so the thinning to
    equal depth draws as before; the corrected effect, and so the verdict, is at equal depth and detection.
    The same estimate at equal depth alone (``detail["equal_depth_effect"]``) and each group's detection of
    each gene, uncorrected, at equal depth and after (``detail["detection"]``), go to the dossier. Without a
    thinning correction (no counts) nothing is equalized. `reference` (the engine's ``norm_pearson`` for
    the line beside the metric) is computed on the same data, raw and corrected
    (``detail["_reference"]``, for ``report.coupling_reference``)."""
    if sesoi_scale not in SESOI_SCALES:
        raise ValueError(f"sesoi_scale must be one of {SESOI_SCALES}, not {sesoi_scale!r}")
    rng = np.random.default_rng(seed)
    claim_genes = [str(g) for g in (claim_genes or ())]
    # change 17': the detection's thinning and the dossier's extra tests draw from streams of their own
    det_rng = np.random.default_rng([int(seed), 17]) if claim_genes else None
    dossier_rng = np.random.default_rng([int(seed), 17, 1]) if claim_genes else None
    within = list(within)
    sd = _restrict(data, group_col, groups)
    design = assess_design(sd, group_col=group_col, groups=groups, within=within,
                           replicate_col=replicate_col, estimand=estimand,
                           min_replicates=min_replicates, spikein_prefix=spikein_prefix,
                           qc_imbalanced=qc_imbalanced, _sd=sd)
    totals, rep_design, tier = design["totals"], design["rep_design"], design["tier"]

    cell_raw = _cell_effect(metric, sd, group_col, groups)
    eq_info, cc, per_set, cell_corr_vals = None, None, [], []
    det_info, per_set_depth, ref_set = None, [], []
    equalizes = design["identifiable"] and design["correction"] in THINNING
    if equalizes:
        for i in range(max(1, n_equalize)):
            eq, info = thin_to_match(sd, group_col, groups, totals=totals, within=within, rng=rng)
            eq_info = eq_info or info
            if i == 0:
                cc = capture_check(sd, eq, group_col, groups, replicate_col)
            if claim_genes:
                if rep_design is not None:
                    per_set_depth.append(_rep_values(metric, eq, group_col, groups, replicate_col, rep_design))
                det = equalize_detection(eq, claim_genes, group_col, groups, within=within, rng=det_rng)
                det_info = det_info or det
            v = _cell_effect(metric, eq, group_col, groups)
            if v is not None:
                cell_corr_vals.append(v)
            if rep_design is not None:
                per_set.append(_rep_values(metric, eq, group_col, groups, replicate_col, rep_design))
                if reference is not None:
                    ref_set.append(_rep_values(reference, eq, group_col, groups, replicate_col, rep_design))
            del eq
    else:
        if cell_raw is not None:
            cell_corr_vals = [cell_raw]
    depth_ratio = eq_info["depth_ratio"] if eq_info else 1.0
    if eq_info is not None:
        eq_info["capture_check"] = cc
    design["equalization"] = eq_info
    cell_corr = float(np.mean(cell_corr_vals)) if cell_corr_vals else None

    detail = dict(cell_level_raw=cell_raw, cell_level_corrected=cell_corr,
                  correction=design["correction"], estimand=estimand, sesoi=sesoi, alpha=alpha)

    if not design["identifiable"]:
        return Assessment("NOT_ESTIMABLE", design["reason"], [], detail), public_design(design)

    if rep_design is None:
        detail.update(effect=cell_corr, raw_effect=cell_raw)
        return Assessment("INCONCLUSIVE", design["tier_reason"], [], detail), public_design(design)

    # ---- 3. replicate values and the graded test ----------------------------------------
    raw_vals = _rep_values(metric, sd, group_col, groups, replicate_col, rep_design)
    ref_raw = _rep_values(reference, sd, group_col, groups, replicate_col, rep_design) if reference is not None else None
    if not per_set:
        per_set = [raw_vals]
    corr_vals = _mean_over_dicts(per_set, rep_design["reps"])
    corr_depth = _mean_over_dicts(per_set_depth, rep_design["reps"]) if per_set_depth else None
    ref_corr = _mean_over_dicts(ref_set, rep_design["reps"]) if ref_set else ref_raw
    # a replicate where the metric gives no value (not measured) leaves both tests and is listed
    raw_vals, corr_vals, unmeasured, n_before, n_after = _measured(raw_vals, corr_vals, rep_design, groups)
    kept = {r for r, v in corr_vals.items() if v is not None}
    if corr_depth is not None:
        corr_depth = {r: (v if r in kept else None) for r, v in corr_depth.items()}

    inf = _infer(corr_vals, rep_design, groups, tier, alpha, n_perm, rng, max_exact, min_replicates)
    inf_raw = _infer(raw_vals, rep_design, groups, tier, alpha, n_perm, rng, max_exact, min_replicates)
    inf_depth = (_infer(corr_depth, rep_design, groups, tier, alpha, n_perm, dossier_rng, max_exact, min_replicates)
                 if corr_depth is not None else None)
    mixture = None
    detail["unmeasured"] = unmeasured
    detection = (dict(genes=claim_genes, cells="each stratum" if within else "all the comparison's cells",
                      records=_detection_records(sd, claim_genes, group_col, groups, within, det_info))
                 if det_info is not None else None)
    reference_values = (dict(kind="replicates", design=rep_design, metric=corr_vals, reference=ref_corr,
                             raw_metric=raw_vals, raw_reference=ref_raw, corrected=equalizes)
                        if reference is not None else None)
    if domain is not None and not domain.get("error"):
        mixture = dict(corrected=_summary(inf, alpha), raw=_summary(inf_raw, alpha),
                       replicate_values=[dict(replicate=r, raw=raw_vals.get(r), corrected=corr_vals.get(r))
                                         for r in rep_design["reps"]],
                       equalization=eq_info, unmeasured=unmeasured)
        if claim_genes:
            mixture.update(equal_depth=_summary(inf_depth, alpha), detection=detection)
        dres = _domain_effect(metric, sd, domain, group_col=group_col, groups=groups, within=within,
                              replicate_col=replicate_col, correction=design["correction"], tier=tier, alpha=alpha,
                              n_perm=n_perm, n_equalize=n_equalize, rng=rng, max_exact=max_exact,
                              min_replicates=min_replicates, spikein_prefix=spikein_prefix,
                              claim_genes=claim_genes, det_rng=det_rng, dossier_rng=dossier_rng, reference=reference)
        detail["mixture"] = mixture
        detail["domain"] = {k: v for k, v in dres.items() if k not in ("inf", "inf_raw", "inf_depth", "eq_info",
                                                                         "detection", "_reference")}
        detail["domain"]["rule"] = domain["rule"]
        if dres["empty"]:
            detail.update(domain_empty=True, effect=None, raw_effect=None)
            reason = ("the claim's domain is empty: no cell type with enough cells and replicates has "
                      + ("the claim's genes detected" if domain["rule"]["kind"] != "global" else "both groups")
                      + f" by the rule ({domain['rule']['text']}); the whole mixture's estimate is reported beside")
            return Assessment("INCONCLUSIVE", reason, [], detail), public_design(design)
        if dres.get("unmeasured_insufficient"):
            detail.update(unmeasured_insufficient=True, effect=None, raw_effect=None)
            n_dom = sum(1 for u in dres["units"] if u["in_domain"])
            n_off = sum(1 for u in dres["unmeasured"] if u["in_domain"])
            return Assessment("NOT_ESTIMABLE", _unmeasured_reason(
                n_off, n_dom, dres["n_after"], min_replicates, "units (replicate, type) of the domain"),
                [], detail), public_design(design)
        inf, inf_raw, inf_depth = dres["inf"], dres["inf_raw"], dres["inf_depth"]
        detection, reference_values = dres["detection"], dres["_reference"]
        if dres["eq_info"] is not None:
            eq_info = dict(dres["eq_info"], capture_check=dres["capture_check"])
            depth_ratio = eq_info["depth_ratio"]
            design["equalization"] = eq_info
    elif unmeasured and n_after < min_replicates <= n_before:
        detail.update(unmeasured_insufficient=True, effect=None, raw_effect=None)
        return Assessment("NOT_ESTIMABLE", _unmeasured_reason(
            len(unmeasured), len(rep_design["reps"]), n_after, min_replicates, "replicates"),
            [], detail), public_design(design)
    if tier != "insufficient" and inf["tier"] == "insufficient":
        design["tier_reason"] = "; ".join(inf["notes"]) or "insufficient replication"
    tier = inf["tier"] if tier != "insufficient" else tier
    design["tier"] = tier
    design["notes"] += inf["notes"]

    detail.update(
        effect=inf["est"], raw_effect=inf_raw["est"], se=inf["se"], df=inf["df"],
        t_method=inf["t_method"], ci=inf["ci"], ci_level=1 - alpha, ci_2alpha=inf["ci_2alpha"],
        p_t=inf["p_t"], p_perm=inf["perm"].p, perm_method=inf["perm"].method,
        perm_n=inf["perm"].n, perm_mc_se=inf["perm"].mc_se,
        min_attainable_p=inf["perm"].min_attainable_p, test=inf["test"], tier=tier,
        n_replicates=inf["n_counts"], design_kind=rep_design["kind"],
        raw_detected=inf_raw["detected"],
        p_raw=(inf_raw["perm"].p if inf_raw["tier"] == "permutation" else inf_raw["p_t"]),
        replicate_values=[dict(replicate=r, raw=raw_vals.get(r), corrected=corr_vals.get(r),
                               stratum=(rep_design["strata"][r].split("\x1f")
                                        if rep_design["strata"][r] else None),
                               group=(str(groups[0]) if groups[0] in rep_design["membership"][r]
                                      else str(groups[1])) if rep_design["kind"] == "nested" else "paired")
                          for r in rep_design["reps"]],
        estimate_over=("the domain's cell types, equally weighted" if mixture is not None else "the replicates"),
    )
    if claim_genes:
        detail.update(detection=detection, equal_depth_effect=_summary(inf_depth, alpha),
                      corrected_for=("depth and the claim genes' detection" if equalizes else "nothing (no counts)"))
    if reference_values is not None:
        detail["_reference"] = reference_values
    raw_est, est = inf_raw["est"], inf["est"]
    retained = est / raw_est if (raw_est and np.isfinite(raw_est) and raw_est != 0) else float("nan")
    detail["retained"] = float(retained) if np.isfinite(retained) else float("nan")

    # ---- 4. the SESOI in the metric's units, the floor, power ------------------------------
    floor = design_floor(inf["se"], inf["df"], alpha)
    floor_raw = design_floor(inf_raw["se"], inf_raw["df"], alpha)
    lam, lam_info = attenuation_lambda(attenuation_half, depth_ratio, sesoi_scale)
    lam_raw, _ = attenuation_lambda(attenuation_half, 1.0, sesoi_scale)
    if sesoi is not None:
        eff_sesoi, eff_sesoi_raw, source = lam * sesoi, lam_raw * sesoi, "declared"
    else:
        eff_sesoi = floor if np.isfinite(floor) else None
        eff_sesoi_raw = floor_raw if np.isfinite(floor_raw) else None
        source = "floor"
    mde_v = mde(inf["se"], inf["df"], alpha, power) if tier != "insufficient" else float("nan")
    design["power"] = dict(mde=mde_v, power=power, sesoi=sesoi, lambda_=lam, effective_sesoi=eff_sesoi,
                           attenuation=lam_info, floor=floor, floor_power=FLOOR_POWER, sesoi_source=source,
                           sesoi_label=sesoi_label)
    design["underpowered"] = (bool(np.isfinite(mde_v) and mde_v > eff_sesoi)
                              if sesoi is not None and eff_sesoi is not None else None)
    design["sesoi_below_floor"] = (bool(np.isfinite(floor) and eff_sesoi < floor)
                                   if sesoi is not None and eff_sesoi is not None else None)
    detail.update(effective_sesoi=eff_sesoi, sesoi_source=source, floor=floor)

    flags = ["PARAMETRIC_ONLY"] if tier == "parametric" else []
    if tier == "insufficient":
        return Assessment("INCONCLUSIVE", design["tier_reason"], flags, detail), public_design(design)

    lo2, hi2 = inf["ci_2alpha"]
    equivalent = (eff_sesoi is not None and eff_sesoi > 0 and np.isfinite(lo2) and np.isfinite(hi2)
                  and lo2 > -eff_sesoi and hi2 < eff_sesoi)
    detail["tost_equivalent"] = bool(equivalent) if eff_sesoi is not None else None
    # "explained by depth" says the raw difference is technical: with a SESOI that also needs the
    # corrected effect shown smaller than the SESOI (TOST), else a real effect the test missed would
    # be called depth (decided 2026-10-08, journal D7); without a SESOI only the retained share speaks
    depth_like = bool(design["correction"] in THINNING and inf_raw["detected"] and not inf["detected"]
                      and np.isfinite(retained) and abs(retained) < 0.5)
    explained = depth_like and (sesoi is None or bool(equivalent))
    detail["explained_by_depth"] = explained
    # a detected raw difference whose sign the correction reverses: the estimate then hinges on
    # how exactly the correction removed a technical difference larger than the effect
    reversed_ = bool(design["correction"] in THINNING and inf_raw["detected"] and inf["detected"]
                     and np.isfinite(raw_est) and np.isfinite(est) and np.sign(est) != np.sign(raw_est))
    detail["reversed_by_correction"] = reversed_
    # a correction can only remove a conclusion, never create one (v0.4 change 6): the corrected
    # effect counts as detected only where the raw test detects it too, in the same direction
    created = bool(design["correction"] in THINNING and inf["detected"] and not reversed_
                   and not (inf_raw["detected"] and np.isfinite(raw_est) and np.sign(raw_est) == np.sign(est)))
    detail["detected_only_after_correction"] = created
    # the TOST on the uncorrected data, at their own depth (decided 2026-10-12): an absence that
    # only the corrected data show is marked, not blocked
    lo2r, hi2r = inf_raw["ci_2alpha"]
    equivalent_raw = (eff_sesoi_raw is not None and eff_sesoi_raw > 0 and np.isfinite(lo2r) and np.isfinite(hi2r)
                      and lo2r > -eff_sesoi_raw and hi2r < eff_sesoi_raw)
    detail.update(tost_equivalent_raw=bool(equivalent_raw) if eff_sesoi_raw is not None else None,
                  effective_sesoi_raw=eff_sesoi_raw, ci_2alpha_raw=inf_raw["ci_2alpha"])

    ci = inf["ci"]
    p_txt = (f"p={inf['perm'].p:.3g} ({inf['perm'].method} permutation)" if tier == "permutation"
             else f"p={inf['p_t']:.3g} (t, parametric only)")
    over = " within the domain's types, equally weighted" if mixture is not None else ""
    est_txt = f"{est:+.4g} [{ci[0]:+.4g}, {ci[1]:+.4g}]{over}"
    if explained:
        what = "depth" if design["correction"] == "depth_thinning" else "capture"
        return Assessment(
            "INCONCLUSIVE",
            f"the raw difference {raw_est:+.4g} is explained by {what}: at equal {what} it is "
            f"{est_txt} ({retained:.0%} retained), {p_txt}",
            flags, detail), public_design(design)
    if reversed_:
        what = "depth" if design["correction"] == "depth_thinning" else "capture"
        return Assessment(
            "INCONCLUSIVE",
            f"the correction reversed the sign: raw {raw_est:+.4g} (detected) -> {est_txt} at equal "
            f"{what}, {p_txt}; the technical difference is larger than the effect, so its direction "
            "depends on how exactly the correction removed it", flags, detail), public_design(design)
    if created:
        what = "depth" if design["correction"] == "depth_thinning" else "capture"
        return Assessment(
            "INCONCLUSIVE",
            f"effect {est_txt} at equal {what}, {p_txt}, but the raw difference {raw_est:+.4g} is not "
            f"detected at the replicate level in that direction (p={detail['p_raw']:.3g}): a correction "
            "can only remove a conclusion, never create one", flags, detail), public_design(design)
    if inf["detected"]:
        return Assessment("DETECTED", f"effect {est_txt}, {p_txt}", flags, detail), public_design(design)
    if equivalent:
        if source == "floor":
            margin = (f"±{eff_sesoi:.4g}, the design's floor: no SESOI was pre-registered, so absence is "
                      f"confirmed only for effects of at least {eff_sesoi:.4g}")
        elif sesoi_label:
            margin = f"±{eff_sesoi:.4g} ({sesoi_label}; TOST)"
        elif sesoi_scale == "observed":
            margin = f"±{eff_sesoi:.4g} (SESOI {sesoi:g} on the metric's scale; TOST)"
        else:
            margin = f"±{eff_sesoi:.4g} (SESOI {sesoi:g} x attenuation {lam:.2f}; TOST)"
        reason = f"effect {est_txt}: the {1 - 2 * alpha:.0%} interval lies within {margin}"
        if source == "floor":
            flags = flags + ["NO_SESOI_FLOOR"]
        if design["correction"] in THINNING and not equivalent_raw:
            what = "depth" if design["correction"] == "depth_thinning" else "capture"
            flags = flags + ["ABSENCE_ONLY_AFTER_CORRECTION"]
            reason += (f"; absence shown only after the {what} correction: on the uncorrected data the "
                       f"{1 - 2 * alpha:.0%} interval [{lo2r:+.4g}, {hi2r:+.4g}] is not within "
                       f"±{eff_sesoi_raw:.4g}")
        return Assessment("NO_DETECTABLE_EFFECT", reason, flags, detail), public_design(design)
    why = (f"neither detected nor equivalent within ±{eff_sesoi:.4g}" + (
        " (the design's floor; no SESOI was pre-registered)" if source == "floor" else "")
           if eff_sesoi is not None else "no SESOI pre-registered and no floor computable, so equivalence cannot be tested")
    if design["correction"] in THINNING and np.isfinite(retained) and abs(retained) < 0.5:
        what = "depth" if design["correction"] == "depth_thinning" else "capture"
        if inf_raw["detected"]:
            why += (f"; the raw difference {raw_est:+.4g} (detected, p={detail['p_raw']:.3g}) shrinks "
                    f"to {retained:.0%} at equal {what}, but the corrected effect is not shown smaller "
                    f"than the SESOI, so {what} does not explain it")
        else:
            why += (f"; the raw difference {raw_est:+.4g} (not detected at the replicate level, "
                    f"p={detail['p_raw']:.3g}) shrinks to {retained:.0%} at equal {what}")
    return Assessment("INCONCLUSIVE", f"effect {est_txt}, {p_txt}; {why}", flags, detail), public_design(design)
