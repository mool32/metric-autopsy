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
   not detected, less than half retained). ``reversed_by_correction`` marks a detected raw
   difference whose sign the correction reverses, also detected: INCONCLUSIVE, because the
   direction then depends on how exactly a technical difference larger than the effect was
   removed. A raw difference that is not itself detected
   cannot be "explained": the effect is then INCONCLUSIVE, and the reason reports how much
   of it the correction left. (Testing raw minus corrected per replicate instead would be
   anti-conservative: the untouched group's component is zero by construction.)
5. **Power.** The SESOI is declared on the construct scale. GATE 0's attenuation under
   depth halving gives a reliability-model estimate of how much of a construct-scale
   difference survives at the analysed depth (lambda). The design is UNDERPOWERED when the
   minimum detectable effect exceeds lambda x SESOI. The same attenuated SESOI is used for TOST.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np

from .core import Assessment, SimpleData
from .equalize import looks_like_counts, spikein_columns, thin_to_match
from .gates import _as_simple, _safe_call
from . import qc as _qc
from .stats import mde, nested_t, paired_t, permutation_test_nested, permutation_test_paired

THINNING = ("depth_thinning", "capture_thinning")


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


def attenuation_lambda(a_half: float | None, depth_ratio: float = 1.0) -> tuple[float, dict]:
    """Fraction of a construct-scale difference retained at the analysed depth.

    Reliability model: measurement noise variance ∝ 1/depth, so reliability at depth d is
    1 / (1 + k/d). If halving depth removes a fraction ``a`` of the signal, then
    k/d = a / (1 - 2a). Equalization further scales the thinned group's depth by
    ``depth_ratio``. Returns 1 when no attenuation was measured (nothing to correct).
    """
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
    if tier == "insufficient":
        design["tier_reason"] = (f"{n_eff} replicate(s) per group (need >= "
                                 f"{max(3, int(min_replicates or 3))}): no effect verdict")
    return design


def public_design(design: dict) -> dict:
    """The design dict without internal arrays/objects (for reports)."""
    return {k: v for k, v in design.items() if k not in ("totals", "rep_design")}


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
    seed: int = 0,
    max_exact: int = 20000,
) -> tuple[Assessment, dict]:
    """Estimate the between-group effect; return (effect assessment, design info)."""
    rng = np.random.default_rng(seed)
    within = list(within)
    sd = _restrict(data, group_col, groups)
    design = assess_design(sd, group_col=group_col, groups=groups, within=within,
                           replicate_col=replicate_col, estimand=estimand,
                           min_replicates=min_replicates, spikein_prefix=spikein_prefix,
                           qc_imbalanced=qc_imbalanced, _sd=sd)
    totals, rep_design, tier = design["totals"], design["rep_design"], design["tier"]

    eq_sets, eq_info = [sd], None
    if design["identifiable"] and design["correction"] in THINNING:
        eq_sets = []
        for _ in range(max(1, n_equalize)):
            eq, info = thin_to_match(sd, group_col, groups, totals=totals, within=within, rng=rng)
            eq_sets.append(eq)
            eq_info = eq_info or info
    depth_ratio = eq_info["depth_ratio"] if eq_info else 1.0
    design["equalization"] = eq_info

    cell_raw = _cell_effect(metric, sd, group_col, groups)
    cell_corr_vals = [v for v in (_cell_effect(metric, e, group_col, groups) for e in eq_sets)
                      if v is not None]
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
    per_set = [_rep_values(metric, e, group_col, groups, replicate_col, rep_design) for e in eq_sets]
    corr_vals = {}
    for r in rep_design["reps"]:
        xs = [p[r] for p in per_set if p.get(r) is not None]
        corr_vals[r] = float(np.mean(xs)) if xs else None

    inf = _infer(corr_vals, rep_design, groups, tier, alpha, n_perm, rng, max_exact, min_replicates)
    inf_raw = _infer(raw_vals, rep_design, groups, tier, alpha, n_perm, rng, max_exact, min_replicates)
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
    )
    raw_est, est = inf_raw["est"], inf["est"]
    retained = est / raw_est if (raw_est and np.isfinite(raw_est) and raw_est != 0) else float("nan")
    detail["retained"] = float(retained) if np.isfinite(retained) else float("nan")

    # ---- 4. power against the (attenuated) SESOI ----------------------------------------
    lam, lam_info = attenuation_lambda(attenuation_half, depth_ratio)
    eff_sesoi = lam * sesoi if sesoi is not None else None
    mde_v = mde(inf["se"], inf["df"], alpha, power) if tier != "insufficient" else float("nan")
    design["power"] = dict(mde=mde_v, power=power, sesoi=sesoi, lambda_=lam,
                           effective_sesoi=eff_sesoi, attenuation=lam_info)
    design["underpowered"] = (bool(np.isfinite(mde_v) and mde_v > eff_sesoi)
                              if eff_sesoi is not None else None)
    detail["effective_sesoi"] = eff_sesoi

    flags = ["PARAMETRIC_ONLY"] if tier == "parametric" else []
    if tier == "insufficient":
        return Assessment("INCONCLUSIVE", design["tier_reason"], flags, detail), public_design(design)

    lo2, hi2 = inf["ci_2alpha"]
    equivalent = (eff_sesoi is not None and eff_sesoi > 0 and np.isfinite(lo2) and np.isfinite(hi2)
                  and lo2 > -eff_sesoi and hi2 < eff_sesoi)
    detail["tost_equivalent"] = bool(equivalent) if eff_sesoi is not None else None
    explained = bool(design["correction"] in THINNING and inf_raw["detected"] and not inf["detected"]
                     and np.isfinite(retained) and abs(retained) < 0.5)
    detail["explained_by_depth"] = explained
    # a detected raw difference whose sign the correction reverses: the estimate then hinges on
    # how exactly the correction removed a technical difference larger than the effect
    reversed_ = bool(design["correction"] in THINNING and inf_raw["detected"] and inf["detected"]
                     and np.isfinite(raw_est) and np.isfinite(est) and np.sign(est) != np.sign(raw_est))
    detail["reversed_by_correction"] = reversed_

    ci = inf["ci"]
    p_txt = (f"p={inf['perm'].p:.3g} ({inf['perm'].method} permutation)" if tier == "permutation"
             else f"p={inf['p_t']:.3g} (t, parametric only)")
    est_txt = f"{est:+.4g} [{ci[0]:+.4g}, {ci[1]:+.4g}]"
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
    if inf["detected"]:
        return Assessment("DETECTED", f"effect {est_txt}, {p_txt}", flags, detail), public_design(design)
    if equivalent:
        return Assessment(
            "NO_DETECTABLE_EFFECT",
            f"effect {est_txt}: the {1 - 2 * alpha:.0%} interval lies within ±{eff_sesoi:.4g} "
            f"(SESOI {sesoi:g} x attenuation {lam:.2f}; TOST)", flags, detail), public_design(design)
    why = ("no SESOI pre-registered, so equivalence cannot be tested" if sesoi is None
           else f"neither detected nor equivalent within ±{eff_sesoi:.4g}")
    if design["correction"] in THINNING and np.isfinite(retained) and abs(retained) < 0.5:
        what = "depth" if design["correction"] == "depth_thinning" else "capture"
        why += (f"; the raw difference {raw_est:+.4g} (not detected at the replicate level, "
                f"p={detail['p_raw']:.3g}) shrinks to {retained:.0%} at equal {what}")
    return Assessment("INCONCLUSIVE", f"effect {est_txt}, {p_txt}; {why}", flags, detail), public_design(design)
