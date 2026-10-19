"""Assemble a gate-by-gate autopsy and a verdict from four independent fields.

The verdict is no longer "the first blocking gate". It is decided by ``decide`` from four
fields that answer different questions:

* ``metric_validity`` — does the metric respond to its construct and resist nuisance?
  PASS needs demonstrated response (a positive control that beats an empirical null, or
  an injected signal); without one it is UNTESTED. A metric that never varies is
  DEGENERATE. A nuisance that *biases* it (GATE 0) or a failed control (GATE 5) is FAIL.
  Attenuation is not a validity question and is not recorded here.
* ``design_adequacy`` — is the comparison identifiable, and can it detect the SESOI? It
  covers QC balance (GATE 1 is a diagnostic), the estimand-dependent correction,
  replication and power. GATE 0's attenuation enters here, as the power check against the
  attenuated SESOI (a SESOI pre-registered on the observed scale is not attenuated). Statuses:
  ADEQUATE, CORRECTED, INSUFFICIENT_REPLICATION,
  UNIDENTIFIABLE; flags: UNDERPOWERED, PARAMETRIC_ONLY, ATTENUATION, POWER_NOT_ASSESSED.
* ``effect`` — the corrected difference at the replicate level: DETECTED,
  NO_DETECTABLE_EFFECT (TOST against the SESOI in the metric's units: GATE 4's response s_m to a change
  of the pre-registered construct-scale SESOI, v0.4 change 3; without a SESOI the design's floor),
  INCONCLUSIVE, NOT_ESTIMABLE. With cell-type labels the difference is combined within the claim's
  domain, the types weighted equally, and the whole mixture's is reported beside (v0.4).
* ``replication`` — REPLICATED / NOT_REPLICATED / INCONCLUSIVE / NOT_RUN (GATE 6).

No rescue language: a FAIL in an earlier field is never softened by a later one. The verdict follows
the first FAIL; since v0.4 the effect, the domain's, the mixture's and the composition are estimated
after it too and reported as the dossier, marked as not deciding the verdict (the owner's decision
of 2026-10-18). Every difference in the dossier is groups[0] - groups[1], the effect test's convention
(the owner's decision of 2026-10-19). For a metric of structure whose groups differ in detection, the
difference their unequal responses to the same coupling would make at an equal coupling (the
differential sensitivity, decided 2026-10-19) is subtracted at its bound before the effect can support
the claim (``differential_sensitivity``).
Judgment gates (4, 7) are pending unless the pre-registration says otherwise, so a
SUPPORTED verdict is never reached by default. The Python API, the CLI and the MCP server
all call this one function.
"""
from __future__ import annotations

import functools
import json
import os
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

from .core import Assessment, GateResult, GateStatus
from . import gates as _g
from . import provenance as _prov

DEFAULT_PREREG = dict(judgment_pending=True, alpha=0.05, power=0.8, min_replicates=3,
                      signal_direction="increase", spikein_prefix="ERCC-",
                      positive_control_dose=2.0, bias_tolerance=0.5)
# The claim's direction: the change of the metric from groups[0] to groups[1]. Mandatory for
# SUPPORTED; "two-sided" declares a non-directional claim, which the verdict marks as such.
DIRECTIONS = ("increase", "decrease", "two-sided")

NOT_RUN = "NOT_RUN"


def _metric_name(metric) -> str:
    """Unwrap functools.partial so partial-bound metrics keep their function name."""
    f = metric
    while isinstance(f, functools.partial):
        f = f.func
    return getattr(f, "__name__", "metric")


def _cell(s) -> str:
    """Escape a dynamic value for safe interpolation into a markdown table/list."""
    return (str(s).replace("\\", "\\\\").replace("|", "\\|")
            .replace("\r", " ").replace("\n", "<br>"))


def normalize_prereg(prereg: dict | None) -> dict:
    """Pre-registration with defaults filled in. ``judgment_resolved: True`` is accepted as
    the inverse of ``judgment_pending``."""
    p = dict(DEFAULT_PREREG)
    p.update(prereg or {})
    if "judgment_resolved" in (prereg or {}):
        p["judgment_pending"] = not bool(prereg["judgment_resolved"])
    if p.get("direction") is not None and p["direction"] not in DIRECTIONS:
        raise ValueError(f"direction must be one of {DIRECTIONS} (the change from groups[0] to "
                         f"groups[1]), not {p['direction']!r}")
    from .effect import SESOI_SCALES
    if p.get("sesoi_scale", "construct") not in SESOI_SCALES:
        raise ValueError(f"sesoi_scale must be one of {SESOI_SCALES}, not {p['sesoi_scale']!r}")
    return p


def delta_min_of(prereg: dict) -> float | None:
    """v0.3's delta_min, the smallest response to an injected signal that matters: pre-registered as
    ``delta_min``, by default ``DELTA_MIN_FRACTION`` (0.5) x SESOI. Since v0.4 (change 3, decided
    2026-10-17) GATE 4 and GATE 5 decide only whether a metric is blind and ``run_autopsy`` does not use
    it; the function stays for reading v0.3's pre-registrations."""
    if prereg.get("delta_min") is not None:
        return float(prereg["delta_min"])
    if prereg.get("sesoi") is not None:
        return _g.DELTA_MIN_FRACTION * float(prereg["sesoi"])
    return None


def _na(reason: str) -> Assessment:
    return Assessment(NOT_RUN, reason)


# The cause of every verdict, one code per branch of `decide_cause` (in its order). The label
# (the verdict's text before " — ") and the cause together say why a claim is or is not
# supported: e.g. NOT SUPPORTED with metric_invalid_gate4 (the metric is shown blind), with
# explained_by_depth or with opposite_direction.
CAUSES = ("degenerate_metric", "metric_invalid_gate0", "metric_invalid_gate5", "metric_invalid_gate4",
          "metric_invalid",
          "unidentifiable", "unmeasured", "effect_not_evaluated", "empty_domain", "explained_by_depth",
          "insufficient_replication",
          "absence_untested_metric", "no_detectable_effect", "detected_only_after_correction",
          "effect_inconclusive",
          "detected_untested_metric", "bias_unsized", "no_direction", "opposite_direction",
          "differential_sensitivity", "not_replicated", "effect_below_sesoi", "judgment_pending", "replicated",
          "provisional")
# The probes of a metric of structure (a pair; a probe that couples a set would join them): for these the
# groups' sensitivity to the same coupling is compared (the owner's decision of 2026-10-19); a metric of
# level (a module's expression) is not
STRUCTURE_PROBES = ("coupling",)


def decide_cause(a: "Autopsy") -> tuple[str, str]:
    """The single verdict rule (shared by the Python API, the CLI and the MCP server): the
    verdict and its cause (one of ``CAUSES``)."""
    mv = a.metric_validity or _na("not evaluated")
    da = a.design_adequacy or _na("not evaluated")
    ef = a.effect or _na("not evaluated")
    rp = a.replication or _na("no independent dataset")
    notes = []
    if "PARAMETRIC_ONLY" in ef.flags:
        notes.append("parametric only — too few replicates for the permutation test to reach alpha")
    if "UNDERPOWERED" in da.flags:
        notes.append("underpowered relative to the SESOI")
    dom = ef.detail.get("domain") if ef.detail else None
    if dom and not dom.get("empty"):
        notes.append(f"judged within the domain: {', '.join(dom.get('types', []))}, the types weighted equally; "
                     "the whole mixture is reported beside")
    for comp in (getattr(a, "composition", None) or {}).get("columns", {}).values():
        if comp.get("differs"):  # a claim about a mixture only with its difference of composition
            notes.append(f"the groups' composition in {comp['column']} differs (total variation distance "
                         f"{comp['tvd']:.3g}, p {comp['tvd_p']:.3g})")
    tail = f" [{'; '.join(notes)}]" if notes else ""

    if mv.status == "DEGENERATE":
        return f"DEGENERATE METRIC — {mv.reason}", "degenerate_metric"
    if mv.status == "FAIL":
        gate = mv.detail.get("failed_gate")
        return (f"NOT SUPPORTED — metric invalid: {mv.reason}",
                f"metric_invalid_gate{gate}" if gate in (0, 4, 5) else "metric_invalid")
    if da.status == "UNIDENTIFIABLE":
        return f"UNIDENTIFIABLE — {da.reason}", "unidentifiable"
    if ef.detail.get("unmeasured_insufficient"):
        return f"UNIDENTIFIABLE — {ef.reason}", "unmeasured"
    if ef.status == NOT_RUN:
        return f"INCONCLUSIVE — effect not evaluated: {ef.reason}", "effect_not_evaluated"
    if ef.detail.get("domain_empty"):
        return f"INCONCLUSIVE — {ef.reason}", "empty_domain"
    if ef.detail.get("explained_by_depth"):
        return f"NOT SUPPORTED — {ef.reason}{tail}", "explained_by_depth"
    if da.status == "INSUFFICIENT_REPLICATION":
        return f"INCONCLUSIVE — insufficient replication: {da.reason}", "insufficient_replication"
    if ef.status == "NO_DETECTABLE_EFFECT":
        if mv.status != "PASS":
            return ("INCONCLUSIVE — no detectable effect within the SESOI, but the metric's "
                    "response to signal is untested, so absence cannot be claimed" + tail,
                    "absence_untested_metric")
        return f"NO DETECTABLE EFFECT — {ef.reason}{tail}", "no_detectable_effect"
    if ef.detail.get("detected_only_after_correction"):
        return f"INCONCLUSIVE — {ef.reason}{tail}", "detected_only_after_correction"
    if ef.status != "DETECTED":
        return f"INCONCLUSIVE — {ef.reason}{tail}", "effect_inconclusive"
    if mv.status != "PASS":
        return ("INCONCLUSIVE — effect detected, but the metric's response to signal is untested "
                "(supply a positive control or an injected signal)" + tail, "detected_untested_metric")
    if "BIAS_UNSIZED" in mv.flags:
        return ("INCONCLUSIVE — effect detected, but a nuisance bias of the metric "
                f"({', '.join(mv.detail.get('unsized_bias', []))}) cannot be sized against the claim "
                "without a SESOI in the metric's units: declare a SESOI on the construct's scale with a probe "
                "(signal_test), which translates it, or on the metric's scale (sesoi_scale='observed'), and the "
                "bias is judged against bias_tolerance x SESOI" + tail, "bias_unsized")
    groups = [str(g) for g in (a.params or {}).get("groups") or ("groups[0]", "groups[1]")]
    direction = (a.prereg or {}).get("direction")
    if direction is None:
        return ("INCONCLUSIVE — effect detected, but no direction was pre-registered (declare "
                f"'increase' or 'decrease' from {groups[0]} to {groups[1]}, or 'two-sided')" + tail,
                "no_direction")
    effect = ef.detail.get("effect")
    if direction != "two-sided" and effect is not None and np.isfinite(effect):
        observed = "decrease" if effect > 0 else "increase"  # the effect is groups[0] - groups[1]
        if observed != direction:
            article = {"increase": "an", "decrease": "a"}
            return (f"NOT SUPPORTED — the effect is in the direction opposite to the pre-registered "
                    f"one: the metric shows {article[observed]} {observed} from {groups[0]} to "
                    f"{groups[1]}, the claim was {article[direction]} {direction}" + tail,
                    "opposite_direction")
    if direction == "two-sided":
        notes.append("non-directional claim")
        tail = f" [{'; '.join(notes)}]"
    ds = ef.detail.get("differential_sensitivity") or {}
    if ds.get("unsized"):
        return ("INCONCLUSIVE — effect detected, but the groups' sensitivity to the same coupling could not be "
                f"compared ({ds.get('error')}), so a difference it alone would make cannot be sized" + tail,
                "differential_sensitivity")
    if ds.get("applies"):
        adj = ds.get("adjusted") or {}
        what = (f"at an equal coupling the groups' unequal responses to it would make a difference of "
                f"{ds['artifact']:+.4g} [{ds['artifact_ci'][0]:+.4g}, {ds['artifact_ci'][1]:+.4g}] "
                f"({ds['rich']} responds {ds['k_rich']:+.4g}, {ds['poor']} {ds['k_poor']:+.4g}, ratio "
                f"{ds['ratio']:.3g}), beyond ±{ds['tolerance']:.4g}; the effect less its bound "
                f"{ds['offset']:+.4g} is {_adjusted_text(adj)}")
        if not adj.get("passes"):
            return f"NOT SUPPORTED — differential sensitivity: {what}{tail}", "differential_sensitivity"
        notes.append(f"judged on the effect less the differential sensitivity's bound {ds['offset']:+.4g}: {what}")
        tail = f" [{'; '.join(notes)}]"
    if rp.status == "NOT_REPLICATED":
        return f"NOT SUPPORTED — the effect did not replicate: {rp.reason}{tail}", "not_replicated"
    if ef.detail.get("tost_equivalent") and ef.detail.get("sesoi_source", "declared") == "declared":
        # detected, in the claimed direction, by the raw and the corrected test, and shown smaller
        # than the SESOI by the same report's TOST: not an effect that matters (v0.4 change 14)
        lo2, hi2 = ef.detail.get("ci_2alpha") or (float("nan"), float("nan"))
        alpha = float((a.prereg or {}).get("alpha", 0.05))
        return (f"EFFECT BELOW SESOI — {ef.reason}; its {1 - 2 * alpha:.0%} interval [{lo2:+.4g}, {hi2:+.4g}] "
                f"lies within ±{ef.detail.get('effective_sesoi'):.4g}, the effective SESOI: detected, but "
                f"smaller than the smallest effect that matters{tail}", "effect_below_sesoi")
    if a.judgment_pending:
        return (f"INCONCLUSIVE — effect detected; judgment gates 4 and 7 are unresolved{tail}",
                "judgment_pending")
    if rp.status == "REPLICATED":
        return f"SUPPORTED — replicated{tail}", "replicated"
    return f"SUPPORTED (provisional until replicated){tail}", "provisional"


def decide(a: "Autopsy") -> str:
    """The verdict of ``decide_cause`` (the single verdict rule)."""
    return decide_cause(a)[0]


def _adjusted_text(adj: dict) -> str:
    if not adj:
        return "not estimated"
    if adj.get("error"):
        return f"not estimated ({adj['error']})"
    e, ci = adj.get("effect"), adj.get("ci") or (float("nan"), float("nan"))
    if e is None or not np.isfinite(e):
        return f"not estimated ({adj.get('status')})"
    return (f"{e:+.4g} [{ci[0]:+.4g}, {ci[1]:+.4g}], p {adj['p']:.3g} ({adj.get('test')}): "
            + ("detected in the claim's direction" if adj.get("passes") else "not detected in the claim's direction"))


def detection_gap(tech: dict | None, cells: str) -> dict:
    """The premise of the differential sensitivity (the owner's decision of 2026-10-19): the groups differ in
    detection on the claim's cells — the genes detected per cell, the claim's own genes left out
    (``composition.technical_covariates`` with `exclude`), with a permutation p below alpha in a stratum or
    over the strata. A test that cannot be made leaves the premise open (the conservative side)."""
    if not tech or tech.get("error"):
        return dict(differs=True, cells=cells, rows=[], error=(tech or {}).get("error") or "not assessed")
    rows = [r for r in tech["rows"] if r["cells"] == cells and r["covariate"] == "genes per cell"]
    sig = [r for r in rows if r.get("difference") is not None and r["p"] < tech["alpha"]]
    return dict(differs=bool(sig), cells=cells, rows=rows, excluded=tech.get("excluded") or [], alpha=tech["alpha"])


def differential_sensitivity(d: dict | None, effect: Assessment | None, prereg: dict,
                             sesoi_metric: float | None, gap: dict | None = None) -> dict:
    """The verdict's use of GATE 4's differential sensitivity (``gates._differential``; the owner's decision
    of 2026-10-19), for a metric of structure. A, the difference the groups' unequal responses would make at
    an equal coupling (groups[0] - groups[1]), applies where the groups differ in detection (`gap`,
    ``detection_gap``: GATE 0's premise fails only there), A lies on the claim's side (the declared
    direction; for a two-sided claim the observed effect's) and its bound on that side exceeds
    bias_tolerance x the SESOI in the metric's units (0 without one): the effect is then judged less that
    bound (`offset`, in the effect's convention), and the verdict falls where the effect less it is not
    detected. Without a difference in detection A is reported only: the groups' responses then differ by
    their biology too (a coupling that changes with the genes' variance changes the metric's reliability:
    on ``tests/test_v03_verdict.py``'s clean data the check would remove a true effect,
    ``validation/exploratory/differential_sensitivity.log``, section 3), which is not an artifact. An A that
    cannot be sized where the groups differ in detection leaves the effect unsized."""
    if d is None:
        return {}
    rec = dict(d)
    gap = gap if gap is not None else dict(differs=True, rows=[], error="not assessed")
    rec.update(detection=gap, applies=False, offset=None, unsized=False)
    if d.get("error"):
        rec["unsized"] = bool(gap.get("differs"))
        return rec
    side = {"decrease": 1.0, "increase": -1.0}.get(prereg.get("direction"))
    if side is None:
        e = (effect.detail or {}).get("effect") if effect is not None else None
        side = float(np.sign(e)) if e is not None and np.isfinite(e) and e != 0 else None
    tol = float(prereg.get("bias_tolerance", DEFAULT_PREREG["bias_tolerance"])) * float(sesoi_metric) if sesoi_metric else 0.0
    A, (lo, hi) = d["artifact"], d["artifact_ci"]
    rec.update(tolerance=tol, side=side)
    if side is None:
        rec["in_direction"] = None
        return rec
    if not (np.isfinite(A) and np.isfinite(lo) and np.isfinite(hi)):
        rec.update(unsized=bool(gap.get("differs")),
                   error="the interval of the difference at an equal coupling could not be computed")
        return rec
    bound = max(side * lo, side * hi)
    rec.update(in_direction=bool(side * A > 0), bound=float(bound))
    if gap.get("differs") and side * A > 0 and bound > tol:
        rec.update(applies=True, offset=float(side * bound))
    return rec


def describe_differential(ds: dict) -> list:
    """Lines of the dossier for the differential sensitivity (``differential_sensitivity``)."""
    if not ds:
        return []
    if ds.get("error") and "responses" not in ds:
        return [f"not sized: {ds['error']}"]

    def iv(x, lo_hi, f="+.4g"):
        lo, hi = lo_hi
        return f"{x:{f}} [{lo:{f}}, {hi:{f}}]" if np.isfinite(lo) and np.isfinite(hi) else f"{x:{f}} [—]"
    change = ds.get("construct_change")
    out = [f"the same change of the construct ({f'{change:g}' if change is not None else 'the declared probe'}) "
           "injected into each group on the claim's units: " + "; ".join(f"{g} responds {iv(v['response'], v['ci'])} ({v['n_units']} units)"
                                            for g, v in ds["responses"].items())
           + f"; the ratio k_{ds['poor']} / k_{ds['rich']} = {iv(ds['ratio'], ds['ratio_ci'], '.3g')}",
           f"the metric minus its null in {ds['rich']}, the more sensitive group: {iv(ds['signal'], ds['signal_ci'])} "
           f"({ds['null_draws']} draws of GATE 0's null of the claim's genes per unit)",
           f"at an equal coupling the groups would differ by (metric − null in {ds['rich']}) × (1 − k_{ds['poor']} / "
           f"k_{ds['rich']}) = {iv(ds['artifact'], ds['artifact_ci'])}"]
    gap = ds.get("detection") or {}
    if gap:
        rows = "; ".join(f"{r['stratum']} {r['difference']:+.4g} (p {r['p']:.3g})" for r in gap.get("rows", [])
                         if r.get("difference") is not None)
        out.append(f"the groups' genes detected per cell on {gap.get('cells', 'the cells')}, the claim's genes left out: "
                   + (rows or f"not compared ({gap.get('error', 'no row')})") + " — "
                   + ("they differ in detection" if gap.get("differs") and not gap.get("error") else
                      "not tested, so taken to differ" if gap.get("differs") else
                      "no difference in detection shown, so the groups' responses may differ by their biology: reported only"))
    if ds.get("unsized"):
        out.append(f"not sized: {ds.get('error')}")
    elif not gap.get("differs", True):
        pass
    elif ds.get("side") is None:
        out.append("reported only: the claim has no direction and the effect none to compare with")
    elif ds.get("applies"):
        out.append(f"on the claim's side and beyond ±{ds['tolerance']:.4g} (bias_tolerance × the SESOI in the metric's "
                   f"units): the effect is judged less the bound {ds['offset']:+.4g}"
                   + (f": {_adjusted_text(ds['adjusted'])}" if "adjusted" in ds else
                      " (the effect is not detected, so the verdict does not change)"))
    elif ds.get("in_direction"):
        out.append(f"on the claim's side, its bound within ±{ds['tolerance']:.4g} (bias_tolerance × the SESOI in the "
                   "metric's units): reported only")
    else:
        out.append("on the other side of the claim's direction: reported only")
    return out


def capture_summary(effect: Assessment | None, design: Assessment | None, groups) -> dict:
    """The correction's line of the dossier (the owner's decision of 2026-10-19): from ``capture_check``, the
    groups' relative gap in genes detected per cell, groups[0] / groups[1] - 1, on the uncorrected data and
    at equal depth, on the domain's cells and on all cells; `widened` where equal depth made it larger."""
    ed = (effect.detail if effect is not None else None) or {}
    dd = (design.detail if design is not None else None) or {}
    sets = []
    dom = ed.get("domain") or {}
    if dom.get("capture_check"):
        sets.append(("the domain's cells", dom["capture_check"]))
        mix = ((ed.get("mixture") or {}).get("equalization") or {}).get("capture_check")
        if mix:
            sets.append(("all cells", mix))
    elif (dd.get("equalization") or {}).get("capture_check"):
        sets.append(("all cells", dd["equalization"]["capture_check"]))
    rows = []
    for where, cc in sets:
        r, e = cc["raw"]["relative_difference"], cc["equalized"]["relative_difference"]
        if not all(np.isfinite(x) for x in (r["genes"], e["genes"])):
            continue
        rows.append(dict(cells=where, genes_raw=float(r["genes"]), genes_equal_depth=float(e["genes"]),
                         depth_raw=float(r["depth"]), depth_equal_depth=float(e["depth"]),
                         widened=bool(abs(e["genes"]) > abs(r["genes"]))))
    if not rows:
        return {}
    return dict(groups=[str(g) for g in groups], rows=rows, widened=any(x["widened"] for x in rows))


def describe_capture(cap: dict) -> list:
    g0, g1 = cap["groups"]
    out = []
    for r in cap["rows"]:
        out.append(f"{r['cells']}: genes detected per cell, {g0} relative to {g1}, {r['genes_raw']:+.1%} uncorrected, "
                   f"{r['genes_equal_depth']:+.1%} at equal depth (the depth {r['depth_raw']:+.1%}, then "
                   f"{r['depth_equal_depth']:+.1%})"
                   + (" — WARNING: the correction widened the groups' gap in detection, so a metric that reads "
                      "detection is compared at a larger gap than the data's" if r["widened"] else ""))
    return out


def _shifted(metric, group_col: str, group, offset: float):
    """`metric` less `offset` on data whose cells all belong to `group` (groups[0]): every difference
    groups[0] - groups[1] the effect field forms (per cell set, replicate, unit, raw and at equal depth) then
    falls by `offset`, and its test becomes the test of the effect against `offset`."""
    def shifted(d):
        v = metric(d)
        g = np.asarray(d.obs[group_col])
        return v - offset if v is not None and len(g) and bool((g == group).all()) else v
    shifted.__name__ = getattr(metric, "__name__", "metric")
    return shifted


def _adjusted_record(adj: Assessment, side: float) -> dict:
    d = adj.detail or {}
    e = d.get("effect")
    p = d.get("p_perm") if d.get("tier") == "permutation" else d.get("p_t")
    passes = bool(adj.status == "DETECTED" and e is not None and np.isfinite(e) and side * e > 0)
    return dict(effect=e, ci=d.get("ci"), p=p, test=d.get("test"), status=adj.status, passes=passes,
                reason=adj.reason)


@dataclass
class Autopsy:
    metric_name: str
    results: list[GateResult] = field(default_factory=list)
    prereg: dict = field(default_factory=dict)
    metric_validity: Assessment | None = None
    design_adequacy: Assessment | None = None
    effect: Assessment | None = None
    replication: Assessment | None = None
    params: dict = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)
    # v0.4: the composition of the comparison and the claim per type (change 7), a second measure of
    # a pair's coupling beside the metric (change 15) and the check of specificity (change 13);
    # reported, not gates
    composition: dict = field(default_factory=dict)
    reference: dict = field(default_factory=dict)
    specificity: dict = field(default_factory=dict)
    # the groups' technical covariates per replicate, counts and genes per cell (a line of the dossier)
    technical: dict = field(default_factory=dict)
    # how the correction changed the groups' gap in genes detected per cell (capture_check; a line of the
    # dossier, with a warning where it grew: the owner's decision of 2026-10-19)
    capture: dict = field(default_factory=dict)

    @property
    def judgment_pending(self) -> bool:
        return bool(self.prereg.get("judgment_pending", True))

    @property
    def verdict(self) -> str:
        return decide(self)

    @property
    def cause(self) -> str:
        """The cause code of the verdict (``CAUSES``)."""
        return decide_cause(self)[1]

    # ------------------------------------------------------------------ output
    def fields(self) -> dict:
        return dict(metric_validity=self.metric_validity, design_adequacy=self.design_adequacy,
                    effect=self.effect, replication=self.replication)

    def to_markdown(self) -> str:
        lines = [f"# Metric autopsy — {_cell(self.metric_name)}", ""]
        prereg_items = {k: v for k, v in self.prereg.items() if k not in DEFAULT_PREREG or k == "estimand"}
        if prereg_items:
            lines += ["## Pre-registration", ""]
            for k, v in prereg_items.items():
                lines.append(f"- **{_cell(k)}:** {_cell(v)}")
            lines.append("")
        lines += ["## Verdict", "", f"**{self.verdict}**", "", f"Cause: `{self.cause}`", ""]
        lines += ["## Assessment", "", "| Field | Status | Reason |", "|---|---|---|"]
        for name, asmt in self.fields().items():
            if asmt is None:
                continue
            flags = f" ({', '.join(asmt.flags)})" if asmt.flags else ""
            lines.append(f"| {name} | **{asmt.status}**{_cell(flags)} | {_cell(asmt.reason)} |")
        lines += ["", "## Gates", "", "| Gate | Name | Status | Finding |", "|---|---|---|---|"]
        for r in sorted(self.results, key=lambda r: r.gate):
            lines.append(f"| {r.gate} | {_cell(r.name)} | **{r.status.value}** | {_cell(r.message)} |")
        g4 = next((r for r in self.results if r.gate == 4), None)
        if g4 is not None and g4.detail.get("dose_response"):
            lines += ["", "## The probe's dose and the metric's response (the SESOI in the metric's units)", "",
                      "| Construct change | Dose | Response | 95% CI |", "|---|---|---|---|"]
            for pnt in g4.detail["dose_response"]:
                dose = pnt.get("dose") or {}
                dtxt = (f"strength {dose['strength']:.4g}" if "strength" in dose else
                        f"fold 1/{dose['fold']:.4g}" if "fold" in dose else "—")
                lo, hi = pnt["ci"]
                lines.append(f"| {pnt['construct_change']:g} | {_cell(dtxt)} | {pnt['response']:+.4g} | "
                             f"[{lo:+.4g}, {hi:+.4g}] |")
        g0n, g1n = [str(g) for g in (self.params.get("groups") or ["groups[0]", "groups[1]"])[:2]]
        lines += ["", "## The dossier (reported beside the verdict, not gates)", "",
                  f"Every difference below is {_cell(g0n)} − {_cell(g1n)} (groups[0] − groups[1], the effect test's "
                  f"convention): positive where {_cell(g0n)} is higher, so a decrease from {_cell(g0n)} to "
                  f"{_cell(g1n)} is a positive difference."]
        ds = ((self.effect.detail if self.effect is not None else None) or {}).get("differential_sensitivity")
        if ds:
            lines += ["", "## Differential sensitivity (a metric of structure: each group's response to the same "
                          "coupling)", ""]
            for ln in describe_differential(ds):
                lines.append(f"- {_cell(ln)}")
        if self.capture.get("rows"):
            lines += ["", "## The correction and the groups' detection (reported, not a gate)", ""]
            for ln in describe_capture(self.capture):
                lines.append(f"- {_cell(ln)}")
        dom = (self.effect.detail.get("domain") if self.effect is not None else None) or {}
        if dom:
            from .domain import describe as _ddesc
            lines += ["", "## The domain (the verdict is judged within it)", ""]
            for ln in _ddesc(self.composition.get("domain") or {}) if self.composition.get("domain") else []:
                lines.append(f"- {_cell(ln)}")
            mix = self.effect.detail.get("mixture") or {}
            out = dom.get("outside") or {}

            def _fx(x):
                if not x or x.get("effect") is None or not np.isfinite(x["effect"]):
                    return "—"
                return f"{x['effect']:+.4g} (p {x['p']:.3g}, {x['tier']})"
            ef = self.effect.detail
            g0, g1 = (self.params.get("groups") or ["groups[0]", "groups[1]"])[:2]
            p_dom = ef.get("p_perm") if ef.get("tier") == "permutation" else ef.get("p_t")
            dom_corr = _fx(dict(effect=ef.get("effect"), p=p_dom, tier=ef.get("tier")))
            dom_raw = _fx(dict(effect=ef.get("raw_effect"), p=ef.get("p_raw"), tier=ef.get("tier")))
            lines += ["", f"| Cells | {_cell(g0)} − {_cell(g1)}, corrected | raw |", "|---|---|---|",
                      f"| the domain ({_cell(', '.join(dom.get('types', [])) or 'none')}) | {_cell(dom_corr)} | "
                      f"{_cell(dom_raw)} |"]
            if out:
                lines.append(f"| outside the domain ({len(out.get('types', []))} types) | {_cell(_fx(out.get('corrected')))} | "
                             f"{_cell(_fx(out.get('raw')))} |")
            if mix:
                lines.append(f"| the whole mixture | {_cell(_fx(mix.get('corrected')))} | {_cell(_fx(mix.get('raw')))} |")
            per = dom.get("per_type") or []
            if per:
                lines += ["", "| Type | In the domain | Detection | The genes' counts per cell (log1p CP10k) | "
                              "Effect (corrected) | Effect (raw) |", "|---|---|---|---|---|---|"]
                for r in per:
                    det = r.get("detection") or {}
                    dtxt = "; ".join(f"{g}: " + ", ".join(f"{k} {v:.0%}" for k, v in d.items()) for g, d in det.items())
                    itxt = "; ".join(f"{g}: {v:.3g}" for g, v in (r.get("index") or {}).items()
                                     if v is not None and np.isfinite(v)) or "—"
                    lines.append(f"| {_cell(r['type'])} | {'yes' if r['in_domain'] else 'no'} | {_cell(dtxt)} | "
                                 f"{_cell(itxt)} | {_cell(_fx(r.get('corrected')))} | {_cell(_fx(r.get('raw')))} |")
        ef_detail = (self.effect.detail if self.effect is not None else None) or {}
        unmeasured = list(ef_detail.get("unmeasured") or [])
        unmeasured_units = [u for u in ((ef_detail.get("domain") or {}).get("unmeasured") or []) if u.get("in_domain")]
        if unmeasured or unmeasured_units:
            lines += ["", "## Left out: the metric gives no value (not measured)", ""]
            for u in unmeasured:
                lines.append(f"- replicate {_cell(u['replicate'])}: uncorrected "
                             f"{'measured' if u.get('raw') is not None else 'no value'}, at equal depth "
                             f"{'measured' if u.get('corrected') is not None else 'no value'}")
            for u in unmeasured_units:
                lines.append(f"- unit ({_cell(u['replicate'])}, {_cell(u['type'])}) of the domain: uncorrected "
                             f"{'measured' if u.get('raw') is not None else 'no value'}, at equal depth "
                             f"{'measured' if u.get('corrected') is not None else 'no value'}")
        if self.composition:
            from .composition import describe as _describe
            lines += ["", "## Composition (reported, not a gate)", ""]
            for comp in self.composition.get("columns", {}).values():
                lines.append(f"- {_cell(_describe(comp))}")
        if self.technical:
            from .composition import describe_technical as _tdesc
            lines += ["", "## Technical covariates per replicate (reported, not a gate)", ""]
            for ln in _tdesc(self.technical):
                lines.append(f"- {_cell(ln)}")
        if self.reference.get("rows"):
            ref = self.reference
            lines += ["", f"## {_cell(ref['measure'])} of {_cell(' and '.join(ref['genes']))} beside the metric "
                          "(reported, not a gate)", "",
                      f"| Stratum | {_cell(self.metric_name)}, {_cell(ref['groups'][0])} − {_cell(ref['groups'][1])} | "
                      f"{_cell(ref['measure'])}, {_cell(ref['groups'][0])} − {_cell(ref['groups'][1])} |", "|---|---|---|"]
            for row in ref["rows"]:
                cells = [f"{x['difference']:+.4g} (p {x['p']:.3g})" if x.get("difference") is not None else "—"
                         for x in (row["metric"], row["reference"])]
                lines.append(f"| {_cell(row['stratum'])} | {_cell(cells[0])} | {_cell(cells[1])} |")
        if self.specificity:
            from .specificity import describe as _sdesc
            lines += ["", "## Specificity (reported, not a gate)", ""]
            for ln in _sdesc(self.specificity):
                lines.append(f"- {_cell(ln)}")
        prov = self.provenance
        if prov:
            lines += ["", "## Provenance", ""]
            for key in ("data_sha256", "prereg_sha256", "claim_id"):
                if prov.get(key):
                    lines.append(f"- {key}: `{prov[key][:16]}…`")
            env = prov.get("environment", {})
            if env:
                lines.append(f"- metric-autopsy {env.get('metric_autopsy')}, seed {self.params.get('seed')}")
            log = prov.get("log") or {}
            if log.get("path"):
                lines.append(f"- run log `{_cell(log['path'])}`: attempt {log.get('attempt')} for this claim")
        lines.append("")
        return "\n".join(lines)

    def to_dict(self) -> dict:
        def asm(a):
            return None if a is None else dict(status=a.status, reason=a.reason, flags=list(a.flags),
                                               detail=_slim(a.detail))
        return _prov.jsonable(dict(
            schema=_prov.SCHEMA, metric=self.metric_name, verdict=self.verdict, cause=self.cause,
            fields={k: asm(v) for k, v in self.fields().items()},
            gates=[dict(gate=r.gate, name=r.name, status=r.status.value, message=r.message,
                        detail=_slim(r.detail)) for r in sorted(self.results, key=lambda r: r.gate)],
            prereg=self.prereg, params=self.params, provenance=self.provenance,
            composition={k: v for k, v in self.composition.items() if k != "domain"} | (
                {"domain": _public_domain(self.composition["domain"])} if self.composition.get("domain") else {}),
            reference=self.reference, specificity=self.specificity, technical=self.technical, capture=self.capture,
            differences=(lambda g: f"{g[0]} - {g[1]}")([str(x) for x in (self.params.get("groups") or
                                                                         ["groups[0]", "groups[1]"])[:2]]),
        ))

    def to_json(self, **kw) -> str:
        kw.setdefault("indent", 2)
        return json.dumps(self.to_dict(), allow_nan=False, **kw)

    def save_json(self, path) -> None:
        with open(path, "w") as fh:
            fh.write(self.to_json())


def _public_domain(dom: dict) -> dict:
    """The domain without its internal design record."""
    return {k: v for k, v in dom.items() if not k.startswith("_")}


def _slim(detail: dict) -> dict:
    """Drop bulky per-cell payloads (e.g. GATE 3 scatter) from serialized details."""
    out = {}
    for k, v in (detail or {}).items():
        if k == "scatter":
            out[k] = "omitted"
        elif k == "mask":
            continue
        else:
            out[k] = v
    return out


# --------------------------------------------------------------------------- #
# field assembly
# --------------------------------------------------------------------------- #
def _controls_tie(metric, pair_metric, data, gene_pair) -> tuple[bool, str, dict]:
    """GATE 5 runs `pair_metric` on the control pairs, so it is evidence about `metric` only if
    `pair_metric` bound to `gene_pair` *is* `metric`. Both are black boxes: compare their values on
    the data. A mismatch would let a valid pair metric certify a blind or random one."""
    if gene_pair is None:
        return False, ("no gene_pair, so pair_metric cannot be checked against the metric: the "
                       "controls would test another metric"), {}
    v = _g._safe_call(lambda: metric(data))
    w = _g._safe_call(lambda: pair_metric(data, gene_a=gene_pair[0], gene_b=gene_pair[1]))
    detail = dict(metric_value=v, pair_metric_value=w)
    if v is not None and w is not None and np.isclose(v, w, rtol=1e-6, atol=1e-12):
        return True, "", detail
    fmt = (lambda x: "no finite value" if x is None else f"{x:.6g}")
    return False, (f"pair_metric on {gene_pair[0]}–{gene_pair[1]} gives {fmt(w)} and the metric gives "
                   f"{fmt(v)}: the controls would test another metric (bind the metric from "
                   "pair_metric, and seed it if it is stochastic)"), detail


_BIAS_FLAGS = {"corrected": "DEPTH_BIAS_CORRECTED", "below tolerance": "BIAS_BELOW_TOLERANCE",
               "unsized": "BIAS_UNSIZED"}


def _metric_validity(g0: GateResult, g5: GateResult | None, g4: GateResult | None) -> Assessment:
    lvl = g0.detail.get("level_shifts", {}) if g0 else {}
    handling = g0.detail.get("bias_handling", {}) if g0 else {}
    detail = dict(gate0=g0.status.value if g0 else None, level_shifts=lvl,
                  gate5=g5.status.value if g5 else None, gate4=g4.status.value if g4 else None,
                  bias_handling=handling,
                  unsized_bias=sorted(k for k, h in handling.items() if h == "unsized"),
                  failed_gate=None)
    flags = (["LEVEL_SHIFT"] if lvl else []) + sorted(
        {_BIAS_FLAGS[h] for h in handling.values() if h in _BIAS_FLAGS})
    bias_txt = (g0.message if g0 is not None and g0.status == GateStatus.WARN
                else "no nuisance bias")
    if g0 is None or g0.status == GateStatus.SKIP:
        base = Assessment("UNTESTED", "nuisance invariance could not be evaluated (GATE 0 skipped)", [], detail)
    elif g0.status == GateStatus.DEGENERATE:
        return Assessment("DEGENERATE", g0.message, [], detail)
    elif g0.status == GateStatus.FAIL:
        return Assessment("FAIL", g0.message, [], {**detail, "failed_gate": 0})
    else:
        base = None
    if g5 is not None and g5.status == GateStatus.FAIL:
        return Assessment("FAIL", f"controls: {g5.message}", [], {**detail, "failed_gate": 5})
    if g4 is not None and g4.status == GateStatus.FAIL:
        return Assessment("FAIL", f"injected signal: {g4.message}", [], {**detail, "failed_gate": 4})
    if base is not None:
        return base
    # GATE 4 judges the response on the analysed construct itself. Where it ran and could show
    # neither a response nor blindness, a positive control on another pair does not stand in
    # for it (decided 2026-10-08, JOURNAL.md D5).
    g4_untested = g4 is not None and g4.status == GateStatus.WARN
    evidence = []
    if g4 is not None and g4.status == GateStatus.PASS:
        evidence.append(g4.message)
    elif not g4_untested:
        if g5 is not None and g5.status == GateStatus.PASS:
            evidence.append("positive control beats the empirical null in every stratum")
        elif g5 is not None and g5.status == GateStatus.WARN and g5.detail.get("pos_demonstrated"):
            evidence.append(f"positive control beats the empirical null in {g5.detail['n_pos_fires']}/"
                            f"{len(g5.detail['rows'])} strata")
    lvl_txt = (" Level shift (reported, not failed): "
               + ", ".join(f"{k} {v:+.0%}" for k, v in lvl.items()) + " of the effect scale.") if lvl else ""
    if evidence:
        return Assessment("PASS", bias_txt + "; " + "; ".join(evidence) + "." + lvl_txt,
                          flags, detail)
    pc_fires = g5 is not None and g5.detail.get("pos_demonstrated")
    why = ((f"GATE 4: {g4.message}" + ("; the positive control responds, but on another pair"
                                       if pc_fires else ""))
           if g4_untested else
           "the positive control did not beat its null in any stratum"
           if g5 is not None and g5.status == GateStatus.WARN else
           f"controls not run: {g5.message}" if g5 is not None and g5.status == GateStatus.SKIP else
           "supply a positive control pair or an injected signal")
    return Assessment("UNTESTED", f"{bias_txt}, but no demonstrated response to signal ({why})."
                      + lvl_txt, flags, detail)


def _design_adequacy(g1: GateResult, design: dict, prereg: dict, attenuation: dict | None = None,
                     sesoi_info: dict | None = None) -> Assessment:
    flags, notes = [], list(design.get("notes", []))
    attenuation = dict(attenuation or {})
    sesoi_info = dict(sesoi_info or {})
    tier = design.get("tier")
    if tier == "parametric":
        flags.append("PARAMETRIC_ONLY")
    if design.get("underpowered"):
        flags.append("UNDERPOWERED")
    elif design.get("underpowered") is None and prereg.get("sesoi") is None and tier not in (None, "insufficient"):
        flags.append("POWER_NOT_ASSESSED")
    if design.get("sesoi_below_floor"):
        flags.append("SESOI_BELOW_FLOOR")
    if sesoi_info.get("not_translated"):
        flags.append("SESOI_NOT_TRANSLATED")
        notes.append(sesoi_info["note"])
    if attenuation:
        flags.append("ATTENUATION")
    if design.get("composite_units"):
        flags.append("COMPOSITE_UNITS")
    detail = dict(qc=g1.status.value, correction=design.get("correction"), tier=tier,
                  replication=design.get("replication"), power=design.get("power"),
                  attenuation=attenuation, equalization=design.get("equalization"), notes=notes,
                  composite_units=design.get("composite_units"), sesoi=sesoi_info or None)
    if g1.status == GateStatus.STOP:
        return Assessment("UNIDENTIFIABLE", g1.message, flags, detail)
    if not design.get("identifiable", True):
        return Assessment("UNIDENTIFIABLE", design.get("reason", "comparison not identifiable"), flags, detail)
    if tier == "insufficient":
        return Assessment("INSUFFICIENT_REPLICATION", design.get("tier_reason", "insufficient replication"),
                          flags, detail)
    power = design.get("power") or {}
    att = ""
    if attenuation:
        att = "; attenuation (GATE 0) " + ", ".join(f"{k} −{v:.0%}" for k, v in attenuation.items())
    pw = ""
    floor = power.get("floor")
    eff = power.get("effective_sesoi")
    if power.get("sesoi_source") == "floor" and floor is not None and np.isfinite(floor):
        pw = (f"; no SESOI in the metric's units: the design's floor {floor:.4g} is the margin of absence, so "
              f"absence is confirmed only for effects of at least {floor:.4g}")
    elif eff is not None and np.isfinite(power.get("mde", np.nan)):
        label = power.get("sesoi_label") or (
            f"SESOI {power['sesoi']:g} on the metric's scale" if prereg.get("sesoi_scale") == "observed" else
            f"SESOI {power['sesoi']:g} x λ {power['lambda_']:.2f}")
        pw = f"; MDE {power['mde']:.4g} vs the SESOI in the metric's units {eff:.4g} ({label})"
        if floor is not None and np.isfinite(floor):
            pw += f"; the design's floor {floor:.4g}"
            if design.get("sesoi_below_floor"):
                pw += (": the SESOI lies below it, so this design cannot show the absence of a change that "
                       "matters")
    rep = design.get("replication") or {}
    rep_txt = f"{rep.get('kind')} replicates {rep.get('counts')}" if rep else ""
    if g1.status == GateStatus.WARN and design.get("correction") in ("depth_thinning", "capture_thinning"):
        return Assessment("CORRECTED", f"{g1.message.split(' — ')[0]}; removed by {design['correction']}; "
                          f"{rep_txt}{att}{pw}", flags, detail)
    return Assessment("ADEQUATE", f"QC parity {g1.status.value}; correction {design.get('correction')}; "
                      f"{rep_txt}{att}{pw}", flags, detail)


def _gate2_result(effect: Assessment, design: dict) -> GateResult:
    d = effect.detail
    if not design.get("identifiable", True):
        return GateResult(2, "Estimand-dependent correction", GateStatus.STOP, design.get("reason", ""), d)
    raw, est = d.get("raw_effect"), d.get("effect")
    txt = (f"{design.get('correction')}: raw {raw:+.4g} -> corrected {est:+.4g}"
           if isinstance(raw, float) and isinstance(est, float) and np.isfinite(raw) and np.isfinite(est)
           else f"{design.get('correction')}: {effect.reason}")
    if d.get("explained_by_depth"):
        return GateResult(2, "Estimand-dependent correction", GateStatus.WARN,
                          txt + " — the raw difference is explained by the technical difference", d)
    if d.get("detected_only_after_correction"):
        return GateResult(2, "Estimand-dependent correction", GateStatus.WARN,
                          txt + " — only the corrected effect is detected: a correction can only remove a "
                          "conclusion, never create one", d)
    return GateResult(2, "Estimand-dependent correction", GateStatus.PASS, txt, d)


def _replication(g6: GateResult | None) -> Assessment:
    if g6 is None:
        return Assessment(NOT_RUN, "no independent dataset supplied")
    status = g6.detail.get("replication", NOT_RUN)
    return Assessment(status, g6.message, list(g6.detail.get("effect_flags", [])), g6.detail)


TECHNICAL_DRAWS = 3  # GATE 0's null draws averaged per unit in the groups' technical difference (change 12)


def technical_difference(metric, data, genes, *, group_col: str, groups, within=(), replicate_col: str,
                         estimand, margin: float, alpha: float = 0.05, power: float = 0.8, n_perm: int = 1000,
                         min_replicates: int = 3, spikein_prefix: str = "ERCC-", domain: dict | None = None,
                         seed: int = 0, n_draws: int = TECHNICAL_DRAWS) -> dict:
    """v0.4 change 12: the groups' technical difference — the claim's metric on GATE 0's null (each claim
    gene's counts redrawn independently from cells at least as deep within the unit, ``gates.claim_null``:
    no gene-gene biology, each group's own technology), averaged over `n_draws` draws per unit, estimated
    between the groups as the effect is (``effect.estimate_effect``: the same units, domain, correction
    and test). Equivalent to no difference where its 90% interval lies within ±`margin` (TOST); a null
    the metric cannot measure (fewer replicates than the test needs) is not equivalent."""
    from .effect import estimate_effect

    rng = np.random.default_rng(seed)

    def null_metric(d):
        sd = _g._as_simple(d)
        vals = [v for v in (_g._safe_call(lambda: metric(_g.claim_null(sd, genes, rng))) for _ in range(n_draws))
                if v is not None]
        return float(np.mean(vals)) if vals else float("nan")

    try:
        eff, _design = estimate_effect(
            null_metric, data, group_col=group_col, groups=groups, within=within, replicate_col=replicate_col,
            estimand=estimand, sesoi=margin, alpha=alpha, power=power, min_replicates=min_replicates,
            n_perm=n_perm, spikein_prefix=spikein_prefix, sesoi_scale="observed", seed=seed, domain=domain)
    except Exception as e:  # never silent: the test is recorded as not done, and the bias keeps blocking
        return dict(equivalent=False, margin=margin, alpha=alpha, error=f"{type(e).__name__}: {e}")
    d = eff.detail
    return dict(equivalent=bool(d.get("tost_equivalent")), margin=float(margin), alpha=alpha,
                effect=d.get("effect"), ci=d.get("ci"), ci_2alpha=d.get("ci_2alpha"), status=eff.status,
                reason=eff.reason, test=d.get("test"), tier=d.get("tier"), estimate_over=d.get("estimate_over"),
                n_replicates=d.get("n_replicates"), draws=n_draws,
                error=None if d.get("effect") is not None else eff.reason)


def coupling_reference(metric, data, gene_pair, *, group_col: str, groups, replicate_col: str, within=(),
                       alpha: float = 0.05, n_perm: int = 1000, seed: int = 0, domain: dict | None = None) -> dict:
    """The engine's own ``norm_pearson`` of the pair beside the user's metric, on the uncorrected data
    at the replicate level, by stratum and over the strata (v0.4 change 15, decided 2026-10-17: on
    preprint section 4 ``mi_3bin`` fell where ``norm_pearson`` rose). Not a gate: always shown. With a
    `domain` both are combined within it, the types weighted equally, as the verdict is. The differences
    are groups[0] - groups[1], the effect test's convention, as every line of the dossier (the owner's
    decision of 2026-10-19; v0.4.0.dev0 gave groups[1] - groups[0] here), each with the permutation (or,
    below its reach, the t) p."""
    from functools import partial
    from .effect import _infer, _replicate_design, _rep_values, _restrict, _units_inference
    from .metrics import norm_pearson

    ref = partial(norm_pearson, gene_a=gene_pair[0], gene_b=gene_pair[1])
    rng = np.random.default_rng(seed)
    if domain is not None and not domain.get("error"):
        from .domain import build_units, unit_values
        units = build_units(domain, data, group_col=group_col, groups=groups, replicate_col=replicate_col)
        if units is None:
            return dict(error="the domain holds no type with enough cells and replicates")
        vals = {k: unit_values(m, data, units) for k, m in (("metric", metric), ("reference", ref))}

        def summary_u(v, sub):
            inf = _units_inference(sub, v, groups, "permutation", alpha, n_perm, rng)
            p = inf["perm"].p if inf["tier"] == "permutation" else inf["p_t"]
            return dict(difference=float(inf["est"]) if np.isfinite(inf["est"]) else None, p=float(p),
                        tier=inf["tier"], n=inf["n_counts"])
        rows = []
        rs = units.rep_stratum[units.unit_rep]
        for si, st in enumerate(units.strata) if len(units.strata) > 1 else []:
            sel = rs == si
            if sel.any():
                sub = units.subset(sel)
                rows.append(dict(stratum=" ".join(st.split("\x1f")) or "all",
                                 **{k: summary_u(v[sel], sub) for k, v in vals.items()}))
        rows.append(dict(stratum="over the strata" if len(units.strata) > 1 else "all",
                         **{k: summary_u(v, units) for k, v in vals.items()}))
        return dict(measure="norm_pearson (the engine's)", genes=list(gene_pair), groups=[str(g) for g in groups],
                    data="uncorrected", over="the domain's types, equally weighted", rows=rows)
    sd = _restrict(data, group_col, groups)
    design = _replicate_design(sd.obs, group_col, groups, replicate_col, list(within))
    if design["kind"] == "mixed":
        return dict(error=design["reason"])
    values = {k: _rep_values(m, sd, group_col, groups, replicate_col, design) for k, m in (("metric", metric), ("reference", ref))}

    def summary(vals, sub):
        inf = _infer(vals, sub, groups, "permutation", alpha, n_perm, rng)
        p = inf["perm"].p if inf["tier"] == "permutation" else inf["p_t"]
        return dict(difference=float(inf["est"]) if np.isfinite(inf["est"]) else None, p=float(p),
                    tier=inf["tier"], n=inf["n_counts"])
    rows = []
    strata = list(dict.fromkeys(design["strata"][r] for r in design["reps"]))
    for st in strata if len(strata) > 1 else []:
        reps = [r for r in design["reps"] if design["strata"][r] == st]
        sub = dict(design, reps=reps, strata={r: "" for r in reps})
        rows.append(dict(stratum=" ".join(st.split("\x1f")) or "all",
                         **{k: summary({r: v[r] for r in reps}, sub) for k, v in values.items()}))
    rows.append(dict(stratum="over the strata" if len(strata) > 1 else "all",
                     **{k: summary(v, design) for k, v in values.items()}))
    return dict(measure="norm_pearson (the engine's)", genes=list(gene_pair), groups=[str(g) for g in groups],
                data="uncorrected", over="the replicates", rows=rows)


def _sesoi_in_metric_units(prereg: dict, g4: GateResult | None, signal_test) -> dict:
    """The SESOI in the metric's units (v0.4 change 3): GATE 4's response s_m to a change of the
    construct of the SESOI's size; the declared SESOI when it is on the metric's scale; none (the
    effect field then uses the design's floor) without a SESOI or without a probe that translates it."""
    sesoi, scale = prereg.get("sesoi"), prereg.get("sesoi_scale", "construct")
    if sesoi is None:
        return dict(value=None, source="none", label=None,
                    note="no SESOI pre-registered: the design's floor is the margin of absence")
    if scale == "observed":
        return dict(value=float(sesoi), source="declared, the metric's scale",
                    label=f"SESOI {float(sesoi):g} on the metric's scale")
    s_m = g4.detail.get("sesoi_metric") if g4 is not None else None
    if s_m is not None:
        lo, hi = g4.detail.get("sesoi_metric_ci") or (float("nan"), float("nan"))
        return dict(value=float(s_m), source="s_m", ci=(lo, hi), construct=float(sesoi),
                    label=(f"s_m, GATE 4's response to a change of {float(sesoi):g} on the construct's scale, "
                           f"95% CI [{lo:+.4g}, {hi:+.4g}] over the replicates"))
    why = ("no probe (signal_test)" if signal_test is None else
           "the probe's response was not shown in the injection's direction" if g4 is not None and g4.detail.get("dose")
           else "no rule translates this probe's dose" if getattr(signal_test, "kind", None) not in ("coupling", "module")
           else "GATE 4 could not measure the claim's statistic")
    return dict(value=None, source="not translated", label=None, not_translated=True,
                note=(f"the SESOI {float(sesoi):g} is on the construct's scale and was not translated into the "
                      f"metric's units ({why}): the design's floor is the margin of absence; supply a coupling or "
                      "module probe, or declare the SESOI on the metric's scale (sesoi_scale='observed')"))


# --------------------------------------------------------------------------- #
# orchestration
# --------------------------------------------------------------------------- #
def run_autopsy(
    metric,
    data,
    *,
    group_col: str,
    groups: tuple,
    gene_pair: tuple | None = None,
    within: Sequence[str] = (),
    pair_metric=None,
    pos_pair: tuple | None = None,
    neg_pair: tuple | None = None,
    data2=None,
    prereg: dict | None = None,
    include_matrix_perturbations: bool | None = None,
    stop_on_first_fail: bool = True,
    replicate_col: str | None = None,
    signal_test=None,
    seed: int = 0,
    n_perm: int = 1000,
    log_path=None,
    composition: Sequence[str] = (),
    gene_metric=None,
    specificity: bool = True,
    specificity_pool=None,
) -> Autopsy:
    """Run the gates and assemble the four-field verdict.

    `metric(data) -> float` is the bound metric. Controls (GATE 5) need the unbound
    `pair_metric(data, *, gene_a, gene_b)`, `gene_pair` and both control pairs; they count only if
    `pair_metric` on `gene_pair` gives the metric's value on the data (else GATE 5 is SKIP: the
    controls would test another metric). `replicate_col` names the biological replicate (mouse, donor,
    plate): without it there is no effect verdict. `signal_test` is an ``injected_signal`` probe.
    `prereg` carries ``estimand`` ('composition' | 'content'), ``sesoi`` (on the construct's scale, a
    change of correlation for a coupling probe or of natural-log expression for a module probe, or on
    the metric's scale with ``sesoi_scale="observed"``), ``min_replicates``, ``alpha``,
    ``judgment_pending``, optionally ``known_interactions`` (pairs) and ``known_to_differ`` (genes) for
    the check of specificity, and free-text commitments. Its hash is recorded, and a run with a
    pre-registration is appended to the run log (`log_path`, else ``$METRIC_AUTOPSY_LOG``, else
    ``metric_autopsy_runs.jsonl``; ``"off"`` disables). `stop_on_first_fail` is kept for compatibility and
    recorded: since v0.4 the effect and every report line are estimated after a FAIL too (the dossier,
    flag ``DOSSIER_ONLY``), and the verdict follows the first FAIL either way.

    v0.4: GATE 4 runs first, and with a construct-scale SESOI its probe is a change of the construct of
    that size, whose response s_m is the SESOI in the metric's units for GATE 0's tolerance, the TOST,
    EFFECT BELOW SESOI and the comparison with the design's floor (change 3 and 4). `composition` names
    obs columns of an annotation, the cell type first: their shares are compared between the groups,
    and the claim is judged within its domain, the types where its genes are detected
    (``domain.claim_domain``), the types weighted equally; GATE 4 runs on the domain's cells, and the
    whole mixture's estimate is reported beside the difference of composition. For a claim about genes,
    the report carries the check of specificity among matched controls (change 13): `pair_metric` (a
    pair) or `gene_metric(data, *, genes)` (a set) gives the metric of other genes, drawn from every gene
    of `specificity_pool` (the same cells with every gene; default `data`). For a pair the engine's own
    ``norm_pearson`` is reported beside the metric (change 15). `delta_min` in a pre-registration is
    recorded and not used.

    For a coupling probe (a metric of structure, ``STRUCTURE_PROBES``) GATE 4 also injects the same change
    of the construct into the other group, and where the groups differ in detection (genes per cell, the
    claim's genes left out, ``detection_gap``) the effect is judged less the bound of the difference their
    unequal responses would make at an equal coupling (``differential_sensitivity``; the owner's decision
    of 2026-10-19), with the same estimate and test (``_shifted``). The dossier carries it always, and the
    correction's change of the groups' gap in genes per cell (``capture_summary``).
    """
    name = _metric_name(metric)
    preregistered = bool(prereg)
    signal_direction = (prereg or {}).get("signal_direction")  # explicit, before the defaults
    prereg = normalize_prereg(prereg)
    groups = tuple(groups)
    within = list(within)
    alpha = float(prereg["alpha"])
    sesoi = prereg.get("sesoi")
    sesoi_scale = prereg.get("sesoi_scale", "construct")
    results: list[GateResult] = []
    composition = [composition] if isinstance(composition, str) else list(composition or ())
    params = dict(metric=name, group_col=group_col, groups=list(groups), within=within,
                  gene_pair=list(gene_pair) if gene_pair else None,
                  pair_metric=_metric_name(pair_metric) if pair_metric else None,
                  gene_metric=_metric_name(gene_metric) if gene_metric else None,
                  pos_pair=list(pos_pair) if pos_pair else None,
                  neg_pair=list(neg_pair) if neg_pair else None, replicate_col=replicate_col,
                  signal_test=getattr(signal_test, "description", None) if signal_test else None,
                  data2=data2 is not None, include_matrix_perturbations=include_matrix_perturbations,
                  stop_on_first_fail=stop_on_first_fail, seed=seed, n_perm=n_perm, composition=composition,
                  specificity=bool(specificity),
                  specificity_pool=list(getattr(specificity_pool.X, "shape", ())) if specificity_pool is not None else None)
    if prereg.get("delta_min") is not None:
        params["delta_min_ignored"] = (f"delta_min {prereg['delta_min']} is recorded and not used: since v0.4 GATE 4 "
                                       "decides only whether the metric is blind (change 3)")

    if include_matrix_perturbations is None:
        include_matrix_perturbations = gene_pair is None

    # every random draw is keyed by the data and the claim's canonical form (v0.4 change 1): a
    # mirrored wording gets the same draws, two claims on the same data different ones
    data_sha = _prov.sha256_data(data, _hashed_columns(group_col, within, replicate_col))
    claim_genes = list(gene_pair or getattr(signal_test, "genes", None) or ())
    claim = _prov.canonical_claim(name, claim_genes, groups, prereg.get("direction"))
    params["claim"] = claim

    def seed_of(stream: str) -> int:
        return _prov.claim_seed(data_sha, claim, seed, stream)

    from .effect import _cell_effect, equalizing_correction, estimate_effect

    # ---- the claim's domain (cell types where its genes are detected) --------------------------
    kind = "pair" if gene_pair else ("set" if claim_genes else "global")
    domain = None
    if composition and replicate_col is not None:
        from .domain import claim_domain
        domain = claim_domain(data, type_col=composition[0], group_col=group_col, groups=groups,
                              replicate_col=replicate_col, within=within, genes=claim_genes, kind=kind,
                              min_replicates=int(prereg["min_replicates"]))
        if domain.get("error"):
            domain = None
    cells_data = data
    if domain is not None and domain["n_domain_types"]:
        cells_data = data[_g._domain_mask(data.obs, domain, within)]

    # ---- GATE 4 first: with a construct-scale SESOI its response is the SESOI in the metric's units
    g4 = None
    if signal_test is not None:
        construct = (float(sesoi) if sesoi is not None and sesoi_scale == "construct"
                     and getattr(signal_test, "kind", None) in ("coupling", "module") else None)
        g4 = _g.gate4_signal_response(metric, data, signal_test,
                                      direction=signal_direction or getattr(signal_test, "direction", None),
                                      alpha=alpha, seed=seed_of("gate4"),
                                      group_col=group_col if replicate_col else None, groups=groups,
                                      replicate_col=replicate_col, within=within,
                                      domain=domain if domain is not None and domain["n_domain_types"] else None,
                                      construct_change=construct,
                                      differential=getattr(signal_test, "kind", None) in STRUCTURE_PROBES)
    sesoi_info = _sesoi_in_metric_units(prereg, g4, signal_test)
    sesoi_metric = sesoi_info["value"]

    # ---- GATE 0 on the claim's cells, its tolerance in the metric's units ----------------------
    raw_cell = _cell_effect(metric, cells_data, group_col, groups)
    scale_candidates = [abs(v) for v in (raw_cell, sesoi_metric) if v is not None and np.isfinite(v) and v != 0]
    effect_scale = max(scale_candidates) if scale_candidates else None
    # a depth bias is not blocking when the declared correction removes depth between groups
    corrected = (("depth_downsample",) if equalizing_correction(
        data, prereg.get("estimand"), prereg.get("spikein_prefix", "ERCC-")) else ())
    g0 = _g.gate0_independence(metric, cells_data, protect_genes=gene_pair or (),
                               include_matrix_perturbations=include_matrix_perturbations,
                               seed=seed_of("gate0"), effect_scale=effect_scale, sesoi=sesoi_metric, corrected=corrected,
                               bias_tolerance=float(prereg["bias_tolerance"]))
    # change 12: a blocking bias blocks only where the groups' technical difference, the metric on GATE 0's
    # null estimated between them as the effect is, cannot be shown within GATE 0's tolerance
    if (g0.status == GateStatus.FAIL and replicate_col is not None
            and any(r.get("handling") == "blocking" for r in g0.detail.get("responses", {}).values())):
        between = technical_difference(
            metric, data, claim_genes, group_col=group_col, groups=groups, within=within,
            replicate_col=replicate_col, estimand=prereg.get("estimand"), margin=g0.detail["bias_tolerance_abs"],
            alpha=alpha, power=float(prereg["power"]), n_perm=n_perm, min_replicates=int(prereg["min_replicates"]),
            spikein_prefix=prereg.get("spikein_prefix", "ERCC-"),
            domain=domain if domain is not None and domain["n_domain_types"] else None, seed=seed_of("gate0:between"))
        g0 = _g.gate0_between_groups(g0, between)
    results.append(g0)
    g5 = None
    if pair_metric is not None and pos_pair is not None and neg_pair is not None:
        tied, why_not, tie = _controls_tie(metric, pair_metric, data, gene_pair)
        if tied:
            g5 = _g.gate5_controls(pair_metric, data, pos_pair, neg_pair, within=within, alpha=alpha,
                                   exclude=gene_pair or (), seed=seed_of("gate5"),
                                   pos_dose=float(prereg["positive_control_dose"]),
                                   direction=prereg.get("signal_direction", "increase"))
            g5.detail.update(controls_tied=True, **tie)
        else:
            g5 = GateResult(5, "Controls", GateStatus.SKIP, why_not,
                            detail=dict(controls_tied=False, **tie))
        results.append(g5)
    if g4 is not None:
        results.append(g4)
    mv = _metric_validity(g0, g5, g4)

    g1 = _g.gate1_qc_parity(data, group_col, groups, within=within, alpha=alpha, seed=seed_of("gate1"))
    results.append(g1)

    att_half = None
    resp = g0.detail.get("responses", {}).get("depth_downsample", {})
    if resp.get("classification") == "attenuation" and "note" not in resp:
        att_half = resp.get("signal_loss")

    effect_sesoi = sesoi_metric
    effect_scale_name = "observed" if sesoi_metric is not None else sesoi_scale
    # the full dossier (v0.4, the owner's decision of 2026-10-18): the effect within the domain, the
    # mixture's, the types outside the domain, the composition and every report line are estimated
    # whatever metric validity says; the verdict still follows the first FAIL (``decide``)
    effect_args = dict(
        group_col=group_col, groups=groups, within=within,
        replicate_col=replicate_col, estimand=prereg.get("estimand"), sesoi=effect_sesoi, alpha=alpha,
        power=float(prereg["power"]), min_replicates=int(prereg["min_replicates"]),
        n_perm=n_perm, spikein_prefix=prereg.get("spikein_prefix", "ERCC-"),
        qc_imbalanced=g1.status == GateStatus.WARN, attenuation_half=att_half,
        sesoi_scale=effect_scale_name, seed=seed_of("effect"), domain=domain, sesoi_label=sesoi_info.get("label"))
    eff, design = estimate_effect(metric, data, **effect_args)
    # the differential sensitivity of a metric of structure (the owner's decision of 2026-10-19): where the
    # groups' unequal responses to the same coupling would make a difference on the claim's side beyond the
    # tolerance, the effect is judged less its bound — the same estimate and test on values shifted by it
    if g4 is not None and "differential" in g4.detail:
        dom_used = domain is not None and domain["n_domain_types"]
        try:
            from .composition import technical_covariates
            gap_tech = technical_covariates(
                data, group_col, groups, replicate_col, within, domain=domain if dom_used else None,
                pool=specificity_pool, alpha=alpha, n_perm=n_perm, rng=np.random.default_rng(seed_of("differential")),
                exclude=claim_genes)
        except Exception as e:  # the premise is then not tested, and stays open
            gap_tech = dict(error=f"{type(e).__name__}: {e}")
        gap = detection_gap(gap_tech, "the domain's cells" if dom_used else "all cells")
        ds = differential_sensitivity(g4.detail["differential"], eff, prereg, sesoi_metric, gap)
        if ds.get("applies") and eff.status == "DETECTED":
            try:
                adj, _ = estimate_effect(_shifted(metric, group_col, groups[0], ds["offset"]), data, **effect_args)
                ds["adjusted"] = _adjusted_record(adj, ds["side"])
            except Exception as e:  # never silent: the effect is then not judged less the bound, and fails
                ds["adjusted"] = dict(error=f"{type(e).__name__}: {e}", passes=False)
        eff.detail["differential_sensitivity"] = ds
    if mv.status in ("FAIL", "DEGENERATE"):
        eff.detail["dossier_only"] = True
        eff = Assessment(eff.status, f"(the dossier: estimated after metric validity's {mv.status}, which "
                                     f"decides the verdict) {eff.reason}", list(eff.flags) + ["DOSSIER_ONLY"],
                         eff.detail)
    results.append(_gate2_result(eff, design))
    da = _design_adequacy(g1, design, prereg, g0.detail.get("attenuation"), sesoi_info)

    if gene_pair is not None:
        results.append(_g.gate3_raw_visibility(data, gene_pair[0], gene_pair[1], group_col, groups))

    g6 = None
    if data2 is not None and eff.status != NOT_RUN:
        g6 = _g.gate6_replication(
            metric, data2, group_col, groups, within=within, replicate_col=replicate_col,
            estimand=prereg.get("estimand"), sesoi=effect_sesoi,
            primary_effect=eff.detail.get("effect"), alpha=alpha,
            min_replicates=int(prereg["min_replicates"]), n_perm=n_perm,
            attenuation_half=att_half, seed=seed_of("gate6"), sesoi_scale=effect_scale_name)
        results.append(g6)
    rp = _replication(g6)

    autopsy = Autopsy(name, results, prereg, mv, da, eff, rp, params)
    autopsy.capture = capture_summary(eff, da, groups)
    if composition and replicate_col is not None:
        from .composition import composition_difference
        autopsy.composition = dict(
            columns={c: composition_difference(data.obs, c, group_col, groups, replicate_col, within, alpha=alpha,
                                               n_perm=n_perm, rng=np.random.default_rng(seed_of(f"composition:{c}")))
                     for c in composition})
        if domain is not None:
            autopsy.composition["domain"] = domain
    if replicate_col is not None:
        from .composition import technical_covariates
        try:
            autopsy.technical = technical_covariates(
                data, group_col, groups, replicate_col, within, domain=domain, pool=specificity_pool, alpha=alpha,
                n_perm=n_perm, rng=np.random.default_rng(seed_of("technical")))
        except Exception as e:  # a report line: reported, never silent, never blocking
            autopsy.technical = dict(error=f"{type(e).__name__}: {e}")
    if gene_pair is not None and replicate_col is not None:
        autopsy.reference = coupling_reference(metric, data, gene_pair, group_col=group_col, groups=groups,
                                               replicate_col=replicate_col, within=within, alpha=alpha,
                                               n_perm=n_perm, seed=seed_of("reference"),
                                               domain=domain if domain is not None and domain["n_domain_types"] else None)
    if specificity and replicate_col is not None and claim_genes:
        from . import specificity as _spec
        metric_of = (None if (kind == "pair" and pair_metric is None) or (kind == "set" and gene_metric is None) else
                     (lambda d, gs: pair_metric(d, gene_a=gs[0], gene_b=gs[1])) if kind == "pair" else
                     (lambda d, gs: gene_metric(d, genes=list(gs))))
        if metric_of is None:
            autopsy.specificity = dict(applicable=True, error=(
                "the metric as a function of its genes is needed: pair_metric for a pair, gene_metric for a set"))
        else:
            try:
                autopsy.specificity = _spec.check(
                    metric_of, data, claim_genes=claim_genes, kind=kind, group_col=group_col, groups=groups,
                    replicate_col=replicate_col, within=within,
                    domain=domain if domain is not None and domain["n_domain_types"] else None,
                    direction=prereg.get("direction"), alpha=alpha, rng=np.random.default_rng(seed_of("specificity")),
                    n_perm=n_perm, pool=specificity_pool, known_genes=prereg.get("known_to_differ") or (),
                    known_pairs=prereg.get("known_interactions") or (), metric=metric)
            except Exception as e:  # a report line: reported, never silent, never blocking
                autopsy.specificity = dict(applicable=True, error=f"{type(e).__name__}: {e}")
    elif specificity and replicate_col is not None and kind == "global":
        from . import specificity as _spec
        autopsy.specificity = _spec.check(None, data, claim_genes=(), kind="global", group_col=group_col,
                                          groups=groups, replicate_col=replicate_col)
    if log_path is None and not preregistered and "METRIC_AUTOPSY_LOG" not in os.environ:
        log_path = "off"  # no pre-registration: nothing to protect, nothing logged by default
    autopsy.provenance = _provenance(autopsy, data, data2, group_col, within, replicate_col, log_path, data_sha)
    return autopsy


def _hashed_columns(group_col, within, replicate_col) -> list:
    return [group_col, *within] + ([replicate_col] if replicate_col else []) + ["n_genes_by_counts", "total_counts"]


def _provenance(a: Autopsy, data, data2, group_col, within, replicate_col, log_path, data_sha=None) -> dict:
    cols = _hashed_columns(group_col, within, replicate_col)
    data_sha = data_sha or _prov.sha256_data(data, cols)
    prereg_sha = _prov.sha256_json(a.prereg)
    claim = _prov.sha256_json(dict(prereg=prereg_sha, data=data_sha, metric=a.metric_name,
                                   group_col=group_col, groups=a.params.get("groups")))
    prov = dict(schema=_prov.SCHEMA, timestamp_utc=_prov.utc_now(), environment=_prov.environment(),
                data_sha256=data_sha, data_shape=list(getattr(data.X, "shape", ())),
                data2_sha256=_prov.sha256_data(data2, cols) if data2 is not None else None,
                prereg_sha256=prereg_sha, claim_id=claim, log=None)
    path = _prov.resolve_log_path(log_path, _prov.DEFAULT_LOG)
    if path is not None:
        previous = [r for r in _prov.read_log(path) if r.get("claim_id") == claim]
        prov["log"] = dict(path=str(path), attempt=len(previous) + 1,
                           previous=[dict(timestamp_utc=r.get("timestamp_utc"), verdict=r.get("verdict"),
                                          params_sha256=r.get("params_sha256")) for r in previous])
        a.provenance = prov
        _prov.append_log(path, dict(timestamp_utc=prov["timestamp_utc"], claim_id=claim,
                                    metric=a.metric_name, data_sha256=data_sha, prereg_sha256=prereg_sha,
                                    params=a.params, params_sha256=_prov.sha256_json(a.params),
                                    verdict=decide(a),
                                    fields={k: (v.status if v else None) for k, v in a.fields().items()},
                                    engine=prov["environment"].get("metric_autopsy")))
    return prov
