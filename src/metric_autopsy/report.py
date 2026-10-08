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
  attenuated SESOI. Statuses: ADEQUATE, CORRECTED, INSUFFICIENT_REPLICATION,
  UNIDENTIFIABLE; flags: UNDERPOWERED, PARAMETRIC_ONLY, ATTENUATION, POWER_NOT_ASSESSED.
* ``effect`` — the corrected difference at the replicate level: DETECTED,
  NO_DETECTABLE_EFFECT (TOST against the pre-registered SESOI), INCONCLUSIVE, NOT_ESTIMABLE.
* ``replication`` — REPLICATED / NOT_REPLICATED / INCONCLUSIVE / NOT_RUN (GATE 6).

No rescue language: a FAIL in an earlier field is never softened by a later one, and with
``stop_on_first_fail`` an invalid metric stops the analysis before the effect is estimated.
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
    return p


def delta_min_of(prereg: dict) -> float | None:
    """delta_min, the smallest response to an injected signal that matters (on the metric's
    scale): pre-registered as ``delta_min``, by default ``DELTA_MIN_FRACTION`` (0.5) x SESOI.
    Without either, GATE 4 and GATE 5 cannot show that a metric is blind (decided 2026-10-08)."""
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
          "unidentifiable", "effect_not_evaluated", "explained_by_depth", "insufficient_replication",
          "absence_untested_metric", "no_detectable_effect", "effect_inconclusive",
          "detected_untested_metric", "bias_unsized", "no_direction", "opposite_direction",
          "not_replicated", "judgment_pending", "replicated", "provisional")


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
    tail = f" [{'; '.join(notes)}]" if notes else ""

    if mv.status == "DEGENERATE":
        return f"DEGENERATE METRIC — {mv.reason}", "degenerate_metric"
    if mv.status == "FAIL":
        gate = mv.detail.get("failed_gate")
        return (f"NOT SUPPORTED — metric invalid: {mv.reason}",
                f"metric_invalid_gate{gate}" if gate in (0, 4, 5) else "metric_invalid")
    if da.status == "UNIDENTIFIABLE":
        return f"UNIDENTIFIABLE — {da.reason}", "unidentifiable"
    if ef.status == NOT_RUN:
        return f"INCONCLUSIVE — effect not evaluated: {ef.reason}", "effect_not_evaluated"
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
    if ef.status != "DETECTED":
        return f"INCONCLUSIVE — {ef.reason}{tail}", "effect_inconclusive"
    if mv.status != "PASS":
        return ("INCONCLUSIVE — effect detected, but the metric's response to signal is untested "
                "(supply a positive control or an injected signal)" + tail, "detected_untested_metric")
    if "BIAS_UNSIZED" in mv.flags:
        return ("INCONCLUSIVE — effect detected, but a nuisance bias of the metric "
                f"({', '.join(mv.detail.get('unsized_bias', []))}) cannot be sized against the claim "
                "without a pre-registered SESOI: declare a SESOI (the smallest effect that matters, on "
                "the metric's scale) in the pre-registration, and the bias is judged against "
                "bias_tolerance x SESOI" + tail, "bias_unsized")
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
    if rp.status == "NOT_REPLICATED":
        return f"NOT SUPPORTED — the effect did not replicate: {rp.reason}{tail}", "not_replicated"
    if a.judgment_pending:
        return (f"INCONCLUSIVE — effect detected; judgment gates 4 and 7 are unresolved{tail}",
                "judgment_pending")
    if rp.status == "REPLICATED":
        return f"SUPPORTED — replicated{tail}", "replicated"
    return f"SUPPORTED (provisional until replicated){tail}", "provisional"


def decide(a: "Autopsy") -> str:
    """The verdict of ``decide_cause`` (the single verdict rule)."""
    return decide_cause(a)[0]


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
        ))

    def to_json(self, **kw) -> str:
        kw.setdefault("indent", 2)
        return json.dumps(self.to_dict(), allow_nan=False, **kw)

    def save_json(self, path) -> None:
        with open(path, "w") as fh:
            fh.write(self.to_json())


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


def _design_adequacy(g1: GateResult, design: dict, prereg: dict, attenuation: dict | None = None) -> Assessment:
    flags, notes = [], list(design.get("notes", []))
    attenuation = dict(attenuation or {})
    tier = design.get("tier")
    if tier == "parametric":
        flags.append("PARAMETRIC_ONLY")
    if design.get("underpowered"):
        flags.append("UNDERPOWERED")
    elif design.get("underpowered") is None and prereg.get("sesoi") is None and tier not in (None, "insufficient"):
        flags.append("POWER_NOT_ASSESSED")
    if attenuation:
        flags.append("ATTENUATION")
    detail = dict(qc=g1.status.value, correction=design.get("correction"), tier=tier,
                  replication=design.get("replication"), power=design.get("power"),
                  attenuation=attenuation, equalization=design.get("equalization"), notes=notes)
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
        att = ("; attenuation (GATE 0) " + ", ".join(f"{k} −{v:.0%}" for k, v in attenuation.items())
               + (f" -> λ={power['lambda_']:.2f} at the analysed depth" if power.get("lambda_") is not None else ""))
    pw = ""
    if power.get("effective_sesoi") is not None and np.isfinite(power.get("mde", np.nan)):
        pw = (f"; MDE {power['mde']:.4g} vs attenuated SESOI {power['effective_sesoi']:.4g} "
              f"(SESOI {power['sesoi']:g} x λ {power['lambda_']:.2f})")
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
    return GateResult(2, "Estimand-dependent correction", GateStatus.PASS, txt, d)


def _replication(g6: GateResult | None) -> Assessment:
    if g6 is None:
        return Assessment(NOT_RUN, "no independent dataset supplied")
    status = g6.detail.get("replication", NOT_RUN)
    return Assessment(status, g6.message, list(g6.detail.get("effect_flags", [])), g6.detail)


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
) -> Autopsy:
    """Run the gates and assemble the four-field verdict.

    `metric(data) -> float` is the bound metric. Controls (GATE 5) need the unbound
    `pair_metric(data, *, gene_a, gene_b)`, `gene_pair` and both control pairs; they count only if
    `pair_metric` on `gene_pair` gives the metric's value on the data (else GATE 5 is SKIP: the
    controls would test another metric). `replicate_col` names the
    biological replicate (mouse, donor, plate): without it there is no effect verdict.
    `signal_test` is an ``injected_signal`` constructor result. `prereg` carries
    ``estimand`` ('composition' | 'content'), ``sesoi`` (construct scale), ``min_replicates``,
    ``alpha``, ``judgment_pending`` and free-text commitments. Its hash is recorded, and a run
    with a pre-registration is appended to the run log (`log_path`, else
    ``$METRIC_AUTOPSY_LOG``, else ``metric_autopsy_runs.jsonl``; ``"off"`` disables).
    With `stop_on_first_fail`, an invalid or degenerate metric stops the analysis before the
    effect is estimated.
    """
    name = _metric_name(metric)
    preregistered = bool(prereg)
    prereg = normalize_prereg(prereg)
    groups = tuple(groups)
    within = list(within)
    alpha = float(prereg["alpha"])
    sesoi = prereg.get("sesoi")
    delta_min = delta_min_of(prereg)
    results: list[GateResult] = []
    params = dict(metric=name, group_col=group_col, groups=list(groups), within=within,
                  gene_pair=list(gene_pair) if gene_pair else None,
                  pair_metric=_metric_name(pair_metric) if pair_metric else None,
                  pos_pair=list(pos_pair) if pos_pair else None,
                  neg_pair=list(neg_pair) if neg_pair else None, replicate_col=replicate_col,
                  signal_test=getattr(signal_test, "description", None) if signal_test else None,
                  data2=data2 is not None, include_matrix_perturbations=include_matrix_perturbations,
                  stop_on_first_fail=stop_on_first_fail, seed=seed, n_perm=n_perm)

    if include_matrix_perturbations is None:
        include_matrix_perturbations = gene_pair is None

    from .effect import _cell_effect, assess_design, equalizing_correction, estimate_effect, public_design

    raw_cell = _cell_effect(metric, data, group_col, groups)
    scale_candidates = [abs(v) for v in (raw_cell, sesoi) if v is not None and np.isfinite(v) and v != 0]
    effect_scale = max(scale_candidates) if scale_candidates else None
    # a depth bias is not blocking when the declared correction removes depth between groups
    corrected = (("depth_downsample",) if equalizing_correction(
        data, prereg.get("estimand"), prereg.get("spikein_prefix", "ERCC-")) else ())

    g0 = _g.gate0_independence(metric, data, protect_genes=gene_pair or (),
                               include_matrix_perturbations=include_matrix_perturbations,
                               seed=seed, effect_scale=effect_scale, sesoi=sesoi, corrected=corrected,
                               bias_tolerance=float(prereg["bias_tolerance"]))
    results.append(g0)
    g5 = None
    if pair_metric is not None and pos_pair is not None and neg_pair is not None:
        tied, why_not, tie = _controls_tie(metric, pair_metric, data, gene_pair)
        if tied:
            g5 = _g.gate5_controls(pair_metric, data, pos_pair, neg_pair, within=within, alpha=alpha,
                                   exclude=gene_pair or (), seed=seed,
                                   pos_dose=float(prereg["positive_control_dose"]), delta_min=delta_min,
                                   direction=prereg.get("signal_direction", "increase"))
            g5.detail.update(controls_tied=True, **tie)
        else:
            g5 = GateResult(5, "Controls", GateStatus.SKIP, why_not,
                            detail=dict(controls_tied=False, **tie))
        results.append(g5)
    g4 = None
    if signal_test is not None:
        g4 = _g.gate4_signal_response(metric, data, signal_test,
                                      direction=prereg.get("signal_direction", "increase"),
                                      delta_min=delta_min, alpha=alpha, seed=seed)
        results.append(g4)
    mv = _metric_validity(g0, g5, g4)

    g1 = _g.gate1_qc_parity(data, group_col, groups, within=within, alpha=alpha, seed=seed)
    results.append(g1)

    att_half = None
    resp = g0.detail.get("responses", {}).get("depth_downsample", {})
    if resp.get("classification") == "attenuation" and "note" not in resp:
        att_half = resp.get("signal_loss")

    if stop_on_first_fail and mv.status in ("FAIL", "DEGENERATE"):
        eff = Assessment(NOT_RUN, f"metric_validity is {mv.status}; the effect is not estimated "
                                  "(stop_on_first_fail)")
        # design adequacy does not depend on the metric: still assess identifiability and replication
        design = public_design(assess_design(
            data, group_col=group_col, groups=groups, within=within, replicate_col=replicate_col,
            estimand=prereg.get("estimand"), min_replicates=int(prereg["min_replicates"]),
            spikein_prefix=prereg.get("spikein_prefix", "ERCC-"),
            qc_imbalanced=g1.status == GateStatus.WARN))
    else:
        eff, design = estimate_effect(
            metric, data, group_col=group_col, groups=groups, within=within,
            replicate_col=replicate_col, estimand=prereg.get("estimand"), sesoi=sesoi, alpha=alpha,
            power=float(prereg["power"]), min_replicates=int(prereg["min_replicates"]),
            n_perm=n_perm, spikein_prefix=prereg.get("spikein_prefix", "ERCC-"),
            qc_imbalanced=g1.status == GateStatus.WARN, attenuation_half=att_half, seed=seed)
        results.append(_gate2_result(eff, design))
    da = _design_adequacy(g1, design, prereg, g0.detail.get("attenuation"))

    if gene_pair is not None:
        results.append(_g.gate3_raw_visibility(data, gene_pair[0], gene_pair[1], group_col, groups))

    g6 = None
    if data2 is not None and eff.status != NOT_RUN:
        g6 = _g.gate6_replication(
            metric, data2, group_col, groups, within=within, replicate_col=replicate_col,
            estimand=prereg.get("estimand"), sesoi=sesoi,
            primary_effect=eff.detail.get("effect"), alpha=alpha,
            min_replicates=int(prereg["min_replicates"]), n_perm=n_perm,
            attenuation_half=att_half, seed=seed)
        results.append(g6)
    rp = _replication(g6)

    autopsy = Autopsy(name, results, prereg, mv, da, eff, rp, params)
    if log_path is None and not preregistered and "METRIC_AUTOPSY_LOG" not in os.environ:
        log_path = "off"  # no pre-registration: nothing to protect, nothing logged by default
    autopsy.provenance = _provenance(autopsy, data, data2, group_col, within, replicate_col, log_path)
    return autopsy


def _provenance(a: Autopsy, data, data2, group_col, within, replicate_col, log_path) -> dict:
    cols = [group_col, *within] + ([replicate_col] if replicate_col else []) + ["n_genes_by_counts", "total_counts"]
    data_sha = _prov.sha256_data(data, cols)
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
