"""Assemble a gate-by-gate autopsy and a single verdict.

The verdict rule is deliberately blunt and matches the research flow's "no rescue
language": the autopsy is decided by the *first blocking gate*. A metric that fails GATE 0
is dead regardless of what GATE 5 says. If every auto gate passes but judgment gates
(4, 7) are unresolved, the verdict is INCONCLUSIVE pending the analyst's ruling — never a
premature PASS.
"""
from __future__ import annotations

import functools
from dataclasses import dataclass, field
from typing import Sequence

from .core import GateResult, GateStatus
from . import gates as _g


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


@dataclass
class Autopsy:
    metric_name: str
    results: list[GateResult] = field(default_factory=list)
    prereg: dict = field(default_factory=dict)

    @property
    def verdict(self) -> str:
        blocking = [r for r in self.results if r.blocking]
        if blocking:
            first = min(blocking, key=lambda r: r.gate)
            return f"FAIL — died at GATE {first.gate} ({first.name})"
        judged = [r for r in self.results if r.status == GateStatus.JUDGMENT]
        pending = self.prereg.get("judgment_pending", bool(judged))
        if pending:
            return "INCONCLUSIVE — auto gates clear; judgment gates (4, 7) unresolved"
        passed = sum(r.status == GateStatus.PASS for r in self.results)
        return f"PASS — cleared {passed} auto gates; treat as provisional until replicated"

    def to_markdown(self) -> str:
        lines = [f"# Metric autopsy — {_cell(self.metric_name)}", ""]
        prereg_items = {k: v for k, v in self.prereg.items() if k != "judgment_pending"}
        if prereg_items:
            lines += ["## Pre-registration", ""]
            for k, v in prereg_items.items():
                lines.append(f"- **{_cell(k)}:** {_cell(v)}")
            lines.append("")
        lines += ["## Gates", "", "| Gate | Name | Status | Finding |", "|---|---|---|---|"]
        for r in sorted(self.results, key=lambda r: r.gate):
            lines.append(f"| {r.gate} | {_cell(r.name)} | **{r.status.value}** | {_cell(r.message)} |")
        lines += ["", "## Verdict", "", f"**{self.verdict}**", ""]
        return "\n".join(lines)


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
) -> Autopsy:
    """Run the auto gates in order and collect them into an Autopsy.

    `metric(data) -> float` is the bound metric (params fixed). Pairwise gates (3, 5) need
    the unbound `pair_metric(data, *, gene_a, gene_b)` plus the relevant gene pairs.
    With `stop_on_first_fail` the sequence halts at the first blocking gate — including the
    pairwise gates — mirroring the manual protocol.
    """
    name = _metric_name(metric)
    prereg = prereg or {}
    results: list[GateResult] = []

    # Whole-matrix metrics (no bound gene pair) are probed with the gene_subsample
    # perturbation in GATE 0; pairwise metrics are not (it is irrelevant to two fixed genes).
    if include_matrix_perturbations is None:
        include_matrix_perturbations = gene_pair is None

    def _record(r: GateResult) -> bool:
        results.append(r)
        return not (stop_on_first_fail and r.blocking)

    if not _record(_g.gate0_independence(
            metric, data, protect_genes=gene_pair or (),
            include_matrix_perturbations=include_matrix_perturbations)):
        return Autopsy(name, results, prereg)
    if not _record(_g.gate1_qc_parity(data, group_col, groups, within=within)):
        return Autopsy(name, results, prereg)
    if not _record(_g.gate2_ngenes_matching(metric, data, group_col, groups)):
        return Autopsy(name, results, prereg)

    if gene_pair is not None:
        if not _record(_g.gate3_raw_visibility(data, gene_pair[0], gene_pair[1], group_col, groups)):
            return Autopsy(name, results, prereg)
    if pair_metric is not None and pos_pair is not None and neg_pair is not None:
        if not _record(_g.gate5_controls(pair_metric, data, pos_pair, neg_pair, within=within)):
            return Autopsy(name, results, prereg)
    if data2 is not None:
        if not _record(_g.gate6_replication(metric, data2, group_col, groups, within=within)):
            return Autopsy(name, results, prereg)

    return Autopsy(name, results, prereg)
