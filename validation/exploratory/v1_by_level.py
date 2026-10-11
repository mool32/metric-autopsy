"""Validation v1's error rates by expression level, from its scores.json. Not a criterion.

The per-condition, per-level and per-truth rates are pre-registered secondary outcomes (v1.md
section 5: every rate "per condition, per expression level and per truth class"; section 7). This
script reads the published `scores.json` (results/panel-v1) and prints them first; then, apart,
what was decided after the results: the sums of SUPPORTED and the answers against the direction,
the pooling of the pure nulls, singling out N8 and N3 at f = 0.4, N8 at its high-level rate on
all its cards, and the comparison with the development probe p16. Nothing here changes v1's
result: the criteria are those of `validation/prereg/v1.md`.

    python validation/exploratory/v1_by_level.py --scores scores.json [--pilot pilot.json] \
        [--p16-log validation/probes/p16_variable_capture_n8.log] \
        [--by-truth-log validation/exploratory/v1_nulls_by_truth.log]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from decimal import ROUND_HALF_UP, Decimal
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from metric_autopsy.stats import clopper_pearson  # noqa: E402

SUP, OPP = "SUPPORTED", "NOT SUPPORTED: opposite direction"
REFUSAL = "REFUSAL: nuisance bias (GATE 0)"
S1_N, S1_MAX = 790, 29  # S1: at most 29 false SUPPORTED of 790 datasets per key null condition
LEVELS = ("high", "medium", "low")


def pct(x) -> str:
    """A percentage to one decimal, rounded half up (18/160 is 11.3%, as in v1.md and p16's log)."""
    d = Decimal(x.numerator) / Decimal(x.denominator) if isinstance(x, Fraction) else Decimal(repr(float(x)))
    return str((d * 100).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))


def fmt(k: int, n: int) -> str:
    lo, hi = clopper_pearson(k, n)
    return f"{k}/{n} = {pct(Fraction(k, n))}% (95% CI {pct(lo)}-{pct(hi)}%)"


def binom_cdf(k: int, n: int, p: float) -> float:
    """P(X <= k) for X ~ Binomial(n, p), summed exactly in log space."""
    def logpmf(i):
        return math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * math.log(p) + (n - i) * math.log1p(-p)
    return float(sum(math.exp(logpmf(i)) for i in range(k + 1)))


def count(entry: dict, *outcomes: str) -> int:
    return sum(int(entry["outcomes"][o]["k"]) for o in outcomes)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--scores", type=Path, required=True)
    ap.add_argument("--pilot", type=Path, default=None, help="pilot.json, for the real SESOI and dose against p16's")
    ap.add_argument("--p16-log", type=Path, default=None)
    ap.add_argument("--by-truth-log", type=Path, default=None, help="v1_nulls_by_truth.log, for N8 on a valid metric")
    args = ap.parse_args(argv)
    raw = args.scores.read_bytes()
    d = json.loads(raw)
    pl, g4 = d["per_level"], d["per_gate"]["gate4"]
    out = ["# Validation v1 by expression level, from scores.json: not a criterion",
           f"# scores.json sha256 {hashlib.sha256(raw).hexdigest()} (key {d['key']})", "",
           "## Pre-registered secondary outcomes: rates per condition, level and truth (v1.md sections 5 and 7)"]

    # the cards by the metric's truth and level: per_gate.gate4 tabulates GATE 4's outcome by
    # truth:level for every card with a report (the run had no crash, so for all of them)
    truths = sorted({k.split(":")[0] for k in g4})
    table = {(t, lev): sum(g4.get(f"{t}:{lev}", {}).values()) for t in truths for lev in LEVELS}
    out.append("cards by the metric's truth (the oracle) and level:")
    out.append("  " + "truth".ljust(10) + "".join(lev.rjust(8) for lev in LEVELS) + "   total")
    for t in truths:
        row = [table[(t, lev)] for lev in LEVELS]
        out.append("  " + t.ljust(10) + "".join(str(x).rjust(8) for x in row) + str(sum(row)).rjust(8))
    tot = [sum(table[(t, lev)] for t in truths) for lev in LEVELS]
    out.append("  " + "all".ljust(10) + "".join(str(x).rjust(8) for x in tot) + str(sum(tot)).rjust(8))
    valid = {lev: table.get(("valid", lev), 0) for lev in LEVELS}
    out.append(f"valid metric: {valid['high']} cards at the high level, {valid['medium']} at medium, {valid['low']} at low")

    # every SUPPORTED and every answer against the direction, by level
    for lev in LEVELS:
        es = [v for k, v in pl.items() if k.endswith(":" + lev)]
        n = sum(int(v["n"]) for v in es)
        out.append(f"  {lev}: SUPPORTED {sum(count(v, SUP) for v in es)} of {n} cards, "
                   f"against the direction {sum(count(v, OPP) for v in es)}")

    def high(key):
        return pl[f"{key}:high"]

    n8, n3 = high("N8:beta(2,2)"), high("N3:f=0.4")
    n8_all = sum(int(v["n"]) for k, v in pl.items() if k.startswith("N8:"))
    pure = [("N1", high("N1:null")), ("N5", high("N5:sham")), ("N7", high("N7:mice"))]
    out.append(f"N8, high level ({n8['n']} of its {n8_all} cards): false SUPPORTED {fmt(count(n8, SUP), n8['n'])} "
               f"(nominal 2.5%)")
    out.append(f"N3 f=0.4, high level: false SUPPORTED {fmt(count(n3, SUP), n3['n'])}")
    for name, e in pure:
        out.append(f"{name}, high level: false SUPPORTED {fmt(count(e, SUP), e['n'])}")
    n7 = high("N7:mice")
    out.append(f"GATE 0's refusals, N7 high level: {fmt(count(n7, REFUSAL), n7['n'])}")

    out += ["", "## Decided after the results",
            "the sums of SUPPORTED and the answers against the direction, high level (nominal 5%):"]
    for name, e in (("N8", n8), ("N3 f=0.4", n3), *pure):
        out.append(f"  {name}: {fmt(count(e, SUP, OPP), e['n'])}")
    k_pure, n_pure = sum(count(e, SUP) for _, e in pure), sum(int(e["n"]) for _, e in pure)
    out.append(f"the pure nulls N1, N5 and N7 pooled, high level: false SUPPORTED {fmt(k_pure, n_pure)}")
    out.append("N8 and N3 at f = 0.4 were singled out after the results, as the two largest high-level rates "
               "(the twelve null cells with Holm's correction: v1_null_cells.log)")
    p = count(n8, SUP) / n8["n"]
    out.append(f"N8 at its high-level rate {100 * p:.2f}% on all {S1_N} cards: {S1_N * p:.1f} false SUPPORTED expected, "
               f"S1 allows {S1_MAX}; P(S1 passes on N8) = {binom_cdf(S1_MAX, S1_N, p):.2f} "
               "(replaced in the text by the criteria's resolution on the cards at risk, v1_null_cells.log)")

    if args.p16_log and args.p16_log.exists():
        text = args.p16_log.read_text()
        m = re.search(r"N8 beta\(2,2\).*?\n\s*false SUPPORTED by the engine (\d+)/(\d+).*?pairs by level (\{[^}]*\})",
                      text, flags=re.S)
        if m:
            k16, n16 = int(m.group(1)), int(m.group(2))
            levels16 = json.loads(m.group(3).replace("'", '"'))
            v16 = levels16.get("high", 0) + levels16.get("medium", 0)  # dry_pilot: valid at high and medium
            import inspect
            sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "validation" / "prereg"))
            import simulate
            dflt = inspect.signature(simulate.dry_pilot).parameters
            real = ""
            if args.pilot and args.pilot.exists():
                pj = json.loads(args.pilot.read_text())
                real = (f", against the real B1's high level {pj['sesoi']['B1']['high']:g} and "
                        f"{pj['saturation_dose']['B1']['high']:g}")
            out += ["the comparison with the development probe p16 (simulated B1; simulate.dry_pilot's stand-in truth: "
                    "valid at the high and medium levels, blind at the low;",
                    f"  its SESOI {dflt['sesoi'].default:g} and dose {dflt['dose'].default:g} at every level{real}):",
                    f"  p16 N8: false SUPPORTED {fmt(k16, n16)} over all its datasets ({levels16}); "
                    "it does not record on which they fell"]
            k8_all = sum(count(v, SUP) for k, v in pl.items() if k.startswith("N8:"))
            ratios = [f"if all {k16} fell on its {v16} datasets with a valid metric: {fmt(k16, v16)}, "
                      f"{(k16 / v16) / p:.1f} times N8's high-level rate"]
            if args.by_truth_log and args.by_truth_log.exists():
                mv = re.search(r"N8:.*?valid\s+high\s+: false SUPPORTED (\d+)/(\d+)", args.by_truth_log.read_text(), flags=re.S)
                if mv:
                    kv, nv = int(mv.group(1)), int(mv.group(2))
                    ratios[0] += f" and {(k16 / v16) / (kv / nv):.1f} times N8's {kv}/{nv} on its cards with a valid metric"
            ratios.append(f"all levels against all levels: {pct(Fraction(k16, n16))}% against N8's "
                          f"{k8_all}/{n8_all} = {pct(Fraction(k8_all, n8_all))}%, {(k16 / n16) / (k8_all / n8_all):.1f} times")
            out += [f"  {r}" for r in ratios]
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
