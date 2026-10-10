"""Validation v1's error rates by expression level, from its scores.json. Exploratory, not a criterion.

v1's criteria count errors over every card of a condition. On the real backgrounds the gene
pairs of the low and medium expression levels are blind to the panel's metric, and the engine
says "metric invalid" there. So a false SUPPORTED can occur only at the high level. This script
reads the published `scores.json` (results/panel-v1) and prints the rates at the high level,
where the errors could occur. Nothing here changes v1's result: the criteria are those of
`validation/prereg/v1.md`, and this is a secondary reading of them, decided after the results.

    python validation/exploratory/v1_by_level.py --scores scores.json \
        [--p16-log validation/probes/p16_variable_capture_n8.log]
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
    ap.add_argument("--p16-log", type=Path, default=None)
    args = ap.parse_args(argv)
    raw = args.scores.read_bytes()
    d = json.loads(raw)
    pl, g4 = d["per_level"], d["per_gate"]["gate4"]
    out = ["# Validation v1 by expression level: exploratory, not a criterion",
           f"# scores.json sha256 {hashlib.sha256(raw).hexdigest()} (key {d['key']})", ""]

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

    n8 = high("N8:beta(2,2)")
    n8_all = sum(int(v["n"]) for k, v in pl.items() if k.startswith("N8:"))
    out += ["", f"N8: {n8['n']} of its {n8_all} cards at the high level, where a false SUPPORTED can occur",
            f"N8, high level: false SUPPORTED {fmt(count(n8, SUP), n8['n'])} (nominal 2.5%); "
            f"SUPPORTED or against the direction {fmt(count(n8, SUP, OPP), n8['n'])} (nominal 5%)"]
    n3 = high("N3:f=0.4")
    out.append(f"N3 f=0.4, high level: false SUPPORTED {fmt(count(n3, SUP), n3['n'])}; "
               f"SUPPORTED or against the direction {fmt(count(n3, SUP, OPP), n3['n'])}")
    pure = [("N1", high("N1:null")), ("N5", high("N5:sham")), ("N7", high("N7:mice"))]
    for name, e in pure:
        out.append(f"{name}, high level: SUPPORTED or against the direction {fmt(count(e, SUP, OPP), e['n'])}")
    k_pure, n_pure = sum(count(e, SUP) for _, e in pure), sum(int(e["n"]) for _, e in pure)
    out.append(f"the pure nulls N1, N5 and N7, high level: false SUPPORTED {fmt(k_pure, n_pure)}")

    p = count(n8, SUP) / n8["n"]
    out += ["", f"N8 at its high-level rate {100 * p:.2f}% on all {S1_N} cards: {S1_N * p:.1f} false SUPPORTED expected, "
            f"S1 allows {S1_MAX}; P(S1 passes on N8) = {binom_cdf(S1_MAX, S1_N, p):.2f}"]
    n7 = high("N7:mice")
    out.append(f"GATE 0's refusals, N7 high level: {fmt(count(n7, REFUSAL), n7['n'])}")

    if args.p16_log and args.p16_log.exists():
        m = re.search(r"N8 beta\(2,2\).*?\n\s*false SUPPORTED by the engine (\d+)/(\d+)", args.p16_log.read_text())
        if m:
            k16, n16 = int(m.group(1)), int(m.group(2))
            out.append(f"dev probe p16 (simulated B1, informative pairs at every level): N8 false SUPPORTED "
                       f"{fmt(k16, n16)}; {(k16 / n16) / p:.1f} times the real high-level rate")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
