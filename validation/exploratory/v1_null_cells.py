"""Validation v1's false-SUPPORTED cells read on the cards at risk, the twelve null cells of the
high expression level, the design effects, and the S4 part of the prediction on N8. Read after
the results; it changes nothing in v1's result.

Each section says what was pre-registered and what was decided after the results:

1. The criteria's resolution on the cards at risk (decided after the results; the pilot showed the
   mix before the key: oc_pilot.log, "valid metric on 2453"). A false SUPPORTED needs a metric
   that responds: on a blind metric GATE 4 stops the analysis first (0 of the 537 blind cards of
   N8 and 0 of the 6,988 medium- and low-level cards were SUPPORTED). A card is at risk if its
   metric is valid or ambiguous by the oracle's truth. For every cell of S1 and S2: its cards, the
   cards at risk, the limit, and on the cards at risk the rate of false SUPPORTED at which the
   cell fails with probability 0.5 and 0.95, and the probability that it passes at 5% and 7.5%.
   On a useless metric (N6a, N6c), and on the real effects' cards where SUPPORTED is an error (all
   with a blind metric: on a valid or ambiguous one SUPPORTED is allowed for a real effect), a false
   SUPPORTED needs GATE 4 or 5 to pass the metric first; no card is at risk, and the same numbers
   are given for that compound error over all cards.
2. The twelve null cells at the high level: k/n and the 95% Clopper-Pearson interval are the
   pre-registered per-condition, per-level rates (v1.md section 5); the one-sided exact binomial
   p against the nominal 2.5% and Holm's correction over the twelve are decided after the results.
3. The design effects (pre-registered, v1.md sections 5 and 7), from scores.json.
4. The S4 part of the prediction on N8 (v1.md section 6: "the answers against the direction on N8
   may fail S4's cell on the nulls as well").

    python validation/exploratory/v1_null_cells.py --results RESULTS --pilot pilot.json --scores scores.json

RESULTS is results/panel-v1 with its reports unpacked (as for v1_nulls_by_truth.py).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation" / "prereg"))
sys.path.insert(0, str(ROOT / "validation" / "exploratory"))
import oc  # noqa: E402
import panel as P  # noqa: E402
import score  # noqa: E402
from v1_by_level import fmt  # noqa: E402  (k/n = rate (95% Clopper-Pearson CI), rounded half up)

NOMINAL = 0.025  # a false SUPPORTED's designed size, alpha/2
RISK = ("valid", "ambiguous")
USELESS = ("useless", "constant")
NULL_CELLS = (("N1", "null"), ("N2", "c=0.5"), ("N2", "c=0.9"), ("N2", "c=0.7"), ("N2", "c=0.3"),
              ("N3", "f=0.1"), ("N3", "f=0.2"), ("N3", "f=0.4"), ("N4", "3v3"), ("N5", "sham"),
              ("N7", "mice"), ("N8", "beta(2,2)"))


def binom_sf(k: int, n: int, p: float) -> float:
    """P(K > k) for K ~ Binomial(n, p)."""
    from scipy import stats
    return float(stats.binom.sf(k, n, p))


def rate_at(n: int, limit: int, prob: float) -> float:
    """The smallest rate p on n cards at which more than `limit` errors occur with probability
    >= `prob` (bisection on the exact binomial)."""
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if binom_sf(limit, n, mid) >= prob:
            hi = mid
        else:
            lo = mid
    return hi


def holm(ps: list[float]) -> list[float]:
    """Holm's step-down adjusted p values, in the input order."""
    m = len(ps)
    order = sorted(range(m), key=lambda i: ps[i])
    adj, run = [0.0] * m, 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (m - rank) * ps[i]))
        adj[i] = run
    return adj


def resolution_line(name: str, rs: list[dict], limit: int) -> str:
    n, k = len(rs), sum(r["sup_error"] for r in rs)
    by = {t: sum(r["truth"] == t for r in rs) for t in ("valid", "ambiguous", "blind", "useless", "constant")}
    risk = [r for r in rs if r["truth"] in RISK]
    head = (f"  {name:7s} n {n:5d}, false SUPPORTED {k:3d}, limit {limit:3d} | valid {by['valid']:4d}, "
            f"ambiguous {by['ambiguous']:3d}, blind {by['blind']:4d}")
    if not risk:
        useless = sum(by[t] for t in USELESS)
        what = "a useless metric" if useless == n else "a blind metric on every card"
        return (head + f", useless {useless:4d} | at risk 0 ({what}: a false SUPPORTED needs GATE 4 or 5 to pass "
                f"the metric first); for that error over all {n} cards: the limit is {100 * limit / n:.1f}% of them; "
                f"fails with P 0.5 at {100 * rate_at(n, limit, 0.5):.1f}%, with P 0.95 at "
                f"{100 * rate_at(n, limit, 0.95):.1f}%")
    nr, kr = len(risk), sum(r["sup_error"] for r in risk)
    return (head + f" | at risk {nr} (false SUPPORTED {kr}); the limit is {100 * limit / nr:.1f}% of them; "
            f"fails with P 0.5 at {100 * rate_at(nr, limit, 0.5):.1f}%, with P 0.95 at "
            f"{100 * rate_at(nr, limit, 0.95):.1f}%; passes at 5% with P {1 - binom_sf(limit, nr, 0.05):.3f}, "
            f"at 7.5% with P {1 - binom_sf(limit, nr, 0.075):.3f}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--results", type=Path, required=True, help="results/panel-v1, its reports unpacked")
    ap.add_argument("--pilot", type=Path, required=True)
    ap.add_argument("--scores", type=Path, required=True)
    args = ap.parse_args(argv)
    manifest = json.loads((args.results / "manifest.json").read_text())
    pilot_raw, scores_raw = args.pilot.read_bytes(), args.scores.read_bytes()
    pilot, scores = json.loads(pilot_raw), json.loads(scores_raw)
    key = manifest["key"]
    reports, _ = score.read_results(args.results, P.check_key(key))
    entries = P.assign(key, pilot.get("dropped", ()), pilot.get("pool_size"))
    rows = score.card_rows(entries, reports, pilot)
    crit = scores["criteria"]
    out = ["# Validation v1: the false-SUPPORTED cells on the cards at risk, the twelve high-level null cells,",
           "# the design effects and the S4 part of the prediction on N8 (read after the results; v1 unchanged)",
           f"# key {key}; pilot.json sha256 {hashlib.sha256(pilot_raw).hexdigest()}; "
           f"scores.json sha256 {hashlib.sha256(scores_raw).hexdigest()}; {len(rows)} cards", ""]

    out += ["## 1. The criteria's resolution on the cards at risk (decided after the results)",
            "A card is at risk if its metric is valid or ambiguous (the oracle's truth); the rates below are rates",
            "of false SUPPORTED among the cards at risk, with no error on the others.", "S1 (limit per key null):"]
    for c in P.CONDITIONS:
        if c.key:
            rs = [r for r in rows if r["condition"] == c.name and r["variant"] == c.key_variant]
            out.append(resolution_line(c.name, rs, int(crit["S1"]["conditions"][c.name]["max_allowed"])))
    out.append("S2 (limit per group of the cards outside S1):")
    for cell, rule in oc.error_cells(rows).items():
        if cell.startswith("S2:"):
            rs = [rows[i] for i in rule["members"]]
            out.append(resolution_line(cell[3:], rs, int(crit["S2"]["cells"][cell]["max_allowed"])))
    blind_n8 = [r for r in rows if r["condition"] == "N8" and r["truth"] == "blind"]
    blind_ml = [r for r in rows if r["level"] in ("medium", "low") and r["truth"] == "blind"]
    out += [f"on blind metrics: N8 {sum(r['outcome'] == P.SUPPORTED for r in blind_n8)} SUPPORTED of {len(blind_n8)}; "
            f"medium and low levels {sum(r['outcome'] == P.SUPPORTED for r in blind_ml)} of {len(blind_ml)}", ""]

    out += ["## 2. The twelve null cells at the high level",
            "k/n and the 95% CI: pre-registered (v1.md section 5); the one-sided p against 2.5% and Holm's correction",
            "over the twelve: decided after the results. Only the cards on which a false SUPPORTED is an error count.",
            f"  {'cell':18s} {'false SUPPORTED (95% CI)':34s} {'valid/ambiguous/blind':22s} {'p (one-sided)':>14s} "
            f"{'Holm':>8s}"]
    cells = []
    for cond, var in NULL_CELLS:
        rs = [r for r in rows if r["condition"] == cond and r["variant"] == var and r["level"] == "high"
              and r["sup_error_possible"]]
        k, n = sum(r["sup_error"] for r in rs), len(rs)
        p = binom_sf(k - 1, n, NOMINAL)  # P(K >= k)
        mix = "/".join(str(sum(r["truth"] == t for r in rs)) for t in ("valid", "ambiguous", "blind"))
        cells.append((f"{cond} {var}", k, n, p, mix))
    for (name, k, n, p, mix), h in zip(cells, holm([c[3] for c in cells])):
        out.append(f"  {name:18s} {fmt(k, n):34s} {mix:22s} {p:14.4f} {h:8.3f}")
    n8v = [r for r in rows if r["condition"] == "N8" and r["level"] == "high" and r["truth"] == "valid"]
    out += [f"  N8 on its high-level cards with a valid metric (decided after the results): "
            f"{fmt(sum(r['sup_error'] for r in n8v), len(n8v))}", ""]

    out.append("## 3. The design effects (pre-registered: v1.md sections 5 and 7; from scores.json)")
    sd = scores["secondary"]["shared_donors"]
    for cell, v in sorted(sd.items(), key=lambda kv: -kv[1]["deff"]):
        if v["deff"] > 1.5:
            out.append(f"  {cell:16s} design effect {v['deff']:.2f}: rate {v['rate']:.3f}, 95% CI at n/deff "
                       f"{v['ci95'][0]:.3f}-{v['ci95'][1]:.3f} (n {v['n']})")
    sup = {c: v for c, v in sd.items() if c.startswith(("S1:", "S2:"))}
    worst = max(sup, key=lambda c: sup[c]["deff"])
    out.append(f"  the false-SUPPORTED cells (S1, S2): largest design effect {sup[worst]['deff']:.4f} ({worst})")
    for c, v in scores["secondary"]["s1_at_realized_deff"]["conditions"].items():
        out.append(f"  S1 {c} at its design effect {v['deff']:.4f}: P(pass | sound) {v['p_pass_sound']:.4f}, "
                   f"P(pass | doubled) {v['p_pass_doubled']:.4f}")
    out += ["  (N2's design effect is 1 + 2e-16; oc.error_rule_deff then takes a beta-binomial in which scipy",
            "   overflows and reads 1.0 / 1.0; at a design effect of 1 the values are N1's)", ""]

    out.append("## 4. The S4 part of the prediction on N8 (v1.md section 6)")
    s4 = crit["S4"]["cells"]["S4:null:opposite"]
    n8 = [r for r in rows if r["condition"] == "N8" and r["opp_null_counted"]]
    out += [f"  S4's null cell, NOT SUPPORTED against the direction: {s4['k']} of {s4['n']}, limit {s4['max_allowed']} "
            f"(passed {s4['passed']})",
            f"  of which on N8: {sum(r['opp_null'] for r in n8)} of {len(n8)} N8 cards in the cell "
            f"({sum(r['opp_null'] for r in n8 if r['level'] == 'high')} at the high level)"]
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
