"""A separate outcome "effect below the SESOI" applied after the fact to validation v1's published
reports (exploratory, post hoc; the third finding of the real-data pilot, v0.4 spec change 14). v1's
result and the engine are unchanged.

In the pilot (validation/exploratory/real_pilot.log, claims 3a2 and 3c3) the engine gave SUPPORTED
where its own TOST showed the corrected effect equivalent to zero within the SESOI. Under the rule
such a SUPPORTED, a detected effect whose 90% interval lies inside +-SESOI x attenuation (the
effective SESOI; `fields.effect.detail.tost_equivalent`), becomes EFFECT BELOW SESOI. Counted on
v1's cards:

1. how many SUPPORTED the rule changes: on null data (each one a false SUPPORTED) and on the real
   effects, where the true change |Δ*| is at least the SESOI (a loss) or smaller (the new outcome
   describes it), by condition and level, with S1's key cells;
2. the same together with change 6's rule (SUPPORTED stands only where the raw and the corrected
   test detect the effect in one direction; validation/exploratory/v1_correction_rule.py).

    python validation/exploratory/v1_below_sesoi.py --results RESULTS --pilot pilot.json > v1_below_sesoi.log
"""
from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation" / "prereg"))
import numpy as np  # noqa: E402

import panel as P  # noqa: E402
import score  # noqa: E402

RISK = ("valid", "ambiguous")


def spread(xs) -> str:
    xs = [x for x in xs if x is not None and np.isfinite(x)]
    if not xs:
        return "none"
    return f"{min(xs):.2f} to {max(xs):.2f}, median {statistics.median(xs):.2f}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--pilot", type=Path, required=True)
    args = ap.parse_args(argv)
    manifest = json.loads((args.results / "manifest.json").read_text())
    pilot_raw = args.pilot.read_bytes()
    pilot = json.loads(pilot_raw)
    key = P.check_key(manifest["key"])
    reports, _ = score.read_results(args.results, key)
    entries = P.assign(key, pilot.get("dropped", ()), pilot.get("pool_size"))
    rows = score.card_rows(entries, reports, pilot)
    conds = P.conditions()
    out = ["# The outcome 'effect below the SESOI' on validation v1's reports (exploratory, post hoc; v1 and the "
           "engine unchanged)",
           f"# key {key}; pilot.json sha256 {hashlib.sha256(pilot_raw).hexdigest()}; {len(rows)} cards", ""]

    sup = []
    for r in rows:
        if r["outcome"] != P.SUPPORTED:
            continue
        det = json.loads((args.results / "reports" / f"{r['id']}.json").read_text())["fields"]["effect"]["detail"]
        raw_det, raw, est = bool(det.get("raw_detected")), det.get("raw_effect"), det.get("effect")
        c = conds[r["condition"]]
        ratio_true = None
        if c.data == "effect":
            ratio_true = abs(float(P.delta_of(c.name, r["variant"], r["pair"], pilot)["value"])) / P.sesoi_of(c, r["pair"], pilot)
        eff_sesoi = det.get("effective_sesoi")
        sup.append(dict(r, data=c.data, below=bool(det.get("tost_equivalent")),
                        kept6=bool(raw_det and raw is not None and est is not None and np.sign(raw) == np.sign(est)),
                        ratio_true=ratio_true,
                        ratio_est=(abs(est) / eff_sesoi if est is not None and eff_sesoi else None)))

    def line(label, rs):
        n, below = len(rs), sum(s["below"] for s in rs)
        return f"  {label}: {n} -> {n - below} SUPPORTED ({below} become EFFECT BELOW SESOI)"

    null_sup = [s for s in sup if s["data"] == "null"]
    eff_sup = [s for s in sup if s["data"] == "effect"]
    out += ["## 1. What the rule changes (a SUPPORTED whose own TOST shows the effect within the effective SESOI)",
            f"SUPPORTED cards: {len(sup)}; with tost_equivalent true: {sum(s['below'] for s in sup)}", "",
            "False SUPPORTED on null data (every one is an error):", line("all null cards", null_sup),
            line("metric at risk", [s for s in null_sup if s["truth"] in RISK]),
            line("blind or useless metric", [s for s in null_sup if s["truth"] not in RISK]),
            "  by condition, variant and level:"]
    for k in sorted({(s["condition"], s["variant"], s["level"]) for s in null_sup}):
        out.append("  " + line(" ".join(k), [s for s in null_sup if (s["condition"], s["variant"], s["level"]) == k]))
    out.append("  S1's key cells (false SUPPORTED of 790, limit 29):")
    for c in P.CONDITIONS:
        if c.key:
            out.append("  " + line(c.name, [s for s in null_sup if s["condition"] == c.name and s["variant"] == c.key_variant]))
    out.append("  |corrected effect| / effective SESOI of those the rule changes: "
               + spread([s["ratio_est"] for s in null_sup if s["below"]]))
    big = [s for s in eff_sup if s["ratio_true"] >= 1]
    small = [s for s in eff_sup if s["ratio_true"] < 1]
    out += ["", "True SUPPORTED on the real effects:", line("all real effects", eff_sup),
            line("|Δ*| >= SESOI (each change is a loss)", big),
            line("|Δ*| < SESOI (the new outcome describes the effect)", small),
            line("valid metric and |Δ*| >= SESOI", [s for s in big if s["truth"] == "valid"]),
            "  by condition and level:"]
    for k in sorted({(s["condition"], s["level"]) for s in eff_sup}):
        out.append("  " + line(" ".join(k), [s for s in eff_sup if (s["condition"], s["level"]) == k]))
    out += ["  |Δ*| / SESOI of those the rule changes: " + spread([s["ratio_true"] for s in eff_sup if s["below"]]),
            "  |corrected effect| / effective SESOI of those the rule changes: "
            + spread([s["ratio_est"] for s in eff_sup if s["below"]]), ""]

    out += ["## 2. Together with change 6's rule (SUPPORTED stands only where the raw and the corrected test "
            "detect the effect in one direction)"]
    for label, rs in (("null data", null_sup), ("null data, metric at risk", [s for s in null_sup if s["truth"] in RISK]),
                      ("real effects", eff_sup), ("real effects, |Δ*| >= SESOI", big)):
        n6 = sum(not s["kept6"] for s in rs)
        nb = sum(s["below"] for s in rs)
        both = sum(s["below"] and not s["kept6"] for s in rs)
        left = sum(s["kept6"] and not s["below"] for s in rs)
        out.append(f"  {label}: SUPPORTED {len(rs)}; change 6 removes {n6}, the new outcome {nb}, both {both}; "
                   f"SUPPORTED left under both rules {left}")
    out.append("  the outcomes of the changed cards under both rules: " + ", ".join(
        f"{k} {v}" for k, v in sorted(Counter(
            ("INCONCLUSIVE" if not s["kept6"] else "EFFECT BELOW SESOI") + f" ({s['data']})"
            for s in sup if s["below"] or not s["kept6"]).items()))
        + " (where both apply, change 6 is applied first: the raw test does not support the effect)")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
