#!/usr/bin/env python3
"""The tallies of the diagnostics on the second run's data (validation/exploratory/
real_pilot_diagnostics_plan.md), from their merged output (run/real_pilot_3.json on the branch
results/real-pilot-3). Exploratory, after the run; it changes no outcome and no rule. It prints, by
check and claim, the outcomes and the details behind them that the merged log does not list: the
controls of (ii) at least as far as the claim, in its direction and in absolute size, with the
Clopper-Pearson 95% interval of their share (how much the draw of 100 controls can move the
outcome), the detection of the drawn controls, (iii)'s values per mouse, and (i)'s thinning.

    git show origin/results/real-pilot-3:run/real_pilot_3.json > real_pilot_3.json
    python validation/exploratory/real_pilot_3_summary.py real_pilot_3.json > real_pilot_3_summary.log
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

from metric_autopsy.stats import clopper_pearson

ALPHA = 0.05


def fmt(x, nd=4) -> str:
    return "nan" if x is None else f"{x:+.{nd}g}"


def main(path: str) -> int:
    raw = Path(path).read_bytes()
    d = json.loads(raw)
    print(f"# The diagnostics' tallies (from {Path(path).name}, sha256 {hashlib.sha256(raw).hexdigest()})")
    print(f"engine {d['engine']['version']}; python {d['python']}; numpy {d['numpy']}; STRING "
          + "; ".join(f"{k} {v}" for k, v in d["string_sha256"].items()))
    print("parts: " + "; ".join(f"{k} {v.get('seconds')} s, complete {v.get('complete')}" for k, v in d["parts"].items()))
    ids = list(d["outcomes"]["i"])
    print("\n## Outcomes (the merged log's table)\n")
    print("| check | " + " | ".join(ids) + " | separates |")
    print("|---|" + "---|" * (len(ids) + 1))
    for chk, row in d["outcomes"].items():
        print(f"| ({chk}) | " + " | ".join(row[i] for i in ids) + f" | {d['separates'][chk]} |")
    print(f"\nthe prediction holds: {d['prediction']['holds']}")

    print("\n## (ii) the controls at least as far as the claim")
    sel = {"1a": d["controls"]["case1"]["pairs"], "1a-males": d["controls"]["case1"]["pairs"],
           "1b": d["controls"]["case1"]["pairs"], "2a": d["controls"]["case2"]["2a"], "2b": d["controls"]["case2"]["2b"]}
    for cid in ids:
        c, s = d["claims"][cid], sel[cid]
        r = c["ii"]
        deltas, ctl = r["deltas"], s["controls"]
        far = [(x, g) for x, g in zip(deltas, ctl)
               if (x <= c["delta"] if c["direction"] == "decrease" else x >= c["delta"])]
        both = [(x, g) for x, g in zip(deltas, ctl) if abs(x) >= abs(c["delta"])]
        lo, hi = clopper_pearson(len(far), r["n"])
        print(f"\n### {cid}: {c['metric']}({', '.join(c['genes'])}), {c['direction']}, Delta {fmt(c['delta'])}")
        print(f"- in the claimed direction: {len(far)} of {r['n']} controls at least as far (p {r['p_one_sided']:.4g}; "
              f"the share {len(far) / r['n']:.3f}, 95% CI {lo:.3f} to {hi:.3f}: the outcome {r['outcome']!r} "
              + ("holds over the whole interval" if (hi < ALPHA) or (lo >= ALPHA) else
                 "is within the interval's reach of the 0.05 line") + ")")
        for x, g in sorted(far, key=lambda t: t[0], reverse=c["direction"] == "increase"):
            print(f"    {fmt(x)}  {' + '.join(g)}")
        lo2, hi2 = clopper_pearson(len(both), r["n"])
        print(f"- in absolute size: {len(both)} of {r['n']} at least as large (two-sided p {r['p_two_sided']:.4g}; "
              f"95% CI of the share {lo2:.3f} to {hi2:.3f})"
              + (": " + ", ".join(f"{' + '.join(g)} {fmt(x)}" for x, g in sorted(both, key=lambda t: -abs(t[0])))
                 if both and len(both) <= 15 else ""))
        det = s["controls_detection"]
        print("- the drawn controls' detection, by the claim gene they match: " + "; ".join(
            f"{g} {s['claim_detection'][g]:.4f}: {min(x[k] for x in det):.4f} to {max(x[k] for x in det):.4f}"
            for k, g in enumerate(c["genes"])))

    print("\n## (iii) per mouse: the claim's value, the response at the largest dose, the change at full removal")
    for cid in ids:
        c = d["claims"][cid]
        r = c["iii"]
        print(f"\n### {cid}: {r['outcome']}" + (f" ({', '.join(r['failed'])})" if r["failed"] else ""))
        for st, s in r["strata"].items():
            print(f"- {st}: h {s['higher']}, l {s['lower']}; |Delta| {s['abs_delta']:.4g}; (a) removal in h "
                  f"{fmt(s['removal_higher'])}, |Delta| / its size {s['ratio_a']:.3g}; (b) the largest mean response in "
                  f"l {fmt(s['max_response_lower'])} at {s['max_dose']}, |Delta| / it {s['ratio_b']:.3g}")
        top = max(next(iter(r["per_mouse"].values()))["dose"], key=float)
        for m, v in sorted(r["per_mouse"].items()):
            print(f"    {m}: value {v['value']:.5f}; response at {top} {fmt(v['dose'][top], 5)}; removal "
                  f"{fmt(v['removal'], 5)}")

    print("\n## (i) the engine without GATE 4 and GATE 5")
    for cid in ids:
        e = d["claims"][cid]["i"]
        f = e["effect"]
        print(f"- {cid}: {e['outcome']} ({e['why']}); raw {fmt(f['raw_effect'])} (p {f['p_raw']:.4g}), corrected "
              f"{fmt(f['effect'])} (p {f['p_perm']:.4g}, {f['status']}), retained {f['retained']:.4f}; thinning "
              + "; ".join(f"{q['stratum']}: {q['thinned_group']} x {q['depth_ratio']:.4f}" for q in e["equalization"])
              + f"; verdict {e['label']} [{e['cause']}]; {e['seconds']} s")

    print("\n## (iv) the two measures by stratum")
    for cid in ids:
        r = d["claims"][cid]["iv"]
        print(f"- {cid} ({r['measures'][0]} against {r['measures'][1]}): {r['outcome']}; " + "; ".join(
            f"{st} {fmt(s['own']['delta'])} (p {s['own']['p']:.4g}) against {fmt(s['other']['delta'])} "
            f"(p {s['other']['p']:.4g})" for st, s in r["strata"].items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
