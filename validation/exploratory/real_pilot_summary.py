#!/usr/bin/env python3
"""The tallies of the pilot on real data (validation/exploratory/real_pilot_plan.md), from its merged
output (run/real_pilot.json on the branch results/real-pilot, or results/real-pilot-2 for the
second run). Exploratory; it reads nothing else and computes no new statistic: the counts of the
conclusions of (a), (b) and (d), where (d) changes (b) and by which gate, the effects found against
their SESOI (the SESOI used: the design's floor where the claim's SESOI is "floor"), each claim's
per-donor values and, in the second run, its floor and per-stratum tests, and the time of the
analyses.

    git show origin/results/real-pilot:run/real_pilot.json > real_pilot.json
    python validation/exploratory/real_pilot_summary.py real_pilot.json > real_pilot_summary.log
"""
from __future__ import annotations

import hashlib
import json
import statistics
import sys
from collections import Counter
from pathlib import Path


def tally(claims, key):
    return dict(sorted(Counter(c.get(key[0], {}).get(key[1]) for c in claims).items(), key=lambda x: str(x[0])))


def main(path: str) -> int:
    raw = Path(path).read_bytes()
    d = json.loads(raw)
    print(f"# The pilot's tallies (from {Path(path).name}, sha256 {hashlib.sha256(raw).hexdigest()})")
    print(f"engine {d['engine']['version']} ({d['engine']['path']}); python {d['python']}; numpy {d['numpy']}")
    print("parts: " + "; ".join(f"{k} {v.get('seconds')} s, complete {v.get('complete')}" for k, v in d["parts"].items()))
    for name, case in d["cases"].items():
        print(f"\n## {name}: {case.get('status')}" + (f" ({case.get('reason') or case.get('error')})"
                                                     if case.get("status") != "run" else ""))
        claims = [c for c in case.get("claims", []) if "skipped" not in c]
        if not claims:
            continue
        print(f"claims: {len(claims)}; cells {case.get('cells')}, genes {case.get('genes')}")
        for an in ("a", "b", "d"):
            print(f"({an}) conclusions: {tally(claims, (an, 'conclusion'))}")
        print(f"(d) under the known issue: {tally(claims, ('d', 'known_issue_conclusion'))}")
        print(f"(d) verdicts: {dict(Counter(c['d'].get('label') for c in claims))}; causes: "
              f"{dict(Counter(c['d'].get('cause') for c in claims))}")
        print(f"raw_detected among (d)'s SUPPORTED: "
              f"{[c['d'].get('raw_detected') for c in claims if c['d'].get('label') == 'SUPPORTED']}")
        ab = sum(c["a"].get("conclusion") != c["b"].get("conclusion") for c in claims)
        print(f"(a) differs from (b) on {ab} of {len(claims)}")
        changed = [c for c in claims if c.get("d_vs_b") != "same"]
        print(f"(d) changes (b) on {len(changed)} of {len(claims)}:")
        for c in changed:
            print(f"  {c['id']}: (b) {c['b'].get('conclusion')} -> (d) {c['d'].get('conclusion')} "
                  f"[{c['d'].get('cause')}]; {c['d_vs_b']}")
        known = [c for c in claims if c.get("d_known_vs_b") != c.get("d_vs_b")]
        print(f"the known issue changes the comparison on {len(known)} of {len(claims)}")
        effects = [c for c in claims if "effect" in (c["b"].get("conclusion"), c["d"].get("conclusion"))]
        if effects:
            print("effects found by (b) or (d), against the SESOI (Delta second - first; (d) corrected, first - second):")
            for c in effects:
                b, e = c["b"], c["d"].get("effect", {})
                use = b.get("sesoi_used", c["sesoi"])
                print(f"  {c['id']} {c['metric']}({', '.join(c['genes'])}): (b) {b.get('conclusion')}, Delta {b.get('delta'):+.4f}, "
                      f"90% CI [{b['ci90'][0]:+.4f}, {b['ci90'][1]:+.4f}], p {b.get('p'):.4g}; (d) {c['d'].get('label')}, "
                      f"effect {e.get('effect', float('nan')):+.4f}, CI {[round(x, 4) for x in e.get('ci') or []]}, "
                      f"p {e.get('p_perm')}, TOST equivalent {e.get('tost_equivalent')} within +-"
                      f"{'none' if e.get('effective_sesoi') is None else format(e['effective_sesoi'], '.4g')}; "
                      f"SESOI {c['sesoi']}" + (f" ({use:.4g})" if use != c["sesoi"] else "")
                      + f"; |Delta| / SESOI {abs(b.get('delta')) / use:.2f}")
        print("(b)'s per-donor values (the metric per donor, or per donor the first group minus the second where "
              "every donor is in both), min to max:")
        for c in claims:
            v = list((c["b"].get("values") or {}).values())
            if v:
                use = c["b"].get("sesoi_used", c["sesoi"])
                print(f"  {c['id']} {c['metric']}: {min(v):+.4f} to {max(v):+.4f} over {len(v)} donors; SESOI {c['sesoi']}"
                      + (f" ({use:.4g})" if use != c["sesoi"] else ""))
        floors = [c for c in claims if "floor" in c["b"]]
        if floors:
            print("(b)'s design floor (v0.4's change 4, by hand) and the SESOI used:")
            for c in floors:
                b = c["b"]
                print(f"  {c['id']} {c['metric']}: floor {b['floor']:.4g} (SE {b['se']:.4g}, df {b['df']:.4g}); "
                      f"SESOI {c['sesoi']}, used {b['sesoi_used']:.4g}; (b) {b['conclusion']}, (d) {c['d'].get('label')}"
                      + (f" with SESOI {c['d']['prereg']['sesoi']:.4g} ({c['d']['prereg'].get('sesoi_scale', 'construct')} scale)"
                         if c["d"].get("prereg") else ""))
        strata = [c for c in claims if c["b"].get("by_stratum")]
        if strata:
            print("(b) per stratum (the exact permutation within each):")
            for c in strata:
                for st, r in c["b"]["by_stratum"].items():
                    g0, g1 = list(r["n"])
                    print(f"  {c['id']} {c['metric']} {st}: {r['n'][g0]} against {r['n'][g1]}, means {r['mean'][g0]:.4g} "
                          f"against {r['mean'][g1]:.4g}, p {r['p']:.4g} ({r['method']}, {r['n_relabelings']})")
        audit = [c for c in claims if c.get("b_vs_audit")]
        for c in audit:
            v = c["b_vs_audit"]
            print(f"{c['id']}: (b)'s per-mouse values against the audit's section E: "
                  + (v.get("note") or f"{v['compared']} mice, max |diff| {v['max_abs_diff']:.2e}, only here "
                     f"{v['only_here'] or 'none'}, only in the audit {v['only_audit'] or 'none'}"))
        sec = [c["d"].get("seconds") for c in claims if c["d"].get("seconds") is not None]
        by = {}
        for c in claims:
            by.setdefault(c["metric"], []).append(c["d"].get("seconds"))
        print("(d) seconds per claim: " + "; ".join(
            f"{m} n={len(v)} median {statistics.median(v):.0f} (min {min(v):.0f}, max {max(v):.0f})" for m, v in by.items())
            + f"; all {sum(sec):.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
