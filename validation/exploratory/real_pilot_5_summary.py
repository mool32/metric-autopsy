#!/usr/bin/env python3
"""The tallies of the third diagnostics' run (real_pilot_5.json, results/real-pilot-5): per acceptance claim
its verdict and cause, the effect within the domain with its interval, s_m, the domain's types, the mixture's
and the outside's estimates, GATE 4's dose, and the line of specificity; the scores; and (b)'s and (a)'s
numbers. Reads the JSON only; computes nothing new.

    python validation/exploratory/real_pilot_5_summary.py RUN_DIR/real_pilot_5.json > real_pilot_5_summary.log
"""
from __future__ import annotations

import json
import sys


def f(x, n=4):
    if x is None:
        return "-"
    if isinstance(x, (list, tuple)):
        return "[" + ", ".join(f(v, n) for v in x) + "]"
    try:
        return f"{float(x):+.{n}g}"
    except (TypeError, ValueError):
        return str(x)


def claim_lines(r: dict) -> list:
    if "error" in r:
        return [f"### {r['id']}: EXCEPTION {r['error']}"]
    c = r["claim"]
    ef, g4, dom, spec = r.get("effect") or {}, r.get("gate4") or {}, r.get("domain") or {}, r.get("specificity") or {}
    mix = (r.get("mixture") or {}).get("corrected") or {}
    out = (dom.get("outside") or {}).get("corrected") or {}
    lines = [f"### {r['id']}: {c['metric']}({', '.join(c['genes'])}), {c['groups'][0]} -> {c['groups'][1]}, "
             f"{c['direction']}, SESOI {c['sesoi']} (construct)",
             f"- verdict: {r['verdict']} [cause {r['cause']}]",
             "- gates: " + "; ".join(f"{g['gate']} {g['status']}" for g in r["gates"]),
             f"- effect ({ef.get('estimate_over')}): {f(ef.get('effect'))} CI {f(ef.get('ci'))}, raw {f(ef.get('raw_effect'))}; "
             f"test {ef.get('test')}; SESOI used {f(ef.get('effective_sesoi'))} ({ef.get('sesoi_source')}); floor {f(ef.get('floor'))}",
             f"- s_m {f(g4.get('sesoi_metric'))} CI {f(g4.get('sesoi_metric_ci'))}; GATE 4 {g4.get('status')}; dose "
             f"{json.dumps(g4.get('dose'), default=str)[:200]}",
             f"- domain {dom.get('types')}; mixture {f(mix.get('effect'))} (p {f(mix.get('p'), 3)}, detected {mix.get('detected')}); "
             f"outside {(dom.get('outside') or {}).get('types')}: {f(out.get('effect'))} (p {f(out.get('p'), 3)}, "
             f"detected {out.get('detected')})"]
    if spec and not spec.get("error") and spec.get("applicable") is not False:
        d, a, sh = spec.get("in_direction") or {}, spec.get("absolute") or {}, spec.get("shift") or {}
        lines.append(f"- specificity: {d.get('k')} of {d.get('n')} controls at least as far in the direction "
                     f"(share CI {f(d.get('share_ci'), 3)}, {d.get('outcome')}); by |t| {a.get('k')} ({a.get('outcome')}); "
                     f"the controls' median t {f(sh.get('median_t'), 3)} CI {f(sh.get('ci'), 3)}"
                     + (f" — WARNING: {sh['warning']}" if sh.get("warning") else ""))
    elif spec:
        lines.append(f"- specificity: {spec.get('error') or spec.get('reason')}")
    return lines


def main(argv=None) -> int:
    path = (argv or sys.argv[1:])[0]
    res = json.loads(open(path).read())
    sc = res["scores"]
    print("# The tallies of the third diagnostics and v0.4's acceptance tests (real_pilot_5_summary.py)\n")
    c = sc["claims"]
    print("| test | criterion | verdict | met |\n|---|---|---|---|")
    for k in ("A1", "A2", "A3", "A4", "A5", "A7"):
        print(f"| {k} | {c[k]['criterion']} | {c[k].get('label')} | {c[k].get('met')} |")
    print(f"| A8 | strict / lenient | {c['A8']['n']} claims; SUPPORTED {c['A8']['supported'] or 'none'} | "
          f"{c['A8']['strict']} / {c['A8']['lenient']} |\n")
    print("A8 verdicts: " + "; ".join(f"{k} {v}" for k, v in sorted(c["A8"]["labels"].items())) + "\n")
    print(f"(b) the prediction: {sc['b_prediction'].get('score')}; (b4) the statement: {sc['b_statement']}")
    print(f"(a1): {sc['a1']}\n(a2): {sc['a2']}\n")
    for part in ("case2", "case1m", "case1", "marsh", "case3a", "case3b", "case3c"):
        for r in ((res["results"].get(part) or {}).get("claims") or []):
            print("\n".join(claim_lines(r)) + "\n")
    b = res["results"].get("b") or {}
    if "b2" in b:
        e = b["b2"]["F"]["non_fibroblast"]
        print(f"(b2) F young minus old {f(e['delta'])}, exact p {f(e['p'], 3)} ({e['relabelings']} relabelings); "
              f"(b3-i) r {f(b['b3_i']['r'], 3)} p {f(b['b3_i']['p'], 3)}; (b3-0) r {f(b['b3_0']['r'], 3)} p "
              f"{f(b['b3_0']['p'], 3)}; (b3-ii) r {f(b['b3_ii']['r'], 3)} p {f(b['b3_ii']['p'], 3)}")
        for g, v in b["b1"].items():
            print(f"  {g}: detection {f(v['det']['non_fibroblast']['delta'])} (p {f(v['det']['non_fibroblast']['p'], 3)}), "
                  f"level {f(v['lvl']['non_fibroblast']['delta'])} (p {f(v['lvl']['non_fibroblast']['p'], 3)})")
    sep = [(r["source"], r["column"], r["sex"]) for r in b.get("b4_columns", [])
           if r.get("outcome") == "separates the ages (shared within an age)"]
    print(f"(b4) columns that separate the ages: {sep or 'none'}")
    a = res["results"].get("a") or {}
    if a.get("a1"):
        print(f"(a1) domain {a['a1']['domain']}; outside {len(a['a1']['outside'])} types")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
