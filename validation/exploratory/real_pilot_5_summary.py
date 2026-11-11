#!/usr/bin/env python3
"""The tallies of the third diagnostics' run (real_pilot_5.json, results/real-pilot-5): per acceptance claim
its verdict and cause, the effect within the domain with its interval, s_m, the domain's types, the mixture's
and the outside's estimates, GATE 4's dose, and the line of specificity; the scores; and (b)'s and (a)'s
numbers. Reads the JSON only; computes nothing new.

After the run (the same JSON, nothing recomputed) it also prints what the run's JSON holds and its log did
not: (a1)'s detections per type and group; (b4)'s columns that separate the ages or are nested within the
units, its numeric columns, and GEO's records; GATE 0's responses where GATE 0 failed; and the values per
mouse of the Marsh claims with their arm and batch (the engine's reports, in PARTS_DIR).

    python validation/exploratory/real_pilot_5_summary.py RUN_DIR/real_pilot_5.json [PARTS_DIR] > real_pilot_5_summary.log
"""
from __future__ import annotations

import json
import os
import sys

from real_pilot_5_spec import MARSH_SAMPLES


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


def a1_lines(a: dict) -> list:
    """(a1)'s detections per type and group, the types with enough cells and units."""
    a1 = a.get("a1") or {}
    rows = [r for r in a1.get("rows", []) if r.get("qualifies")]
    if not rows:
        return []
    genes = a1["rule"]["genes"]
    groups = list(rows[0]["detection"])
    cols = [f"{g} {grp}" for grp in groups for g in genes]
    lines = [f"## (a1) the detections per type ({a1['rule']['text']})", "",
             "| type | fibroblast-like | " + " | ".join(cols) + " | by the rule |",
             "|---|---|" + "---|" * len(cols) + "---|"]
    for r in rows:
        d = r["detection"]
        vals = [f"{d[grp][g]:.3f}" for grp in groups for g in genes]
        lines.append(f"| {r['type']} | {r['fibroblast_like']} | " + " | ".join(vals) + f" | {r['in_domain']} |")
    return lines + [""]


def b4_lines(b: dict) -> list:
    """(b4)'s columns that separate the ages or are nested within the units, its numeric columns, GEO."""
    rows = b.get("b4_columns") or []
    if not rows:
        return []
    lines = ["## (b4) the metadata", "",
             "Categorical columns that separate the ages or are nested within the units (up to three values each side):"]
    for r in rows:
        if r.get("outcome") in ("separates the ages (shared within an age)", "nested within units"):
            lines.append(f"- {r['source']}, {r['column']} ({r['sex']}): {r['outcome']}; {r['values']} values; "
                         f"young {r['young_values'][:3]}; old {r['old_values'][:3]}")
    lines.append("\nNumeric columns, the units' medians, young minus old, exact p over the relabelings:")
    for r in rows:
        if r.get("kind") == "numeric":
            lines.append(f"- {r['source']}, {r['column']} ({r['sex']}): {f(r.get('young_minus_old'), 5)}, "
                         f"p {f(r.get('p'), 3)} ({r.get('relabelings')} relabelings)")
    lines.append("\nGEO:")
    for acc, g in (b.get("b4_geo") or {}).items():
        if isinstance(g, dict):
            first = ["; ".join(t.split(" | ")[:2]) for t in g.get("titles", [])]
            lines.append(f"- {acc}: {g.get('samples')} samples; case 1's units named: "
                         f"{sum(g.get('units_named', {}).values())}; the series' lines quoted: "
                         f"{g.get('series_lines')}; the first samples: {first}")
    return lines + [""]


def report_of(parts_dir: str, part: str, cid: str):
    p = os.path.join(parts_dir, f"part-{part}", "reports", f"{cid}.json")
    return json.loads(open(p).read()) if os.path.exists(p) else None


def gate0_lines(rep: dict, cid: str) -> list:
    """GATE 0's responses where it failed."""
    g0 = next((g for g in rep.get("gates", []) if g.get("gate") == 0), None)
    if g0 is None or g0.get("status") != "FAIL":
        return []
    d = g0["detail"]
    lines = [f"### {cid}: GATE 0's responses: the claim's value {f(d.get('baseline'))}, the null's centre "
             f"{f(d.get('null_center'))}, the signal above the null {f(d.get('signal'))}; the tolerance "
             f"{f(d.get('bias_tolerance_abs'))} ({d.get('bias_tolerance')} x s_m {f(d.get('sesoi'))})"]
    for k, r in (d.get("responses") or {}).items():
        null = (f"; the null moves {f(r['null_shift'])} (z {f(r.get('null_shift_z'), 3)})"
                if r.get("null_shift") is not None else "")
        lines.append(f"- {k}: the claim's value moves {f(r.get('shift'))} (z {f(r.get('z'), 3)}){null}; "
                     f"{r.get('classification')}, {r.get('handling', 'not blocking')}")
    return lines + [""]


def unit_lines(rep: dict, cid: str) -> list:
    """The values per replicate within the domain (raw / at equal depth), with the Marsh arm and batch."""
    units = ((rep.get("fields", {}).get("effect", {}).get("detail") or {}).get("domain") or {}).get("units") or []
    if not units:
        return []
    lines = [f"### {cid}: the values per mouse within the domain (raw / at equal depth)"]
    for u in units:
        _, arm, batch = MARSH_SAMPLES.get(u["replicate"], ("", "?", "?"))
        lines.append(f"- {u['replicate']} {arm}, {batch}: {f(u.get('raw'))} / {f(u.get('corrected'))}")
    return lines + [""]


def main(argv=None) -> int:
    args = argv or sys.argv[1:]
    path = args[0]
    parts_dir = args[1] if len(args) > 1 else os.path.join(os.path.dirname(os.path.abspath(path)), "parts")
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
    print("\n# What the JSON holds beside the log (the same run, nothing recomputed)\n")
    print("\n".join(a1_lines(a)))
    print("\n".join(b4_lines(b)))
    for part, cid in (("case1m", "A3"), ("case1", "A4")):
        rep = report_of(parts_dir, part, cid)
        if rep is not None:
            print("\n".join(gate0_lines(rep, cid)))
    for cid in ("A5", "A6", "A7", "A5m", "A6m", "A7m"):
        rep = report_of(parts_dir, "marsh", cid)
        if rep is not None:
            print("\n".join(unit_lines(rep, cid)))
    print("## The claims' own donor-level t, as the line of specificity computes it (within the domain, "
          "uncorrected; the engine's reports): in the claimed direction (positive where the claim's "
          "direction holds; a two-sided claim takes its observed one), and as groups[0] - groups[1]\n")
    for part in ("case2", "case1m", "case1", "marsh", "case3a", "case3b", "case3c"):
        for r in ((res["results"].get(part) or {}).get("claims") or []):
            rep = report_of(parts_dir, part, r["id"])
            s = (rep or {}).get("specificity") or {}
            if s.get("t_claim") is not None:
                print(f"- {r['id']}: {f(s.get('t_claim_signed'), 4)}; groups[0] - groups[1] {f(s['t_claim'], 4)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
