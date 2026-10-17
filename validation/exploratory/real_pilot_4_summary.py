#!/usr/bin/env python3
"""The tallies of the second diagnostics on the second run's data (validation/exploratory/
real_pilot_diagnostics_2_plan.md), from their output (run/real_pilot_4.json on the branch
results/real-pilot-4). Exploratory, after the run; it changes no outcome and no rule. It counts, per
sex, the cell types whose difference of MI (and of norm_pearson) falls, rises or stays, the share of
the overall difference that the types keep, the fibroblast-like types against the others, and, per
claim of (ii′), the positions and the quantiles of the controls' t.

    git show origin/results/real-pilot-4:run/real_pilot_4.json > real_pilot_4.json
    python validation/exploratory/real_pilot_4_summary.py real_pilot_4.json > real_pilot_4_summary.log
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

FIBRO = ("fibroblast", "mesenchymal", "stromal")  # the plan's rule for fibroblast-like labels


def fmt(x, nd=4) -> str:
    return "nan" if x is None else f"{x:+.{nd}g}"


def signs(values) -> str:
    neg = sum(v < 0 for v in values)
    pos = sum(v > 0 for v in values)
    return f"{neg} lower, {pos} higher, {len(values) - neg - pos} equal"


def mean(values):
    return sum(values) / len(values) if values else None


def main(path: str) -> int:
    raw = Path(path).read_bytes()
    d = json.loads(raw)
    print(f"# The second diagnostics' tallies (from {Path(path).name}, sha256 {hashlib.sha256(raw).hexdigest()})")
    print(f"engine {d['engine']['version']}; python {d['python']}; numpy {d['numpy']}; settings {d['settings']}")

    print("\n## (v) the composition of cells on §4")
    comp = d["composition"]
    for sx, r in comp["sexes"].items():
        print(f"\n### {sx}: {len(r['units'])} units; MI over all cells, old minus young {fmt(r['delta_all'])}")
        for what in ("cell_type", "tissue"):
            c = r[what]
            print(f"- {what}: total variation distance {c['tvd']:.4f}, exact p {c['tvd_p']:.4g} ({c['relabelings']} "
                  f"relabelings); levels {len(c['by_level'])}, with an exact p below 0.05 (uncorrected): "
                  f"{sum(x['p'] < 0.05 for x in c['by_level'])}")
        for kind in ("mi_3bin", "norm_pearson"):
            w = r["within"][kind]
            types = w["types"]
            deltas = {t: v["delta"] for t, v in types.items()}
            fib = [v for t, v in deltas.items() if any(f in t.lower() for f in FIBRO)]
            other = [v for t, v in deltas.items() if not any(f in t.lower() for f in FIBRO)]
            comb = w["combined"]
            # the share kept is the plan's measure for the males' drop only (the females' overall
            # difference is about zero, so a share of it means nothing)
            kept = (comb["delta"] / r["delta_all"]) if kind == "mi_3bin" and sx == "male" else None
            print(f"- {kind} within {len(types)} types: combined {fmt(comb['delta'])}, exact p {comb['p']:.4g}"
                  + (f", {kept:.1%} of the overall difference" if kept is not None else "")
                  + f"; types {signs(list(deltas.values()))}; with an exact p below 0.05 (uncorrected): "
                    f"{sum(v['p'] < 0.05 for v in types.values())}")
            print(f"    fibroblast-like types ({len(fib)}): mean difference {fmt(mean(fib))}, {signs(fib)}; "
                  f"the other types ({len(other)}): {fmt(mean(other))}, {signs(other)}")
            for t in sorted(t for t in deltas if any(f in t.lower() for f in FIBRO)):
                v = types[t]
                print(f"    {t}: {fmt(v['delta'])} (p {v['p']:.4g}; young {v['n']['young']}, old {v['n']['old']} units)")
        a = r["association"]
        print(f"- MI against the share of fibroblast-like cells: Pearson r {fmt(a['pearson']['r'], 3)} (p "
              f"{a['pearson']['p']:.4g}), Spearman rho {fmt(a['spearman']['r'], 3)} (p {a['spearman']['p']:.4g}), "
              f"{a['orderings']} orderings")
    o = comp.get("outcome", {})
    print(f"\n(v) in males: not detected {o.get('not_detected')}; less than half {o.get('less_than_half')}; composition "
          f"accounts for the drop: {o.get('accounts_for_the_drop')}; the prediction holds: {comp['prediction']['holds']}")

    print("\n## (ii′) specificity on the donor-level t")
    sp = d["specificity"]
    for cid in ("1a-males", "2a", "2b"):
        r = sp[cid]
        print(f"\n### {cid}: {r['metric']}({', '.join(r['genes'])}), {r['first']} -> {r['second']}, {r['direction']}: "
              f"Delta {fmt(r['delta'])}, SE {r['se']:.4g}, t {fmt(r['t'])}; {r['ranked']} of {r['controls']} controls "
              "with a t")
        for key in ("claimed", "absolute"):
            c = r[key]
            print(f"- {key}: k {c['k']} (share 95% CI {c['share_ci'][0]:.4f} to {c['share_ci'][1]:.4f}), p {c['p']:.4g} "
                  f"-> {c['outcome']}")
        print("- the controls' t: " + ", ".join(f"{q} {fmt(v)}" for q, v in r["control_t_quantiles"].items()))
        if "vs_run2_max_abs_diff" in r:
            print(f"- per mouse against the second run's values: max |diff| {r['vs_run2_max_abs_diff']:.3g}")
    pr = sp["prediction"]
    print("\n(ii′) " + "; ".join(f"{k}: {v}" for k, v in pr["parts"].items()) + f" -> the prediction holds: {pr['holds']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
