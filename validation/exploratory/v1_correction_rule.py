"""The rule "a correction can only remove a conclusion, never create one", applied after the fact
to validation v1's published reports (exploratory, post hoc; the owner's decision of 2026-10-11 on
v0.4 change 6). v1's result and the engine are unchanged.

Under the rule SUPPORTED stands only where the engine's raw replicate-level test (the same data
without the correction: `fields.effect.detail.raw_detected`, `raw_effect`) and the corrected test
both detect the effect, in the same direction; every other SUPPORTED becomes INCONCLUSIVE. For a
user of v0.3.0 the same reading is the known issue of the release notes: a SUPPORTED whose
`raw_detected` is false was made by the depth correction.

1. What the rule changes on v1's cards: the false SUPPORTED on null data (all, at risk, by
   condition, and S1's key cells) and the true SUPPORTED on the real effects (all, and the valid
   cards with |Δ*| >= SESOI); for those it removes, what the correction did to the raw difference
   (enlarged, shrank or reversed it), and how they relate to the false "effects" that only (d)
   makes in `v1_analyzers.log` section 5.
2. Where v1's four analyzers differ on the null cards (`v1_analyzers.log`): the false "effect" of
   (b) and (d) on blind metrics and on the cards at risk, so that the share of the reduction that is
   the refusal on blind metrics is a number.
3. Absences the correction made: (d)'s NO DETECTABLE EFFECT where (c), the same engine on the
   uncorrected data, does not reach it (whether the rule should bind absence claims too).
4. An illustration, on v1's cards, of the exchange proposed for v2's primary comparison
   (validation/plans/v2_requirements.md): erroneous definite conclusions (a claim of an effect, in
   either direction, on null data or against the true direction on a real one; a claim of absence
   on a real effect with |Δ*| >= SESOI or on a blind or useless metric) and sensitivity (an effect
   in the true direction on the real effects at risk), for (b), (c), (d) and (d) under the rule.

    python validation/exploratory/v1_correction_rule.py --results RESULTS --pilot pilot.json \
        --cards validation/exploratory/v1_analyzers_cards.jsonl.gz > v1_correction_rule.log
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation" / "prereg"))
sys.path.insert(0, str(ROOT / "validation" / "exploratory"))
import numpy as np  # noqa: E402

import panel as P  # noqa: E402
import score  # noqa: E402
from v1_by_level import fmt  # noqa: E402

RISK = ("valid", "ambiguous")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--pilot", type=Path, required=True)
    ap.add_argument("--cards", type=Path, required=True)
    args = ap.parse_args(argv)
    manifest = json.loads((args.results / "manifest.json").read_text())
    pilot_raw = args.pilot.read_bytes()
    pilot = json.loads(pilot_raw)
    key = P.check_key(manifest["key"])
    reports, _ = score.read_results(args.results, key)
    entries = P.assign(key, pilot.get("dropped", ()), pilot.get("pool_size"))
    rows = score.card_rows(entries, reports, pilot)
    conds = P.conditions()
    cards_raw = gzip.decompress(args.cards.read_bytes())
    mine = {json.loads(line)["card"]: json.loads(line) for line in cards_raw.decode().splitlines() if line.strip()}
    out = ["# The rule 'a correction can only remove a conclusion, never create one' on validation v1's reports",
           "# (exploratory, post hoc; v1 and the engine unchanged)",
           f"# key {key}; pilot.json sha256 {hashlib.sha256(pilot_raw).hexdigest()}; "
           f"per-card analyzers sha256 {hashlib.sha256(cards_raw).hexdigest()}; {len(rows)} cards", ""]

    # every SUPPORTED card with its raw replicate-level test
    sup = []
    for r in rows:
        if r["outcome"] != P.SUPPORTED:
            continue
        det = json.loads((args.results / "reports" / f"{r['id']}.json").read_text())["fields"]["effect"]["detail"]
        raw_det, raw, est = bool(det.get("raw_detected")), det.get("raw_effect"), det.get("effect")
        same = raw_det and raw is not None and est is not None and np.sign(raw) == np.sign(est)
        c = conds[r["condition"]]
        big = None
        if c.data == "effect":
            big = abs(float(P.delta_of(c.name, r["variant"], r["pair"], pilot)["value"])) >= P.sesoi_of(c, r["pair"], pilot)
        sup.append(dict(r, raw_detected=raw_det, kept=bool(same), data=c.data, big=big, raw=raw, est=est,
                        p_raw=det.get("p_raw")))
    raw_opposite = sum(s["raw_detected"] and not s["kept"] for s in sup)

    out += ["## 1. What the rule changes (SUPPORTED stands only where the raw and the corrected test detect it, "
            "in one direction)",
            f"SUPPORTED cards: {len(sup)}; with raw_detected false: {sum(not s['raw_detected'] for s in sup)}; "
            f"raw detected but in the other direction: {raw_opposite}", ""]
    null_sup = [s for s in sup if s["data"] == "null"]
    eff_sup = [s for s in sup if s["data"] == "effect"]

    def gone(rs):
        return sum(not s["kept"] for s in rs)

    out += ["False SUPPORTED on null data (every one is an error):",
            f"  all null cards: {len(null_sup)} -> {len(null_sup) - gone(null_sup)} "
            f"({gone(null_sup)} become INCONCLUSIVE)",
            f"  metric at risk: {sum(s['truth'] in RISK for s in null_sup)} -> "
            f"{sum(s['truth'] in RISK and s['kept'] for s in null_sup)}",
            f"  blind or useless metric: {sum(s['truth'] not in RISK for s in null_sup)} -> "
            f"{sum(s['truth'] not in RISK and s['kept'] for s in null_sup)}"]
    by = Counter((s["condition"], s["variant"], s["level"]) for s in null_sup)
    out.append("  by condition, variant and level (before -> after):")
    for (cname, var, lv), n in sorted(by.items()):
        rs = [s for s in null_sup if (s["condition"], s["variant"], s["level"]) == (cname, var, lv)]
        out.append(f"    {cname} {var} {lv}: {n} -> {n - gone(rs)}")
    out.append("  S1's key cells (false SUPPORTED of 790, limit 29), before -> after:")
    for c in P.CONDITIONS:
        if c.key:
            rs = [s for s in null_sup if s["condition"] == c.name and s["variant"] == c.key_variant]
            out.append(f"    {c.name}: {len(rs)} -> {len(rs) - gone(rs)}")

    def done(s):
        """What the correction did to the raw difference (retained = corrected / raw effect)."""
        if not s["raw"] or s["est"] is None:
            return "no raw effect"
        k = s["est"] / s["raw"]
        return "enlarged it" if k > 1 else ("shrank it" if k > 0 else "reversed its sign")

    def split(rs):
        return ", ".join(f"{k} {v}" for k, v in sorted(Counter(done(s) for s in rs).items()))

    lost = [s for s in null_sup if not s["kept"]]
    out.append(f"  what the correction did to the raw difference of the {len(lost)} (retained = corrected / raw "
               f"effect; enlarged > 1, shrank 0 to 1, reversed < 0): {split(lost)}; by condition: "
               + "; ".join(f"{cn} {split([s for s in lost if s['condition'] == cn])}"
                           for cn in sorted({s['condition'] for s in lost})))
    d_only = [s for s in null_sup if s["truth"] in RISK and mine[s["id"]]["b_call"] != "+"]
    kept_d = [s for s in d_only if s["raw_detected"]]
    out.append(f"  against v1_analyzers.log section 5 (the false 'effects' of (d) at risk that (b) does not make): "
               f"{len(d_only)}; with raw_detected false {sum(not s['raw_detected'] for s in d_only)}, the "
               f"rule's {len(lost)} {'all among them' if all(s in d_only for s in lost) else 'NOT all among them'}; "
               f"raw detected {len(kept_d)}, which the rule keeps: "
               + ", ".join(f"{k} {v}" for k, v in sorted(Counter(
                   f"{s['condition']} {s['variant']}, (b) {mine[s['id']]['b_test']}" for s in kept_d).items()))
               + (f"; (b)'s p on them {min(mine[s['id']]['b_p'] for s in kept_d):.3f} to "
                  f"{max(mine[s['id']]['b_p'] for s in kept_d):.3f}, the engine's raw p "
                  f"{min(float(s['p_raw']) for s in kept_d):.4f} to {max(float(s['p_raw']) for s in kept_d):.4f}"
                  if kept_d else ""))
    out += ["", "True SUPPORTED on the real effects (each one the rule turns into INCONCLUSIVE is a loss):",
            f"  all real effects: {len(eff_sup)} -> {len(eff_sup) - gone(eff_sup)}",
            "  " + ", ".join(f"{c.name} {sum(s['condition'] == c.name for s in eff_sup)} -> "
                             f"{sum(s['condition'] == c.name and s['kept'] for s in eff_sup)}"
                             for c in P.CONDITIONS if c.data == "effect"),
            f"  valid metric and |Δ*| >= SESOI: {sum(s['truth'] == 'valid' and s['big'] for s in eff_sup)} -> "
            f"{sum(s['truth'] == 'valid' and s['big'] and s['kept'] for s in eff_sup)}",
            "  what the correction did to the raw difference of the "
            f"{sum(not s['kept'] for s in eff_sup)}: {split([s for s in eff_sup if not s['kept']])}", ""]

    # 2. the reduction of false 'effect' between (b) and (d), split by the metric's truth
    nulls = [r for r in rows if conds[r["condition"]].data == "null" and conds[r["condition"]].metric == "norm_pearson"]
    out += ["## 2. False 'effect' on the null cards (N1-N5, N7, N8): (b) donors without gates against (d) the engine"]
    tot = {}
    for label, sel in (("all", lambda r: True), ("blind metric", lambda r: r["truth"] == "blind"),
                       ("metric at risk", lambda r: r["truth"] in RISK)):
        rs = [r for r in nulls if sel(r)]
        b = sum(mine[r["id"]]["b_call"] == "+" for r in rs)
        d = sum(r["outcome"] == P.SUPPORTED for r in rs)
        dd = sum(r["outcome"] == P.SUPPORTED and next((s["kept"] for s in sup if s["id"] == r["id"]), False) for r in rs)
        tot[label] = (b, d)
        out.append(f"  {label:<15} n {len(rs):>5}: (b) {fmt(b, len(rs))}; (d) {fmt(d, len(rs))}; "
                   f"(d) under the rule {fmt(dd, len(rs))}")
    red_all = tot["all"][0] - tot["all"][1]
    red_blind = tot["blind metric"][0] - tot["blind metric"][1]
    out.append(f"  the reduction from (b) to (d): {red_all} false effects, {red_blind} of them on blind metrics "
               f"({100 * red_blind / red_all:.0f}%)")
    out.append("  artifacts N2, N3, N8 by level, (b) -> (d):")
    for lv in P.LEVELS:
        rs = [r for r in nulls if r["condition"] in ("N2", "N3", "N8") and r["level"] == lv]
        b = sum(mine[r["id"]]["b_call"] == "+" for r in rs)
        d = sum(r["outcome"] == P.SUPPORTED for r in rs)
        risk = sum(r["truth"] in RISK for r in rs)
        out.append(f"    {lv:<6} n {len(rs):>4} (at risk {risk:>4}): (b) {fmt(b, len(rs))} -> (d) {fmt(d, len(rs))}")
    out.append("")

    # 3. absences the correction made
    nde_d = [r for r in rows if r["outcome"] == P.NDE]
    made_abs = [r for r in nde_d if mine[r["id"]]["c_outcome"] != P.NDE]
    out += ["## 3. NO DETECTABLE EFFECT made by the correction ((d) NDE where (c), uncorrected, is not)",
            f"  (d) NDE {len(nde_d)}; of them (c) not NDE {len(made_abs)}: "
            + ", ".join(f"{k} {v}" for k, v in sorted(Counter(
                f"{conds[r['condition']].data}/{r['truth']}" for r in made_abs).items())),
            "  of these, wrong by v1's rules (NDE outside the card's allowed set): "
            f"{sum(r['error'] for r in made_abs)}", ""]

    # 4. the exchange proposed for v2, illustrated on v1
    kept_sup = {s["id"] for s in sup if s["kept"]}

    def outcome_of(r, an):
        m = mine[r["id"]]
        if an == "b":
            return {"+": "effect", "-": "against", "0": "none"}[m["b_call"]]
        o = m["c_outcome"] if an == "c" else r["outcome"]
        if an == "d*" and o == P.SUPPORTED and r["id"] not in kept_sup:
            o = P.INCONCLUSIVE
        return {P.SUPPORTED: "effect", P.NS_OPPOSITE: "against", P.NDE: "absent"}.get(o, "none")

    def big(r):
        c = conds[r["condition"]]
        return c.data == "effect" and abs(float(P.delta_of(c.name, r["variant"], r["pair"], pilot)["value"])) >= P.sesoi_of(c, r["pair"], pilot)

    def kinds(r, an):
        """(false effect claim, false absence claim, unwarranted absence): an effect claimed in either
        direction on null data or against the true direction on a real one; an absence claimed on a real
        effect with |Δ*| >= SESOI or with a blind or useless metric; an absence claimed on null data with
        a blind or useless metric (true of the data, but not shown by the metric; v1 counted it an error)."""
        o, c = outcome_of(r, an), conds[r["condition"]]
        if c.data == "null":
            return o in ("effect", "against"), False, o == "absent" and r["truth"] not in RISK
        return o == "against", o == "absent" and (big(r) or r["truth"] not in RISK), False

    def wrong(r, an):
        fe, fa, _ = kinds(r, an)
        return fe or fa

    out += ["## 4. The exchange proposed for v2, on v1's cards (an illustration, after the results)",
            "  E: false definite conclusions per 1,000 cards (a false effect claim or a false absence claim); S: effects",
            "  found in the true direction per 1,000 real-effect cards. (b) makes no claim of absence (no TOST); (c) is",
            "  the engine's donor-level test with its TOST and no gates; d* is (d) under the rule of section 1.",
            "  Components (false effect / false absence / unwarranted absence on null data with a blind or useless",
            "  metric, which E leaves out):"]
    for an in ("b", "c", "d", "d*"):
        parts = []
        for label, rs in (("at risk", [r for r in rows if r["truth"] in RISK]),
                          ("blind or useless", [r for r in rows if r["truth"] not in RISK])):
            k = [sum(kinds(r, an)[i] for r in rs) for i in range(3)]
            parts.append(f"{label} {k[0]} / {k[1]} / {k[2]}")
        out.append(f"    {an:<2}: " + "; ".join(parts))
    risk_all = [r for r in rows if r["truth"] in RISK]
    risk_eff = [r for r in risk_all if conds[r["condition"]].data == "effect"]
    other = [r for r in rows if r["truth"] not in RISK]
    vals = {}
    for an in ("b", "c", "d", "d*"):
        e_r = sum(wrong(r, an) for r in risk_all)
        e_o = sum(wrong(r, an) for r in other)
        s_r = sum(outcome_of(r, an) == "effect" for r in risk_eff)
        vals[an] = (1000 * e_r / len(risk_all), 1000 * e_o / len(other), 1000 * s_r / len(risk_eff))
        out.append(f"  {an:<2}: at risk ({len(risk_all)} cards) E {e_r} = {vals[an][0]:.1f}; blind or useless "
                   f"({len(other)}) E {e_o} = {vals[an][1]:.1f}; S on the {len(risk_eff)} real effects at risk "
                   f"{s_r} = {vals[an][2]:.1f}")
    for ref in ("b", "c"):
        for an in ("d", "d*"):
            de, ds = vals[ref][0] - vals[an][0], vals[ref][2] - vals[an][2]
            out.append(f"  {an} against {ref}, at risk: errors removed {de:+.1f} per 1,000 cards, detections lost "
                       f"{ds:+.1f} per 1,000 real-effect cards; the rule 'removed >= lost / 4' "
                       f"{'holds' if de >= ds / 4 else 'fails'}; blind or useless: errors removed "
                       f"{vals[ref][1] - vals[an][1]:+.1f} per 1,000 cards")
    out.append("")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
