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
4. An illustration, on v1's cards, of the exchange of v2's primary comparison
   (validation/plans/v2_requirements.md; the owner's decision of 2026-10-12): per stratum, the false
   definite conclusions the gates remove per null card (r) and per effect card (r1), the true
   detections they lose per effect card (l), and pi* = r / (r + l / 4 - r1), the largest share of
   true effects among the claims at which they pay off, with bootstrap intervals over cards and over
   gene pairs; for (d) and (d) under the rule, against (b) and (c).
5. The gene pairs at risk v2 needs for the interval of pi* over pairs to be within +-0.1 (the owner,
   2026-10-12): v1's pairs resampled K at a time, and the square-root rule from v1's own interval.
6. An illustration for the real claims' second question (validation/plans/real_claims.md): the
   votes between (d) and the donor-level test with its TOST ((c)) on v1's artifacts (as negative
   controls) and real effects (as positive ones), under two rules of counting, per dataset of m
   claims; and the control datasets needed for 13 and 15 deciding votes.

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

    def unit(r):
        """The gene pair a card is drawn on (background, pool index): the cluster of the bootstrap."""
        return f"{conds[r['condition']].background}:{r['pair']}"

    def pi_star(r0, r1, lost):
        """The largest share of true effects at which the gates pay off, (1 - pi) r0 + pi r1 >= pi lost / 4:
        1.0 where they pay off at every share, 0.0 where at none, nan where only above a share."""
        a, b = r0, lost / 4 - r1
        if a + b > 0:
            return min(max(a / (a + b), 0.0), 1.0)
        return 1.0 if a >= 0 else float("nan")

    rng = np.random.default_rng(20261012)
    n_boot = 2000
    out += ["## 4. The exchange for v2, on v1's cards (an illustration, after the results)",
            "  The owner's decision of 2026-10-12: one false conclusion is worth four lost detections. With a share pi",
            "  of true effects among the claims, the gates pay off against a comparator where",
            "  (1 - pi) * r + pi * r1 >= pi * l / 4: r the false conclusions they remove per null card, r1 per effect",
            "  card, l the true detections they lose per effect card; pi* = r / (r + l / 4 - r1), the largest share of",
            "  true effects at which they pay off. A false conclusion: an effect claimed in either direction on null",
            "  data or against the true direction on a real one; an absence claimed on a real effect with",
            "  |Δ*| >= SESOI or with a blind or useless metric. An absence claimed on null data with a blind or useless",
            "  metric is true of the data but not shown by the metric: counted apart as unwarranted. A detection: an",
            "  effect found in the true direction. (b) makes no claim of absence (no TOST); (c) is the engine's",
            "  donor-level test with its TOST and no gates; d* is (d) under the rule of section 1. Intervals: 95%",
            f"  percentile bootstrap, {n_boot} resamples each, over cards (null and effect cards apart) and over gene",
            "  pairs (background, pool index; all of a pair's cards together); the donors that datasets share are",
            "  not resampled."]
    for label, rs in (("at risk", [r for r in rows if r["truth"] in RISK]),
                      ("blind or useless", [r for r in rows if r["truth"] not in RISK])):
        is_null = np.array([conds[r["condition"]].data == "null" for r in rs])
        units = sorted({unit(r) for r in rs})
        uidx = np.array([units.index(unit(r)) for r in rs])
        n0, n1 = int(is_null.sum()), int((~is_null).sum())
        rec = {an: (np.array([wrong(r, an) for r in rs], dtype=float),
                    np.array([not z and outcome_of(r, an) == "effect" for r, z in zip(rs, is_null)], dtype=float))
               for an in ("b", "c", "d", "d*")}
        out.append(f"  {label}: {n0} null cards, {n1} effect cards, {len(units)} gene pairs. False conclusions on "
                   "null cards / on effect cards; detections; unwarranted absences:")
        for an, (e, d) in rec.items():
            out.append(f"    {an:<2}: {int(e[is_null].sum())} / {int(e[~is_null].sum())}; {int(d.sum())}; "
                       f"{sum(kinds(r, an)[2] for r in rs)}")
        i0, i1 = np.flatnonzero(is_null), np.flatnonzero(~is_null)
        for ref in ("b", "c"):
            for an in ("d", "d*"):
                dE, dS = rec[ref][0] - rec[an][0], rec[ref][1] - rec[an][1]
                e0, e1, ls = int(dE[i0].sum()), int(dE[i1].sum()), int(dS[i1].sum())
                est = pi_star(e0 / n0, e1 / n1, ls / n1)
                card = []
                for _ in range(n_boot):
                    j0, j1 = rng.choice(i0, n0), rng.choice(i1, n1)
                    card.append(pi_star(dE[j0].sum() / n0, dE[j1].sum() / n1, dS[j1].sum() / n1))
                agg = np.zeros((len(units), 5))
                np.add.at(agg, uidx, np.column_stack([is_null, ~is_null, dE * is_null, dE * ~is_null, dS]))
                pair = []
                for _ in range(n_boot):
                    t = np.bincount(rng.integers(0, len(units), len(units)), minlength=len(units)) @ agg
                    pair.append(pi_star(t[2] / t[0], t[3] / t[1], t[4] / t[1]) if t[0] and t[1] else float("nan"))

                def ci(v):
                    v = np.asarray(v)
                    lo, hi = np.nanpercentile(v, [2.5, 97.5])
                    nan = int(np.isnan(v).sum())
                    return f"{lo:.3f}-{hi:.3f}" + (f", {nan} resamples undefined" if nan else "")

                pays = e0 + e1 >= ls / 4
                out.append(f"    {an} against {ref}: r {e0}/{n0} = {e0 / n0:.4f}, r1 {e1}/{n1} = {e1 / n1:.4f}, l {ls}/{n1} "
                           f"= {ls / n1:.4f}; pi* {est:.3f} (cards {ci(card)}; gene pairs {ci(pair)}); at v1's share "
                           f"pi = {n1}/{len(rs)} = {n1 / len(rs):.3f}: removed {e0 + e1} against lost / 4 = {ls / 4:.1f}, "
                           f"{'pays off' if pays else 'does not pay off'}")
    # 5. the gene pairs v2 needs for pi* within +-0.1 (the owner, 2026-10-12)
    rs = [r for r in rows if r["truth"] in RISK]
    is_null = np.array([conds[r["condition"]].data == "null" for r in rs])
    units = sorted({unit(r) for r in rs})
    uidx = np.array([units.index(unit(r)) for r in rs])
    rec = {an: (np.array([wrong(r, an) for r in rs], dtype=float),
                np.array([not z and outcome_of(r, an) == "effect" for r, z in zip(rs, is_null)], dtype=float))
           for an in ("b", "c", "d", "d*")}
    grid = (18, 25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 400, 500)
    out += ["## 5. Gene pairs at risk needed for the interval of pi* over pairs to be within +-0.1",
            f"  The cards at risk drawn on K gene pairs resampled with replacement from v1's {len(units)} (a pair's",
            f"  cards together, {n_boot} resamples per K): the half-width of the 95% percentile interval of pi*, and",
            "  the K a square-root rule gives from v1's own half-width (K = 18 * (half-width / 0.1)^2)."]
    for ref, an in (("b", "d"), ("b", "d*"), ("c", "d"), ("c", "d*")):
        dE, dS = rec[ref][0] - rec[an][0], rec[ref][1] - rec[an][1]
        agg = np.zeros((len(units), 5))
        np.add.at(agg, uidx, np.column_stack([is_null, ~is_null, dE * is_null, dE * ~is_null, dS]))
        half, first = {}, None
        for k in grid:
            v = []
            for _ in range(n_boot):
                t = np.bincount(rng.integers(0, len(units), k), minlength=len(units)) @ agg
                v.append(pi_star(t[2] / t[0], t[3] / t[1], t[4] / t[1]) if t[0] and t[1] else float("nan"))
            lo, hi = np.nanpercentile(np.asarray(v), [2.5, 97.5])
            half[k] = (hi - lo) / 2
            if first is None and half[k] <= 0.1:
                first = k
        sq = len(units) * (half[len(units)] / 0.1) ** 2
        out.append(f"    {an} against {ref}: half-width " + ", ".join(f"K {k} {half[k]:.3f}" for k in grid)
                   + f"; the first K within 0.1: {first}; the square-root rule: {sq:.0f}")
    out.append("")

    # 6. the votes a control dataset gives, on v1's cards (an illustration for the real claims' second question)
    def judged(r, an):
        """+1 a conclusion the truth confirms, -1 one it rejects, 0 none (a refusal, INCONCLUSIVE). Null data,
        as a negative control by construction: an effect claimed either way is rejected, an absence confirmed.
        A real effect, as a positive control: the true direction confirmed, the other rejected, an absence
        rejected where |Δ*| >= SESOI or the metric is blind or useless, confirmed otherwise."""
        o, c = outcome_of(r, an), conds[r["condition"]]
        if o == "none":
            return 0
        if c.data == "null":
            return 1 if o == "absent" else -1
        if o == "absent":
            return -1 if (big(r) or r["truth"] not in RISK) else 1
        return 1 if o == "effect" else -1

    def vote(r, effects_only=False):
        """+1 a vote for the engine (d), -1 for the donor-level test with its TOST ((c) here), 0 none. With
        effects_only a difference decides only where one side claims an effect (either way), not where an
        absence meets no conclusion."""
        oc_, od = outcome_of(r, "c"), outcome_of(r, "d")
        if oc_ == od or (effects_only and {oc_, od} <= {"absent", "none"}):
            return 0
        return int(np.sign(judged(r, "d") - judged(r, "c")))

    controls = (("negative (N2, N3, N8: artifacts, no effect)", [r for r in rows if r["condition"] in ("N2", "N3", "N8")]),
                ("positive (E1, E2, E3: real effects)", [r for r in rows if conds[r["condition"]].data == "effect"]))
    out += ["## 6. Votes a control dataset gives, on v1's cards (an illustration, after the results)",
            "  A vote compares (d) with (c), which stands in for the real track's donor-level test with its TOST,",
            "  where their conclusions differ: for (d) if it avoids a conclusion the truth rejects or finds one it",
            "  confirms, for the test if the reverse, none where both are equally right or wrong. The null cards",
            "  of the artifacts stand in for negative controls (truth by construction: an absence is true), the real",
            "  effects for positive ones. Two rules: every difference of conclusion decides (as decided), or only",
            "  one where a side claims an effect (an absence against no conclusion does not decide). A dataset",
            "  with m claims is m cards drawn at random from the stratum (v1's cards come from different datasets,",
            "  so the claims of one real dataset may agree more); its vote is the side of the majority of its card",
            "  votes, none on a tie. 20000 datasets per m."]
    q = {}
    for rule, eo in (("every difference decides", False), ("effect claims only", True)):
        out.append(f"  rule: {rule}")
        votes = {}
        for label, rs in controls:
            v = np.array([vote(r, eo) for r in rs])
            votes[label] = v
            out.append(f"    {label}: {len(rs)} cards; per card for (d) {int((v > 0).sum())}, for the test "
                       f"{int((v < 0).sum())}, none {int((v == 0).sum())}")
            for m in (1, 8, 24):
                tot = v[rng.integers(0, len(v), (20000, m))].sum(axis=1)
                dec = tot != 0
                q[(rule, label, m)] = (dec.mean(), (tot > 0).sum() / max(dec.sum(), 1))
                out.append(f"      m = {m:>2}: a deciding vote in {dec.mean():.3f} of the datasets; for (d) in "
                           f"{(tot > 0).sum() / max(dec.sum(), 1):.3f} of those")
        for m in (8, 24):
            (qn, pn), (qp, pp) = q[(rule, controls[0][0], m)], q[(rule, controls[1][0], m)]
            out.append(f"    half negative and half positive datasets, m = {m}: for (d) in "
                       f"{(qn * pn + qp * pp) / (qn + qp):.3f} of the deciding votes")
    out += ["  Datasets needed for 13 and 15 deciding votes, half negative and half positive controls, after the",
            "  votes of N papers of which a share f change the conclusion and half of those meet a deciding verdict:"]
    for rule in ("every difference decides", "effect claims only"):
        out.append(f"  rule: {rule}")
        for m in (1, 8, 24):
            qbar = (q[(rule, controls[0][0], m)][0] + q[(rule, controls[1][0], m)][0]) / 2
            parts = []
            for npap, f in ((20, 0.2), (20, 0.3), (30, 0.2), (30, 0.3)):
                pv = npap * f * 0.5
                parts.append(f"N {npap}, f {f}: papers {pv:.1f} votes, controls "
                             + " / ".join(str(int(np.ceil(max(t - pv, 0) / qbar))) for t in (13, 15)))
            out.append(f"    m = {m:>2} (a deciding vote in {qbar:.3f} of control datasets): " + "; ".join(parts))
    out.append("")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
