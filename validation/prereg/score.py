"""Scoring of the confirmatory validation (validation/prereg/v1.md, sections 5 and 6).

    python validation/prereg/score.py --key KEY --results RESULTS --pilot pilot.json --out scores.json

Run by the workflow right after the blind run's results are committed. The key (the drand
round's randomness in the results' key record) re-derives the condition, variant and gene pair
of every dataset ID (``panel.assign``). Each claim card's report is reduced to its outcome, the
label and cause of the verdict (``panel.outcome``), and compared with the card's allowed
outcomes (``panel.card_allowed``: by the truth about the metric on the pair and about the data,
one rule for every condition, without the effect verdicts the engine's rules exclude on the card). The primary outcomes and the criteria S1-S7 follow, every rate with
its two-sided 95% Clopper-Pearson interval (exact: given the background every dataset is an
independent draw, so the counts are binomial), and for every criterion the design effect of
datasets that share donors (``overlap_interval``: how much the outcomes depend on the donors a
dataset draws, so how far a rate carries beyond the background's donors); a design effect above
1.5 is reported as the pre-registered limitation on that reading. A verdict the engine's
deterministic rules cannot give there (``oc.S7_OUTCOMES``: a label that does not match its cause,
a cause the panel does not expect, UNIDENTIFIABLE; DEGENERATE METRIC on a metric that varies; an
effect verdict where its rules exclude one) is a rule violation: an error (P2) that fails S7a, none
allowed. An engine exception or a missing report is a crash (no verdict): counted in S7b, at most
the principle's threshold at a nominal 0.1% (the fifth round), and every crash is listed with its
traceback in scores.json ("crashes"). Every other wrong outcome is counted in its cell of S2, S4,
S5 or S6 (``oc.error_cells``). Thresholds are computed with the rules of ``oc.py`` on the realized
cards; S3's strata are judged at the tiers pilot.json fixed before the key ("s3_rules"). A cell
without cards holds (no card on which its error can occur); S1, S7a and S7b always have cards; S3
fails unless every stratum is judged (what it shows is not shown).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import stats

import oc
import panel as P

DEFF_LIMIT = 1.5  # a criterion's design effect above this is a pre-registered limitation
# the flag of a card's row that counts the outcome of a cell, by the cell's field (oc.FIELDS)
FLAG = {"sup": "sup_error", "inv": "false_invalid", "nde": "nde_error", "opp": "opp_error", "depth": "depth_error",
        "oppnull": "opp_null"}


def cp(k: float, n: float, alpha: float = 0.05) -> tuple[float, float]:
    """Clopper-Pearson; with a non-integer (effective) k and n it is the Korn-Graubard interval."""
    lo = 0.0 if k <= 0 else float(stats.beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k >= n else float(stats.beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


def rate(k: int, n: int) -> dict:
    lo, hi = cp(k, n) if n else (float("nan"), float("nan"))
    return dict(k=int(k), n=int(n), rate=k / n if n else float("nan"), ci95=[lo, hi])


def overlap_interval(y, donor_sets: list, alpha: float = 0.05) -> dict:
    """A rate's interval that allows for datasets sharing donors (chosen before the key).

    Given the background the datasets are independent draws and the Clopper-Pearson interval is
    exact for the rate on it (the third review); but datasets share donors and their outcomes can
    depend on the donors drawn, which bounds how far the rate carries to other donors of the same
    kind (donor-disjoint subsets would leave about 11 datasets per condition). This interval is for
    that reading, secondary. Working model: the
    covariance of two datasets' outcomes grows with the number of donors they share,
    Cov(y_j, y_k) = beta * |S_j & S_k| (an additive donor effect). beta is the moment estimate
    over all pairs of datasets, truncated at 0; the variance of the rate is then
    (sum of squared residuals + beta * sum of overlaps) / n^2, its ratio to the binomial variance
    the design effect, and the interval Clopper-Pearson at the effective sample size n / deff
    (Korn and Graubard 1998)."""
    y = np.asarray(y, float)
    n = len(y)
    if n == 0:
        return dict(rate=float("nan"), n=0, deff=float("nan"), n_eff=0.0, beta=0.0, ci95=[float("nan")] * 2)
    p = float(y.mean())
    donors = sorted({d for s in donor_sets for d in s})
    at = {d: i for i, d in enumerate(donors)}
    M = np.zeros((n, len(donors)), dtype=np.float32)
    for j, s in enumerate(donor_sets):
        M[j, [at[d] for d in s]] = 1.0
    O = M @ M.T
    np.fill_diagonal(O, 0.0)
    e = y - p
    den = float((O.astype(np.float64) ** 2).sum())
    beta = max(0.0, float(e @ (O @ e)) / den) if den > 0 else 0.0
    var = (float(e @ e) + beta * float(O.sum())) / n ** 2
    vb = p * (1 - p) / n
    deff = max(1.0, var / vb) if vb > 0 else 1.0
    n_eff = n / deff
    lo, hi = cp(p * n_eff, n_eff, alpha)
    return dict(rate=p, n=n, deff=deff, n_eff=n_eff, beta=beta, ci95=[lo, hi])


def card_rows(entries: list[dict], reports: dict, pilot: dict) -> list[dict]:
    """One row per claim card: its condition, variant, level and truth, its outcome and how it is
    scored (with ``oc.card_model``: where each error can occur, the sound validator's rates and the
    rule nominals). `reports`: {card id: (verdict text, cause[, gate record])}; a card without one
    (no report) or with a dict (the engine's exception: its error and traceback, ``read_results``)
    is a crash."""
    conds = P.conditions()
    rows = []
    for e in entries:
        c = conds[e["condition"]]
        pair = int(e["pair"])
        level = P.pool_level(c, pair, pilot)
        truth = P.metric_truth(c, e["variant"], pair, pilot)
        model = oc.card_model(c, e["variant"], pair, pilot)
        for k, cid in enumerate(P.card_ids(e)):
            ok = P.card_allowed(c, e["variant"], pair, pilot, k)  # without what the engine's rules exclude
            good = ok & P.DEFINITE
            kind = oc.card_kind(model, c, k)
            rec = reports.get(cid)
            if rec is None or isinstance(rec, dict):  # no report, or the engine's exception: no verdict
                o, gates, crash = P.CRASH, {}, dict(rec or {"error": "no report", "traceback": []})
            else:
                gates, crash = (rec[2] if len(rec) > 2 else {}), None
                o = P.outcome(rec[0], rec[1])
            rows.append(dict(kind, id=cid, dataset=e["id"], variant=e["variant"], level=level,
                             pair=pair, truth=truth, outcome=o,
                             sup_error=o == P.SUPPORTED and P.SUPPORTED not in ok,
                             false_invalid=truth == "valid" and o == P.NS_INVALID,
                             nde_error=o == P.NDE and P.NDE not in ok,
                             opp_error=o == P.NS_OPPOSITE and P.NS_OPPOSITE not in ok,
                             opp_null=o == P.NS_OPPOSITE and kind["opp_null_counted"],
                             depth_error=o == P.NS_DEPTH and P.NS_DEPTH not in ok,
                             rule_violation=(o in oc.S7_OUTCOMES or (o == P.DEGENERATE and P.DEGENERATE not in ok)
                                             or o in P.excluded(c, k)),
                             crashed=o == P.CRASH, crash=crash,
                             error=o not in ok, correct_definite=o in good, definite=o in P.DEFINITE,
                             gates=gates))
    return rows


def gate_record(rep: dict) -> dict:
    """What P5 reads from a report: GATE 0's refusal, GATE 1's flag, GATE 4's outcome, GATE 5's
    status and negative control, and the effect's diagnosis and retained share."""
    g = {r.get("gate"): r for r in rep.get("gates", [])}
    ef = (rep.get("fields") or {}).get("effect") or {}
    rows5 = (g.get(5, {}).get("detail") or {}).get("rows") or []
    return dict(gate0=g.get(0, {}).get("status"), gate1=g.get(1, {}).get("status"),
                gate4=(g.get(4, {}).get("detail") or {}).get("outcome") or g.get(4, {}).get("status"),
                gate5=g.get(5, {}).get("status"), gate5_negative_fails=any(not r.get("neg_ok", True) for r in rows5),
                explained_by_depth=bool((ef.get("detail") or {}).get("explained_by_depth")),
                retained=(ef.get("detail") or {}).get("retained"))


def per_gate(rows: list[dict]) -> dict:
    """P5: per-gate sensitivity and specificity (v1.md, section 5)."""
    def share(rs, fn):
        return rate(sum(bool(fn(r)) for r in rs), len(rs))

    def by(keyf, rs):
        out = {}
        for r in rs:
            out.setdefault(keyf(r), []).append(r)
        return out
    known = [r for r in rows if r["gates"]]
    out = dict(
        gate1_flags={k: share(rs, lambda r: r["gates"].get("gate1") in ("WARN", "STOP"))
                     for k, rs in by(lambda r: f"{r['condition']}:{r['variant']}",
                                     [r for r in known if r["condition"] in ("N1", "N2", "N3")]).items()},
        explained_by_depth={k: share(rs, lambda r: r["gates"].get("explained_by_depth"))
                            for k, rs in by(lambda r: f"{r['variant']}:{r['level']}",
                                            [r for r in known if r["condition"] == "N2"]).items()},
        gate4={k: {o: sum(r["gates"].get("gate4") == o for r in rs) for o in ("PASS", "FAIL", "UNTESTED", "SKIP")}
               for k, rs in by(lambda r: f"{r['truth']}:{r['level']}", known).items()},
        gate5_negative_control_fails=share([r for r in known if r["condition"] == "N1"],
                                           lambda r: r["gates"].get("gate5_negative_fails")),
        gate0_refusals={k: share(rs, lambda r: r["outcome"] == P.REFUSAL)
                        for k, rs in by(lambda r: r["truth"], rows).items()},
        retained={k: dict(n=len(v), median=float(np.median(v)), q25=float(np.quantile(v, 0.25)),
                          q75=float(np.quantile(v, 0.75)))
                  for k, rs in by(lambda r: f"{r['condition']}:{r['variant']}",
                                  [r for r in known if r["condition"] in ("E1", "E2", "E3")]).items()
                  for v in [[float(r["gates"]["retained"]) for r in rs
                             if r["gates"].get("retained") is not None and np.isfinite(r["gates"]["retained"])]] if v})
    return out


def score(entries: list[dict], reports: dict, pilot: dict, donors: dict | None = None, sims: int = 20_000) -> dict:
    """entries: panel.assign(key, dropped, pool sizes) (or a subset); reports: {card id: (verdict,
    cause)}; donors: {dataset id: [donor ids]} from the manifest (for the design effects); sims: the
    simulations of a sound validator for the joint probability (``oc.criteria_rules``)."""
    rows = card_rows(entries, reports, pilot)
    out = dict(per_condition={}, per_level={}, per_truth={}, criteria={}, secondary={})

    def summary(rs):
        return dict(n=len(rs), errors=rate(sum(r["error"] for r in rs), len(rs)),
                    false_supported=rate(sum(r["sup_error"] for r in rs), sum(r["sup_error_possible"] for r in rs)),
                    false_invalid=rate(sum(r["false_invalid"] for r in rs), sum(r["valid"] for r in rs)),
                    false_nde=rate(sum(r["nde_error"] for r in rs), sum(r["nde_error_possible"] for r in rs)),
                    false_opposite=rate(sum(r["opp_error"] for r in rs), sum(r["opp_error_possible"] for r in rs)),
                    false_depth=rate(sum(r["depth_error"] for r in rs), sum(r["depth_error_possible"] for r in rs)),
                    opposite_on_null=rate(sum(r["opp_null"] for r in rs), sum(r["opp_null_counted"] for r in rs)),
                    rule_violations=rate(sum(r["rule_violation"] for r in rs), len(rs)),
                    crashes=rate(sum(r["crashed"] for r in rs), len(rs)),
                    correct_definite=rate(sum(r["correct_definite"] for r in rs if r["establishable"]),
                                          sum(r["establishable"] for r in rs)),
                    outcomes={o: rate(sum(r["outcome"] == o for r in rs), len(rs)) for o in P.OUTCOMES})

    groups = {}
    for r in rows:
        groups.setdefault(("per_condition", f"{r['condition']}:{r['variant']}"), []).append(r)
        groups.setdefault(("per_level", f"{r['condition']}:{r['variant']}:{r['level']}"), []).append(r)
        groups.setdefault(("per_truth", r["truth"]), []).append(r)
    for (where, key), rs in groups.items():
        out[where][key] = summary(rs)
    out["pooled"] = summary(rows)  # P1-P4 and P6-P7 pooled over every card (reported)
    crit = out["criteria"]
    members = {}
    fixed = (pilot.get("s3_rules") or {}).get("tiers")
    rules = oc.criteria_rules(rows, sims=sims, s3_tiers=fixed)
    # S1: false SUPPORTED at every key null condition's key variant; passes only if all pass
    s1 = {}
    for c in P.CONDITIONS:
        if not c.key:
            continue
        rs = [r for r in rows if r["condition"] == c.name and r["variant"] == c.key_variant]
        k = sum(r["sup_error"] for r in rs)
        kmax = rules["S1"][c.name]["max_allowed"]
        s1[c.name] = dict(rate(k, len(rs)), max_allowed=kmax, passed=bool(rs) and k <= kmax)
        members[f"S1:{c.name}"] = (rs, "sup_error")
    crit["S1"] = dict(conditions=s1, passed=bool(s1) and all(v["passed"] for v in s1.values()))

    # S2, S4, S5, S6: every wrong outcome in its cell (oc.error_cells): a false SUPPORTED per group of
    # the cards outside S1; a false NOT SUPPORTED against the direction or explained by depth, a false
    # "metric invalid" and a false NO DETECTABLE EFFECT, each per stratum. A criterion holds if every
    # cell holds; a cell without cards holds (no card on which its error can occur).
    for name in oc.ERROR_CRITERIA:
        cells = {}
        for cell, rule in rules["cells"].items():
            if cell.split(":")[0] != name:
                continue
            rs = [rows[i] for i in rule["members"]]
            k = sum(r[FLAG[rule["field"]]] for r in rs)
            cells[cell] = dict(rate(k, len(rs)), outcome=rule["outcome"], nominal=rule["nominal"],
                               rule_nominal=rule["rule_nominal"], raised=rule["raised"], judged=rule["judged"],
                               max_allowed=rule["max_allowed"], p_pass_sound=rule["p_pass_sound"],
                               p_pass_doubled=rule["p_pass_doubled"], p_pass_measured=rule.get("p_pass_measured"),
                               passed=(k <= rule["max_allowed"]) if rule["judged"] else True)
            members[cell] = (rs, FLAG[rule["field"]])
        crit[name] = dict(cells=cells, passed=all(v["passed"] for v in cells.values()))
    # S3: correct definite outcomes on establishable cards, per stratum (real effects, nulls, blind
    # or useless metrics; not N3 and N8), each at its mean measured decisiveness (at most 0.85) and
    # at the tier fixed before the key (oc.criteria_rules). S3 holds if every stratum is judged and
    # holds: a stratum without a tier is decisiveness not shown, and one with a valid metric where
    # GATE 0 refused more often than oc.GATE0_BOUND before the key fails (the fifth round).
    s3, over = {}, oc.refusal_bound_failures(pilot)  # GATE 0 above its bound before the key: S3 fails there
    for st in oc.STRATA:
        rs = [r for r in rows if r["establishable"] and r["stratum"] == st]
        k = sum(r["correct_definite"] for r in rs)
        rule = rules["S3"][st]
        need = rule["min_required"]
        s3[st] = dict(rate(k, len(rs)), nominal=rule["nominal"], measured=rule["measured"], tier=rule["tier"],
                      min_required=need, judged=need is not None,
                      passed=need is not None and k >= need and st not in over,
                      gate0_refusal_bound_failed=st in over,
                      tiers=rule["tiers"], definite_any=rate(sum(r["definite"] for r in rs), len(rs)))
        members[f"S3:{st}"] = (rs, "correct_definite")
    crit["S3"] = dict(strata=s3, tiers_fixed_before_the_key=fixed is not None, gate0_refusal_bound_failed=over,
                      passed=all(s3[st]["passed"] for st in oc.STRATA))
    # S7a: no verdict the engine's deterministic rules cannot give, on any card
    k = sum(r["rule_violation"] for r in rows)
    crit["S7a"] = dict(rate(k, len(rows)), max_allowed=0, passed=bool(rows) and k == 0)
    members["S7a"] = (rows, "rule_violation")
    # S7b: crashes (an engine exception or no report) at most the threshold of the nominal 0.1%
    k, rule = sum(r["crashed"] for r in rows), rules["S7b"]
    crit["S7b"] = dict(rate(k, len(rows)), nominal=rule["nominal"], rule_nominal=rule["rule_nominal"],
                       raised=rule["raised"], judged=rule["judged"], max_allowed=rule["max_allowed"],
                       p_pass_sound=rule["p_pass_sound"], p_pass_doubled=rule["p_pass_doubled"],
                       passed=bool(rows) and k <= rule["max_allowed"])
    members["S7b"] = (rows, "crashed")
    out["crashes"] = [dict(card=r["id"], dataset=r["dataset"], condition=r["condition"], variant=r["variant"],
                           **r["crash"]) for r in rows if r["crashed"]]
    out["passed"] = all(crit[s]["passed"] for s in oc.CRITERIA)
    out["n_cards"] = len(rows)
    out["joint_pass_probability_sound"] = rules["joint"]
    out["joint_pass_probability_measured"] = rules["joint_measured"]
    out["per_gate"] = per_gate(rows)
    # the design effect of datasets sharing donors, for every primary criterion
    limitations = []
    if donors:
        de = {}
        for name, (rs, field) in members.items():
            de[name] = overlap_interval([r[field] for r in rs], [donors.get(r["dataset"], []) for r in rs])
            if de[name]["deff"] > DEFF_LIMIT:
                limitations.append(f"{name}: design effect {de[name]['deff']:.2f} > {DEFF_LIMIT} (shared donors: "
                                   "the rate is read beyond the background's donors with care)")
        out["secondary"]["shared_donors"] = de
        # S1's operating characteristics at the realized design effects (v1.md section 4): exact
        # given the background; read beyond its donors, the results state them at these effects
        per = {}
        for c, v in s1.items():
            d = de.get(f"S1:{c}", {}).get("deff", 1.0)
            d = 1.0 if not np.isfinite(d) else d
            _, ps, pd = oc.error_rule_deff(v["n"], oc.E_SUPPORTED, d) if v["n"] else (-1, float("nan"), float("nan"))
            per[c] = dict(deff=d, p_pass_sound=ps, p_pass_doubled=pd)
        out["secondary"]["s1_at_realized_deff"] = dict(
            conditions=per, p_all_pass_sound=float(np.prod([v["p_pass_sound"] for v in per.values()])) if per else float("nan"))
    out["limitations"] = limitations
    return out


def read_results(results: Path, key: str | None = None) -> tuple[dict, dict]:
    """(verdict, cause, gate record) by card id — for an engine exception its error and traceback
    (a dict), for a missing report nothing — and donors by dataset id, from the blind run's merged
    results; with `key`, the results must be the run of that key (the third review: the key given
    to the scoring was never compared with the manifest's)."""
    manifest = json.loads((results / "manifest.json").read_text())
    if key is not None and manifest.get("key") != key:
        raise SystemExit("the key given is not the key of these results (manifest.json)")
    reports, donors = {}, {}
    for d in manifest["datasets"]:
        donors[d["id"]] = d.get("donors", [])
        for c in d["cards"]:
            rep_path = results / "reports" / f"{c['id']}.json"
            if not rep_path.exists():
                continue  # a crash: no report
            rep = json.loads(rep_path.read_text())
            reports[c["id"]] = (dict(error=rep.get("error") or "an empty report", traceback=rep.get("traceback") or [])
                                if "error" in rep or not rep else
                                (rep.get("verdict"), rep.get("cause"), gate_record(rep)))
    return reports, donors


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--key", required=True, help="the key: the drand round's randomness (64 hex characters)")
    p.add_argument("--results", required=True)
    p.add_argument("--pilot", required=True)
    p.add_argument("--out", default="scores.json")
    args = p.parse_args(argv)
    pilot = json.loads(Path(args.pilot).read_text())
    if not (pilot.get("s3_rules") or {}).get("tiers"):  # the fourth review: chosen after the key otherwise
        raise SystemExit("pilot.json has no S3 tiers fixed before the key (oc.py --pilot ... --write-judged)")
    reports, donors = read_results(Path(args.results), P.check_key(args.key))
    entries = P.assign(args.key, pilot.get("dropped", ()), pilot.get("pool_size"))
    res = score(entries, reports, pilot, donors)
    res["key"] = args.key
    Path(args.out).write_text(json.dumps(res, indent=1))
    for s in oc.CRITERIA:
        print(s, "PASS" if res["criteria"][s]["passed"] else "FAIL")
        for cell, v in (res["criteria"][s].get("cells") or {}).items():
            if not v["passed"]:
                print(f"  {cell}: {v['k']}/{v['n']}, at most {v['max_allowed']} allowed")
    for st, v in res["criteria"]["S3"]["strata"].items():
        print(f"  S3 {st}: {v['k']}/{v['n']} correct definite"
              + (f", at least {v['min_required']} required ({v['tier']} tier)" if v["judged"] else
                 " — not judged: S3 fails")
              + (f"; GATE 0 refused above {oc.GATE0_BOUND:g} before the key: S3 fails" if v["gate0_refusal_bound_failed"] else ""))
    v = res["criteria"]["S7b"]
    print(f"  S7b: {v['k']} crashes of {v['n']} cards, at most {v['max_allowed']} allowed")
    for c in res["crashes"]:
        print(f"  crash {c['card']} ({c['condition']}:{c['variant']}): {c['error']}")
        for line in c.get("traceback", []):
            print(f"      {line}")
    for lim in res["limitations"]:
        print("limitation:", lim)
    print("validation", "PASSES" if res["passed"] else "FAILS")


if __name__ == "__main__":
    main()
