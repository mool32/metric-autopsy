"""Scoring of the confirmatory validation (validation/prereg/v1.md, sections 5 and 6).

    python validation/prereg/score.py --key KEY --results RESULTS --pilot pilot.json --out scores.json

Run by the workflow right after the blind run's results are committed. The key (the drand
round's randomness in the results' key record) re-derives the condition, variant and gene pair
of every dataset ID (``panel.assign``). Each claim card's report is reduced to its outcome, the
label and cause of the verdict (``panel.outcome``), and compared with the card's allowed
outcomes (``panel.allowed``: by the truth about the metric on the pair and about the data, one
rule for every condition). The primary outcomes and the criteria S1-S7 follow, every rate with
its two-sided 95% Clopper-Pearson interval, and for every criterion the design effect of
datasets that share donors (``overlap_interval``); a design effect above 1.5 is reported as the
pre-registered limitation. A missing report, an engine error or an unexpected verdict counts as
an error (S4) and fails S7. Thresholds are computed with the rules of ``oc.py`` on the realized
cards; S3's strata are judged at the tiers pilot.json fixed before the key ("s3_rules"). A
criterion without cards: S2, S5 and S6 hold (no card on which their error can occur); S1, S4 and
S7 always have cards; S3 fails unless every stratum is judged (what it shows is not shown).
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

    Datasets drawn from one background share donors, so their verdicts can be correlated, and
    donor-disjoint subsets would leave about 11 datasets per condition. Working model: the
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
    rule nominals). `reports`: {card id: (verdict text, cause[, gate record])} (None for a missing
    report or an engine error)."""
    conds = P.conditions()
    rows = []
    for e in entries:
        c = conds[e["condition"]]
        pair = int(e["pair"])
        level = P.pool_level(c, pair, pilot)
        ok = P.allowed(c, e["variant"], pair, pilot)
        good = P.definite(c, e["variant"], pair, pilot)
        truth = P.metric_truth(c, e["variant"], pair, pilot)
        model = oc.card_model(c, e["variant"], pair, pilot)
        for cid in P.card_ids(e):
            rec = reports.get(cid) or (None, None)
            verdict, cause = rec[0], rec[1]
            gates = rec[2] if len(rec) > 2 else {}
            o = P.outcome(verdict, cause)
            rows.append(dict(model, id=cid, dataset=e["id"], variant=e["variant"], level=level, pair=pair,
                             truth=truth, outcome=o,
                             sup_error=o == P.SUPPORTED and P.SUPPORTED not in ok,
                             false_invalid=truth == "valid" and o == P.NS_INVALID,
                             nde_error=o == P.NDE and P.NDE not in ok,
                             engine_error=o in (P.ERROR, P.OTHER),
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
                    engine_errors=rate(sum(r["engine_error"] for r in rs), len(rs)),
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

    def error_criterion(name, rs, field, must_have_cards=False):
        rule = rules[name]
        k = sum(r[field] for r in rs)
        passed = (k <= rule["max_allowed"]) if rule["judged"] else (not must_have_cards or bool(rs))
        crit[name] = dict(rate(k, len(rs)), nominal=rule["nominal"], rule_nominal=rule["rule_nominal"],
                          raised=rule["raised"], max_allowed=rule["max_allowed"], judged=rule["judged"],
                          p_pass_sound=rule["p_pass_sound"], p_pass_doubled=rule["p_pass_doubled"], passed=passed)
        members[name] = (rs, field)
    # S2: false SUPPORTED on every card where SUPPORTED is an error, at their mean rule nominal
    error_criterion("S2", [r for r in rows if r["sup_error_possible"]], "sup_error")
    # S3: correct definite outcomes on establishable cards, per stratum (real effects, nulls, blind
    # or useless metrics; not N3 and N8), each at its mean measured decisiveness (at most 0.85) and
    # at the tier fixed before the key (oc.criteria_rules). S3 holds if every stratum is judged and
    # holds: a stratum without a tier is decisiveness not shown.
    s3 = {}
    for st in oc.STRATA:
        rs = [r for r in rows if r["establishable"] and r["stratum"] == st]
        k = sum(r["correct_definite"] for r in rs)
        rule = rules["S3"][st]
        need = rule["min_required"]
        s3[st] = dict(rate(k, len(rs)), nominal=rule["nominal"], measured=rule["measured"], tier=rule["tier"],
                      min_required=need, judged=need is not None, passed=need is not None and k >= need,
                      tiers=rule["tiers"], definite_any=rate(sum(r["definite"] for r in rs), len(rs)))
        members[f"S3:{st}"] = (rs, "correct_definite")
    crit["S3"] = dict(strata=s3, tiers_fixed_before_the_key=fixed is not None,
                      passed=all(s3[st]["passed"] for st in oc.STRATA))
    # S4: outside the allowed outcomes, all cards, against the mean of the cards' rule nominals
    error_criterion("S4", rows, "error", must_have_cards=True)
    # S5: false "metric invalid" (GATE 4/5) on cards whose metric is valid by the truth
    error_criterion("S5", [r for r in rows if r["valid"]], "false_invalid")
    # S6: false NO DETECTABLE EFFECT where it is an error (a blind or useless metric; a real effect
    # at or above the SESOI)
    error_criterion("S6", [r for r in rows if r["nde_error_possible"]], "nde_error")
    # S7: no engine error, missing report or unexpected verdict on any card
    k = sum(r["engine_error"] for r in rows)
    crit["S7"] = dict(rate(k, len(rows)), max_allowed=0, passed=bool(rows) and k == 0)
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
                limitations.append(f"{name}: design effect {de[name]['deff']:.2f} > {DEFF_LIMIT} (shared donors)")
        out["secondary"]["shared_donors"] = de
        # S1's operating characteristics at the realized design effects (v1.md section 4): the
        # principle holds for independent datasets; with shared donors the results state them
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


def read_results(results: Path) -> tuple[dict, dict]:
    """(verdict, cause) by card id and donors by dataset id, from the blind run's merged results."""
    manifest = json.loads((results / "manifest.json").read_text())
    reports, donors = {}, {}
    for d in manifest["datasets"]:
        donors[d["id"]] = d.get("donors", [])
        for c in d["cards"]:
            rep_path = results / "reports" / f"{c['id']}.json"
            rep = json.loads(rep_path.read_text()) if rep_path.exists() else {}
            reports[c["id"]] = (None if "error" in rep or not rep else
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
    reports, donors = read_results(Path(args.results))
    entries = P.assign(args.key, pilot.get("dropped", ()), pilot.get("pool_size"))
    res = score(entries, reports, pilot, donors)
    res["key"] = args.key
    Path(args.out).write_text(json.dumps(res, indent=1))
    for s in oc.CRITERIA:
        print(s, "PASS" if res["criteria"][s]["passed"] else "FAIL")
    for st, v in res["criteria"]["S3"]["strata"].items():
        print(f"  S3 {st}: {v['k']}/{v['n']} correct definite"
              + (f", at least {v['min_required']} required ({v['tier']} tier)" if v["judged"] else
                 " — not judged: S3 fails"))
    for lim in res["limitations"]:
        print("limitation:", lim)
    print("validation", "PASSES" if res["passed"] else "FAILS")


if __name__ == "__main__":
    main()
