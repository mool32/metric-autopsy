"""Scoring of the confirmatory validation (validation/prereg/v1.md, sections 5 and 6).

    python validation/prereg/score.py --key KEY --results RESULTS --pilot pilot.json --out scores.json

Run by the workflow right after the blind run's results are committed. The key (the drand
round's randomness in the results' key record) re-derives the condition, variant and gene pair
of every dataset ID (``panel.assign``). Each claim card's report is reduced to its outcome, the
label and cause of the verdict (``panel.outcome``), and compared with the card's allowed
outcomes (``panel.allowed``: by the truth about the metric on the pair and about the data, one
rule for every condition). The primary outcomes and the criteria S1-S5 follow, every rate with
its two-sided 95% Clopper-Pearson interval, and for every criterion the design effect of
datasets that share donors (``overlap_interval``); a design effect above 1.5 is reported as the
pre-registered limitation. A missing report or an engine error counts as an error.
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
    scored. `reports`: {card id: (verdict text, cause[, gate record])} (None for a missing report
    or an engine error)."""
    conds = P.conditions()
    est = pilot.get("establishable", {})
    rows = []
    for e in entries:
        c = conds[e["condition"]]
        pair = int(e["pair"])
        level = P.pool_level(c, pair, pilot)
        ok = P.allowed(c, e["variant"], pair, pilot)
        good = P.definite(c, e["variant"], pair, pilot)
        truth = P.metric_truth(c, pair, pilot)
        establishable = c.oracle and bool(est.get(f"{c.name}:{e['variant']}:{pair}", {}).get("establishable"))
        nominal = oc.nominal_error(ok)
        for cid in P.card_ids(e):
            rec = reports.get(cid) or (None, None)
            verdict, cause = rec[0], rec[1]
            gates = rec[2] if len(rec) > 2 else {}
            o = P.outcome(verdict, cause)
            rows.append(dict(id=cid, dataset=e["id"], condition=c.name, variant=e["variant"], level=level,
                             pair=pair, truth=truth, outcome=o,
                             key=c.key and e["variant"] == c.key_variant,
                             sup_error_possible=P.SUPPORTED not in ok, sup_error=o == P.SUPPORTED and P.SUPPORTED not in ok,
                             invalid_error_possible=P.NS_INVALID not in ok, nde_error_possible=P.NDE not in ok,
                             valid=truth == "valid", false_invalid=truth == "valid" and o == P.NS_INVALID,
                             error=o not in ok, establishable=establishable, correct_definite=o in good,
                             definite=o in P.DEFINITE, nominal=nominal, gates=gates))
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


def score(entries: list[dict], reports: dict, pilot: dict, donors: dict | None = None) -> dict:
    """entries: panel.assign(key, dropped, pool sizes) (or a subset); reports: {card id: (verdict,
    cause)}; donors: {dataset id: [donor ids]} from the manifest (for the design effects)."""
    rows = card_rows(entries, reports, pilot)
    out = dict(per_condition={}, per_level={}, per_truth={}, criteria={}, secondary={})

    def summary(rs):
        return dict(n=len(rs), errors=rate(sum(r["error"] for r in rs), len(rs)),
                    false_supported=rate(sum(r["sup_error"] for r in rs), sum(r["sup_error_possible"] for r in rs)),
                    false_invalid=rate(sum(r["false_invalid"] for r in rs), sum(r["valid"] for r in rs)),
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
    # S1: false SUPPORTED at every key null condition's key variant; passes only if all pass
    s1 = {}
    for c in P.CONDITIONS:
        if not c.key:
            continue
        rs = [r for r in rows if r["condition"] == c.name and r["variant"] == c.key_variant]
        k = sum(r["sup_error"] for r in rs)
        kmax = oc.error_rule(len(rs), oc.E_SUPPORTED)[0] if rs else -1
        s1[c.name] = dict(rate(k, len(rs)), max_allowed=kmax, passed=bool(rs) and k <= kmax)
        members[f"S1:{c.name}"] = (rs, "sup_error")
    crit["S1"] = dict(conditions=s1, passed=bool(s1) and all(v["passed"] for v in s1.values()))
    # S2: false SUPPORTED on every card where SUPPORTED is an error
    rs = [r for r in rows if r["sup_error_possible"]]
    k = sum(r["sup_error"] for r in rs)
    kmax = oc.error_rule(len(rs), oc.E_SUPPORTED)[0] if rs else -1
    crit["S2"] = dict(rate(k, len(rs)), nominal=oc.E_SUPPORTED, max_allowed=kmax, passed=bool(rs) and k <= kmax)
    members["S2"] = (rs, "sup_error")
    # S3: correct definite outcomes on establishable cards
    rs = [r for r in rows if r["establishable"]]
    k = sum(r["correct_definite"] for r in rs)
    kmin = oc.decisiveness_rule(len(rs), oc.D_NOMINAL)[0] if rs else 1
    crit["S3"] = dict(rate(k, len(rs)), nominal=oc.D_NOMINAL, min_required=kmin, passed=bool(rs) and k >= kmin,
                      definite_any=rate(sum(r["definite"] for r in rs), len(rs)))
    members["S3"] = (rs, "correct_definite")
    # S4: outside the allowed outcomes, all cards, against the mean of the cards' nominal rates
    k = sum(r["error"] for r in rows)
    e_bar = float(np.mean([r["nominal"] for r in rows])) if rows else 0.0
    kmax = oc.error_rule(len(rows), e_bar)[0] if rows else -1
    crit["S4"] = dict(rate(k, len(rows)), nominal=e_bar, max_allowed=kmax, passed=bool(rows) and k <= kmax)
    members["S4"] = (rows, "error")
    # S5: false "metric invalid" (GATE 4/5) on cards whose metric is valid by the truth
    rs = [r for r in rows if r["valid"]]
    k = sum(r["false_invalid"] for r in rs)
    kmax = oc.error_rule(len(rs), oc.E_INVALID)[0] if rs else -1
    crit["S5"] = dict(rate(k, len(rs)), nominal=oc.E_INVALID, max_allowed=kmax, passed=bool(rs) and k <= kmax)
    members["S5"] = (rs, "false_invalid")
    out["passed"] = all(crit[s]["passed"] for s in ("S1", "S2", "S3", "S4", "S5"))
    out["n_cards"] = len(rows)
    out["joint_pass_probability_sound"] = oc.joint_pass_probability(rows)
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
    for s in ("S1", "S2", "S3", "S4", "S5"):
        print(s, "PASS" if res["criteria"][s]["passed"] else "FAIL")
    for lim in res["limitations"]:
        print("limitation:", lim)
    print("validation", "PASSES" if res["passed"] else "FAILS")


if __name__ == "__main__":
    main()
