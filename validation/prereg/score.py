"""Scoring of the confirmatory validation (validation/prereg/v1.md, sections 5 and 6).

    python validation/prereg/score.py --key KEY --results RESULTS --pilot pilot.json --out scores.json

Run only after the blind run's results (reports, run log, manifest) are committed with their
sha256 (section 8). The key re-derives the condition, variant and gene pair of every dataset
ID (``panel.assign``); each claim card's verdict is reduced to its label and compared with the
design's allowed set (``panel.allowed``: for the real effects NO DETECTABLE EFFECT is allowed
only where |Δ*| < SESOI); the primary outcomes P1-P3 and the criteria S1-S4 follow, every rate
with its two-sided 95% Clopper-Pearson interval. A missing report or an engine error counts as a
verdict outside the allowed set. Secondary: every rate per expression level, and intervals that
allow for datasets sharing donors (``overlap_interval``).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import stats

import oc
import panel as P

LABELS = P.LABELS
NOMINAL_ERROR = 0.025       # false SUPPORTED with directional claims (alpha / 2); outside the set
NOMINAL_DECISIVENESS = 0.85  # correct definite verdicts on establishable cards


def label(verdict: str | None) -> str:
    """The verdict's label: the text before ' — ', without bracketed qualifiers.
    'SUPPORTED (provisional until replicated)' and 'SUPPORTED — replicated' are SUPPORTED."""
    if not verdict:
        return "ERROR"
    head = verdict.split(" — ")[0].split(" [")[0].split(" (")[0].strip()
    return head if head in LABELS else "ERROR"


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


def score(entries: list[dict], verdicts: dict, pilot: dict, donors: dict | None = None) -> dict:
    """entries: panel.assign(key) (or a subset); verdicts: {card id: verdict text or None};
    donors: {dataset id: [donor ids]} from the manifest (for the overlap intervals)."""
    conds = P.conditions()
    est = pilot.get("establishable", {})
    rows = []
    for e in entries:
        c = conds[e["condition"]]
        level = P.level_of_pair(e["pair"])
        allowed = P.allowed(c, e["variant"], e["pair"], pilot)
        good = P.definite(c, e["variant"], e["pair"], pilot)
        establishable = c.oracle and bool(est.get(f"{c.name}:{e['variant']}:{level}", {}).get("establishable"))
        for cid in P.card_ids(e):
            lab = label(verdicts.get(cid))
            rows.append(dict(id=cid, dataset=e["id"], condition=c.name, variant=e["variant"], level=level,
                             label=lab, null=c.null, key=c.key and e["variant"] == c.key_variant,
                             outside=lab not in allowed, establishable=establishable,
                             correct_definite=lab in good, supported=lab == "SUPPORTED"))
    out = dict(per_condition={}, per_level={}, criteria={}, secondary={})

    def summary(rs):
        return dict(n=len(rs), outside=rate(sum(r["outside"] for r in rs), len(rs)),
                    false_supported=(rate(sum(r["supported"] for r in rs), len(rs)) if rs and rs[0]["null"] else None),
                    correct_definite=(rate(sum(r["correct_definite"] for r in rs), len(rs))
                                      if rs and rs[0]["establishable"] else None),
                    labels={lab: sum(r["label"] == lab for r in rs) for lab in (*LABELS, "ERROR")})

    groups = {}
    for r in rows:
        groups.setdefault(f"{r['condition']}:{r['variant']}", []).append(r)
        groups.setdefault(f"{r['condition']}:{r['variant']}:{r['level']}", []).append(r)
    for key, rs in groups.items():
        (out["per_level"] if key.count(":") == 2 else out["per_condition"])[key] = summary(rs)
    nulls = [r for r in rows if r["null"]]
    estab = [r for r in rows if r["establishable"]]
    crit = out["criteria"]
    # S1: false SUPPORTED at every key null condition's key variant; passes only if all pass
    s1 = {}
    for c in P.CONDITIONS:
        if not c.key:
            continue
        rs = [r for r in rows if r["condition"] == c.name and r["variant"] == c.key_variant]
        k = sum(r["supported"] for r in rs)
        kmax = oc.error_rule(len(rs), NOMINAL_ERROR)[0] if rs else -1
        s1[c.name] = dict(rate(k, len(rs)), max_allowed=kmax, passed=bool(rs) and k <= kmax)
    crit["S1"] = dict(conditions=s1, passed=bool(s1) and all(v["passed"] for v in s1.values()))
    k = sum(r["supported"] for r in nulls)
    kmax = oc.error_rule(len(nulls), NOMINAL_ERROR)[0] if nulls else -1
    crit["S2"] = dict(rate(k, len(nulls)), max_allowed=kmax, passed=bool(nulls) and k <= kmax)
    k = sum(r["correct_definite"] for r in estab)
    kmin = oc.decisiveness_rule(len(estab), NOMINAL_DECISIVENESS)[0] if estab else 1
    crit["S3"] = dict(rate(k, len(estab)), min_required=kmin, passed=bool(estab) and k >= kmin,
                      definite_any=rate(sum(r["label"] not in ("INCONCLUSIVE", "ERROR") for r in estab), len(estab)))
    k = sum(r["outside"] for r in rows)
    kmax = oc.error_rule(len(rows), NOMINAL_ERROR)[0]
    crit["S4"] = dict(rate(k, len(rows)), max_allowed=kmax, passed=k <= kmax)
    out["passed"] = all(crit[s]["passed"] for s in ("S1", "S2", "S3", "S4"))
    out["n_cards"] = len(rows)
    # secondary: the primary rates per level, and with intervals that allow for shared donors
    for level in P.LEVELS:
        rl = [r for r in rows if r["level"] == level]
        out["secondary"][level] = dict(
            false_supported=rate(sum(r["supported"] for r in rl if r["null"]), sum(r["null"] for r in rl)),
            outside=rate(sum(r["outside"] for r in rl), len(rl)),
            correct_definite=rate(sum(r["correct_definite"] for r in rl if r["establishable"]),
                                  sum(r["establishable"] for r in rl)))
    if donors:
        def ov(rs, field):
            return overlap_interval([r[field] for r in rs], [donors.get(r["dataset"], []) for r in rs])
        out["secondary"]["shared_donors"] = dict(
            S1={c: ov([r for r in rows if r["condition"] == c and r["key"]], "supported") for c in s1},
            S2=ov(nulls, "supported"), S3=ov(estab, "correct_definite"), S4=ov(rows, "outside"))
    return out


def read_results(results: Path) -> tuple[dict, dict]:
    """Verdicts by card id and donors by dataset id, from the blind run's merged results."""
    manifest = json.loads((results / "manifest.json").read_text())
    verdicts, donors = {}, {}
    for d in manifest["datasets"]:
        donors[d["id"]] = d.get("donors", [])
        for c in d["cards"]:
            rep_path = results / "reports" / f"{c['id']}.json"
            rep = json.loads(rep_path.read_text()) if rep_path.exists() else {}
            verdicts[c["id"]] = rep.get("verdict") if "error" not in rep else None
    return verdicts, donors


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--key", required=True, help="the revealed key (32 hex characters)")
    p.add_argument("--results", required=True)
    p.add_argument("--pilot", required=True)
    p.add_argument("--out", default="scores.json")
    args = p.parse_args(argv)
    pilot = json.loads(Path(args.pilot).read_text())
    verdicts, donors = read_results(Path(args.results))
    res = score(P.assign(args.key, pilot.get("dropped", ())), verdicts, pilot, donors)
    res["key_sha256"] = P.key_commitment(args.key)
    Path(args.out).write_text(json.dumps(res, indent=1))
    for s in ("S1", "S2", "S3", "S4"):
        print(s, "PASS" if res["criteria"][s]["passed"] else "FAIL")
    print("validation", "PASSES" if res["passed"] else "FAILS")


if __name__ == "__main__":
    main()
