"""Scoring of the confirmatory validation (validation/prereg/v1.md, sections 5 and 6).

    python validation/prereg/score.py --key-seed SEED --panel DIR --reports REPORTS --pilot pilot.json

Run only after the engine's reports and run log are committed (section 8, step 4). The key seed
re-derives the condition of every dataset ID (``panel.assign``); each claim card's verdict is
reduced to its label and compared with the condition's allowed set; the primary outcomes
P1-P3 and the criteria S1-S4 follow, every rate with its two-sided 95% Clopper-Pearson
interval. A missing report or an engine error counts as a verdict outside the allowed set.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scipy import stats

import oc
import panel as P

LABELS = ("SUPPORTED", "NOT SUPPORTED", "NO DETECTABLE EFFECT", "INCONCLUSIVE", "UNIDENTIFIABLE",
          "DEGENERATE METRIC")


def label(verdict: str | None) -> str:
    """The verdict's label: the text before ' — ', without bracketed qualifiers.
    'SUPPORTED (provisional until replicated)' and 'SUPPORTED — replicated' are SUPPORTED."""
    if not verdict:
        return "ERROR"
    head = verdict.split(" — ")[0].split(" [")[0].split(" (")[0].strip()
    return head if head in LABELS else "ERROR"


def cp(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    lo = 0.0 if k == 0 else float(stats.beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(stats.beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi


def rate(k: int, n: int) -> dict:
    lo, hi = cp(k, n) if n else (float("nan"), float("nan"))
    return dict(k=k, n=n, rate=k / n if n else float("nan"), ci95=[lo, hi])


def score(entries: list[dict], verdicts: dict, pilot: dict) -> dict:
    """entries: panel.assign(key) (or a subset); verdicts: {card id: verdict text or None}."""
    conds = P.conditions()
    est = pilot.get("establishable", {})
    rows = []
    for e in entries:
        c = conds[e["condition"]]
        ids = [e["id"]] if c.cards == 1 else [f"{e['id']}a", f"{e['id']}b"]
        for cid in ids:
            lab = label(verdicts.get(cid))
            key = f"{c.name}:{e['variant']}"
            rows.append(dict(id=cid, condition=c.name, variant=e["variant"], label=lab,
                             null=c.null, key=c.key and e["variant"] == c.key_variant,
                             outside=lab not in c.allowed,
                             establishable=bool(est.get(key, {}).get("establishable")) and c.definite is not None,
                             definite=c.definite))
    out = dict(per_condition={}, criteria={})
    groups = {}
    for r in rows:
        groups.setdefault(f"{r['condition']}:{r['variant']}", []).append(r)
    for key, rs in groups.items():
        out["per_condition"][key] = dict(
            n=len(rs), outside=rate(sum(r["outside"] for r in rs), len(rs)),
            false_supported=rate(sum(r["label"] == "SUPPORTED" for r in rs), len(rs)) if rs[0]["null"] else None,
            definite=(rate(sum(r["label"] != "INCONCLUSIVE" and r["label"] != "ERROR" for r in rs), len(rs))
                      if rs[0]["establishable"] else None),
            correct_definite=(rate(sum(r["label"] == r["definite"] for r in rs), len(rs))
                              if rs[0]["establishable"] else None),
            labels={lab: sum(r["label"] == lab for r in rs) for lab in (*LABELS, "ERROR")})
    nulls = [r for r in rows if r["null"]]
    estab = [r for r in rows if r["establishable"]]
    crit = out["criteria"]
    # S1: false SUPPORTED at each key null condition's key variant
    s1 = {}
    for c in P.CONDITIONS:
        if not c.key:
            continue
        rs = [r for r in rows if r["condition"] == c.name and r["variant"] == c.key_variant]
        k = sum(r["label"] == "SUPPORTED" for r in rs)
        kmax = oc.error_rule(len(rs), 0.025)[0]
        s1[c.name] = dict(rate(k, len(rs)), max_allowed=kmax, passed=k <= kmax)
    crit["S1"] = dict(conditions=s1, passed=all(v["passed"] for v in s1.values()))
    k = sum(r["label"] == "SUPPORTED" for r in nulls)
    kmax = oc.error_rule(len(nulls), 0.025)[0] if nulls else -1
    crit["S2"] = dict(rate(k, len(nulls)), max_allowed=kmax, passed=bool(nulls) and k <= kmax)
    k = sum(r["label"] not in ("INCONCLUSIVE", "ERROR") for r in estab)
    kmin = oc.decisiveness_rule(len(estab), 0.85)[0] if estab else 1
    crit["S3"] = dict(rate(k, len(estab)), min_required=kmin, passed=bool(estab) and k >= kmin,
                      correct=rate(sum(r["label"] == r["definite"] for r in estab), len(estab)))
    k = sum(r["outside"] for r in rows)
    kmax = oc.error_rule(len(rows), 0.025)[0]
    crit["S4"] = dict(rate(k, len(rows)), max_allowed=kmax, passed=k <= kmax)
    out["passed"] = all(crit[s]["passed"] for s in ("S1", "S2", "S3", "S4"))
    out["n_cards"] = len(rows)
    return out


def read_verdicts(reports: Path) -> dict:
    out = {}
    for f in sorted(reports.glob("*.json")):
        rep = json.loads(f.read_text())
        out[f.stem] = rep.get("verdict") if "error" not in rep else None
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--key-seed", type=int, required=True)
    p.add_argument("--reports", required=True)
    p.add_argument("--pilot", required=True)
    p.add_argument("--out", default="scores.json")
    args = p.parse_args(argv)
    pilot = json.loads(Path(args.pilot).read_text())
    res = score(P.assign(args.key_seed), read_verdicts(Path(args.reports)), pilot)
    Path(args.out).write_text(json.dumps(res, indent=1))
    for s in ("S1", "S2", "S3", "S4"):
        print(s, "PASS" if res["criteria"][s]["passed"] else "FAIL")
    print("validation", "PASSES" if res["passed"] else "FAILS")


if __name__ == "__main__":
    main()
