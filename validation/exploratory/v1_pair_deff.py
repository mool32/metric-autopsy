"""Two more numbers for validation v1, from its published reports (exploratory; planned in
v1_analyzers_plan.md before the run): v1's sensitivity on the real effects, and the design effect of
cards that share a gene pair, next to v1's design effect of datasets that share donors.

1. Sensitivity: the share SUPPORTED among the cards of a real effect with a valid metric and
   |Δ*| >= SESOI (the cells of S6's effect stratum), and among all with a valid metric; per
   condition and variant; with every other outcome of those cards.
2. Design effects for every cell of v1's criteria (S1 per key null, the cells of S2, S4, S5 and S6,
   S3's strata and conditions, S7a, S7b), each cell's cards and outcome as v1's scoring forms them
   (score.score):
   - by donor: v1's estimator (score.overlap_interval: Cov(y_j, y_k) = beta x the donors two
     datasets share); it reproduces scores.json's "shared_donors";
   - by gene pair: the same estimator with the gene pair (background, pool index) as the one unit
     two cards share or not (a one-way cluster design effect);
   - both: Cov = beta_donor x (shared donors) + beta_pair x (same pair), fitted by least squares on
     the residuals' cross-products over all pairs of cards, each beta truncated at 0 (the other then
     refitted alone).
   The useless metrics (N6a-c) do not read their pair; their pair is only a label.
3. With --cards (the per-card results of v1_analyzers.py): the same three design effects for the
   main rates of the four analyzers.

    python validation/exploratory/v1_pair_deff.py --results RESULTS --pilot pilot.json \
        --scores scores.json [--cards cards.jsonl] > v1_pair_deff.log
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation" / "prereg"))
sys.path.insert(0, str(ROOT / "validation" / "exploratory"))
import numpy as np  # noqa: E402

import oc  # noqa: E402
import panel as P  # noqa: E402
import score  # noqa: E402
from v1_by_level import fmt  # noqa: E402


def pair_unit(r: dict) -> str:
    return f"{P.conditions()[r['condition']].background}:{r['pair']}"


def _units(sets: list) -> np.ndarray:
    """The cards x units incidence matrix (float64) of each card's set of units."""
    names = sorted({u for st in sets for u in st})
    at = {u: i for i, u in enumerate(names)}
    M = np.zeros((len(sets), len(names)))
    for j, st in enumerate(sets):
        M[j, [at[u] for u in st]] = 1.0
    return M


def _moments(M: np.ndarray, e: np.ndarray) -> tuple[float, float, float]:
    """For O = M M^T without its diagonal: (sum of O, sum of O^2, e' O e), through the small
    matrix M'M (units x units), so that no cards x cards matrix is formed."""
    S = M.sum(axis=1)
    G = M.T @ M
    tot = M.sum(axis=0)
    return (float(tot @ tot - S.sum()), float((G * G).sum() - (S * S).sum()),
            float((M.T @ e) @ (M.T @ e) - (e * e * S).sum()))


def overlap_deff(y, sets: list) -> dict:
    """score.overlap_interval's design effect, Cov(y_j, y_k) = beta x |S_j & S_k|, computed
    through M'M (the same estimator; checked against it on small cells in main)."""
    y = np.asarray(y, float)
    n = len(y)
    p = float(y.mean())
    vb = p * (1 - p) / n
    if vb == 0:
        return dict(deff=1.0, beta=0.0)
    e = y - p
    s_o, s_o2, eoe = _moments(_units(sets), e)
    beta = max(0.0, eoe / s_o2) if s_o2 > 0 else 0.0
    var = (float(e @ e) + beta * s_o) / n ** 2
    return dict(deff=max(1.0, var / vb), beta=beta)


def joint_deff(y, donor_sets: list, units: list) -> dict:
    """Cov(y_j, y_k) = b_d |S_j & S_k| + b_p [unit_j == unit_k] for j != k, least squares on the
    residuals' cross-products, each b >= 0; the design effect of the rate."""
    y = np.asarray(y, float)
    n = len(y)
    p = float(y.mean())
    vb = p * (1 - p) / n
    if vb == 0:
        return dict(deff=1.0, beta_donor=0.0, beta_pair=0.0)
    e = y - p
    Md, Mp = _units(donor_sets), _units([[u] for u in units])
    sd, a11, r1 = _moments(Md, e)
    sp, a22, r2 = _moments(Mp, e)
    a12 = 0.0  # sum over j != k of the shared donors of two cards of one pair
    S = Md.sum(axis=1)
    for c in range(Mp.shape[1]):
        rows = Mp[:, c] > 0
        t = Md[rows].sum(axis=0)
        a12 += float(t @ t - S[rows].sum())
    det = a11 * a22 - a12 * a12
    bd, bp = ((r1 * a22 - r2 * a12) / det, (r2 * a11 - r1 * a12) / det) if det > 0 else (0.0, 0.0)
    if bd < 0 or bp < 0:  # truncate the negative one at 0 and refit the other alone
        cand = [(max(0.0, r1 / a11) if a11 > 0 else 0.0, 0.0), (0.0, max(0.0, r2 / a22) if a22 > 0 else 0.0)]

        def sse(b):
            return -2 * (b[0] * r1 + b[1] * r2) + b[0] ** 2 * a11 + b[1] ** 2 * a22 + 2 * b[0] * b[1] * a12
        bd, bp = min(cand, key=sse)
    var = (float(e @ e) + bd * sd + bp * sp) / n ** 2
    return dict(deff=max(1.0, var / vb), beta_donor=float(bd), beta_pair=float(bp))


def deffs(rs: list, flag, donors: dict) -> tuple[dict, dict, dict]:
    y = [float(flag(r)) for r in rs]
    dsets = [donors.get(r["dataset"], []) for r in rs]
    units = [pair_unit(r) for r in rs]
    return overlap_deff(y, dsets), overlap_deff(y, [[u] for u in units]), joint_deff(y, dsets, units)


def members_of(rows: list, pilot: dict) -> dict:
    """Every criterion cell's cards and outcome field, as score.score forms them."""
    fixed = (pilot.get("s3_rules") or {}).get("tiers")
    rules = oc.criteria_rules(rows, sims=200, s3_tiers=fixed)  # the cells' members do not depend on sims
    members = {}
    for c in P.CONDITIONS:
        if c.key:
            members[f"S1:{c.name}"] = ([r for r in rows if r["condition"] == c.name and r["variant"] == c.key_variant],
                                       "sup_error")
    for cell, rule in rules["cells"].items():
        if cell.split(":")[0] in oc.ERROR_CRITERIA:
            members[cell] = ([rows[i] for i in rule["members"]], score.FLAG[rule["field"]])
    for st in oc.STRATA:
        members[f"S3:{st}"] = ([r for r in rows if r["establishable"] and r["stratum"] == st], "correct_definite")
    for name, rule in rules["S3c"].items():
        members[f"S3c:{name}"] = ([rows[i] for i in rule["members"]], "correct_definite")
    members["S7a"] = (rows, "rule_violation")
    members["S7b"] = (rows, "crashed")
    return members


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--pilot", type=Path, required=True)
    ap.add_argument("--scores", type=Path, required=True)
    ap.add_argument("--cards", type=Path)
    args = ap.parse_args(argv)
    manifest = json.loads((args.results / "manifest.json").read_text())
    pilot_raw, scores_raw = args.pilot.read_bytes(), args.scores.read_bytes()
    pilot, scores = json.loads(pilot_raw), json.loads(scores_raw)
    key = P.check_key(manifest["key"])
    reports, donors = score.read_results(args.results, key)
    entries = P.assign(key, pilot.get("dropped", ()), pilot.get("pool_size"))
    rows = score.card_rows(entries, reports, pilot)
    conds = P.conditions()
    out = ["# Validation v1: its sensitivity, and the design effect of cards sharing a gene pair (exploratory;",
           "# planned before the run in v1_analyzers_plan.md; v1 unchanged)",
           f"# key {key}; pilot.json sha256 {hashlib.sha256(pilot_raw).hexdigest()}; "
           f"scores.json sha256 {hashlib.sha256(scores_raw).hexdigest()}; {len(rows)} cards", ""]

    # 1. sensitivity
    eff = [r for r in rows if conds[r["condition"]].data == "effect"]
    for r in eff:
        r["delta_star"] = float(P.delta_of(r["condition"], r["variant"], r["pair"], pilot)["value"])
        r["big"] = abs(r["delta_star"]) >= P.sesoi_of(conds[r["condition"]], r["pair"], pilot)
    main_cell = [r for r in eff if r["truth"] == "valid" and r["big"]]
    valid = [r for r in eff if r["truth"] == "valid"]
    out += ["## 1. Sensitivity: SUPPORTED on the real effects",
            f"  a valid metric and |Δ*| >= SESOI: SUPPORTED {fmt(sum(r['outcome'] == P.SUPPORTED for r in main_cell), len(main_cell))}",
            f"  a valid metric, any Δ*:          SUPPORTED {fmt(sum(r['outcome'] == P.SUPPORTED for r in valid), len(valid))}",
            "  the outcomes of the main cell: " + "; ".join(
                f"{o} {sum(r['outcome'] == o for r in main_cell)}" for o in P.OUTCOMES
                if any(r["outcome"] == o for r in main_cell)),
            "  by condition and variant (main cell):"]
    for c in P.CONDITIONS:
        if c.data != "effect":
            continue
        for v, _ in c.variants:
            rs = [r for r in main_cell if r["condition"] == c.name and r["variant"] == v]
            if rs:
                out.append(f"    {c.name} {v:<10} SUPPORTED {fmt(sum(r['outcome'] == P.SUPPORTED for r in rs), len(rs))}")
    by_d, by_p, both = deffs(main_cell, lambda r: r["outcome"] == P.SUPPORTED, donors)
    out += [f"  design effect of the main cell: donors {by_d['deff']:.2f}, pair {by_p['deff']:.2f}, both {both['deff']:.2f}",
            f"  its {len(main_cell)} cards share {len({pair_unit(r) for r in main_cell})} gene pairs", ""]

    # the low-memory estimator against v1's own (score.overlap_interval) on the cells of up to 1,500 cards
    checks = []
    for name, (rs, field) in members_of(rows, pilot).items():
        if 0 < len(rs) <= 1500:
            y = [float(r[field]) for r in rs]
            for sets in ([donors.get(r["dataset"], []) for r in rs], [[pair_unit(r)] for r in rs]):
                a, b = overlap_deff(y, sets)["deff"], score.overlap_interval(y, sets)["deff"]
                checks.append(abs(a - b) / b)
    out += [f"(the design effects below through M'M: on {len(checks)} cell and unit choices of up to 1,500 cards they",
            f"equal score.overlap_interval's within a relative {max(checks):.1e})", ""]

    # 2. design effects of every criterion cell
    out += ["## 2. Design effects of v1's criterion cells: by donor (v1's), by gene pair, and both",
            "  cell                        n     k | by donor (scores.json) | by pair | both (b_donor, b_pair) | pairs"]
    published = (scores.get("secondary") or {}).get("shared_donors") or {}
    worst = []
    for name, (rs, field) in members_of(rows, pilot).items():
        if not rs:
            continue
        by_d, by_p, both = deffs(rs, lambda r: r[field], donors)
        pub = published.get(name, {}).get("deff")
        k = sum(bool(r[field]) for r in rs)
        out.append(f"  {name:<26} {len(rs):>5} {k:>5} | {by_d['deff']:>6.2f} ({pub if pub is None else f'{pub:.2f}'}) "
                   f"| {by_p['deff']:>6.2f} | {both['deff']:>6.2f} ({both['beta_donor']:.2e}, {both['beta_pair']:.2e}) "
                   f"| {len({pair_unit(r) for r in rs})}")
        worst.append((by_p["deff"], name))
    worst.sort(reverse=True)
    out += ["  the largest design effects by pair: " + ", ".join(f"{n} {d:.2f}" for d, n in worst[:8])]
    cells = members_of(rows, pilot)
    out.append("  where the counted outcomes fall, in the cells whose design effect by pair exceeds 2 (k/n per pair):")
    for d, name in worst:
        if d <= 2:
            break
        rs, field = cells[name]
        per = {}
        for r in rs:
            per.setdefault(pair_unit(r), []).append(bool(r[field]))
        if name.startswith("S3"):  # a success counted: the pairs where it is rarest
            top = sorted(per.items(), key=lambda kv: (sum(kv[1]) / len(kv[1]), kv[0]))[:3]
            which = "the three lowest rates"
        else:
            top = sorted(per.items(), key=lambda kv: (-sum(kv[1]), kv[0]))[:3]
            which = "the three with the most"
        k = sum(bool(r[field]) for r in rs)
        out.append(f"    {name:<18} k {k} of {len(rs)} on {len(per)} pairs; {which}: " + ", ".join(
            f"{u} {sum(v)}/{len(v)}" for u, v in top))
    out.append("")

    # S1's operating characteristics at these design effects (as scores.json's s1_at_realized_deff
    # does for the donors' one: oc.error_rule_deff, a beta-binomial with deff times the variance)
    out += ["## 2b. S1 at the realized design effects: P(the cell passes) at the nominal 2.5% and at a doubled 5%",
            "  (the principle wants P(pass | doubled) <= 0.05; given the background the count is binomial, deff 1)"]
    for name, (rs, field) in members_of(rows, pilot).items():
        if not name.startswith("S1:"):
            continue
        by_d, by_p, both = deffs(rs, lambda r: r[field], donors)
        cells = []
        for label, d in (("deff 1", 1.0), ("donor", by_d["deff"]), ("pair", by_p["deff"]), ("both", both["deff"])):
            d = 1.0 if d < 1 + 1e-6 else d  # oc.error_rule_deff's beta-binomial overflows at a deff just above 1
            _, ps, pd = oc.error_rule_deff(len(rs), oc.E_SUPPORTED, d)
            cells.append(f"{label} {d:.2f}: {ps:.3f} / {pd:.3f}")
        out.append(f"  {name:<8} " + " | ".join(cells))
    out.append("")

    # 3. the analyzers' main rates
    if args.cards:
        mine = {}
        for line in args.cards.read_text().splitlines():
            if line.strip():
                m = json.loads(line)
                mine[m["card"]] = m
        byid = {r["id"]: r for r in rows}
        joined = [dict(byid[cid], **m) for cid, m in mine.items() if m.get("data_ok")]

        def call(r, an):
            if an in ("a", "b"):
                return r[f"{an}_call"]
            o = r["c_outcome"] if an == "c" else r["outcome"]
            return "+" if o == P.SUPPORTED else "-" if o == P.NS_OPPOSITE else "0"

        nulls = [r for r in joined if conds[r["condition"]].data == "null" and conds[r["condition"]].metric == "norm_pearson"]
        risk = [r for r in nulls if r["truth"] in ("valid", "ambiguous")]
        effs = [r for r in joined if conds[r["condition"]].data == "effect"]
        for r in effs:
            r["big"] = abs(float(P.delta_of(r["condition"], r["variant"], r["pair"], pilot)["value"])) >= \
                P.sesoi_of(conds[r["condition"]], r["pair"], pilot)
        sens = [r for r in effs if r["truth"] == "valid" and r["big"]]
        out += ["## 3. Design effects of the four analyzers' main rates (by donor, by pair, both)"]
        for label, rs in (("false 'effect', all nulls", nulls), ("false 'effect', nulls at risk", risk),
                          ("'effect', real, valid, |Δ*| >= SESOI", sens)):
            for an in ("a", "b", "c", "d"):
                y = lambda r: call(r, an) == "+"  # noqa: E731
                by_d, by_p, both = deffs(rs, y, donors)
                out.append(f"  {label:<38} {an}: {fmt(sum(y(r) for r in rs), len(rs))}; deff donor {by_d['deff']:.2f}, "
                           f"pair {by_p['deff']:.2f}, both {both['deff']:.2f}")
        out.append("")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
