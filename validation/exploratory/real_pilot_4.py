#!/usr/bin/env python3
"""The second diagnostics on the data of the pilot's second run
(validation/exploratory/real_pilot_diagnostics_2_plan.md). Exploratory, development data only.

- (v) the composition of cells on preprint section 4: per unit the shares of cell types and tissues;
  their difference between the groups per sex (the total variation distance and per type, exact
  permutations over the units); mi_3bin(Smad3, Col1a1), and norm_pearson beside it, within each cell
  type with at least 20 cells in at least 3 units of each group, and combined with equal weights; a
  unit's MI against its share of fibroblast-like cells;
- (ii′) specificity again, for 1a-males, 2a and 2b: the donor-level t of the claim among 1,000
  controls matched on mean and detection over the full set of genes (real_pilot_4_extract.py), in
  the claimed direction and by absolute value, with the Clopper-Pearson interval of each share.

    python validation/exploratory/real_pilot_4.py --data RUN2_DATA --full EXTRACT_OUT --run2 real_pilot.json --out OUT

It writes OUT/real_pilot_4.json and OUT/real_pilot_4.log. The engine (the tag v0.3.0-prereg's, as
installed by the workflow) supplies the metrics and the statistics only; run_autopsy does not run.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import platform
import sys
import time
import traceback
from functools import partial
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "prereg"))
sys.path.insert(0, str(HERE))
import frozen  # noqa: E402  (standard library only)

NUMERICS = frozen.pin_numerics()  # v1's numerical path, set before numpy loads

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import metric_autopsy as MA  # noqa: E402
from metric_autopsy import SimpleData, metrics, stats  # noqa: E402
from metric_autopsy.stats import clopper_pearson  # noqa: E402
from anchors import log_cp10k_mean  # noqa: E402  (R1's metric)
import real_pilot_4_controls as C  # noqa: E402

ALPHA = 0.05
MIN_CELLS, MIN_UNITS = 20, 3
FIBRO = ("fibroblast", "mesenchymal", "stromal")
DONOR_PERMS, DONOR_MAX_EXACT = 9999, 20000
REST = "__rest_of_counts__"
GROUPS1 = ("young", "old")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fmt(x, nd=4) -> str:
    return "nan" if x is None or not np.isfinite(x) else f"{x:+.{nd}g}"


def fp(x, nd=4) -> str:
    return "nan" if x is None or not np.isfinite(x) else f"{x:.{nd}g}"


# --------------------------------------------------------------------------- #
# measures on reduced matrices: the named genes' counts and the rest of each cell's counts
# --------------------------------------------------------------------------- #
def reduced(X, cols, tot, rows=None) -> np.ndarray:
    sub = X[:, cols] if rows is None else X[np.ix_(rows, cols)]
    t = tot if rows is None else tot[rows]
    rest = t - sub.sum(axis=1, dtype=np.int64)
    return np.column_stack([sub.astype(np.float64), rest.astype(np.float64)])


def as_data(M, names) -> SimpleData:
    return SimpleData(M, pd.DataFrame(index=np.arange(M.shape[0])), [*names, REST])


def measure(kind: str, names):
    names = list(names)
    if kind == "mi_3bin":
        return partial(metrics.mi_3bin, gene_a=names[0], gene_b=names[1])
    if kind == "norm_pearson":
        return partial(metrics.norm_pearson, gene_a=names[0], gene_b=names[1])
    return partial(log_cp10k_mean, cols=list(range(len(names))))


def value(kind, X, genes, names, tot, rows) -> float:
    return float(measure(kind, names)(as_data(reduced(X, [genes.index(g) for g in names], tot, rows), names)))


def run2_values(run2: dict | None, case: str, claim: str) -> dict | None:
    """The second run's per-unit values of a claim's analysis (b), keyed by unit, or None."""
    if run2 is None:
        return None
    return next((x.get("b", {}).get("values") for x in run2.get("cases", {}).get(case, {}).get("claims", [])
                 if x.get("id") == claim), None)


def relabelings(n: int, n_second: int) -> np.ndarray:
    """Every assignment of n units to the second group of size n_second (a boolean matrix)."""
    out = []
    for comb in itertools.combinations(range(n), n_second):
        r = np.zeros(n, bool)
        r[list(comb)] = True
        out.append(r)
    return np.array(out)


def ranks(x) -> np.ndarray:
    return pd.Series(np.asarray(x, float)).rank(method="average").to_numpy()


def corr_perm(x, y) -> dict:
    """Pearson r and Spearman rho of x and y with the exact permutation p over every ordering of y."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    perms = np.array(list(itertools.permutations(range(n))), dtype=np.int16)
    out = {}
    for name, (a, b) in (("pearson", (x, y)), ("spearman", (ranks(x), ranks(y)))):
        az, bz = a - a.mean(), b - b.mean()
        den = np.sqrt((az ** 2).sum() * (bz ** 2).sum())
        if den == 0:
            out[name] = dict(r=float("nan"), p=float("nan"))
            continue
        r_obs = float((az * bz).sum() / den)
        r_all = (bz[perms] * az).sum(axis=1) / den
        tol = 1e-12
        out[name] = dict(r=r_obs, p=float(np.mean(np.abs(r_all) >= abs(r_obs) - tol)))
    out["orderings"] = int(len(perms))
    return out


# --------------------------------------------------------------------------- #
# (v) the composition of cells
# --------------------------------------------------------------------------- #
def composition(X, genes, obs, run2_values, log) -> dict:
    tot = X.sum(axis=1, dtype=np.int64)
    pair = list(C.PAIR)
    unit = obs["mouse"].astype(str).to_numpy()
    sex = obs["sex"].astype(str).to_numpy()
    age = obs["age"].astype(str).to_numpy()
    ctype = obs["cell_type"].astype(str).to_numpy()
    tissue = obs["tissue"].astype(str).to_numpy()
    fib_labels = sorted({t for t in set(ctype) if any(w in t.lower() for w in FIBRO)})
    fib = np.isin(ctype, fib_labels)
    log("## (v) the composition of cells on §4", "",
        f"fibroblast-like labels (the Census's cell_type containing {', '.join(repr(w) for w in FIBRO)}): "
        + "; ".join(f"{t} ({int((ctype == t).sum())} cells)" for t in fib_labels))
    res = dict(fibroblast_like_labels={t: int((ctype == t).sum()) for t in fib_labels}, sexes={})
    for sx in sorted(set(sex)):
        m_sex = sex == sx
        units = sorted(set(unit[m_sex]))
        grp = {u: str(next(iter(set(age[m_sex & (unit == u)])))) for u in units}
        if any(g not in GROUPS1 for g in grp.values()):
            raise ValueError(f"{sx}: groups {sorted(set(grp.values()))}")
        old = np.array([grp[u] == "old" for u in units])
        R = relabelings(len(units), int(old.sum()))
        rec = dict(units=units, groups=grp)
        # per-unit values over all cells, and the shares
        rows_u = {u: np.where(m_sex & (unit == u))[0] for u in units}
        mi_all = {u: value("mi_3bin", X, genes, pair, tot, r) for u, r in rows_u.items()}
        pe_all = {u: value("norm_pearson", X, genes, pair, tot, r) for u, r in rows_u.items()}
        common = [u for u in units if run2_values and u in run2_values]
        diff_vs_run2 = max(abs(mi_all[u] - run2_values[u]) for u in common) if common else float("nan")
        delta_all = float(np.mean([mi_all[u] for u in units if grp[u] == "old"]) - np.mean([mi_all[u] for u in units if grp[u] == "young"]))
        rec.update(mi_all=mi_all, pearson_all=pe_all, delta_all=delta_all, mi_vs_run2_max_abs_diff=diff_vs_run2,
                   mi_vs_run2_units=len(common))
        for what, lab in (("cell_type", ctype), ("tissue", tissue)):
            tab = pd.crosstab(unit[m_sex], lab[m_sex], normalize="index").reindex(units).fillna(0.0)
            S = tab.to_numpy()
            d_obs = S[old].mean(0) - S[~old].mean(0)
            d_all = np.array([S[r].mean(0) - S[~r].mean(0) for r in R])
            tvd_obs = 0.5 * float(np.abs(d_obs).sum())
            tvd_all = 0.5 * np.abs(d_all).sum(axis=1)
            p_tvd = float(np.mean(tvd_all >= tvd_obs - 1e-12))
            p_each = (np.abs(d_all) >= np.abs(d_obs) - 1e-12).mean(axis=0)
            rows = sorted(zip(tab.columns, S[~old].mean(0), S[old].mean(0), d_obs, p_each), key=lambda t: -abs(t[3]))
            rec[what] = dict(tvd=tvd_obs, tvd_p=p_tvd, relabelings=int(len(R)),
                             by_level=[dict(level=str(a), young=float(b), old=float(c), diff=float(dd), p=float(pp))
                                       for a, b, c, dd, pp in rows],
                             per_unit={u: {str(k): float(v) for k, v in tab.loc[u].items() if v > 0} for u in units})
        # MI and Pearson within cell types
        types = sorted(set(ctype[m_sex]))
        qual = {}
        for t in types:
            n_cells = {u: int((m_sex & (unit == u) & (ctype == t)).sum()) for u in units}
            ok = [u for u in units if n_cells[u] >= MIN_CELLS]
            if sum(grp[u] == "young" for u in ok) >= MIN_UNITS and sum(grp[u] == "old" for u in ok) >= MIN_UNITS:
                qual[t] = ok
        within = {}
        for kind in ("mi_3bin", "norm_pearson"):
            per_type = {}
            for t, ok in qual.items():
                vals = {u: value(kind, X, genes, pair, tot, np.where(m_sex & (unit == u) & (ctype == t))[0]) for u in ok}
                y = np.array([vals[u] for u in ok])
                g = np.array([1 if grp[u] == "young" else 0 for u in ok])
                d = float(y[g == 0].mean() - y[g == 1].mean())
                pr = stats.permutation_test_nested(y, g, np.zeros(len(y), int), n_perm=DONOR_PERMS,
                                                   rng=np.random.default_rng(C.SEED), max_exact=DONOR_MAX_EXACT)
                per_type[t] = dict(values=vals, delta=d, p=float(pr.p), method=pr.method, n_relabelings=int(pr.n),
                                   n={"young": int((g == 1).sum()), "old": int((g == 0).sum())})
            # combined, equal weights; every relabeling of the sex's units
            idx = {u: i for i, u in enumerate(units)}
            obs_d = float(np.mean([v["delta"] for v in per_type.values()])) if per_type else float("nan")
            perm_d = []
            for r in R:
                ds = []
                for t, v in per_type.items():
                    ok = qual[t]
                    yo = [v["values"][u] for u in ok if r[idx[u]]]
                    yy = [v["values"][u] for u in ok if not r[idx[u]]]
                    if yo and yy:
                        ds.append(np.mean(yo) - np.mean(yy))
                perm_d.append(np.mean(ds) if ds else np.nan)
            perm_d = np.array(perm_d)
            p_comb = float(np.mean(np.abs(perm_d[np.isfinite(perm_d)]) >= abs(obs_d) - 1e-12)) if per_type else float("nan")
            within[kind] = dict(types=per_type, combined=dict(delta=obs_d, p=p_comb, types=len(per_type),
                                                               relabelings=int(len(R))))
        rec["within"] = within
        # the association with the share of fibroblast-like cells
        share = {u: float(fib[rows_u[u]].mean()) for u in units}
        rec["fibroblast_like_share"] = share
        rec["association"] = corr_perm([mi_all[u] for u in units], [share[u] for u in units])
        res["sexes"][sx] = rec
        log(*lines_composition(sx, rec))
    m = res["sexes"].get("male")
    if m:
        c = m["within"]["mi_3bin"]["combined"]
        not_detected = bool(np.isfinite(c["p"]) and c["p"] >= ALPHA)
        less_than_half = bool(np.isfinite(c["delta"]) and c["delta"] > 0.5 * m["delta_all"])
        res["outcome"] = dict(not_detected=not_detected, less_than_half=less_than_half,
                              accounts_for_the_drop=not_detected and less_than_half)
        res["prediction"] = dict(text="внутри типов падение MI у самцов исчезает", holds=not_detected and less_than_half)
        log("", f"(v) outcome in males: Delta within types {fmt(c['delta'])} against the overall {fmt(m['delta_all'])}; not "
            f"detected (p {fp(c['p'])} >= 0.05): {not_detected}; less than half the overall drop: {less_than_half} -> "
            f"composition {'accounts' if res['outcome']['accounts_for_the_drop'] else 'does not account'} for the drop",
            f"the prediction («{res['prediction']['text']}») {'holds' if res['prediction']['holds'] else 'does not hold'}")
    return res


def lines_composition(sx, rec) -> list:
    out = ["", f"### {sx}: units {len(rec['units'])} (young {sum(g == 'young' for g in rec['groups'].values())}, old "
           f"{sum(g == 'old' for g in rec['groups'].values())}); MI over all cells: old minus young {fmt(rec['delta_all'])} "
           f"(against the second run's values of {rec['mi_vs_run2_units']} units, max |diff| {rec['mi_vs_run2_max_abs_diff']:.1e})"]
    for what in ("cell_type", "tissue"):
        r = rec[what]
        out.append(f"- {what}: total variation distance between the groups' mean compositions {r['tvd']:.4f}, exact p "
                   f"{fp(r['tvd_p'])} ({r['relabelings']} relabelings); the largest differences (young -> old, p):")
        for row in r["by_level"][:12]:
            out.append(f"    {row['level']}: {row['young']:.4f} -> {row['old']:.4f} ({fmt(row['diff'])}, p {fp(row['p'])})")
    for kind in ("mi_3bin", "norm_pearson"):
        w = rec["within"][kind]
        c = w["combined"]
        out.append(f"- {kind} within cell types ({c['types']} types with >= {MIN_CELLS} cells in >= {MIN_UNITS} units of each "
                   f"group): combined Delta (old minus young, equal weights) {fmt(c['delta'])}, exact p {fp(c['p'])} "
                   f"({c['relabelings']} relabelings)")
        for t, v in sorted(w["types"].items(), key=lambda kv: kv[1]["delta"]):
            ys = [v["values"][u] for u in v["values"] if rec["groups"][u] == "young"]
            os_ = [v["values"][u] for u in v["values"] if rec["groups"][u] == "old"]
            out.append(f"    {t}: young {np.mean(ys):.4g} ({v['n']['young']} units) -> old {np.mean(os_):.4g} "
                       f"({v['n']['old']}), Delta {fmt(v['delta'])}, p {fp(v['p'])}")
    a = rec["association"]
    out.append(f"- MI over all cells against the share of fibroblast-like cells, across the {len(rec['units'])} units: "
               f"Pearson r {fmt(a['pearson']['r'], 3)} (p {fp(a['pearson']['p'])}), Spearman rho {fmt(a['spearman']['r'], 3)} "
               f"(p {fp(a['spearman']['p'])}), exact over {a['orderings']} orderings")
    for u in rec["units"]:
        out.append(f"    {u} ({rec['groups'][u]}): MI {rec['mi_all'][u]:.5f}, Pearson {rec['pearson_all'][u]:+.4f}, "
                   f"fibroblast-like {rec['fibroblast_like_share'][u]:.4f}")
    return out


# --------------------------------------------------------------------------- #
# (ii′) specificity on the donor-level t
# --------------------------------------------------------------------------- #
def t_of(values: dict, mice, second: set) -> tuple[float, float, float]:
    """Delta (second minus first), its SE and t, by stats.nested_t in one stratum."""
    y = np.array([values[m] for m in mice], float)
    g = np.array([0 if m in second else 1 for m in mice])
    te = stats.nested_t(y, g, np.zeros(len(y), int))
    d, se = -float(te.est), float(te.se)
    return d, se, (d / se if np.isfinite(se) and se > 0 else float("nan"))


def outcome3(k: int, n: int) -> tuple[str, float, float]:
    lo, hi = clopper_pearson(k, n)
    return ("specific" if hi < ALPHA else "not specific" if lo >= ALPHA else "undecided"), float(lo), float(hi)


def specificity(full: Path, sel: dict, run2: dict | None, log) -> dict:
    res = {}
    for case, claims in (("case1", ("1a-males",)), ("case2", ("2a", "2b"))):
        z = np.load(full / f"full_{case}.npz", allow_pickle=False)
        X = np.asarray(z["X"], dtype=np.int32)
        genes = [str(g) for g in z["genes"]]
        obs = pd.DataFrame({k[4:]: z[k] for k in z.files if k.startswith("obs_")})
        tot = np.asarray(z["tot"], dtype=np.int64) if "tot" in z.files else X.sum(axis=1, dtype=np.int64)
        mice_col = obs["mouse"].astype(str).to_numpy()
        mice = sorted(set(mice_col))
        rows = {m: np.where(mice_col == m)[0] for m in mice}
        for cid in claims:
            if cid == "1a-males":
                kind, names, first, second_g, direction = "mi_3bin", list(C.PAIR), "young", "old", "decrease"
                grp = {m: str(next(iter(set(obs["age"].astype(str).to_numpy()[rows[m]])))) for m in mice}
            else:
                kind = "log_cp10k_mean"
                names = ["Xist"] if cid == "2a" else list(C.Y_GENES)
                first, second_g, direction = "female", "male", ("decrease" if cid == "2a" else "increase")
                grp = {m: str(next(iter(set(obs["sex"].astype(str).to_numpy()[rows[m]])))) for m in mice}
            second = {m for m in mice if grp[m] == second_g}
            vals = {m: value(kind, X, genes, names, tot, rows[m]) for m in mice}
            d, se, t = t_of(vals, mice, second)
            ctl = sel["cases"][case]["selection"][cid]["controls"]
            tc = []
            for c in ctl:
                v = {m: value(kind, X, genes, c, tot, rows[m]) for m in mice}
                tc.append(t_of(v, mice, second))
            ts = np.array([x[2] for x in tc])
            ok = np.isfinite(ts)
            n = int(ok.sum())
            tol = 1e-12 * max(1.0, abs(t))
            far1 = (ts[ok] <= t + tol) if direction == "decrease" else (ts[ok] >= t - tol)
            far2 = np.abs(ts[ok]) >= abs(t) - tol
            k1, k2 = int(far1.sum()), int(far2.sum())
            o1, lo1, hi1 = outcome3(k1, n)
            o2, lo2, hi2 = outcome3(k2, n)
            idx = np.where(ok)[0]
            rec = dict(metric=kind, genes=names, direction=direction, first=first, second=second_g, values=vals,
                       delta=d, se=se, t=t, controls=len(ctl), ranked=n,
                       undefined=[ctl[i] for i in np.where(~ok)[0]],
                       claimed=dict(k=k1, p=(1 + k1) / (1 + n), share_ci=[lo1, hi1], outcome=o1,
                                    at_least_as_far=sorted([(float(ts[i]), ctl[i]) for i, f in zip(idx, far1) if f],
                                                           reverse=direction == "increase")),
                       absolute=dict(k=k2, p=(1 + k2) / (1 + n), share_ci=[lo2, hi2], outcome=o2,
                                     at_least_as_large=sorted([(float(ts[i]), ctl[i]) for i, f in zip(idx, far2) if f],
                                                              key=lambda x: -abs(x[0]))),
                       control_t_quantiles=dict(zip(("5%", "25%", "50%", "75%", "95%"),
                                                    np.quantile(ts[ok], [0.05, 0.25, 0.5, 0.75, 0.95]).tolist()))
                       if n else {},
                       matched=sel["cases"][case]["selection"][cid]["matched"])
            theirs = run2_values(run2, case, {"1a-males": "1a"}.get(cid, cid))
            if theirs is not None:
                common = sorted(set(theirs) & set(vals))
                rec["run2_values"] = {m: theirs[m] for m in common}
                rec["vs_run2_max_abs_diff"] = float(max(abs(theirs[m] - vals[m]) for m in common)) if common else float("nan")
            res[cid] = rec
            log(*lines_specificity(cid, rec))
    pred = {"2a specific in the claimed direction": res["2a"]["claimed"]["outcome"] == "specific",
            "2a specific by absolute value": res["2a"]["absolute"]["outcome"] == "specific",
            "2b specific in the claimed direction": res["2b"]["claimed"]["outcome"] == "specific",
            "2b specific by absolute value": res["2b"]["absolute"]["outcome"] == "specific",
            "1a-males not specific in the claimed direction": res["1a-males"]["claimed"]["outcome"] == "not specific"}
    res["prediction"] = dict(text="Xist и Y-гены — крайние в обе стороны, §4 по самцам — типичный", parts=pred,
                             holds=all(pred.values()))
    log("", "(ii′) the prediction («" + res["prediction"]["text"] + "»):",
        *[f"- {k}: {'yes' if v else 'no'}" for k, v in pred.items()],
        f"- the prediction {'holds' if res['prediction']['holds'] else 'does not hold'}")
    return res


def lines_specificity(cid, r) -> list:
    out = ["", f"### {cid}: {r['metric']}({', '.join(r['genes'])}), {r['first']} -> {r['second']}, {r['direction']}: Delta "
           f"{fmt(r['delta'])}, SE {fp(r['se'])}, t {fmt(r['t'])}"
           + (f"; the second run's per-mouse values: max |diff| {r['vs_run2_max_abs_diff']:.2e}" if "vs_run2_max_abs_diff" in r else "")]
    for key, label in (("claimed", "in the claimed direction"), ("absolute", "by absolute value")):
        c = r[key]
        lst = c.get("at_least_as_far", c.get("at_least_as_large", []))
        out.append(f"- {label}: {c['k']} of {r['ranked']} controls at least as far (p {fp(c['p'])}; the share's 95% CI "
                   f"{c['share_ci'][0]:.4f} to {c['share_ci'][1]:.4f}) -> {c['outcome']}")
        for tval, g in lst[:25]:
            out.append(f"    {fmt(tval)}  {' + '.join(g)}")
        if len(lst) > 25:
            out.append(f"    ... and {len(lst) - 25} more (in the JSON)")
    q = r["control_t_quantiles"]
    out.append("- the controls' t: " + ", ".join(f"{k} {fmt(v)}" for k, v in q.items())
               + f"; without a t (SE 0): {len(r['undefined'])}")
    out.append("- matched: " + "; ".join(f"{g}: mean {v['mean']:.4g} (controls {v['mean_range'][0]:.4g} to "
                                          f"{v['mean_range'][1]:.4g}), detection {v['detection']:.4f} ({v['detection_range'][0]:.4f} "
                                          f"to {v['detection_range'][1]:.4f})" for g, v in r["matched"].items()))
    return out


# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--data", type=Path, required=True, help="the second run's select.json and case files")
    p.add_argument("--full", type=Path, required=True, help="real_pilot_4_extract.py's output")
    p.add_argument("--run2", type=Path, help="the second run's real_pilot.json (its per-mouse values, compared)")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    lines = []

    def log(*msg):
        for m in msg:
            print(m, flush=True)
            lines.append(str(m))
    spec_bytes = (args.data / "select.json").read_bytes()
    spec = json.loads(spec_bytes)
    sel_bytes = (args.full / "selection.json").read_bytes()
    sel = json.loads(sel_bytes)
    run2 = json.loads(args.run2.read_text()) if args.run2 else None
    tag = "SYNTHETIC " if spec.get("synthetic") else ""
    out = dict(synthetic=spec.get("synthetic"), engine=dict(version=MA.__version__, path=str(Path(MA.__file__).resolve().parent)),
               python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, numerics=NUMERICS,
               select_sha256=hashlib.sha256(spec_bytes).hexdigest(), selection_sha256=hashlib.sha256(sel_bytes).hexdigest(),
               run2_sha256=sha256_file(args.run2) if args.run2 else None,
               full_files={f: sha256_file(args.full / f) for f in ("full_case1.npz", "full_case2.npz") if (args.full / f).exists()},
               settings=dict(alpha=ALPHA, min_cells=MIN_CELLS, min_units=MIN_UNITS, fibroblast_like=FIBRO, k=C.K, seed=C.SEED))
    log(f"# {tag}The second diagnostics on the second run's data (validation/exploratory/real_pilot_diagnostics_2_plan.md)",
        "", f"engine {MA.__version__} from {Path(MA.__file__).resolve().parent}; python {platform.python_version()}; numpy "
        f"{np.__version__}; pandas {pd.__version__}; numerics {NUMERICS}",
        f"select.json sha256 {out['select_sha256']}; the extraction's selection.json sha256 {out['selection_sha256']}; "
        + "; ".join(f"{k} {v}" for k, v in out["full_files"].items())
        + (f"; the second run's real_pilot.json sha256 {out['run2_sha256']}" if args.run2 else ""), "")
    t0 = time.time()

    def save():
        out["seconds"] = round(time.time() - t0, 1)
        (args.out / "real_pilot_4.json").write_text(json.dumps(out, indent=1, default=str))
        (args.out / "real_pilot_4.log").write_text("\n".join(lines) + "\n")
    try:
        case = spec["cases"]["case1"]
        path = args.data / case["data"]["file"]
        if sha256_file(path) != case["data"]["sha256"]:
            raise ValueError(f"{path.name}: sha256 differs from select.json's")
        z = np.load(path, allow_pickle=False)
        X = np.asarray(z["X"], dtype=np.int32)
        genes = [str(g) for g in z["genes"]]
        obs = pd.DataFrame({k[4:]: z[k] for k in z.files if k.startswith("obs_")})
        out["composition"] = composition(X, genes, obs, run2_values(run2, "case1", "1a"), log)
    except Exception as exc:  # reported as it is
        out["composition"] = dict(error=repr(exc), traceback=traceback.format_exc().splitlines())
        log(f"(v) EXCEPTION {exc!r}", *out["composition"]["traceback"][-8:])
    save()
    try:
        if "error" in sel:
            raise ValueError(f"the extraction failed: {sel['error']}")
        out["specificity"] = specificity(args.full, sel, run2, log)
    except Exception as exc:  # reported as it is
        out["specificity"] = dict(error=repr(exc), traceback=traceback.format_exc().splitlines())
        log(f"(ii′) EXCEPTION {exc!r}", *out["specificity"]["traceback"][-8:])
    log("", f"{tag}the second diagnostics: {time.time() - t0:.1f} s")
    save()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
