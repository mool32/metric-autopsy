#!/usr/bin/env python3
"""The analyses of the third diagnostics and of v0.4's acceptance tests
(validation/exploratory/real_pilot_diagnostics_3_plan.md). Exploratory, for development data only.

Parts (the workflow runs them as parallel jobs):

- ``b``: §4's technical counts in the males of case 1: (b1) the detection and level of the plan's
  fibroblast genes per (mouse, type), (b2) the index of foreign counts, combined over the
  non-fibroblast types with the exact relabelings of the 9 units, (b3) its association with MI, and
  (b4) the metadata (the Census's obs, the source h5ad's obs, GEO's records);
- ``a``: the rule of the domain on §4 in males and on Marsh et al. 2022's three contrasts;
- ``case2``, ``case1m``, ``case1``, ``marsh``, ``case3a``, ``case3b``, ``case3c``: the acceptance tests,
  ``run_autopsy`` of the engine installed (v0.4.0.dev0 of claude/v0.4) on each claim;
- ``--merge``: the parts' outputs into real_pilot_5.json and real_pilot_5.log, with every criterion and
  prediction scored as the plan fixes it.

    python validation/exploratory/real_pilot_5.py --inputs INPUTS --full FULL --part PART --out OUT
    python validation/exploratory/real_pilot_5.py --merge PARTS_DIR --out OUT [--synthetic]

INPUTS holds data/ (the case files, select.json, select1.json) and geo/ (GEO's SOFT records); FULL the
extraction's outputs (real_pilot_5_extract.py). Every exception is caught per analysis and reported
with its traceback; a part that raises still writes what it has.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
import platform
import re
import sys
import time
import traceback
from functools import partial
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "prereg"))
import frozen  # noqa: E402  (standard library only)

NUMERICS = frozen.pin_numerics()  # v1's numerical path, set before numpy loads

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import metric_autopsy as MA  # noqa: E402  (the engine of claude/v0.4, installed by the workflow)
from metric_autopsy import SimpleData, injected_signal, metrics, run_autopsy  # noqa: E402
from metric_autopsy.domain import DOMAIN_DETECTION, claim_domain  # noqa: E402
from anchors import log_cp10k_mean  # noqa: E402  (R1's metric, the same function)
import real_pilot_5_spec as S  # noqa: E402

N_GENES = 2000  # D1's cut


# --------------------------------------------------------------------------- #
# loading
# --------------------------------------------------------------------------- #
def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_npz(path: Path):
    """(X, genes, obs): X dense as stored, or CSR for a sparse file (panel.save_npz)."""
    z = np.load(path, allow_pickle=False)
    obs = pd.DataFrame({k[4:]: z[k] for k in z.files if k.startswith("obs_")})
    if "X_data" in z.files:
        from scipy import sparse
        X = sparse.csr_matrix((z["X_data"], z["X_indices"], z["X_indptr"]), shape=tuple(z["X_shape"]))
    else:
        X = z["X"]
    return X, [str(g) for g in z["genes"]], obs


class Pool:
    """The check of specificity's pool: the same cells as the claim's data, every gene (sparse)."""

    def __init__(self, X, var_names, obs):
        self.X, self.var_names, self.obs = X, list(var_names), obs


def load_case(inp: Path, full: Path, name: str, synthetic: bool, with_pool: bool = True, rows=None):
    """The case file (D1's 2,000 genes, as SimpleData) and its full set of genes (a Pool), the same
    cells in the same order; with `rows` (a boolean mask on the obs), only those cells."""
    path = inp / "data" / f"{name}.npz"
    sha = sha256_file(path)
    if not synthetic and sha != S.CASE_SHA256[path.name]:
        raise ValueError(f"{path.name}: sha256 {sha}, the plan has {S.CASE_SHA256[path.name]}")
    X, genes, obs = load_npz(path)
    keep = np.ones(len(obs), bool) if rows is None else np.asarray(rows(obs), bool)
    X, obs = X[keep], obs[keep].astype(str).reset_index(drop=True)
    data = SimpleData(np.asarray(X, dtype=np.float64), obs, genes)
    del X
    if not with_pool:
        return data, sha
    XF, genesF, obsF = load_npz(full / f"full_{name}.npz")
    if not np.array_equal(obsF["soma_joinid"].to_numpy(np.int64)[keep], obs["soma_joinid"].to_numpy(np.int64)):
        raise ValueError(f"full_{name}.npz: its cells are not {path.name}'s, in order")
    XF = XF.tocsr()[np.where(keep)[0]] if rows is not None else XF.tocsr()
    return data, Pool(XF, genesF, obs), sha


def d1_cut(X, genes, named) -> tuple[np.ndarray, list]:
    """D1's cut (anchors._load): the named genes, then the most expressed by mean count (ties by column),
    N_GENES in all, in column order."""
    n = X.shape[0]
    mean = np.asarray(X.sum(axis=0, dtype=np.int64)).ravel() / max(n, 1)
    first = sorted(j for j, g in enumerate(genes) if g in set(named))
    order = np.argsort(-mean, kind="mergesort")
    taken = set(first)
    keep = sorted(list(dict.fromkeys([*first, *[int(j) for j in order if j not in taken]]))[:max(N_GENES, len(first))])
    sub = X[:, keep]
    return (sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)).astype(np.float64), [genes[j] for j in keep]


def load_marsh(full: Path):
    """Marsh et al.'s cells: D1's cut as SimpleData, and the full matrix as the pool."""
    X, genes, obs = load_npz(full / "marsh.npz")
    obs = obs.astype(str).reset_index(drop=True)
    dense, kept = d1_cut(X, genes, S.MARSH_PAIR)
    return SimpleData(dense, obs, kept), Pool(X.tocsr(), genes, obs)


def subset(data: SimpleData, pool: Pool, mask):
    mask = np.asarray(mask, bool)
    obs = data.obs[mask].reset_index(drop=True)
    return SimpleData(data.X[mask], obs, data.var_names), Pool(pool.X[np.where(mask)[0]], pool.var_names, obs)


# --------------------------------------------------------------------------- #
# small statistics
# --------------------------------------------------------------------------- #
def pearson(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def combined_delta(V: np.ndarray, young: np.ndarray) -> tuple[float, np.ndarray]:
    """V: types x units (NaN where a unit has fewer than the minimum cells of the type). Per type the young
    units' mean minus the old units' (each unit weighted equally); the mean over the types where both
    groups have a unit, equal weights."""
    with np.errstate(invalid="ignore"):
        fy = np.isfinite(V) & young[None, :]
        fo = np.isfinite(V) & ~young[None, :]
        my = np.where(fy, V, 0.0).sum(axis=1) / np.maximum(fy.sum(axis=1), 1)
        mo = np.where(fo, V, 0.0).sum(axis=1) / np.maximum(fo.sum(axis=1), 1)
    ok = (fy.sum(axis=1) > 0) & (fo.sum(axis=1) > 0)
    d = np.where(ok, my - mo, np.nan)
    return (float(np.nanmean(d)) if ok.any() else float("nan")), d


def exact_combined(V: np.ndarray, young: np.ndarray) -> dict:
    """The combined difference and its exact two-sided p over every relabeling of the units into the
    groups' sizes (each type's difference recomputed from its own units)."""
    obs, per = combined_delta(V, young)
    n, k = V.shape[1], int(young.sum())
    stats = []
    for pick in itertools.combinations(range(n), k):
        lab = np.zeros(n, bool)
        lab[list(pick)] = True
        stats.append(combined_delta(V, lab)[0])
    stats = np.asarray(stats)
    ok = np.isfinite(stats)
    p = float(np.mean(np.abs(stats[ok]) >= abs(obs) - 1e-12)) if ok.any() and np.isfinite(obs) else float("nan")
    return dict(delta=obs, p=p, relabelings=int(len(stats)), min_p=float(1 / max(ok.sum(), 1)) * 2,
                per_type=[None if not np.isfinite(x) else float(x) for x in per])


def perm_r(x, y, cells, rng, n_perm: int) -> dict:
    """Pearson's r of x and y after each cell's mean is removed from both, and its two-sided p from
    permutations of y among the units of each cell."""
    x, y, cells = np.asarray(x, float), np.asarray(y, float), np.asarray(cells)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y, cells = x[ok], y[ok], cells[ok]
    groups = [np.where(cells == c)[0] for c in pd.unique(cells)]
    xc, yc = x.copy(), y.copy()
    for g in groups:
        xc[g] -= x[g].mean()
        yc[g] -= y[g].mean()
    r = pearson(xc, yc)
    if not np.isfinite(r):
        return dict(r=r, p=float("nan"), n=int(ok.sum()), cells=len(groups), permutations=0)
    hits = 0
    yp = yc.copy()
    for _ in range(n_perm):
        for g in groups:
            if len(g) > 1:
                yp[g] = yc[rng.permutation(g)]
        rp = pearson(xc, yp)
        hits += int(np.isfinite(rp) and abs(rp) >= abs(r) - 1e-12)
    return dict(r=r, p=float((1 + hits) / (1 + n_perm)), n=int(ok.sum()), cells=len(groups), permutations=n_perm)


def fib_like(label: str) -> bool:
    return any(w in str(label).lower() for w in S.FIB_WORDS)


# --------------------------------------------------------------------------- #
# (b) §4's technical counts
# --------------------------------------------------------------------------- #
def b_tables(data: SimpleData, pool: Pool, log) -> dict:
    """(b1)-(b3) on the males of case 1."""
    obs = data.obs
    male = (obs["sex"] == "male").to_numpy()
    d, p = subset(data, pool, male)
    o = d.obs
    unit, age, ctype = (o[c].astype(str).to_numpy() for c in ("mouse", "age", "cell_type"))
    units = sorted(set(unit))
    young_u = np.array([bool((age[unit == u] == "young").all()) for u in units])
    tot = np.asarray(p.X.sum(axis=1), dtype=float).ravel()
    tot = np.where(tot > 0, tot, 1.0)
    col = {g: j for j, g in enumerate(p.var_names)}
    present = [g for g in S.FIB_GENES if g in col]
    missing = [g for g in S.FIB_GENES if g not in col]

    def counts(genes, rows):
        js = [col[g] for g in genes if g in col]
        if not js:
            return np.zeros(len(rows))
        return np.asarray(p.X[rows][:, js].sum(axis=1), dtype=float).ravel()

    ia, ib = d.var_names.index(S.PAIR[0]), d.var_names.index(S.PAIR[1])
    types, rows_ut = [], {}
    for t in sorted(set(ctype)):
        n_y = sum((((unit == u) & (ctype == t)).sum() >= S.MIN_UNIT_CELLS) for u, yy in zip(units, young_u) if yy)
        n_o = sum((((unit == u) & (ctype == t)).sum() >= S.MIN_UNIT_CELLS) for u, yy in zip(units, young_u) if not yy)
        if n_y >= S.MIN_UNITS and n_o >= S.MIN_UNITS:
            types.append(t)
            for u in units:
                r = np.where((unit == u) & (ctype == t))[0]
                if len(r) >= S.MIN_UNIT_CELLS:
                    rows_ut[(u, t)] = r
    nonfib = [t for t in types if not fib_like(t)]
    fibt = [t for t in types if fib_like(t)]
    log(f"(b) males: {len(units)} units ({int(young_u.sum())} young, {int((~young_u).sum())} old); {len(types)} types "
        f"with >= {S.MIN_UNIT_CELLS} cells in >= {S.MIN_UNITS} units of each group: {len(fibt)} fibroblast-like "
        f"({', '.join(fibt) or 'none'}), {len(nonfib)} others; the list's genes not measured: {missing or 'none'}")

    recs = {}
    for (u, t), r in rows_ut.items():
        rec = dict(unit=u, type=t, age="young" if young_u[units.index(u)] else "old", cells=int(len(r)),
                   fibroblast_like=fib_like(t))
        for g in present:
            c = np.asarray(p.X[r][:, col[g]].todense()).ravel()
            rec[f"det:{g}"] = float(np.mean(c > 0))
            rec[f"lvl:{g}"] = float(np.mean(np.log1p(c / tot[r] * 1e4)))
        rec["F"] = float(np.mean(np.log1p(counts(S.INDEX_GENES, r) / tot[r] * 1e4)))
        rec["F_plus"] = float(np.mean(np.log1p(counts(S.INDEX_PLUS, r) / tot[r] * 1e4)))
        sd = SimpleData(d.X[r][:, [ia, ib]], o.iloc[r], [S.PAIR[0], S.PAIR[1]])
        rec["MI"] = float(metrics.mi_3bin(sd, gene_a=S.PAIR[0], gene_b=S.PAIR[1]))
        recs[(u, t)] = rec

    def matrix(key, ts):
        return np.array([[recs[(u, t)][key] if (u, t) in recs else np.nan for u in units] for t in ts], dtype=float)

    out = dict(units=units, young=[u for u, y in zip(units, young_u) if y], types=types, fibroblast_like=fibt,
               non_fibroblast=nonfib, genes_missing=missing, index_genes=list(S.INDEX_GENES),
               index_plus=list(S.INDEX_PLUS), unit_rows=list(recs.values()))
    b1 = {}
    for g in present:
        b1[g] = {}
        for what in ("det", "lvl"):
            b1[g][what] = dict(non_fibroblast=exact_combined(matrix(f"{what}:{g}", nonfib), young_u),
                               fibroblast_like=exact_combined(matrix(f"{what}:{g}", fibt), young_u) if fibt else None)
    out["b1"] = b1
    out["b2"] = {k: dict(non_fibroblast=exact_combined(matrix(k, nonfib), young_u),
                         fibroblast_like=exact_combined(matrix(k, fibt), young_u) if fibt else None)
                 for k in ("F", "F_plus")}
    out["mi_combined"] = dict(non_fibroblast=exact_combined(matrix("MI", nonfib), young_u),
                              fibroblast_like=exact_combined(matrix("MI", fibt), young_u) if fibt else None)
    rng = np.random.default_rng(S.SEED)
    nf = [r for r in recs.values() if not r["fibroblast_like"]]
    mi = [r["MI"] for r in nf]
    f = [r["F"] for r in nf]
    out["b3_i"] = perm_r(mi, f, [f"{r['type']}|{r['age']}" for r in nf], rng, S.N_PERM)
    out["b3_0"] = perm_r(mi, f, [r["type"] for r in nf], rng, S.N_PERM)
    dmi = np.array(out["mi_combined"]["non_fibroblast"]["per_type"], dtype=float)
    dF = np.array(out["b2"]["F"]["non_fibroblast"]["per_type"], dtype=float)
    out["b3_ii"] = perm_r(dmi, dF, ["all"] * len(dmi), rng, S.N_PERM)
    out["per_type"] = [dict(type=t, fibroblast_like=fib_like(t),
                            units_young=int(sum((u, t) in recs and recs[(u, t)]["age"] == "young" for u in units)),
                            units_old=int(sum((u, t) in recs and recs[(u, t)]["age"] == "old" for u in units)),
                            d_MI=(out["mi_combined"]["non_fibroblast" if not fib_like(t) else "fibroblast_like"]["per_type"]
                                  [(nonfib if not fib_like(t) else fibt).index(t)]),
                            d_F=(out["b2"]["F"]["non_fibroblast" if not fib_like(t) else "fibroblast_like"]["per_type"]
                                 [(nonfib if not fib_like(t) else fibt).index(t)]))
                       for t in types]
    return out


def column_check(frame: pd.DataFrame, unit_col: str, sex_of: dict, age_of: dict, exclude, log, what: str) -> list:
    """Each column against the ages, within each sex: does it separate the young units from the old? The
    units' sex and age are case 1's (`sex_of`, `age_of`), so a frame without those columns is checked too."""
    rows = []
    fr = frame[frame[unit_col].astype(str).isin(set(age_of))].copy()
    fr["_unit"] = fr[unit_col].astype(str)
    ids = set()
    for c in fr.columns:
        if c == "_unit":
            continue
        g = fr.groupby("_unit")[c].nunique(dropna=False)
        back = fr.groupby(c)["_unit"].nunique(dropna=False) if fr[c].notna().any() else pd.Series(dtype=int)
        if len(g) and (g == 1).all() and len(back) and (back == 1).all() and len(back) == len(g):
            ids.add(c)  # one value per unit, one unit per value: an id of the unit
    for c in fr.columns:
        if c in ("_unit", unit_col) or c in exclude or c.lower() in ("age", "sex", "sex_ontology_term_id") or c in ids:
            continue
        for sex in sorted(set(sex_of.values())):
            units = sorted(u for u in fr["_unit"].unique() if sex_of.get(u) == sex)
            s = fr[fr["_unit"].isin(units)]
            ages = {u: age_of[u] for u in units}
            if sum(a == "young" for a in ages.values()) < 1 or sum(a == "old" for a in ages.values()) < 1:
                continue
            vals = s[c]
            num = pd.to_numeric(vals, errors="coerce")
            if num.notna().all() and vals.nunique() > 10:
                med = s.assign(_v=num).groupby("_unit")["_v"].median()
                young = np.array([ages[u] == "young" for u in med.index.astype(str)])
                ex = exact_combined(med.to_numpy(float)[None, :], young)
                rows.append(dict(source=what, column=c, sex=sex, kind="numeric", unit_medians=med.round(6).to_dict(),
                                 young_minus_old=ex["delta"], p=ex["p"], relabelings=ex["relabelings"]))
                continue
            by = s.assign(_v=vals.astype(str)).groupby("_unit")["_v"].agg(lambda x: sorted(set(x)))
            yv = {v for u, vs in by.items() if ages[str(u)] == "young" for v in vs}
            ov = {v for u, vs in by.items() if ages[str(u)] == "old" for v in vs}
            share = pd.Series([v for vs in by for v in vs]).value_counts()
            n_vals = len(yv | ov)
            disjoint = not (yv & ov)
            shared = bool((share >= 2).any())
            outcome = ("one value" if n_vals < 2 else "separates the ages (shared within an age)" if disjoint and shared
                       else "nested within units" if disjoint else "shared between the ages")
            rows.append(dict(source=what, column=c, sex=sex, kind="categorical", values=n_vals,
                             young_values=sorted(yv)[:12], old_values=sorted(ov)[:12], outcome=outcome))
    sep = [r for r in rows if r.get("outcome") == "separates the ages (shared within an age)"]
    log(f"(b4) {what}: {len(rows)} column x sex checks; the unit's ids left out: {sorted(ids) or 'none'}; columns that "
        f"separate the ages: {sorted({(r['column'], r['sex']) for r in sep}) or 'none'}")
    return rows


def unit_pattern(uid: str) -> list:
    """The id as a whole token with '_', '-' or '.' between its parts, also with its sex before its number."""
    parts = uid.split("_")
    orders = [parts] + ([[parts[0], parts[2], parts[1]]] if len(parts) == 3 else [])
    sep = r"[-_.]"
    return [re.compile(r"(?<![A-Za-z0-9])" + sep.join(map(re.escape, o)) + r"(?![A-Za-z0-9])", re.I) for o in orders]


def soft_text(path: Path) -> tuple[dict, list]:
    """GSM -> its !Sample_ lines' values joined; and the !Series_ lines."""
    samples, series, cur = {}, [], None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                cur = line.split("=", 1)[1].strip()
                samples[cur] = []
            elif line.startswith("^"):
                cur = None
            elif line.startswith("!Series_"):
                series.append(line)
            elif cur is not None and line.startswith("!Sample_") and " = " in line:
                samples[cur].append(line.split(" = ", 1)[1])
    return {k: " | ".join(v) for k, v in samples.items()}, series


def geo_check(geo: Path, log) -> dict:
    out = {}
    for gse in (S.GEO_TM, S.GEO_TMS):
        path = geo / f"{gse}_family.soft.gz"
        if not path.exists():
            out[gse] = dict(error="missing")
            continue
        samples, series = soft_text(path)
        hits = {}
        for uid in (*S.YOUNG_UNITS, *S.OLD_UNITS):
            pats = unit_pattern(uid)
            hits[uid] = sorted(g for g, txt in samples.items() if any(pt.search(txt) for pt in pats))
        quote = [s for s in series if s.startswith(("!Series_summary", "!Series_overall_design"))
                 and re.search(r"3[- ]?month|\b3m\b|tabula muris|previously", s, re.I)]
        out[gse] = dict(samples=len(samples), sha256=sha256_file(path), units_named={u: len(v) for u, v in hits.items()},
                        examples={u: v[:3] for u, v in hits.items() if v}, series_lines=quote[:12],
                        titles=[samples[g][:160] for g in list(samples)[:3]])
        log(f"(b4) GEO {gse}: {len(samples)} samples; units named: "
            + ", ".join(f"{u} {len(v)}" for u, v in hits.items()) + f"; series lines quoted: {len(quote)}")
    tm = out.get(S.GEO_TM, {}).get("units_named")
    if not tm:
        score = "undecided (GSE109774's records not read)"
    elif not any(tm.values()):
        score = "undecided: GSE109774's records name none of case 1's units"
    elif all(tm[u] for u in S.YOUNG_UNITS) and not any(tm[u] for u in S.OLD_UNITS):
        score = "holds: GSE109774's records name every 3-month unit of case 1 and no old one"
    else:
        score = ("does not hold: GSE109774's records name "
                 + ", ".join(f"{u}" for u in (*S.YOUNG_UNITS, *S.OLD_UNITS) if tm[u])
                 + (f"; missing 3-month units: {[u for u in S.YOUNG_UNITS if not tm[u]]}" if not all(
                     tm[u] for u in S.YOUNG_UNITS) else ""))
    out["statement_score"] = score
    return out


EXCLUDE_COLUMNS = {"soma_joinid", "observation_joinid", "donor_id", "development_stage",
                   "development_stage_ontology_term_id", "cell_type", "cell_type_ontology_term_id", "tissue",
                   "tissue_ontology_term_id", "tissue_general", "tissue_general_ontology_term_id", "index"}


def b_part(args, log) -> dict:
    males = (lambda o: o["sex"].astype(str).to_numpy() == "male")
    data, pool, sha = load_case(args.inputs, args.full, "case1", args.synthetic, rows=males)
    out = dict(case1_sha256=sha)
    t0 = time.time()
    try:
        out.update(b_tables(data, pool, log))
    except Exception as exc:
        out["tables_error"] = dict(error=repr(exc), traceback=traceback.format_exc().splitlines())
        log(f"EXCEPTION in (b1)-(b3): {exc!r}")
    out["seconds_tables"] = round(time.time() - t0, 1)
    del data, pool
    meta = []
    z = np.load(args.inputs / "data" / "case1.npz", allow_pickle=False)  # every unit's sex and age: the obs only
    o = pd.DataFrame({k[4:]: z[k] for k in z.files if k.startswith("obs_")}).astype(str)
    sex_of = dict(zip(o["mouse"], o["sex"]))
    age_of = dict(zip(o["mouse"], o["age"]))
    try:
        co = pd.read_csv(args.full / "census_obs_case1.csv.gz", low_memory=False)
        co = co.set_index("soma_joinid").loc[o["soma_joinid"].astype(np.int64)].reset_index()
        co["_case_unit"] = o["mouse"].to_numpy()
        meta += column_check(co, "_case_unit", sex_of, age_of, EXCLUDE_COLUMNS, log, "the Census's obs")
    except Exception as exc:
        out["census_obs_error"] = repr(exc)
        log(f"EXCEPTION in (b4), the Census's obs: {exc!r}")
    try:
        so = pd.read_csv(args.full / "source_obs_case1.csv.gz", low_memory=False)
        ucol = next(c for c in ("donor_id", "mouse.id", "mouse_id") if c in so.columns)
        so = so[so[ucol].astype(str).isin(set(age_of))]
        meta += column_check(so, ucol, sex_of, age_of, EXCLUDE_COLUMNS, log, "the source h5ad's obs")
        out["source_unit_column"] = ucol
        out["source_index_examples"] = {a: so[so[ucol].astype(str).map(age_of) == a].iloc[:3, 0].astype(str).tolist()
                                        for a in ("young", "old")}
    except Exception as exc:
        out["source_obs_error"] = repr(exc)
        log(f"EXCEPTION in (b4), the source h5ad's obs: {exc!r}")
    out["b4_columns"] = meta
    try:
        out["b4_geo"] = geo_check(args.inputs / "geo", log)
    except Exception as exc:
        out["b4_geo"] = dict(error=repr(exc), statement_score="undecided (an exception)")
        log(f"EXCEPTION in (b4), GEO: {exc!r}")
    return out


# --------------------------------------------------------------------------- #
# (a) the rule of the domain
# --------------------------------------------------------------------------- #
def domain_rows(dom: dict) -> list:
    out = []
    for s in dom.get("strata", []):
        for r in s["types"]:
            out.append(dict(stratum=s["stratum"], type=r["type"], cells=r["cells"], qualifies=r["qualifies"],
                            detection=r["detection"], in_domain=r["in_domain"], absent_from=r["absent_from"],
                            fibroblast_like=fib_like(r["type"])))
    return out


def a_part(args, log) -> dict:
    out = {}
    data, sha = load_case(args.inputs, args.full, "case1", args.synthetic, with_pool=False)
    male = (data.obs["sex"] == "male").to_numpy()
    d = SimpleData(data.X[male], data.obs[male].reset_index(drop=True), data.var_names)
    del data
    dom = claim_domain(d, type_col="cell_type", group_col="age", groups=("young", "old"), replicate_col="mouse",
                       genes=list(S.PAIR), kind="pair")
    types = sorted({t for v in dom["domain"].values() for t in v})
    ok = bool(types) and all(fib_like(t) for t in types)
    out["a1"] = dict(rule=dom["rule"], domain=types, outside=sorted({t for v in dom["outside"].values() for t in v}),
                     rows=domain_rows(dom), score=("holds" if ok else "does not hold") + f": the domain is {types or 'empty'}")
    log(f"(a1) §4 in males: the domain {types or 'empty'}; {out['a1']['score']}")
    mdata, mpool = load_marsh(args.full)
    arm = mdata.obs["arm"].to_numpy()
    out["a2"] = {}
    for c, (groups, direction, why) in S.CONTRASTS.items():
        mask = np.isin(arm, groups)
        dd, _ = subset(mdata, mpool, mask)
        dom = claim_domain(dd, type_col="cell_type", group_col="arm", groups=groups, replicate_col="mouse",
                           genes=list(S.MARSH_PAIR), kind="pair")
        row = domain_rows(dom)[0] if domain_rows(dom) else None
        cols = [dd.var_names.index(g) for g in S.MARSH_PAIR]
        pooled = {g: float(np.mean(dd.X[:, j] > 0)) for g, j in zip(S.MARSH_PAIR, cols)}
        alt = dict(both_groups=bool(row and all(all(row["detection"][g][n] >= DOMAIN_DETECTION for n in S.MARSH_PAIR)
                                                for g in groups)),
                   pooled_cells=bool(all(v >= DOMAIN_DETECTION for v in pooled.values())))
        in_dom = bool(row and row["in_domain"] and row["qualifies"])
        out["a2"][c] = dict(groups=list(groups), row=row, pooled_detection=pooled, in_domain=in_dom,
                            alternatives_not_used=alt)
        log(f"(a2) Marsh {c} {groups[0]} -> {groups[1]}: in the domain {in_dom}; detection "
            + (", ".join(f"{g}: " + ", ".join(f"{n} {row['detection'][g][n]:.3f}" for n in S.MARSH_PAIR) for g in groups)
               if row else "none") + f"; under 'both groups' {alt['both_groups']}, 'pooled cells' {alt['pooled_cells']}")
    out["a2_score"] = ("holds" if out["a2"]["M1"]["in_domain"] else "does not hold") + ": M1's type " + (
        "is" if out["a2"]["M1"]["in_domain"] else "is not") + " in its domain"
    return out


# --------------------------------------------------------------------------- #
# the acceptance tests
# --------------------------------------------------------------------------- #
def lcm(data, *, genes):
    """R1's metric on the named genes: the mean over cells of log1p CP10k of their summed counts."""
    names = list(map(str, data.var_names))
    return log_cp10k_mean(data, cols=[names.index(g) for g in genes])


def bind(claim: dict):
    """(metric, pair_metric, gene_metric, gene_pair)."""
    m, genes = claim["metric"], list(claim["genes"])
    if m in ("mi_3bin", "norm_pearson"):
        fn = getattr(metrics, m)
        return partial(fn, gene_a=genes[0], gene_b=genes[1]), fn, None, (genes[0], genes[1])
    if m == "log_cp10k_mean":
        return partial(lcm, genes=genes), None, lcm, None
    raise ValueError(f"no binding for {m!r}")


def probe(claim: dict):
    if claim["probe"] == "coupling":
        return injected_signal.coupling(*claim["genes"])
    return injected_signal.module(list(claim["genes"]))


def label_of(verdict: str) -> str:
    return verdict.split(" — ")[0].split(" [")[0].split(" (")[0].strip()


def summarize(rep: dict) -> dict:
    """What the log and the scores read from an engine's report."""
    ef = (rep["fields"].get("effect") or {})
    det = ef.get("detail") or {}
    dom = det.get("domain") or {}
    g4 = next((g for g in rep["gates"] if g["gate"] == 4), None)
    g4d = (g4 or {}).get("detail") or {}
    spec = rep.get("specificity") or {}
    keep = ("effect", "raw_effect", "ci", "p_perm", "p_raw", "p_t", "raw_detected", "tost_equivalent", "sesoi",
            "effective_sesoi", "sesoi_source", "floor", "tier", "test", "n_replicates", "estimate_over",
            "detected_only_after_correction", "explained_by_depth", "domain_empty")
    return dict(
        verdict=rep["verdict"], label=label_of(rep["verdict"]), cause=rep["cause"],
        fields={k: (v and dict(status=v["status"], reason=v["reason"], flags=v["flags"])) for k, v in rep["fields"].items()},
        gates=[dict(gate=g["gate"], name=g["name"], status=g["status"], message=g["message"]) for g in rep["gates"]],
        effect={k: det.get(k) for k in keep if k in det},
        domain=dict(types=dom.get("types"), empty=dom.get("empty"), n_units=dom.get("n_units"),
                    outside=dom.get("outside"), per_type=dom.get("per_type")) if dom else None,
        mixture=det.get("mixture") and {k: det["mixture"].get(k) for k in ("corrected", "raw")},
        gate4=dict(status=(g4 or {}).get("status"), sesoi_metric=g4d.get("sesoi_metric"),
                   sesoi_metric_ci=g4d.get("sesoi_metric_ci"), dose=g4d.get("dose"),
                   dose_response=g4d.get("dose_response")) if g4 else None,
        specificity={k: spec.get(k) for k in ("applicable", "error", "t_claim", "n_controls", "in_direction", "absolute",
                                              "nearest", "shift", "pool")} | (
            {"at_least_as_far": (spec.get("at_least_as_far") or [])[:40]} if spec else {}),
        reference=rep.get("reference"),
        composition={c: {k: v.get(k) for k in ("tvd", "tvd_p", "differs")}
                     for c, v in ((rep.get("composition") or {}).get("columns") or {}).items()})


def accept_run(claim: dict, data: SimpleData, pool: Pool, out_dir: Path, log) -> dict:
    metric, pair_metric, gene_metric, gene_pair = bind(claim)
    prereg = dict(estimand="composition", direction=claim["direction"], sesoi=float(claim["sesoi"]),
                  judgment_pending=False)
    if claim.get("known_to_differ"):
        prereg["known_to_differ"] = list(S.SEX_GENES)
    ctl = bool(claim.get("controls"))
    t0 = time.time()
    a = run_autopsy(metric, data, group_col=claim["group_col"], groups=tuple(claim["groups"]), gene_pair=gene_pair,
                    within=list(claim["within"]), pair_metric=pair_metric,
                    pos_pair=S.S4_POS if ctl else None, neg_pair=S.S4_NEG if ctl else None, prereg=prereg,
                    replicate_col=claim["replicate"], signal_test=probe(claim), seed=0, log_path="off",
                    composition=list(claim["composition"]), gene_metric=gene_metric, specificity=True,
                    specificity_pool=pool)
    a.metric_name = claim["metric"]
    rep = a.to_dict()
    try:
        text = json.dumps(rep, indent=1, allow_nan=False)
    except ValueError:  # a NaN in the report: kept as it is, not lost
        text = json.dumps(rep, indent=1, default=str)
    (out_dir / "reports").mkdir(parents=True, exist_ok=True)
    (out_dir / "reports" / f"{claim['id']}.json").write_text(text)
    rec = summarize(json.loads(text))
    rec.update(id=claim["id"], claim={k: v for k, v in claim.items() if k != "part"}, prereg=prereg,
               report=f"reports/{claim['id']}.json", report_sha256=hashlib.sha256(text.encode()).hexdigest(),
               seconds=round(time.time() - t0, 1), cells=int(data.n_obs))
    ef = rec["effect"]
    log(f"- {claim['id']}: {rec['verdict']} [cause {rec['cause']}] ({rec['seconds']} s)",
        f"    effect {ef.get('effect')} CI {ef.get('ci')} over {ef.get('estimate_over')}; s_m "
        f"{(rec['gate4'] or {}).get('sesoi_metric')}; domain {(rec['domain'] or {}).get('types')}",
        "    gates: " + "; ".join(f"{g['gate']} {g['status']}" for g in rec["gates"]))
    return rec


def case3_claims(inp: Path, prefix: str) -> list:
    spec1 = json.loads((inp / "data" / "select1.json").read_text())
    out = []
    for c in spec1["cases"]["case3"]["claims"]:
        if not c["id"].startswith(prefix):
            continue
        out.append(dict(id=f"A8-{c['id']}", part=f"case{prefix}", case="case3", metric=c["metric"],
                        genes=list(c["genes"]), groups=("A", "B"), group_col="group", replicate="donor", within=[],
                        direction=c["direction"], sesoi=float(c["sesoi"]),
                        probe=c["signal"]["kind"], composition=["cell_type"], criterion="EFFECT BELOW SESOI or none"))
    return out


def accept_part(args, part: str, log) -> dict:
    out = dict(claims=[])
    if part == "case3a" or part == "case3b" or part == "case3c":
        claims = case3_claims(args.inputs, part[4:])
    else:
        claims = [c for c in S.ACCEPT if c["part"] == part]
    if part == "marsh":
        data, pool = load_marsh(args.full)
        out["marsh_caveat"] = S.MARSH_CAVEAT
        out["parametric_only_note"] = ("3 mice per arm: the exact permutation over C(6, 3) = 20 relabelings cannot "
                                       "reach alpha 0.05 (its smallest two-sided p is 0.1); the replicate-level t decides")
    else:
        name = {"case2": "case2", "case1m": "case1", "case1": "case1"}.get(part, "case3")
        males = (lambda o: o["sex"].astype(str).to_numpy() == "male") if part == "case1m" else None
        data, pool, sha = load_case(args.inputs, args.full, name, args.synthetic, rows=males)
        out["case_sha256"] = sha
    for claim in claims:
        d, p = data, pool
        if claim["case"] == "marsh":
            d, p = subset(data, pool, np.isin(data.obs["arm"].to_numpy(), claim["groups"]))
        try:
            out["claims"].append(accept_run(claim, d, p, args.out, log))
        except Exception as exc:
            out["claims"].append(dict(id=claim["id"], claim=claim, error=repr(exc),
                                      traceback=traceback.format_exc().splitlines()))
            log(f"- {claim['id']}: EXCEPTION {exc!r}")
    return out


# --------------------------------------------------------------------------- #
# the scores
# --------------------------------------------------------------------------- #
def score_claims(recs: dict) -> dict:
    def lab(i):
        r = recs.get(i)
        return None if r is None or "error" in r else r["label"]
    out = {}
    for i in ("A1", "A2", "A5"):
        out[i] = dict(criterion="SUPPORTED", label=lab(i), met=lab(i) == "SUPPORTED" if lab(i) else None)
    for i in ("A4", "A7"):
        out[i] = dict(criterion="not SUPPORTED", label=lab(i), met=(lab(i) != "SUPPORTED") if lab(i) else None)
    r = recs.get("A3")
    if r is None or "error" in r:
        out["A3"] = dict(criterion="not SUPPORTED; no effect in the domain; the drop in the mixture from outside", met=None)
    else:
        f = r["fields"].get("effect") or {}
        dom = r.get("domain") or {}
        mix = ((r.get("mixture") or {}).get("corrected")) or {}
        outside = ((dom.get("outside") or {}).get("corrected")) or {}
        c1 = r["label"] != "SUPPORTED"
        c2 = bool(dom.get("types")) and not dom.get("empty") and f.get("status") != "DETECTED"
        c3 = bool(mix.get("detected") and (mix.get("effect") or 0) > 0 and outside.get("detected")
                  and (outside.get("effect") or 0) > 0)
        out["A3"] = dict(criterion="not SUPPORTED; no effect in the domain; the drop in the mixture from outside",
                         label=r["label"], not_supported=c1, no_effect_in_domain=c2, drop_outside=c3,
                         mixture=mix, outside=outside, met=bool(c1 and c2 and c3))
    a8 = {k: v for k, v in recs.items() if k.startswith("A8-")}
    labels = {k: (v.get("label") if "error" not in v else "EXCEPTION") for k, v in a8.items()}
    out["A8"] = dict(n=len(a8), labels=labels,
                     strict=bool(a8) and all(lb in ("EFFECT BELOW SESOI", "NO DETECTABLE EFFECT") for lb in labels.values()),
                     lenient=bool(a8) and all(lb != "SUPPORTED" and lb != "EXCEPTION" for lb in labels.values()),
                     supported=sorted(k for k, lb in labels.items() if lb == "SUPPORTED"))
    return out


def score_b(b: dict) -> dict:
    if "b2" not in b:
        return dict(score="not scored: (b1)-(b3) raised", part1=None, part2=None)
    f = b["b2"]["F"]["non_fibroblast"]
    r = b["b3_i"]
    p1 = bool(np.isfinite(f["delta"]) and f["delta"] > 0 and f["p"] < S.ALPHA)
    p2 = bool(np.isfinite(r["r"]) and r["r"] > 0 and r["p"] < S.ALPHA)
    score = "holds" if p1 and p2 else "partly holds" if (p1 or p2) else "does not hold"
    return dict(score=score, part1=p1, part2=p2, F=f, b3_i=r)


# --------------------------------------------------------------------------- #
# the log of a part, and the merge
# --------------------------------------------------------------------------- #
def fmt(x, n=4):
    return "nan" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:+.{n}g}"


def b_lines(b: dict) -> list:
    if "b2" not in b:
        return ["(b) not computed: " + str(b.get("tables_error", {}).get("error"))]
    out = ["### (b1) detection and level of the list's genes, young minus old, combined over the non-fibroblast types",
           "| gene | detection Δ | p | level Δ | p | fibroblast-like: detection Δ | level Δ |", "|---|---|---|---|---|---|---|"]
    for g, v in b["b1"].items():
        nf, fl = v["det"]["non_fibroblast"], v["det"]["fibroblast_like"] or {}
        ln, ll = v["lvl"]["non_fibroblast"], v["lvl"]["fibroblast_like"] or {}
        out.append(f"| {g} | {fmt(nf['delta'])} | {nf['p']:.3g} | {fmt(ln['delta'])} | {ln['p']:.3g} | "
                   f"{fmt(fl.get('delta'))} | {fmt(ll.get('delta'))} |")
    for k, name in (("F", "F (the index genes)"), ("F_plus", "F⁺ (with Col1a1)")):
        e = b["b2"][k]["non_fibroblast"]
        out.append(f"(b2) {name}: young minus old {fmt(e['delta'])}, exact p {e['p']:.4g} over {e['relabelings']} "
                   "relabelings (non-fibroblast types)")
    e = b["mi_combined"]["non_fibroblast"]
    out.append(f"MI within the non-fibroblast types: young minus old {fmt(e['delta'])}, exact p {e['p']:.4g}")
    for k, name in (("b3_i", "(b3-i) within type and age, scored"), ("b3_0", "(b3-0) within type, reported"),
                    ("b3_ii", "(b3-ii) across types, reported")):
        r = b[k]
        out.append(f"{name}: r {fmt(r['r'], 3)}, p {r['p']:.4g} ({r['n']} values, {r['cells']} cells, "
                   f"{r['permutations']} permutations)")
    out += ["", "| type | fibroblast-like | units young/old | Δ MI | Δ F |", "|---|---|---|---|---|"]
    for t in b["per_type"]:
        out.append(f"| {t['type']} | {t['fibroblast_like']} | {t['units_young']}/{t['units_old']} | {fmt(t['d_MI'])} | "
                   f"{fmt(t['d_F'])} |")
    return out


def merge(parts_dir: Path, out_dir: Path, synthetic: bool) -> int:
    parts, logs = {}, []
    for part in S.PARTS:
        f = parts_dir / f"part-{part}" / f"{part}.json"
        if f.exists():
            parts[part] = json.loads(f.read_text())
            lg = parts_dir / f"part-{part}" / f"{part}.log"
            logs += lg.read_text().splitlines() if lg.exists() else []
        else:
            parts[part] = dict(missing=True)
    recs = {}
    for part in ("case2", "case1m", "case1", "marsh", "case3a", "case3b", "case3c"):
        for r in parts[part].get("result", {}).get("claims", []):
            recs[r["id"]] = r
    scores = dict(claims=score_claims(recs))
    b = parts["b"].get("result", {})
    a = parts["a"].get("result", {})
    scores["b_prediction"] = score_b(b) if b else dict(score="not scored: part b missing")
    scores["b_statement"] = (b.get("b4_geo") or {}).get("statement_score", "not scored")
    scores["a1"] = (a.get("a1") or {}).get("score", "not scored")
    scores["a2"] = a.get("a2_score", "not scored")
    tag = "SYNTHETIC " if synthetic else ""
    res = dict(synthetic=synthetic, parts={k: {kk: vv for kk, vv in v.items() if kk != "result"} for k, v in parts.items()},
               results={k: v.get("result") for k, v in parts.items()}, scores=scores)
    (out_dir / "real_pilot_5.json").write_text(json.dumps(res, indent=1, default=str))
    lines = [f"# {tag}The third diagnostics and v0.4's acceptance tests (validation/exploratory/"
             "real_pilot_diagnostics_3_plan.md)", "", "## The scores", ""]
    c = scores["claims"]
    lines += ["| test | criterion | verdict | met |", "|---|---|---|---|"]
    for k in ("A1", "A2", "A3", "A4", "A5", "A7"):
        v = c[k]
        lines.append(f"| {k} | {v['criterion']} | {v.get('label')} | {v.get('met')} |")
    lines.append(f"| A8 | strict: each EFFECT BELOW SESOI or NO DETECTABLE EFFECT; lenient: none SUPPORTED | "
                 f"{c['A8']['n']} claims, SUPPORTED: {c['A8']['supported'] or 'none'} | strict {c['A8']['strict']}, "
                 f"lenient {c['A8']['lenient']} |")
    if "no_effect_in_domain" in c["A3"]:
        lines.append(f"A3 in parts: not SUPPORTED {c['A3']['not_supported']}; no effect in the domain "
                     f"{c['A3']['no_effect_in_domain']}; the drop in the mixture and outside the domain {c['A3']['drop_outside']}")
    lines += ["", f"(b) the owner's prediction: {scores['b_prediction'].get('score')}",
              f"(b4) the owner's statement on Tabula Muris: {scores['b_statement']}",
              f"(a1) §4 in males: {scores['a1']}", f"(a2) Marsh: {scores['a2']}", "", "## The claims", ""]
    for k, r in recs.items():
        lines.append(f"- {k}: {r.get('verdict', r.get('error'))} [cause {r.get('cause')}]")
    lines += ["", "## (b)", ""] + (b_lines(b) if b else ["part b missing"])
    lines += ["", "## The parts' logs", ""] + logs
    (out_dir / "real_pilot_5.log").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:40]))
    return 0


# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--inputs", type=Path)
    p.add_argument("--full", type=Path)
    p.add_argument("--part", choices=S.PARTS)
    p.add_argument("--merge", type=Path, help="the directory of the parts' outputs (part-<name>/)")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--synthetic", action="store_true")
    args = p.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    if args.merge:
        return merge(args.merge, args.out, args.synthetic)
    lines = []

    def log(*msg):
        for m in msg:
            print(m, flush=True)
            lines.append(str(m))
    tag = "SYNTHETIC " if args.synthetic else ""
    log(f"# {tag}part {args.part} (validation/exploratory/real_pilot_diagnostics_3_plan.md)", "",
        f"engine {MA.__version__} from {Path(MA.__file__).resolve().parent}; python {platform.python_version()}; numpy "
        f"{np.__version__}; pandas {pd.__version__}; numerics {NUMERICS}", "")
    t0, status = time.time(), 0
    try:
        if args.part == "b":
            result = b_part(args, log)
        elif args.part == "a":
            result = a_part(args, log)
        else:
            result = accept_part(args, args.part, log)
    except Exception as exc:
        result = dict(error=repr(exc), traceback=traceback.format_exc().splitlines())
        log(f"EXCEPTION {exc!r}", *result["traceback"])
        status = 1
    if args.part == "b" and "b2" in result:
        log("", *b_lines(result))
    out = dict(part=args.part, synthetic=bool(args.synthetic), engine=MA.__version__, seconds=round(time.time() - t0, 1),
               result=result)
    (args.out / f"{args.part}.json").write_text(json.dumps(out, indent=1, default=str))
    log("", f"{tag}part {args.part}: {out['seconds']} s")
    (args.out / f"{args.part}.log").write_text("\n".join(lines) + "\n")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
