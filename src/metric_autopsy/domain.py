"""The domain of a claim: the cell types where it is judged (v0.4, decided 2026-10-18).

When the data carry cell-type labels (the first column of ``run_autopsy(composition=...)``), the
verdict comes from the claim's estimate combined over the types of its domain, the types weighted
equally (``combined``), and the estimate over the whole mixture goes into the report beside the
groups' difference of composition. The rule, fixed from principle before the third diagnostics
(``validation/exploratory/real_pilot_diagnostics_3_plan.md``), not from section 4's data:

* a type is **in the domain of a claim about a pair** where, in at least one of the two groups, both
  genes are detected in at least ``DOMAIN_DETECTION`` (10%) of the type's cells, the lower bound of
  validation v1's expression levels (no metric was tested below it). Detection in a group is the mean
  over its replicates, each weighted equally, of the share of the replicate's cells of the type with a
  count of the gene. *In at least one group*, because a claim of change may be the very appearance of
  the coupling: requiring both groups would remove a type where the change itself moves a gene across
  the threshold, and pooling the cells would make the rule depend on the groups' numbers of cells;
* for a claim about a set of genes (a module, or one gene), where the set (a count of any of its
  genes) is detected so in at least one group;
* for a claim that names no genes, every type;
* and only where it can be estimated: at least ``MIN_CELLS`` cells in at least ``MIN_REPLICATES``
  replicates of each group (in a paired design, replicates with that many cells in both groups).

The rule is applied within each stratum of ``within``. A type absent from every replicate of one group
cannot be compared there (change 7's UNIDENTIFIABLE) and is listed. Nothing here reads the metric.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np

from .combined import Units
from .composition import MIN_CELLS, MIN_REPLICATES

DOMAIN_DETECTION = 0.10


def _gene_columns(data, genes) -> np.ndarray:
    from .core import as_dense, unique_col_index
    names = list(map(str, data.var_names))
    cols = [unique_col_index(names, str(g)) for g in genes]
    X = data.X
    sub = X[:, cols]
    return np.asarray(as_dense(sub), dtype=float).reshape(X.shape[0], len(cols))


def claim_domain(data, *, type_col: str, group_col: str, groups: Sequence, replicate_col: str,
                 within: Sequence[str] = (), genes: Sequence[str] = (), kind: str = "global",
                 min_cells: int = MIN_CELLS, min_replicates: int = MIN_REPLICATES,
                 threshold: float = DOMAIN_DETECTION) -> dict:
    """The claim's domain: per stratum, every type with its cells, its replicates holding at least
    `min_cells` of them, the detection of the claim's genes in each group, whether it qualifies for an
    estimate and whether it is in the domain. `kind` is "pair", "set" or "global"."""
    from .effect import _replicate_design

    if kind not in ("pair", "set", "global"):
        raise ValueError(f"kind is 'pair', 'set' or 'global', not {kind!r}")
    obs = data.obs
    grp_all = np.asarray(obs[group_col])
    keep = (grp_all == groups[0]) | (grp_all == groups[1])
    sub_obs = obs[keep]
    design = _replicate_design(sub_obs, group_col, groups, replicate_col, list(within))
    rule = dict(threshold=threshold, kind=kind, genes=[str(g) for g in genes], min_cells=min_cells,
                min_replicates=min_replicates,
                text=("types where, in at least one group, " + (
                    "both genes are" if kind == "pair" else "the set of genes is" if kind == "set" else "")
                    + f" detected in at least {threshold:.0%} of the type's cells" if kind != "global"
                    else "every type (the claim names no genes)")
                + f"; estimable with at least {min_cells} cells in at least {min_replicates} replicates of each group")
    if design["kind"] == "mixed":
        return dict(column=type_col, rule=rule, error=design["reason"])
    lab = np.asarray(sub_obs[type_col]).astype(object)
    missing = np.array([v is None or (isinstance(v, float) and np.isnan(v)) for v in lab])
    lab = np.array([str(v) for v in lab], dtype=object)
    grp = grp_all[keep]
    rep = np.asarray(sub_obs[replicate_col]).astype(str)
    G = _gene_columns(data, genes)[keep] if kind != "global" and len(genes) else None
    if G is not None:
        detected = (G > 0) if kind == "pair" else (G.sum(axis=1, keepdims=True) > 0)
        names = [str(g) for g in genes] if kind == "pair" else ["set"]
    need = max(int(min_replicates), MIN_REPLICATES)
    g0, g1 = (str(x) for x in groups)
    strata_labels = list(dict.fromkeys(design["strata"][r] for r in design["reps"]))
    out_strata = []
    for st in strata_labels:
        reps_st = [r for r in design["reps"] if design["strata"][r] == st]
        in_st = np.isin(rep, reps_st) & ~missing
        rows = []
        for t in sorted(set(lab[in_st].tolist())):
            mt = in_st & (lab == t)
            row = dict(type=t, cells={}, replicates={}, detection={})
            units = {}
            for gname in (g0, g1):
                mg = mt & (grp.astype(str) == gname)
                counts = {r: int((mg & (rep == r)).sum()) for r in reps_st}
                row["cells"][gname] = int(mg.sum())
                ok = sorted(r for r, c in counts.items() if c >= min_cells)
                units[gname] = ok
                if G is not None:
                    det = [detected[mg & (rep == r)].mean(axis=0) for r in ok]
                    row["detection"][gname] = ({n: float(np.mean([d[i] for d in det])) for i, n in enumerate(names)}
                                               if det else {n: float("nan") for n in names})
            if design["kind"] == "paired":
                both = sorted(set(units[g0]) & set(units[g1]))
                row["replicates"] = dict(pairs=both)
                row["qualifies"] = len(both) >= need
            else:
                row["replicates"] = {g: units[g] for g in (g0, g1)}
                row["qualifies"] = all(len(units[g]) >= need for g in (g0, g1))
            row["absent_from"] = next((g for g in (g0, g1) if row["cells"][g] == 0), None)
            if kind == "global":
                row["in_domain"] = True
            else:
                row["in_domain"] = any(
                    all(np.isfinite(row["detection"][g][n]) and row["detection"][g][n] >= threshold for n in names)
                    for g in (g0, g1))
            rows.append(row)
        out_strata.append(dict(stratum=" ".join(st.split("\x1f")) if st else "all", key=st, types=rows))
    domain = {s["key"]: [r["type"] for r in s["types"] if r["qualifies"] and r["in_domain"]] for s in out_strata}
    outside = {s["key"]: [r["type"] for r in s["types"] if r["qualifies"] and not r["in_domain"]] for s in out_strata}
    return dict(column=type_col, rule=rule, design=design["kind"], strata=out_strata, domain=domain, outside=outside,
                n_domain_types=len({t for v in domain.values() for t in v}),
                unlabelled_cells=int(missing.sum()), _design=design)


def build_units(dom: dict, data, *, group_col: str, groups: Sequence, replicate_col: str,
                which: str = "domain", types: dict | None = None) -> Units | None:
    """The units (replicate, type) of `which` ("domain", "outside", or "all" qualifying types; or the
    explicit {stratum key: [types]} of `types`), with the cell indices of each unit in `data` (any data
    holding the claim's cells with the same obs columns: the cells are found by their labels)."""
    design = dom["_design"]
    obs = data.obs
    grp = np.asarray(obs[group_col]).astype(str)
    rep = np.asarray(obs[replicate_col]).astype(str)
    raw = np.asarray(obs[dom["column"]]).astype(object)
    missing = np.array([v is None or (isinstance(v, float) and np.isnan(v)) for v in raw])
    lab = np.array([str(v) for v in raw], dtype=object)
    idx_all = np.arange(len(obs))
    g0, g1 = (str(x) for x in groups)
    if types is None:
        if which == "domain":
            types = dom["domain"]
        elif which == "outside":
            types = dom["outside"]
        else:
            types = {k: sorted(set(dom["domain"].get(k, [])) | set(dom["outside"].get(k, []))) for k in dom["domain"]}
    all_types = sorted({t for v in types.values() for t in v})
    if not all_types:
        return None
    t_index = {t: i for i, t in enumerate(all_types)}
    reps = list(design["reps"])
    r_index = {r: i for i, r in enumerate(reps)}
    strata = list(dict.fromkeys(design["strata"][r] for r in reps))
    s_index = {s: i for i, s in enumerate(strata)}
    rep_first = np.array([1 if groups[0] in design["membership"][r] else 0 for r in reps]) if design["kind"] == "nested" \
        else np.zeros(len(reps), int)
    rep_stratum = np.array([s_index[design["strata"][r]] for r in reps])
    u_rep, u_type, cells = [], [], []
    by_stratum = {s["key"]: s for s in dom["strata"]}
    for st, ts in types.items():
        rows = {r["type"]: r for r in by_stratum[st]["types"]}
        for t in ts:
            row = rows[t]
            mt = (lab == t) & ~missing
            if design["kind"] == "paired":
                for r in row["replicates"]["pairs"]:
                    m = mt & (rep == r)
                    cells.append((idx_all[m & (grp == g0)], idx_all[m & (grp == g1)]))
                    u_rep.append(r_index[r]), u_type.append(t_index[t])
            else:
                for gname in (g0, g1):
                    for r in row["replicates"][gname]:
                        m = mt & (rep == r) & (grp == gname)
                        cells.append(idx_all[m])
                        u_rep.append(r_index[r]), u_type.append(t_index[t])
    return Units(design["kind"], reps, rep_first, rep_stratum, strata, all_types,
                 np.array(u_rep, int), np.array(u_type, int), cells)


def unit_values(metric, data, units: Units) -> np.ndarray:
    """The metric on each unit's cells (a paired unit: its first group's value minus its second's),
    NaN where the metric gives no finite value."""
    from .gates import _safe_call
    out = np.full(units.n_units, np.nan)
    for i, c in enumerate(units.cells):
        if units.kind == "paired":
            a, b = (_safe_call(lambda c=x: metric(data[c])) for x in c)
            if a is not None and b is not None:
                out[i] = a - b
        else:
            v = _safe_call(lambda: metric(data[c]))
            if v is not None:
                out[i] = v
    return out


def describe(dom: dict) -> list:
    """Lines for the report: per stratum, the domain's types and the types left outside it."""
    if dom.get("error"):
        return [f"the domain could not be built: {dom['error']}"]
    lines = [f"rule: {dom['rule']['text']}"]
    for s in dom["strata"]:
        inside = [r["type"] for r in s["types"] if r["qualifies"] and r["in_domain"]]
        outside = [r["type"] for r in s["types"] if r["qualifies"] and not r["in_domain"]]
        small = [r["type"] for r in s["types"] if not r["qualifies"] and r["in_domain"]]
        lines.append(f"{s['stratum']}: domain {', '.join(inside) or 'none'}; outside it {len(outside)} type(s)"
                     + (f"; in the domain by detection but not estimable: {', '.join(small)}" if small else ""))
    return lines
