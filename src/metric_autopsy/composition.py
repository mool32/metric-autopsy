"""The composition of a comparison in cell types and tissues (v0.4 change 7, the decided minimum).

v1 used one cell type per background, and v0.3 stratifies only by the factors a user declares, so
nothing in v0.3 names a difference of composition: on preprint section 4's setup the old females are
the mammary gland alone against 22 tissues in the young (``validation/exploratory/real_pilot_2.log``).
The decided minimum:

(i)   a diagnostic: the shares of each level of an annotation (a cell type, a tissue) per replicate,
      compared between the groups at the replicate level (``composition_difference``);
(ii)  a type absent from every replicate of one group cannot be compared, and is listed;
(iii) estimates per type;
(iv)  a claim about a mixture only together with its difference of composition (the verdict's note).

On section 4 in males the drop of MI stayed within the cell types (68.7% of it, exact p 0.01587;
``validation/exploratory/real_pilot_4.log``). The rule decided on 2026-10-18 is the claim's domain
(``domain``): the verdict follows the estimate combined over the types where the claim's genes are
expressed, the types weighted equally, and the whole mixture's estimate is reported beside the
difference of composition. This module reports; it removes nothing.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np

from .stats import nested_relabelings, relabeling_p

MIN_CELLS = 20       # a replicate counts for a type with at least this many of its cells
MIN_REPLICATES = 3   # a type is judged with at least this many such replicates in each group


def _labels(obs, col) -> np.ndarray:
    return np.asarray(obs[col]).astype(str)


def composition_difference(obs, col: str, group_col: str, groups: Sequence, replicate_col: str,
                           within: Sequence[str] = (), *, alpha: float = 0.05, n_perm: int = 1000,
                           rng=None, max_exact: int = 20000) -> dict:
    """The groups' composition in the levels of `col`, compared at the replicate level.

    Per replicate, the share of its cells (of the two groups) in each level. Per level, the mean share
    in each group (each replicate weighted equally) and the difference from groups[0] to groups[1],
    with its two-sided permutation p; as one number, the total variation distance between the groups'
    mean compositions, TVD = 1/2 sum |difference|, with its permutation p. The relabelings are those of
    the effect test: the group labels of the replicates within the strata of `within` (every one when
    there are at most `max_exact`, else `n_perm` draws); for replicates that hold both groups (a paired
    design) the per-replicate differences change sign. Levels absent from every replicate of one group
    are listed (``absent``)."""
    from .effect import _replicate_design

    rng = rng if rng is not None else np.random.default_rng(0)
    grp = np.asarray(obs[group_col])
    keep = (grp == groups[0]) | (grp == groups[1])
    obs = obs[keep]
    grp = grp[keep]
    if len(obs) == 0:
        return dict(column=col, error="no cells of the two groups")
    design = _replicate_design(obs, group_col, groups, replicate_col, within)
    if design["kind"] == "mixed":
        return dict(column=col, error=design["reason"])
    lab = _labels(obs, col)
    rep = np.asarray(obs[replicate_col]).astype(str)
    levels = sorted(set(lab.tolist()))
    index = {v: i for i, v in enumerate(levels)}
    codes = np.array([index[v] for v in lab])

    def shares(mask) -> np.ndarray:
        n = int(mask.sum())
        return np.bincount(codes[mask], minlength=len(levels)) / n if n else np.full(len(levels), np.nan)

    reps = design["reps"]
    present = {str(g): np.zeros(len(levels), bool) for g in groups}
    if design["kind"] == "nested":
        S = np.array([shares(rep == r) for r in reps])
        first = np.array([1 if groups[0] in design["membership"][r] else 0 for r in reps])
        for r, row in zip(reps, S):
            present[str(groups[0] if groups[0] in design["membership"][r] else groups[1])] |= row > 0
        st_labels = {st: i for i, st in enumerate(dict.fromkeys(design["strata"][r] for r in reps))}
        strata = np.array([st_labels[design["strata"][r]] for r in reps])
        relab, method = nested_relabelings(first, strata, n_perm=n_perm, rng=rng, max_exact=max_exact)

        def diff(labels: np.ndarray) -> np.ndarray:  # second minus first, precision weights over strata
            num, wsum = np.zeros(len(levels)), 0.0
            for st in np.unique(strata):
                m = strata == st
                a, b = m & (labels == 1), m & (labels == 0)
                if a.any() and b.any():
                    w = 1.0 / (1.0 / a.sum() + 1.0 / b.sum())
                    num += w * (S[b].mean(axis=0) - S[a].mean(axis=0))
                    wsum += w
            return num / wsum if wsum else np.full(len(levels), np.nan)
        observed = diff(first)
        null = np.array([diff(x) for x in relab])
        mean_first, mean_second = S[first == 1].mean(axis=0), S[first == 0].mean(axis=0)
        n_counts = {str(groups[0]): int(first.sum()), str(groups[1]): int((1 - first).sum())}
    else:
        D, A, B = [], [], []
        for r in reps:
            a, b = shares((rep == r) & (grp == groups[0])), shares((rep == r) & (grp == groups[1]))
            present[str(groups[0])] |= a > 0
            present[str(groups[1])] |= b > 0
            A.append(a), B.append(b), D.append(b - a)
        D, A, B = np.array(D), np.array(A), np.array(B)
        observed = D.mean(axis=0)
        n = len(reps)
        if 2 ** n <= max_exact:
            import itertools
            signs, method = np.array(list(itertools.product((1.0, -1.0), repeat=n))), "exact"
        else:
            signs, method = rng.choice((1.0, -1.0), size=(n_perm, n)), "monte-carlo"
        null = signs @ D / n
        mean_first, mean_second = A.mean(axis=0), B.mean(axis=0)
        n_counts = dict(pairs=n)
    tvd_obs = 0.5 * float(np.nansum(np.abs(observed)))
    tvd_null = 0.5 * np.nansum(np.abs(null), axis=1)
    p_tvd = relabeling_p(tvd_obs, tvd_null, method)
    rows = []
    for j, v in enumerate(levels):
        rows.append(dict(level=v, first=float(mean_first[j]), second=float(mean_second[j]),
                         difference=float(observed[j]),
                         p=relabeling_p(abs(observed[j]), np.abs(null[:, j]), method)))
    rows.sort(key=lambda x: -abs(x["difference"]))
    absent = {v: g for g in present for j, v in enumerate(levels)
              if not present[g][j] and any(present[h][j] for h in present if h != g)}
    return dict(column=col, groups=[str(g) for g in groups], design=design["kind"], replicates=n_counts,
                tvd=tvd_obs, tvd_p=p_tvd, method=method, n_relabelings=int(len(null)),
                differs=bool(p_tvd < alpha), levels=rows, absent=absent)


def type_replicates(obs, type_col: str, group_col: str, groups: Sequence, replicate_col: str,
                    min_cells: int = MIN_CELLS) -> dict:
    """Per type, per group: the replicates holding at least `min_cells` cells of it, and whether the
    type has any cell in a replicate of the group."""
    lab = _labels(obs, type_col)
    grp = np.asarray(obs[group_col])
    rep = np.asarray(obs[replicate_col]).astype(str)
    out = {}
    for t in sorted(set(lab[(grp == groups[0]) | (grp == groups[1])].tolist())):
        row = {}
        for g in groups:
            m = (lab == t) & (grp == g)
            counts = {r: int((m & (rep == r)).sum()) for r in dict.fromkeys(rep[m])}
            row[str(g)] = dict(any=bool(m.any()), cells=int(m.sum()),
                               replicates=sorted(r for r, c in counts.items() if c >= min_cells))
        out[t] = row
    return out


TECHNICAL = ("counts per cell", "genes per cell")


def _domain_cells(obs, dom: dict, replicate_col: str) -> np.ndarray:
    """The cells of the domain's types, each in its own stratum (the domain's per-stratum lists)."""
    rep = np.asarray(obs[replicate_col]).astype(str)
    lab = _labels(obs, dom["column"])
    strata = dom["_design"]["strata"]
    key = np.array([strata.get(r, None) for r in rep], dtype=object)
    out = np.zeros(len(obs), bool)
    for st, types in dom["domain"].items():
        out |= (key == st) & np.isin(lab, list(types))
    return out


def technical_covariates(data, group_col: str, groups: Sequence, replicate_col: str, within: Sequence[str] = (), *,
                         domain: dict | None = None, pool=None, alpha: float = 0.05, n_perm: int = 1000, rng=None,
                         max_exact: int = 20000) -> dict:
    """The groups' technical covariates per replicate (a line of the dossier, the owner's decision of
    2026-10-18): per cell its total count and its number of genes with a count, over every gene of
    `pool` (the same cells with every gene) or else of `data`; per replicate (per replicate and group in
    a paired design) the median over its cells; the difference from groups[0] to groups[1]
    (groups[0] − groups[1], as the effect) of the replicates' means within each stratum of `within`,
    and over the strata with the effect test's precision weights, each with its two-sided permutation
    p over the effect test's relabelings (every one when there are at most `max_exact`). On all the
    comparison's cells, and on the cells of the domain's types when a domain is given. GATE 1 compares
    the cells' median numbers of genes with a ratio threshold of 1.5; on preprint section 4's setup it
    passed at 1.32 while the young males' units detected 1,298 more genes per cell (exact p 0.016 over
    126 relabelings, the third diagnostics, ``validation/exploratory/real_pilot_5_summary.log``)."""
    from .effect import _replicate_design
    from .stats import permutation_test_nested, permutation_test_paired

    rng = rng if rng is not None else np.random.default_rng(0)
    obs = data.obs
    src = pool if pool is not None else data
    n = len(obs)
    if int(getattr(src.X, "shape", (0,))[0]) != n:
        return dict(error=f"the pool has {getattr(src.X, 'shape', ('?',))[0]} cells, the data {n}")
    grp = np.asarray(obs[group_col])
    keep = (grp == groups[0]) | (grp == groups[1])
    design = _replicate_design(obs[keep], group_col, groups, replicate_col, list(within))
    if design["kind"] == "mixed":
        return dict(error=design["reason"])
    X = src.X
    per_cell = {TECHNICAL[0]: np.asarray(X.sum(axis=1), dtype=float).ravel(),
                TECHNICAL[1]: np.asarray((X > 0).sum(axis=1), dtype=float).ravel()}
    rep = np.asarray(obs[replicate_col]).astype(str)
    cell_sets = {"all cells": keep}
    if domain is not None and not domain.get("error") and domain.get("n_domain_types"):
        cell_sets["the domain's cells"] = keep & _domain_cells(obs, domain, replicate_col)
    reps = design["reps"]
    strata = design["strata"]
    st_order = list(dict.fromkeys(strata[r] for r in reps))
    g0, g1 = groups[0], groups[1]
    rows = []
    for cells, cm in cell_sets.items():
        for cov, v in per_cell.items():
            if design["kind"] == "nested":
                vals, lab, st = [], [], []
                for r in reps:
                    m = cm & (rep == r)
                    if m.any():
                        vals.append(float(np.median(v[m])))
                        lab.append(1 if g0 in design["membership"][r] else 0)
                        st.append(strata[r])
                y, g, s = np.asarray(vals), np.asarray(lab, int), np.asarray(st, dtype=object)
                parts = [(name, s == key) for name, key in
                         ([(" ".join(k.split("\x1f")) or "all", k) for k in st_order] if len(st_order) > 1 else [])]
                parts.append(("combined" if len(st_order) > 1 else "all", np.ones(len(y), bool)))
                for name, m in parts:
                    if not ((g[m] == 1).any() and (g[m] == 0).any()):
                        rows.append(dict(cells=cells, covariate=cov, stratum=name, difference=None,
                                         reason="a group without replicates"))
                        continue
                    res = permutation_test_nested(y[m], g[m], s[m], n_perm=n_perm, rng=rng, max_exact=max_exact)
                    rows.append(dict(cells=cells, covariate=cov, stratum=name,
                                     first=float(y[m][g[m] == 1].mean()), second=float(y[m][g[m] == 0].mean()),
                                     difference=float(res.observed), p=float(res.p), method=res.method,
                                     n_relabelings=int(res.n), replicates={str(g0): int((g[m] == 1).sum()),
                                                                           str(g1): int((g[m] == 0).sum())}))
            else:
                d, st = [], []
                for r in reps:
                    a, b = cm & (rep == r) & (grp == g0), cm & (rep == r) & (grp == g1)
                    if a.any() and b.any():
                        d.append((float(np.median(v[a])), float(np.median(v[b]))))
                        st.append(strata[r])
                D, s = np.asarray(d, float).reshape(-1, 2), np.asarray(st, dtype=object)
                parts = [(name, s == key) for name, key in
                         ([(" ".join(k.split("\x1f")) or "all", k) for k in st_order] if len(st_order) > 1 else [])]
                parts.append(("combined" if len(st_order) > 1 else "all", np.ones(len(D), bool)))
                for name, m in parts:
                    if not m.any():
                        rows.append(dict(cells=cells, covariate=cov, stratum=name, difference=None,
                                         reason="no replicate with cells of both groups"))
                        continue
                    res = permutation_test_paired(D[m, 0] - D[m, 1], n_perm=n_perm, rng=rng, max_exact=max_exact)
                    rows.append(dict(cells=cells, covariate=cov, stratum=name, first=float(D[m, 0].mean()),
                                     second=float(D[m, 1].mean()), difference=float(res.observed), p=float(res.p),
                                     method=res.method, n_relabelings=int(res.n), replicates=dict(pairs=int(m.sum()))))
    return dict(groups=[str(g0), str(g1)], design=design["kind"], source="pool" if pool is not None else "data",
                alpha=alpha, rows=rows,
                differs=[dict(cells=r["cells"], covariate=r["covariate"], stratum=r["stratum"]) for r in rows
                         if r.get("difference") is not None and r["p"] < alpha])


def describe_technical(tech: dict) -> list:
    """Lines for the report: per cell set, covariate and stratum, the difference with its p."""
    if not tech:
        return []
    if tech.get("error"):
        return [f"not assessed ({tech['error']})"]
    g0, g1 = tech["groups"]
    out = []
    for r in tech["rows"]:
        head = f"{r['cells']}, {r['covariate']}, {r['stratum']}"
        if r.get("difference") is None:
            out.append(f"{head}: not compared ({r['reason']})")
            continue
        out.append(f"{head}: {g0} {r['first']:.4g}, {g1} {r['second']:.4g}, {g0} − {g1} {r['difference']:+.4g} "
                   f"({r['method']} p {r['p']:.3g}, {r['n_relabelings']} relabelings; the replicates' medians per cell"
                   f"{', over every gene of the pool' if tech['source'] == 'pool' else ''})")
    return out


def describe(comp: dict, top: int = 5) -> str:
    """One line: the TVD with its p and the largest differences."""
    if comp.get("error"):
        return f"{comp['column']}: not assessed ({comp['error']})"
    g0, g1 = comp["groups"]
    lv = "; ".join(f"{r['level']} {r['first']:.3f} -> {r['second']:.3f}" for r in comp["levels"][:top])
    absent = ("; absent from every replicate of one group: "
              + ", ".join(f"{v} (no {g})" for v, g in sorted(comp["absent"].items()))) if comp["absent"] else ""
    return (f"{comp['column']}: total variation distance {comp['tvd']:.3f} between the groups' mean compositions "
            f"({comp['method']} p {comp['tvd_p']:.3g}, {comp['n_relabelings']} relabelings); the largest "
            f"differences from {g0} to {g1}: {lv}{absent}")
