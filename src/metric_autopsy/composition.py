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
