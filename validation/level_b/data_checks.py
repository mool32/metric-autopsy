#!/usr/bin/env python3
"""Level B, step S5 (`plan_s5.md`): the checks of a kept claim's data, independent of the engine.

Each function takes what a dataset's loader gives — the counts (cells x genes, dense or scipy.sparse), the gene
symbols and the cells' metadata (a pandas frame) — and returns plain numbers and words for the list (S7). Nothing
here reads a claim's metric; nothing imports the engine.

1. `counts_check`: the counts are non-negative integers (rule 3: raw counts or the authors' counts).
2. `donors_check`: the donors per group of the claim's comparison, nested or paired, ids that name several donors,
   and rule 4 (`rule4`: the smallest p of an exact two-sided permutation test below 0.05).
3. `sex_check`: where both sexes are among the donors, each donor's sex from the X-inactivation gene and the Y
   genes (`controls.sex_call`), a donor whose cells carry both in many cells, and a call against its label.
4. `genes_check`: the claim's genes in the matrix (by symbol, case ignored).
5. `source_check`: the condition against the data's source (12.3): fails where the groups share no source,
   balanced where every source holds both groups, partly otherwise.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from controls import SEX_GENES, sex_call  # noqa: E402  the same call as the controls'

ALPHA = 0.05
SPLIT_SHARE = 0.20  # a donor whose cells carry the X-inactivation gene and the Y genes each in >= 20% of them


def rule4(design: str, n1: int, n2: int | None = None) -> dict:
    """Rule 4 of the protocol: an exact two-sided permutation test over donors can reach p < 0.05. Nested: the
    smallest p is 2 / C(n1 + n2, n1) (4 against 4: 2/70); paired: 2 / 2^n (6 pairs: 2/64)."""
    if design == "paired":
        relabelings = 2 ** n1
    else:
        relabelings = math.comb(n1 + (n2 or 0), n1)
    p = 2 / relabelings if relabelings else 1.0
    return dict(design=design, relabelings=relabelings, min_p=min(p, 1.0), passes=p < ALPHA)


def counts_check(X) -> dict:
    v = X.data if hasattr(X, "data") and not isinstance(X, np.ndarray) else np.asarray(X)
    v = np.asarray(v)
    return dict(nonnegative=bool(np.all(v >= 0)), integers=bool(np.all(v == np.round(v))),
                raw=bool(np.all(v >= 0) and np.all(v == np.round(v))), max=float(v.max()) if v.size else 0.0)


def donors_check(obs: pd.DataFrame, donor_col: str, group_col: str, groups: tuple) -> dict:
    """The donors of each group; paired where every donor is in both groups, nested where each is in one,
    mixed otherwise; ids naming several donors ('/', '+', ',' or ';' in the id)."""
    o = obs[obs[group_col].astype(str).isin([str(g) for g in groups])]
    donor = o[donor_col].astype(str)
    of = {str(g): sorted(set(donor[o[group_col].astype(str) == str(g)])) for g in groups}
    a, b = (set(of[str(g)]) for g in groups)
    composite = sorted(d for d in a | b if any(c in d for c in "/+,;"))
    if a and a == b:
        design, n1, n2 = "paired", len(a), None
    elif not (a & b):
        design, n1, n2 = "nested", len(a), len(b)
    else:
        design, n1, n2 = "mixed", len(a), len(b)
    cells = o.groupby([group_col, donor_col]).size()
    res = dict(design=design, donors={k: len(v) for k, v in of.items()}, composite=composite,
               cells_per_donor={str(k): int(v) for k, v in cells.items()})
    res["rule4"] = rule4(design, n1, n2) if design != "mixed" else dict(design="mixed", passes=False)
    return res


def sex_check(X, genes: list, obs: pd.DataFrame, donor_col: str, organism: str, label_col: str | None = None) -> dict:
    """Each donor's shares of cells with an X-inactivation count and with a Y-gene count, its call, whether its
    cells carry both in many of them (`SPLIT_SHARE`: two individuals under one id), and its label."""
    xg, yg = SEX_GENES[organism]
    upper = {str(g).upper(): j for j, g in enumerate(genes)}
    jx = upper.get(xg.upper())
    jy = [upper[g.upper()] for g in yg if g.upper() in upper]
    rows = []
    for d, idx in obs.groupby(donor_col).indices.items():
        sub = X[idx]
        fx = float((np.asarray(sub[:, [jx]].sum(axis=1)).ravel() > 0).mean()) if jx is not None else float("nan")
        fy = float((np.asarray(sub[:, jy].sum(axis=1)).ravel() > 0).mean()) if jy else float("nan")
        label = str(obs[label_col].iloc[idx[0]]) if label_col else "-"
        call = sex_call(fx, fy)
        rows.append(dict(donor=str(d), cells=len(idx), xist=fx, y=fy, call=call, label=label,
                         split=bool(fx >= SPLIT_SHARE and fy >= SPLIT_SHARE),
                         against_label=bool(label in ("female", "male") and call in ("female", "male") and label != call)))
    calls = {r["call"] for r in rows}
    both = {"female", "male"} <= calls or {"female", "male"} <= {r["label"] for r in rows}
    return dict(both_sexes=both, donors=rows, fails=[r["donor"] for r in rows if both and (r["split"] or r["against_label"])])


def genes_check(genes: list, claim_genes: list) -> dict:
    have = {str(g).upper() for g in genes}
    missing = [g for g in claim_genes if str(g).upper() not in have]
    return dict(present=[g for g in claim_genes if g not in missing], missing=missing, passes=not missing)


def source_check(obs: pd.DataFrame, group_col: str, groups: tuple, source_col: str) -> dict:
    """The cells' sources (a batch, run, lane, platform or study) per group: fails where no source holds both
    groups, balanced where every source does, partly otherwise; with the shares of each source per group."""
    o = obs[obs[group_col].astype(str).isin([str(g) for g in groups])]
    tab = pd.crosstab(o[source_col].astype(str), o[group_col].astype(str))
    both = [s for s in tab.index if all(tab.loc[s, str(g)] > 0 for g in groups if str(g) in tab.columns)
            and all(str(g) in tab.columns for g in groups)]
    if not both:
        verdict = "fails"
    elif len(both) == len(tab.index):
        verdict = "balanced"
    else:
        verdict = "partly"
    shares = {str(g): {s: float(tab.loc[s, str(g)] / tab[str(g)].sum()) for s in tab.index}
              for g in groups if str(g) in tab.columns}
    return dict(verdict=verdict, sources=list(map(str, tab.index)), shared=both, shares=shares)
