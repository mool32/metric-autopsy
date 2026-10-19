"""The check of specificity: the claim among matched controls (v0.4 change 13; a line of the report,
not a condition of SUPPORTED, decided 2026-10-18).

For a claim about given genes (a pair, or a set), K = 1,000 controls are drawn in the claim's own data
before the claim's metric is read on them, from every gene the data hold (not a cut of the most
expressed), matched to each claim gene on its mean count and its detection in the claim's cells: the
rule of the pilot's second diagnostics (``validation/exploratory/real_pilot_4_controls.py``). A pair's
controls are pairs of one gene among the 100 nearest to each claim gene; a single gene's, genes among
its 2,000 nearest; a set's, sets of one gene among the 100 nearest to each of its genes. The controls
leave out the claim's genes, the genes the user names as known to differ by the claim's factor
(``known_to_differ``) and, for pairs, the pairs of the user's optional table of known interactions
(``known_interactions``). The draw is keyed by the claim (``provenance.claim_seed``).

The statistic is the claim's donor-level t (the claim's estimate over its standard error, within the
domain when the claim has one, on the uncorrected data), the same for every control on the same cells.
The line reports:

* how many of the K controls lie at least as far as the claim, in the claimed direction and by absolute
  value, with p = (1 + k) / (1 + K) and the Clopper-Pearson interval of k / K; where that interval
  crosses alpha the outcome is "undecided" (decided 2026-10-17);
* the nearest controls (the 10 whose t lies closest to the claim's), and every control at least as far;
* the controls' common shift: the median of their t in the claim's direction, with an interval from the
  relabelings of the replicates (the median recomputed under each; a shift model), and a warning when it
  lies wholly on the claim's side of 0: the claim's change may be part of a change shared by many genes.

A global claim (GCL and the like) is not checked: a global shift is unspecific by construction.

The metric is read on each unit's counts of the claim's and the controls' genes plus one column of the
cells' other counts (``REST``), so a metric that reads its genes and the cells' totals sees the data
exactly; a metric that reads other genes is not supported (the line checks the claim's own value).
"""
from __future__ import annotations

from typing import Callable, Sequence

import numpy as np

from .combined import Estimator, Units
from .core import SimpleData, as_dense
from .stats import clopper_pearson

K_CONTROLS = 1000
NEAREST_PAIR, NEAREST_GENE = 100, 2000
N_NEAREST = 10
REST = "__other_counts__"


def _gene_stats(X, rows) -> tuple[np.ndarray, np.ndarray]:
    """The mean count per cell and the share of cells with a count, per gene, over `rows`."""
    sub = X[rows]
    n = sub.shape[0]
    if hasattr(sub, "tocsr"):
        sub = sub.tocsr()
        g = sub.shape[1]
        total = np.bincount(sub.indices, weights=np.asarray(sub.data, dtype=np.float64), minlength=g)
        nz = np.bincount(sub.indices[np.asarray(sub.data) > 0], minlength=g)
    else:
        sub = np.asarray(sub, dtype=np.float64)
        total, nz = sub.sum(axis=0), (sub > 0).sum(axis=0)
    return total / max(n, 1), nz / max(n, 1)


def _distances(mean, det, cand, j0) -> np.ndarray:
    cand = np.asarray(cand)
    lm = np.log10(mean[cand])
    s_m, s_d = float(np.std(lm)), float(np.std(det[cand]))
    s_m, s_d = (s_m if s_m > 0 else 1.0), (s_d if s_d > 0 else 1.0)
    return np.sqrt(((lm - np.log10(mean[j0])) / s_m) ** 2 + ((det[cand] - det[j0]) / s_d) ** 2)


def _nearest(mean, det, cand, j0, k) -> list:
    cand = sorted(cand)
    order = np.argsort(_distances(mean, det, cand, j0), kind="mergesort")
    return [cand[i] for i in order[:k]]


def draw_controls(names, mean, det, claim_genes, kind: str, rng, *, k=K_CONTROLS, known_genes=(),
                  known_pairs=()) -> dict:
    """The controls, by the rule above; `names`, `mean` and `det` over every gene of the pool."""
    names = [str(g) for g in names]
    seen = {}
    for g in names:
        seen[g] = seen.get(g, 0) + 1
    index = {g: j for j, g in enumerate(names)}
    claim_genes = [str(g) for g in claim_genes]
    exclude = set(claim_genes) | {str(g) for g in known_genes}
    cand = [j for j, g in enumerate(names) if mean[j] > 0 and seen[g] == 1 and g not in exclude]
    js = [index[g] for g in claim_genes]
    known = {frozenset((str(a), str(b))) for a, b in known_pairs}
    if kind == "pair":
        A, B = _nearest(mean, det, cand, js[0], NEAREST_PAIR), _nearest(mean, det, cand, js[1], NEAREST_PAIR)
        eligible, keys, n_known = [], set(), 0
        for a in sorted(A):
            for b in sorted(B):
                key = (min(a, b), max(a, b))
                if a == b or key in keys:
                    continue
                keys.add(key)
                if frozenset((names[a], names[b])) in known:
                    n_known += 1
                    continue
                eligible.append((a, b))
        pick = rng.choice(len(eligible), size=min(k, len(eligible)), replace=False) if eligible else []
        controls = [eligible[i] for i in pick]
        lists = (A, B)
        extra = dict(eligible=len(eligible), known_pairs_excluded=n_known)
    elif len(js) == 1:
        near = sorted(_nearest(mean, det, cand, js[0], NEAREST_GENE))
        pick = rng.choice(len(near), size=min(k, len(near)), replace=False) if near else []
        controls = [(near[i],) for i in pick]
        lists = (near,)
        extra = dict(eligible=len(near))
    else:
        lists = [sorted(_nearest(mean, det, cand, j, NEAREST_PAIR)) for j in js]
        controls, keys, tries = [], set(), 0
        while len(controls) < k and tries < 200 * k and all(lists):
            tries += 1
            s = tuple(L[int(rng.integers(len(L)))] for L in lists)
            key = frozenset(s)
            if len(key) < len(s) or key in keys:
                continue
            keys.add(key)
            controls.append(s)
        extra = dict(draws=tries)
    matched = {g: dict(mean=float(mean[j]), detection=float(det[j]),
                       mean_range=[float(mean[L].min()), float(mean[L].max())] if len(L) else None,
                       detection_range=[float(det[L].min()), float(det[L].max())] if len(L) else None)
               for g, j, L in zip(claim_genes, js, [np.asarray(x, int) for x in lists])}
    return dict(controls=[[names[j] for j in c] for c in controls], candidates=len(cand), matched=matched, **extra)


def _units_on_rows(data, group_col, groups, replicate_col, within, domain) -> Units | None:
    """The claim's units with cell indices into `data`'s rows (the pool shares them)."""
    if domain is not None:
        from .domain import build_units
        return build_units(domain, data, group_col=group_col, groups=groups, replicate_col=replicate_col)
    from .effect import _replicate_design
    grp = np.asarray(data.obs[group_col])
    keep = (grp == groups[0]) | (grp == groups[1])
    design = _replicate_design(data.obs[keep], group_col, groups, replicate_col, list(within))
    if design["kind"] == "mixed":
        return None
    rep = np.asarray(data.obs[replicate_col]).astype(str)
    reps = list(design["reps"])
    strata = list(dict.fromkeys(design["strata"][r] for r in reps))
    s_index = {s: i for i, s in enumerate(strata)}
    first = (np.array([1 if groups[0] in design["membership"][r] else 0 for r in reps])
             if design["kind"] == "nested" else np.zeros(len(reps), int))
    cells = []
    for r in reps:
        m = keep & (rep == r)
        cells.append((np.where(m & (grp == groups[0]))[0], np.where(m & (grp == groups[1]))[0])
                     if design["kind"] == "paired" else np.where(m)[0])
    return Units(design["kind"], reps, first, np.array([s_index[design["strata"][r]] for r in reps]), strata,
                 ["all cells"], np.arange(len(reps)), np.zeros(len(reps), int), cells)


def _reduced(pool_X, rows, cols, names, obs) -> SimpleData:
    """The cells `rows` with the columns `cols` and one column of their other counts."""
    sub = pool_X[rows]
    tot = np.asarray(sub.sum(axis=1), dtype=float).ravel()
    sel = np.asarray(as_dense(sub[:, cols]), dtype=float).reshape(len(rows), len(cols))
    rest = np.clip(tot - sel.sum(axis=1), 0.0, None)
    return SimpleData(np.column_stack([sel, rest]), obs.iloc[rows], list(names) + [REST])


def _values(metric_of, pool_X, units, gene_sets, col_of, names_sel, obs) -> np.ndarray:
    """(len(gene_sets), U): the metric of each gene set on each unit's reduced cells."""
    from .gates import _safe_call
    cols = [col_of[g] for g in names_sel]
    V = np.full((len(gene_sets), units.n_units), np.nan)
    for u, c in enumerate(units.cells):
        parts = c if units.kind == "paired" else (c,)
        subs = [_reduced(pool_X, np.asarray(p, int), cols, names_sel, obs) for p in parts]
        for i, gs in enumerate(gene_sets):
            vals = [_safe_call(lambda s=s: metric_of(s, gs)) for s in subs]
            if all(v is not None for v in vals):
                V[i, u] = vals[0] - vals[1] if units.kind == "paired" else vals[0]
    return V


def check(metric_of: Callable, data, *, claim_genes: Sequence[str], kind: str, group_col: str, groups,
          replicate_col: str, within: Sequence[str] = (), domain: dict | None = None, direction: str | None = None,
          alpha: float = 0.05, rng=None, n_perm: int = 1000, max_exact: int = 20000, k: int = K_CONTROLS,
          pool=None, known_genes: Sequence[str] = (), known_pairs: Sequence = (), metric=None) -> dict:
    """The line of the report (see the module's docstring). `metric_of(data, genes)` is the metric as a
    function of its genes; `pool` holds the same cells as `data` with every gene (default `data`)."""
    rng = rng if rng is not None else np.random.default_rng(0)
    if kind == "global" or not claim_genes:
        return dict(applicable=False, reason="a global claim names no genes: a shift shared by all genes is "
                                             "unspecific by construction, so specificity is not checked")
    pool = data if pool is None else pool
    if pool.X.shape[0] != len(data.obs):
        return dict(applicable=True, error=f"the pool holds {pool.X.shape[0]} cells, the data {len(data.obs)}")
    units = _units_on_rows(data, group_col, groups, replicate_col, within, domain)
    if units is None or units.n_units == 0:
        return dict(applicable=True, error="the claim has no units (replicates, or the domain's types)")
    rows = np.unique(np.concatenate([np.concatenate(c) if units.kind == "paired" else c for c in units.cells]))
    names = [str(g) for g in pool.var_names]
    mean, det = _gene_stats(pool.X, rows)
    drawn = draw_controls(names, mean, det, claim_genes, kind, rng, k=k, known_genes=known_genes,
                          known_pairs=known_pairs)
    controls = drawn.pop("controls")
    if not controls:
        return dict(applicable=True, error="no control could be drawn", **drawn)
    claim_genes = [str(g) for g in claim_genes]
    names_sel = list(dict.fromkeys(claim_genes + [g for c in controls for g in c]))
    col_of = {g: j for j, g in enumerate(names)}
    V = _values(metric_of, pool.X, units, [claim_genes] + controls, col_of, names_sel, data.obs)
    tie = None
    if pool is data and metric is not None:
        from .gates import _safe_call
        u0 = units.cells[0][0] if units.kind == "paired" else units.cells[0]
        full = _safe_call(lambda: metric(data[np.asarray(u0, int)]))
        red = _safe_call(lambda: metric_of(_reduced(pool.X, np.asarray(u0, int), [col_of[g] for g in names_sel],
                                                    names_sel, data.obs), claim_genes))
        tie = dict(full=full, reduced=red)
        if full is None or red is None or not np.isclose(full, red, rtol=1e-6, atol=1e-9):
            return dict(applicable=True, error=(f"the metric on a unit's reduced data ({red}) differs from the metric on "
                                                f"its cells ({full}): it reads more than its genes and the cells' totals"),
                        tie=tie)
    est = Estimator(units)
    D, se, df = est.estimate(V)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = D / se
    t_claim, t_c = float(t[0]), t[1:]
    ok = np.isfinite(t_c)
    if not np.isfinite(t_claim) or not ok.any():
        return dict(applicable=True, error="the claim's or the controls' t could not be computed", **drawn)
    if direction in ("increase", "decrease"):
        sign = 1.0 if direction == "decrease" else -1.0  # the estimate is groups[0] - groups[1]
        side = f"in the claimed direction ({direction} from {groups[0]} to {groups[1]})"
    else:
        sign = 1.0 if t_claim >= 0 else -1.0
        side = "in the claim's observed direction (a two-sided claim)"
    ts, tc = sign * t_claim, sign * t_c[ok]
    names_c = [controls[i] for i in np.where(ok)[0]]
    K = int(ok.sum())
    k_dir = int(np.sum(tc >= ts))
    k_abs = int(np.sum(np.abs(tc) >= abs(ts)))

    def outcome(kk):
        lo, hi = clopper_pearson(kk, K, alpha)
        p = (1 + kk) / (1 + K)
        what = "specific" if hi < alpha else "not specific" if lo > alpha else "undecided"
        return dict(k=kk, n=K, share=kk / K, share_ci=(lo, hi), p=p, outcome=what)
    near = np.argsort(np.abs(tc - ts), kind="mergesort")[:N_NEAREST]
    as_far = np.where(tc >= ts)[0]
    shift = _shift(est, V[1:][ok], sign, alpha, n_perm, rng, max_exact)
    return dict(applicable=True, statistic="the donor-level t of the claim's estimate"
                + (" within the domain's types" if domain is not None else ""), data="uncorrected", side=side,
                groups=[str(groups[0]), str(groups[1])], sign=sign,
                t_claim=t_claim, t_claim_signed=ts, k_requested=k, n_controls=K,
                in_direction=outcome(k_dir), absolute=outcome(k_abs),
                nearest=[dict(genes=names_c[i], t=float(sign * tc[i]), t_signed=float(tc[i]),
                              at_least_as_far=bool(tc[i] >= ts)) for i in near],
                at_least_as_far=[dict(genes=names_c[i], t_signed=float(tc[i])) for i in as_far[np.argsort(-tc[as_far])]],
                shift=shift, pool=dict(genes=len(names), same_as_data=pool is data), tie=tie, **drawn)


def _shift(est: Estimator, Vc: np.ndarray, sign: float, alpha: float, n_perm: int, rng, max_exact: int) -> dict:
    """The controls' median t in the claim's direction, with its interval over the relabelings of the
    replicates (a shift model: the observed median minus the relabelings' quantiles) and the one-sided
    p of a shift in the claim's direction."""
    D, se, _ = est.estimate(Vc)
    with np.errstate(invalid="ignore", divide="ignore"):
        m_obs = float(np.nanmedian(sign * D / se))
    labs, method, _ = est.relabelings(n_perm=n_perm, rng=rng, max_exact=max_exact)
    null = []
    for lab in labs:
        d, s, _ = est.estimate(Vc, lab)
        with np.errstate(invalid="ignore", divide="ignore"):
            null.append(np.nanmedian(sign * d / s))
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    if not len(null) or not np.isfinite(m_obs):
        return dict(median_t=m_obs, error="no relabeling gave a median")
    q_lo, q_hi = np.quantile(null, [alpha / 2.0, 1.0 - alpha / 2.0])
    ci = (m_obs - q_hi, m_obs - q_lo)
    count = int(np.sum(null >= m_obs - 1e-12 * max(1.0, abs(m_obs))))
    p = count / len(null) if method == "exact" else (1 + count) / (1 + len(null))
    return dict(median_t=m_obs, ci=ci, ci_level=1.0 - alpha, p_one_sided=float(p), method=method,
                n_relabelings=int(len(null)), significant=bool(ci[0] > 0),
                warning=("the controls shift together in the claim's direction: the claim's change may be part of a "
                         "change shared by many genes") if ci[0] > 0 else None)


def describe(line: dict) -> list:
    """Lines of the report."""
    if not line:
        return []
    if not line.get("applicable", True):
        return [line["reason"]]
    if line.get("error"):
        return [f"not computed: {line['error']}"]
    a, b = line["in_direction"], line["absolute"]
    # the t values in the dossier's convention, groups[0] - groups[1] (the owner's decision of 2026-10-19);
    # the counts stay in the claim's direction, which `side` names
    sign = float(line.get("sign", 1.0))
    g0, g1 = line.get("groups") or ("groups[0]", "groups[1]")
    claim_sign = "positive" if sign > 0 else "negative"
    out = [f"{line['statistic']} ({line['data']} data), {g0} − {g1}: t = {sign * line['t_claim_signed']:+.3g}; "
           f"{line['side']}, a {claim_sign} t",
           f"{a['k']} of {a['n']} matched controls at least as far in the claim's direction (p {a['p']:.3g}, share's "
           f"interval {a['share_ci'][0]:.3g}–{a['share_ci'][1]:.3g}): {a['outcome']}; by absolute value {b['k']} "
           f"(p {b['p']:.3g}): {b['outcome']}",
           "nearest controls: " + "; ".join(f"{'–'.join(c['genes'])} {sign * c['t_signed']:+.3g}"
                                           + (" (at least as far)" if c["at_least_as_far"] else "")
                                           for c in line["nearest"])]
    sh = line.get("shift") or {}
    if "ci" in sh:
        lo, hi = sorted((sign * sh["ci"][0], sign * sh["ci"][1]))
        out.append(f"the controls' common shift: median t {sign * sh['median_t']:+.3g} ({sh['ci_level']:.0%} interval over "
                   f"{sh['n_relabelings']} relabelings [{lo:+.3g}, {hi:+.3g}], one-sided p in the claim's direction "
                   f"{sh['p_one_sided']:.3g})" + (f" — WARNING: {sh['warning']}" if sh.get("warning") else ""))
    return out
