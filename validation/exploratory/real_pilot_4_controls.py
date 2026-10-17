"""The rule that draws the controls of the second diagnostics' check (ii′)
(validation/exploratory/real_pilot_diagnostics_2_plan.md): genes matched to a claim gene on mean and
detection over the full set of genes, K = 1,000 per claim. numpy and pandas only, no engine: the
extraction (real_pilot_4_extract.py, in the selection's environment) draws the controls, and the
analysis (real_pilot_4.py) reads them.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np

SEED = 20261018
K = 1000
NEAREST_PAIR, NEAREST_GENE = 100, 2000
STRING_MIN_SCORE = 400
STRING_FILES = ("10090.protein.info.v12.0.txt.gz", "10090.protein.links.v12.0.txt.gz")
STRING_SHA256 = {"10090.protein.info.v12.0.txt.gz": "ededc5cefddea16013c4d8f7d5a1080efea22d4b3f1f090c27019eb1d5209c08",
                 "10090.protein.links.v12.0.txt.gz": "803775e1d1a47864561e41944bf1f988094a836991e7cf13e4229233fc8cff3d"}
PAIR = ("Smad3", "Col1a1")
Y_GENES = ("Ddx3y", "Eif2s3y", "Kdm5d", "Uty")
S4_CONTROLS = ("Actb", "Gapdh", "Tst", "Lrrc42")
SEX_GENES = ("Xist", "Tsix", "Jpx", "Ftx",                                    # the first diagnostics' list
             "Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Zfy1", "Zfy2", "Sry", "Usp9y", "Ube1y1",
             "Ddx3x", "Eif2s3x", "Kdm5c", "Kdm6a", "Zfx", "Usp9x", "Uba1",
             "Mid1", "Erdr1", "Sts")


def gene_stats(X) -> tuple[np.ndarray, np.ndarray]:
    """The mean count per cell and the share of cells with a count, per gene, of a CSR or dense matrix."""
    n = X.shape[0]
    if hasattr(X, "indices"):
        X = X.tocsr()
        g = X.shape[1]
        total = np.bincount(X.indices, weights=X.data.astype(np.float64), minlength=g)
        nz = np.bincount(X.indices[X.data > 0], minlength=g)
    else:
        total = np.asarray(X, dtype=np.float64).sum(axis=0)
        nz = (np.asarray(X) > 0).sum(axis=0)
    return total / n, nz / n


def pool(symbols, mean, det, exclude, mapped=None) -> list:
    """The candidates: a count in the claim's cells, a symbol that names one feature, not excluded,
    and (for pairs) a STRING preferred name."""
    symbols = list(symbols)
    seen = {}
    for s in symbols:
        seen[s] = seen.get(s, 0) + 1
    ex = set(exclude)
    return [j for j, s in enumerate(symbols)
            if mean[j] > 0 and seen[s] == 1 and s not in ex and (mapped is None or s in mapped)]


def distances(mean, det, cand, j0) -> np.ndarray:
    """sqrt(((log10 m - log10 m0) / s_m)^2 + ((d - d0) / s_d)^2) from gene j0 to every candidate, the
    SDs over the candidates."""
    cand = np.asarray(cand)
    lm = np.log10(mean[cand])
    s_m, s_d = float(np.std(lm)), float(np.std(det[cand]))
    s_m, s_d = (s_m if s_m > 0 else 1.0), (s_d if s_d > 0 else 1.0)
    return np.sqrt(((lm - np.log10(mean[j0])) / s_m) ** 2 + ((det[cand] - det[j0]) / s_d) ** 2)


def nearest(mean, det, cand, j0, k) -> list:
    """The k candidates nearest to gene j0, ties by the features' order."""
    cand = sorted(cand)
    order = np.argsort(distances(mean, det, cand, j0), kind="mergesort")
    return [cand[i] for i in order[:k]]


def string_names(info_path: Path) -> dict:
    """{preferred_name: [protein ids]} of STRING's protein.info."""
    out = {}
    with gzip.open(info_path, "rt", encoding="utf-8") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        if head[:2] != ["#string_protein_id", "preferred_name"]:
            raise ValueError(f"{info_path.name}: unexpected header {head[:2]}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out.setdefault(f[1], []).append(f[0])
    return out


def string_links(links_path: Path, proteins: set, min_score: int = STRING_MIN_SCORE) -> set:
    """The links of combined score >= min_score among `proteins`."""
    out = set()
    with gzip.open(links_path, "rt", encoding="utf-8") as fh:
        head = fh.readline().split()
        if head != ["protein1", "protein2", "combined_score"]:
            raise ValueError(f"{links_path.name}: unexpected header {head}")
        for line in fh:
            p1, p2, s = line.split()
            if p1 in proteins and p2 in proteins and int(s) >= min_score:
                out.add(frozenset((p1, p2)))
    return out


def draw_pairs(symbols, mean, det, claim_genes, names: dict, links_path: Path, k=K, nearest_k=NEAREST_PAIR) -> dict:
    """1a-males: one gene from each claim gene's nearest, two genes, no STRING link >= 400 between
    their proteins, each unordered pair once; k drawn without replacement."""
    symbols = list(symbols)
    ia, ib = (symbols.index(g) for g in claim_genes)
    cand = pool(symbols, mean, det, set(claim_genes) | set(S4_CONTROLS) | set(SEX_GENES), mapped=names)
    A, B = nearest(mean, det, cand, ia, nearest_k), nearest(mean, det, cand, ib, nearest_k)
    linked = string_links(links_path, {p for j in set(A) | set(B) for p in names[symbols[j]]})
    eligible, seen, n_linked = [], set(), 0
    for a in sorted(A):
        for b in sorted(B):
            key = (min(a, b), max(a, b))
            if a == b or key in seen:
                continue
            seen.add(key)
            if any(frozenset((p, q)) in linked for p in names[symbols[a]] for q in names[symbols[b]]):
                n_linked += 1
                continue
            eligible.append((a, b))
    rng = np.random.default_rng(SEED)
    pick = rng.choice(len(eligible), size=min(k, len(eligible)), replace=False) if eligible else []
    pairs = [eligible[i] for i in pick]
    return dict(controls=[[symbols[a], symbols[b]] for a, b in pairs], columns=[[int(a), int(b)] for a, b in pairs],
                eligible=len(eligible), linked_excluded=n_linked, candidates=len(cand),
                nearest={claim_genes[0]: [symbols[j] for j in A], claim_genes[1]: [symbols[j] for j in B]},
                matched=_matched(mean, det, claim_genes, symbols, (A, B)))


def draw_genes(symbols, mean, det, gene, k=K, nearest_k=NEAREST_GENE) -> dict:
    """2a: k genes drawn without replacement from the nearest_k candidates nearest to the claim gene."""
    symbols = list(symbols)
    cand = pool(symbols, mean, det, {gene} | set(SEX_GENES))
    near = sorted(nearest(mean, det, cand, symbols.index(gene), nearest_k))
    rng = np.random.default_rng(SEED)
    pick = rng.choice(len(near), size=min(k, len(near)), replace=False)
    return dict(controls=[[symbols[near[i]]] for i in pick], columns=[[int(near[i])] for i in pick],
                candidates=len(cand), matched=_matched(mean, det, [gene], symbols, (near,)))


def draw_sets(symbols, mean, det, claim_genes, k=K, nearest_k=NEAREST_PAIR, max_tries=1_000_000) -> dict:
    """2b: sets of one gene from each claim gene's nearest candidates, all different, each set once."""
    symbols = list(symbols)
    cand = pool(symbols, mean, det, set(claim_genes) | set(SEX_GENES))
    lists = [sorted(nearest(mean, det, cand, symbols.index(g), nearest_k)) for g in claim_genes]
    rng = np.random.default_rng(SEED)
    sets, seen, tries = [], set(), 0
    while len(sets) < k and tries < max_tries:
        tries += 1
        s = [L[int(rng.integers(len(L)))] for L in lists]
        key = frozenset(s)
        if len(key) < len(s) or key in seen:
            continue
        seen.add(key)
        sets.append(s)
    return dict(controls=[[symbols[j] for j in s] for s in sets], columns=[[int(j) for j in s] for s in sets],
                candidates=len(cand), draws=tries, matched=_matched(mean, det, claim_genes, symbols, lists))


def _matched(mean, det, claim_genes, symbols, lists) -> dict:
    """The claim genes' mean and detection against the range of their matched candidates."""
    out = {}
    for g, L in zip(claim_genes, lists):
        j = symbols.index(g)
        out[g] = dict(mean=float(mean[j]), detection=float(det[j]),
                      mean_range=[float(mean[L].min()), float(mean[L].max())],
                      detection_range=[float(det[L].min()), float(det[L].max())])
    return out
