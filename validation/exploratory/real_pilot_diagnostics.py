#!/usr/bin/env python3
"""Diagnostics on the data of the pilot's second run (validation/exploratory/real_pilot_diagnostics_plan.md).

Exploratory, on development data only. Four checks, each asked whether it separates the drop of MI on
preprint section 4 from the sex differences, on the claims 1a (mi_3bin of Smad3 and Col1a1, young ->
old within sex), 1a-males (the same in males), 1b (norm_pearson of the pair), 2a (Xist, female ->
male) and 2b (the four Y genes):

- (i) the engine of the tag v0.3.0-prereg with GATE 4 and GATE 5 off: the raw and the corrected effect;
- (ii) specificity: the claim's donor-level difference among 100 controls of similar detection without
  known interactions (STRING v12.0 for the pairs), chosen by the plan's rule;
- (iii) plausibility: per mouse, the metric's response to the maximal dose of the construct and to its
  full removal, against the observed difference, by stratum;
- (iv) agreement of measures: the sign of the change by two measures of the same construct, by stratum.

    python validation/exploratory/real_pilot_diagnostics.py --data DATA --string STRING --out OUT --part case1
    python validation/exploratory/real_pilot_diagnostics.py --merge PART_DIR... --out OUT

DATA holds the second run's select.json, case1.npz and case2.npz (their sha256 are checked against
select.json); STRING holds 10090.protein.info.v12.0.txt.gz and 10090.protein.links.v12.0.txt.gz. A part
writes OUT/<part>.json, OUT/<part>.log and OUT/reports/<claim>.json (the engine's reports of (i)), after
every step; the merge writes OUT/real_pilot_3.json and OUT/real_pilot_3.log with the outcomes by check
and claim, whether each check separates, and the owner's prediction scored. --fabricate-string writes
STRING-like files for the genes of a fabricated case (the dry run).

(ii) to (iv) compute a metric on the claim's genes and one more column, the rest of the cell's counts:
the metrics read only the named genes and the row sum, so the values are the full matrix's exactly
(checked on the claim per mouse, and against the second run's (b) where its JSON is given).
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
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

import metric_autopsy as MA  # noqa: E402  (the engine: the tag's, installed by the workflow)
from metric_autopsy import SimpleData, injected_signal, metrics, run_autopsy, stats  # noqa: E402
from anchors import log_cp10k_mean  # noqa: E402  (R1's metric)
from real_pilot import design_floor, label_of  # noqa: E402  (the second run's floor and verdict label)

ALPHA = 0.05
SEED = 20261017
N_CONTROLS = 100
NEAREST_PAIR, NEAREST_GENE = 50, 200
STRING_MIN_SCORE = 400
STRING_FILES = ("10090.protein.info.v12.0.txt.gz", "10090.protein.links.v12.0.txt.gz")
STRENGTHS = (1.0, 3.0, 10.0, 50.0)   # (iii): coupling, the largest response over these
FOLDS = (2.0, 10.0, 100.0)           # (iii): the module probe in every cell, the largest response
K_DRAWS = 20                         # (iii): injections and shams per dose and mouse; permutations
DONOR_PERMS, DONOR_MAX_EXACT = 9999, 20000
REST = "__rest_of_counts__"
PAIR = ("Smad3", "Col1a1")
Y_GENES = ("Ddx3y", "Eif2s3y", "Kdm5d", "Uty")
S4_CONTROLS = ("Actb", "Gapdh", "Tst", "Lrrc42")
SEX_GENES = ("Xist", "Tsix", "Jpx", "Ftx",                                    # the X-inactivation centre
             "Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Zfy1", "Zfy2", "Sry", "Usp9y", "Ube1y1",  # Y genes
             "Ddx3x", "Eif2s3x", "Kdm5c", "Kdm6a", "Zfx", "Usp9x", "Uba1",     # their X homologues
             "Mid1", "Erdr1", "Sts")                                           # pseudoautosomal
CLAIMS = (
    dict(id="1a", case="case1", metric="mi_3bin", genes=list(PAIR), direction="decrease", males=False,
         sesoi="floor", sesoi_scale="observed", other="norm_pearson"),
    dict(id="1a-males", case="case1", metric="mi_3bin", genes=list(PAIR), direction="decrease", males=True,
         sesoi="floor", sesoi_scale="observed", other="norm_pearson"),
    dict(id="1b", case="case1", metric="norm_pearson", genes=list(PAIR), direction="decrease", males=False,
         sesoi=0.1, other="mi_3bin"),
    dict(id="2a", case="case2", metric="log_cp10k_mean", genes=["Xist"], direction="decrease", males=False,
         sesoi=0.5, other="detection"),
    dict(id="2b", case="case2", metric="log_cp10k_mean", genes=list(Y_GENES), direction="increase", males=False,
         sesoi=0.5, other="detection"),
)
SECTION4, SEX = ("1a", "1a-males"), ("2a", "2b")
CHECKS = ("i", "ii", "iii", "iv")
FLAGS, PASSES, NA = "flags", "does not flag", "not applicable"


# --------------------------------------------------------------------------- #
# the data: a case, a claim's cells and mice, reduced matrices
# --------------------------------------------------------------------------- #
def load_case(data_dir: Path, case: dict):
    """The case's counts (int32), genes and obs, its sha256 checked against select.json."""
    path = data_dir / case["data"]["file"]
    got = sha256_file(path)
    if got != case["data"]["sha256"]:
        raise ValueError(f"{path.name}: sha256 {got}, select.json has {case['data']['sha256']}")
    z = np.load(path, allow_pickle=False)
    X = np.asarray(z["X"], dtype=np.int32)
    genes = [str(g) for g in z["genes"]]
    if len(set(genes)) != len(genes):
        raise ValueError(f"{path.name}: gene names are not unique")
    obs = pd.DataFrame({k[4:]: z[k] for k in z.files if k.startswith("obs_")})
    return X, genes, obs, got


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class Design:
    """A claim's cells, mice, groups and strata (the second run's (b)): g is 1 for the first group."""

    def __init__(self, obs: pd.DataFrame, case: dict, claim: dict):
        self.group_col, self.groups, self.rep = case["group_col"], tuple(case["groups"]), case["replicate"]
        self.within = [] if claim["males"] else list(case["within"])
        cells = obs[self.group_col].astype(str).isin(self.groups).to_numpy()
        if claim["males"]:
            cells &= (obs["sex"].astype(str) == "male").to_numpy()
        self.cells = np.where(cells)[0]
        o = obs.iloc[self.cells]
        donor = o[self.rep].astype(str).to_numpy()
        group = o[self.group_col].astype(str).to_numpy()
        strata = (o[self.within].astype(str).agg("|".join, axis=1).to_numpy() if self.within
                  else np.array([""] * len(o)))
        self.mice = sorted(set(donor))
        member = {m: set(group[donor == m]) for m in self.mice}
        if any(len(v) != 1 for v in member.values()):
            raise ValueError("a mouse in both groups: the diagnostics take nested designs only")
        st_of = {m: set(strata[donor == m]) for m in self.mice}
        if any(len(v) != 1 for v in st_of.values()):
            raise ValueError("a mouse in two strata")
        self.labels = {s: i for i, s in enumerate(sorted({next(iter(v)) for v in st_of.values()}))}
        self.g = np.array([1 if self.groups[0] in member[m] else 0 for m in self.mice])
        self.s = np.array([self.labels[next(iter(st_of[m]))] for m in self.mice])
        self.rows = {m: self.cells[donor == m] for m in self.mice}  # row indices into the case's matrix

    def delta(self, values: dict) -> tuple[float, object]:
        """Delta (second - first): over several strata the engine's pooled-stratified estimate."""
        y = np.array([values[m] for m in self.mice], float)
        te = stats.nested_t(y, self.g, self.s)
        return -float(te.est), te

    def by_stratum(self, values: dict, rng=None) -> dict:
        """Per stratum: group means, Delta, the exact permutation within it, the standardized size."""
        y = np.array([values[m] for m in self.mice], float)
        out = {}
        for st, i in self.labels.items():
            m = self.s == i
            a, b = y[m & (self.g == 1)], y[m & (self.g == 0)]
            if not len(a) or not len(b):
                continue
            pr = stats.permutation_test_nested(y[m], self.g[m], np.zeros(int(m.sum()), int), n_perm=DONOR_PERMS,
                                               rng=rng if rng is not None else np.random.default_rng(SEED),
                                               max_exact=DONOR_MAX_EXACT)
            dfree = len(a) + len(b) - 2
            sd = (math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / dfree)
                  if dfree > 0 and len(a) > 1 and len(b) > 1 else float("nan"))
            d = float(b.mean() - a.mean())
            out[st or "all"] = dict(n={self.groups[0]: int(len(a)), self.groups[1]: int(len(b))},
                                    mean={self.groups[0]: float(a.mean()), self.groups[1]: float(b.mean())},
                                    delta=d, p=float(pr.p), method=pr.method, n_relabelings=int(pr.n),
                                    detected=bool(pr.p < ALPHA), sd_pooled=sd,
                                    standardized=(d / sd if np.isfinite(sd) and sd > 0 else float("nan")))
        return out

    def strata_of(self) -> dict:
        """{stratum: (index mask over mice)}."""
        return {st or "all": self.s == i for st, i in self.labels.items()}


def reduced(X: np.ndarray, genes: list, names, tot: np.ndarray, rows: np.ndarray | None = None) -> np.ndarray:
    """The named genes' counts and the rest of each cell's counts (the row sum less the named), as
    float64: the metrics read nothing else, so their values are the full matrix's."""
    cols = [genes.index(g) for g in names]
    sub = X[:, cols] if rows is None else X[np.ix_(rows, cols)]
    t = tot if rows is None else tot[rows]
    rest = t - sub.sum(axis=1, dtype=np.int64)
    return np.column_stack([sub.astype(np.float64), rest.astype(np.float64)])


def as_data(M: np.ndarray, names) -> SimpleData:
    return SimpleData(M, pd.DataFrame(index=np.arange(M.shape[0])), [*names, REST])


def measure_fn(kind: str, names):
    """The measure on a reduced matrix's SimpleData (genes first, then REST)."""
    names = list(names)
    if kind == "mi_3bin":
        return partial(metrics.mi_3bin, gene_a=names[0], gene_b=names[1])
    if kind == "norm_pearson":
        return partial(metrics.norm_pearson, gene_a=names[0], gene_b=names[1])
    if kind == "log_cp10k_mean":
        return partial(log_cp10k_mean, cols=list(range(len(names))))
    if kind == "detection":
        return lambda d: float(np.mean(np.asarray(d.X)[:, :len(names)].sum(axis=1) > 0))
    raise ValueError(kind)


def per_mouse(kind: str, X, genes, names, tot, design: Design) -> dict:
    fn = measure_fn(kind, names)
    return {m: float(fn(as_data(reduced(X, genes, names, tot, r), names))) for m, r in design.rows.items()}


def full_metric(claim: dict, genes: list):
    """The claim's metric on the full matrix, as the engine takes it."""
    g = claim["genes"]
    if claim["metric"] in ("mi_3bin", "norm_pearson"):
        return partial(getattr(metrics, claim["metric"]), gene_a=g[0], gene_b=g[1])
    return partial(log_cp10k_mean, cols=[genes.index(x) for x in g])


# --------------------------------------------------------------------------- #
# STRING v12.0 and the controls of (ii)
# --------------------------------------------------------------------------- #
def string_names(info_path: Path) -> dict:
    """{preferred_name: [protein ids]} of protein.info (tab-separated, a header line)."""
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
    """The links of combined score >= min_score among `proteins` (space-separated, a header line)."""
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


def string_stats(info_path: Path, links_path: Path, probe=("Smad3", "Smad4")) -> dict:
    """What the files hold (the dry run's check of the real download): proteins, names, links >= 400, and
    the link of two named genes."""
    names = string_names(info_path)
    n, n400 = 0, 0
    with gzip.open(links_path, "rt", encoding="utf-8") as fh:
        fh.readline()
        for line in fh:
            n += 1
            n400 += int(line.rsplit(None, 1)[1]) >= STRING_MIN_SCORE
    want = {p for g in probe for p in names.get(g, [])}
    link = string_links(links_path, want) if len(want) >= 2 else set()
    return dict(proteins=sum(len(v) for v in names.values()), names=len(names), links=n, links_400=n400,
                probe=list(probe), probe_mapped=[g for g in probe if g in names], probe_linked=bool(link),
                sha256={p.name: sha256_file(p) for p in (info_path, links_path)})


def nearest(det: np.ndarray, j0: int, cand: list, k: int) -> list:
    """The k candidates nearest to gene j0 by |detection difference|, ties by the file's order."""
    cand = sorted(cand)
    d = np.abs(det[cand] - det[j0])
    return [cand[i] for i in np.argsort(d, kind="mergesort")[:k]]


def candidates(genes: list, exclude, mapped=None) -> list:
    ex = set(exclude)
    return [j for j, g in enumerate(genes) if g not in ex and (mapped is None or g in mapped)]


def control_pairs(det, genes, claim_genes, names: dict, links_path: Path, n=N_CONTROLS, k=NEAREST_PAIR) -> dict:
    """1a, 1a-males, 1b: one gene of each claim gene's k nearest STRING-mapped candidates, two genes,
    no STRING link >= 400 between them, each unordered pair once; n drawn without replacement."""
    ia, ib = (genes.index(g) for g in claim_genes)
    cand = candidates(genes, set(claim_genes) | set(S4_CONTROLS) | set(SEX_GENES), mapped=names)
    A, B = nearest(det, ia, cand, k), nearest(det, ib, cand, k)
    prot = {p for j in set(A) | set(B) for p in names[genes[j]]}
    linked = string_links(links_path, prot)
    eligible, seen, n_linked = [], set(), 0
    for a in sorted(A):
        for b in sorted(B):
            key = (min(a, b), max(a, b))
            if a == b or key in seen:
                continue
            seen.add(key)
            if any(frozenset((p, q)) in linked for p in names[genes[a]] for q in names[genes[b]]):
                n_linked += 1
                continue
            eligible.append((a, b))
    rng = np.random.default_rng(SEED)
    pick = rng.choice(len(eligible), size=min(n, len(eligible)), replace=False) if eligible else []
    pairs = [eligible[i] for i in pick]
    return dict(controls=[[genes[a], genes[b]] for a, b in pairs], eligible=len(eligible), linked_excluded=n_linked,
                controls_detection=[[float(det[a]), float(det[b])] for a, b in pairs],
                candidates=len(cand),
                nearest={claim_genes[0]: [genes[j] for j in A], claim_genes[1]: [genes[j] for j in B]},
                detection_range={claim_genes[0]: [float(det[A].min()), float(det[A].max())],
                                 claim_genes[1]: [float(det[B].min()), float(det[B].max())]},
                claim_detection={g: float(det[genes.index(g)]) for g in claim_genes},
                string_links_among_nearest=len(linked))


def control_genes(det, genes, gene: str, n=N_CONTROLS, k=NEAREST_GENE) -> dict:
    """2a: n genes drawn without replacement from the k candidates nearest to the claim gene."""
    cand = candidates(genes, {gene} | set(SEX_GENES))
    near = sorted(nearest(det, genes.index(gene), cand, k))
    rng = np.random.default_rng(SEED)
    pick = rng.choice(len(near), size=min(n, len(near)), replace=False)
    return dict(controls=[[genes[near[i]]] for i in pick], candidates=len(cand),
                controls_detection=[[float(det[near[i]])] for i in pick],
                detection_range=[float(det[near].min()), float(det[near].max())],
                claim_detection={gene: float(det[genes.index(gene)])})


def control_sets(det, genes, claim_genes, n=N_CONTROLS, k=NEAREST_PAIR, max_tries=200000) -> dict:
    """2b: sets of one gene from each claim gene's k nearest candidates, all different, each set once."""
    cand = candidates(genes, set(claim_genes) | set(SEX_GENES))
    lists = [sorted(nearest(det, genes.index(g), cand, k)) for g in claim_genes]
    rng = np.random.default_rng(SEED)
    sets, seen, tries = [], set(), 0
    while len(sets) < n and tries < max_tries:
        tries += 1
        s = [L[int(rng.integers(len(L)))] for L in lists]
        key = frozenset(s)
        if len(key) < len(s) or key in seen:
            continue
        seen.add(key)
        sets.append(s)
    return dict(controls=[[genes[j] for j in s] for s in sets], candidates=len(cand), draws=tries,
                controls_detection=[[float(det[j]) for j in s] for s in sets],
                detection_range={g: [float(det[L].min()), float(det[L].max())] for g, L in zip(claim_genes, lists)},
                claim_detection={g: float(det[genes.index(g)]) for g in claim_genes})


def position(delta: float, controls, direction: str) -> dict:
    """The claim among the controls: one-sided in the claimed direction, and two-sided."""
    c = np.asarray([x for x in controls if np.isfinite(x)], float)
    tol = 1e-12 * max(1.0, abs(delta))
    if direction == "decrease":
        far = int(np.sum(c <= delta + tol))
    else:
        far = int(np.sum(c >= delta - tol))
    both = int(np.sum(np.abs(c) >= abs(delta) - tol))
    p = (1 + far) / (1 + len(c))
    q = np.quantile(c, [0.05, 0.25, 0.5, 0.75, 0.95]).tolist() if len(c) else []
    return dict(n=int(len(c)), at_least_as_far=far, p_one_sided=float(p), p_two_sided=float((1 + both) / (1 + len(c))),
                quantiles=dict(zip(("5%", "25%", "50%", "75%", "95%"), q)),
                outcome=FLAGS if p >= ALPHA else PASSES)


# --------------------------------------------------------------------------- #
# the checks
# --------------------------------------------------------------------------- #
def outcome_i(det: dict) -> tuple[str, str]:
    """(i)'s outcome from the engine's effect detail."""
    if not det or det.get("raw_detected") is None or det.get("effect") is None:
        return NA, "the effect was not estimated"
    if not det.get("raw_detected"):
        return NA, f"the raw difference is not detected (p_raw {det.get('p_raw'):.4g})"
    raw, est = det.get("raw_effect"), det.get("effect")
    if det.get("explained_by_depth"):
        return FLAGS, "explained by depth (the engine's flag)"
    if det.get("reversed_by_correction"):
        return FLAGS, "the correction reversed the sign"
    detected = det.get("status") == "DETECTED"
    if detected and np.sign(raw) == np.sign(est):
        return PASSES, "the corrected effect is detected with the raw's sign"
    return FLAGS, "the raw difference is detected, the corrected effect is not"


def check_i(claim, case, X, genes, obs, design: Design, floor_te, out_dir: Path) -> dict:
    """The engine with GATE 4 and GATE 5 off (no injected signal, no controls), the effect always estimated."""
    sesoi = claim["sesoi"]
    if sesoi == "floor":
        sesoi = design_floor(floor_te)
        if not np.isfinite(sesoi):
            raise ValueError("no design floor: no finite SE")
    cols = list(dict.fromkeys([design.group_col, design.rep, *design.within]))
    o = obs.iloc[design.cells][cols].astype(str).reset_index(drop=True)
    data = SimpleData(X[design.cells].astype(np.float64), o, list(genes))
    prereg = dict(estimand="composition", direction=claim["direction"], sesoi=float(sesoi), judgment_pending=False)
    if claim.get("sesoi_scale"):
        prereg["sesoi_scale"] = claim["sesoi_scale"]
    pair = tuple(claim["genes"]) if claim["metric"] in ("mi_3bin", "norm_pearson") else None
    a = run_autopsy(full_metric(claim, genes), data, group_col=design.group_col, groups=design.groups, gene_pair=pair,
                    within=design.within, prereg=prereg, replicate_col=design.rep, signal_test=None,
                    stop_on_first_fail=False, seed=0, log_path="off")
    a.metric_name = claim["metric"]
    rep = a.to_dict()
    try:
        text = json.dumps(rep, indent=1, allow_nan=False)
    except ValueError:
        text = json.dumps(rep, indent=1, default=str)
    (out_dir / "reports").mkdir(parents=True, exist_ok=True)
    (out_dir / "reports" / f"{claim['id']}.json").write_text(text)
    eff = rep["fields"].get("effect") or {}
    det = dict(eff.get("detail") or {}, status=eff.get("status"))
    outcome, why = outcome_i(det)
    keep = ("effect", "raw_effect", "ci", "p_perm", "p_raw", "raw_detected", "retained", "explained_by_depth",
            "reversed_by_correction", "tost_equivalent", "sesoi", "effective_sesoi", "test", "n_replicates", "status")
    eq = ((rep["fields"].get("design_adequacy") or {}).get("detail") or {}).get("equalization") or {}
    del data
    return dict(verdict=a.verdict, label=label_of(a.verdict), cause=a.cause, prereg=prereg,
                fields={k: (v and dict(status=v["status"], reason=v["reason"])) for k, v in rep["fields"].items()},
                gates=[dict(gate=g["gate"], status=g["status"], message=g["message"]) for g in rep["gates"]],
                effect={k: det.get(k) for k in keep},
                equalization=[dict(stratum="|".join(map(str, r["stratum"].values())) if isinstance(r.get("stratum"), dict)
                                   else str(r.get("stratum")), thinned_group=r.get("thinned_group"),
                                   depth_ratio=r.get("depth_ratio")) for r in eq.get("strata", [])],
                report=f"reports/{claim['id']}.json", report_sha256=hashlib.sha256(text.encode()).hexdigest(),
                outcome=outcome, why=why)


def responses(kind: str, names, M: np.ndarray, rng_key) -> dict:
    """(iii) on one mouse's reduced matrix: the observed value, the mean response (injected - sham) at
    every dose, and the mean change under the full removal (decoupling, or loss)."""
    fn = measure_fn(kind, names)
    obs_v = float(fn(as_data(M, names)))
    pair = kind in ("mi_3bin", "norm_pearson")
    doses = {}
    for di, dose in enumerate(STRENGTHS if pair else FOLDS):
        inj = (injected_signal.coupling(names[0], names[1], strength=dose) if pair
               else injected_signal.module(list(names), fold=dose, frac=1.0))
        rng = np.random.default_rng([*rng_key, 1, di])
        d = as_data(M, names)
        doses[str(dose)] = float(np.mean([fn(inj(d, rng)) - fn(inj.sham(d, rng)) for _ in range(K_DRAWS)]))
    if pair:
        rng = np.random.default_rng([*rng_key, 2])
        vals = []
        for _ in range(K_DRAWS):
            P_ = M.copy()
            P_[:, 1] = P_[rng.permutation(len(P_)), 1]  # the second gene's counts permuted across the cells
            vals.append(fn(as_data(P_, names)) - obs_v)
        removal = float(np.mean(vals))
    else:
        Z = M.copy()
        Z[:, :len(names)] = 0.0  # the genes' counts lost
        removal = float(fn(as_data(Z, names)) - obs_v)
    return dict(value=obs_v, dose=doses, removal=removal)


def check_iii(claim, design: Design, per: dict) -> dict:
    """Per stratum, h the higher group and l the lower by the claim's metric: (a) the removal lowers h by
    at least |Delta|; (b) the maximal dose (the largest group-mean response over the doses) raises l by
    at least |Delta|."""
    strata = {}
    fails = []
    for st, m in design.strata_of().items():
        mice = [x for x, keep in zip(design.mice, m) if keep]
        g = {x: gi for x, gi in zip(design.mice, design.g)}
        first = [x for x in mice if g[x] == 1]
        second = [x for x in mice if g[x] == 0]
        if not first or not second:
            continue
        mean = {design.groups[0]: float(np.mean([per[x]["value"] for x in first])),
                design.groups[1]: float(np.mean([per[x]["value"] for x in second]))}
        h_name, l_name = ((design.groups[0], design.groups[1]) if mean[design.groups[0]] >= mean[design.groups[1]]
                          else (design.groups[1], design.groups[0]))
        h = first if h_name == design.groups[0] else second
        lo = second if h is first else first
        delta = abs(mean[design.groups[1]] - mean[design.groups[0]])
        removal = float(np.mean([per[x]["removal"] for x in h]))
        by_dose = {k: float(np.mean([per[x]["dose"][k] for x in lo])) for k in per[lo[0]]["dose"]}
        best = max(by_dose, key=lambda k: by_dose[k])
        ok_a, ok_b = -removal >= delta, by_dose[best] >= delta
        if not ok_a:
            fails.append(f"{st}: (a)")
        if not ok_b:
            fails.append(f"{st}: (b)")
        strata[st] = dict(higher=h_name, lower=l_name, mean=mean, abs_delta=delta, removal_higher=removal,
                          dose_lower=by_dose, max_dose=best, max_response_lower=by_dose[best],
                          a_holds=bool(ok_a), b_holds=bool(ok_b),
                          ratio_a=(delta / -removal if removal < 0 else float("inf")),
                          ratio_b=(delta / by_dose[best] if by_dose[best] > 0 else float("inf")))
    return dict(strata=strata, failed=fails, outcome=FLAGS if fails else PASSES)


def check_iv(claim, design: Design, own: dict, other: dict) -> dict:
    """Per stratum, the two measures' Delta, exact p and standardized size: flags where the signs differ
    and at least one is detected."""
    rng = np.random.default_rng([SEED, 4])
    a, b = design.by_stratum(own, rng), design.by_stratum(other, rng)
    strata, fails = {}, []
    for st in a:
        x, y = a[st], b[st]
        opposite = bool(np.sign(x["delta"]) * np.sign(y["delta"]) < 0)
        flag = opposite and (x["detected"] or y["detected"])
        if flag:
            fails.append(st)
        strata[st] = dict(own=x, other=y, opposite=opposite, flagged=bool(flag))
    return dict(measures=[claim["metric"], claim["other"]], strata=strata, failed=fails,
                outcome=FLAGS if fails else PASSES)


# --------------------------------------------------------------------------- #
# a part: one case's claims, every check
# --------------------------------------------------------------------------- #
def run_part(part: str, spec: dict, data_dir: Path, string_dir: Path, out_dir: Path, log, save, run2=None) -> dict:
    case = spec["cases"][part]
    if case.get("status") != "run":
        raise ValueError(f"{part} did not run in the second run")
    t0 = time.time()
    X, genes, obs, sha = load_case(data_dir, case)
    tot = X.sum(axis=1, dtype=np.int64)
    det = (X > 0).mean(axis=0)
    res = dict(data=dict(file=case["data"]["file"], sha256=sha, cells=int(X.shape[0]), genes=int(X.shape[1])),
               claims={}, controls={})
    log(f"## {part}: {X.shape[0]} cells x {X.shape[1]} genes, {case['data']['file']} sha256 {sha} (as select.json); "
        f"groups {case['groups'][0]!r} -> {case['groups'][1]!r}, unit {case['replicate']} ({time.time() - t0:.1f} s)")
    claims = [c for c in CLAIMS if c["case"] == part]
    designs = {c["id"]: Design(obs, case, c) for c in claims}
    for c in claims:
        d = designs[c["id"]]
        log(f"{c['id']}: {len(d.cells)} cells, mice " + "; ".join(
            f"{st or 'all'} {int(((d.s == i) & (d.g == 1)).sum())} -> {int(((d.s == i) & (d.g == 0)).sum())}"
            for st, i in d.labels.items()))

    # the claims' own per-mouse values, checked against the full matrix (and the second run's (b))
    values = {}
    for c in claims:
        d = designs[c["id"]]
        for kind in dict.fromkeys([c["metric"], c["other"]]):
            values[(c["id"], kind)] = per_mouse(kind, X, genes, c["genes"], tot, d)
        fm = full_metric(c, genes)
        diffs = [abs(float(fm(SimpleData(X[r].astype(np.float64), pd.DataFrame(index=np.arange(len(r))), genes)))
                     - values[(c["id"], c["metric"])][m]) for m, r in d.rows.items()]
        rec = dict(id=c["id"], metric=c["metric"], genes=c["genes"], direction=c["direction"], other=c["other"],
                   males=c["males"], within=d.within, values=values[(c["id"], c["metric"])],
                   other_values=values[(c["id"], c["other"])], reduced_vs_full=float(max(diffs)))
        delta, te = d.delta(rec["values"])
        rec.update(delta=delta, se=float(te.se), df=float(te.df), floor=design_floor(te))
        if run2 is not None:
            src = {"1a-males": "1a"}.get(c["id"], c["id"])
            theirs = next((x["b"]["values"] for x in run2["cases"].get(part, {}).get("claims", []) if x["id"] == src), None)
            if theirs is not None:
                common = sorted(set(theirs) & set(rec["values"]))
                rec["vs_run2_b"] = dict(compared=len(common), only_here=sorted(set(rec["values"]) - set(theirs)),
                                        max_abs_diff=float(max(abs(theirs[k] - rec["values"][k]) for k in common))
                                        if common else float("nan"))
        res["claims"][c["id"]] = rec
        log(f"{c['id']}: Delta {delta:+.5g} (SE {te.se:.4g}, df {te.df:.4g}); per-mouse values against the full "
            f"matrix: max |diff| {rec['reduced_vs_full']:.1e}"
            + (f"; against the second run's (b): {rec['vs_run2_b']['compared']} mice, max |diff| "
               f"{rec['vs_run2_b']['max_abs_diff']:.1e}" if "vs_run2_b" in rec else ""))
    save(res)

    # (ii) the controls, chosen by the plan's rule, and the claims' positions among them
    t1 = time.time()
    try:
        if part == "case1":
            names = string_names(string_dir / STRING_FILES[0])
            sel = control_pairs(det, genes, list(PAIR), names, string_dir / STRING_FILES[1])
            sel["string_mapped_genes"] = sum(g in names for g in genes)
            groups_of = {"pairs": sel}
        else:
            groups_of = {"2a": control_genes(det, genes, "Xist"), "2b": control_sets(det, genes, list(Y_GENES))}
        res["controls"] = groups_of
        for key, sel in groups_of.items():
            log(f"(ii) controls {key}: {len(sel['controls'])} drawn from {sel.get('eligible', sel['candidates'])} "
                f"({sel['candidates']} candidates); the claim genes' detection {sel['claim_detection']}; the nearest's "
                f"detection {sel['detection_range']}"
                + (f"; {sel['linked_excluded']} pairs excluded by STRING >= {STRING_MIN_SCORE}; "
                   f"{sel['string_mapped_genes']} of {len(genes)} genes mapped to STRING" if key == "pairs" else ""))
        for c in claims:
            d = designs[c["id"]]
            sel = groups_of["pairs"] if part == "case1" else groups_of[c["id"]]
            deltas, firsts = [], []
            for names_c in sel["controls"]:
                v = per_mouse(c["metric"], X, genes, names_c, tot, d)
                deltas.append(d.delta(v)[0])
                firsts.append(float(np.mean([v[m] for m, gi in zip(d.mice, d.g) if gi == 1])))
            rec = res["claims"][c["id"]]
            pos = position(rec["delta"], deltas, c["direction"])
            first_claim = float(np.mean([rec["values"][m] for m, gi in zip(d.mice, d.g) if gi == 1]))
            rec["ii"] = dict(pos, deltas=deltas, first_group_means=dict(
                claim=first_claim, median=float(np.median(firsts)), min=float(np.min(firsts)), max=float(np.max(firsts))))
            log(*lines_ii(c, rec))
    except Exception as exc:  # reported as it is
        for c in claims:
            res["claims"][c["id"]]["ii"] = dict(error=repr(exc), traceback=traceback.format_exc().splitlines(),
                                                outcome="error")
        log(f"(ii) EXCEPTION {exc!r}", *traceback.format_exc().splitlines()[-6:])
    log(f"(ii): {time.time() - t1:.1f} s")
    save(res)

    # (iii) per mouse, once per measure; (iv) from the per-mouse values
    t1 = time.time()
    share = 2 * (1 - 0.5 * (1 + math.erf(math.log(99) / 50 / math.sqrt(2))))
    log(f"(iii) coupling at strength 50: the shared keep-probability is below 0.01 or above 0.99 in {share:.1%} of "
        "cells (z standard normal)")
    per_cache = {}
    mouse_index = {m: i for i, m in enumerate(sorted(set(obs[case["replicate"]].astype(str))))}
    for c in claims:
        rec = res["claims"][c["id"]]
        d = designs[c["id"]]
        try:
            key = (c["metric"], tuple(c["genes"]))
            if key not in per_cache:
                per_cache[key] = {}
            code = int(hashlib.sha256(f"{c['metric']}:{','.join(c['genes'])}".encode()).hexdigest()[:8], 16)
            for m, r in d.rows.items():
                if m not in per_cache[key]:
                    M = reduced(X, genes, c["genes"], tot, r)
                    per_cache[key][m] = responses(c["metric"], c["genes"], M, (SEED, code, mouse_index[m]))
            per = {m: per_cache[key][m] for m in d.mice}
            for m in d.mice:
                if abs(per[m]["value"] - rec["values"][m]) > 1e-12:
                    raise ValueError(f"(iii)'s value of {m} differs from the per-mouse value")
            rec["iii"] = dict(check_iii(c, d, per), per_mouse=per)
        except Exception as exc:
            rec["iii"] = dict(error=repr(exc), traceback=traceback.format_exc().splitlines(), outcome="error")
        try:
            rec["iv"] = check_iv(c, d, rec["values"], rec["other_values"])
        except Exception as exc:
            rec["iv"] = dict(error=repr(exc), traceback=traceback.format_exc().splitlines(), outcome="error")
        log(*lines_iii(c, rec), *lines_iv(c, rec))
        save(res)
    log(f"(iii) and (iv): {time.time() - t1:.1f} s")

    # (i) the engine, claim by claim (the slow part, last)
    for c in claims:
        t1 = time.time()
        rec = res["claims"][c["id"]]
        d = designs[c["id"]]
        try:
            _, te = d.delta(rec["values"])
            rec["i"] = check_i(c, case, X, genes, obs, d, te, out_dir)
        except Exception as exc:
            rec["i"] = dict(error=repr(exc), traceback=traceback.format_exc().splitlines(), outcome="error")
        rec["i"]["seconds"] = round(time.time() - t1, 1)
        log(*lines_i(c, rec))
        save(res)
    return res


# --------------------------------------------------------------------------- #
# the log
# --------------------------------------------------------------------------- #
def fmt(x, nd=4) -> str:
    return "nan" if x is None or not np.isfinite(x) else f"{x:+.{nd}g}"


def fp(x, nd=4) -> str:
    """A p-value or a size, unsigned."""
    return "nan" if x is None or not np.isfinite(x) else f"{x:.{nd}g}"


def lines_i(c, rec) -> list:
    r = rec.get("i", {})
    if "error" in r:
        return [f"### {c['id']} (i) EXCEPTION {r['error']}", *[f"    {x}" for x in r["traceback"][-6:]]]
    e = r["effect"]
    return [f"### {c['id']} (i) the engine without GATE 4 and GATE 5 ({r['seconds']} s): {r['verdict']} [cause {r['cause']}]",
            f"- raw effect (first - second) {fmt(e.get('raw_effect'))}, p_raw {fp(e.get('p_raw'))}, raw_detected "
            f"{e.get('raw_detected')}; corrected {fmt(e.get('effect'))}, p {fp(e.get('p_perm'))} ({e.get('status')}); "
            f"retained {fmt(e.get('retained'))}; explained_by_depth {e.get('explained_by_depth')}; reversed "
            f"{e.get('reversed_by_correction')}; TOST {e.get('tost_equivalent')} at +-{fp(e.get('effective_sesoi'))} "
            f"(SESOI {fp(r['prereg']['sesoi'])}, {r['prereg'].get('sesoi_scale', 'construct')})",
            "- thinning: " + "; ".join(f"{q['stratum']}: {q['thinned_group']} x {q['depth_ratio']:.4g}"
                                       for q in r["equalization"]),
            "- gates: " + "; ".join(f"{g['gate']} {g['status']}" for g in r["gates"])
            + "; fields: " + "; ".join(f"{k} {v['status']}" for k, v in r["fields"].items() if v),
            f"- (i) {r['outcome']}: {r['why']}"]


def lines_ii(c, rec) -> list:
    r = rec["ii"]
    q = r["quantiles"]
    f = r["first_group_means"]
    return [f"### {c['id']} (ii) specificity: Delta {fmt(rec['delta'])} against {r['n']} controls: "
            f"{r['at_least_as_far']} at least as far in the claimed direction ({c['direction']}), p {r['p_one_sided']:.4g}; "
            f"two-sided p {r['p_two_sided']:.4g}",
            "- the controls' Delta: " + ", ".join(f"{k} {fmt(v)}" for k, v in q.items())
            + f"; their first group's mean {fmt(f['median'])} (median; {fmt(f['min'])} to {fmt(f['max'])}), the claim's "
              f"{fmt(f['claim'])}",
            f"- (ii) {r['outcome']}"]


def lines_iii(c, rec) -> list:
    r = rec.get("iii", {})
    if "error" in r:
        return [f"### {c['id']} (iii) EXCEPTION {r['error']}", *[f"    {x}" for x in r["traceback"][-6:]]]
    out = [f"### {c['id']} (iii) plausibility, by stratum (h the higher group, l the lower):"]
    for st, s in r["strata"].items():
        out.append(f"- {st}: |Delta| {s['abs_delta']:.4g} ({s['higher']} {s['mean'][s['higher']]:.4g} > {s['lower']} "
                   f"{s['mean'][s['lower']]:.4g}); (a) removal in h {fmt(s['removal_higher'])}: "
                   f"{'holds' if s['a_holds'] else 'FAILS'} (|Delta| / its size {s['ratio_a']:.3g}); (b) the dose "
                   f"responses in l " + ", ".join(f"{k}: {fmt(v)}" for k, v in s["dose_lower"].items())
                   + f", the largest at {s['max_dose']}: {'holds' if s['b_holds'] else 'FAILS'} (|Delta| / it "
                     f"{s['ratio_b']:.3g})")
    out.append(f"- (iii) {r['outcome']}" + (f": {', '.join(r['failed'])}" if r["failed"] else ""))
    return out


def lines_iv(c, rec) -> list:
    r = rec.get("iv", {})
    if "error" in r:
        return [f"### {c['id']} (iv) EXCEPTION {r['error']}"]
    out = [f"### {c['id']} (iv) agreement of {r['measures'][0]} and {r['measures'][1]}, by stratum:"]
    for st, s in r["strata"].items():
        x, y = s["own"], s["other"]
        out.append(f"- {st}: {r['measures'][0]} Delta {fmt(x['delta'])} (p {x['p']:.4g}, standardized "
                   f"{fmt(x['standardized'], 3)}); {r['measures'][1]} Delta {fmt(y['delta'])} (p {y['p']:.4g}, "
                   f"standardized {fmt(y['standardized'], 3)})"
                   + ("; opposite signs" if s["opposite"] else "; the same sign")
                   + ("; FLAGGED" if s["flagged"] else ""))
    out.append(f"- (iv) {r['outcome']}")
    return out


# --------------------------------------------------------------------------- #
# the merge: outcomes, "separates", the prediction
# --------------------------------------------------------------------------- #
def outcome_of(rec: dict, check: str) -> str:
    r = rec.get(check)
    return "missing" if r is None else r.get("outcome", "error")


def separates(outcomes: dict) -> str:
    """flags 1a and 1a-males and neither sex claim; partly: one of them and neither sex claim."""
    s4 = [outcomes.get(k) == FLAGS for k in SECTION4]
    sex = [outcomes.get(k) == FLAGS for k in SEX]
    known = all(outcomes.get(k) in (FLAGS, PASSES, NA) for k in SECTION4 + SEX)
    if not known:
        return "undetermined (an exception or a missing claim)"
    if any(sex) or not any(s4):
        return "does not separate"
    return "separates" if all(s4) else "partly"


def prediction(table: dict) -> dict:
    """The owner's: (ii) and (iv) each flag 1a and 1a-males and flag neither 2a nor 2b."""
    parts = {}
    for chk in ("ii", "iv"):
        for k in SECTION4:
            parts[f"({chk}) flags {k}"] = table[chk].get(k) == FLAGS
        for k in SEX:
            parts[f"({chk}) does not flag {k}"] = table[chk].get(k) == PASSES
    return dict(parts=parts, holds=all(parts.values()))


def merge(dirs, out_dir: Path) -> int:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "reports").mkdir(exist_ok=True)
    found = {}
    for d in dirs:
        for f in sorted(Path(d).glob("case*.json")):
            if f.stem in ("case1", "case2"):
                found[f.stem] = (json.loads(f.read_text()), f.with_suffix(".log").read_text()
                                 if f.with_suffix(".log").exists() else "")
                for r in (Path(d) / "reports").glob("*.json"):
                    (out_dir / "reports" / r.name).write_bytes(r.read_bytes())
    claims, text, parts = {}, [], {}  # the parts' logs, then the outcomes
    for part in ("case1", "case2"):
        if part not in found:
            parts[part] = "missing"
            text += [f"## {part}: MISSING (the job did not finish or wrote nothing)", ""]
            continue
        res, log = found[part]
        parts[part] = dict(seconds=res.get("seconds"), complete=res.get("complete"))
        text.append(log)
        if not res.get("complete"):
            text += [f"## {part}: INCOMPLETE (the job stopped; the steps above are those it finished)", ""]
        claims.update(res.get("result", {}).get("claims", {}))
    table = {chk: {c["id"]: outcome_of(claims.get(c["id"], {}), chk) for c in CLAIMS} for chk in CHECKS}
    sep = {chk: separates(table[chk]) for chk in CHECKS}
    pred = prediction(table)
    meta = next(iter(found.values()))[0] if found else {}
    tag = "SYNTHETIC " if meta.get("synthetic") else ""
    start = len(text)
    text += ["## Outcomes by check and claim", "",
             "| check | " + " | ".join(c["id"] for c in CLAIMS) + " | separates §4 from the sex differences |",
             "|---|" + "---|" * (len(CLAIMS) + 1)]
    for chk in CHECKS:
        text.append(f"| ({chk}) | " + " | ".join(table[chk][c["id"]] for c in CLAIMS) + f" | {sep[chk]} |")
    text += ["", "The prediction (the owner, 2026-10-17): «если §4 — артефакт детекции, его отметят проверки (ii) и (iv), "
             "а Xist и Y-гены пройдут»; scored as (ii) and (iv) each flag 1a and 1a-males and flag neither 2a nor 2b:",
             *[f"- {k}: {'yes' if v else 'no'}" for k, v in pred["parts"].items()],
             f"- the prediction {'holds' if pred['holds'] else 'does not hold'} (the part of (iv) on §4 was known before the "
             "run: real_pilot_diagnostics_plan.md, \"Known before the run\")", "",
             "Time per part (s): " + "; ".join(f"{p} {v.get('seconds') if isinstance(v, dict) else v}"
                                               for p, v in parts.items()),
             f"{tag}merged {len(found)} of 2 parts"]
    merged = {k: v for k, v in meta.items() if k not in ("result", "seconds", "part", "complete")}
    merged.update(parts=parts, claims=claims, outcomes=table, separates=sep, prediction=pred,
                  controls={p: found[p][0].get("result", {}).get("controls") for p in found},
                  data={p: found[p][0].get("result", {}).get("data") for p in found})
    (out_dir / "real_pilot_3.json").write_text(json.dumps(merged, indent=1, default=str))
    (out_dir / "real_pilot_3.log").write_text("\n".join(text) + "\n")
    print("\n".join(text[start:]))
    return 0


# --------------------------------------------------------------------------- #
# the dry run: STRING-like files for a fabricated case
# --------------------------------------------------------------------------- #
def fabricate_string(genes, out_dir: Path, seed: int = 7) -> None:
    """protein.info and protein.links in STRING's formats for the genes of a fabricated case: every
    third gene unmapped, a few linked pairs among the rest (some below the threshold)."""
    rng = np.random.default_rng(seed)
    out_dir.mkdir(parents=True, exist_ok=True)
    mapped = [g for i, g in enumerate(genes) if i % 3 != 2 or g in PAIR + S4_CONTROLS]
    ids = {g: f"10090.SYNTHP{i:011d}" for i, g in enumerate(mapped)}
    with gzip.open(out_dir / STRING_FILES[0], "wt", encoding="utf-8") as fh:
        fh.write("#string_protein_id\tpreferred_name\tprotein_size\tannotation\n")
        for g, p in ids.items():
            fh.write(f"{p}\t{g}\t{int(rng.integers(100, 2000))}\tA synthetic protein.\n")
    with gzip.open(out_dir / STRING_FILES[1], "wt", encoding="utf-8") as fh:
        fh.write("protein1 protein2 combined_score\n")
        names = list(ids)
        for _ in range(4 * len(names)):
            a, b = rng.choice(len(names), 2, replace=False)
            s = int(rng.integers(150, 1000))
            fh.write(f"{ids[names[a]]} {ids[names[b]]} {s}\n{ids[names[b]]} {ids[names[a]]} {s}\n")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--data", type=Path, help="the second run's select.json and case files")
    p.add_argument("--string", type=Path, help="the directory of STRING's two files")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--part", choices=("case1", "case2"))
    p.add_argument("--run2", type=Path, help="the second run's real_pilot.json: (b)'s per-mouse values, compared")
    p.add_argument("--merge", nargs="+", type=Path, metavar="PART_DIR")
    p.add_argument("--fabricate-string", type=Path, metavar="DIR", help="STRING-like files for --data's case 1 (dry run)")
    p.add_argument("--string-stats", action="store_true", help="what --string's files hold (the dry run's check)")
    args = p.parse_args(argv)
    if args.merge:
        return merge(args.merge, args.out)
    if args.string_stats:
        st = string_stats(args.string / STRING_FILES[0], args.string / STRING_FILES[1])
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / "string_stats.json").write_text(json.dumps(st, indent=1))
        print(json.dumps(st, indent=1))
        ok = st["probe_mapped"] == list(st["probe"]) and st["probe_linked"] and st["links_400"] > 0
        return 0 if ok else 1
    if args.data is None:
        p.error("--data is required")
    spec_bytes = (args.data / "select.json").read_bytes()
    spec = json.loads(spec_bytes)
    if args.fabricate_string:
        z = np.load(args.data / spec["cases"]["case1"]["data"]["file"], allow_pickle=False)
        fabricate_string([str(g) for g in z["genes"]], args.fabricate_string)
        return 0
    if args.part is None or args.string is None:
        p.error("--part and --string are required")
    args.out.mkdir(parents=True, exist_ok=True)
    lines = []

    def log(*msg):
        for m in msg:
            print(m, flush=True)
            lines.append(str(m))
    tag = "SYNTHETIC " if spec.get("synthetic") else ""
    string_sha = {f: sha256_file(args.string / f) for f in STRING_FILES}
    run2 = json.loads(args.run2.read_text()) if args.run2 else None
    out = dict(synthetic=spec.get("synthetic"), part=args.part, engine=dict(version=MA.__version__,
               path=str(Path(MA.__file__).resolve().parent)), python=platform.python_version(), numpy=np.__version__,
               pandas=pd.__version__, numerics=NUMERICS, select_sha256=hashlib.sha256(spec_bytes).hexdigest(),
               string_sha256=string_sha, run2_sha256=sha256_file(args.run2) if args.run2 else None,
               settings=dict(alpha=ALPHA, seed=SEED, n_controls=N_CONTROLS, nearest_pair=NEAREST_PAIR,
                             nearest_gene=NEAREST_GENE, string_min_score=STRING_MIN_SCORE, strengths=STRENGTHS,
                             folds=FOLDS, k_draws=K_DRAWS, engine_seed=0, estimand="composition"))
    log(f"# {tag}Diagnostics on the second run's data (validation/exploratory/real_pilot_diagnostics_plan.md): {args.part}",
        "", f"engine {MA.__version__} from {Path(MA.__file__).resolve().parent}; python {platform.python_version()}; "
        f"numpy {np.__version__}; pandas {pd.__version__}; numerics {NUMERICS}",
        f"select.json sha256 {out['select_sha256']}; STRING " + "; ".join(f"{k} {v}" for k, v in string_sha.items())
        + (f"; the second run's real_pilot.json sha256 {out['run2_sha256']}" if args.run2 else ""), "")
    t0 = time.time()

    def save(res, complete=False):
        out["result"], out["complete"] = res, complete
        out["seconds"] = round(time.time() - t0, 1)
        (args.out / f"{args.part}.json").write_text(json.dumps(out, indent=1, default=str))
        (args.out / f"{args.part}.log").write_text("\n".join(lines) + "\n")
    try:
        res = run_part(args.part, spec, args.data, args.string, args.out, log, save, run2)
    except Exception as exc:  # reported as it is
        res = dict(status="error", error=repr(exc), traceback=traceback.format_exc().splitlines(), claims={})
        log(f"## {args.part}: EXCEPTION {exc!r}", *res["traceback"])
    log("", f"{tag}{args.part}: {time.time() - t0:.1f} s in all")
    save(res, complete=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
