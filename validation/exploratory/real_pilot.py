#!/usr/bin/env python3
"""The analyses of the pilot on real data (validation/exploratory/real_pilot_plan.md).

Exploratory, for development data only. For every claim of select.json (real_pilot_select.py):

- (a) a test over cells: the claim's metric on all cells of each group, within the claim's strata
  (the mean of the strata's differences); the group labels permuted across the cells 999 times
  within strata; a two-sided p;
- (b) a test over donors, without gates: the metric per donor (per donor and condition where each
  donor is in both groups); the engine's own replicate-level statistics (stats.nested_t and
  permutation_test_nested within strata, or paired_t and the sign flips), exact up to 20,000
  relabelings and else 9,999 random ones; absence by a TOST at the SESOI, the 90% t interval of the
  donor-level difference inside +-SESOI;
- (d) the engine of the tag v0.3.0-prereg: run_autopsy as a user runs the command line (the
  composition estimand, the claim's direction and SESOI on the default construct scale, delta_min
  by default, the judgment gates resolved as in v1's panel and anchor R1, seed 0, run log off), its
  verdict, cause and four fields, raw_detected, and its reading under v0.3.0's known issue
  (SUPPORTED with raw_detected false reads INCONCLUSIVE).

The conclusions are compared: an effect in the claimed direction (either, for a two-sided claim),
one against it, an absence, or none (a refusal). Where (d) differs from (b), the gate is named by
the engine's cause.

    python validation/exploratory/real_pilot.py --data DATA --out OUT

It writes OUT/real_pilot.json, OUT/real_pilot.log and OUT/reports/<claim>.json (the engine's
reports). (a) uses vectorized kernels of the metrics for its permutations; each kernel is checked
against the engine's metric (R1's metric for the means) on the observed groups and donors, and the
check is reported. Every exception is caught per analysis and reported with its traceback.

The second run (the plan's section "Deviation: a second run of cases 1 and 2") adds three things.
(b) reports the design's floor of every claim, v0.4's change 4 computed by hand (design_floor); a
claim whose SESOI is "floor" takes it as its SESOI, in (b)'s TOST and in (d) on the observed scale,
computed before (d) runs. Where the design has several strata, (b) also tests each stratum on its
own (the exact permutation within the stratum). Where select.json carries the flagship audit's
per-mouse table, (b)'s per-mouse values of mi_3bin(Smad3, Col1a1) are compared with it.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
import traceback
from functools import partial
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "prereg"))
import frozen  # noqa: E402  (standard library only)

NUMERICS = frozen.pin_numerics()  # v1's numerical path, set before numpy loads

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import metric_autopsy as MA  # noqa: E402  (the engine: the tag's, installed by the workflow)
from metric_autopsy import SimpleData, injected_signal, metrics, run_autopsy, stats  # noqa: E402
from anchors import log_cp10k_mean  # noqa: E402  (R1's metric, the same function)

ALPHA = 0.05
FLOOR_POWER = 0.95  # the design's floor: absence confirmed by the TOST with this power when there is no effect
CELL_PERMS = 999
DONOR_PERMS, DONOR_MAX_EXACT = 9999, 20000
TOST_LEVEL = 0.90
SEED = 20261015
KERNEL_TOL = 1e-9
CAUSE_GATE = {  # where (d) parts from (b): the gate or field that the engine's cause names
    "degenerate_metric": "metric validity (a degenerate metric)",
    "metric_invalid_gate0": "GATE 0 (the metric's bias under nuisance perturbations)",
    "metric_invalid_gate4": "GATE 4 (no response to the injected signal)",
    "metric_invalid_gate5": "GATE 5 (the controls)",
    "metric_invalid": "metric validity",
    "unidentifiable": "design adequacy (the estimand is unidentifiable)",
    "effect_not_evaluated": "the effect field (not evaluated)",
    "explained_by_depth": "GATE 2 (the correction of depth: explained by depth)",
    "insufficient_replication": "design adequacy (insufficient replication)",
    "absence_untested_metric": "metric validity (the response to signal untested: absence refused)",
    "no_detectable_effect": "the effect field (equivalence at the effective SESOI, after the correction of depth)",
    "effect_inconclusive": "the effect field (after the correction of depth: neither detected nor equivalent)",
    "detected_untested_metric": "metric validity (the response to signal untested)",
    "bias_unsized": "GATE 0 (a bias without a SESOI)",
    "no_direction": "the pre-registration (no direction)",
    "opposite_direction": "the effect field (the opposite direction)",
    "not_replicated": "GATE 6 (replication)",
    "judgment_pending": "the judgment gates",
    "replicated": "", "provisional": "",
}


# --------------------------------------------------------------------------- #
# the claim's metric: the engine's function and a vectorized kernel for (a)
# --------------------------------------------------------------------------- #
def bind(claim: dict, genes: list):
    """The claim's metric as the engine takes it, its unbound pair metric (GATE 5) and gene pair."""
    kind, g = claim["metric"], claim["genes"]
    if kind in ("mi_3bin", "norm_pearson"):
        fn = getattr(metrics, kind)
        return partial(fn, gene_a=g[0], gene_b=g[1]), fn, (g[0], g[1])
    if kind == "log_cp10k_mean":
        return partial(log_cp10k_mean, cols=[genes.index(x) for x in g]), None, None
    raise ValueError(f"unknown metric {kind!r}")


def _bins(x: np.ndarray) -> np.ndarray:
    """metrics._to_3bin: 0 not detected, 1 at or below the median of the detected values, 2 above."""
    out = np.zeros(len(x), dtype=int)
    nz = x > 0
    if nz.any():
        med = np.median(x[nz])
        out[nz] = np.where(x[nz] > med, 2, 1)
    return out


def _mi(a: np.ndarray, b: np.ndarray) -> float:
    """metrics.mi_3bin's sum, in its order, on the bins of a subset."""
    n = len(a)
    if n == 0:
        return 0.0
    mi = 0.0
    for ia in (0, 1, 2):
        pa = np.mean(a == ia)
        if pa == 0:
            continue
        for ib in (0, 1, 2):
            pb = np.mean(b == ib)
            if pb == 0:
                continue
            pab = np.mean((a == ia) & (b == ib))
            if pab > 0:
                mi += pab * np.log(pab / (pa * pb))
    return float(mi)


class Kernel:
    """The claim's metric on any subset of the cells, from per-cell quantities computed once."""

    def __init__(self, X: np.ndarray, genes: list, claim: dict):
        self.kind = claim["metric"]
        cols = [genes.index(g) for g in claim["genes"]]
        tot = X.sum(axis=1)
        tot = np.where(tot == 0, 1.0, tot)
        if self.kind == "mi_3bin":
            self.xa, self.xb = X[:, cols[0]], X[:, cols[1]]
        elif self.kind == "norm_pearson":
            self.both = (X[:, cols[0]] > 0) & (X[:, cols[1]] > 0)
            self.a = np.log1p(X[:, cols[0]] / tot * 1e4)
            self.b = np.log1p(X[:, cols[1]] / tot * 1e4)
        else:
            self.s = np.log1p(X[:, cols].sum(axis=1) / tot * 1e4)

    def value(self, mask: np.ndarray) -> float:
        if self.kind == "mi_3bin":
            return _mi(_bins(self.xa[mask]), _bins(self.xb[mask]))
        if self.kind == "log_cp10k_mean":
            return float(self.s[mask].mean())
        m = mask & self.both
        if m.sum() < 3:
            return 0.0
        a, b = self.a[m], self.b[m]
        if a.std() == 0 or b.std() == 0:
            return 0.0
        return float(np.corrcoef(a, b)[0, 1])

    def values(self, masks: np.ndarray) -> np.ndarray:
        """The metric on every row of a boolean matrix of cell subsets."""
        if self.kind == "mi_3bin":
            return np.array([self.value(m) for m in masks])
        def dot(M, v):  # by einsum's own loops, not BLAS (see real_pilot_select.pooled_pairs)
            return np.einsum("pc,c->p", M, v, optimize=False)
        if self.kind == "log_cp10k_mean":
            M = masks.astype(np.float64)
            return dot(M, self.s) / M.sum(axis=1)
        M = (masks & self.both).astype(np.float64)
        n = M.sum(axis=1)
        sa, sb = dot(M, self.a), dot(M, self.b)
        saa, sbb, sab = dot(M, self.a ** 2), dot(M, self.b ** 2), dot(M, self.a * self.b)
        va, vb = n * saa - sa ** 2, n * sbb - sb ** 2
        ok = (n >= 3) & (va > 1e-12 * np.maximum(n * saa, 1e-300)) & (vb > 1e-12 * np.maximum(n * sbb, 1e-300))
        r = np.zeros(len(n))
        r[ok] = (n[ok] * sab[ok] - sa[ok] * sb[ok]) / np.sqrt(va[ok] * vb[ok])
        return r


# --------------------------------------------------------------------------- #
# (a), (b), (d)
# --------------------------------------------------------------------------- #
def call(delta: float, p: float, direction: str) -> str:
    """'effect' in the claimed direction (either, for a two-sided claim), 'against', or 'none'."""
    if not (np.isfinite(p) and p < ALPHA) or not np.isfinite(delta) or delta == 0:
        return "none"
    if direction == "two-sided":
        return "effect"
    return "effect" if (delta > 0) == (direction == "increase") else "against"


def cells_test(kern: Kernel, group: np.ndarray, strata: np.ndarray, groups, direction, rng) -> dict:
    """(a): Delta = the mean over strata of metric(second group) - metric(first group); the labels
    permuted across the cells within each stratum (sizes kept); p = (1 + #{|null| >= |Delta|}) / 1000."""
    g0, g1 = groups
    parts, null = [], []
    for st in sorted(set(strata)):
        m = strata == st
        a, b = m & (group == g0), m & (group == g1)
        if not a.any() or not b.any():
            continue
        parts.append(kern.value(b) - kern.value(a))
        idx = np.where(m)[0]
        nb = int(b.sum())
        keys = rng.random((CELL_PERMS, len(idx)))
        pick = np.argpartition(keys, nb - 1, axis=1)[:, :nb]
        pb = np.zeros((CELL_PERMS, len(group)), dtype=bool)
        np.put_along_axis(pb, idx[pick], True, axis=1)
        pa = np.zeros_like(pb)
        pa[:, idx] = True
        pa &= ~pb
        null.append(kern.values(pb) - kern.values(pa))
    if not parts:
        return dict(delta=float("nan"), p=float("nan"), strata=0, conclusion="none")
    delta = float(np.mean(parts))
    null = np.mean(np.vstack(null), axis=0)
    tol = 1e-12 * max(1.0, abs(delta))
    p = (1 + int(np.sum(np.abs(null) >= abs(delta) - tol))) / (1 + CELL_PERMS)
    return dict(delta=delta, p=p, strata=len(parts), per_stratum=[float(x) for x in parts],
                conclusion=call(delta, p, direction))


def design_floor(te) -> float:
    """The design's floor (v0.4's change 4, by hand: v0.3.0 has no such rule): the smallest SESOI at
    which the TOST (the 90% t interval inside +-SESOI) confirms absence with power 0.95 when there is
    no effect, from the replicate-level difference's SE and df: SE x (t_0.95,df + t_0.975,df)."""
    if not (np.isfinite(te.se) and te.se > 0 and te.df > 0):
        return float("nan")
    return float(te.se * (stats.t_ppf(0.5 + TOST_LEVEL / 2, te.df) + stats.t_ppf(0.5 + FLOOR_POWER / 2, te.df)))


def stratum_tests(y, g, s, labels, groups, rng) -> dict:
    """Each stratum on its own: the exact permutation of the group labels within it (the engine's
    permutation_test_nested on the stratum's donors), the group means and Delta (second - first)."""
    g0, g1 = groups
    out = {}
    for st, i in labels.items():
        m = s == i
        a, b = y[m & (g == 1)], y[m & (g == 0)]
        if not len(a) or not len(b):
            continue
        pr = stats.permutation_test_nested(y[m], g[m], np.zeros(int(m.sum()), int), n_perm=DONOR_PERMS, rng=rng,
                                           max_exact=DONOR_MAX_EXACT)
        out[st] = dict(n={str(g0): len(a), str(g1): len(b)}, mean={str(g0): float(a.mean()), str(g1): float(b.mean())},
                       delta=float(b.mean() - a.mean()), p=float(pr.p), method=pr.method, n_relabelings=int(pr.n),
                       min_attainable_p=float(pr.min_attainable_p))
    return out


def donors_test(metric, data, group, donor, strata, groups, direction, sesoi, rng) -> dict:
    """(b): one value per donor (or per donor and group where every donor is in both groups); the
    engine's replicate-level statistics on the uncorrected data, without gates; TOST at the SESOI
    (at the design's floor where the SESOI is "floor"); each stratum on its own where there are
    several."""
    g0, g1 = groups
    reps = sorted(set(donor))
    member = {r: set(group[donor == r]) & {g0, g1} for r in reps}
    if all(len(m) == 1 for m in member.values()):
        design = "nested"
        vals = {r: float(metric(data[donor == r])) for r in reps}
        st_of = {r: next(iter(set(strata[donor == r]))) for r in reps}
        labels = {s: i for i, s in enumerate(sorted(set(st_of.values())))}
        y = np.array([vals[r] for r in reps])
        g = np.array([1 if g0 in member[r] else 0 for r in reps])
        s = np.array([labels[st_of[r]] for r in reps])
        te = stats.nested_t(y, g, s)
        perm = stats.permutation_test_nested(y, g, s, n_perm=DONOR_PERMS, rng=rng, max_exact=DONOR_MAX_EXACT)
        counts = {str(g0): int(g.sum()), str(g1): int((1 - g).sum())}
        per_stratum = {st: {str(g0): int(((s == i) & (g == 1)).sum()), str(g1): int(((s == i) & (g == 0)).sum())}
                       for st, i in labels.items()}
        by_stratum = stratum_tests(y, g, s, labels, groups, rng) if len(labels) > 1 else {}
    elif all(len(m) == 2 for m in member.values()):
        design = "paired"
        vals = {r: float(metric(data[(donor == r) & (group == g0)])) - float(metric(data[(donor == r) & (group == g1)]))
                for r in reps}
        d = np.array([vals[r] for r in reps])
        te = stats.paired_t(d)
        perm = stats.permutation_test_paired(d, n_perm=DONOR_PERMS, rng=rng, max_exact=DONOR_MAX_EXACT)
        counts, per_stratum, by_stratum = dict(pairs=len(d)), {}, {}
    else:
        return dict(design="mixed", conclusion="none", note="donors partly in both groups")
    delta = -float(te.est)  # the engine's estimate is first - second; Delta is second - first
    lo, hi = te.ci(TOST_LEVEL)
    ci = (-hi, -lo)
    floor = design_floor(te)
    use = floor if sesoi == "floor" else float(sesoi)
    tost = bool(np.isfinite(lo) and np.isfinite(hi) and np.isfinite(use) and -use < lo and hi < use)
    concl = call(delta, perm.p, direction)
    if concl == "none" and tost:
        concl = "absence"
    return dict(design=design, delta=delta, se=float(te.se), df=float(te.df), t_method=te.method, p=float(perm.p),
                perm_method=perm.method, perm_n=int(perm.n), min_attainable_p=float(perm.min_attainable_p),
                ci90=[float(ci[0]), float(ci[1])], floor=floor, sesoi_used=use, tost=tost, counts=counts,
                per_stratum=per_stratum, by_stratum=by_stratum, values={r: vals[r] for r in reps}, conclusion=concl)


def label_of(verdict: str) -> str:
    return verdict.split(" — ")[0].split(" [")[0].split(" (")[0].strip()


def d_conclusion(label: str, cause: str) -> str:
    if label == "SUPPORTED":
        return "effect"
    if cause == "opposite_direction":
        return "against"
    if label == "NO DETECTABLE EFFECT":
        return "absence"
    return "none"


def engine_run(claim: dict, case: dict, X: np.ndarray, obs: pd.DataFrame, genes: list, out_dir: Path) -> dict:
    """(d): run_autopsy with the case's settings; its verdict, cause, fields and the known-issue reading."""
    metric, pair_metric, gene_pair = bind(claim, genes)
    sig = claim["signal"]
    if sig["kind"] == "coupling":
        signal = injected_signal.coupling(*claim["genes"], strength=sig["strength"])
    else:
        signal = injected_signal.module(list(claim["genes"]), fold=sig["fold"], frac=sig["frac"])
    ctl = claim.get("controls")
    cols = list(dict.fromkeys([case["group_col"], case["replicate"], *case["within"]]))
    data = SimpleData(X.copy(), obs[cols].astype(str).reset_index(drop=True), list(genes))
    prereg = dict(estimand="composition", direction=claim["direction"], sesoi=float(claim["sesoi"]),
                  judgment_pending=False)
    if claim.get("sesoi_scale"):
        prereg["sesoi_scale"] = claim["sesoi_scale"]
    a = run_autopsy(metric, data, group_col=case["group_col"], groups=tuple(case["groups"]), gene_pair=gene_pair,
                    within=list(case["within"]), pair_metric=pair_metric if ctl else None,
                    pos_pair=tuple(ctl["pos"]) if ctl else None, neg_pair=tuple(ctl["neg"]) if ctl else None,
                    prereg=prereg, replicate_col=case["replicate"], signal_test=signal, seed=0, log_path="off")
    a.metric_name = claim["metric"]
    rep = a.to_dict()
    try:
        text = json.dumps(rep, indent=1, allow_nan=False)
    except ValueError:  # a NaN in the report: kept as it is, not lost
        text = json.dumps(rep, indent=1, default=str)
    (out_dir / "reports" / f"{claim['id']}.json").write_text(text)
    ef = (rep["fields"].get("effect") or {}).get("detail") or {}
    label = label_of(a.verdict)
    raw = ef.get("raw_detected")
    known = "INCONCLUSIVE" if (label == "SUPPORTED" and raw is False) else label
    keep = ("effect", "raw_effect", "ci", "p_perm", "p_raw", "p_t", "raw_detected", "tost_equivalent", "sesoi",
            "effective_sesoi", "explained_by_depth", "retained", "reversed_by_correction", "tier", "test",
            "correction", "n_replicates", "design_kind")
    return dict(verdict=a.verdict, label=label, cause=a.cause, conclusion=d_conclusion(label, a.cause),
                known_issue_label=known, known_issue_conclusion=d_conclusion(known, a.cause if known == label else ""),
                raw_detected=raw, fields={k: (v and dict(status=v["status"], reason=v["reason"], flags=v["flags"]))
                                          for k, v in rep["fields"].items()},
                effect={k: ef.get(k) for k in keep if k in ef},
                gates=[dict(gate=g["gate"], name=g["name"], status=g["status"], message=g["message"]) for g in rep["gates"]],
                prereg=prereg, report=f"reports/{claim['id']}.json",
                report_sha256=hashlib.sha256(text.encode()).hexdigest())


def gate_of(d: dict) -> str:
    """The gate or field behind (d)'s conclusion where it differs from (b)'s."""
    cause = d.get("cause")
    if "error" in d:
        return "an exception in the engine"
    if cause in ("replicated", "provisional"):
        return ("GATE 2 (the correction of depth: the effect is detected at equal depth and not before it, "
                "raw_detected false)" if d.get("raw_detected") is False
                else "the effect field (the effect detected at equal depth)")
    return CAUSE_GATE.get(cause) or str(cause)


# --------------------------------------------------------------------------- #
def kernel_check(kern: Kernel, metric, data, group, donor, strata, groups) -> float:
    """The largest difference between the kernel and the engine's metric on the observed groups (per
    stratum) and on each donor's cells (per group where donors are in both)."""
    diffs = []
    for st in sorted(set(strata)):
        for g in groups:
            m = (strata == st) & (group == g)
            if m.any():
                diffs.append(abs(kern.value(m) - float(metric(data[m]))))
    for r in sorted(set(donor)):
        for g in groups:
            m = (donor == r) & (group == g)
            if m.any():
                diffs.append(abs(kern.value(m) - float(metric(data[m]))))
    return float(max(diffs)) if diffs else float("nan")


def run_case(name: str, case: dict, data_dir: Path, out_dir: Path, log, on_claim=None) -> dict:
    path = data_dir / case["data"]["file"]
    got = hashlib.sha256(path.read_bytes()).hexdigest()
    if got != case["data"]["sha256"]:
        raise ValueError(f"{path.name}: sha256 {got}, select.json has {case['data']['sha256']}")
    z = np.load(path, allow_pickle=False)
    X = z["X"].astype(np.float64)
    genes = [str(g) for g in z["genes"]]
    obs = pd.DataFrame({k[4:]: z[k] for k in z.files if k.startswith("obs_")})
    gc, groups, rep = case["group_col"], tuple(case["groups"]), case["replicate"]
    group = obs[gc].astype(str).to_numpy()
    donor = obs[rep].astype(str).to_numpy()
    strata = (obs[case["within"]].astype(str).agg("|".join, axis=1).to_numpy() if case["within"]
              else np.array([""] * len(obs)))
    data = SimpleData(X, obs.astype(str), genes)
    design = obs.groupby([*case["within"], gc])[rep].nunique().rename("donors").reset_index()
    log(f"## {name}: {len(obs)} cells x {len(genes)} genes; groups {groups[0]!r} -> {groups[1]!r}"
        + (f" within {case['within']}" if case["within"] else "") + f"; the unit: {rep}")
    log("Donors per group" + (" and stratum" if case["within"] else "") + ":",
        *[f"  {', '.join(f'{c}={row[c]}' for c in design.columns)}" for _, row in design.iterrows()])
    records = []
    for claim in case["claims"]:
        rec = dict(id=claim["id"], metric=claim["metric"], genes=claim["genes"], direction=claim["direction"],
                   sesoi=claim["sesoi"], signal=claim["signal"], controls=claim.get("controls"))
        for k in ("sesoi_scale", "sesoi_note", "seed_id"):
            if claim.get(k):
                rec[k] = claim[k]
        if claim.get("skip"):
            rec["skipped"] = claim["skip"]
            records.append(rec)
            log(f"{claim['id']}: skipped ({claim['skip']})")
            continue
        metric, _, _ = bind(claim, genes)
        kern = Kernel(X, genes, claim)
        t0 = time.time()
        try:
            rec["kernel_max_diff"] = kernel_check(kern, metric, data, group, donor, strata, groups)
            rec["kernel_ok"] = bool(rec["kernel_max_diff"] <= KERNEL_TOL)
            if "pooled_r" in claim:
                eng = float(metric(data))
                rec["pooled_r_engine"], rec["pooled_r_diff"] = eng, abs(eng - claim["pooled_r"])
        except Exception as exc:
            rec["kernel_error"] = dict(error=repr(exc), traceback=traceback.format_exc().splitlines())
        rec["seconds_check"] = round(time.time() - t0, 2)
        cs = claim_seed(claim.get("seed_id") or claim["id"])  # a SESOI variant shares its claim's permutations

        def engine_claim(claim=claim, rec=rec):  # the floor, from (b), before the engine runs the claim
            if claim["sesoi"] != "floor":
                return engine_run(claim, case, X, obs, genes, out_dir)
            floor = rec["b"].get("floor", float("nan"))
            if not np.isfinite(floor):
                raise ValueError("no design floor: (b) has no finite SE")
            return engine_run(dict(claim, sesoi=floor), case, X, obs, genes, out_dir)
        for an, fn in (("a", lambda: cells_test(kern, group, strata, groups, claim["direction"],
                                                np.random.default_rng([SEED, cs, 1]))),
                       ("b", lambda: donors_test(metric, data, group, donor, strata, groups, claim["direction"],
                                                 claim["sesoi"], np.random.default_rng([SEED, cs, 2]))),
                       ("d", engine_claim)):
            t0 = time.time()
            try:
                rec[an] = fn()
            except Exception as exc:  # reported as it is
                rec[an] = dict(error=repr(exc), traceback=traceback.format_exc().splitlines(), conclusion="error")
            rec[an]["seconds"] = round(time.time() - t0, 2)
        rec["sesoi_used"] = rec["b"].get("sesoi_used", claim["sesoi"])
        if case.get("audit_mice") and claim["metric"] == "mi_3bin" and "values" in rec["b"]:
            rec["b_vs_audit"] = against_audit(rec["b"]["values"], case["audit_mice"], claim["genes"])
        b, d, dk = rec["b"].get("conclusion"), rec["d"].get("conclusion"), rec["d"].get("known_issue_conclusion")
        rec["d_vs_b"] = "same" if d == b else f"changed by {gate_of(rec['d'])}"
        rec["d_known_vs_b"] = ("same" if dk == b else
                               "changed by v0.3.0's known issue (SUPPORTED with raw_detected false)" if d == b else
                               rec["d_vs_b"])
        records.append(rec)
        log(*claim_lines(rec))
        if on_claim is not None:  # the outputs so far, written after every claim (a job may be stopped)
            on_claim(dict(cells=int(len(obs)), genes=int(len(genes)), claims=records, data_sha256=got, status="running"))
    return dict(cells=int(len(obs)), genes=int(len(genes)), design=json.loads(design.to_json(orient="records")),
                claims=records, data_sha256=got)


def against_audit(values: dict, audit_mice: list, genes) -> dict:
    """(b)'s per-mouse values of mi_3bin against the flagship audit's (section E: Smad3, Col1a1, the
    same cells by its own selection)."""
    if list(genes) != ["Smad3", "Col1a1"]:
        return dict(note="the audit's table is of mi_3bin(Smad3, Col1a1) only")
    theirs = {str(m["donor_id"]): float(m["mi_3bin"]) for m in audit_mice}
    both = sorted(set(values) & set(theirs))
    diffs = [abs(values[k] - theirs[k]) for k in both]
    return dict(compared=len(both), only_here=sorted(set(values) - set(theirs)), only_audit=sorted(set(theirs) - set(values)),
                max_abs_diff=float(max(diffs)) if diffs else float("nan"))


def fmt(x, nd=4) -> str:
    return "nan" if x is None or not np.isfinite(x) else f"{x:+.{nd}f}"


def claim_lines(rec: dict) -> list:
    a, b, d = rec.get("a", {}), rec.get("b", {}), rec.get("d", {})
    out = [f"### {rec['id']}: {rec['metric']}({', '.join(rec['genes'])}), {rec['direction']}, SESOI {rec['sesoi']}"]
    if "kernel_max_diff" in rec:
        out.append(f"- kernel against the engine's metric: max |diff| {rec['kernel_max_diff']:.2e} "
                   f"({'ok' if rec['kernel_ok'] else 'MISMATCH'})"
                   + (f"; pooled r {rec['pooled_r_engine']:+.6f} (|diff| from the selection {rec['pooled_r_diff']:.1e})"
                      if "pooled_r_engine" in rec else ""))
    if "kernel_error" in rec:
        out.append(f"- kernel check: EXCEPTION {rec['kernel_error']['error']}")
    out.append(f"- (a) cells: Delta {fmt(a.get('delta'))}, p {a.get('p', float('nan')):.4g} over {a.get('strata', 0)} "
               f"strata -> {a.get('conclusion')} ({a.get('seconds')} s)" if "error" not in a else
               f"- (a) EXCEPTION {a['error']}")
    if "error" in b:
        out.append(f"- (b) EXCEPTION {b['error']}")
    elif b.get("design") == "mixed":
        out.append("- (b) donors: mixed design, no test")
    else:
        use = b.get("sesoi_used", rec["sesoi"])
        out.append(f"- (b) donors ({b['design']}, {b['counts']}): Delta {fmt(b['delta'])}, p {b['p']:.4g} "
                   f"({b['perm_method']}, {b['perm_n']}; min attainable {b['min_attainable_p']:.3g}), 90% CI "
                   f"[{fmt(b['ci90'][0])}, {fmt(b['ci90'][1])}], TOST at +-{use:.4g}"
                   + (" (the design's floor)" if rec["sesoi"] == "floor" else "") + f": {b['tost']} -> "
                   f"{b['conclusion']} ({b['seconds']} s)")
        if "floor" in b:
            out.append(f"    the design's floor (v0.4's change 4, by hand): {b['floor']:.4g} = SE {b['se']:.4g} x "
                       f"(t_0.95 + t_0.975 at df {b['df']:.4g}); the SESOI used {use:.4g}")
        for st, r in (b.get("by_stratum") or {}).items():
            g0, g1 = list(r["n"])
            out.append(f"    stratum {st}: {r['n'][g0]} against {r['n'][g1]} donors, means {r['mean'][g0]:.4g} against "
                       f"{r['mean'][g1]:.4g}, Delta {fmt(r['delta'])}, p {r['p']:.4g} ({r['method']}, {r['n_relabelings']}; "
                       f"min attainable {r['min_attainable_p']:.3g})")
        if rec.get("b_vs_audit"):
            v = rec["b_vs_audit"]
            out.append("    against the audit's per-mouse mi_3bin (section E): " + (
                v["note"] if "note" in v else f"{v['compared']} mice compared, max |diff| {v['max_abs_diff']:.2e}; "
                f"only here {v['only_here'] or 'none'}, only in the audit {v['only_audit'] or 'none'}"))
    if "error" in d:
        out.append(f"- (d) EXCEPTION {d['error']}")
        out += [f"    {ln}" for ln in d["traceback"][-6:]]
    else:
        e = d.get("effect", {})
        out.append(f"- (d) engine: {d['verdict']} [cause {d['cause']}] -> {d['conclusion']} ({d['seconds']} s)")
        out.append("    fields: " + "; ".join(f"{k} {v['status']}" for k, v in d["fields"].items() if v)
                   + f"; raw_detected {d['raw_detected']}; corrected effect (first - second) {fmt(e.get('effect'))}, "
                   f"raw {fmt(e.get('raw_effect'))}; effective SESOI {e.get('effective_sesoi')}")
        out.append("    gates: " + "; ".join(f"{g['gate']} {g['status']}" for g in d["gates"]))
        out.append(f"    under the known issue: {d['known_issue_label']} -> {d['known_issue_conclusion']}")
    out.append(f"- (d) against (b): {rec.get('d_vs_b')}; under the known issue: {rec.get('d_known_vs_b')}")
    return out


def summary_lines(results: dict) -> list:
    out = ["## Summary", "",
           "| case | claim | metric | (a) | (b) | (d) verdict | cause | (d) | (d) known issue | (d) against (b) | s (a/b/d) |",
           "|---|---|---|---|---|---|---|---|---|---|---|"]
    for name, res in results.items():
        for r in res.get("claims", []):
            if "skipped" in r:
                out.append(f"| {name} | {r['id']} | {r['metric']} | skipped: {r['skipped']} |||||||||")
                continue
            a, b, d = r.get("a", {}), r.get("b", {}), r.get("d", {})
            out.append(f"| {name} | {r['id']} | {r['metric']}({', '.join(r['genes'])}) | {a.get('conclusion')} | "
                       f"{b.get('conclusion')} | {d.get('label', 'error')} | {d.get('cause', '')} | {d.get('conclusion')} | "
                       f"{d.get('known_issue_conclusion')} | {r.get('d_vs_b')} | {a.get('seconds')}/{b.get('seconds')}/"
                       f"{d.get('seconds')} |")
    return out


PARTS = {"case1": ("case1", ""), "case2": ("case2", ""), "case3a": ("case3", "3a"), "case3b": ("case3", "3b"),
         "case3c": ("case3", "3c"),  # the workflow's parallel jobs: a case, or a case's claims by prefix
         "case1a": ("case1", ("1a",)), "case1f": ("case1", ("1a-floor",)),  # the second run: case 1's claims
         "case1b": ("case1", ("1b",))}                                       # one by one (a tuple: exact ids)
FIRST_RUN_PARTS = ("case1", "case2", "case3a", "case3b", "case3c")  # what a merge expects by default


def claim_seed(claim_id: str) -> int:
    """The seed of a claim's permutations: the first 32 bits of its id's sha256 (the same in any part)."""
    return int(hashlib.sha256(claim_id.encode()).hexdigest()[:8], 16)


def merge(dirs, out_dir: Path, expect=None) -> int:
    """The parts' outputs (each a directory with <part>.json, <part>.log and reports/) into one
    real_pilot.json and real_pilot.log; a part of `expect` (by default the first run's five) that is
    missing is recorded as missing."""
    expect = list(expect or FIRST_RUN_PARTS)
    (out_dir / "reports").mkdir(parents=True, exist_ok=True)
    found = {}
    for d in dirs:
        for f in sorted(Path(d).glob("case*.json")):
            if f.stem in PARTS:
                found[f.stem] = (json.loads(f.read_text()), (f.with_suffix(".log")).read_text()
                                 if f.with_suffix(".log").exists() else "")
                for r in (Path(d) / "reports").glob("*.json"):
                    (out_dir / "reports" / r.name).write_bytes(r.read_bytes())
    meta = next(iter(found.values()))[0] if found else {}
    merged = {k: v for k, v in meta.items() if k not in ("cases", "seconds", "part")}
    merged.update(parts={}, cases={})
    text = []
    for part in expect:
        case_name = PARTS[part][0]
        if part not in found:
            merged["parts"][part] = dict(status="missing (the job did not finish or wrote nothing)")
            text += [f"## {part}: MISSING (the job did not finish or wrote nothing)", ""]
            continue
        res, log = found[part]
        merged["parts"][part] = dict(seconds=res.get("seconds"), complete=res.get("complete"))
        text.append(log)
        if not res.get("complete"):
            text += [f"## {part}: INCOMPLETE (the job stopped; the claims above are those it finished)", ""]
        case = res["cases"].get(case_name, {})
        into = merged["cases"].setdefault(case_name, {k: v for k, v in case.items() if k != "claims"})
        into.setdefault("claims", []).extend(case.get("claims", []))
        if case.get("status") == "error":
            into["status"], into["error"] = "error", case.get("error")
    tag = "SYNTHETIC " if merged.get("synthetic") else ""
    text += summary_lines(merged["cases"])
    text += ["", "Time per part (s): " + "; ".join(f"{p} {v.get('seconds', v.get('status'))}"
                                                  for p, v in merged["parts"].items())]
    text.append(f"{tag}merged {len(set(found) & set(expect))} of {len(expect)} parts")
    (out_dir / "real_pilot.json").write_text(json.dumps(merged, indent=1, default=str))
    (out_dir / "real_pilot.log").write_text("\n".join(text) + "\n")
    print("\n".join(text[-(len(summary_lines(merged['cases'])) + 3):]))
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--data", type=Path, help="real_pilot_select.py's output directory")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--part", choices=sorted(PARTS), help="only this part (the workflow runs the parts in parallel)")
    p.add_argument("--merge", nargs="+", type=Path, metavar="PART_DIR", help="merge the parts' outputs into OUT")
    p.add_argument("--expect", nargs="+", choices=sorted(PARTS),
                   help="the parts the merge expects (default: the first run's five)")
    args = p.parse_args(argv)
    if args.merge:
        return merge(args.merge, args.out, args.expect)
    if args.data is None:
        p.error("--data is required")
    (args.out / "reports").mkdir(parents=True, exist_ok=True)
    spec_bytes = (args.data / "select.json").read_bytes()
    spec = json.loads(spec_bytes)
    lines = []

    def log(*msg):
        for m in msg:
            print(m, flush=True)
            lines.append(str(m))
    tag = "SYNTHETIC " if spec.get("synthetic") else ""
    t0 = time.time()
    second = ", second run (the plan's section \"Deviation: a second run of cases 1 and 2\")" if spec.get("run") == "second" else ""
    log(f"# {tag}The pilot's analyses (validation/exploratory/real_pilot_plan.md){second}" + (f": {args.part}" if args.part else ""),
        "", f"engine {MA.__version__} from {Path(MA.__file__).resolve().parent}; python {platform.python_version()}; "
        f"numpy {np.__version__}; pandas {pd.__version__}; numerics {NUMERICS}",
        f"select.json sha256 {hashlib.sha256(spec_bytes).hexdigest()}; Census release {spec.get('release')}", "")
    out = dict(synthetic=spec.get("synthetic"), run=spec.get("run", "first"), part=args.part,
               engine=dict(version=MA.__version__, path=str(Path(MA.__file__).resolve().parent)),
               python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, numerics=NUMERICS,
               select_sha256=hashlib.sha256(spec_bytes).hexdigest(), settings=dict(
                   alpha=ALPHA, cell_perms=CELL_PERMS, donor_perms=DONOR_PERMS, donor_max_exact=DONOR_MAX_EXACT,
                   tost_level=TOST_LEVEL, floor_power=FLOOR_POWER, seed=SEED, engine_seed=0, judgment_pending=False,
                   estimand="composition"),
               cases={})
    want_case, prefix = PARTS[args.part] if args.part else (None, "")
    picks = (lambda cid: cid in prefix) if isinstance(prefix, tuple) else (lambda cid: cid.startswith(prefix))
    stem = args.part or "real_pilot"

    def save(complete: bool):
        out["complete"] = complete
        (args.out / f"{stem}.json").write_text(json.dumps(out, indent=1, default=str))
        (args.out / f"{stem}.log").write_text("\n".join(lines) + "\n")
    for name, case in spec["cases"].items():
        if want_case and name != want_case:
            continue
        t1 = time.time()
        if case.get("status") != "run":
            out["cases"][name] = dict(status=case.get("status"), reason=case.get("reason") or case.get("error"))
            log(f"## {name}: not run ({out['cases'][name]['reason']})", "")
            continue
        case = dict(case, claims=[c for c in case["claims"] if picks(c["id"])])

        def partial_save(res, name=name):
            out["cases"][name] = res
            save(False)
        try:
            res = run_case(name, case, args.data, args.out, log, on_claim=partial_save)
            res["status"] = "run"
        except Exception as exc:  # reported as it is
            res = dict(status="error", error=repr(exc), traceback=traceback.format_exc().splitlines())
            log(f"## {name}: EXCEPTION {exc!r}", *res["traceback"])
        res["seconds"] = round(time.time() - t1, 1)
        out["cases"][name] = res
        log(f"{name}{' ' + ', '.join(prefix) if isinstance(prefix, tuple) else ' ' + prefix if prefix else ''}: "
            f"{res['seconds']} s", "")
    out["seconds"] = round(time.time() - t0, 1)
    if not args.part:
        log(*summary_lines(out["cases"]))
    log("", f"{tag}analyses{' (' + args.part + ')' if args.part else ''}: {out['seconds']} s in all")
    save(True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
