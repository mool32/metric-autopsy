"""What v0.3's gates add over a donor-level test: four analyzers on every claim card of validation
v1 (exploratory; the plan, committed before the run: v1_analyzers_plan.md).

v1's datasets are development data now (v1.md section 6). Nothing here changes v1's result or the
engine: the engine's functions are imported from the tag's src/, unchanged.

The analyzers, on the same data and claim cards:

(a) cells: the claim's metric on all cells of each group, a two-sided permutation test that
    relabels cells (999 relabelings, the group sizes kept); donors ignored;
(b) donors, without gates: the metric on each donor's cells, a two-sided permutation test over
    donors of the difference of the groups' means (exact up to 20,000 relabelings, else 9,999
    random ones; N5, whose donors have cells in both groups: the donor's difference B - A and all
    sign flips); numpy only, no engine;
(c) the engine with GATE 0, 2, 4 and 5 off: GATE 1, the engine's replicate design, graded rule,
    power and TOST on the uncorrected data, and its verdict rule with the metric's validity taken
    as passed;
(d) the full engine: the published reports.

    python validation/exploratory/v1_analyzers.py run --results RESULTS --compact COMPACT \
        --pilot pilot.json --out cards.jsonl [--workers 4] [--datasets N]
    python validation/exploratory/v1_analyzers.py report --results RESULTS --pilot pilot.json \
        --cards cards.jsonl --table v1_analyzers_cards.csv.gz > v1_analyzers.log

RESULTS is results/panel-v1 with its reports unpacked (manifest.json, reports/); COMPACT the run's
compact backgrounds (compact.tar unpacked); pilot.json the run tag's. `run` rebuilds every dataset
from the key (panel.build), checks it against the manifest's data_sha256 and writes one JSON line
per card; `report` joins them with the published reports and prints the tables.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import sys
from collections import Counter
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation" / "prereg"))
sys.path.insert(0, str(ROOT / "validation" / "exploratory"))
sys.path.insert(0, str(ROOT / "src"))
import frozen  # noqa: E402  (standard library only)

frozen.pin_numerics()  # one BLAS thread per process, set before numpy loads, as the run did

import numpy as np  # noqa: E402

import panel as P  # noqa: E402

ALPHA = 0.05
CELL_PERMS = 999
DONOR_MAX_EXACT = 20_000
DONOR_MC = 9_999
GROUPS = ("A", "B")

# --------------------------------------------------------------------------- #
# the claim's metric, computed without the engine (for (a) and (b))
# --------------------------------------------------------------------------- #


def _seed(card_id: str) -> int:
    """The engine's seed of a card (run_panel._seed): the first 32 bits of its ID's sha256."""
    return int(hashlib.sha256(card_id.encode()).hexdigest()[:8], 16)


class Metric:
    """The claim's metric on any subset of a dataset's cells, from per-cell quantities computed
    once: log1p(CP10k) of the pair's genes and their co-detection (norm_pearson), the mean
    log1p(CP10k) of the score's genes (score), a stream of normal draws (random), 0 (constant)."""

    def __init__(self, X: np.ndarray, genes: list, card: dict, stream: int):
        self.kind = card["metric"]
        self.rng = np.random.default_rng([_seed(card["id"]), stream])
        tot = X.sum(axis=1)
        tot = np.where(tot == 0, 1.0, tot)
        if self.kind == "norm_pearson":
            ia, ib = genes.index(card["gene_pair"][0]), genes.index(card["gene_pair"][1])
            self.both = (X[:, ia] > 0) & (X[:, ib] > 0)
            self.a = np.log1p(X[:, ia] / tot * 1e4)
            self.b = np.log1p(X[:, ib] / tot * 1e4)
        elif self.kind == "score":
            cols = [genes.index(g) for g in card["score_genes"]]
            self.s = np.log1p(X[:, cols] / tot[:, None] * 1e4).mean(axis=1)

    def value(self, mask: np.ndarray) -> float:
        if self.kind == "constant":
            return 0.0
        if self.kind == "random":
            return float(self.rng.normal())
        if self.kind == "score":
            return float(self.s[mask].mean())
        m = mask & self.both
        if m.sum() < 3:
            return 0.0
        a, b = self.a[m], self.b[m]
        if a.std() == 0 or b.std() == 0:
            return 0.0
        return float(np.corrcoef(a, b)[0, 1])

    def values(self, masks: np.ndarray) -> np.ndarray:
        """The metric on every row of a boolean matrix of cell subsets, vectorized."""
        k = masks.shape[0]
        if self.kind == "constant":
            return np.zeros(k)
        if self.kind == "random":
            return self.rng.normal(size=k)
        if self.kind == "score":
            M = masks.astype(np.float64)
            return (M @ self.s) / M.sum(axis=1)
        M = (masks & self.both).astype(np.float64)
        n = M.sum(axis=1)
        sa, sb = M @ self.a, M @ self.b
        saa, sbb, sab = M @ (self.a ** 2), M @ (self.b ** 2), M @ (self.a * self.b)
        va, vb = n * saa - sa ** 2, n * sbb - sb ** 2
        ok = (n >= 3) & (va > 1e-12 * np.maximum(n * saa, 1e-300)) & (vb > 1e-12 * np.maximum(n * sbb, 1e-300))
        r = np.zeros(k)
        r[ok] = (n[ok] * sab[ok] - sa[ok] * sb[ok]) / np.sqrt(va[ok] * vb[ok])
        return r


def call(delta: float, p: float, direction: str) -> str:
    """'+' an effect in the claimed direction, '-' against it, '0' not significant (p >= alpha).
    delta is the change from A to B."""
    if not (np.isfinite(p) and p < ALPHA) or delta == 0:
        return "0"
    return "+" if (delta > 0) == (direction == "increase") else "-"


def two_sided_p(null: np.ndarray, obs: float, plus_one: bool) -> float:
    tol = 1e-12 * max(1.0, abs(obs))
    hits = int(np.sum(np.abs(null) >= abs(obs) - tol))
    return (1 + hits) / (1 + len(null)) if plus_one else hits / len(null)


def cells_test(metric: Metric, group: np.ndarray, rng) -> dict:
    """(a): the metric on all cells of each group; relabel cells at random (sizes kept)."""
    inB = group == "B"
    obs = metric.value(inB) - metric.value(~inB)
    n, nb = len(group), int(inB.sum())
    keys = rng.random((CELL_PERMS, n))
    perm_b = np.zeros((CELL_PERMS, n), dtype=bool)
    np.put_along_axis(perm_b, np.argpartition(keys, nb - 1, axis=1)[:, :nb], True, axis=1)
    null = metric.values(perm_b) - metric.values(~perm_b)
    return dict(delta=float(obs), p=two_sided_p(null, obs, plus_one=True))


def donors_test(metric: Metric, group: np.ndarray, donor: np.ndarray, rng) -> dict:
    """(b): one value per donor; a permutation over donors (paired where donors hold both groups)."""
    donors = list(dict.fromkeys(donor))
    member = {d: set(group[donor == d]) for d in donors}
    if all(len(m) == 2 for m in member.values()):  # N5: every donor has cells in both groups
        d = np.array([metric.value((donor == x) & (group == "B")) - metric.value((donor == x) & (group == "A"))
                      for x in donors])
        signs = np.array(list(itertools.product((1.0, -1.0), repeat=len(d))))
        null = (signs @ d) / len(d)
        return dict(delta=float(d.mean()), p=two_sided_p(null, float(d.mean()), plus_one=False),
                    design="paired", n_donors=len(d), test=f"exact sign flips ({len(signs)})")
    if not all(len(m) == 1 for m in member.values()):
        return dict(delta=float("nan"), p=float("nan"), design="mixed", n_donors=len(donors), test="none")
    v = np.array([metric.value(donor == x) for x in donors])
    isb = np.array([next(iter(member[x])) == "B" for x in donors])
    nb, n = int(isb.sum()), len(v)
    obs = float(v[isb].mean() - v[~isb].mean())
    if comb(n, nb) <= DONOR_MAX_EXACT:
        M = np.zeros((comb(n, nb), n))
        for i, c in enumerate(itertools.combinations(range(n), nb)):
            M[i, list(c)] = 1.0
        exact, test = True, f"exact ({comb(n, nb)} relabelings)"
    else:
        M = np.zeros((DONOR_MC, n))
        for i in range(DONOR_MC):
            M[i, rng.choice(n, size=nb, replace=False)] = 1.0
        exact, test = False, f"Monte Carlo ({DONOR_MC} relabelings)"
    null = (M @ v) / nb - ((1 - M) @ v) / (n - nb)
    return dict(delta=obs, p=two_sided_p(null, obs, plus_one=not exact), design="nested",
                n_donors=n, test=test)


# --------------------------------------------------------------------------- #
# (c): the engine with GATE 0, 2, 4 and 5 off
# --------------------------------------------------------------------------- #


def engine_off(X: np.ndarray, obs, genes: list, card: dict) -> dict:
    """GATE 1, the engine's design, graded rule, power and TOST on the uncorrected data, and its
    verdict rule with the metric's validity taken as passed. Mirrors effect.estimate_effect without
    the correction (no thinning, so no explained-by-depth or reversed-by-correction) and
    report.run_autopsy without GATE 0, 4 and 5."""
    from metric_autopsy import effect as E
    from metric_autopsy import gates as G
    from metric_autopsy import report as Rp
    from metric_autopsy.core import Assessment, GateStatus, SimpleData
    from metric_autopsy.stats import mde

    import run_panel as R
    metric, kw = R.make_run_args(card)
    seed, groups, rep_col = kw["seed"], tuple(card["groups"]), card["replicate_col"]
    prereg = Rp.normalize_prereg(card["prereg"])
    alpha, sesoi = float(prereg["alpha"]), prereg.get("sesoi")
    data = SimpleData(X, obs, genes)
    g1 = G.gate1_qc_parity(data, card["group_col"], groups, within=[], alpha=alpha, seed=seed)
    sd = E._restrict(data, card["group_col"], groups)
    design = E.assess_design(sd, group_col=card["group_col"], groups=groups, within=[], replicate_col=rep_col,
                             estimand=prereg.get("estimand"), min_replicates=int(prereg["min_replicates"]),
                             qc_imbalanced=g1.status == GateStatus.WARN, _sd=sd)
    design.update(correction="none (GATE 2 off)", totals=None)
    rep_design, tier = design["rep_design"], design["tier"]
    detail = dict(correction=design["correction"], estimand=prereg.get("estimand"), sesoi=sesoi, alpha=alpha)
    p_test, est = float("nan"), float("nan")
    if not design["identifiable"]:
        eff = Assessment("NOT_ESTIMABLE", design["reason"], [], detail)
    elif rep_design is None:
        eff = Assessment("INCONCLUSIVE", design["tier_reason"], [], detail)
    else:
        rng = np.random.default_rng(seed)
        vals = E._rep_values(metric, sd, card["group_col"], groups, rep_col, rep_design)
        inf = E._infer(vals, rep_design, groups, tier, alpha, 1000, rng, 20000, int(prereg["min_replicates"]))
        if tier != "insufficient" and inf["tier"] == "insufficient":
            design["tier_reason"] = "; ".join(inf["notes"]) or "insufficient replication"
        tier = inf["tier"] if tier != "insufficient" else tier
        design["tier"] = tier
        lam, lam_info = E.attenuation_lambda(None, 1.0, prereg.get("sesoi_scale", "construct"))
        eff_sesoi = lam * sesoi if sesoi is not None else None
        mde_v = mde(inf["se"], inf["df"], alpha, float(prereg["power"])) if tier != "insufficient" else float("nan")
        design["power"] = dict(mde=mde_v, power=float(prereg["power"]), sesoi=sesoi, lambda_=lam,
                               effective_sesoi=eff_sesoi, attenuation=lam_info)
        design["underpowered"] = (bool(np.isfinite(mde_v) and mde_v > eff_sesoi) if eff_sesoi is not None else None)
        est = float(inf["est"])
        p_test = float(inf["perm"].p if tier == "permutation" else inf["p_t"])
        detail.update(effect=inf["est"], p=p_test, test=inf["test"], tier=tier, perm_method=inf["perm"].method)
        flags = ["PARAMETRIC_ONLY"] if tier == "parametric" else []
        lo2, hi2 = inf["ci_2alpha"]
        equivalent = (eff_sesoi is not None and eff_sesoi > 0 and np.isfinite(lo2) and np.isfinite(hi2)
                      and lo2 > -eff_sesoi and hi2 < eff_sesoi)
        if tier == "insufficient":
            eff = Assessment("INCONCLUSIVE", design["tier_reason"], flags, detail)
        elif inf["detected"]:
            eff = Assessment("DETECTED", "detected (uncorrected)", flags, detail)
        elif equivalent:
            eff = Assessment("NO_DETECTABLE_EFFECT", "TOST within the SESOI (uncorrected)", flags, detail)
        else:
            eff = Assessment("INCONCLUSIVE", "neither detected nor equivalent (uncorrected)", flags, detail)
    da = Rp._design_adequacy(g1, design, prereg, None)
    mv = Assessment("PASS", "taken as passed: GATE 0, 4 and 5 are off (analyzer c)")
    a = Rp.Autopsy(card["metric"], [g1], prereg, mv, da, eff, Rp._replication(None), dict(groups=list(groups)))
    verdict, cause = Rp.decide_cause(a)
    return dict(outcome=P.outcome(verdict, cause), cause=cause, effect=est, p=p_test,
                tier=design.get("tier"), qc=g1.status.value)


# --------------------------------------------------------------------------- #
# run: rebuild every dataset and apply (a), (b), (c)
# --------------------------------------------------------------------------- #
STATE: dict = {}


def analyze(entry: dict) -> list[dict]:
    X, obs, genes, cards = P.build(entry, STATE["bgs"], STATE["pilot"])
    ok = P.dataset_sha256(X, obs, genes) == STATE["want"][entry["id"]]["data_sha256"]
    group, donor = np.asarray(obs["group"]), np.asarray(obs["donor"]).astype(str)
    out, shared = [], None
    for card in cards:
        rec = dict(card=card["id"], dataset=entry["id"], data_ok=bool(ok), metric=card["metric"],
                   direction=card["prereg"]["direction"], replicate_col=card["replicate_col"])
        if ok:
            if shared is None:  # (a) and (b) read the dataset, not the card: one answer per dataset
                first = cards[0]
                a = cells_test(Metric(X, genes, first, 1), group, np.random.default_rng([_seed(first["id"]), 2]))
                b = donors_test(Metric(X, genes, first, 3), group, donor, np.random.default_rng([_seed(first["id"]), 4]))
                shared = dict(a_delta=a["delta"], a_p=a["p"], b_delta=b["delta"], b_p=b["p"],
                              b_design=b["design"], b_donors=b["n_donors"], b_test=b["test"])
            rec.update(shared)
            rec.update(a_call=call(rec["a_delta"], rec["a_p"], rec["direction"]),
                       b_call=call(rec["b_delta"], rec["b_p"], rec["direction"]))
            c = engine_off(X, obs, genes, card)
            rec.update(c_outcome=c["outcome"], c_cause=c["cause"], c_effect=c["effect"], c_p=c["p"],
                       c_tier=c["tier"], c_qc=c["qc"])
        out.append(rec)
    return out


def run(args) -> int:
    import multiprocessing as mp

    import blind
    manifest = json.loads((args.results / "manifest.json").read_text())
    pilot = json.loads(args.pilot.read_text())
    if manifest["pilot_sha256"] != blind.pilot_sha256(pilot):
        raise SystemExit("the pilot differs from the run's")
    bgs, info = blind.load_prepared(args.compact, None)
    if json.dumps(info, sort_keys=True) != json.dumps(manifest["backgrounds"], sort_keys=True):
        raise SystemExit("the compact backgrounds differ from the run's")
    want = {d["id"]: d for d in manifest["datasets"]}
    entries = [e for e in P.assign(P.check_key(manifest["key"]), pilot.get("dropped", ()), pilot.get("pool_size"))
               if e["id"] in want]
    if args.datasets:
        entries = entries[:args.datasets]
    done = set()
    if args.out.exists():
        done = {json.loads(line)["dataset"] for line in args.out.read_text().splitlines() if line.strip()}
    todo = [e for e in entries if e["id"] not in done]
    STATE.update(bgs=bgs, pilot=pilot, want=want)
    print(f"{len(entries)} datasets, {len(done)} done before, {len(todo)} to run on {args.workers} workers",
          file=sys.stderr)
    ctx = mp.get_context("fork")
    with open(args.out, "a") as fh, ctx.Pool(args.workers) as pool:
        for i, recs in enumerate(pool.imap_unordered(analyze, todo, chunksize=4), 1):
            for r in recs:
                fh.write(json.dumps(r) + "\n")
            if i % 250 == 0:
                fh.flush()
                print(f"{i}/{len(todo)}", file=sys.stderr)
    return 0


# --------------------------------------------------------------------------- #
# report: join with the published reports and print the tables
# --------------------------------------------------------------------------- #
ANALYZERS = ("a", "b", "c", "d")
NAMES = {"a": "(a) cells", "b": "(b) donors, no gates", "c": "(c) engine, GATE 0/2/4/5 off", "d": "(d) full engine"}


def effect_call(row: dict, an: str) -> str:
    """'+', '-' or '0' for every analyzer: (c) and (d) SUPPORTED is '+', NOT SUPPORTED (opposite
    direction) '-', anything else '0'."""
    if an in ("a", "b"):
        return row[f"{an}_call"]
    o = row[f"{an}_outcome"]
    return "+" if o == P.SUPPORTED else "-" if o == P.NS_OPPOSITE else "0"


def report(args) -> int:
    import score
    from v1_by_level import fmt
    manifest = json.loads((args.results / "manifest.json").read_text())
    pilot_raw = args.pilot.read_bytes()
    pilot = json.loads(pilot_raw)
    key = P.check_key(manifest["key"])
    reports, _ = score.read_results(args.results, key)
    entries = P.assign(key, pilot.get("dropped", ()), pilot.get("pool_size"))
    published = {r["id"]: r for r in score.card_rows(entries, reports, pilot)}
    by_entry = {e["id"]: e for e in entries}
    cards_raw = args.cards.read_bytes()
    mine = [json.loads(line) for line in cards_raw.decode().splitlines() if line.strip()]
    conds = P.conditions()
    rows = []
    for m in mine:
        pub = published[m["card"]]
        e = by_entry[m["dataset"]]
        c = conds[e["condition"]]
        d_star = sesoi = None
        if c.data == "effect":
            d_star = float(P.delta_of(c.name, e["variant"], e["pair"], pilot)["value"])
        sesoi = P.sesoi_of(c, e["pair"], pilot)
        rows.append(dict(m, condition=c.name, variant=e["variant"], level=pub["level"], pair=int(e["pair"]),
                         background=c.background, truth=pub["truth"], data=c.data, delta_star=d_star, sesoi=sesoi,
                         big=(d_star is not None and abs(d_star) >= sesoi), d_outcome=pub["outcome"],
                         establishable=pub.get("establishable")))
    ok = [r for r in rows if r["data_ok"]]
    out = ["# What the gates add over a donor-level test: four analyzers on validation v1's cards",
           "# (exploratory; plan committed before the run at d3db111: validation/exploratory/v1_analyzers_plan.md)",
           f"# key {key}; pilot.json sha256 {hashlib.sha256(pilot_raw).hexdigest()}; "
           f"per-card results sha256 {hashlib.sha256(cards_raw).hexdigest()}",
           f"# cards {len(rows)} of {len(published)}; datasets rebuilt {len({r['dataset'] for r in rows})}, "
           f"identical to the manifest's {len({r['dataset'] for r in ok})}", ""]
    for an in ANALYZERS:
        out.append(f"{NAMES[an]}")
    out += ["Changes after the plan (none to the analyzers, the sample or the outcomes; the code was tested on the",
            "first 20 and 60 datasets, and the summary's grouping (section 0) was added after a look at the first",
            "3,200 cards): section 4 leaves out N6a, whose random metric cannot agree by construction, and compares",
            "Monte Carlo p values by their detection only; '-' marks an outcome an analyzer cannot give.", ""]
    out += ["'effect': (a), (b) p < 0.05 in the claimed direction; (c), (d) SUPPORTED. 'against': p < 0.05 the other",
            "way; (c), (d) NOT SUPPORTED, opposite direction. Rates are k/n (95% Clopper-Pearson interval).", ""]

    def line(label, rs, what="+"):
        cells = []
        for an in ANALYZERS:
            k = sum(effect_call(r, an) in what for r in rs)
            cells.append(f"{an} {fmt(k, len(rs))}" if rs else f"{an} -")
        return f"  {label:<26} n {len(rs):>5} | " + " | ".join(cells)

    null_conds = [c for c in P.CONDITIONS if c.data == "null" and c.metric == "norm_pearson"]
    nulls = [r for r in ok if r["data"] == "null" and r["metric"] == "norm_pearson"]
    eff = [r for r in ok if r["data"] == "effect"]
    risk = ("valid", "ambiguous")
    pure, art = ("N1", "N5", "N7"), ("N2", "N3", "N8")

    def said(r, an, what):
        """'+'/'-' as effect_call; 'absent': (c), (d) NO DETECTABLE EFFECT; 'silent': (a), (b) not
        significant; 'flag': (d) metric invalid, DEGENERATE METRIC or GATE 0's refusal."""
        if what in ("+", "-", "+-"):
            return effect_call(r, an) in what
        if what == "absent":
            return an in ("c", "d") and r[f"{an}_outcome"] == P.NDE
        if what == "silent":
            return an in ("a", "b") and r[f"{an}_call"] == "0"
        if what == "flag":
            return an == "d" and r["d_outcome"] in (P.NS_INVALID, P.DEGENERATE, P.REFUSAL)
        raise ValueError(what)

    applies = {"silent": ("a", "b"), "absent": ("c", "d"), "flag": ("d",)}

    def srow(label, rs, what):
        cells = [f"{an} {fmt(sum(said(r, an, what) for r in rs), len(rs))}"
                 if rs and an in applies.get(what, ANALYZERS) else f"{an} -" for an in ANALYZERS]
        return f"  {label:<44} n {len(rs):>5} | " + " | ".join(cells)

    out += ["## 0. Summary",
            "False 'effect' on null data (nominal 2.5%):",
            srow("pure nulls N1, N5, N7, metric at risk", [r for r in nulls if r["condition"] in pure and r["truth"] in risk], "+"),
            srow("pure nulls N1, N5, N7, blind metric", [r for r in nulls if r["condition"] in pure and r["truth"] == "blind"], "+"),
            srow("artifacts N2, N3, N8, metric at risk", [r for r in nulls if r["condition"] in art and r["truth"] in risk], "+"),
            srow("artifacts N2, N3, N8, blind metric", [r for r in nulls if r["condition"] in art and r["truth"] == "blind"], "+"),
            srow("N4, 3 against 3 donors (both cards)", [r for r in nulls if r["condition"] == "N4"], "+"),
            "Useless metrics N6a-c (only 'metric invalid' or DEGENERATE METRIC is right; '-': the analyzer cannot say it):",
            srow("'effect'", [r for r in ok if r["condition"] in ("N6a", "N6b", "N6c")], "+"),
            srow("NO DETECTABLE EFFECT claimed", [r for r in ok if r["condition"] in ("N6a", "N6b", "N6c")], "absent"),
            srow("flagged (metric invalid, DEGENERATE)", [r for r in ok if r["condition"] in ("N6a", "N6b", "N6c")], "flag"),
            "Real effects:",
            srow("'effect', valid metric and |Δ*| >= SESOI", [r for r in eff if r["truth"] == "valid" and r["big"]], "+"),
            srow("'effect', ambiguous metric", [r for r in eff if r["truth"] == "ambiguous"], "+"),
            srow("blind metric: 'effect'", [r for r in eff if r["truth"] == "blind"], "+"),
            srow("blind metric: not significant ((a), (b))", [r for r in eff if r["truth"] == "blind"], "silent"),
            srow("blind metric: NO DETECTABLE EFFECT ((c), (d))", [r for r in eff if r["truth"] == "blind"], "absent"),
            srow("blind metric: flagged ((d))", [r for r in eff if r["truth"] == "blind"], "flag"),
            ""]
    out += ["## 1. False 'effect' on null data (N1-N5, N7, N8; nominal 2.5%)", "By condition and variant, all levels:"]
    for c in null_conds:
        for v, _ in c.variants:
            out.append(line(f"{c.name} {v}", [r for r in ok if r["condition"] == c.name and r["variant"] == v]))
    out.append("By expression level, all null conditions:")
    for lv in P.LEVELS:
        out.append(line(lv, [r for r in nulls if r["level"] == lv]))
    out.append("By condition and variant, the high level only:")
    for c in null_conds:
        for v, _ in c.variants:
            out.append(line(f"{c.name} {v} high", [r for r in nulls if r["condition"] == c.name
                                                    and r["variant"] == v and r["level"] == "high"]))
    out.append("By the metric's truth (at risk: valid or ambiguous):")
    out.append(line("at risk", [r for r in nulls if r["truth"] in ("valid", "ambiguous")]))
    out.append(line("blind", [r for r in nulls if r["truth"] == "blind"]))
    out.append("By condition, the cards at risk only:")
    for c in null_conds:
        out.append(line(f"{c.name} at risk", [r for r in nulls if r["condition"] == c.name
                                                and r["truth"] in ("valid", "ambiguous")]))
    out.append("'effect' or 'against' (any false detection; nominal 5%):")
    out.append(line("all nulls", nulls, "+-"))
    out.append(line("at risk", [r for r in nulls if r["truth"] in ("valid", "ambiguous")], "+-"))
    for c in null_conds:
        out.append(line(f"{c.name}", [r for r in nulls if r["condition"] == c.name], "+-"))
    out.append("")

    out += ["## 2. Real effects: the share 'effect' (sensitivity)",
            "The main cell: a valid metric and |Δ*| >= SESOI."]
    out.append(line("valid, |Δ*| >= SESOI", [r for r in eff if r["truth"] == "valid" and r["big"]]))
    out.append(line("valid, |Δ*| < SESOI", [r for r in eff if r["truth"] == "valid" and not r["big"]]))
    out.append(line("ambiguous", [r for r in eff if r["truth"] == "ambiguous"]))
    out.append(line("blind", [r for r in eff if r["truth"] == "blind"]))
    out.append("By condition and variant, a valid metric with |Δ*| >= SESOI:")
    for c in P.CONDITIONS:
        if c.data != "effect":
            continue
        for v, _ in c.variants:
            out.append(line(f"{c.name} {v}", [r for r in eff if r["condition"] == c.name and r["variant"] == v
                                                and r["truth"] == "valid" and r["big"]]))
    out.append("By condition, all truths (and 'against'):")
    for c in P.CONDITIONS:
        if c.data == "effect":
            rs = [r for r in eff if r["condition"] == c.name]
            out.append(line(f"{c.name} effect", rs))
            out.append(line(f"{c.name} against", rs, "-"))
    out.append("By level, all real effects:")
    for lv in P.LEVELS:
        out.append(line(lv, [r for r in eff if r["level"] == lv]))
    out.append("")

    out += ["## 3. Blind and useless metrics: every analyzer's outcomes"]

    def dist(label, rs):
        out.append(f"  {label} (n {len(rs)}):")
        if not rs:
            out.append("    none")
            return
        for an in ANALYZERS:
            if an in ("a", "b"):
                cnt = Counter({"+": "effect", "-": "against", "0": "not significant"}[r[f"{an}_call"]] for r in rs)
            else:
                cnt = Counter(r[f"{an}_outcome"] for r in rs)
            out.append(f"    {an}: " + "; ".join(f"{k} {fmt(v, len(rs))}" for k, v in cnt.most_common()))

    dist("null data, blind metric", [r for r in nulls if r["truth"] == "blind"])
    dist("real effects, blind metric", [r for r in eff if r["truth"] == "blind"])
    for name in ("N6a", "N6b", "N6c"):
        dist(f"{name} ({conds[name].metric})", [r for r in ok if r["condition"] == name])
    dist("null data, at risk", [r for r in nulls if r["truth"] in ("valid", "ambiguous")])
    dist("real effects, valid, |Δ*| >= SESOI", [r for r in eff if r["truth"] == "valid" and r["big"]])
    out.append("")

    out += ["## 4. (c) against the engine's own raw replicate-level test in the published reports",
            "Where (d) estimated the effect, its report holds the raw test at the replicate level (raw_detected,",
            "p_raw; the same data without the correction). (c) runs that test; with an exact permutation they must agree."]
    agree = disagree = mc = mc_same = 0
    examples = []
    for r in ok:
        if r["metric"] == "random":  # a new draw at every call: its values depend on the calls before it
            continue
        rep = json.loads((args.results / "reports" / f"{r['card']}.json").read_text())
        det = (rep.get("fields") or {}).get("effect", {}).get("detail", {})
        if "p_raw" not in det or det.get("p_raw") is None or not np.isfinite(r["c_p"]):
            continue
        if det.get("perm_method") != "exact" and det.get("tier") == "permutation":
            mc += 1
            mc_same += bool(det.get("raw_detected")) == bool(r["c_p"] < ALPHA)
            continue
        same = abs(float(det["p_raw"]) - r["c_p"]) <= 1e-9 * max(1.0, abs(r["c_p"]))
        agree += same
        disagree += not same
        if not same and len(examples) < 5:
            examples.append(f"{r['card']} {r['condition']}: report p_raw {det['p_raw']}, (c) {r['c_p']}")
    out.append(f"  exact permutation or t interval: {agree} agree, {disagree} differ; Monte Carlo permutation (N7): "
               f"{mc}, the same detection at alpha on {mc_same} (their p differ by the Monte Carlo error: other draws)")
    out.append("  (N6a's random metric is not compared: it draws a new number at every call, so its values depend on")
    out.append("  the calls made before the effect, and (d) made GATE 0's, 4's and 5's first)")
    out += [f"  {x}" for x in examples]
    out.append("")
    print("\n".join(out))

    if args.table:
        cols = ["card", "dataset", "condition", "variant", "background", "level", "pair", "truth", "data", "metric",
                "direction", "delta_star", "sesoi", "a_delta", "a_p", "a_call", "b_delta", "b_p", "b_call", "b_test",
                "c_outcome", "c_cause", "c_effect", "c_p", "c_tier", "d_outcome", "data_ok"]
        with gzip.open(args.table, "wt", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(cols)
            for r in sorted(rows, key=lambda r: r["card"]):
                w.writerow(["" if r.get(k) is None else (f"{r[k]:.10g}" if isinstance(r.get(k), float) else r[k])
                            for k in cols])
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--results", type=Path, required=True)
    r.add_argument("--compact", type=Path, required=True)
    r.add_argument("--pilot", type=Path, required=True)
    r.add_argument("--out", type=Path, required=True)
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--datasets", type=int)
    s = sub.add_parser("report")
    s.add_argument("--results", type=Path, required=True)
    s.add_argument("--pilot", type=Path, required=True)
    s.add_argument("--cards", type=Path, required=True)
    s.add_argument("--table", type=Path)
    args = ap.parse_args(argv)
    return run(args) if args.cmd == "run" else report(args)


if __name__ == "__main__":
    raise SystemExit(main())
