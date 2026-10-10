"""The real-data anchors R1-R3 (validation/prereg/v1.md, section 3.2): single datasets whose truth
is public, run with the frozen engine after the panel, reported one by one (not pooled into rates).

* R1 (B2, mice): log-normalized Xist (mean log1p CP10k over a mouse's cells) from female to male,
  a decrease; the Y-gene score (the same on the summed counts of Ddx3y, Eif2s3y, Kdm5d, Uty), an
  increase; within age (development stage), mouse as replicate, composition estimand, SESOI 0.5 on
  the log1p-CP10k scale (delta_min 0.25), the genes raised 2-fold in every cell as GATE 4's module
  signal (a weaker injection, 2-fold in 30% of cells, moves such a level metric by only 0.1-0.2,
  below delta_min). Allowed: SUPPORTED where the engine's replicate rule gives an effect verdict
  (at least 3 mice of each sex: effect._tier counts the replicates per group over the strata;
  3 is the parametric tier, 4 or more the permutation tier), else INCONCLUSIVE. The sham: female
  mice split in two by the public seed (Xist, two-sided): NO DETECTABLE EFFECT or INCONCLUSIVE.
* R2 (B3, mESC sorted by DNA content): R2a the G2M score (mean log1p CP10k of the G2M genes,
  matched to the mouse symbols case-insensitively) from G2M to G1, a decrease, composition; R2b the
  log total endogenous counts, content estimand with the ERCC spike-ins, a decrease; R2c R2b with
  the ERCC rows removed. Each phase is one capture batch (the replicate): R2a and R2b INCONCLUSIVE
  (insufficient replication), R2c UNIDENTIFIABLE.
* R3 (B4): dropped with B4, whose deposit's counts are not raw (section 3.1).

As in the panel, GATE 0's block for a nuisance bias is a refusal: allowed for every anchor, never
definite.

Writes anchors.json (each claim: the verdict, its cause, the four fields, the allowed set, whether
the verdict is in it) and anchors.md.

    python validation/prereg/anchors.py --backgrounds backgrounds.json --data-dir DATA --out anchors.json
"""
from __future__ import annotations

import argparse
import json
from functools import partial
from pathlib import Path

import frozen  # standard library only

frozen.pin_numerics()  # one numerical path, as run_panel.py: set before numpy loads

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import panel as P  # noqa: E402

Y_GENES = ("Ddx3y", "Eif2s3y", "Kdm5d", "Uty")
SESOI = 0.5  # on the log1p-CP10k scale
MODULE_FOLD, MODULE_FRAC = 2.0, 1.0  # GATE 4's signal: the claimed genes 2-fold in every cell


def _load(spec_path: Path, name: str, data_dir: Path, named=()):
    """A background's counts (dense), unique gene symbols and obs. A sparse file is cut to the
    panel's N_GENES genes by the panel's rule, `named` first (DEVIATIONS.md, D1): SimpleData takes
    a dense matrix, and B2's 53,384 genes dense in float64 would take about 8 GB."""
    spec = json.loads(Path(spec_path).read_text())
    if name not in spec:
        return None
    rec = P.resolve_background(spec_path, name, spec[name], data_dir)
    z = np.load(rec["path"], allow_pickle=False)
    from scipy import sparse
    genes = P.unique_names(z["genes"])  # unique symbols, as the panel's backgrounds
    obs = pd.DataFrame({k[4:]: z[k] for k in z.files if k.startswith("obs_")})
    if "X_data" not in z.files:
        return z["X"], genes, obs
    X = sparse.csr_matrix((z["X_data"], z["X_indices"], z["X_indptr"]), shape=tuple(z["X_shape"]))
    # the panel's rule (panel.py, where it plans a background): the named genes, then the most
    # expressed by mean count over the cells (ties by column), N_GENES in all, in column order
    want = {str(g).upper() for g in named}
    first = sorted(j for j, g in enumerate(genes) if g.upper() in want)
    mean, _ = P._col_stats(X)
    order = np.argsort(-mean, kind="mergesort")
    keep = list(dict.fromkeys([*first, *[int(j) for j in order if j not in set(first)]]))
    keep = sorted(keep[:max(P.N_GENES, len(first))])
    return X[:, keep].toarray(), [genes[j] for j in keep], obs


def _cols(genes, wanted):
    up = {g.upper(): j for j, g in enumerate(genes)}
    return [up[w.upper()] for w in wanted if w.upper() in up]


def log_cp10k_mean(data, *, cols):
    """Mean over cells of log1p(CP10k) of the summed counts of `cols`."""
    from metric_autopsy.core import as_dense
    X = as_dense(data.X)
    tot = X.sum(axis=1)
    tot = np.where(tot > 0, tot, 1.0)
    return float(np.mean(np.log1p(X[:, cols].sum(axis=1) / tot * 1e4)))


def log_total(data, *, endogenous):
    from metric_autopsy.core import as_dense
    X = as_dense(data.X)
    return float(np.mean(np.log1p(X[:, endogenous].sum(axis=1))))


def _run(metric, data, groups, within, replicate, prereg, genes_signal, group_col):
    from metric_autopsy import injected_signal, run_autopsy
    a = run_autopsy(metric, data, group_col=group_col, groups=groups, within=within, replicate_col=replicate,
                    signal_test=(injected_signal.module(genes_signal, fold=MODULE_FOLD, frac=MODULE_FRAC)
                                 if genes_signal else None),
                    prereg=dict(prereg, judgment_pending=False, sesoi=SESOI), log_path="off", seed=0)
    return a


def _record(name, a, allowed, note=""):
    label = (a.verdict.split(" — ")[0].split(" [")[0].split(" (")[0]).strip()
    return dict(claim=name, verdict=a.verdict, cause=a.cause, label=label, allowed=list(allowed),
                in_allowed=label in allowed or a.cause == "metric_invalid_gate0", note=note,
                fields={k: (v.status if v else None) for k, v in a.fields().items()})


def r1_allowed(obs) -> tuple:
    """R1's allowed outcome by the engine's replicate rule (effect._tier counts the replicates of
    each group over all strata): SUPPORTED with at least 3 mice of each sex, else INCONCLUSIVE."""
    per_sex = obs.groupby("sex")["mouse"].nunique()
    enough = min(int(per_sex.get("female", 0)), int(per_sex.get("male", 0))) >= 3
    return ("SUPPORTED",) if enough else ("INCONCLUSIVE",)


def r2_allowed(obs) -> tuple:
    """R2a's and R2b's allowed outcomes by the design (v1.md 3.2; the third review found the
    second branch unimplemented): SUPPORTED or INCONCLUSIVE with at least 3 independent captures
    (batches) in each of G1 and G2M, else INCONCLUSIVE (one capture per phase: no replicates)."""
    per_phase = obs.groupby("phase")["batch"].nunique()
    enough = min(int(per_phase.get("G1", 0)), int(per_phase.get("G2M", 0))) >= 3
    return ("SUPPORTED", "INCONCLUSIVE") if enough else ("INCONCLUSIVE",)


def r1(spec, data_dir) -> list[dict]:
    from metric_autopsy import SimpleData
    got = _load(spec, "B2", data_dir, named=("Xist", *Y_GENES))
    if got is None:
        return [dict(claim="R1", skipped="B2 was dropped")]
    X, genes, obs = got
    obs = obs.rename(columns={"donor_id": "mouse", "development_stage": "age"})
    data = SimpleData(X, obs[["sex", "age", "mouse"]].astype(str), genes)
    sexes = obs.assign(n=1).groupby(["age", "sex"])["mouse"].nunique()
    allowed = r1_allowed(obs)
    out = []
    xist = _cols(genes, ["Xist"])
    ygenes = _cols(genes, Y_GENES)
    if xist:
        a = _run(partial(log_cp10k_mean, cols=xist), data, ("female", "male"), ["age"], "mouse",
                 dict(estimand="composition", direction="decrease"), [genes[j] for j in xist], "sex")
        out.append(_record("R1 Xist", a, allowed, f"mice per age and sex: {sexes.to_dict()}"))
    if ygenes:
        a = _run(partial(log_cp10k_mean, cols=ygenes), data, ("female", "male"), ["age"], "mouse",
                 dict(estimand="composition", direction="increase"), [genes[j] for j in ygenes], "sex")
        out.append(_record("R1 Y genes", a, allowed))
    female = obs["sex"].astype(str) == "female"
    if xist and female.any():
        mice = sorted(obs.loc[female, "mouse"].astype(str).unique())
        rng = np.random.default_rng(P.PLAN_SEED)
        half = set(rng.permutation(mice)[: len(mice) // 2])
        sham_obs = obs.loc[female, ["mouse", "age"]].astype(str).copy()
        sham_obs["arm"] = np.where(sham_obs["mouse"].isin(half), "A", "B")
        sham = SimpleData(X[np.where(female)[0]], sham_obs.reset_index(drop=True), genes)
        a = _run(partial(log_cp10k_mean, cols=xist), sham, ("A", "B"), ["age"], "mouse",
                 dict(estimand="composition", direction="two-sided"), [genes[j] for j in xist], "arm")
        out.append(_record("R1 sham (female vs female)", a, ("NO DETECTABLE EFFECT", "INCONCLUSIVE")))
    return out


def r2(spec, data_dir) -> list[dict]:
    from metric_autopsy import SimpleData
    got = _load(spec, "B3", data_dir)
    if got is None:
        return [dict(claim="R2", skipped="B3 was dropped")]
    X, genes, obs = got
    obs = obs[["phase", "batch"]].astype(str)
    keep = np.isin(obs["phase"].to_numpy(), ["G1", "G2M"])
    X, obs = X[keep], obs[keep].reset_index(drop=True)
    data = SimpleData(X, obs, genes)
    g2m = _cols(genes, P.G2M_GENES)
    endo = [j for j, g in enumerate(genes) if not g.upper().startswith("ERCC")]
    out, allowed = [], r2_allowed(obs)
    a = _run(partial(log_cp10k_mean, cols=g2m), data, ("G2M", "G1"), [], "batch",
             dict(estimand="composition", direction="decrease"), [genes[j] for j in g2m], "phase")
    out.append(_record("R2a G2M score", a, allowed, f"{len(g2m)} G2M genes matched"))
    a = _run(partial(log_total, endogenous=endo), data, ("G2M", "G1"), [], "batch",
             dict(estimand="content", direction="decrease", spikein_prefix="ERCC"), None, "phase")
    out.append(_record("R2b total RNA, ERCC present", a, allowed))
    no_ercc = SimpleData(X[:, endo], obs, [genes[j] for j in endo])
    a = _run(partial(log_total, endogenous=list(range(len(endo)))), no_ercc, ("G2M", "G1"), [], "batch",
             dict(estimand="content", direction="decrease", spikein_prefix="ERCC"), None, "phase")
    out.append(_record("R2c total RNA, ERCC removed", a, ("UNIDENTIFIABLE",)))
    return out


def r3(spec, data_dir) -> list[dict]:
    """B4 is never selected (its deposit's counts are not raw; select_backgrounds.fetch_b4)."""
    return [dict(claim="R3", skipped="B4 was dropped (section 3.1)")]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--backgrounds", type=Path, required=True)
    p.add_argument("--data-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, default=Path("anchors.json"))
    args = p.parse_args(argv)
    res = []
    for fn in (r1, r2, r3):
        try:
            res += fn(args.backgrounds, args.data_dir)
        except Exception as exc:  # an anchor that cannot run is reported, not hidden
            res.append(dict(claim=fn.__name__.upper(), error=repr(exc)))
    args.out.write_text(json.dumps(res, indent=1, default=str))
    lines = ["# Anchors R1-R3 (validation/prereg/v1.md, section 3.2)", "", "| Claim | Verdict | Cause | Allowed | In allowed |",
             "|---|---|---|---|---|"]
    for r in res:
        lines.append(f"| {r['claim']} | {r.get('verdict', r.get('skipped', r.get('error')))} | {r.get('cause', '')} | "
                     f"{', '.join(r.get('allowed', []))} | {r.get('in_allowed', '')} |")
    args.out.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
