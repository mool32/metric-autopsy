#!/usr/bin/env python3
"""Flagship-data audit for preprint §4 (Tabula Muris Senis FACS). Validation-plan step 1.

The audit checks §4 *as computed*, so it is pinned to the v0.1.1 engine: run it from a
worktree of the v0.1.1 tag, not from a checkout with reworked verdict logic.

    git worktree add ../metric-autopsy-v0.1.1 v0.1.1
    pip install -e ../metric-autopsy-v0.1.1 cellxgene-census
    python validation/flagship_audit/audit_tms.py --out validation/flagship_audit/out

Network: census.cellxgene.cziscience.com and the Census S3 bucket
(cellxgene-census-public-us-west-2.s3.us-west-2.amazonaws.com). ``--source-h5ad`` also
downloads each TMS FACS source .h5ad, several GB, to look for ERCC spike-ins in the
submitted files. ``--dry-run-synthetic`` runs the whole analysis on a fabricated
TMS-shaped dataset so the code path can be checked offline. Its output is labelled
SYNTHETIC and says nothing about TMS.

What it answers. Each answer goes to ``audit_report.md`` and ``audit_results.json``:

  A. the raw ``development_stage`` values, with cells and unique mice (``donor_id``) per
     sex x stage;
  B. how ``download_data.py``'s age regex maps those strings (does "20m" merge several
     ages?), and the age encoded in the TMS mouse id;
  C. whether ERCC spike-ins are present: in the Census var, and in the source .h5ad;
  D. §4 as published, recomputed with the v0.1.1 engine at the cell level;
  E. §4 with the mouse as the unit: per-mouse QC and mi_3bin, an exact permutation over
     mice within sex, and tissue composition per group;
  F. same-data positive controls. Xist and the Y genes, female vs male within age: the
     v0.1.1 verdict, per-mouse values, and a per-mouse sex-label sanity check.

It writes no conclusions about §4 beyond the numbers. Which claims still hold is for
the human-written report.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import importlib.util
import itertools
import json
import platform
import re
import subprocess
import sys
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
PINNED_ENGINE = "0.1.1"
DOWNLOAD_MODULE = REPO / "examples" / "mi_coupling_tms" / "download_data.py"

PAIR = ("Smad3", "Col1a1")
POS_PAIR = ("Actb", "Gapdh")
NEG_PAIR = ("Tst", "Lrrc42")               # the expression-matched negative control in §4
SEX_GENES = ("Xist", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty")
GENES = PAIR + POS_PAIR + NEG_PAIR + SEX_GENES
Y_GENES = SEX_GENES[1:]

# Numbers as stated in preprint §4 (paper/manuscript.md), for side-by-side comparison.
PAPER = {
    "n_cells": 110824,
    "median_nnz": {("male", "3m"): 2799, ("male", "18m"): 1882, ("male", "20m"): 1701,
                   ("female", "3m"): 2495, ("female", "18m"): 2291, ("female", "20m"): 3413},
    "n_cells_female_20m": 728,
    "gate1_male_ratio": 1.65,
    "gate0_dropout_rel": 0.42, "gate0_dropout_z": 63.2,
    "gate2_pooled_retained": 0.43, "gate2_male_retained": 0.30,
    "neg_control_mi": 0.038,
    "zero_frac_young_old": (0.92, 0.97),
}


# --------------------------------------------------------------------------- #
# loading
# --------------------------------------------------------------------------- #
def _load_download_module():
    """examples/mi_coupling_tms/download_data.py, imported by path (it is not a package)."""
    spec = importlib.util.spec_from_file_location("ma_download_data", DOWNLOAD_MODULE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_census(census_version: str, cache: Path):
    """obs for every TMS FACS Smart-seq2 primary cell + the audit genes, cached locally."""
    import anndata
    import cellxgene_census

    dl = _load_download_module()
    cache.mkdir(parents=True, exist_ok=True)
    h5 = cache / f"tms_facs_audit_genes_{census_version}.h5ad"
    meta_path = cache / f"census_meta_{census_version}.json"
    if h5.exists() and meta_path.exists():
        return anndata.read_h5ad(h5), json.loads(meta_path.read_text())

    obs_cols = ["assay", "dataset_id", "cell_type", "tissue", "tissue_general", "sex",
                "development_stage", "development_stage_ontology_term_id", "donor_id",
                "is_primary_data", "nnz", "raw_sum"]
    with cellxgene_census.open_soma(census_version=census_version) as census:
        summary = census["census_info"]["summary"].read().concat().to_pandas()
        ids = dl._resolve_tms_dataset_ids(census)
        datasets = census["census_info"]["datasets"].read().concat().to_pandas()
        id_list = ", ".join(f"'{i}'" for i in ids)
        value_filter = (f"dataset_id in [{id_list}] and assay == 'Smart-seq2' "
                        "and is_primary_data == True")
        gene_list = ", ".join(f"'{g}'" for g in GENES)
        adata = cellxgene_census.get_anndata(
            census, organism="Mus musculus", obs_value_filter=value_filter,
            var_value_filter=f"feature_name in [{gene_list}]",
            column_names={"obs": obs_cols, "var": ["feature_id", "feature_name"]},
        )
        var = census["census_data"]["mus_musculus"].ms["RNA"].var.read(
            column_names=["feature_id", "feature_name"]).concat().to_pandas()
    meta = {
        "census_version": census_version,
        "census_summary": dict(zip(summary["label"].astype(str), summary["value"].astype(str))),
        "tms_dataset_ids": list(ids),
        "tms_dataset_titles": datasets.set_index("dataset_id").loc[list(ids), "dataset_title"].tolist(),
        "census_var_n": int(len(var)),
        "census_var_ercc": int(var["feature_id"].astype(str).str.startswith("ERCC").sum()
                               + var["feature_name"].astype(str).str.startswith("ERCC").sum()),
    }
    # Census var_names default to a RangeIndex; gene symbols are needed for name lookup.
    # (download_data.py does not set them, so --gene-a Smad3 would not resolve on its output.)
    meta["census_var_names_were_symbols"] = bool(
        set(map(str, adata.var_names)) & set(GENES))
    adata.var_names = adata.var["feature_name"].astype(str).values
    adata.obs_names = adata.obs_names.astype(str)
    adata.write_h5ad(h5)
    meta_path.write_text(json.dumps(meta, indent=2))
    return adata, meta


def check_source_h5ad(dataset_ids, census_version: str, cache: Path) -> list[dict]:
    """Download each source .h5ad and look for ERCC spike-ins (feature ids/names, biotype)."""
    import anndata
    import cellxgene_census

    out = []
    for ds in dataset_ids:
        path = cache / f"source_{ds}.h5ad"
        if not path.exists():
            cellxgene_census.download_source_h5ad(ds, to_path=str(path), census_version=census_version)
        a = anndata.read_h5ad(path, backed="r")
        var = a.var.copy()
        names = pd.Series(a.var_names.astype(str))
        is_ercc = names.str.startswith("ERCC").values
        for col in ("feature_name", "gene_symbols", "feature_id"):
            if col in var:
                is_ercc |= var[col].astype(str).str.startswith("ERCC").values
        if "feature_biotype" in var:
            is_ercc |= (var["feature_biotype"].astype(str) == "spike-in").values
        rec = dict(dataset_id=ds, n_vars=int(a.n_vars), n_ercc_vars=int(is_ercc.sum()),
                   has_raw=a.raw is not None)
        if is_ercc.any():
            # per-cell ERCC totals in row chunks (backed mode)
            cols = np.where(is_ercc)[0]
            totals = []
            for start in range(0, a.n_obs, 20000):
                chunk = a.X[start:start + 20000]
                chunk = chunk[:, cols]
                totals.append(np.asarray(chunk.sum(axis=1)).ravel())
            tot = np.concatenate(totals)
            obs = a.obs[[c for c in ("donor_id", "sex", "development_stage") if c in a.obs]].copy()
            obs["ercc_total"] = tot
            rec["frac_cells_with_ercc"] = float(np.mean(tot > 0))
            rec["ercc_by_group"] = (obs.groupby([c for c in ("sex", "development_stage") if c in obs],
                                                observed=True)["ercc_total"]
                                    .median().reset_index().to_dict("records"))
        a.file.close()
        out.append(rec)
    return out


# --------------------------------------------------------------------------- #
# analysis (engine-version pinned; works on any AnnData-like object)
# --------------------------------------------------------------------------- #
def donor_age_months(donor_id: str):
    """TMS mouse ids start with the age in months (e.g. '24_58_M', '3-M-8')."""
    m = re.match(r"\s*(\d+)", str(donor_id))
    return int(m.group(1)) if m else None


def section_a_b(obs: pd.DataFrame, age_token) -> dict:
    obs = obs.copy()
    obs["development_stage"] = obs["development_stage"].astype(str)
    obs["donor_id"] = obs["donor_id"].astype(str)
    stage = (obs.groupby(["development_stage", "sex"], observed=True)
             .agg(n_cells=("donor_id", "size"), n_mice=("donor_id", "nunique")).reset_index())
    mapping = {s: age_token(s) for s in sorted(obs["development_stage"].unique())}
    inverse = {}
    for raw, tok in mapping.items():
        inverse.setdefault(tok, []).append(raw)
    collisions = {tok: raws for tok, raws in inverse.items() if len(raws) > 1}
    obs["age"] = obs["development_stage"].map(mapping)
    obs["donor_age_months"] = obs["donor_id"].map(donor_age_months)
    donors = (obs.groupby("donor_id", observed=True)
              .agg(sex=("sex", "first"), age_token=("age", "first"),
                   stages=("development_stage", lambda s: sorted(set(s))),
                   donor_age_months=("donor_age_months", "first"),
                   n_cells=("sex", "size"), n_tissues=("tissue", "nunique"))
              .reset_index())
    tok_months = donors["age_token"].str.extract(r"^(\d+)m$")[0].astype(float)
    donors["token_matches_donor_age"] = tok_months == donors["donor_age_months"]
    groups = (donors.groupby(["sex", "age_token"], observed=True)
              .agg(n_mice=("donor_id", "size"), n_cells=("n_cells", "sum"),
                   mice=("donor_id", lambda s: sorted(s)),
                   donor_ages=("donor_age_months", lambda s: sorted(set(s.dropna().astype(int)))))
              .reset_index())
    return dict(stage_table=stage, regex_mapping=mapping, regex_collisions=collisions,
                donors=donors, groups=groups, obs=obs)


def _dense(x):
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x, dtype=float)


def _per_cell_lognorm(data, genes):
    """log1p(CP10k) using the library size from obs (the object holds only the audit genes).

    Caveat: v0.1.1's GATE 0 perturbs X but not obs, so on this gene-subset object GATE 0 is
    not meaningful for library-normalized metrics; GATE 1 (data-level) is what section F
    needs.
    """
    from metric_autopsy.core import gene_column
    lib = np.asarray(data.obs["total_counts"], dtype=float)
    lib[lib <= 0] = 1.0
    X = np.column_stack([gene_column(data, g) for g in genes])
    return np.log1p(X / lib[:, None] * 1e4)


def make_sex_metrics():
    def mean_lognorm_xist(data):
        return float(np.mean(_per_cell_lognorm(data, ["Xist"])))

    def y_score(data):
        return float(np.mean(_per_cell_lognorm(data, Y_GENES)))
    return mean_lognorm_xist, y_score


def exact_permutation_p(a, b):
    """Two-sided exact permutation p for a difference in means (mice as units)."""
    a, b = list(a), list(b)
    allv = np.array(a + b, dtype=float)
    n, na = len(allv), len(a)
    if na == 0 or n - na == 0:
        return float("nan"), float("nan")
    obs = abs(np.mean(a) - np.mean(b))
    diffs = []
    for comb in itertools.combinations(range(n), na):
        mask = np.zeros(n, bool)
        mask[list(comb)] = True
        diffs.append(abs(allv[mask].mean() - allv[~mask].mean()))
    diffs = np.array(diffs)
    total = len(diffs)
    p = float(np.mean(diffs >= obs - 1e-12))
    # equal group sizes: each split and its mirror give the same |difference|
    return p, (2.0 if na == n - na else 1.0) / total


def section_d(adata, groups) -> dict:
    """§4 recomputed at the cell level with the pinned engine."""
    from metric_autopsy import metrics, run_autopsy, gate2_ngenes_matching, gate3_raw_visibility

    sub = adata[np.isin(np.asarray(adata.obs["age"]), groups)].copy()
    mi = partial(metrics.mi_3bin, gene_a=PAIR[0], gene_b=PAIR[1])
    a = run_autopsy(mi, sub, group_col="age", groups=tuple(groups), within=["sex"],
                    gene_pair=PAIR, pair_metric=metrics.mi_3bin, pos_pair=POS_PAIR,
                    neg_pair=NEG_PAIR, stop_on_first_fail=False)
    res = {r.gate: r for r in a.results}
    g0 = res[0].detail.get("responses", {}).get("extra_dropout", {})
    g1_rows = res[1].detail.get("table", [])
    male = sub[np.asarray(sub.obs["sex"]) == "male"].copy()
    g2_male = gate2_ngenes_matching(mi, male, "age", tuple(groups))
    g3 = gate3_raw_visibility(sub, PAIR[0], PAIR[1], "age", tuple(groups))
    neg_val = float(metrics.mi_3bin(sub, gene_a=NEG_PAIR[0], gene_b=NEG_PAIR[1]))
    return dict(
        verdict=a.verdict, markdown=a.to_markdown(),
        gate0_dropout_rel=g0.get("rel_change"), gate0_dropout_z=g0.get("z"),
        gate1_table=g1_rows,
        gate2_pooled=res[2].detail, gate2_pooled_message=res[2].message,
        gate2_male=g2_male.detail, gate2_male_message=g2_male.message,
        gate3_zero_frac=(g3.detail["group_a"]["zero_frac"], g3.detail["group_b"]["zero_frac"]),
        neg_control_mi=neg_val,
        median_nnz=(sub.obs.groupby(["sex", "age"], observed=True)["n_genes_by_counts"]
                    .median().to_dict()),
    )


def section_e(adata, groups) -> dict:
    """Mouse as the unit: per-mouse QC and mi_3bin; exact permutation within sex."""
    from metric_autopsy import metrics

    sub = adata[np.isin(np.asarray(adata.obs["age"]), groups)].copy()
    rows = []
    for mouse in pd.unique(np.asarray(sub.obs["donor_id"]).astype(str)):
        m = sub[np.asarray(sub.obs["donor_id"]).astype(str) == mouse]
        rows.append(dict(
            donor_id=mouse, sex=str(m.obs["sex"].iloc[0]), age=str(m.obs["age"].iloc[0]),
            donor_age_months=donor_age_months(mouse), n_cells=int(m.n_obs),
            n_tissues=int(m.obs["tissue"].nunique()),
            median_nnz=float(np.median(m.obs["n_genes_by_counts"])),
            mi_3bin=float(metrics.mi_3bin(m, gene_a=PAIR[0], gene_b=PAIR[1])),
        ))
    mice = pd.DataFrame(rows)
    tests = {}
    for sex in sorted(mice["sex"].unique()):
        ms = mice[mice["sex"] == sex]
        young = ms[ms["age"] == groups[0]]
        old = ms[ms["age"] == groups[1]]
        for col in ("median_nnz", "mi_3bin"):
            p, pmin = exact_permutation_p(young[col], old[col])
            tests[f"{sex}:{col}"] = dict(
                n_young=int(len(young)), n_old=int(len(old)),
                mean_young=float(young[col].mean()) if len(young) else None,
                mean_old=float(old[col].mean()) if len(old) else None,
                p_exact=p, min_attainable_p=pmin)
    tissue = pd.crosstab(np.asarray(sub.obs["tissue"]).astype(str),
                         [np.asarray(sub.obs["sex"]).astype(str), np.asarray(sub.obs["age"]).astype(str)],
                         rownames=["tissue"], colnames=["sex", "age"])
    tissue.columns = [f"{a}:{b}" for a, b in tissue.columns]
    return dict(mice=mice, tests=tests, tissue_cells=tissue)


def section_f(adata) -> dict:
    """Xist / Y genes, female vs male within age: v0.1.1 verdict + per-mouse sanity check."""
    from metric_autopsy import gate1_qc_parity, gate2_ngenes_matching

    mean_lognorm_xist, y_score = make_sex_metrics()
    sub = adata[np.isin(np.asarray(adata.obs["sex"]).astype(str), ["female", "male"])].copy()
    out = {}
    # GATE 0 is not evaluable here: v0.1.1 perturbs X but not the obs library size this
    # gene-subset object relies on. GATE 1 is data-level (metric-independent), GATE 2 uses
    # the metric on unperturbed data; the v0.1.1 verdict is the first blocking gate among them.
    g1 = gate1_qc_parity(sub, "sex", ("female", "male"), within=["age"])
    for name, metric in (("Xist", mean_lognorm_xist), ("Y_score", y_score)):
        g2 = gate2_ngenes_matching(metric, sub, "sex", ("female", "male"))
        first = next((g for g in (g1, g2) if g.blocking), None)
        verdict = (f"FAIL — died at GATE {first.gate} ({first.name})" if first
                   else "no blocking gate among GATES 1-2")
        markdown = "\n".join([
            "| Gate | Status | Finding |", "|---|---|---|",
            "| 0 | not evaluable | gene-subset object: perturbations do not update the obs library size |",
            f"| 1 | {g1.status.value} | {g1.message} |",
            f"| 2 | {g2.status.value} | {g2.message} |", "",
            f"**v0.1.1 verdict (GATES 1-2): {verdict}**"])
        out[name] = dict(verdict=verdict, markdown=markdown, gate1_table=g1.detail.get("table", []),
                         gate2=g2.detail)
    rows = []
    for mouse in pd.unique(np.asarray(sub.obs["donor_id"]).astype(str)):
        m = sub[np.asarray(sub.obs["donor_id"]).astype(str) == mouse]
        rows.append(dict(donor_id=mouse, sex=str(m.obs["sex"].iloc[0]), age=str(m.obs["age"].iloc[0]),
                         n_cells=int(m.n_obs), xist=mean_lognorm_xist(m), y_score=y_score(m)))
    mice = pd.DataFrame(rows)
    # sex-label sanity: a female mouse must have Xist >> Y, a male the opposite
    mice["label_consistent"] = np.where(mice["sex"] == "female",
                                        mice["xist"] > mice["y_score"], mice["y_score"] > mice["xist"])
    out["mice"] = mice
    return out


# --------------------------------------------------------------------------- #
# reporting
# --------------------------------------------------------------------------- #
def _md_table(df: pd.DataFrame, max_rows=60) -> str:
    df = df.head(max_rows)
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        cells = [f"{v:.4g}" if isinstance(v, float) else str(v) for v in r.values]
        lines.append("| " + " | ".join(c.replace("|", "\\|") for c in cells) + " |")
    return "\n".join(lines)


def _jsonable(o):
    if isinstance(o, pd.DataFrame):
        return o.to_dict("records")
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if hasattr(o, "value") and hasattr(o, "name"):  # enums
        return o.value
    return o


def write_report(out: Path, meta, ab, c, d, e, f, groups, synthetic: bool):
    out.mkdir(parents=True, exist_ok=True)
    banner = ("> **SYNTHETIC DRY RUN — fabricated TMS-shaped data. Nothing here describes Tabula "
              "Muris Senis.**\n\n" if synthetic else "")
    L = [f"# Flagship-data audit (§4, TMS FACS) — engine v{meta.get('engine_version')}", "", banner]
    L += ["## Provenance", "", "```json", json.dumps(_jsonable(meta), indent=2), "```", ""]
    L += ["## A. development_stage × sex (cells, unique mice)", "", _md_table(ab["stage_table"]), ""]
    L += ["## B. Age regex (download_data.py) and mouse ids", "",
          "Mapping raw → token:", "", "```json", json.dumps(ab["regex_mapping"], indent=2), "```", "",
          f"Collisions (one token, several raw strings): `{json.dumps(ab['regex_collisions'])}`", "",
          "Groups as the §4 analysis sees them (sex × token):", "",
          _md_table(ab["groups"][["sex", "age_token", "n_mice", "n_cells", "donor_ages", "mice"]]), "",
          "Mice whose id-encoded age disagrees with their token:", "",
          _md_table(ab["donors"][~ab["donors"]["token_matches_donor_age"]]
                    [["donor_id", "sex", "age_token", "donor_age_months", "stages", "n_cells"]]), ""]
    L += ["## C. ERCC spike-ins", "", "```json", json.dumps(_jsonable(c), indent=2), "```", ""]
    L += [f"## D. §4 recomputed with the pinned engine (cell level; groups {groups})", "",
          d["markdown"], "",
          "| quantity | §4 says | recomputed |", "|---|---|---|",
          f"| GATE 0 dropout rel / z | {PAPER['gate0_dropout_rel']} / {PAPER['gate0_dropout_z']} | "
          f"{d['gate0_dropout_rel']} / {d['gate0_dropout_z']} |",
          f"| GATE 2 pooled retained | {PAPER['gate2_pooled_retained']} | "
          f"{d['gate2_pooled'].get('retained_frac')} ({d['gate2_pooled_message']}) |",
          f"| GATE 2 male retained | {PAPER['gate2_male_retained']} | "
          f"{d['gate2_male'].get('retained_frac')} ({d['gate2_male_message']}) |",
          f"| negative control MI (Tst-Lrrc42) | {PAPER['neg_control_mi']} | {d['neg_control_mi']:.4g} |",
          f"| zero-fraction young/old | {PAPER['zero_frac_young_old']} | "
          f"{tuple(round(x, 3) for x in d['gate3_zero_frac'])} |", "",
          "Median nnz by sex × age (cells): " + json.dumps({f"{k[0]}:{k[1]}": v for k, v in
                                                            d["median_nnz"].items()}), "",
          "GATE 1 table:", "", _md_table(pd.DataFrame(d["gate1_table"])), ""]
    L += ["## E. Mouse as the unit", "", _md_table(e["mice"]), "",
          "Exact permutation over mice within sex (young vs old):", "",
          _md_table(pd.DataFrame(e["tests"]).T.reset_index().rename(columns={"index": "test"})), "",
          "Cells per tissue × (sex, age):", "", _md_table(e["tissue_cells"].reset_index()), ""]
    L += ["## F. Same-data positive controls (female vs male, within age)", ""]
    for name in ("Xist", "Y_score"):
        L += [f"### {name}", "", f[name]["markdown"], ""]
    L += ["Per-mouse sex markers:", "", _md_table(f["mice"]), ""]
    (out / "audit_report.md").write_text("\n".join(L))
    (out / "audit_results.json").write_text(json.dumps(_jsonable(dict(
        meta=meta, A=ab["stage_table"], B=dict(mapping=ab["regex_mapping"],
                                                collisions=ab["regex_collisions"],
                                                donors=ab["donors"], groups=ab["groups"]),
        C=c, D={k: v for k, v in d.items() if k != "markdown"},
        E=dict(mice=e["mice"], tests=e["tests"]),
        F={k: (v if k == "mice" else {kk: vv for kk, vv in v.items() if kk != "markdown"})
           for k, v in f.items()},
        paper=dict((k if not isinstance(k, tuple) else ":".join(k), v) for k, v in PAPER.items()),
    )), indent=2, default=str))


# --------------------------------------------------------------------------- #
# provenance (self-contained: the pinned v0.1.1 engine has no provenance module)
# --------------------------------------------------------------------------- #
OBS_HASH_COLS = ("donor_id", "sex", "development_stage", "tissue", "assay", "age",
                 "n_genes_by_counts", "total_counts")


def _sha256_file(path: Path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_anndata(adata, obs_columns=OBS_HASH_COLS) -> str:
    """Content hash of X (sparse or dense), var_names and the obs columns the audit reads —
    the same recipe as metric_autopsy.provenance.sha256_data in v0.3."""
    h = hashlib.sha256()
    X = adata.X
    if hasattr(X, "tocsr"):
        csr = X.tocsr()
        h.update(f"csr{csr.shape}{csr.dtype}".encode())
        for arr in (csr.indptr, csr.indices, csr.data):
            h.update(np.ascontiguousarray(arr).tobytes())
    else:
        arr = np.ascontiguousarray(np.asarray(X))
        h.update(f"dense{arr.shape}{arr.dtype}".encode())
        h.update(arr.tobytes())
    h.update("\x1f".join(map(str, adata.var_names)).encode())
    for col in sorted(c for c in obs_columns if c in adata.obs.columns):
        h.update(col.encode() + b"\x1e")
        h.update("\x1f".join(map(str, np.asarray(adata.obs[col]))).encode())
    return h.hexdigest()


def provenance(adata, args, download_module_path, cached_h5ad=None) -> dict:
    def git(*cmd):
        try:
            return subprocess.run(["git", *cmd], cwd=REPO, capture_output=True, text=True).stdout.strip()
        except Exception:
            return None
    versions = dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__)
    for mod in ("anndata", "scipy", "cellxgene_census"):
        try:
            from importlib import metadata
            versions[mod] = metadata.version(mod.replace("_", "-"))
        except Exception:
            versions[mod] = None
    import metric_autopsy  # the imported module, not the (possibly stale) editable-install metadata
    versions["metric_autopsy"] = metric_autopsy.__version__
    params = {k: v for k, v in vars(args).items() if k not in ("out", "cache")}
    return dict(
        timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        git_commit=git("rev-parse", "HEAD"), git_dirty=bool(git("status", "--porcelain")),
        script_sha256=_sha256_file(Path(__file__)),
        download_module_sha256=_sha256_file(download_module_path),
        data_sha256=sha256_anndata(adata), data_hash_obs_columns=list(OBS_HASH_COLS),
        cached_h5ad_sha256=_sha256_file(cached_h5ad) if cached_h5ad and cached_h5ad.exists() else None,
        params=params, params_sha256=hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest(),
        versions=versions,
    )


# --------------------------------------------------------------------------- #
# synthetic dry run (offline check of the code path; NOT TMS)
# --------------------------------------------------------------------------- #
def synthetic_tms(seed=0):
    import anndata
    rng = np.random.default_rng(seed)
    design = [("3_1_M", "male", "3-month-old stage"), ("3_2_M", "male", "3-month-old stage"),
              ("3_3_F", "female", "3-month-old stage"), ("3_4_F", "female", "3-month-old stage"),
              ("18_5_M", "male", "18-month-old stage"), ("18_6_F", "female", "18-month-old stage"),
              ("21_7_F", "female", "20 month-old stage and over"),
              ("24_8_M", "male", "20 month-old stage and over"),
              ("24_9_M", "male", "20 month-old stage and over")]
    tissues = ["Lung", "Liver", "Heart", "Marrow"]
    obs_rows, Xs = [], []
    for donor, sex, stage in design:
        n = int(rng.integers(150, 300))
        eff = 0.5 if donor.startswith("24") else 1.0
        lam = np.exp(rng.normal(0.5, 0.8, (n, len(GENES)))) * eff
        if sex == "female":
            lam[:, GENES.index("Xist")] *= 30
            lam[:, [GENES.index(g) for g in Y_GENES]] *= 0.0
        else:
            lam[:, GENES.index("Xist")] *= 0.0
            lam[:, [GENES.index(g) for g in Y_GENES]] *= 8
        Xs.append(rng.poisson(lam).astype(float))
        nnz = rng.normal(2500 * eff, 400, n).clip(200)
        for k in range(n):
            obs_rows.append(dict(donor_id=donor, sex=sex, development_stage=stage,
                                 tissue=tissues[k % len(tissues)], assay="Smart-seq2",
                                 nnz=nnz[k], raw_sum=nnz[k] * 300))
    obs = pd.DataFrame(obs_rows)
    obs.index = [f"cell{i}" for i in range(len(obs))]
    ad = anndata.AnnData(X=np.vstack(Xs), obs=obs, var=pd.DataFrame(index=list(GENES)))
    meta = dict(census_version="SYNTHETIC", census_var_ercc=0, tms_dataset_ids=["synthetic"])
    return ad, meta


# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "out"))
    ap.add_argument("--cache", default=str(REPO / "data" / "raw" / "audit_cache"))
    ap.add_argument("--census-version", default="stable")
    ap.add_argument("--young", default="3m")
    ap.add_argument("--old", default="20m", help="the token §4 used for 'old'")
    ap.add_argument("--source-h5ad", action="store_true",
                    help="also download the source .h5ad files and look for ERCC spike-ins")
    ap.add_argument("--dry-run-synthetic", action="store_true")
    ap.add_argument("--allow-engine-version", action="store_true",
                    help=f"run even if metric_autopsy is not v{PINNED_ENGINE}")
    args = ap.parse_args(argv)

    import metric_autopsy
    if metric_autopsy.__version__ != PINNED_ENGINE and not args.allow_engine_version:
        sys.exit(f"metric_autopsy {metric_autopsy.__version__} is installed; this audit is pinned to "
                 f"v{PINNED_ENGINE} (see the module docstring for the worktree recipe)")

    dl = _load_download_module()
    if args.dry_run_synthetic:
        adata, meta = synthetic_tms()
        out = Path(args.out) / "DRY_RUN_SYNTHETIC"
    else:
        adata, meta = load_census(args.census_version, Path(args.cache))
        out = Path(args.out)
    meta["engine_version"] = metric_autopsy.__version__
    meta["n_cells"] = int(adata.n_obs)

    # the same QC-column mapping download_data.py applies (nnz -> n_genes_by_counts, ...)
    adata = dl._make_metric_autopsy_ready(adata)
    ab = section_a_b(adata.obs, _age_token_from(dl))
    adata.obs["age"] = ab["obs"]["age"].values
    groups = (args.young, args.old)

    c = dict(census_var_ercc=meta.get("census_var_ercc"))
    if args.source_h5ad and not args.dry_run_synthetic:
        c["source_h5ad"] = check_source_h5ad(meta["tms_dataset_ids"], args.census_version, Path(args.cache))

    cached = (None if args.dry_run_synthetic else
              Path(args.cache) / f"tms_facs_audit_genes_{args.census_version}.h5ad")
    meta["provenance"] = provenance(adata, args, DOWNLOAD_MODULE, cached)

    d = section_d(adata, groups)
    e = section_e(adata, groups)
    f = section_f(adata)
    write_report(out, meta, ab, c, d, e, f, groups, synthetic=args.dry_run_synthetic)
    print(f"wrote {out / 'audit_report.md'} and {out / 'audit_results.json'}")


def _age_token_from(dl):
    """download_data.py defines _age_token inside _make_metric_autopsy_ready; rebuild the same
    regex here (kept byte-identical to download_data.py:62-65) and assert it matches."""
    def _age_token(s):
        s = str(s)
        m = re.search(r"(\d+)\s*-?\s*(month|year)", s)
        return f"{m.group(1)}{'m' if m.group(2) == 'month' else 'y'}" if m else s
    import inspect
    src = inspect.getsource(dl._make_metric_autopsy_ready)
    assert r're.search(r"(\d+)\s*-?\s*(month|year)", s)' in src, "download_data.py age regex changed"
    return _age_token


if __name__ == "__main__":
    main()
