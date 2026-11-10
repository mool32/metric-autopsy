#!/usr/bin/env python3
"""Level B, step S6 (`plan_s6.md`): the controls, selected by the rules registered before any was selected.

Run once in GitHub Actions (`.github/workflows/level-b-controls.yml`, the label `level-b-controls`); this
environment does not read the Census. Nothing of the engine is imported or run: the rules count cells and
detections, and the negative controls' claims come from v1's rule of pairs (`panel.plan_background`, unchanged).

* **Negative controls** (`plan_s6.md` section 1): the CELLxGENE Census's `stable` release, primary cells of
  tissues (`tissue_type == 'tissue'`, `JOURNAL.md` J5), human and mouse. A candidate is a collection's
  tissue and cell type in two technical conditions (assay | suspension type) with at least `MIN_DONORS`
  donors of at least `MIN_CELLS` cells in each, a donor's disease, sex and stage the same in both; ranked by
  its donors, then the median of its cells per donor and condition, then its ids; one per collection (its
  first), the first `N_NEGATIVE` collections whose extraction has raw counts and v1's pool of pairs, over at
  most `MAX_TRIED` collections. The claims: every pair of the pool, a set per level of the high and medium
  levels' pairs, each claim's direction by a fair coin (`SEED_DIRECTIONS`, keyed by the control's number).
* **The positive control** (section 2): Kang et al. 2018, GEO GSE96583's second batch, fetched here; its
  structure checked; in CD14+ monocytes the twelve interferon-stimulated genes `ISGS`, eligible where detected
  in at least 10% of the cells of each condition and induced in every donor; a seeded matching into disjoint
  pairs (`SEED_MATCHING`) and one set.

    python validation/level_b/controls.py --out OUT              # the run
    python validation/level_b/controls.py --out OUT --synthetic  # the dry run: a fabricated Census and GEO

OUT gets controls.log (the report), selection.json, candidates.tsv, tried.tsv, claims.json, joinids/,
kang/ (the files' requests and the structure check), and data/ (the extractions, as npz; the workflow keeps
them as an artifact, the branch gets their SHA256). With --synthetic every output says SYNTHETIC.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import platform
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "prereg"))
import panel as P  # noqa: E402  v1's rule of pairs, npz I/O and gene names: numpy and pandas, no engine
from select_backgrounds import _known, count_cells, raw_counts  # noqa: E402  v1's selection helpers

# --------------------------------------------------------------------------- #
# the plan's numbers (plan_s6.md, committed in e857481 before any selection)
# --------------------------------------------------------------------------- #
MIN_DONORS = 6                    # pairs: the engine's exact paired test reaches p < 0.05 at 6 (2/2^6)
MIN_CELLS = P.CELLS_PER_DONOR     # 200, v1's cells per donor, in each condition
MAX_DONORS, MAX_CELLS = 200, 400  # v1's extraction caps (select_backgrounds.MAX_DONORS, MAX_CELLS)
MAX_TRIED = 10                    # collections tried for raw counts and the pool
N_NEGATIVE = 3
SEED_CELLS, SEED_DIRECTIONS, SEED_MATCHING = 20261110, 20261111, 20261112
SET_LEVELS = ("high", "medium")
ISGS = ("ISG15", "IFI6", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "OAS1", "OAS2", "RSAD2", "IFI44L", "ISG20")
DOMAIN_DETECTION = 0.10           # the engine's domain rule for a pair (metric_autopsy.domain; a test checks it)
KANG_TYPE = "CD14+ Monocytes"
MIN_KANG_DONORS = 6
SEX_GENES = {"homo_sapiens": ("XIST", ("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "ZFY", "USP9Y")),
             "mus_musculus": ("Xist", ("Ddx3y", "Eif2s3y", "Kdm5d", "Uty"))}
SEX_SHARE, SEX_RATIO = 0.05, 10.0  # a sex is called where the other's genes are in <= 5% of the cells and its
#                                    own in at least 10 times that share and 5% (the real pilot's contrast call)
ORGANISMS = {"homo_sapiens": "Homo sapiens", "mus_musculus": "Mus musculus"}
BASE_FILTER = "is_primary_data == True and tissue_type == 'tissue'"
COUNTED = ["dataset_id", "tissue", "cell_type", "donor_id", "assay", "suspension_type", "disease", "sex",
           "development_stage"]
# not data the engine was developed or validated on, nor the positive control (plan_s6.md section 1)
EXCLUDED_DATASETS = {"37a17b78-4864-4a42-b67b-31c00962795a": "v1's background B1",
                     "49e4ffcc-5444-406d-bdee-577127404ba8": "v1's background B2"}
EXCLUDED_NAMES = {"tabula muris": "Tabula Muris (Senis): the real pilot",
                  "sea-ad": "SEA-AD: the real pilot",
                  "multiplexed droplet single-cell rna-sequencing using natural genetic variation":
                      "Kang et al. 2018: the positive control"}
EXCLUDED_DOIS = {"10.1038/s41586-020-2496-1": "Tabula Muris Senis: the real pilot",
                 "10.1038/s41593-024-01774-5": "SEA-AD: the real pilot",
                 "10.1038/s41593-022-01022-8": "Marsh et al. 2022: the real pilot",
                 "10.1038/nbt.4042": "Kang et al. 2018: the positive control"}
GEO = "https://ftp.ncbi.nlm.nih.gov/geo"
KANG_FILES = {  # name: (url, what)
    "filelist.txt": (f"{GEO}/series/GSE96nnn/GSE96583/suppl/filelist.txt", "the series' file list"),
    "genes.tsv.gz": (f"{GEO}/series/GSE96nnn/GSE96583/suppl/GSE96583_batch2.genes.tsv.gz", "batch 2's genes"),
    "cells.tsv.gz": (f"{GEO}/series/GSE96nnn/GSE96583/suppl/GSE96583_batch2.total.tsne.df.tsv.gz",
                     "the authors' table of batch 2's cells"),
    "ctrl.mtx.gz": (f"{GEO}/samples/GSM2560nnn/GSM2560248/suppl/GSM2560248_2.1.mtx.gz", "lane 2.1's counts"),
    "ctrl.barcodes.tsv.gz": (f"{GEO}/samples/GSM2560nnn/GSM2560248/suppl/GSM2560248_barcodes.tsv.gz",
                             "lane 2.1's barcodes"),
    "stim.mtx.gz": (f"{GEO}/samples/GSM2560nnn/GSM2560249/suppl/GSM2560249_2.2.mtx.gz", "lane 2.2's counts"),
    "stim.barcodes.tsv.gz": (f"{GEO}/samples/GSM2560nnn/GSM2560249/suppl/GSM2560249_barcodes.tsv.gz",
                             "lane 2.2's barcodes"),
}
USER_AGENT = "metric-autopsy-validation/0.4 (level B controls; python urllib)"


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def units_ok(donor: pd.Series) -> np.ndarray:
    """A donor id that names one donor: present (v1's placeholders excluded) and not composite ('/')."""
    return _known(donor) & ~donor.astype(str).str.contains("/", regex=False).to_numpy()


# --------------------------------------------------------------------------- #
# the negative controls' rule, on cell counts (pure: tested on synthetic tables)
# --------------------------------------------------------------------------- #
def exclusions(datasets: pd.DataFrame) -> dict:
    """{collection_id: why it is excluded} for the collections of EXCLUDED_DATASETS, EXCLUDED_NAMES (in the
    collection's name or a dataset's title) and EXCLUDED_DOIS (the collection's DOI or a citation)."""
    out = {}
    for _, d in datasets.iterrows():
        text = f"{d.get('collection_name', '')} | {d.get('dataset_title', '')}".lower()
        cite = f"{d.get('collection_doi', '')} | {d.get('citation', '')}".lower()
        why = ([EXCLUDED_DATASETS[str(d["dataset_id"])]] if str(d["dataset_id"]) in EXCLUDED_DATASETS else [])
        why += [w for n, w in EXCLUDED_NAMES.items() if n in text] + [w for x, w in EXCLUDED_DOIS.items() if x in cite]
        if why:
            c = str(d["collection_id"])
            out[c] = "; ".join(sorted(set(out.get(c, "").split("; ")) - {""} | set(why)))
    return out


def unit_table(counts: pd.DataFrame, datasets: pd.DataFrame) -> pd.DataFrame:
    """One row per collection, organism, tissue, cell type, condition and donor: its cells, and the commonest
    disease, sex and development stage of them (jointly, by cells; ties by the values in order) with that
    combination's cells (the unit's cells: a donor's cells of another stage are another time point)."""
    coll = dict(zip(datasets["dataset_id"].astype(str), datasets["collection_id"].astype(str)))
    c = counts[counts["dataset_id"].astype(str).isin(coll)].copy()
    c = c[units_ok(c["donor_id"])]
    c["collection_id"] = c["dataset_id"].astype(str).map(coll)
    c["condition"] = c["assay"].astype(str) + " | " + c["suspension_type"].astype(str)
    key = ["organism", "collection_id", "tissue", "cell_type", "condition", "donor_id"]
    meta = ["disease", "sex", "development_stage"]
    per = c.groupby(key + meta, sort=True)["cells"].sum().reset_index()
    per = per.sort_values(key + ["cells"] + meta, ascending=[True] * len(key) + [False] + [True] * 3, kind="mergesort")
    top = per.drop_duplicates(key).rename(columns={"cells": "unit_cells"})
    total = per.groupby(key)["cells"].sum().rename("all_cells").reset_index()
    return top.merge(total, on=key)


def negative_candidates(counts: pd.DataFrame, datasets: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Every qualifying candidate (plan_s6.md section 1), ranked; and its donors. A donor qualifies where its
    units in both conditions have at least MIN_CELLS cells of their commonest disease, sex and stage, and these
    are the same in both. `counts`: cells per COUNTED and organism (strings). The excluded collections are
    kept, with their reason, and ranked after none: `rank` counts the others only."""
    u = unit_table(counts, datasets)
    u = u[u["unit_cells"] >= MIN_CELLS]
    on = ["organism", "collection_id", "tissue", "cell_type", "donor_id"]
    m = u.merge(u, on=on, suffixes=("_a", "_b"))
    m = m[m["condition_a"] < m["condition_b"]].copy()
    m["same"] = ((m["disease_a"] == m["disease_b"]) & (m["sex_a"] == m["sex_b"])
                 & (m["development_stage_a"] == m["development_stage_b"]))
    g = ["organism", "collection_id", "tissue", "cell_type", "condition_a", "condition_b"]
    cols = g + ["donors", "median_cells", "mismatched", "excluded", "rank"]
    if m.empty:
        return pd.DataFrame(columns=cols), m
    rows = []
    for k, grp in m.groupby(g, sort=True):
        ok = grp[grp["same"]]
        cells = np.concatenate([ok["unit_cells_a"].to_numpy(), ok["unit_cells_b"].to_numpy()])
        rows.append(dict(zip(g, k), donors=int(ok["donor_id"].nunique()),
                         median_cells=float(np.median(cells)) if len(cells) else 0.0,
                         mismatched=int((~grp["same"]).sum())))
    cand = pd.DataFrame(rows)
    cand = cand[cand["donors"] >= MIN_DONORS].copy()
    excl = exclusions(datasets)
    cand["excluded"] = cand["collection_id"].map(excl).fillna("")
    cand = cand.sort_values(["donors", "median_cells", *g], ascending=[False, False] + [True] * len(g),
                            kind="mergesort").reset_index(drop=True)
    rank, ranks = 0, []
    for e in cand["excluded"]:
        rank += not e
        ranks.append(rank if not e else 0)
    cand["rank"] = ranks
    return cand[cols], m[m["same"]]


def first_per_collection(cand: pd.DataFrame) -> pd.DataFrame:
    """The ranked candidates, excluded ones left out, one per collection: its first."""
    ok = cand[cand["excluded"] == ""]
    return ok.drop_duplicates("collection_id", keep="first").reset_index(drop=True)


def sample_units(cells: pd.DataFrame, donors: list, seed: int = SEED_CELLS) -> np.ndarray:
    """The joinids to extract: at most MAX_DONORS of the qualifying donors (a sample of their sorted ids where
    there are more), then at most MAX_CELLS cells of each donor and condition (a sample of its sorted joinids
    where there are more), one generator with the public seed, in that order."""
    rng = np.random.default_rng(seed)
    donors = sorted(map(str, donors))
    if len(donors) > MAX_DONORS:
        donors = sorted(rng.choice(donors, size=MAX_DONORS, replace=False).tolist())
    out = []
    for d in donors:
        for cond in sorted(set(cells.loc[cells["donor_id"] == d, "condition"])):
            ids = np.sort(cells.loc[(cells["donor_id"] == d) & (cells["condition"] == cond), "soma_joinid"]
                          .to_numpy(np.int64))
            if len(ids) > MAX_CELLS:
                ids = np.sort(rng.choice(ids, size=MAX_CELLS, replace=False))
            out.append(ids)
    return np.sort(np.concatenate(out)) if out else np.zeros(0, np.int64)


def negative_claims(plan: dict, number: int) -> list:
    """The claims of negative control `number` from v1's plan: every pool pair, then one set per level of
    SET_LEVELS (the genes of its pairs, in the pool's order); each claim's direction by a fair coin from a
    generator keyed by the seed and the control's number, in that order."""
    rng = np.random.default_rng([SEED_DIRECTIONS, number])
    claims = []
    for pe in plan["pool"]:
        claims.append(dict(id=f"N{number}p{pe['index'] + 1:02d}", kind="pair", level=pe["level"], genes=list(pe["pair"]),
                           pooled_r=pe["r"], gate5=dict(negative=list(pe["neg_pair"]), positive=list(pe["pos_pair"]))))
    for lvl in SET_LEVELS:
        genes = [g for pe in plan["pool"] if pe["level"] == lvl for g in pe["pair"]]
        if genes:
            claims.append(dict(id=f"N{number}s{lvl[0]}", kind="set", level=lvl, genes=genes))
    for c in claims:
        c["higher_in"] = "A" if rng.random() < 0.5 else "B"
    return claims


def sex_call(f_x: float, f_y: float) -> str:
    """female, male or unclear from the shares of cells with an X-inactivation count and with a Y-gene count."""
    if f_y <= SEX_SHARE and f_x >= max(SEX_SHARE, SEX_RATIO * f_y):
        return "female"
    if f_x <= SEX_SHARE and f_y >= max(SEX_SHARE, SEX_RATIO * f_x):
        return "male"
    return "unclear"


def sex_shares(X, genes: list, rows: np.ndarray, organism: str) -> tuple[float, float]:
    """The shares of `rows` (cells) with a count of the X-inactivation gene and with a count of any Y gene."""
    xg, yg = SEX_GENES[organism]
    upper = {g.upper(): j for j, g in enumerate(genes)}
    sub = X[rows]
    fx = float((np.asarray(sub[:, upper[xg.upper()]].todense()).ravel() > 0).mean()) if xg.upper() in upper else float("nan")
    ys = [upper[g.upper()] for g in yg if g.upper() in upper]
    fy = (float((np.asarray(sub[:, ys].sum(axis=1)).ravel() > 0).mean()) if ys else float("nan"))
    return fx, fy


# --------------------------------------------------------------------------- #
# the Census, and a fabricated one for the dry run
# --------------------------------------------------------------------------- #
class Census:
    """The CELLxGENE Census's `stable` release (cellxgene_census; validation/prereg/requirements-select.txt)."""

    synthetic = False

    def __init__(self):
        import cellxgene_census
        self._cc = cellxgene_census
        desc = cellxgene_census.get_census_version_description("stable")
        self.release = desc.get("release_build") or "stable"
        self._census = cellxgene_census.open_soma(census_version=self.release)

    def close(self):
        self._census.close()

    def datasets(self) -> pd.DataFrame:
        return self._census["census_info"]["datasets"].read().concat().to_pandas()

    def count(self, organism: str) -> pd.DataFrame:
        obs = self._census["census_data"][organism].obs
        return count_cells(obs.read(value_filter=BASE_FILTER, column_names=COUNTED), COUNTED)

    def cells(self, organism: str, dataset_ids: list, tissue: str, cell_type: str) -> pd.DataFrame:
        flt = (BASE_FILTER + " and dataset_id in [" + ", ".join(repr(str(d)) for d in dataset_ids) + "]"
               + f" and tissue == {tissue!r} and cell_type == {cell_type!r}")
        obs = self._census["census_data"][organism].obs
        return obs.read(value_filter=flt, column_names=["soma_joinid", *COUNTED]).concat().to_pandas()

    def extract(self, organism: str, joinids: np.ndarray) -> tuple:
        from scipy import sparse
        ad = self._cc.get_anndata(self._census, ORGANISMS[organism], X_name="raw", obs_coords=np.asarray(joinids),
                                  obs_column_names=["soma_joinid", *COUNTED],
                                  var_column_names=["soma_joinid", "feature_id", "feature_name"])
        o = np.argsort(ad.obs["soma_joinid"].to_numpy(), kind="mergesort")
        X = sparse.csr_matrix(ad.X)[o]
        return X, ad.obs.iloc[o].reset_index(drop=True), list(ad.var["feature_name"].astype(str))


class SyntheticCensus:
    """A fabricated Census for the dry run (SYNTHETIC, nothing real): cells and nuclei of one human tissue
    from 7 donors (the first choice), 10x 3' v2 and v3 of another from 6 (one of whose donors is relabelled
    in one condition, which leaves 5: it does not qualify), a mouse collection by two chemistries (6 mice),
    a collection the rules exclude by its DOI with more donors than any, and one whose counts are not
    integers (it fails at its extraction). Counts from v1's simulated backgrounds (`simulate.py`): three
    levels of expression with coupled pairs, so that v1's rule of pairs finds its pool."""

    synthetic = True
    release = "SYNTHETIC"

    def __init__(self, seed: int = 11):
        import simulate
        self._sim = simulate
        self.rng = np.random.default_rng(seed)
        self._datasets, self._obs = [], {"homo_sapiens": [], "mus_musculus": []}
        human = "homo_sapiens"
        self._add("heart-cells", "c-heart", "Synthetic heart: cells and nuclei", "10.0000/syn.heart", human,
                  "heart left ventricle", "fibroblast", "10x 3' v3", "cell", [f"H{k}" for k in range(1, 8)], 230)
        self._add("heart-nuclei", "c-heart", "Synthetic heart: cells and nuclei", "10.0000/syn.heart", human,
                  "heart left ventricle", "fibroblast", "10x 3' v3", "nucleus", [f"H{k}" for k in range(1, 8)], 240)
        self._add("lung-v2", "c-lung", "Synthetic lung by two chemistries", "10.0000/syn.lung", human, "lung",
                  "alveolar macrophage", "10x 3' v2", "cell", [f"L{k}" for k in range(1, 7)], 210)
        self._add("lung-v3", "c-lung", "Synthetic lung by two chemistries", "10.0000/syn.lung", human, "lung",
                  "alveolar macrophage", "10x 3' v3", "cell", [f"L{k}" for k in range(1, 7)], 210,
                  stage_of={"L3": "61-year-old stage"})
        self._add("kidney-a", "c-kidney", "Synthetic kidney", "10.0000/syn.kidney", human, "kidney",
                  "kidney epithelial cell", "10x 3' v2", "cell", [f"K{k}" for k in range(1, 7)], 220)
        self._add("kidney-b", "c-kidney", "Synthetic kidney", "10.0000/syn.kidney", human, "kidney",
                  "kidney epithelial cell", "10x 3' v2", "nucleus", [f"K{k}" for k in range(1, 7)], 220)
        self._add("seaad-like", "c-excl", "SEA-AD: a synthetic stand-in", "10.1038/s41593-024-01774-5", human,
                  "middle temporal gyrus", "astrocyte", "10x 3' v3", "nucleus", [f"S{k}" for k in range(1, 13)], 210)
        self._add("seaad-like-2", "c-excl", "SEA-AD: a synthetic stand-in", "10.1038/s41593-024-01774-5", human,
                  "middle temporal gyrus", "astrocyte", "10x multiome", "nucleus", [f"S{k}" for k in range(1, 13)], 210)
        self._add("mouse-v2", "c-mouse", "Synthetic mouse spleen", "10.0000/syn.mouse", "mus_musculus", "spleen",
                  "B cell", "10x 3' v2", "cell", [f"M{k}" for k in range(1, 7)], 205)
        self._add("mouse-v3", "c-mouse", "Synthetic mouse spleen", "10.0000/syn.mouse", "mus_musculus", "spleen",
                  "B cell", "10x 3' v3", "cell", [f"M{k}" for k in range(1, 7)], 205)
        self._datasets = pd.DataFrame(self._datasets)
        self._datasets["soma_joinid"] = np.arange(len(self._datasets))
        self._X, self._genes = {}, {}
        for org, rows in self._obs.items():
            obs = pd.DataFrame(rows)
            obs.insert(0, "soma_joinid", np.arange(len(obs), dtype=np.int64))
            self._obs[org] = obs
            self._X[org], self._genes[org] = self._counts(obs, org)

    def _add(self, ds, coll, name, doi, org, tissue, ct, assay, susp, donors, n, stage_of=None):
        self._datasets.append(dict(dataset_id=ds, collection_id=coll, collection_name=name, collection_doi=doi,
                                   dataset_title=f"{name}: {ds}", citation=f"Publication: https://doi.org/{doi}"))
        for k, d in enumerate(donors):
            sex = "female" if k % 2 == 0 else "male"
            stage = (stage_of or {}).get(d, "60-year-old stage" if org == "homo_sapiens" else "3-month-old stage")
            self._obs[org] += [dict(dataset_id=ds, assay=assay, suspension_type=susp, tissue=tissue, tissue_type="tissue",
                                    cell_type=ct, sex=sex, development_stage=stage, donor_id=d, disease="normal",
                                    is_primary_data=True)] * n
        # a few cells of an organoid of the same donors: never counted (tissue_type)
        self._obs[org] += [dict(dataset_id=ds, assay=assay, suspension_type=susp, tissue=tissue, tissue_type="organoid",
                                cell_type=ct, sex="female", development_stage="unknown", donor_id=donors[0],
                                disease="normal", is_primary_data=True)] * 3

    def _counts(self, obs: pd.DataFrame, org: str):
        from scipy import sparse
        n_genes = 420
        bg = self._sim.simulated_background("syn", donors=1, cells=len(obs), n_genes=n_genes,
                                            seed=int(self.rng.integers(1 << 30)), plan=False)
        X = np.asarray(bg.X, dtype=float)
        genes = list(bg.genes)
        xg, yg = SEX_GENES[org]
        female = (obs["sex"] == "female").to_numpy()
        sexcols = np.zeros((len(obs), 1 + len(yg)))
        sexcols[:, 0] = self.rng.poisson(3.0 * female)
        for j in range(len(yg)):
            sexcols[:, 1 + j] = self.rng.poisson(1.0 * ~female)
        X = np.hstack([X, sexcols])
        genes += [xg, *yg]
        bad = (obs["dataset_id"] == "kidney-a").to_numpy()
        X[bad] = X[bad] * 0.5  # not integers: the kidney fails at its extraction
        return sparse.csr_matrix(X), genes

    def close(self):
        pass

    def datasets(self) -> pd.DataFrame:
        return self._datasets.copy()

    def count(self, organism: str) -> pd.DataFrame:
        obs = self._obs[organism]
        sel = obs[(obs["is_primary_data"]) & (obs["tissue_type"] == "tissue")]
        return count_cells([sel[COUNTED]], COUNTED)

    def cells(self, organism, dataset_ids, tissue, cell_type) -> pd.DataFrame:
        obs = self._obs[organism]
        m = (obs["is_primary_data"] & (obs["tissue_type"] == "tissue") & obs["dataset_id"].isin(list(dataset_ids))
             & (obs["tissue"] == tissue) & (obs["cell_type"] == cell_type))
        return obs.loc[m, ["soma_joinid", *COUNTED]].reset_index(drop=True)

    def extract(self, organism, joinids) -> tuple:
        ids = np.sort(np.asarray(joinids, dtype=np.int64))
        obs = self._obs[organism].set_index("soma_joinid").loc[ids].reset_index()
        return self._X[organism][ids], obs[["soma_joinid", *COUNTED]], list(self._genes[organism])


# --------------------------------------------------------------------------- #
# the negative controls
# --------------------------------------------------------------------------- #
def select_negatives(census, out: Path, log) -> dict:
    datasets = census.datasets()
    parts = []
    for org in ORGANISMS:
        t0 = time.time()
        c = census.count(org)
        c.insert(0, "organism", org)
        parts.append(c)
        log(f"{org}: {int(c['cells'].sum()):,} cells counted ({BASE_FILTER}) in {time.time() - t0:.0f} s")
    counts = pd.concat(parts, ignore_index=True)
    cand, donors = negative_candidates(counts, datasets)
    names = datasets.groupby(datasets["collection_id"].astype(str)).agg(
        collection_name=("collection_name", "first"), collection_doi=("collection_doi", "first"))
    cand = cand.merge(names, left_on="collection_id", right_index=True, how="left")
    cand.to_csv(out / "candidates.tsv", sep="\t", index=False)
    log("")
    log(f"Qualifying candidates (>= {MIN_DONORS} donors with >= {MIN_CELLS} cells of one disease, sex and stage in "
        f"each of two technical conditions): {len(cand)}, of them excluded {int((cand['excluded'] != '').sum())} "
        f"(candidates.tsv). Excluded collections: "
        + ("; ".join(f"{c} ({w})" for c, w in sorted(set(zip(cand.loc[cand['excluded'] != '', 'collection_name'],
                                                                 cand.loc[cand['excluded'] != '', 'excluded'])))) or "none"))
    firsts = first_per_collection(cand)
    log(f"Collections with a qualifying candidate, not excluded: {len(firsts)}; tried in order, at most {MAX_TRIED}.")
    tried, chosen = [], []
    (out / "data").mkdir(exist_ok=True)
    (out / "joinids").mkdir(exist_ok=True)
    for i, ch in firsts.head(MAX_TRIED).iterrows():
        if len(chosen) == N_NEGATIVE:
            break
        org = ch["organism"]
        ds_ids = sorted(datasets.loc[datasets["collection_id"].astype(str) == ch["collection_id"], "dataset_id"].astype(str))
        q = donors[(donors["organism"] == org) & (donors["collection_id"] == ch["collection_id"])
                   & (donors["tissue"] == ch["tissue"]) & (donors["cell_type"] == ch["cell_type"])
                   & (donors["condition_a"] == ch["condition_a"]) & (donors["condition_b"] == ch["condition_b"])]
        meta = {r["donor_id"]: (r["disease_a"], r["sex_a"], r["development_stage_a"]) for _, r in q.iterrows()}
        cells = census.cells(org, ds_ids, ch["tissue"], ch["cell_type"])
        cells = cells.assign(condition=cells["assay"].astype(str) + " | " + cells["suspension_type"].astype(str),
                             donor_id=cells["donor_id"].astype(str))
        keep = (cells["donor_id"].isin(meta) & cells["condition"].isin([ch["condition_a"], ch["condition_b"]])).to_numpy()
        keep = keep & np.array([meta.get(d) == (str(a), str(b), str(c)) for d, a, b, c in
                                zip(cells["donor_id"], cells["disease"], cells["sex"], cells["development_stage"])],
                               dtype=bool)
        cells = cells[keep]
        ids = sample_units(cells, sorted(meta))
        X, obs, genes = census.extract(org, ids)
        row = dict(order=i + 1, rank=int(ch["rank"]), collection_id=ch["collection_id"],
                   collection_name=ch["collection_name"], organism=org, tissue=ch["tissue"], cell_type=ch["cell_type"],
                   condition_a=ch["condition_a"], condition_b=ch["condition_b"], donors=int(ch["donors"]),
                   cells=int(len(ids)))
        if not raw_counts(X):
            tried.append(dict(row, outcome="fails: the counts are not non-negative integers"))
            log(f"Tried {i + 1}: {ch['collection_name']!r}: the counts are not raw integers")
            continue
        number = len(chosen) + 1
        obs = obs.assign(condition=obs["assay"].astype(str) + " | " + obs["suspension_type"].astype(str))
        obs["group"] = np.where(obs["condition"] == ch["condition_a"], "A", "B")
        obs["unit"] = obs["donor_id"].astype(str) + " | " + obs["group"]
        path = out / "data" / f"N{number}.npz"
        P.save_npz(path, X, obs[["soma_joinid", "unit", "group", *COUNTED, "condition"]].astype(
            {c: str for c in ["unit", "group", *COUNTED, "condition"]}), genes)
        try:
            bg = P.load_background(f"N{number}", dict(path=str(path), donor="unit", counts="X"))
            plan = P.plan_background(bg)
        except ValueError as exc:
            path.unlink()
            tried.append(dict(row, outcome=f"fails: no pool of pairs ({exc})"))
            log(f"Tried {i + 1}: {ch['collection_name']!r}: no pool of pairs ({exc})")
            continue
        claims = negative_claims(plan, number)
        sexes, names = [], P.unique_names(genes)
        for unit in sorted(set(obs["unit"])):
            rows = np.where(obs["unit"].to_numpy() == unit)[0]
            fx, fy = sex_shares(X, names, rows, org)
            sexes.append(dict(unit=unit, cells=len(rows), label=str(obs["sex"].iloc[rows[0]]), xist=fx, y=fy,
                              call=sex_call(fx, fy)))
        by_donor = {}
        for s in sexes:
            by_donor.setdefault(s["unit"].rsplit(" | ", 1)[0], set()).add(s["call"])
        conflicts = sorted(d for d, cs in by_donor.items() if {"female", "male"} <= cs)
        (out / "joinids" / f"N{number}.txt").write_text("\n".join(map(str, ids)) + "\n")
        chosen.append(dict(id=f"N{number}", **row, release=census.release, datasets=ds_ids,
                           npz_sha256=P.sha256(path), genes=len(genes), joinids=f"joinids/N{number}.txt",
                           groups=dict(A=ch["condition_a"], B=ch["condition_b"]),
                           donor_meta={d: dict(disease=m[0], sex=m[1], development_stage=m[2]) for d, m in sorted(meta.items())},
                           units=sexes, sex_conflicts=conflicts,
                           pool=dict(pairs_per_level=plan["pairs_per_level"], fewer_pairs_because=plan["fewer_pairs_because"],
                                     levels=plan["levels"], positive_control=plan["positive_control"]),
                           claims=claims))
        tried.append(dict(row, outcome=f"chosen: N{number}"))
        log(f"Tried {i + 1}: {ch['collection_name']!r} ({ch['collection_id']}), {org}, {ch['tissue']!r}, "
            f"{ch['cell_type']!r}: A = {ch['condition_a']!r}, B = {ch['condition_b']!r}; {int(ch['donors'])} donors, "
            f"{len(ids)} cells; chosen as N{number}: {plan['pairs_per_level']} pairs per level, {len(claims)} claims; "
            f"sex calls in conflict between the conditions: {conflicts or 'none'}")
    pd.DataFrame(tried).to_csv(out / "tried.tsv", sep="\t", index=False)
    if len(chosen) < N_NEGATIVE:
        log(f"Only {len(chosen)} of {N_NEGATIVE} negative controls: nothing is put in their place (plan_s6.md).")
    return dict(candidates=int(len(cand)), excluded=int((cand["excluded"] != "").sum()), collections=int(len(firsts)),
                tried=len(tried), chosen=chosen)


# --------------------------------------------------------------------------- #
# the positive control: Kang et al. 2018
# --------------------------------------------------------------------------- #
class Fetch:
    """GET with up to 4 attempts on a network error or HTTP 5xx; every attempt recorded."""

    def __init__(self, log):
        self.log, self.rows = log, []

    def get(self, url: str) -> bytes:
        for attempt in range(4):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
                with urllib.request.urlopen(req, timeout=300) as r:
                    body = r.read()
                self.rows.append(dict(url=url, status=200, bytes=len(body), sha256=sha256_bytes(body), attempt=attempt + 1))
                return body
            except urllib.error.HTTPError as e:
                self.rows.append(dict(url=url, status=e.code, bytes=0, sha256="-", attempt=attempt + 1))
                if e.code < 500:
                    raise
            except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
                self.rows.append(dict(url=url, status=0, bytes=0, sha256="-", attempt=attempt + 1))
                self.log(f"  attempt {attempt + 1} failed ({type(e).__name__}): {url}")
            if attempt < 3:
                time.sleep(2 ** (attempt + 1))
        raise RuntimeError(f"four attempts failed: {url}")


def read_mtx(raw: bytes):
    """A gzipped MatrixMarket coordinate file as CSR (as stored): its first line that is not a comment gives
    the rows, the columns and the entries, each later line one entry (1-based row, column, value)."""
    from scipy import sparse
    df = pd.read_csv(io.BytesIO(gzip.decompress(raw)), sep=r"\s+", comment="%", header=None, engine="c",
                     dtype=np.float64)
    n_rows, n_cols, nnz = (int(x) for x in df.iloc[0, :3])
    data = df.iloc[1:, :3].to_numpy()
    if len(data) != nnz:
        raise ValueError(f"the matrix declares {nnz} entries and holds {len(data)}")
    return sparse.csr_matrix((data[:, 2], (data[:, 0].astype(np.int64) - 1, data[:, 1].astype(np.int64) - 1)),
                             shape=(n_rows, n_cols))


def read_lines(raw: bytes) -> list:
    return [ln for ln in gzip.decompress(raw).decode().splitlines() if ln.strip()]


def kang_cells(table_raw: bytes) -> pd.DataFrame:
    """The authors' table: a header of 7 names and rows of 8 fields, the first the cell's name."""
    lines = read_lines(table_raw)
    head = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[1:]]
    if any(len(r) != len(head) + 1 for r in rows):
        raise ValueError("the table's rows do not have one field more than its header")
    df = pd.DataFrame([r[1:] for r in rows], columns=head)
    df.insert(0, "name", [r[0] for r in rows])
    return df


def kang_lane_names(ctrl: list, stim: list) -> list:
    """The table's names of the stimulated lane's barcodes: a barcode also in the control lane carries '1'."""
    seen = set(ctrl)
    return [b + "1" if b in seen else b for b in stim]


def positive_kang(files: dict, out: Path, log) -> dict:
    """`files`: {name of KANG_FILES: bytes}. The structure check, then the claims (plan_s6.md section 2)."""
    from scipy import sparse
    genes_rows = [ln.split("\t") for ln in read_lines(files["genes.tsv.gz"])]
    symbols = P.unique_names([r[1] if len(r) > 1 else r[0] for r in genes_rows])
    table = kang_cells(files["cells.tsv.gz"]).set_index("name")
    lanes = dict(ctrl=read_lines(files["ctrl.barcodes.tsv.gz"]), stim=read_lines(files["stim.barcodes.tsv.gz"]))
    names = dict(ctrl=lanes["ctrl"], stim=kang_lane_names(lanes["ctrl"], lanes["stim"]))
    if not table.index.is_unique or any(len(set(v)) != len(v) for v in lanes.values()):
        raise ValueError("a cell's name is not unique in the authors' table or in a lane's barcodes")
    blocks, obs_parts, problems = [], [], []
    for lane in ("ctrl", "stim"):
        M = read_mtx(files[f"{lane}.mtx.gz"])
        if M.shape == (len(symbols), len(names[lane])):
            M = M.T.tocsr()
        elif M.shape != (len(names[lane]), len(symbols)):
            raise ValueError(f"{lane}: a matrix of {M.shape} for {len(symbols)} genes and {len(names[lane])} barcodes")
        missing = [n for n in names[lane] if n not in table.index]
        if missing:
            raise ValueError(f"{lane}: {len(missing)} barcodes not in the authors' table, e.g. {missing[:3]}")
        o = table.loc[names[lane]].reset_index()
        wrong = int((o["stim"] != lane).sum())
        if wrong:
            problems.append(f"{lane}: {wrong} of its barcodes are labelled otherwise in the table")
        blocks.append(M)
        obs_parts.append(o.assign(lane=lane))
    X = sparse.vstack(blocks, format="csr")
    obs = pd.concat(obs_parts, ignore_index=True)
    obs["ind"] = obs["ind"].astype(str)
    integers = raw_counts(X)
    singlet = (obs["multiplets"] == "singlet").to_numpy()
    rows = []
    for (d, cond), g in obs[singlet].groupby(["ind", "stim"], sort=True):
        idx = g.index.to_numpy()
        fx, fy = sex_shares(X, symbols, idx, "homo_sapiens")
        rows.append(dict(donor=d, condition=cond, singlets=len(idx),
                         monocytes=int((g["cell"] == KANG_TYPE).sum()), xist=fx, y=fy, call=sex_call(fx, fy)))
    structure = pd.DataFrame(rows)
    (out / "kang").mkdir(exist_ok=True)
    structure.to_csv(out / "kang" / "structure.tsv", sep="\t", index=False)
    calls = structure.groupby("donor")["call"].agg(lambda s: set(s))
    conflicts = sorted(d for d, cs in calls.items() if {"female", "male"} <= cs)
    both = sorted(d for d, g in structure.groupby("donor") if set(g["condition"]) >= {"ctrl", "stim"}
                  and (g["monocytes"] > 0).all())
    res = dict(source="Kang et al. 2018, GEO GSE96583, batch 2 (GSM2560248 control, GSM2560249 stimulated)",
               cells=int(len(obs)), singlets=int(singlet.sum()), integers=integers, label_problems=problems,
               donors_in_both=both, sex_conflicts=conflicts, type=KANG_TYPE)
    log(f"Kang: {len(obs)} cells ({int(singlet.sum())} singlets); counts integers: {integers}; "
        f"donors with {KANG_TYPE} in both conditions: {len(both)}; sex calls in conflict: {conflicts or 'none'}; "
        f"label problems: {problems or 'none'}")
    if not integers or problems or conflicts or len(both) < MIN_KANG_DONORS:
        res.update(passed=False, claims=[])
        log("Kang fails its structure check: no positive control (plan_s6.md).")
        return res
    mono = np.where(singlet & (obs["cell"] == KANG_TYPE).to_numpy() & obs["ind"].isin(both).to_numpy())[0]
    rng = np.random.default_rng(SEED_CELLS)
    keep = []
    for (d, cond), g in obs.iloc[mono].groupby(["ind", "stim"], sort=True):
        idx = np.sort(g.index.to_numpy())
        if len(idx) > MAX_CELLS:
            idx = np.sort(rng.choice(idx, size=MAX_CELLS, replace=False))
        keep.append(idx)
    keep = np.sort(np.concatenate(keep))
    Xm, om = X[keep], obs.iloc[keep].reset_index(drop=True)
    tot = np.asarray(Xm.sum(axis=1)).ravel().astype(float)
    tot = np.where(tot == 0, 1.0, tot)
    col = {g: j for j, g in enumerate(symbols)}
    genes = []
    for g in ISGS:
        if g not in col:
            genes.append(dict(gene=g, present=False, eligible=False))
            continue
        v = np.asarray(Xm[:, col[g]].todense()).ravel()
        det, induced = {}, []
        for cond in ("ctrl", "stim"):
            det[cond] = float(np.mean([(v[(om["ind"] == d).to_numpy() & (om["stim"] == cond).to_numpy()] > 0).mean()
                                       for d in both]))
        for d in both:
            lv = np.log1p(v / tot * 1e4)
            m = {c: lv[(om["ind"] == d).to_numpy() & (om["stim"] == c).to_numpy()].mean() for c in ("ctrl", "stim")}
            induced.append(bool(m["stim"] > m["ctrl"]))
        ok = det["ctrl"] >= DOMAIN_DETECTION and det["stim"] >= DOMAIN_DETECTION and all(induced)
        genes.append(dict(gene=g, present=True, detection_ctrl=det["ctrl"], detection_stim=det["stim"],
                          induced_in=int(sum(induced)), donors=len(both), eligible=ok))
    eligible = [g["gene"] for g in genes if g["eligible"]]
    perm = np.random.default_rng(SEED_MATCHING).permutation(len(eligible))
    order = [eligible[k] for k in perm]
    claims = [dict(id=f"K1p{k // 2 + 1}", kind="pair", genes=[order[k], order[k + 1]], higher_in="stim")
              for k in range(0, len(order) - 1, 2)]
    if len(eligible) >= 4:
        claims.append(dict(id="K1s", kind="set", genes=list(eligible), higher_in="stim"))
    (out / "data").mkdir(exist_ok=True)
    om_save = om[["name", "ind", "stim", "cell", "multiplets", "lane"]].astype(str)
    P.save_npz(out / "data" / "K1.npz", Xm, om_save, symbols)
    res.update(passed=True, genes=genes, eligible=eligible, unmatched=[order[-1]] if len(order) % 2 else [],
               cells_kept=int(len(keep)), npz_sha256=P.sha256(out / "data" / "K1.npz"), claims=claims,
               groups=dict(A="stim", B="ctrl"))
    log(f"Kang: {len(keep)} {KANG_TYPE} kept (at most {MAX_CELLS} per donor and condition); eligible genes "
        f"{len(eligible)} of {len(ISGS)}: {eligible}; claims: " + "; ".join(f"{c['id']} {'-'.join(c['genes'])}" if c['kind'] == 'pair'
                                                                              else f"{c['id']} set of {len(c['genes'])}"
                                                                              for c in claims))
    return res


def synthetic_kang_files(seed: int = 5) -> dict:
    """GEO's files for the dry run, fabricated in their real formats (SYNTHETIC): 8 donors, both conditions,
    monocytes and T cells; the ISGs induced with a shared per-cell factor in the stimulated cells; IFIT2 below
    the detection rule in the control cells; one doublet per lane; a stimulated barcode that is also a control's."""
    rng = np.random.default_rng(seed)
    filler = [f"GENE{i:03d}" for i in range(60)]
    symbols = list(ISGS) + ["XIST", "RPS4Y1", "DDX3Y"] + filler
    donors = ["101", "1015", "1016", "1039", "107", "1244", "1256", "1488"]
    files, table = {}, ["tsne1\ttsne2\tind\tstim\tcluster\tcell\tmultiplets"]
    ctrl_barcodes = []
    for lane in ("ctrl", "stim"):
        cols, barcodes = [], []
        for k, d in enumerate(donors):
            for t, (ct, n) in enumerate(((KANG_TYPE, 50), ("CD4 T cells", 20))):
                for i in range(n):
                    bc = f"{lane.upper()[:1]}{k:02d}T{t}{i:04d}-1"
                    if lane == "stim" and k == 0 and ct == KANG_TYPE and i == 0:
                        bc = ctrl_barcodes[0]  # the same barcode in both lanes
                    lib = np.exp(rng.normal(0, 0.3))
                    z = rng.normal(0, 1)
                    lam = np.full(len(symbols), 0.0)
                    base = 0.6 if lane == "ctrl" else 6.0
                    lam[:len(ISGS)] = lib * base * np.exp(0.6 * z * (lane == "stim"))
                    if lane == "ctrl":
                        lam[ISGS.index("IFIT2")] = 0.02
                    lam[len(ISGS)] = lib * 3.0 * (k % 2 == 0)
                    lam[len(ISGS) + 1:len(ISGS) + 3] = lib * 2.0 * (k % 2 == 1)
                    lam[len(ISGS) + 3:] = lib * rng.uniform(0.5, 3.0, len(filler))
                    cols.append(rng.poisson(lam))
                    barcodes.append(bc)
                    name = bc + "1" if lane == "stim" and bc in set(ctrl_barcodes) else bc
                    mult = "doublet" if i == 1 and k == 0 else "singlet"
                    table.append(f"{name}\t0.1\t0.2\t{d}\t{lane}\t1\t{ct}\t{mult}")
        if lane == "ctrl":
            ctrl_barcodes = list(barcodes)
        M = np.array(cols).T  # genes x cells, as 10x writes it
        r, c = np.nonzero(M)
        mtx = "%%MatrixMarket matrix coordinate integer general\n%\n" + f"{M.shape[0]} {M.shape[1]} {len(r)}\n"
        mtx += "".join(f"{i + 1} {j + 1} {M[i, j]}\n" for i, j in zip(r, c))
        files[f"{lane}.mtx.gz"] = gzip.compress(mtx.encode(), mtime=0)
        files[f"{lane}.barcodes.tsv.gz"] = gzip.compress(("\n".join(barcodes) + "\n").encode(), mtime=0)
    files["genes.tsv.gz"] = gzip.compress("".join(f"ENSG{i:011d}\t{s}\n" for i, s in enumerate(symbols)).encode(), mtime=0)
    files["cells.tsv.gz"] = gzip.compress(("\n".join(table) + "\n").encode(), mtime=0)
    files["filelist.txt"] = b"#Archive/File\tName\tTime\tSize\tType\n(SYNTHETIC)\n"
    return files


# --------------------------------------------------------------------------- #
# the run
# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--synthetic", action="store_true", help="the dry run: a fabricated Census and GEO files")
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    tag = "SYNTHETIC: " if a.synthetic else ""
    log(f"# {tag}Level B, S6: the controls (validation/level_b/controls.py; plan_s6.md)")
    log(f"python {platform.python_version()}; numpy {np.__version__}; pandas {pd.__version__}; "
        f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    census = SyntheticCensus() if a.synthetic else Census()
    log(f"Census release: {census.release}")
    log("")
    log("## Negative controls")
    try:
        neg = select_negatives(census, a.out, log)
    finally:
        census.close()
    log("")
    log("## The positive control")
    fetch = Fetch(log)
    if a.synthetic:
        files = synthetic_kang_files()
        for name, body in files.items():
            fetch.rows.append(dict(url=f"SYNTHETIC:{name}", status=200, bytes=len(body), sha256=sha256_bytes(body), attempt=1))
    else:
        files = {name: fetch.get(url) for name, (url, _) in KANG_FILES.items()}
    (a.out / "kang").mkdir(exist_ok=True)
    pd.DataFrame(fetch.rows).to_csv(a.out / "kang" / "requests.tsv", sep="\t", index=False)
    log("Files: " + "; ".join(f"{n} {len(b):,} bytes sha256 {sha256_bytes(b)}" for n, b in files.items()))
    pos = positive_kang(files, a.out, log)
    claims = dict(synthetic=a.synthetic, release=census.release,
                  negative=[dict({k: v for k, v in c.items() if k not in ("units", "donor_meta")}) for c in neg["chosen"]],
                  positive=[pos] if pos.get("passed") else [], positive_failed=[] if pos.get("passed") else [pos])
    (a.out / "claims.json").write_text(json.dumps(claims, indent=1, default=str) + "\n")
    sel = dict(synthetic=a.synthetic, release=census.release, plan="validation/level_b/plan_s6.md",
               rules=dict(min_donors=MIN_DONORS, min_cells=MIN_CELLS, max_donors=MAX_DONORS, max_cells=MAX_CELLS,
                          max_tried=MAX_TRIED, n_negative=N_NEGATIVE, seeds=[SEED_CELLS, SEED_DIRECTIONS, SEED_MATCHING],
                          base_filter=BASE_FILTER, excluded_datasets=EXCLUDED_DATASETS, excluded_names=EXCLUDED_NAMES,
                          excluded_dois=EXCLUDED_DOIS, isgs=list(ISGS), domain_detection=DOMAIN_DETECTION),
               negative=dict({k: v for k, v in neg.items() if k != "chosen"},
                             chosen=[dict(c, claims=len(c["claims"])) for c in neg["chosen"]]),
               positive={k: v for k, v in pos.items() if k != "claims"})
    (a.out / "selection.json").write_text(json.dumps(sel, indent=1, default=str) + "\n")
    log("")
    log(f"Negative controls: {len(neg['chosen'])}; positive: {1 if pos.get('passed') else 0}. claims.json, "
        "selection.json, candidates.tsv, tried.tsv, joinids/, kang/; the extractions in data/.")
    (a.out / "controls.log").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
