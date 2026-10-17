#!/usr/bin/env python3
"""The extraction of the pilot on real data (validation/exploratory/real_pilot_plan.md).

Exploratory, for development data only. It reads the CELLxGENE Census, release 2025-11-08 (v1's),
in GitHub Actions (.github/workflows/real-pilot.yml), applies each case's rules before any metric
is computed, and writes each case's cells for real_pilot.py. It does not import the engine: the
rules count cells and detections, case 2's structure check reads the sex genes, and case 3's
claims are chosen on the pooled cells of both conditions, before any analysis.

    python validation/exploratory/real_pilot_select.py --out-dir DATA               # the Census
    python validation/exploratory/real_pilot_select.py --out-dir DATA --synthetic   # offline

In DATA it writes case1.npz, case2.npz and case3.npz for the cases that run (counts dense int32
over 2,000 genes, the gene symbols, the obs; panel.save_npz) and select.json (every choice, the
facts that decided it, the claims); the log goes to stdout. With --synthetic the same code runs on
a fabricated Census (SyntheticCensus): every output says SYNTHETIC and says nothing about real data.

Loading, as deviation D1 loads a background: one cell type; at most 500 cells per donor and
condition (where there are more, a sample of the sorted joinids drawn with the seed 20261015); the
case's named genes, then the most expressed by mean count over the extracted cells (ties by column
order), 2,000 in all, in column order, among the genes the Census marks as measured in every
dataset the case reads. Donor ids that name two animals ("/", as 3_38_F/3_39_F) and missing ids
are excluded. Only the donors that qualify under the case's rule are extracted (case 3: at most 24
of them, a sample of the sorted ids drawn with the seed where more qualify, for the runners' time
and memory).

An exception in a case is printed with its traceback, and the script exits with 1 after the other
cases, so that the analyses do not start; a case that does not qualify by its rule is recorded and
does not run.

The second run (--second; the plan's section "Deviation: a second run of cases 1 and 2", decided by
the owner on 2026-10-16 and written before it) extracts cases 1 and 2 only, with two changes. Case 1
is preprint section 4's setup as the flagship audit reproduced it (validation/flagship_audit,
section E): every primary Smart-seq2 TMS FACS cell at the 3-month stage and the "20-month-old stage
and over" stage, both sexes, every cell type and tissue, no cap per mouse, and the donor ids as the
Census gives them for units, composite ids included; with --audit-results, its units and their
cells are checked against the audit's per-mouse table. Case 2 calls a mouse's sex by the contrast of
the sex genes (sex_call, "contrast"), not by the plan's thresholds.
"""
from __future__ import annotations

import argparse
import itertools
import json
import platform
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "prereg"))
import panel as P  # noqa: E402  (unique_names, _col_stats, save_npz, sha256: numpy and pandas, no engine)
from select_backgrounds import _known, count_cells  # noqa: E402  (v1's selection: placeholder ids, chunked counts)

RELEASE = "2025-11-08"                # v1's Census release (the audit's "stable" build)
SEED = 20261015                       # the samples of cells
MAX_CELLS = 500                       # per donor and condition
N_GENES = P.N_GENES                   # 2,000, as D1
MIN_MOUSE_CELLS, MIN_GROUP_MICE = 20, 3          # cases 1 and 2
NEG_MIN_DONORS, NEG_MIN_CELLS = 4, 50            # case 3
NEG_MAX_DONORS = 24                   # case 3: at most 24 donors (B1's minimum in v1), a sample where more qualify
MAX_RANKS = 5                         # case 3: candidates tried, in order, for raw counts (as v1)
PAIR = ("Smad3", "Col1a1")
POS_PAIR, NEG_PAIR = ("Actb", "Gapdh"), ("Tst", "Lrrc42")
Y_GENES = ("Ddx3y", "Eif2s3y", "Kdm5d", "Uty")
SEX_GENES = ("Xist",) + Y_GENES
TMS_GENES = PAIR + POS_PAIR + NEG_PAIR + SEX_GENES
YOUNG, OLD = "3-month-old stage", "20-month-old stage and over"
CALL_MIN, CALL_MAX = 0.50, 0.05       # case 2: a female has Xist in >= 50% of cells and a Y gene in <= 5%
CALL_SHARE, CALL_RATIO = 0.05, 10.0   # case 2, second run: the other sex's genes in <= 5% of cells, the own in >= 10x that share and >= 5%
CHUNK_CELLS = 10000                   # case 1, second run: the cells read per call to the Census
AUDIT_AGE = {"3m": "young", "20m": "old"}  # the audit's age tokens (flagship_audit/audit_tms.py)
# Case 3 leaves out, by collection name, title or DOI: Tabula Muris Senis and Tabula Muris (their
# mice by FACS and by droplet are a candidate truth for v2), and the confirmatory set's Kang et al.
# 2018 and Levy et al. 2020.
EXCLUDED_NAMES = ("tabula muris",
                  "multiplexed droplet single-cell rna-sequencing using natural genetic variation",
                  "age-related loss of gene-to-gene transcriptional coordination")
EXCLUDED_DOIS = ("10.1038/nbt.4042", "10.1038/s42255-020-00304-4")
ORGANISMS = {"mus_musculus": "Mus musculus", "homo_sapiens": "Homo sapiens"}
CASE3_COLS = ["dataset_id", "donor_id", "assay", "suspension_type", "tissue", "cell_type", "sex",
              "development_stage"]
OBS_COLS = ["dataset_id", "donor_id", "assay", "suspension_type", "tissue", "cell_type", "sex",
            "development_stage"]
CASE3_PAIRS, CASE3_TOP, CASE3_MID, CASE3_MEANS = 8, (0, 100), (500, 600), 4
SESOI_PAIR, SESOI_MEAN = 0.1, 0.5
COUPLING_STRENGTH, MODULE_FOLD, MODULE_FRAC = 1.0, 2.0, 1.0


# --------------------------------------------------------------------------- #
# filters and small helpers
# --------------------------------------------------------------------------- #
def soma_filter(flt: dict) -> str:
    """A SOMA value filter from {column: value, list of values, or bool}."""
    parts = []
    for col, val in flt.items():
        if isinstance(val, bool):
            parts.append(f"{col} == {val}")
        elif isinstance(val, (list, tuple)):
            parts.append(f"{col} in [" + ", ".join(repr(str(v)) for v in val) + "]")
        else:
            parts.append(f"{col} == {str(val)!r}")
    return " and ".join(parts)


def apply_filter(df: pd.DataFrame, flt: dict) -> np.ndarray:
    """The rows of `df` that `soma_filter(flt)` selects (the synthetic Census)."""
    keep = np.ones(len(df), bool)
    for col, val in flt.items():
        v = df[col]
        if isinstance(val, bool):
            keep &= v.to_numpy() == val
        elif isinstance(val, (list, tuple)):
            keep &= v.astype(str).isin([str(x) for x in val]).to_numpy()
        else:
            keep &= (v.astype(str) == str(val)).to_numpy()
    return keep


def units(donor: pd.Series) -> np.ndarray:
    """True where a donor id names one donor: present (not a placeholder) and not composite."""
    composite = donor.astype(str).str.contains("/", regex=False).to_numpy()
    return _known(donor) & ~composite


def sample_ids(ids_by_unit: dict, max_cells: int, seed: int = SEED) -> np.ndarray:
    """At most `max_cells` joinids of each unit (donor, or donor and condition), a sample of its
    sorted joinids drawn with `seed` where there are more; the units in sorted order."""
    rng = np.random.default_rng(seed)
    out = []
    for u in sorted(ids_by_unit):
        ids = np.sort(np.asarray(ids_by_unit[u], dtype=np.int64))
        if len(ids) > max_cells:
            ids = np.sort(rng.choice(ids, size=max_cells, replace=False))
        out.append(ids)
    return np.sort(np.concatenate(out)) if out else np.zeros(0, np.int64)


def raw_counts(X) -> bool:
    v = X.data if hasattr(X, "data") and not isinstance(X, np.ndarray) else np.asarray(X)
    return bool(v.size == 0 or (np.all(v >= 0) and np.all(v == np.round(v))))


def cut_genes(X, names: list, named=()) -> tuple[np.ndarray, list, list]:
    """D1's cut (anchors._load): the named genes, then the most expressed by mean count over the
    cells (ties by column order), N_GENES in all, in column order; dense int32. Returns the matrix,
    the kept symbols (unique) and the named genes not found."""
    names = P.unique_names(names)
    want = {str(g).upper() for g in named}
    first = sorted(j for j, g in enumerate(names) if g.upper() in want)
    found = {names[j].upper() for j in first}
    missing = [g for g in named if g.upper() not in found]
    mean, _ = P._col_stats(X)
    order = np.argsort(-mean, kind="mergesort")
    taken = set(first)
    keep = list(dict.fromkeys([*first, *[int(j) for j in order if j not in taken]]))
    keep = sorted(keep[:max(N_GENES, len(first))])
    sub = X[:, keep]
    dense = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
    return dense.astype(np.int32), [names[j] for j in keep], missing


def md_table(df: pd.DataFrame, cols=None, max_rows: int = 15) -> list[str]:
    cols = list(cols or df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, row in df.head(max_rows).iterrows():
        out.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
    if len(df) > max_rows:
        out.append(f"| ... {len(df) - max_rows} more rows | " + " | ".join("" for _ in cols[1:]) + " |")
    return out


# --------------------------------------------------------------------------- #
# the Census, and a fabricated one for the dry run
# --------------------------------------------------------------------------- #
class Census:
    """The CELLxGENE Census at RELEASE (cellxgene_census; requirements-select.txt)."""

    synthetic = False

    def __init__(self):
        import cellxgene_census
        self._cc = cellxgene_census
        self._census = cellxgene_census.open_soma(census_version=RELEASE)

    def close(self):
        self._census.close()

    def summary(self) -> dict:
        s = self._census["census_info"]["summary"].read().concat().to_pandas()
        return dict(zip(s["label"].astype(str), s["value"].astype(str)))

    def datasets(self) -> pd.DataFrame:
        return self._census["census_info"]["datasets"].read().concat().to_pandas()

    def _obs(self, organism: str):
        return self._census["census_data"][organism].obs

    def obs_columns(self, organism: str) -> list:
        return list(self._obs(organism).schema.names)

    def obs_batches(self, organism: str, flt: dict, cols):
        return self._obs(organism).read(value_filter=soma_filter(flt), column_names=list(cols))

    def obs_table(self, organism: str, flt: dict, cols) -> pd.DataFrame:
        return self.obs_batches(organism, flt, ["soma_joinid", *cols]).concat().to_pandas()

    def var(self, organism: str) -> pd.DataFrame:
        return (self._census["census_data"][organism].ms["RNA"].var
                .read(column_names=["soma_joinid", "feature_id", "feature_name"]).concat().to_pandas())

    def presence(self, organism: str):
        """Datasets (rows, by the datasets table's soma_joinid) x genes (columns, by var soma_joinid)."""
        return self._cc.get_presence_matrix(self._census, ORGANISMS[organism], "RNA").tocsr()

    def counts(self, organism: str, *, flt=None, obs_joinids=None, var_joinids=None, var_names=None,
               obs_cols=()):
        """Raw counts (CSR, cells by joinid, genes by joinid), obs and var of the cells that `flt`
        or `obs_joinids` select, over the genes `var_joinids` or `var_names` select (else all)."""
        from scipy import sparse
        kw = dict(X_name="raw", obs_column_names=["soma_joinid", *obs_cols],
                  var_column_names=["soma_joinid", "feature_id", "feature_name"])
        if flt is not None:
            kw["obs_value_filter"] = soma_filter(flt)
        if obs_joinids is not None:
            kw["obs_coords"] = np.asarray(obs_joinids, dtype=np.int64)
        if var_joinids is not None:
            kw["var_coords"] = np.asarray(var_joinids, dtype=np.int64)
        if var_names is not None:
            kw["var_value_filter"] = soma_filter(dict(feature_name=list(var_names)))
        ad = self._cc.get_anndata(self._census, ORGANISMS[organism], **kw)
        o = np.argsort(ad.obs["soma_joinid"].to_numpy(), kind="mergesort")
        v = np.argsort(ad.var["soma_joinid"].to_numpy(), kind="mergesort")
        X = sparse.csr_matrix(ad.X)[o][:, v]
        return (X, ad.obs.iloc[o].reset_index(drop=True), ad.var.iloc[v].reset_index(drop=True))


class SyntheticCensus:
    """A fabricated Census with the same interface, for the dry run: TMS-shaped FACS tissues with
    the audit's design of mice (composite ids, the 21-month females of the mammary gland only, one
    mouse whose sex genes contradict its label), a collection of cells and nuclei from the same
    donors, and collections that the rules exclude (by name, by DOI, by disease). Nothing is real."""

    synthetic = True
    TMS_TITLE = "{} - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - {}"

    def __init__(self, seed: int = 7):
        self.rng = np.random.default_rng(seed)
        self._datasets = []
        self._obs, self._X, self._var = {}, {}, {}
        self._build_mouse()
        self._build_human()
        self._datasets = pd.DataFrame(self._datasets)
        self._datasets["soma_joinid"] = np.arange(len(self._datasets))

    def _dataset(self, dataset_id, collection_id, collection_name, doi, title):
        self._datasets.append(dict(dataset_id=dataset_id, collection_id=collection_id,
                                   collection_name=collection_name, collection_doi=doi,
                                   dataset_title=title, citation=f"Publication: https://doi.org/{doi}"))

    @staticmethod
    def _frame(rows):
        df = pd.DataFrame(rows)
        df.insert(0, "soma_joinid", np.arange(len(df), dtype=np.int64))
        return df

    def _build_mouse(self):
        rng = self.rng
        tissues = {"Lung": ("lung", ["fibroblast of lung", "endothelial cell", "B cell"]),
                   "Mammary gland": ("mammary gland", ["stromal cell", "endothelial cell", "basal cell"]),
                   "Limb muscle": ("limb muscle", ["skeletal muscle satellite cell", "endothelial cell",
                                                   "mesenchymal stem cell"])}
        mice = [("3_8_M", "male", YOUNG), ("3_9_M", "male", YOUNG), ("3_10_M", "male", YOUNG),
                ("3_11_M", "male", YOUNG), ("3_10_M/3_11_M", "male", YOUNG),
                ("24_58_M", "male", OLD), ("24_59_M", "male", OLD), ("24_60_M", "male", OLD),
                ("24_61_M", "male", OLD), ("3_38_F", "female", YOUNG), ("3_39_F", "female", YOUNG),
                ("3_56_F", "female", YOUNG), ("3_57_F", "female", YOUNG), ("3_58_F", "female", YOUNG),
                ("3_38_F/3_39_F", "female", YOUNG), ("21_53_F", "female", OLD),
                ("21_54_F", "female", OLD), ("21_55_F", "female", OLD),
                ("18_45_M", "male", "18-month-old stage"), ("18_47_F", "female", "18-month-old stage")]
        rows = []
        for t, (tissue, types) in tissues.items():
            ds = f"tms-facs-{tissue.replace(' ', '-')}"
            self._dataset(ds, "tms", "Tabula Muris Senis", "10.1038/s41586-020-2496-1",
                          self.TMS_TITLE.format(t, "Smart-seq2"))
            for donor, sex, stage in mice:
                if sex == "male" and tissue == "mammary gland":
                    continue
                if stage == OLD and sex == "female" and tissue != "mammary gland":
                    continue  # as in TMS: the 21-month females are of the mammary gland only
                for ct in types:
                    n = 8 if "/" in donor else int(rng.integers(14, 70))
                    rows += [dict(dataset_id=ds, assay="Smart-seq2", suspension_type="cell", tissue=tissue,
                                  tissue_type="tissue", cell_type=ct, sex=sex, development_stage=stage,
                                  donor_id=donor, disease="normal", is_primary_data=True)] * n
        # an aggregate copy of some lung cells (not primary) and the same mice by droplet
        self._dataset("tms-facs-all", "tms", "Tabula Muris Senis", "10.1038/s41586-020-2496-1",
                      self.TMS_TITLE.format("All", "Smart-seq2"))
        self._dataset("tms-droplet-lung", "tms", "Tabula Muris Senis", "10.1038/s41586-020-2496-1",
                      self.TMS_TITLE.format("Lung", "10x 3' v2"))
        lung = [r for r in rows if r["tissue"] == "lung"][:200]
        rows += [dict(r, dataset_id="tms-facs-all", is_primary_data=False) for r in lung]
        rows += [dict(r, dataset_id="tms-droplet-lung", assay="10x 3' v2") for r in
                 [r for r in rows if r["dataset_id"] == "tms-facs-lung"]]
        obs = self._frame(rows)
        genes = list(TMS_GENES) + ["Dup1", "Dup1"] + [f"Gm{i:04d}" for i in range(787)]
        X = self._counts(obs, genes, mouse=True)
        self._obs["mus_musculus"], self._X["mus_musculus"] = obs, X
        self._var["mus_musculus"] = pd.DataFrame(dict(soma_joinid=np.arange(len(genes)),
                                                      feature_id=[f"ENSMUSG{i:011d}" for i in range(len(genes))],
                                                      feature_name=genes))

    def _build_human(self):
        rng = self.rng
        rows = []

        def add(ds, assay, susp, tissue, ct, donor, sex, stage, n, disease="normal"):
            rows.extend([dict(dataset_id=ds, assay=assay, suspension_type=susp, tissue=tissue, tissue_type="tissue",
                              cell_type=ct, sex=sex, development_stage=stage, donor_id=donor, disease=disease,
                              is_primary_data=True)] * n)
        self._dataset("heart-cells", "col-heart", "Synthetic cells and nuclei of the adult heart",
                      "10.0000/synthetic.heart", "Heart: cells")
        self._dataset("heart-nuclei", "col-heart", "Synthetic cells and nuclei of the adult heart",
                      "10.0000/synthetic.heart", "Heart: nuclei")
        for k in range(1, 8):
            sex, stage = ("female" if k % 2 else "male"), f"{40 + k}-year-old stage"
            add("heart-cells", "10x 3' v2", "cell", "heart left ventricle", "fibroblast", f"D{k}", sex, stage,
                int(rng.integers(55, 90)))
            add("heart-cells", "10x 3' v2", "cell", "heart left ventricle", "cardiac muscle cell", f"D{k}", sex,
                stage, int(rng.integers(20, 40)))
            if k < 7:  # D7: cells only
                add("heart-nuclei", "10x 3' v3", "nucleus", "heart left ventricle", "fibroblast", f"D{k}", sex,
                    stage, int(rng.integers(55, 90)))
                add("heart-nuclei", "10x 3' v3", "nucleus", "heart left ventricle", "cardiac muscle cell", f"D{k}",
                    sex, stage, int(rng.integers(60, 90)))
        self._dataset("kang", "col-kang", "Multiplexed droplet single-cell RNA-sequencing using natural genetic "
                      "variation", "10.1038/nbt.4042", "PBMC")
        for k in range(1, 9):
            for assay in ("10x 3' v1", "10x 3' v2"):
                add("kang", assay, "cell", "blood", "CD14-positive monocyte", f"K{k}", "male", "adult stage",
                    int(rng.integers(60, 90)))
        self._dataset("covid", "col-covid", "Synthetic COVID-19 blood", "10.0000/synthetic.covid", "Blood")
        for k in range(1, 11):
            for assay in ("10x 3' v3", "10x 5' v2"):
                add("covid", assay, "cell", "blood", "T cell", f"C{k}", "female", "adult stage", 70,
                    disease="COVID-19")
        self._dataset("lung-3p", "col-lung", "Synthetic lung by two chemistries", "10.0000/synthetic.lung", "3'")
        self._dataset("lung-5p", "col-lung", "Synthetic lung by two chemistries", "10.0000/synthetic.lung", "5'")
        for k in range(1, 6):
            add("lung-3p", "10x 3' v3", "cell", "lung", "alveolar macrophage", f"L{k}", "male", "adult stage", 60)
            add("lung-5p", "10x 5' v2", "cell", "lung", "alveolar macrophage", f"L{k}", "male", "adult stage", 60)
        obs = self._frame(rows)
        genes = [f"GENE{i:04d}" for i in range(797)] + ["MT-CO1", "MT-ND1", "MT-CYB"]
        X = self._counts(obs, genes, mouse=False)
        absent = genes.index("GENE0005")  # not measured in the nuclei's dataset: the presence rule drops it
        X[(obs["dataset_id"] == "heart-nuclei").to_numpy(), absent] = 0
        self._obs["homo_sapiens"], self._X["homo_sapiens"] = obs, X
        self._var["homo_sapiens"] = pd.DataFrame(dict(soma_joinid=np.arange(len(genes)),
                                                      feature_id=[f"ENSG{i:011d}" for i in range(len(genes))],
                                                      feature_name=genes))

    def _counts(self, obs: pd.DataFrame, genes: list, mouse: bool) -> np.ndarray:
        rng = self.rng
        n, g = len(obs), len(genes)
        base = np.exp(rng.normal(-1.2, 1.6, g))
        lib = np.exp(rng.normal(0.0, 0.35, n))
        if mouse:
            old_male = ((obs["development_stage"] == OLD) & (obs["sex"] == "male")).to_numpy()
            lib[old_male] *= 0.6  # the depth difference of TMS's old males
        else:
            nuc = (obs["suspension_type"] == "nucleus").to_numpy()
            lib[nuc] *= 0.5
        lam = np.outer(lib, base)
        for ct in obs["cell_type"].unique():  # a module per cell type
            m = (obs["cell_type"] == ct).to_numpy()
            mod = rng.choice(g, size=25, replace=False)
            lam[np.ix_(m, mod)] *= 3.0
        if not mouse:
            lam[np.ix_((obs["suspension_type"] == "nucleus").to_numpy(), np.arange(30))] *= 0.3
            lam[:, g - 3:] *= 40.0  # mitochondrial genes, high
        else:
            ix = {s: genes.index(s) for s in TMS_GENES}
            z = rng.normal(0, 1, n)
            fib = obs["cell_type"].isin(["fibroblast of lung", "stromal cell", "mesenchymal stem cell"]).to_numpy()
            for s in PAIR:
                lam[:, ix[s]] = lib * (2.0 * fib + 0.3) * np.exp(0.8 * z * fib)
            hk = rng.normal(0, 1, n)
            lam[:, ix["Actb"]], lam[:, ix["Gapdh"]] = lib * 60 * np.exp(0.5 * hk), lib * 40 * np.exp(0.5 * hk)
            lam[:, ix["Tst"]], lam[:, ix["Lrrc42"]] = lib * 1.0, lib * 1.0
            female = (obs["sex"] == "female").to_numpy()
            lam[:, ix["Xist"]] = lib * 15.0 * female
            low = (obs["donor_id"] == "3_38_F").to_numpy()  # Xist in about 40% of its cells, as the real 3_38_F
            lam[low, ix["Xist"]] = lib[low] * 0.5
            for s in Y_GENES:
                lam[:, ix[s]] = lib * 3.0 * ~female
            odd = (obs["donor_id"] == "3_57_F").to_numpy()  # labelled female; half its cells carry Y genes
            half = odd & (rng.random(n) < 0.5)
            for s in Y_GENES:
                lam[half, ix[s]] = lib[half] * 3.0
        return rng.poisson(lam).astype(np.int32)

    def close(self):
        pass

    def summary(self) -> dict:
        return dict(census_build_date="SYNTHETIC", census_schema_version="SYNTHETIC")

    def datasets(self) -> pd.DataFrame:
        return self._datasets.copy()

    def obs_columns(self, organism):
        return list(self._obs[organism].columns)

    def obs_batches(self, organism, flt, cols, chunk: int = 1000):
        obs = self._obs[organism]
        sel = obs[apply_filter(obs, flt)]
        for i in range(0, len(sel), chunk):  # categorical columns, as the Census's Arrow tables give
            yield sel.iloc[i:i + chunk][list(cols)].astype("category")

    def obs_table(self, organism, flt, cols) -> pd.DataFrame:
        obs = self._obs[organism]
        return obs[apply_filter(obs, flt)][["soma_joinid", *cols]].reset_index(drop=True)

    def var(self, organism):
        return self._var[organism].copy()

    def presence(self, organism):
        from scipy import sparse
        var, ds = self._var[organism], self._datasets
        M = np.ones((len(ds), len(var)), dtype=np.int8)
        if organism == "homo_sapiens":
            row = int(ds.loc[ds["dataset_id"] == "heart-nuclei", "soma_joinid"].iloc[0])
            M[row, list(var["feature_name"]).index("GENE0005")] = 0
        return sparse.csr_matrix(M)

    def counts(self, organism, *, flt=None, obs_joinids=None, var_joinids=None, var_names=None, obs_cols=()):
        from scipy import sparse
        obs, var, X = self._obs[organism], self._var[organism], self._X[organism]
        rows = np.ones(len(obs), bool) if flt is None else apply_filter(obs, flt)
        if obs_joinids is not None:
            rows &= obs["soma_joinid"].isin(np.asarray(obs_joinids)).to_numpy()
        cols = np.ones(len(var), bool)
        if var_joinids is not None:
            cols &= var["soma_joinid"].isin(np.asarray(var_joinids)).to_numpy()
        if var_names is not None:
            cols &= var["feature_name"].isin(list(var_names)).to_numpy()
        r, c = np.where(rows)[0], np.where(cols)[0]
        Xs = sparse.csr_matrix(X[np.ix_(r, c)].astype(np.float32))  # the Census's raw X is float32
        return (Xs, obs.iloc[r][["soma_joinid", *obs_cols]].reset_index(drop=True),
                var.iloc[c].reset_index(drop=True))


# --------------------------------------------------------------------------- #
# Tabula Muris Senis (FACS): cases 1 and 2
# --------------------------------------------------------------------------- #
def tms_dataset_ids(datasets: pd.DataFrame) -> list:
    """The TMS FACS datasets as the audit resolves them (examples/mi_coupling_tms/download_data.py)."""
    is_tms = datasets["collection_name"].str.contains("Tabula Muris Senis", case=False, na=False)
    is_facs = datasets["dataset_title"].str.contains("facs|smart-seq", case=False, na=False)
    hit = datasets[is_tms & is_facs]
    if hit.empty:
        hit = datasets[is_tms]
    return hit["dataset_id"].tolist()


def tms_flt(ids) -> dict:
    return dict(dataset_id=list(ids), assay="Smart-seq2", is_primary_data=True)


def read_tms(census, ids) -> pd.DataFrame:
    """Every TMS FACS primary Smart-seq2 cell: its obs and the counts of the named genes (g_<gene>)."""
    X, obs, var = census.counts("mus_musculus", flt=tms_flt(ids), var_names=TMS_GENES, obs_cols=OBS_COLS)
    tab = obs.copy()
    for c in OBS_COLS:
        tab[c] = tab[c].astype(str)
    X = X.tocsc()
    names = list(var["feature_name"].astype(str))
    for g in TMS_GENES:
        j = [k for k, n in enumerate(names) if n == g]
        tab[f"g_{g}"] = np.asarray(X[:, j].sum(axis=1)).ravel() if j else np.zeros(len(tab))
    tab.attrs["missing"] = [g for g in TMS_GENES if g not in names]
    tab.attrs["raw"] = raw_counts(X)
    return tab


def _age(stage: pd.Series) -> np.ndarray:
    return np.where(stage == YOUNG, "young", np.where(stage == OLD, "old", ""))


def case1_rule(tms: pd.DataFrame):
    """Case 1's cell type. A (type, sex) is eligible where >= 3 young and >= 3 old mice of the sex have
    >= 20 cells of the type (the qualifying mice). Among the eligible types, the one with the most
    cells in which Smad3 and Col1a1 are both detected, over the qualifying mice of its eligible sexes;
    ties: more such cells in all, then the type's name. Returns (choice or None, table)."""
    t = tms[units(tms["donor_id"]) & tms["development_stage"].isin([YOUNG, OLD]).to_numpy()].copy()
    t["age"] = _age(t["development_stage"])
    t["unit"] = t["sex"] + "|" + t["age"] + "|" + t["donor_id"]
    per = t.groupby(["cell_type", "sex", "age", "donor_id", "unit"]).size().rename("cells").reset_index()
    q = per[per["cells"] >= MIN_MOUSE_CELLS]
    n = q.groupby(["cell_type", "sex", "age"])["donor_id"].nunique().unstack("age", fill_value=0)
    for a in ("young", "old"):
        if a not in n:
            n[a] = 0
    elig = n[(n["young"] >= MIN_GROUP_MICE) & (n["old"] >= MIN_GROUP_MICE)].reset_index()
    both = ((t[f"g_{PAIR[0]}"] > 0) & (t[f"g_{PAIR[1]}"] > 0)).to_numpy()
    rows = []
    for ct, g in elig.groupby("cell_type", sort=True):
        sexes = sorted(g["sex"])
        mice = q[(q["cell_type"] == ct) & q["sex"].isin(sexes)]
        cells = ((t["cell_type"] == ct) & t["unit"].isin(set(mice["unit"]))).to_numpy()
        per_sex = "; ".join(f"{s}: {int(r['young'])} young, {int(r['old'])} old" for s, r in g.set_index("sex").iterrows())
        rows.append(dict(cell_type=ct, sexes=",".join(sexes), mice=per_sex, cells=int(cells.sum()),
                         codetected=int((cells & both).sum()), _units=sorted(mice["unit"])))
    table = pd.DataFrame(rows, columns=["cell_type", "sexes", "mice", "cells", "codetected", "_units"])
    if table.empty:
        return None, table
    table = table.sort_values(["codetected", "cells", "cell_type"], ascending=[False, False, True],
                              kind="mergesort").reset_index(drop=True)
    top = table.iloc[0]
    choice = dict(cell_type=top["cell_type"], sexes=top["sexes"].split(","), units=list(top["_units"]),
                  cells=int(top["cells"]), codetected=int(top["codetected"]))
    return choice, table


def sex_call(xs: float, ys: float, rule: str = "threshold") -> str:
    """The sex that a mouse's shares of cells with Xist (xs) and with a Y gene (ys) call, or "neither".
    "threshold" (the plan): female where xs >= 0.50 and ys <= 0.05, male the reverse. "contrast" (the
    second run, fixed before it): female where ys <= 0.05 and xs is at least ten times ys and at least
    0.05; male the reverse."""
    if rule == "threshold":
        return ("female" if xs >= CALL_MIN and ys <= CALL_MAX else
                "male" if ys >= CALL_MIN and xs <= CALL_MAX else "neither")
    if rule == "contrast":
        def over(own, other):
            return other <= CALL_SHARE and own >= CALL_RATIO * other and own >= CALL_SHARE
        return "female" if over(xs, ys) else "male" if over(ys, xs) else "neither"
    raise ValueError(f"unknown rule {rule!r}")


def case2_structure(tms: pd.DataFrame, rule: str = "threshold") -> pd.DataFrame:
    """Case 2's structure, per mouse of the 3-month stratum (all its cells): cells, the share with
    Xist > 0, the share with any Y gene > 0, the annotated sex, the tissues, the call by the sex
    genes (sex_call with `rule`) and whether the mouse is kept (a call that matches the annotation)."""
    t = tms[units(tms["donor_id"]) & (tms["development_stage"] == YOUNG).to_numpy()].copy()
    t["xist"] = t["g_Xist"] > 0
    t["y"] = t[[f"g_{g}" for g in Y_GENES]].sum(axis=1) > 0
    rows = []
    for donor, g in t.groupby("donor_id", sort=True):
        sexes = sorted(set(g["sex"]))
        xs, ys = float(g["xist"].mean()), float(g["y"].mean())
        call = sex_call(xs, ys, rule)
        label = sexes[0] if len(sexes) == 1 else "/".join(sexes)
        if call == "neither":
            why = "neither call"
        elif call != label:
            why = f"call {call} contradicts the label {label}"
        else:
            why = ""
        rows.append(dict(mouse=donor, sex=label, cells=len(g), xist_share=round(xs, 4), y_share=round(ys, 4),
                         call=call, kept=not why, excluded_because=why,
                         tissues=",".join(sorted(set(g["tissue"])))))
    return pd.DataFrame(rows, columns=["mouse", "sex", "cells", "xist_share", "y_share", "call", "kept",
                                       "excluded_because", "tissues"])


def case1_audit(tms: pd.DataFrame, audit_mice=None):
    """Case 1 of the second run: section 4's setup as the flagship audit reproduced it
    (flagship_audit/audit_tms.py, section E): every cell of `tms` (TMS FACS, primary, Smart-seq2) at
    the 3-month stage or the "20-month-old stage and over" stage, both sexes, every cell type and
    tissue, no cap; the units are the donor ids as the Census gives them, composite ids included.
    Returns the cells, the table of units and, with the audit's per-mouse table (audit_results.json,
    E.mice), whether each unit's sex, age and cells are the audit's (None without it)."""
    t = tms[tms["development_stage"].isin([YOUNG, OLD]).to_numpy()].copy()
    t["age"] = _age(t["development_stage"])
    per = (t.groupby("donor_id", sort=True)
           .agg(sex=("sex", lambda v: "/".join(sorted(set(v)))), age=("age", lambda v: "/".join(sorted(set(v)))),
                cells=("soma_joinid", "size"), tissues=("tissue", "nunique"))
           .reset_index())
    per["composite"] = per["donor_id"].astype(str).str.contains("/", regex=False)
    per = per.sort_values(["sex", "age", "donor_id"], ascending=[True, False, True], kind="mergesort").reset_index(drop=True)
    check = None
    if audit_mice is not None:
        want = {str(m["donor_id"]): (str(m["sex"]), AUDIT_AGE.get(str(m["age"]), str(m["age"])), int(m["n_cells"]))
                for m in audit_mice}
        got = {str(r.donor_id): (str(r.sex), str(r.age), int(r.cells)) for r in per.itertuples()}
        diff = sorted(set(want) ^ set(got)) + sorted(k for k in set(want) & set(got) if want[k] != got[k])
        check = dict(same=not diff, units=len(got), audit_units=len(want),
                     differences=[dict(donor=k, audit=want.get(k), here=got.get(k)) for k in diff])
    return t, per, check


def case2_rule(tms: pd.DataFrame, structure: pd.DataFrame):
    """Case 2's cell type over the kept mice: the type with >= 20 cells in every kept mouse and the
    most cells; if there is none, the type with >= 20 cells in the most mice (at least 3 of each sex),
    then the most cells; ties: the type's name. Returns (choice or None, table, reason)."""
    kept = structure[structure["kept"]]
    per_sex = kept.groupby("sex")["mouse"].nunique()
    if min(int(per_sex.get("female", 0)), int(per_sex.get("male", 0))) < MIN_GROUP_MICE:
        return None, pd.DataFrame(), f"fewer than {MIN_GROUP_MICE} kept mice of a sex: {per_sex.to_dict()}"
    sex_of = dict(zip(kept["mouse"], kept["sex"]))
    t = tms[tms["donor_id"].isin(set(kept["mouse"])) & (tms["development_stage"] == YOUNG).to_numpy()]
    per = t.groupby(["cell_type", "donor_id"]).size().rename("cells").reset_index()
    rows = []
    for ct, g in per.groupby("cell_type", sort=True):
        ok = g[g["cells"] >= MIN_MOUSE_CELLS]
        sexes = pd.Series([sex_of[m] for m in ok["donor_id"]], dtype=object)
        rows.append(dict(cell_type=ct, mice=len(ok), female=int((sexes == "female").sum()),
                         male=int((sexes == "male").sum()), every=len(ok) == len(kept),
                         cells=int(ok["cells"].sum()), _units=sorted(ok["donor_id"])))
    table = pd.DataFrame(rows, columns=["cell_type", "mice", "female", "male", "every", "cells", "_units"])
    full = table[table["every"]].sort_values(["cells", "cell_type"], ascending=[False, True], kind="mergesort")
    if len(full):
        top, rule = full.iloc[0], "at least 20 cells in every kept mouse, the most cells"
    else:
        part = table[(table["female"] >= MIN_GROUP_MICE) & (table["male"] >= MIN_GROUP_MICE)]
        part = part.sort_values(["mice", "cells", "cell_type"], ascending=[False, False, True], kind="mergesort")
        if part.empty:
            return None, table, "no type has at least 20 cells in 3 mice of each sex"
        top, rule = part.iloc[0], "no type in every kept mouse: the most mice with at least 20 cells, then the most cells"
    table = table.sort_values(["every", "mice", "cells", "cell_type"], ascending=[False, False, False, True],
                              kind="mergesort").reset_index(drop=True)
    choice = dict(cell_type=top["cell_type"], units=list(top["_units"]), cells=int(top["cells"]), rule=rule)
    return choice, table, ""


# --------------------------------------------------------------------------- #
# case 3: the same donors in two technical conditions
# --------------------------------------------------------------------------- #
def excluded_reason(datasets: pd.DataFrame) -> pd.Series:
    """Per dataset_id, why case 3 leaves it out ("" when it does not): the collection name, the
    dataset title or the citation names an excluded study, or a DOI is an excluded one."""
    text = (datasets["collection_name"].astype(str) + " | " + datasets["dataset_title"].astype(str)).str.lower()
    cite = (datasets["collection_doi"].astype(str) + " | "
            + (datasets["citation"].astype(str) if "citation" in datasets else "")).str.lower()
    out = []
    for t, c in zip(text, cite):
        why = [f"name '{n}'" for n in EXCLUDED_NAMES if n in t] + [f"DOI {d}" for d in EXCLUDED_DOIS if d in c]
        out.append("; ".join(why))
    return pd.Series(out, index=datasets["dataset_id"].astype(str).to_numpy())


def case3_candidates(counts: pd.DataFrame, datasets: pd.DataFrame):
    """Every qualifying set: one collection, tissue and cell type, two conditions (assay | suspension
    type) and the donors (the same id within the collection, one sex and stage, the same in both
    conditions) with >= 50 cells in each condition; >= 4 such donors. Sorted by the rule: most
    donors, then most cells (the donors' cells in both conditions), then the ids in order. Excluded
    collections are listed with their reason and never chosen. Returns (candidates, donor rows)."""
    ds = datasets.assign(dataset_id=datasets["dataset_id"].astype(str)).set_index("dataset_id")
    c = counts[counts["dataset_id"].isin(ds.index)].copy()
    c = c[units(c["donor_id"])]
    c["collection_id"] = c["dataset_id"].map(ds["collection_id"]).astype(str)
    c["condition"] = c["assay"] + " | " + c["suspension_type"]
    c["sexstage"] = c["sex"] + " | " + c["development_stage"]
    key = ["organism", "collection_id", "tissue", "cell_type", "condition", "donor_id"]
    per = (c.groupby(key).agg(cells=("cells", "sum"), n_sexstage=("sexstage", "nunique"),
                              sexstage=("sexstage", "first"))
           .reset_index())
    per = per[(per["cells"] >= NEG_MIN_CELLS) & (per["n_sexstage"] == 1)]
    on = ["organism", "collection_id", "tissue", "cell_type", "donor_id"]
    m = per.merge(per, on=on, suffixes=("_1", "_2"))
    m = m[(m["condition_1"] < m["condition_2"]) & (m["sexstage_1"] == m["sexstage_2"])].copy()
    m["pair_cells"] = m["cells_1"] + m["cells_2"]
    g = ["organism", "collection_id", "tissue", "cell_type", "condition_1", "condition_2"]
    cand = m.groupby(g).agg(donors=("donor_id", "nunique"), cells=("pair_cells", "sum")).reset_index()
    cand = cand[cand["donors"] >= NEG_MIN_DONORS]
    reason = excluded_reason(datasets)
    coll = datasets.assign(reason=reason.to_numpy()).groupby("collection_id").agg(
        collection_name=("collection_name", "first"), collection_doi=("collection_doi", "first"),
        excluded=("reason", lambda s: "; ".join(sorted({r for r in s if r}))))
    cand = cand.merge(coll, left_on="collection_id", right_index=True, how="left")
    cand = cand.sort_values(["donors", "cells", *g], ascending=[False, False] + [True] * len(g),
                            kind="mergesort").reset_index(drop=True)
    return cand, m


def pooled_pairs(X: np.ndarray, genes: list, idx, k: int) -> list:
    """The `k` most correlated pairs of the genes `idx`, without a shared gene, taken greedily (by r,
    then the columns in order): r is norm_pearson (log1p CP10k with the cell's sum over the kept
    genes, the cells where both genes are detected, Pearson; 0 below 3 such cells or without
    variance) over the pooled cells. Returns [(gene_a, gene_b, r)]."""
    tot = X.sum(axis=1).astype(float)
    tot = np.where(tot == 0, 1.0, tot)
    sub = X[:, idx].astype(float)
    D = (sub > 0).astype(float)
    L = np.log1p(sub / tot[:, None] * 1e4) * D
    # sums over the cells by einsum's own loops, not BLAS: OpenBLAS's dgemm of these transposed
    # products crashed on a GitHub runner with v1's OPENBLAS_CORETYPE=Haswell (the dry run of 0cd81db)
    N = np.einsum("ci,cj->ij", D, D, optimize=False)
    SA = np.einsum("ci,cj->ij", L, D, optimize=False)
    SAA = np.einsum("ci,cj->ij", L * L, D, optimize=False)
    SAB = np.einsum("ci,cj->ij", L, L, optimize=False)
    va = N * SAA - SA ** 2
    vb = va.T
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (N * SAB - SA * SA.T) / np.sqrt(va * vb)
    tiny_a = va <= 1e-12 * np.maximum(N * SAA, 1e-300)
    r = np.where((N >= 3) & ~tiny_a & ~tiny_a.T & np.isfinite(r), np.clip(r, -1.0, 1.0), 0.0)  # as np.corrcoef
    pairs = sorted(((float(r[i, j]), i, j) for i, j in itertools.combinations(range(len(idx)), 2)),
                   key=lambda x: (-x[0], x[1], x[2]))
    used, out = set(), []
    for val, i, j in pairs:
        if i in used or j in used:
            continue
        used |= {i, j}
        out.append((genes[idx[i]], genes[idx[j]], val))
        if len(out) == k:
            break
    return out


def case3_claims(X: np.ndarray, genes: list) -> tuple[list, dict]:
    """Case 3's 20 two-sided claims on the pooled cells: (3a) norm_pearson on the 8 most correlated
    pairs among the 100 most expressed genes (by mean CP10k), (3b) the same among the genes ranked
    501 to 600, (3c) the mean log1p CP10k of each of the 4 most expressed genes."""
    tot = X.sum(axis=1).astype(float)
    tot = np.where(tot == 0, 1.0, tot)
    mean_cp10k = (X / tot[:, None] * 1e4).mean(axis=0)
    order = [int(j) for j in np.argsort(-mean_cp10k, kind="mergesort")]
    claims, facts = [], dict(top_genes=[genes[j] for j in order[:10]])
    for tag, (lo, hi) in (("3a", CASE3_TOP), ("3b", CASE3_MID)):
        idx = order[lo:hi]
        if len(idx) < hi - lo:
            facts[f"{tag}_missing"] = f"only {len(order)} genes: ranks {lo + 1}-{hi} incomplete"
            if len(idx) < 2:
                continue
        for k, (a, b, r) in enumerate(pooled_pairs(X, genes, idx, CASE3_PAIRS), 1):
            claims.append(dict(id=f"{tag}{k}", metric="norm_pearson", genes=[a, b], pooled_r=r,
                               direction="two-sided", sesoi=SESOI_PAIR,
                               signal=dict(kind="coupling", strength=COUPLING_STRENGTH)))
    for k, j in enumerate(order[:CASE3_MEANS], 1):
        claims.append(dict(id=f"3c{k}", metric="log_cp10k_mean", genes=[genes[j]],
                           mean_cp10k=float(mean_cp10k[j]), direction="two-sided", sesoi=SESOI_MEAN,
                           signal=dict(kind="module", fold=MODULE_FOLD, frac=MODULE_FRAC)))
    return claims, facts


# --------------------------------------------------------------------------- #
# extraction
# --------------------------------------------------------------------------- #
def present_genes(census, organism: str, datasets: pd.DataFrame, ids) -> np.ndarray:
    """The var joinids that the presence matrix marks as measured in every dataset of `ids`."""
    rows = datasets.set_index(datasets["dataset_id"].astype(str)).loc[list(ids), "soma_joinid"].to_numpy()
    M = census.presence(organism)[rows]
    every = np.asarray((M > 0).sum(axis=0)).ravel() == len(rows)
    return np.where(every)[0]


def read_blocks(census, organism: str, joinids: np.ndarray, present: np.ndarray, chunk: int):
    """The counts of `joinids` over the genes `present`, read `chunk` cells at a time (sorted joinids
    in sorted blocks) and stacked: int32 CSR (the counts must be raw), obs, var. None for the
    counts where a block is not raw."""
    from scipy import sparse
    blocks = np.array_split(np.asarray(joinids, dtype=np.int64), max(1, int(np.ceil(len(joinids) / chunk))))
    Xs, obs, var = [], [], None
    for b in blocks:
        Xb, ob, vb = census.counts(organism, obs_joinids=b, var_joinids=present, obs_cols=OBS_COLS)
        if Xb.shape[0] != len(b) or not np.array_equal(ob["soma_joinid"].to_numpy(np.int64), b):
            raise ValueError(f"{len(b)} cells asked, {Xb.shape[0]} read, or not those")
        if var is None:
            var = vb
        elif not np.array_equal(vb["soma_joinid"].to_numpy(), var["soma_joinid"].to_numpy()):
            raise ValueError("the genes read differ between blocks")
        if not raw_counts(Xb):
            return None, None, var
        Xs.append(sparse.csr_matrix(Xb).astype(np.int32))
        obs.append(ob)
    return sparse.vstack(Xs, format="csr"), pd.concat(obs, ignore_index=True), var


def extract(census, organism: str, datasets: pd.DataFrame, joinids: np.ndarray, dataset_ids, named=(),
            chunk: int | None = None) -> dict:
    """The cells' counts over the genes measured in every dataset they come from (`dataset_ids`), cut
    by D1's rule; with `chunk`, read that many cells at a time (read_blocks)."""
    ids = sorted(set(map(str, dataset_ids)))
    var = census.var(organism)
    present = present_genes(census, organism, datasets, ids)
    if len(present) and present.max() >= len(var):
        raise ValueError("the presence matrix has more genes than the var table")
    if chunk:
        X, obs, v = read_blocks(census, organism, joinids, present, chunk)
        if X is None:
            return dict(raw=False, cells=int(len(joinids)), datasets=ids)
    else:
        X, obs, v = census.counts(organism, obs_joinids=joinids, var_joinids=present, obs_cols=OBS_COLS)
    if X.shape[0] != len(joinids):
        raise ValueError(f"{len(joinids)} cells asked, {X.shape[0]} read")
    if sorted(set(obs["dataset_id"].astype(str))) != ids:
        raise ValueError("the cells read come from other datasets than their obs said")
    if not raw_counts(X):
        return dict(raw=False, cells=int(X.shape[0]), datasets=ids)
    if not chunk:
        X = X.astype(np.int64)
    dense, genes, missing = cut_genes(X, list(v["feature_name"].astype(str)), named)
    for c in OBS_COLS:
        obs[c] = obs[c].astype(str)
    return dict(raw=True, X=dense, genes=genes, obs=obs, missing=missing, datasets=ids,
                measured_genes=int(len(present)), var_genes=int(len(var)))


def write_case(out_dir: Path, name: str, rec: dict, obs: pd.DataFrame) -> dict:
    path = out_dir / f"{name}.npz"
    P.save_npz(path, rec["X"], obs, rec["genes"])
    return dict(file=path.name, sha256=P.sha256(path), cells=int(rec["X"].shape[0]),
                genes=int(rec["X"].shape[1]), datasets=rec["datasets"], measured_genes=rec["measured_genes"],
                var_genes=rec["var_genes"], missing_named=rec["missing"])


def per_unit_table(obs: pd.DataFrame, cols) -> list[str]:
    n = obs.groupby(list(cols)).size().rename("cells").reset_index()
    return md_table(n, max_rows=60)


# --------------------------------------------------------------------------- #
# the cases
# --------------------------------------------------------------------------- #
def run_case1(census, tms, datasets, out_dir, log) -> dict:
    log("## Case 1: MI coupling on Tabula Muris Senis (FACS), young (3 months) against old "
        f"({OLD!r}), within sex")
    choice, table = case1_rule(tms)
    log("Eligible cell types (a sex is eligible with >= 3 young and >= 3 old mice of >= 20 cells), by the rule:")
    log(*md_table(table, ["cell_type", "sexes", "mice", "cells", "codetected"]))
    if choice is None:
        log("No cell type is eligible: case 1 does not run.")
        return dict(status="not run", reason="no eligible cell type")
    t = tms[(tms["cell_type"] == choice["cell_type"]).to_numpy()].copy()
    t["age"] = _age(t["development_stage"])
    t["unit"] = t["sex"] + "|" + t["age"] + "|" + t["donor_id"]
    t = t[t["unit"].isin(set(choice["units"]))]
    ids = sample_ids({u: g["soma_joinid"].to_numpy() for u, g in t.groupby("unit")}, MAX_CELLS)
    log(f"Chosen: {choice['cell_type']!r} within {', '.join(choice['sexes'])}; {len(choice['units'])} mice, "
        f"{choice['cells']} cells ({choice['codetected']} with Smad3 and Col1a1 both detected); "
        f"extracted {len(ids)} cells (at most {MAX_CELLS} per mouse)")
    rec = extract(census, "mus_musculus", datasets, ids, t.loc[t["soma_joinid"].isin(ids), "dataset_id"],
                  named=PAIR + POS_PAIR + NEG_PAIR)
    if not rec["raw"]:
        log("The counts are not raw (non-negative integers): case 1 does not run.")
        return dict(status="not run", reason="counts not raw", choice=_public(choice))
    obs = rec["obs"]
    obs = pd.DataFrame(dict(mouse=obs["donor_id"], age=_age(obs["development_stage"]), sex=obs["sex"],
                            tissue=obs["tissue"], cell_type=obs["cell_type"], dataset_id=obs["dataset_id"],
                            development_stage=obs["development_stage"], soma_joinid=obs["soma_joinid"].astype(np.int64)))
    log("Cells per mouse (extracted):", *per_unit_table(obs, ["sex", "age", "mouse"]))
    log("Tissues of the extracted cells, per sex and age:",
        *per_unit_table(obs, ["sex", "age", "tissue"]))
    f = write_case(out_dir, "case1", rec, obs)
    log(f"case1.npz: {f['cells']} cells x {f['genes']} genes ({f['measured_genes']} of {f['var_genes']} genes "
        f"measured in all {len(f['datasets'])} datasets read); named genes missing: {f['missing_named'] or 'none'}; "
        f"sha256 {f['sha256']}")
    missing = set(f["missing_named"])
    claims = []
    for cid, metric in (("1a", "mi_3bin"), ("1b", "norm_pearson")):
        claims.append(dict(id=cid, metric=metric, genes=list(PAIR), direction="decrease", sesoi=0.1,
                           signal=dict(kind="coupling", strength=COUPLING_STRENGTH),
                           controls=dict(pos=list(POS_PAIR), neg=list(NEG_PAIR)),
                           skip=(f"named genes missing: {sorted(missing & set(PAIR + POS_PAIR + NEG_PAIR))}"
                                 if missing & set(PAIR + POS_PAIR + NEG_PAIR) else "")))
    return dict(status="run", choice=_public(choice), data=f, group_col="age", groups=["young", "old"],
                within=["sex"], replicate="mouse", claims=claims, eligible=_records(table))


def run_case1_audit(census, tms, datasets, out_dir, log, audit=None, chunk: int = CHUNK_CELLS) -> dict:
    """Case 1 of the second run: section 4's setup as the flagship audit reproduced it (case1_audit)."""
    log("## Case 1, second run: section 4's setup as the flagship audit reproduced it: every TMS FACS cell at "
        f"{YOUNG!r} and {OLD!r}, both sexes, every cell type and tissue, no cap per mouse; the units are the donor "
        "ids as the Census gives them, composite ids included; young against old within sex")
    t, per, check = case1_audit(tms, audit["E"]["mice"] if audit else None)
    log("Units (donor ids) by sex and age:", *md_table(per, max_rows=60))
    n = per.groupby(["sex", "age"]).size()
    log("Units per sex, young against old: " + "; ".join(
        f"{sx} {int(n.get((sx, 'young'), 0))} against {int(n.get((sx, 'old'), 0))}" for sx in sorted(set(per["sex"]))))
    if check is not None:
        log(f"The units and their cells are the audit's section E (audit_results.json, E.mice): {check['same']} "
            f"({check['units']} units here, {check['audit_units']} there)"
            + "".join(f"; {d['donor']}: audit {d['audit']}, here {d['here']}" for d in check["differences"]))
    mixed = per[per["sex"].str.contains("/", regex=False) | per["age"].str.contains("/", regex=False)]
    if len(mixed):
        log(f"Donor ids with two sexes or two ages: {list(mixed['donor_id'])}: case 1 does not run.")
        return dict(status="not run", reason="a donor id with two sexes or ages", units=_records(per), audit_check=check)
    ids = np.sort(t["soma_joinid"].to_numpy(np.int64))
    log(f"Extracting {len(ids)} cells of {t['dataset_id'].nunique()} datasets, {chunk} at a time")
    rec = extract(census, "mus_musculus", datasets, ids, t["dataset_id"], named=PAIR + POS_PAIR + NEG_PAIR, chunk=chunk)
    if not rec["raw"]:
        log("The counts are not raw (non-negative integers): case 1 does not run.")
        return dict(status="not run", reason="counts not raw", units=_records(per), audit_check=check)
    o = rec["obs"]
    obs = pd.DataFrame(dict(mouse=o["donor_id"], age=_age(o["development_stage"]), sex=o["sex"], tissue=o["tissue"],
                            cell_type=o["cell_type"], dataset_id=o["dataset_id"], development_stage=o["development_stage"],
                            soma_joinid=o["soma_joinid"].astype(np.int64)))
    log("Cells per mouse (extracted):", *per_unit_table(obs, ["sex", "age", "mouse"]))
    log("Tissues of the extracted cells, per sex and age:", *md_table(
        obs.groupby(["sex", "age", "tissue"]).size().rename("cells").reset_index(), max_rows=100))
    f = write_case(out_dir, "case1", rec, obs)
    log(f"case1.npz: {f['cells']} cells x {f['genes']} genes ({f['measured_genes']} of {f['var_genes']} genes "
        f"measured in all {len(f['datasets'])} datasets read); named genes missing: {f['missing_named'] or 'none'}; "
        f"sha256 {f['sha256']}")
    missing = set(f["missing_named"]) & set(PAIR + POS_PAIR + NEG_PAIR)
    skip = f"named genes missing: {sorted(missing)}" if missing else ""
    base = dict(genes=list(PAIR), direction="decrease", signal=dict(kind="coupling", strength=COUPLING_STRENGTH),
                controls=dict(pos=list(POS_PAIR), neg=list(NEG_PAIR)), skip=skip)
    claims = [dict(base, id="1a", metric="mi_3bin", sesoi=0.1,
                   sesoi_note="AGENTS.md's example, wrong in scale for mi_3bin"),
              dict(base, id="1a-floor", metric="mi_3bin", sesoi="floor", sesoi_scale="observed", seed_id="1a",
                   sesoi_note="the design's floor by v0.4's change 4, computed by real_pilot.py from (b) before (d)"),
              dict(base, id="1b", metric="norm_pearson", sesoi=0.1)]
    out = dict(status="run", run="second", setup="section 4 as the flagship audit reproduced it", data=f,
               group_col="age", groups=["young", "old"], within=["sex"], replicate="mouse", claims=claims,
               units=_records(per), audit_check=check)
    if audit:
        out["audit_mice"] = [{k: m[k] for k in ("donor_id", "sex", "age", "n_cells", "median_nnz", "mi_3bin")}
                             for m in audit["E"]["mice"]]
        out["audit_tests"] = audit["E"]["tests"]
    return out


def run_case2(census, tms, datasets, out_dir, log, rule: str = "threshold") -> dict:
    log("## Case 2: Xist and the Y genes on Tabula Muris Senis (FACS), female against male, 3 months")
    st = case2_structure(tms, rule)
    how = (f"female: Xist in >= {CALL_MIN:.0%} of cells and a Y gene in <= {CALL_MAX:.0%}; male: the reverse"
           if rule == "threshold" else
           f"by contrast: female where a Y gene is in <= {CALL_SHARE:.0%} of cells and Xist in >= {CALL_RATIO:g} times "
           f"that share and >= {CALL_SHARE:.0%}; male: the reverse")
    log(f"Structure per mouse of the 3-month stratum (all its cells), checked before any metric ({how}):")
    log(*md_table(st, max_rows=60))
    choice, table, why = case2_rule(tms, st)
    if len(table):
        log("Cell types over the kept mice:", *md_table(table, ["cell_type", "mice", "female", "male", "every", "cells"]))
    if choice is None:
        log(f"Case 2 does not run: {why}.")
        return dict(status="not run", reason=why, structure=_records(st), sex_rule=rule)
    t = tms[(tms["cell_type"] == choice["cell_type"]).to_numpy() & tms["donor_id"].isin(set(choice["units"])).to_numpy()
            & (tms["development_stage"] == YOUNG).to_numpy()]
    ids = sample_ids({u: g["soma_joinid"].to_numpy() for u, g in t.groupby("donor_id")}, MAX_CELLS)
    log(f"Chosen: {choice['cell_type']!r} ({choice['rule']}); {len(choice['units'])} mice, {choice['cells']} cells; "
        f"extracted {len(ids)} cells (at most {MAX_CELLS} per mouse)")
    rec = extract(census, "mus_musculus", datasets, ids, t.loc[t["soma_joinid"].isin(ids), "dataset_id"],
                  named=SEX_GENES)
    if not rec["raw"]:
        log("The counts are not raw (non-negative integers): case 2 does not run.")
        return dict(status="not run", reason="counts not raw", choice=_public(choice), structure=_records(st))
    obs = rec["obs"]
    obs = pd.DataFrame(dict(mouse=obs["donor_id"], sex=obs["sex"], age=obs["development_stage"], tissue=obs["tissue"],
                            cell_type=obs["cell_type"], dataset_id=obs["dataset_id"],
                            soma_joinid=obs["soma_joinid"].astype(np.int64)))
    log("Cells per mouse (extracted):", *per_unit_table(obs, ["sex", "mouse"]))
    f = write_case(out_dir, "case2", rec, obs)
    log(f"case2.npz: {f['cells']} cells x {f['genes']} genes ({f['measured_genes']} of {f['var_genes']} genes "
        f"measured in all {len(f['datasets'])} datasets read); named genes missing: {f['missing_named'] or 'none'}; "
        f"sha256 {f['sha256']}")
    missing = set(f["missing_named"])
    claims = [dict(id="2a", metric="log_cp10k_mean", genes=["Xist"], direction="decrease", sesoi=0.5,
                   signal=dict(kind="module", fold=MODULE_FOLD, frac=MODULE_FRAC),
                   skip="Xist missing" if "Xist" in missing else ""),
              dict(id="2b", metric="log_cp10k_mean", genes=[g for g in Y_GENES if g not in missing],
                   direction="increase", sesoi=0.5, signal=dict(kind="module", fold=MODULE_FOLD, frac=MODULE_FRAC),
                   skip="no Y gene" if set(Y_GENES) <= missing else "")]
    return dict(status="run", choice=_public(choice), data=f, group_col="sex", groups=["female", "male"],
                within=[], replicate="mouse", claims=claims, structure=_records(st), types=_records(table), sex_rule=rule)


def run_case3(census, datasets, out_dir, log) -> dict:
    log("## Case 3: a negative control, the same donors in two technical conditions")
    counts, filters = [], {}
    for organism in ORGANISMS:
        flt = dict(is_primary_data=True, disease="normal")
        if "tissue_type" in census.obs_columns(organism):
            flt["tissue_type"] = "tissue"
        filters[organism] = soma_filter(flt)
        t0 = time.time()
        c = count_cells(census.obs_batches(organism, flt, CASE3_COLS), CASE3_COLS)
        c.insert(0, "organism", organism)
        counts.append(c)
        log(f"{organism}: {int(c['cells'].sum()):,} cells counted ({filters[organism]}) in {time.time() - t0:.0f} s")
    counts = pd.concat(counts, ignore_index=True)
    cand, donors = case3_candidates(counts, datasets)
    log("Qualifying sets (>= 4 donors with >= 50 cells in both conditions), by the rule; excluded ones are not chosen:")
    log(*md_table(cand.assign(excluded=cand["excluded"].replace("", "-")),
                  ["organism", "collection_name", "tissue", "cell_type", "condition_1", "condition_2", "donors",
                   "cells", "excluded"]))
    ok = cand[cand["excluded"] == ""].reset_index(drop=True)
    if ok.empty:
        log("No set qualifies: case 3 does not run.")
        return dict(status="not run", reason="no qualifying set", filters=filters)
    tried = []
    for rank in range(1, min(len(ok), MAX_RANKS) + 1):
        ch = ok.iloc[rank - 1]
        keys = ["organism", "collection_id", "tissue", "cell_type", "condition_1", "condition_2"]
        d = donors[np.logical_and.reduce([(donors[k] == ch[k]).to_numpy() for k in keys])]
        qualifying = sorted(set(d["donor_id"]))
        qual = qualifying
        if len(qual) > NEG_MAX_DONORS:  # a sample of the sorted ids, drawn with the seed
            qual = sorted(np.random.default_rng(SEED).choice(qual, size=NEG_MAX_DONORS, replace=False).tolist())
        coll_ids = datasets.loc[datasets["collection_id"].astype(str) == ch["collection_id"], "dataset_id"].astype(str).tolist()
        (a1, s1), (a2, s2) = ch["condition_1"].split(" | "), ch["condition_2"].split(" | ")
        flt = dict(dataset_id=coll_ids, is_primary_data=True, disease="normal", tissue=ch["tissue"],
                   cell_type=ch["cell_type"], assay=sorted({a1, a2}), suspension_type=sorted({s1, s2}))
        if "tissue_type" in census.obs_columns(ch["organism"]):
            flt["tissue_type"] = "tissue"
        obs = census.obs_table(ch["organism"], flt, OBS_COLS)
        for c in OBS_COLS:
            obs[c] = obs[c].astype(str)
        obs["condition"] = obs["assay"] + " | " + obs["suspension_type"]
        sel = obs[(obs["tissue"] == ch["tissue"]) & (obs["cell_type"] == ch["cell_type"])
                  & obs["condition"].isin([ch["condition_1"], ch["condition_2"]]) & obs["donor_id"].isin(qual)]
        ids = sample_ids({f"{dn}|{cd}": g["soma_joinid"].to_numpy() for (dn, cd), g in sel.groupby(["donor_id", "condition"])},
                         MAX_CELLS)
        log(f"Rank {rank}: {ch['collection_name']!r} ({ch['collection_id']}, DOI {ch['collection_doi']}), "
            f"{ch['organism']}, {ch['tissue']!r}, {ch['cell_type']!r}: A = {ch['condition_1']!r}, "
            f"B = {ch['condition_2']!r}; {len(qualifying)} qualifying donors"
            + (f", {len(qual)} of them sampled (at most {NEG_MAX_DONORS})" if len(qual) < len(qualifying) else "")
            + f"; extracting {len(ids)} cells")
        rec = extract(census, ch["organism"], datasets, ids, sel.loc[sel["soma_joinid"].isin(ids), "dataset_id"])
        if rec["raw"]:
            break
        tried.append(dict(rank=rank, collection_id=ch["collection_id"], reason="counts not raw"))
        log(f"Rank {rank}: the counts are not raw (non-negative integers): the next set is tried")
    else:
        log(f"None of the first {MAX_RANKS} sets has raw counts: case 3 does not run.")
        return dict(status="not run", reason="no set with raw counts", tried=tried, filters=filters)
    o = rec["obs"]
    o["condition"] = o["assay"] + " | " + o["suspension_type"]
    obs = pd.DataFrame(dict(donor=o["donor_id"], group=np.where(o["condition"] == ch["condition_1"], "A", "B"),
                            condition=o["condition"], assay=o["assay"], suspension_type=o["suspension_type"],
                            tissue=o["tissue"], cell_type=o["cell_type"], sex=o["sex"],
                            development_stage=o["development_stage"], dataset_id=o["dataset_id"],
                            soma_joinid=o["soma_joinid"].astype(np.int64)))
    log("Cells per donor and condition (extracted):", *per_unit_table(obs, ["donor", "group", "sex", "development_stage"]))
    f = write_case(out_dir, "case3", rec, obs)
    log(f"case3.npz: {f['cells']} cells x {f['genes']} genes ({f['measured_genes']} of {f['var_genes']} genes "
        f"measured in all {len(f['datasets'])} datasets read); sha256 {f['sha256']}")
    claims, facts = case3_claims(rec["X"], rec["genes"])
    log("Claims, chosen on the pooled cells of both conditions (two-sided; no difference exists by construction):")
    for c in claims:
        extra = f"pooled r {c['pooled_r']:.4f}" if "pooled_r" in c else f"mean CP10k {c['mean_cp10k']:.2f}"
        log(f"  {c['id']}: {c['metric']}({', '.join(c['genes'])}); {extra}")
    for k, v in facts.items():
        log(f"  {k}: {v}")
    choice = {k: (v.item() if hasattr(v, "item") else v) for k, v in ch.items()}
    return dict(status="run", choice=choice, rank=rank, tried=tried, filters=filters, data=f, group_col="group",
                groups=["A", "B"], within=[], replicate="donor", conditions=dict(A=ch["condition_1"], B=ch["condition_2"]),
                donors=qual, qualifying_donors=len(qualifying), claims=claims, facts=facts,
                candidates=_records(cand.head(30)))


def _public(choice: dict) -> dict:
    return {k: v for k, v in choice.items() if not k.startswith("_")}


def _records(df: pd.DataFrame) -> list:
    return json.loads(df.drop(columns=[c for c in df.columns if str(c).startswith("_")]).to_json(orient="records"))


# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--synthetic", action="store_true", help="a fabricated Census (offline): the code path only")
    p.add_argument("--audit-results", type=Path,
                   help="the audit's out/audit_results.json: its TMS dataset ids must be the ones resolved here")
    p.add_argument("--second", action="store_true",
                   help="the second run (the plan's deviation of 2026-10-16): cases 1 and 2 only, case 1 as the audit "
                        "reproduced section 4, case 2 with the sex called by contrast")
    args = p.parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    lines = []

    def log(*msg):
        for m in msg:
            print(m, flush=True)
            lines.append(str(m))
    tag = "SYNTHETIC " if args.synthetic else ""
    t0 = time.time()
    log(f"# {tag}The pilot's extraction (validation/exploratory/real_pilot_plan.md)"
        + (", second run (the plan's section \"Deviation: a second run of cases 1 and 2\")" if args.second else ""), "")
    census = SyntheticCensus() if args.synthetic else Census()
    spec = dict(synthetic=args.synthetic, release="SYNTHETIC" if args.synthetic else RELEASE, seed=SEED,
                rules=dict(max_cells=MAX_CELLS, n_genes=N_GENES, min_mouse_cells=MIN_MOUSE_CELLS,
                           min_group_mice=MIN_GROUP_MICE, neg_min_donors=NEG_MIN_DONORS, neg_min_cells=NEG_MIN_CELLS,
                           neg_max_donors=NEG_MAX_DONORS,
                           excluded_names=list(EXCLUDED_NAMES), excluded_dois=list(EXCLUDED_DOIS)),
                python=platform.python_version(), run="second" if args.second else "first", cases={})
    if args.second:
        spec["rules"].update(sex_call="contrast", call_share=CALL_SHARE, call_ratio=CALL_RATIO, chunk_cells=CHUNK_CELLS,
                             case1="section 4 as the flagship audit reproduced it")
    failed = []
    try:
        spec["census_summary"] = census.summary()
        log(f"{tag}Census release {spec['release']}: {spec['census_summary']}")
        datasets = census.datasets()
        ids = tms_dataset_ids(datasets)
        spec["tms_dataset_ids"] = ids
        audit = json.loads(args.audit_results.read_text()) if args.audit_results else None
        if audit:
            same = sorted(audit["meta"]["tms_dataset_ids"]) == sorted(ids)
            spec["tms_ids_as_audit"] = same
            log(f"TMS FACS datasets: {len(ids)}; the same as the audit's: {same}")
            if not same:
                raise SystemExit("the TMS FACS datasets differ from the audit's")
        t1 = time.time()
        tms = read_tms(census, ids)
        log(f"TMS FACS: {len(tms):,} primary Smart-seq2 cells of {len(ids)} datasets with the named genes read in "
            f"{time.time() - t1:.0f} s; named genes missing: {tms.attrs['missing'] or 'none'}; raw: {tms.attrs['raw']}")
        stages = tms.groupby(["sex", "development_stage"]).agg(cells=("soma_joinid", "size"),
                                                               mice=("donor_id", "nunique")).reset_index()
        log("Cells and donor ids per sex and stage (composite ids included here):", *md_table(stages, max_rows=40))
        if args.second:  # a small block on the fabricated Census, so that the blocks are stacked there too
            cases = (("case1", lambda: run_case1_audit(census, tms, datasets, args.out_dir, log, audit,
                                                       chunk=1000 if args.synthetic else CHUNK_CELLS)),
                     ("case2", lambda: run_case2(census, tms, datasets, args.out_dir, log, rule="contrast")))
        else:
            cases = (("case1", lambda: run_case1(census, tms, datasets, args.out_dir, log)),
                     ("case2", lambda: run_case2(census, tms, datasets, args.out_dir, log)),
                     ("case3", lambda: run_case3(census, datasets, args.out_dir, log)))
        for name, fn in cases:
            t1 = time.time()
            log("")
            try:
                spec["cases"][name] = fn()
            except Exception as exc:  # printed, and the analyses do not start
                tb = traceback.format_exc()
                log(f"{name}: EXCEPTION {exc!r}", *tb.splitlines())
                spec["cases"][name] = dict(status="error", error=repr(exc), traceback=tb.splitlines())
                failed.append(name)
            spec["cases"][name]["seconds"] = round(time.time() - t1, 1)
            log(f"{name}: {spec['cases'][name]['status']} ({spec['cases'][name]['seconds']} s)")
    finally:
        census.close()
    spec["seconds"] = round(time.time() - t0, 1)
    log("", f"{tag}extraction: {spec['seconds']} s in all; failed cases: {failed or 'none'}")
    (args.out_dir / "select.json").write_text(json.dumps(spec, indent=1, default=str))
    (args.out_dir / "real_pilot_select.log").write_text("\n".join(lines) + "\n")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
