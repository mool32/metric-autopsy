#!/usr/bin/env python3
"""Level B, step S5 (`plan_s5.md`): one kept paper's data, read from the authors' deposit, and its checks.

Run in GitHub Actions, once per dataset, by the label `level-b-data` (`.github/workflows/level-b-data.yml`), on the
dataset's spec (`data_specs/<rid>.json`, committed before its run). Nothing of the engine is imported: the checks
are `data_checks.py`'s. The spec names the deposit's files, how to read them, which cells are the claim's, and
which metadata give each cell's group, donor, source and sex label:

    {"rid": "R01225", "organism": "homo_sapiens", "loader": "csv_counts",
     "files": {"counts": "<url>", "meta": "<url>"},
     "cells": {"location": ["PT", "FT"]},               # optional: keep the cells whose metadata have these values
     "groups": {"column": "group", "values": ["tumour", "normal"],
                "rules": [{"when": {"location": ["PT"], "Cell_Type": ["Epithelial cells"]}, "value": "tumour"}, ...]},
     "donor": "orig.ident", "source": null, "sex_label": null,
     "claims": [{"id": "1", "genes": ["SERPINE2", "TFF3"]}, ...]}

A spec may add columns to the cells' metadata from a table the paper publishes beside its deposit (`tables`: the
table's file among `files`, a CSV or an Excel sheet read with the standard library, its key and value columns, and
the cells' key, a column or a regular expression's first group on one), e.g. each sample's patient:

     "tables": [{"as": "patient", "file": "pairs", "sheet": "Sheet1", "key": "Sample_ID", "value": "Patient.ID",
                 "on": {"column": "title", "regex": "(\\d+)$"}}]

The cells of the two groups that no donor names are left out of the checks and listed (`without_donor`). A paper
that kept only the donors with at least n cells in each group (`"min_cells_each_group": n`) has the others left
out and listed (`donors_left_out`), as the paper did.

Loaders (`LOADERS`):
* `csv_counts` — a counts table as CSV (genes by cells or cells by genes, found from the metadata's cell ids) and a
  table of cells as CSV (the cell ids in its first column);
* `tenx_samples` — a GEO archive of per-sample 10x files (barcodes, genes, matrix) and the series matrix, whose
  per-sample titles and characteristics become every cell's metadata; with `loader_args` `{"min_umis": n}` (a
  deposit of Cell Ranger's raw matrices, every barcode) a cell is a barcode with at least n counts, the authors'
  own threshold;
* `mtx_libraries` — one large MatrixMarket file (genes by cells) read in a stream, its barcodes prefixed by their
  library's title, and the series matrices, whose per-library characteristics become the cells' metadata; only the
  X-inactivation and Y genes and the claims' genes are kept, with every cell's total;
* `r_export` — the counts, genes and cells R wrote from the deposit's R object (`r_export.R`, run before this);
  with `loader_args` `{"symbol_column": c}` the genes are named by the object's rowData column c.

    python -I validation/level_b/data_run.py --spec validation/level_b/data_specs/R01225.json --out OUT
    python -I validation/level_b/data_run.py --spec ... --out OUT --files DIR   # the files already in DIR (no fetch)
    python -I validation/level_b/data_run.py --spec ... --out OUT --files-to DIR --fetch-only   # an R object, fetched
                                                                                  # before r_export.R exports it

OUT gets report.json and report.md (the six checks of `plan_s5.md`), files.tsv (each file's URL, size and
SHA256) and the log.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
import tarfile
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import data_checks as D  # noqa: E402
from controls import SEX_GENES  # noqa: E402

USER_AGENT = "metric-autopsy-validation/0.4 (level B data checks; python urllib)"


# --------------------------------------------------------------------------- #
# files
# --------------------------------------------------------------------------- #
def fetch(url: str, path: Path, log) -> dict:
    """One file to `path` (streamed), up to 4 attempts on a network error or HTTP 5xx; its size and SHA256."""
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            h = hashlib.sha256()
            with urllib.request.urlopen(req, timeout=600) as r, path.open("wb") as f:
                while True:
                    block = r.read(1 << 22)
                    if not block:
                        break
                    h.update(block)
                    f.write(block)
            return dict(url=url, file=path.name, bytes=path.stat().st_size, sha256=h.hexdigest(), attempts=attempt + 1)
        except urllib.error.HTTPError as e:
            if e.code < 500:
                raise
            log(f"  attempt {attempt + 1}: HTTP {e.code}: {url}")
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
            log(f"  attempt {attempt + 1}: {type(e).__name__}: {url}")
        time.sleep(2 ** (attempt + 1))
    raise RuntimeError(f"four attempts failed: {url}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def _open(path: Path):
    return gzip.open(path, "rb") if path.suffix == ".gz" else path.open("rb")


# --------------------------------------------------------------------------- #
# metadata from GEO's series matrix
# --------------------------------------------------------------------------- #
def series_samples(path: Path) -> pd.DataFrame:
    """One row per sample of a series matrix: its accession, title, source and each characteristic ("key: value")
    as a column."""
    rows = {}
    with _open(path) as f:
        for raw in f:
            line = raw.decode("utf-8", "replace").rstrip("\n")
            if not line.startswith("!Sample_"):
                continue
            key, *vals = line.split("\t")
            vals = [v.strip('"') for v in vals]
            rows.setdefault(key, []).append(vals)
    acc = rows["!Sample_geo_accession"][0]
    df = pd.DataFrame(dict(gsm=acc, title=rows.get("!Sample_title", [[""] * len(acc)])[0],
                           source=rows.get("!Sample_source_name_ch1", [[""] * len(acc)])[0]))
    for vals in rows.get("!Sample_characteristics_ch1", []):
        for i, v in enumerate(vals):
            if ":" in v:
                k, val = v.split(":", 1)
                df.loc[i, k.strip().lower()] = val.strip()
    return df


# --------------------------------------------------------------------------- #
# a table the paper publishes (its samples' patients, say): CSV, TSV or an Excel sheet
# --------------------------------------------------------------------------- #
_XL = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
       "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
       "p": "http://schemas.openxmlformats.org/package/2006/relationships"}


def _col_index(ref: str) -> int:
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group(0):
        n = n * 26 + ord(ch) - 64
    return n - 1


def xlsx_rows(path: Path, sheet: str | None = None) -> list:
    """One sheet of an Excel workbook (.xlsx), read with the standard library (no Excel reader in the run's
    environment): its rows, each a list over every column the sheet uses (None where a cell is empty), every value
    a string as the workbook stores it ("69", "pat. 5"); rows without a value dropped."""
    import xml.etree.ElementTree as ET
    import zipfile
    m = "{%s}" % _XL["m"]
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        sheets = [(e.get("name"), e.get("{%s}id" % _XL["r"])) for e in wb.iter(m + "sheet")]
        rid = next((r for n, r in sheets if sheet is None or n == sheet), None)
        if rid is None:
            raise ValueError(f"{path.name}: no sheet {sheet!r} (it has {[n for n, _ in sheets]})")
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        target = next(e.get("Target") for e in rels.iter("{%s}Relationship" % _XL["p"]) if e.get("Id") == rid)
        target = target.lstrip("/")
        target = target if target.startswith("xl/") else "xl/" + target
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).iter(m + "si"):
                shared.append("".join(t.text or "" for t in si.iter(m + "t")))
        root = ET.fromstring(z.read(target))
    rows = []
    for row in root.iter(m + "row"):
        vals = {}
        for i, c in enumerate(row.iter(m + "c")):
            j = _col_index(c.get("r")) if c.get("r") else i
            v = c.find(m + "v")
            if c.get("t") == "s" and v is not None:
                vals[j] = shared[int(v.text)]
            elif c.get("t") == "inlineStr":
                vals[j] = "".join(t.text or "" for t in c.iter(m + "t"))
            elif v is not None and v.text is not None:
                vals[j] = v.text
        if vals:
            rows.append(vals)
    width = max((max(r) + 1 for r in rows), default=0)
    return [[r.get(j) for j in range(width)] for r in rows]


def xlsx_sheet(path: Path, sheet: str | None = None) -> pd.DataFrame:
    """One sheet as a table: its first row the header (an unnamed column named by its place), every value a
    string as the workbook stores it, rows without a value dropped."""
    rows = xlsx_rows(path, sheet)
    if not rows:
        return pd.DataFrame()
    header = [h if h is not None else f"column {j + 1}" for j, h in enumerate(rows[0])]
    return pd.DataFrame(rows[1:], columns=header)


def paper_table(path: Path, sheet: str | None = None) -> pd.DataFrame:
    if path.suffix.lower() == ".xlsx":
        return xlsx_sheet(path, sheet)
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def table_column(obs: pd.DataFrame, path: Path, rule: dict) -> pd.Series:
    """Each cell's value from the paper's table: the cell's key (a column, or a regular expression's first group on
    one) looked up in the table's key column; a key the table does not hold gives no value."""
    t = paper_table(path, rule.get("sheet")).dropna(subset=[rule["key"], rule["value"]])
    t = t.assign(_k=t[rule["key"]].astype(str).str.strip(), _v=t[rule["value"]].astype(str).str.strip())
    clash = t.groupby("_k")["_v"].nunique()
    if (clash > 1).any():
        raise ValueError(f"{path.name}: keys with two values {sorted(clash[clash > 1].index)}")
    lut = dict(zip(t["_k"], t["_v"]))
    on = rule["on"]
    keys = (obs[on["column"]].astype(str).str.extract(on["regex"], expand=False) if isinstance(on, dict)
            else obs[on].astype(str))
    return keys.map(lut)


# --------------------------------------------------------------------------- #
# loaders: (X cells x genes CSR, gene symbols, obs) or, for the streamed one, a summary of the kept genes
# --------------------------------------------------------------------------- #
def load_csv_counts(files: dict, args: dict, log):
    from scipy import sparse
    meta = pd.read_csv(files["meta"], index_col=0)
    meta.index = meta.index.astype(str)
    counts = pd.read_csv(files["counts"], index_col=0)
    cols, rows = set(map(str, counts.columns)), set(map(str, counts.index))
    ids = set(meta.index)
    if len(ids & cols) >= len(ids & rows):  # genes by cells
        counts.columns = counts.columns.astype(str)
        counts = counts.T
    counts.index = counts.index.astype(str)
    common = [c for c in meta.index if c in set(counts.index)]
    log(f"csv_counts: {counts.shape[0]} cells by {counts.shape[1]} genes; {len(common)} of the table's {len(meta)} "
        "cells in the counts")
    X = sparse.csr_matrix(counts.loc[common].to_numpy(dtype=np.float64))
    return X, [str(g) for g in counts.columns], meta.loc[common].reset_index().rename(columns={"index": "cell"})


def _read_tenx(barcodes: bytes, genes: bytes, matrix: bytes):
    from scipy import sparse
    from controls import read_mtx
    bc = [ln for ln in gzip.decompress(barcodes).decode().splitlines() if ln.strip()]
    gn = [ln.split("\t") for ln in gzip.decompress(genes).decode().splitlines() if ln.strip()]
    M = read_mtx(matrix)
    if M.shape == (len(gn), len(bc)):
        M = M.T.tocsr()
    return sparse.csr_matrix(M), [g[1] if len(g) > 1 else g[0] for g in gn], bc


def load_tenx_samples(files: dict, args: dict, log):
    """A GEO RAW archive of `GSMxxxx_<name>_{barcodes,genes,matrix}.*.gz` triplets; the cells' metadata are their
    sample's (`series`)."""
    from scipy import sparse
    samples = series_samples(files["series"])
    min_umis = int(args.get("min_umis", 0))
    by_gsm = {}
    with tarfile.open(files["raw"]) as tar:
        for m in tar.getmembers():
            g = re.match(r"(GSM\d+)_.*?(barcodes|genes|features|matrix)\.", m.name)
            if g:
                by_gsm.setdefault(g.group(1), {})["genes" if g.group(2) == "features" else g.group(2)] = \
                    tar.extractfile(m).read()
    Xs, obs, genes = [], [], None
    for gsm in sorted(by_gsm):
        parts = by_gsm[gsm]
        if set(parts) != {"barcodes", "genes", "matrix"}:
            log(f"tenx_samples: {gsm} lacks {sorted({'barcodes', 'genes', 'matrix'} - set(parts))}")
            continue
        X, g, bc = _read_tenx(parts["barcodes"], parts["genes"], parts["matrix"])
        if min_umis:
            tot = np.asarray(X.sum(axis=1)).ravel()
            keep = np.flatnonzero(tot >= min_umis)
            log(f"tenx_samples: {gsm}: {len(keep)} of {len(bc)} barcodes with at least {min_umis} counts")
            X, bc = X[keep], [bc[i] for i in keep]
        if genes is None:
            genes = g
        elif g != genes:
            raise ValueError(f"{gsm}: its genes are not the first sample's")
        Xs.append(X)
        obs.append(pd.DataFrame(dict(cell=[f"{gsm}_{b}" for b in bc], gsm=gsm)))
    obs = pd.concat(obs, ignore_index=True).merge(samples, on="gsm", how="left")
    log(f"tenx_samples: {len(Xs)} samples, {sum(x.shape[0] for x in Xs)} cells, {len(genes)} genes")
    return sparse.vstack(Xs, format="csr"), genes, obs


def load_r_export(files: dict, args: dict, log):
    """The files `r_export.R` wrote into `args["export_dir"]` from the deposit's R object: counts.mtx.gz (genes by
    cells), genes.tsv.gz, cells.csv.gz (the object's metadata, the cells in the counts' order)."""
    from scipy import sparse
    from controls import read_mtx
    d = Path(args["export_dir"])
    M = read_mtx((d / "counts.mtx.gz").read_bytes())
    genes = [ln for ln in gzip.decompress((d / "genes.tsv.gz").read_bytes()).decode().splitlines() if ln.strip()]
    if args.get("symbol_column"):  # the genes' symbols from the object's rowData (`rowdata.csv.gz`), in the same order
        rd = pd.read_csv(d / "rowdata.csv.gz", dtype=str)
        if rd["gene"].tolist() != genes:
            raise ValueError("rowdata.csv.gz's genes are not genes.tsv.gz's")
        genes = rd[args["symbol_column"]].fillna("").tolist()
    obs = pd.read_csv(d / "cells.csv.gz")
    if M.shape == (len(genes), len(obs)):
        M = M.T
    log(f"r_export: {M.shape[0]} cells by {M.shape[1]} genes")
    return sparse.csr_matrix(M), genes, obs


def stream_mtx_libraries(files: dict, args: dict, keep_genes: list, log) -> dict:
    """A large genes-by-cells MatrixMarket file read in chunks: every value checked as a non-negative integer, every
    cell's total, and the values of `keep_genes` only. Returns the kept genes' CSR (cells x kept), the cells'
    totals, whether every value is a count, and the cells' metadata from their libraries (`series`)."""
    from scipy import sparse
    feats = [ln.split("\t") for ln in gzip.decompress(files["features"].read_bytes()).decode().splitlines() if ln.strip()]
    symbols = [f[1] if len(f) > 1 else f[0] for f in feats]
    bcd = [ln for ln in gzip.decompress(files["barcodes"].read_bytes()).decode().splitlines() if ln.strip()]
    want = {g.upper() for g in keep_genes}
    kept = [i for i, s in enumerate(symbols) if s.upper() in want]
    lut = np.full(len(symbols), -1, dtype=np.int64)
    lut[kept] = np.arange(len(kept))
    kept_names = [symbols[i] for i in kept]
    totals = np.zeros(len(bcd))
    rows, cols, vals = [], [], []
    integers, nonneg, n_rows = True, True, 0
    header = None
    reader = pd.read_csv(files["matrix"], sep=" ", comment="%", header=None, chunksize=20_000_000, dtype=np.float64)
    for chunk in reader:
        a = chunk.to_numpy()
        if header is None:
            header, a = a[0], a[1:]
            if int(header[0]) != len(symbols) or int(header[1]) != len(bcd):
                raise ValueError(f"the matrix is {header[:2]}, the features {len(symbols)} and barcodes {len(bcd)}")
        n_rows += len(a)
        v = a[:, 2]
        integers &= bool(np.all(v == np.round(v)))
        nonneg &= bool(np.all(v >= 0))
        g, c = a[:, 0].astype(np.int64) - 1, a[:, 1].astype(np.int64) - 1
        totals += np.bincount(c, weights=v, minlength=len(bcd))
        k = lut[g]
        sel = k >= 0
        if sel.any():
            rows.append(c[sel])
            cols.append(k[sel])
            vals.append(v[sel])
        log(f"  {n_rows:,} entries read")
    if header is not None and n_rows != int(header[2]):
        raise ValueError(f"the matrix declares {int(header[2])} entries and holds {n_rows}")
    K = sparse.csr_matrix((np.concatenate(vals) if vals else [], (np.concatenate(rows) if rows else [],
                           np.concatenate(cols) if cols else [])), shape=(len(bcd), len(kept_names)))
    samples = pd.concat([series_samples(p) for p in files["series"]], ignore_index=True)
    lib = [b.rsplit("_", 1)[0] for b in bcd]
    obs = pd.DataFrame(dict(cell=bcd, library=lib)).merge(samples.rename(columns={"title": "library"}), on="library",
                                                          how="left")
    return dict(K=K, genes=kept_names, totals=totals, raw=integers and nonneg, obs=obs, all_genes=symbols)


LOADERS = dict(csv_counts=load_csv_counts, tenx_samples=load_tenx_samples, r_export=load_r_export)


# --------------------------------------------------------------------------- #
# the spec's cells, groups, donors
# --------------------------------------------------------------------------- #
def _match(obs: pd.DataFrame, when: dict) -> np.ndarray:
    m = np.ones(len(obs), bool)
    for col, vals in when.items():
        if col not in obs.columns:
            raise ValueError(f"the metadata have no column {col!r} (they have {list(obs.columns)[:20]})")
        m &= obs[col].astype(str).isin([str(v) for v in vals]).to_numpy()
    return m


def _labelled(obs: pd.DataFrame, spec: dict):
    """Every cell with its group (a column, or rules in order), donor (a column, or a regular expression's first
    group on a column), source and sex label; and which cells are the spec's (its cells, in one of its groups)."""
    obs = obs.copy()
    keep = _match(obs, spec.get("cells") or {})
    g = spec["groups"]
    if g.get("rules"):
        grp = np.array([None] * len(obs), dtype=object)
        for rule in g["rules"]:
            m = _match(obs, rule["when"]) & (grp == None)  # noqa: E711 - the first rule that matches
            grp[m] = rule["value"]
        obs["group"] = grp
    else:
        obs["group"] = obs[g["column"]].astype(str).map(g.get("map") or (lambda x: x))
    keep &= obs["group"].isin(g["values"]).to_numpy()
    def column(rule):  # a column, or {"column": ..., "regex": ...}: the first group of the pattern on it; a missing
        # value stays missing (pandas 2's astype(str) would make it the string "nan")
        c = obs[rule["column"] if isinstance(rule, dict) else rule]
        c = c.astype(str).where(c.notna())
        return c.str.extract(rule["regex"], expand=False) if isinstance(rule, dict) else c
    obs["donor"] = column(spec["donor"])
    obs["source"] = column(spec["source"]).fillna("unmatched") if spec.get("source") else "no source named"
    obs["sex_label"] = column(spec["sex_label"]).str.lower() if spec.get("sex_label") else "-"
    return obs, keep


def assign(obs: pd.DataFrame, spec: dict) -> pd.DataFrame:
    """The spec's cells whose donor is named, each with its group, donor and source."""
    obs, keep = _labelled(obs, spec)
    return obs[keep & obs["donor"].notna().to_numpy()]


def min_cells_filter(obs: pd.DataFrame, spec: dict):
    """The spec's `min_cells_each_group`: the donors with at least that many cells in each of the two groups, and
    the others (left out, with their cells per group)."""
    n = spec.get("min_cells_each_group")
    if not n:
        return obs, {}
    groups = [str(g) for g in spec["groups"]["values"]]
    ct = obs.groupby(["donor", "group"]).size().unstack(fill_value=0).reindex(columns=groups, fill_value=0)
    ok = set(ct.index[(ct >= n).all(axis=1)])
    out = {str(d): {g: int(ct.loc[d, g]) for g in groups} for d in ct.index if d not in ok}
    return obs[obs["donor"].isin(ok).to_numpy()], out


def without_donor(obs: pd.DataFrame, spec: dict) -> dict:
    """The spec's cells that no donor names (left out of the checks), by sample (its title) where the metadata
    have one, and by group."""
    obs, keep = _labelled(obs, spec)
    o = obs[keep & obs["donor"].isna().to_numpy()]
    by = ["title", "group"] if "title" in o.columns else ["group"]
    return {" / ".join(map(str, k if isinstance(k, tuple) else (k,))): int(v) for k, v in o.groupby(by).size().items()}


def checks(X, genes: list, obs: pd.DataFrame, spec: dict, raw: dict) -> dict:
    groups = tuple(spec["groups"]["values"])
    res = dict(rid=spec["rid"], counts=raw, cells=int(len(obs)),
               cells_per_group={k: int(v) for k, v in obs["group"].value_counts().items()},
               donors=D.donors_check(obs, "donor", "group", groups),
               sex=D.sex_check(X, genes, obs.reset_index(drop=True), "donor", spec["organism"],
                               label_col="sex_label" if spec.get("sex_label") else None),
               genes={c["id"]: D.genes_check(genes, c["genes"]) for c in spec["claims"]},
               source=(D.source_check(obs, "group", groups, "source") if spec.get("source")
                       else dict(verdict="no source named")))
    return res


def add_tables(obs: pd.DataFrame, spec: dict, paths: dict, log) -> None:
    """The spec's `tables`: each adds a column to the cells' metadata from a table the paper publishes."""
    for t in spec.get("tables") or []:
        obs[t["as"]] = table_column(obs, paths[t["file"]], t).to_numpy()
        log(f"table {t['file']}: {t['as']} for {int(obs[t['as']].notna().sum())} of {len(obs)} cells")


def report_md(res: dict, spec: dict, files: list) -> str:
    dn = res["donors"]
    lines = [f"# Level B, S5: {spec['rid']}'s data and their checks (data_run.py; plan_s5.md)", "",
             "Files: " + "; ".join(f"{f['file']} {f['bytes']:,} bytes sha256 {f['sha256']}" for f in files), "",
             f"1. Counts: non-negative integers {res['counts'].get('raw')}.",
             f"2. Donors: {dn['design']}, per group {dn['donors']}; ids naming several donors {dn['composite'] or 'none'}; "
             f"rule 4 {'passes' if dn['rule4'].get('passes') else 'fails'} (smallest p {dn['rule4'].get('min_p', '-')})"
             + (f"; cells no donor names, left out: {res['without_donor']}" if res.get("without_donor") else "")
             + (f"; donors with fewer than {spec.get('min_cells_each_group')} cells in a group, left out as the paper did: "
                f"{res['donors_left_out']}" if res.get("donors_left_out") else "") + ".",
             f"3. Each id one individual: both sexes among the donors {res['sex']['both_sexes']}; failing "
             f"{res['sex']['fails'] or 'none'}.",
             "4. The claims' genes: " + "; ".join(f"{k} {'present' if v['passes'] else 'missing ' + ', '.join(v['missing'])}"
                                                  for k, v in res["genes"].items()) + ".",
             f"5. The source: {res['source']['verdict']}" + (f" (sources {res['source'].get('sources')})"
                                                             if res['source'].get('sources') else "") + ".",
             f"6. The paper's design (from the spec): {spec.get('design', '-')}.",
             "", f"Cells: {res['cells']} ({res['cells_per_group']})."]
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--spec", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--files", type=Path, help="the spec's files already in this directory, by their names")
    ap.add_argument("--files-to", type=Path, help="where to download the spec's files (default OUT/files)")
    ap.add_argument("--fetch-only", action="store_true", help="download the files and stop (an R object before its export)")
    a = ap.parse_args(argv)
    spec = json.loads(a.spec.read_text())
    a.out.mkdir(parents=True, exist_ok=True)
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    log(f"# Level B, S5: {spec['rid']} ({spec.get('accession', '-')}), loader {spec['loader']}")
    fdir = a.files or a.files_to or (a.out / "files")
    fdir.mkdir(parents=True, exist_ok=True)
    got, paths = [], {}
    for key, url in spec["files"].items():
        urls = url if isinstance(url, list) else [url]
        ps = []
        for u in urls:
            p = fdir / u.rsplit("/", 1)[-1]
            rec = (dict(url=u, file=p.name, bytes=p.stat().st_size, sha256=sha256(p), attempts=0) if a.files
                   else fetch(u, p, log))
            got.append(rec)
            ps.append(p)
        paths[key] = ps if isinstance(url, list) else ps[0]
    pd.DataFrame(got).to_csv(a.out / "files.tsv", sep="\t", index=False)
    if a.fetch_only:
        (a.out / "data_run.log").write_text("\n".join(lines) + "\n")
        return 0
    claim_genes = sorted({g for c in spec["claims"] for g in c["genes"]})
    if spec["loader"] == "mtx_libraries":
        xg, yg = SEX_GENES[spec["organism"]]
        s = stream_mtx_libraries(paths, spec.get("loader_args") or {}, [xg, *yg, *claim_genes], log)
        add_tables(s["obs"], spec, paths, log)
        missing = without_donor(s["obs"], spec)
        obs, left = min_cells_filter(assign(s["obs"], spec), spec)
        X = s["K"][obs.index.to_numpy()]
        res = checks(X, s["genes"], obs.reset_index(drop=True), spec, dict(raw=s["raw"]))
        res["genes"] = {c["id"]: D.genes_check(s["all_genes"], c["genes"]) for c in spec["claims"]}
    else:
        X, genes, obs = LOADERS[spec["loader"]](paths, spec.get("loader_args") or {}, log)
        add_tables(obs, spec, paths, log)
        missing = without_donor(obs, spec)
        obs, left = min_cells_filter(assign(obs, spec), spec)
        X = X[obs.index.to_numpy()]
        res = checks(X, genes, obs.reset_index(drop=True), spec, D.counts_check(X))
    res["without_donor"] = missing
    res["donors_left_out"] = left
    res["files"] = got
    (a.out / "report.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    md = report_md(res, spec, got)
    (a.out / "report.md").write_text(md)
    log(md)
    (a.out / "data_run.log").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
