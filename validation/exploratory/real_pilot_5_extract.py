#!/usr/bin/env python3
"""The extraction of the third diagnostics and of v0.4's acceptance tests
(validation/exploratory/real_pilot_diagnostics_3_plan.md). Exploratory, development data; it does not
import the engine (v1's selection environment, requirements-select.txt).

- **The full set of genes** for the cells of cases 1 and 2 (the second run's files) and 3 (the first
  run's), read again from the CELLxGENE Census (release 2025-11-08) by their soma_joinid, 10,000 cells at
  a time, checked raw (real_pilot_select.read_blocks): every gene the Census marks as measured in every
  dataset of the case, symbols that name more than one feature left out. The case files' genes must
  have the same counts. Saved in the case file's order of cells.
- **The metadata of case 1** for check (b4): every obs column of the Census for its cells, the datasets
  table's row of its dataset, and the obs of the dataset's source h5ad (its matrix is never read).
- **Marsh et al. 2022**: GSE152183_RAW.tar (downloaded by the workflow) checked against GEO's file list,
  its 12 samples against the registered table (real_pilot_5_spec.MARSH_SAMPLES) and GEO's SOFT
  records, and merged into one matrix of every barcode (Cell Ranger's called cells), symbols that name
  more than one feature left out.

    python real_pilot_5_extract.py --data INPUTS/data --geo INPUTS/geo --out OUT
    python real_pilot_5_extract.py --data DATA --geo GEO --out OUT --synthetic   # offline, fabricated

With --synthetic the same code runs on a fabricated Census (real_pilot_select's, with the list's
genes added) and fabricated GEO files (--fabricate-geo DIR writes them first): every output says
SYNTHETIC. --check-census opens the Census and counts the cases' measured genes, reading no count.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import platform
import sys
import tarfile
import tempfile
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "prereg"))
import panel as P  # noqa: E402  (save_npz, sha256: numpy and pandas, no engine)
import real_pilot_5_spec as S  # noqa: E402
import real_pilot_select as RS  # noqa: E402  (Census, SyntheticCensus, present_genes, read_blocks)

CHUNK = 10000
ORGANISM = {"case1": "mus_musculus", "case2": "mus_musculus", "case3": "homo_sapiens"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------- #
# the Census, with what this run needs beside real_pilot_select's
# --------------------------------------------------------------------------- #
class Census(RS.Census):
    def obs_rows(self, organism: str, joinids, cols=None) -> pd.DataFrame:
        """Every obs column (or `cols`) of the cells `joinids`, in the order of `joinids`."""
        j = np.asarray(joinids, dtype=np.int64)
        t = self._obs(organism).read(coords=(np.sort(j),), column_names=cols).concat().to_pandas()
        return t.set_index("soma_joinid").loc[j].reset_index()

    def source_obs(self, dataset_id: str, workdir: Path) -> tuple[pd.DataFrame, dict]:
        """The obs of the dataset's source h5ad, as the release names it; X is never read."""
        import anndata
        path = Path(workdir) / f"{dataset_id}.h5ad"
        self._cc.download_source_h5ad(dataset_id, to_path=str(path), census_version=S.RELEASE, progress_bar=False)
        info = dict(file=path.name, bytes=path.stat().st_size, sha256=sha256_file(path))
        ad = anndata.read_h5ad(path, backed="r")
        obs = ad.obs.copy()
        info.update(cells=int(ad.n_obs), genes=int(ad.n_vars), obs_index_name=str(obs.index.name))
        ad.file.close()
        path.unlink()  # 'only metadata': the file goes as soon as its obs is read
        return obs.reset_index(), info


class SyntheticCensus(RS.SyntheticCensus):
    """real_pilot_select's fabricated Census with the list's genes added (fabricated counts, a matrix
    module in the fibroblast-like types and a little everywhere else), a fabricated source h5ad's obs
    and the Census's numeric obs columns. Nothing is real."""

    EXTRA = tuple(g for g in S.FIB_GENES if g not in RS.TMS_GENES)

    def _build_mouse(self):
        super()._build_mouse()
        obs, X = self._obs["mus_musculus"], self._X["mus_musculus"]
        rng = np.random.default_rng(11)  # its own draws: the base counts stay those of real_pilot_select's
        fib = obs["cell_type"].astype(str).str.contains("|".join(S.FIB_WORDS), case=False).to_numpy()
        young = (obs["development_stage"] == RS.YOUNG).to_numpy()
        lib = X.sum(axis=1) / max(X.sum(axis=1).mean(), 1.0)
        lam = np.outer(lib, np.ones(len(self.EXTRA))) * np.where(fib, 20.0, 0.05 + 0.25 * young)[:, None]
        extra = rng.poisson(lam).astype(np.int32)
        self._X["mus_musculus"] = np.hstack([X, extra])
        var = self._var["mus_musculus"]
        n0 = len(var)
        add = pd.DataFrame(dict(soma_joinid=np.arange(n0, n0 + len(self.EXTRA)),
                                feature_id=[f"ENSMUSG9{i:010d}" for i in range(len(self.EXTRA))],
                                feature_name=list(self.EXTRA)))
        self._var["mus_musculus"] = pd.concat([var, add], ignore_index=True)
        obs["raw_sum"] = self._X["mus_musculus"].sum(axis=1).astype(float)
        obs["nnz"] = (self._X["mus_musculus"] > 0).sum(axis=1).astype(np.int64)
        obs["observation_joinid"] = [f"obs{i:06d}" for i in range(len(obs))]

    def obs_rows(self, organism, joinids, cols=None) -> pd.DataFrame:
        obs = self._obs[organism].set_index("soma_joinid")
        t = obs.loc[np.asarray(joinids, dtype=np.int64)].reset_index()
        return t if cols is None else t[list(cols)]

    def source_obs(self, dataset_id, workdir) -> tuple[pd.DataFrame, dict]:
        obs = self._obs["mus_musculus"]
        sel = obs[obs["assay"] == "Smart-seq2"].copy()
        out = pd.DataFrame({"index": [f"A{i % 12 + 1}.B{100 + i // 300:06d}.{d}.1.1"
                                      for i, d in enumerate(sel["donor_id"].astype(str))],
                            "mouse.id": sel["donor_id"].astype(str).to_numpy(),
                            "age": np.where(sel["development_stage"] == RS.YOUNG, "3m", "24m"),
                            "FACS.selection": "Multiple", "method": "facs",
                            "donor_id": sel["donor_id"].astype(str).to_numpy(),
                            "development_stage": sel["development_stage"].astype(str).to_numpy()})
        return out, dict(file="SYNTHETIC", bytes=0, sha256="SYNTHETIC", cells=int(len(out)), genes=0,
                         obs_index_name="index")


# --------------------------------------------------------------------------- #
# the full set of genes for a case's cells
# --------------------------------------------------------------------------- #
def load_case(path: Path, want: str | None):
    got = sha256_file(path)
    if want is not None and got != want:
        raise ValueError(f"{path.name}: sha256 {got}, expected {want}")
    z = np.load(path, allow_pickle=False)
    obs = pd.DataFrame({k[4:]: z[k] for k in z.files if k.startswith("obs_")})
    return np.asarray(z["X"], dtype=np.int32), [str(g) for g in z["genes"]], obs, got


def unique_columns(symbols: list) -> np.ndarray:
    """The columns whose symbol names one feature only (the others are left out)."""
    s = pd.Series(symbols)
    return np.where(~s.duplicated(keep=False).to_numpy())[0]


def full_genes(census, organism, datasets, obs_case: pd.DataFrame, dataset_ids, chunk: int):
    """The case's cells (in the case file's order) x every measured gene with a unique symbol. Read as
    real_pilot_select.read_blocks reads (sorted joinids in sorted blocks, each checked raw), the columns
    kept block by block, so that the full matrix is held about twice at most."""
    from scipy import sparse
    joinids = obs_case["soma_joinid"].to_numpy(np.int64)
    order = np.argsort(joinids, kind="mergesort")
    ids = sorted(set(map(str, dataset_ids)))
    present = RS.present_genes(census, organism, datasets, ids)
    blocks = np.array_split(joinids[order], max(1, int(np.ceil(len(joinids) / chunk))))
    Xs, keep, symbols, var0 = [], None, None, None
    for b in blocks:
        Xb, ob, vb = census.counts(organism, obs_joinids=b, var_joinids=present, obs_cols=RS.OBS_COLS)
        if Xb.shape[0] != len(b) or not np.array_equal(ob["soma_joinid"].to_numpy(np.int64), b):
            raise ValueError(f"{len(b)} cells asked, {Xb.shape[0]} read, or not those")
        if var0 is None:
            var0 = vb["soma_joinid"].to_numpy()
            symbols = vb["feature_name"].astype(str).tolist()
            keep = unique_columns(symbols)
        elif not np.array_equal(vb["soma_joinid"].to_numpy(), var0):
            raise ValueError("the genes read differ between blocks")
        if not RS.raw_counts(Xb):
            raise ValueError("the counts read are not raw (non-negative integers)")
        Xs.append(sparse.csr_matrix(Xb)[:, keep].astype(np.int32))
        del Xb
    X = sparse.vstack(Xs, format="csr")
    del Xs
    inv = np.empty_like(order)
    inv[order] = np.arange(len(order))
    X = X[inv]  # back to the case file's order
    kept = set(keep.tolist())
    left_out = sorted({symbols[j] for j in range(len(symbols)) if j not in kept})
    return X.tocsr(), [symbols[j] for j in keep], dict(measured_genes=int(len(present)),
                                                       kept_genes=int(len(keep)), left_out=left_out)


def same_counts(X_full, genes_full, X_case, genes_case, names) -> dict:
    col = {g: j for j, g in enumerate(genes_full)}
    out = {}
    for g in names:
        if g not in col or g not in genes_case:
            out[g] = None
            continue
        a = np.asarray(X_full[:, col[g]].todense()).ravel().astype(np.int64)
        out[g] = bool(np.array_equal(a, X_case[:, genes_case.index(g)].astype(np.int64)))
    return out


def run_census(args, log, out) -> None:
    census = SyntheticCensus() if args.synthetic else Census()
    out["census_summary"] = census.summary()
    datasets = census.datasets()
    spec2 = json.loads((args.data / "select.json").read_text())
    spec1 = json.loads((args.data / "select1.json").read_text())
    cases = {"case1": spec2["cases"]["case1"], "case2": spec2["cases"]["case2"], "case3": spec1["cases"]["case3"]}
    named = {"case1": [*S.PAIR, *S.FIB_GENES, *S.S4_POS, *S.S4_NEG], "case2": ["Xist", *S.Y_GENES], "case3": []}
    for name, case in cases.items():
        t0 = time.time()
        path = args.data / case["data"]["file"]
        Xc, genes_c, obs_c, sha = load_case(path, None if args.synthetic else S.CASE_SHA256[path.name])
        X, genes, info = full_genes(census, ORGANISM[name], datasets, obs_c, case["data"]["datasets"], args.chunk)
        check = same_counts(X, genes, Xc, genes_c, sorted(set(genes_c) & set(genes))[:50] + named[name])
        bad = {g: v for g, v in check.items() if v is False}
        log(f"## {name}: {X.shape[0]} cells x {X.shape[1]} genes ({info['measured_genes']} measured, "
            f"{len(info['left_out'])} symbols naming several features left out; {time.time() - t0:.0f} s); "
            f"the case file's genes with the same counts: {sum(v is True for v in check.values())} of "
            f"{sum(v is not None for v in check.values())} checked" + (f"; DIFFERENT: {bad}" if bad else ""))
        missing = [g for g in named[name] if g not in genes]
        if missing:
            log(f"   named genes not measured: {missing}")
        if bad:
            raise ValueError(f"{name}: the counts differ from {path.name}'s: {bad}")
        f = args.out / f"full_{name}.npz"
        P.save_npz(f, X, obs_c[["soma_joinid"]], genes)
        out["cases"][name] = dict(case_file=path.name, case_sha256=sha, cells=int(X.shape[0]), **info,
                                  checked=check, named_missing=missing, file=f.name, sha256=sha256_file(f),
                                  seconds=round(time.time() - t0, 1))
        if name == "case1":
            t1 = time.time()
            cols = census.obs_columns("mus_musculus")
            co = census.obs_rows("mus_musculus", obs_c["soma_joinid"].to_numpy(np.int64), cols)
            co.to_csv(args.out / "census_obs_case1.csv.gz", index=False)
            row = datasets[datasets["dataset_id"].astype(str).isin(case["data"]["datasets"])]
            (args.out / "dataset_case1.json").write_text(row.astype(str).to_json(orient="records", indent=1))
            with tempfile.TemporaryDirectory(dir=args.out) as tmp:
                src, sinfo = census.source_obs(S.TMS_DATASET, Path(tmp))
            src.to_csv(args.out / "source_obs_case1.csv.gz", index=False)
            out["metadata_case1"] = dict(census_obs_columns=list(map(str, cols)), census_obs_rows=int(len(co)),
                                         datasets=row["dataset_id"].astype(str).tolist(), source_h5ad=sinfo,
                                         source_obs_columns=list(map(str, src.columns)),
                                         seconds=round(time.time() - t1, 1))
            log(f"   the metadata: {len(cols)} Census obs columns for {len(co)} cells; the source h5ad "
                f"({sinfo['bytes']:,} bytes, sha256 {sinfo['sha256']}): {sinfo['cells']} cells, obs columns "
                f"{list(map(str, src.columns))}")
        del X
    census.close()


# --------------------------------------------------------------------------- #
# Marsh et al. 2022
# --------------------------------------------------------------------------- #
def soft_samples(path: Path) -> dict:
    """GSM -> {title, characteristics {key: value}} from a GEO family SOFT file (text only)."""
    out, cur = {}, None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                cur = out.setdefault(line.split("=", 1)[1].strip(), dict(title=None, ch={}))
            elif line.startswith("^"):
                cur = None
            elif cur is not None and " = " in line:
                k, v = line.split(" = ", 1)
                if k == "!Sample_title":
                    cur["title"] = v.strip()
                elif k.startswith("!Sample_characteristics") and ": " in v:
                    a, b = v.split(": ", 1)
                    cur["ch"][a.strip().lower()] = b.strip()
    return out


def check_marsh_soft(path: Path) -> list:
    """The registered table against GEO's records: the problems found (none expected)."""
    got, bad = soft_samples(path), []
    if set(got) != set(S.MARSH_SAMPLES):
        bad.append(f"GSMs: GEO {sorted(got)}, registered {sorted(S.MARSH_SAMPLES)}")
    for gsm, (title, arm, batch) in S.MARSH_SAMPLES.items():
        g = got.get(gsm)
        if g is None:
            continue
        digestion, inhib = S.ARM_CHARACTERISTICS[arm]
        want = {"digestion type": digestion, "inhibitors": inhib, "batch": batch, "sex": "M"}
        if g["title"] != title:
            bad.append(f"{gsm}: title {g['title']!r}, registered {title!r}")
        for k, v in want.items():
            if g["ch"].get(k) != v:
                bad.append(f"{gsm}: {k} {g['ch'].get(k)!r}, registered {v!r}")
    return bad


def filelist(path: Path) -> dict:
    """GEO's suppl/filelist.txt: name -> size in bytes."""
    out = {}
    for line in path.read_text().splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) >= 4:
            out[parts[1]] = int(parts[3])
    return out


def read_10x(tar: tarfile.TarFile, members: dict, gsm: str, title: str):
    """One sample's Cell Ranger files from the tar's members (read as bytes, nothing written)."""
    from scipy import io as sio
    from scipy import sparse
    stem = f"{gsm}_{title}_"

    def raw(suffix):
        return gzip.decompress(tar.extractfile(members[stem + suffix]).read())
    barcodes = raw("barcodes.tsv.gz").decode().split()
    feats = [line.split("\t") for line in raw("features.tsv.gz").decode().splitlines() if line]
    M = sio.mmread(io.BytesIO(raw("matrix.mtx.gz")), spmatrix=False)
    M = sparse.csr_matrix(M).T.tocsr()  # cells x genes
    if M.shape != (len(barcodes), len(feats)):
        raise ValueError(f"{gsm}: matrix {M.shape}, {len(barcodes)} barcodes and {len(feats)} features")
    if M.nnz and (M.data.min() < 0 or not np.all(M.data == np.round(M.data))):
        raise ValueError(f"{gsm}: the matrix is not raw counts")
    return M.astype(np.int32), barcodes, [f[0] for f in feats], [f[1] if len(f) > 1 else f[0] for f in feats]


def run_marsh(args, log, out) -> None:
    from scipy import sparse
    t0 = time.time()
    bad = check_marsh_soft(args.geo / "GSE152183_family.soft.gz")
    if bad:
        raise ValueError("the registered sample table differs from GEO's records: " + "; ".join(bad))
    log("## Marsh et al. 2022 (GSE152183): the registered table of 12 samples as GEO's records")
    tar_path = args.geo / S.MARSH_TAR
    size, listed = tar_path.stat().st_size, filelist(args.geo / "GSE152183_suppl_filelist.txt")
    if not args.synthetic and size != S.MARSH_TAR_BYTES:
        raise ValueError(f"{S.MARSH_TAR}: {size} bytes, GEO's file list has {S.MARSH_TAR_BYTES}")
    Xs, obs, ids, syms = [], [], None, None
    with tarfile.open(tar_path) as tar:
        members = {m.name: m for m in tar.getmembers() if m.isfile()}
        want = {k: v for k, v in listed.items() if k != S.MARSH_TAR}
        diff = sorted(set(members) ^ set(want)) + sorted(k for k in set(members) & set(want)
                                                        if members[k].size != want[k])
        if diff:
            raise ValueError(f"the tar's files differ from GEO's file list: {diff[:6]}")
        for gsm, (title, arm, batch) in S.MARSH_SAMPLES.items():
            M, barcodes, fid, fsym = read_10x(tar, members, gsm, title)
            if ids is None:
                ids, syms = fid, fsym
            elif fid != ids:
                raise ValueError(f"{gsm}: its features differ from the first sample's")
            Xs.append(M)
            obs.append(pd.DataFrame(dict(mouse=gsm, arm=arm, batch=batch, title=title, barcode=barcodes)))
            log(f"   {gsm} {title}: {arm}, {batch}; {len(barcodes)} barcodes; {int(M.sum()):,} counts")
    X = sparse.vstack(Xs, format="csr")
    obs = pd.concat(obs, ignore_index=True)
    obs["cell_type"] = S.MARSH_TYPE
    keep = unique_columns(syms)
    genes = [syms[j] for j in keep]
    X = X[:, keep].tocsr()
    missing = [g for g in (*S.MARSH_PAIR, *S.MARSH_SET) if g not in genes]
    f = args.out / "marsh.npz"
    P.save_npz(f, X, obs, genes)
    per = obs.groupby(["arm", "batch", "mouse"]).size().rename("barcodes").reset_index()
    out["marsh"] = dict(tar_bytes=size, tar_sha256=sha256_file(tar_path), cells=int(X.shape[0]), features=len(syms),
                        kept_genes=int(len(keep)), left_out=len(syms) - int(len(keep)), set_missing=missing,
                        barcodes=per.to_dict(orient="records"),
                        post_qc={f"{a}|{b}": v for (a, b), v in S.MARSH_POST_QC.items()},
                        file=f.name, sha256=sha256_file(f), seconds=round(time.time() - t0, 1))
    log(f"   merged: {X.shape[0]} cells x {len(genes)} genes ({len(syms) - len(keep)} features whose symbol names "
        f"several left out); the set's genes not found: {missing or 'none'}")


# --------------------------------------------------------------------------- #
# fabricated GEO files for the dry run
# --------------------------------------------------------------------------- #
def fabricate_geo(geo: Path, seed: int = 5) -> None:
    """GSE152183's SOFT, file list and tar, and GSE132042's and GSE109774's SOFT, all fabricated."""
    from scipy import io as sio
    from scipy import sparse
    geo.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    genes = list(dict.fromkeys([*S.MARSH_SET, "Actb", "P2ry12", "Tmem119", "Cx3cr1", "Hexb", "Dup2", "Dup2",
                                *[f"Gm{i:04d}" for i in range(240)]]))
    genes = genes + ["Dup2"]  # a symbol that names two features
    soft = ["^SERIES = GSE152183", "!Series_title = SYNTHETIC"]
    files = []
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as tar:
        for gsm, (title, arm, batch) in S.MARSH_SAMPLES.items():
            digestion, inhib = S.ARM_CHARACTERISTICS[arm]
            soft += [f"^SAMPLE = {gsm}", f"!Sample_title = {title}", f"!Sample_characteristics_ch1 = batch: {batch}",
                     "!Sample_characteristics_ch1 = brain region: Whole Brain",
                     f"!Sample_characteristics_ch1 = digestion type: {digestion}",
                     f"!Sample_characteristics_ch1 = inhibitors: {inhib}", "!Sample_characteristics_ch1 = sex: M"]
            n = int(rng.integers(150, 260))
            lib = np.exp(rng.normal(0, 0.3, n))
            lam = np.outer(lib, np.exp(rng.normal(-0.5, 1.2, len(genes))))
            if arm == "ENZ-NONE":  # the induction: Fos and Jun up and coupled
                z = rng.normal(0, 1, n)
                for g in S.MARSH_PAIR:
                    lam[:, genes.index(g)] = lib * 4.0 * np.exp(0.9 * z)
            M = sparse.csr_matrix(rng.poisson(lam).T.astype(np.int64))  # genes x cells, as Cell Ranger writes
            mtx = io.BytesIO()
            sio.mmwrite(mtx, M, field="integer")
            parts = {"barcodes.tsv.gz": "\n".join(f"{''.join(rng.choice(list('ACGT'), 16))}-1" for _ in range(n)) + "\n",
                     "features.tsv.gz": "\n".join(f"ENSMUSG{i:011d}\t{g}\tGene Expression"
                                                  for i, g in enumerate(genes)) + "\n"}
            for suffix, content in [*parts.items(), ("matrix.mtx.gz", mtx.getvalue())]:
                data = gzip.compress(content.encode() if isinstance(content, str) else content, mtime=0)
                info = tarfile.TarInfo(f"{gsm}_{title}_{suffix}")
                info.size = len(data)
                tar.addfile(info, io.BytesIO(data))
                files.append((info.name, len(data)))
    (geo / S.MARSH_TAR).write_bytes(buf.getvalue())
    lines = ["#Archive/File\tName\tTime\tSize\tType", f"Archive\t{S.MARSH_TAR}\t\t{len(buf.getvalue())}\tTAR"]
    lines += [f"File\t{n}\t\t{s}\tGZ" for n, s in files]
    (geo / "GSE152183_suppl_filelist.txt").write_text("\n".join(lines) + "\n")
    (geo / "GSE152183_family.soft.gz").write_bytes(gzip.compress(("\n".join(soft) + "\n").encode(), mtime=0))
    tms = ["^SERIES = GSE132042", "!Series_title = SYNTHETIC Tabula Muris Senis",
           "!Series_summary = SYNTHETIC: the 3-month data were published as Tabula Muris.",
           "!Series_overall_design = SYNTHETIC: 1, 3, 18, 21, 24 and 30 months."]
    tm = ["^SERIES = GSE109774", "!Series_title = SYNTHETIC Tabula Muris"]
    for i, m in enumerate(("3_8_M", "3_9_M", "3-M-10", "3_11_M", "3_38_F", "3_39_F", "3_56_F", "3_57_F")):
        tm += [f"^SAMPLE = GSM9{i:06d}", f"!Sample_title = SYNTHETIC plate B00{i} {m}",
               f"!Sample_characteristics_ch1 = mouse id: {m}"]
    for name, text in (("GSE132042", tms), ("GSE109774", tm)):
        (geo / f"{name}_family.soft.gz").write_bytes(gzip.compress(("\n".join(text) + "\n").encode(), mtime=0))


# --------------------------------------------------------------------------- #
def check_census(args) -> int:
    """The dry run's check of the real inputs: the Census opens at the release and the cases' datasets have
    measured genes; GEO's SOFT of GSE152183 matches the registered table. No count is read."""
    spec2 = json.loads((args.data / "select.json").read_text())
    spec1 = json.loads((args.data / "select1.json").read_text())
    census = SyntheticCensus() if args.synthetic else Census()
    datasets = census.datasets()
    out = dict(release="SYNTHETIC" if args.synthetic else S.RELEASE, summary=census.summary(), cases={},
               marsh_soft_problems=check_marsh_soft(args.geo / "GSE152183_family.soft.gz"))
    for name, case in (("case1", spec2["cases"]["case1"]), ("case2", spec2["cases"]["case2"]),
                       ("case3", spec1["cases"]["case3"])):
        ids = sorted(set(map(str, case["data"]["datasets"])))
        present = RS.present_genes(census, ORGANISM[name], datasets, ids)
        out["cases"][name] = dict(datasets=ids, measured_genes=int(len(present)))
    out["source_h5ad_uri"] = (census._cc.get_source_h5ad_uri(S.TMS_DATASET, census_version=S.RELEASE)
                              if not args.synthetic else "SYNTHETIC")
    z = np.load(args.data / "case1.npz", allow_pickle=False)  # the obs of five cells: metadata only
    some = z["obs_soma_joinid"][:5].astype(np.int64)
    rows = census.obs_rows("mus_musculus", some)
    out["census_obs_rows"] = dict(cells=int(len(rows)), columns=list(map(str, rows.columns)),
                                  same_cells=bool(np.array_equal(rows["soma_joinid"].to_numpy(np.int64), some)))
    census.close()
    (args.out / "census_check.json").write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps(out, indent=1, default=str))
    ok = (all(c["measured_genes"] > 0 for c in out["cases"].values()) and not out["marsh_soft_problems"]
          and out["census_obs_rows"]["same_cells"])
    return 0 if ok else 1


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--data", type=Path, required=True, help="the case files, select.json (run 2), select1.json (run 1)")
    p.add_argument("--geo", type=Path, required=True, help="GEO's files: GSE152183's tar, SOFT and file list")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--synthetic", action="store_true")
    p.add_argument("--fabricate-geo", action="store_true", help="write fabricated GEO files into --geo first")
    p.add_argument("--chunk", type=int, default=CHUNK)
    p.add_argument("--only", choices=["census", "marsh"], help="one of the two extractions (the dry run's checks)")
    p.add_argument("--check-census", action="store_true")
    args = p.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    if args.fabricate_geo:
        fabricate_geo(args.geo)
    if args.check_census:
        return check_census(args)
    lines = []

    def log(*msg):
        for m in msg:
            print(m, flush=True)
            lines.append(str(m))
    tag = "SYNTHETIC " if args.synthetic else ""
    log(f"# {tag}The third diagnostics' extraction (validation/exploratory/real_pilot_diagnostics_3_plan.md)", "",
        f"python {platform.python_version()}; numpy {np.__version__}; pandas {pd.__version__}; Census release "
        f"{'SYNTHETIC' if args.synthetic else S.RELEASE}", "")
    out = dict(synthetic=bool(args.synthetic), release="SYNTHETIC" if args.synthetic else S.RELEASE,
               python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, cases={})
    t0, status = time.time(), 0
    try:
        if args.only in (None, "marsh"):
            run_marsh(args, log, out)
        if args.only in (None, "census"):
            run_census(args, log, out)
    except Exception as exc:  # reported as it is; the analyses do not start
        out.update(error=repr(exc), traceback=traceback.format_exc().splitlines())
        log(f"EXCEPTION {exc!r}", *out["traceback"])
        status = 1
    out["seconds"] = round(time.time() - t0, 1)
    log("", f"{tag}extraction: {out['seconds']} s")
    (args.out / "extract.json").write_text(json.dumps(out, indent=1, default=str))
    (args.out / "extract.log").write_text("\n".join(lines) + "\n")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
