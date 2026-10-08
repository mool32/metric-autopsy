"""Mechanical selection of the panel's backgrounds (validation/prereg/v1.md, section 3.1).

Run by the validation workflow after the tag ``v0.3.0-prereg``, in GitHub Actions (the agent's
container cannot reach the data hosts). B1 and B2 come from the CELLxGENE Census (its ``stable``
release at the time, recorded): every candidate that meets every criterion is listed with the
facts that decide it, the choice rule picks one, and the chosen cells are extracted with their
raw counts. B3 and B4 are the named candidates, downloaded from their public sources and checked
against the criteria. A background without a qualifying candidate is dropped with the cases that
need it, and the drop is recorded. Nothing is substituted by judgment.

Rules (section 3.1):

* B1 — human, droplet 3' (10x 3' v1/v2/v3), primary cells; one cell type of one dataset with >= 200
  cells in each of >= 24 donors (donor IDs present). Choice: the most such donors; ties: more cells
  per qualifying donor, then the dataset and cell type IDs in order. At most 200 of the qualifying
  donors (a sample with the public seed, when there are more) and at most 400 cells per donor (a
  sample with the public seed) are extracted.
* B2 — mouse, droplet 3', primary cells; one tissue and cell type of one dataset with >= 12 mice of
  >= 200 cells each, both sexes among them. Choice: the most such mice; ties: more cells per mouse,
  then the IDs in order. All qualifying mice, at most 400 cells each.
* B3 — Buettner et al. 2015 (E-MTAB-2805): plate-based, ERCC spike-ins, phase from a DNA-content
  sort; >= 50 cells per phase; raw counts including ERCC.
* B4 — Mahdessian et al. 2021 (GEO GSE146773): U2OS FUCCI, full-length scRNA-seq with raw counts
  and FUCCI intensities for every cell; >= 300 cells; plate IDs. The deposit's counts are RSEM
  expected counts (not integers), so B4 fails the raw-counts criterion unless the deposit changes.

Raw counts (non-negative integers) are a criterion of B1 and B2 that only the extracted matrix
shows: it is checked on the extraction, in the order of the choice rule, over at most the first
MAX_RANKS candidates; the first that has them is chosen.

    python select_backgrounds.py --out-dir DATA --spec backgrounds.json --report selection.md
    python select_backgrounds.py --rehearsal      # before the tag: the code path, no identities printed
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

import panel as P

DROPLET_3P = ("10x 3' v1", "10x 3' v2", "10x 3' v3", "10x 3' transcription profiling")
MIN_CELLS = P.CELLS_PER_DONOR          # cells per donor (mouse)
B1_MIN_DONORS, B2_MIN_MICE = 24, 12
MAX_DONORS, MAX_CELLS = 200, 400       # extraction caps (public-seed samples)
MAX_RANKS = 5                          # candidates tried, in order, for the raw-counts criterion
MISSING_DONOR = {"", "na", "n/a", "nan", "none", "unknown", "not applicable"}
SEED = P.PLAN_SEED


# --------------------------------------------------------------------------- #
# the rules, on cell counts (pure functions: tested on synthetic tables)
# --------------------------------------------------------------------------- #
def _codes_mask(s: pd.Series, hit_of_categories) -> np.ndarray:
    """A boolean per row of a categorical column, computed once per category: the 75 million
    human cells never become strings one by one (a missing value is the string "nan")."""
    cats = pd.Index(s.cat.categories.astype(str))
    hit = np.asarray(hit_of_categories(cats), dtype=bool)
    codes = s.cat.codes.to_numpy()
    missing = bool(np.asarray(hit_of_categories(pd.Index(["nan"])))[0])
    return np.where(codes >= 0, hit[np.maximum(codes, 0)], missing)


def _known(donor: pd.Series) -> np.ndarray:
    """True where the donor ID is present (not empty, missing or a placeholder)."""
    def ok(v: pd.Index):
        return ~v.str.strip().str.lower().isin(MISSING_DONOR)
    if isinstance(donor.dtype, pd.CategoricalDtype):
        return _codes_mask(donor, ok)
    return np.asarray(ok(pd.Index(donor.astype(str))))


def _equals(s: pd.Series, value) -> np.ndarray:
    if isinstance(s.dtype, pd.CategoricalDtype):
        return _codes_mask(s, lambda v: v == str(value))
    return (s.astype(str) == str(value)).to_numpy()


def count_cells(batches, cols: list[str]) -> pd.DataFrame:
    """Cells per combination of `cols`, summed over `batches` (Arrow tables or pandas frames read one
    at a time, e.g. the Census obs in chunks), with the values as strings (missing: "nan"). Each
    batch is grouped on its category codes; only the combinations become strings."""
    cols = list(cols)
    parts = []
    for b in batches:
        df = b.to_pandas() if hasattr(b, "to_pandas") else b
        n = df[cols].groupby(cols, observed=True, dropna=False).size()
        part = n.rename("cells").reset_index()
        for c in cols:
            part[c] = part[c].astype(str)
        parts.append(part)
    if not parts:
        return pd.DataFrame({**{c: pd.Series(dtype=str) for c in cols}, "cells": pd.Series(dtype=int)})
    return pd.concat(parts, ignore_index=True).groupby(cols, sort=True)["cells"].sum().reset_index()


def candidates(counts: pd.DataFrame, keys: list[str], min_donors: int, both_sexes: bool = False) -> pd.DataFrame:
    """Every group of `keys` (e.g. dataset and cell type) with >= `min_donors` donors of >= MIN_CELLS
    cells, with the facts that decide the choice: qualifying donors, cells per qualifying donor,
    (for B2) the sexes among them. `counts` holds the cells per (keys, donor_id[, sex]) from
    `count_cells`; a donor's sex is the commonest among its cells (ties: the first in order).
    Sorted by the choice rule: most donors, then most cells per donor, then the keys in order."""
    counts = counts[_known(counts["donor_id"])]
    g = list(keys) + ["donor_id"]
    per = counts.groupby(g, sort=True)["cells"].sum().reset_index()
    per = per[per["cells"] >= MIN_CELLS]
    if both_sexes:
        top = counts.sort_values(g + ["cells", "sex"], ascending=[True] * len(g) + [False, True], kind="mergesort")
        per = per.merge(top.drop_duplicates(g)[g + ["sex"]], on=g)
    agg = dict(donors=("donor_id", "nunique"), cells=("cells", "sum"))
    if both_sexes:
        agg["sexes"] = ("sex", lambda s: ",".join(sorted(set(map(str, s)))))
    out = per.groupby(keys, observed=True).agg(**agg).reset_index()
    out["cells_per_donor"] = out["cells"] / out["donors"]
    out = out[out["donors"] >= min_donors]
    if both_sexes:
        out = out[out["sexes"].str.split(",").map(lambda x: "female" in x and "male" in x)]
    return out.sort_values(["donors", "cells_per_donor", *keys], ascending=[False, False, *([True] * len(keys))],
                           kind="mergesort").reset_index(drop=True)


def sample_cells(obs: pd.DataFrame, choice: dict, keys: list[str], max_donors: int | None, max_cells: int,
                 seed: int = SEED) -> np.ndarray:
    """The soma_joinids to extract: the chosen group's qualifying donors (at most `max_donors`, a
    public-seed sample of the sorted IDs; None: all) and at most `max_cells` of each one's cells (a
    public-seed sample of its sorted joinids)."""
    rng = np.random.default_rng(seed)
    mask = np.logical_and.reduce([_equals(obs[k], choice[k]) for k in keys]) & _known(obs["donor_id"])
    sel = pd.DataFrame({"soma_joinid": obs["soma_joinid"].to_numpy()[mask],
                        "donor_id": obs["donor_id"][mask].astype(str).to_numpy()})
    ids_of = {d: np.sort(v.to_numpy()) for d, v in sel.groupby("donor_id")["soma_joinid"]}
    donors = sorted(d for d, ids in ids_of.items() if len(ids) >= MIN_CELLS)
    if max_donors is not None and len(donors) > max_donors:
        donors = sorted(rng.choice(donors, size=max_donors, replace=False).tolist())
    out = []
    for d in donors:
        ids = ids_of[d]
        if len(ids) > max_cells:
            ids = np.sort(rng.choice(ids, size=max_cells, replace=False))
        out.append(ids)
    return np.concatenate(out)


# --------------------------------------------------------------------------- #
# the Census: candidates and extraction
# --------------------------------------------------------------------------- #
CENSUS_FILTER = "is_primary_data == True and assay in [" + ", ".join(f'"{a}"' for a in DROPLET_3P) + "]"


def _obs(census, organism: str):
    return census["census_data"][organism.lower().replace(" ", "_")].obs


def _census_counts(census, organism: str, cols: list[str]) -> pd.DataFrame:
    """Cells per (cols) over the organism's primary droplet-3' cells, read in chunks."""
    return count_cells(_obs(census, organism).read(value_filter=CENSUS_FILTER, column_names=list(cols)), cols)


def _census_group(census, organism: str, choice: dict, keys: list[str]) -> pd.DataFrame:
    """The chosen dataset's primary droplet-3' cells (joinid, keys, donor) for the extraction sample."""
    flt = CENSUS_FILTER + f' and dataset_id == "{choice["dataset_id"]}"'
    return (_obs(census, organism).read(value_filter=flt, column_names=["soma_joinid", *keys, "donor_id"])
            .concat().to_pandas())


def raw_counts(X) -> bool:
    """The "raw counts" criterion on an extracted matrix: non-negative integers."""
    v = X.data if hasattr(X, "data") and not isinstance(X, np.ndarray) else np.asarray(X)
    return bool(np.all(v >= 0) and np.all(v == np.round(v)))


def extract(census, organism: str, joinids: np.ndarray, path: Path, obs_columns: list[str]) -> dict:
    """The cells' raw counts (all genes, sparse) with gene symbols and their obs, written by
    ``panel.save_npz`` when they are raw counts; returns the file's sha256 and shape, and whether
    the counts are raw (non-negative integers: else nothing is written)."""
    import cellxgene_census
    ad = cellxgene_census.get_anndata(census, organism, X_name="raw", obs_coords=joinids,
                                      obs_column_names=["soma_joinid", *obs_columns],
                                      var_column_names=["feature_id", "feature_name"])
    order = np.argsort(ad.obs["soma_joinid"].to_numpy(), kind="mergesort")
    ad = ad[order]
    X = ad.X.tocsr()
    if not raw_counts(X):
        return dict(raw=False, cells=int(ad.n_obs), genes=int(ad.n_vars))
    obs = ad.obs[obs_columns].astype(str).reset_index(drop=True)
    P.save_npz(path, X, obs, list(ad.var["feature_name"].astype(str)))
    return dict(raw=True, sha256=P.sha256(path), cells=int(ad.n_obs), genes=int(ad.n_vars))


def _datasets(census) -> pd.DataFrame:
    return census["census_info"]["datasets"].read().concat().to_pandas()


def select_census(out_dir: Path, rehearsal: bool = False) -> tuple[dict, list[str]]:
    """B1 and B2 from the Census stable release. Returns their backgrounds.json entries and the
    report lines (in a rehearsal: only that each step ran, with a hash of the choice, and a small
    extraction to exercise the code path; nothing is kept)."""
    import cellxgene_census
    version = cellxgene_census.get_census_version_description("stable")
    release = version.get("release_build") or "stable"
    spec, lines = {}, [f"Census release: stable = {release} (LTS {version.get('lts', False)})"]
    with cellxgene_census.open_soma(census_version=release) as census:
        titles = _datasets(census).set_index("dataset_id")["dataset_title"].to_dict()
        rules = (("B1", "Homo sapiens", ["dataset_id", "cell_type"], B1_MIN_DONORS, False, MAX_DONORS),
                 ("B2", "Mus musculus", ["dataset_id", "tissue", "cell_type"], B2_MIN_MICE, True, None))
        for name, organism, keys, k_min, sexes, max_donors in rules:
            t0 = time.time()
            counted = keys + ["donor_id"] + (["sex"] if sexes else [])
            counts = _census_counts(census, organism, counted)
            lines.append(f"{name}: {int(counts['cells'].sum()):,} primary droplet-3' cells of {organism} counted "
                         f"in {time.time() - t0:.0f} s")
            cand = candidates(counts, keys, k_min, both_sexes=sexes)
            if cand.empty:
                lines.append(f"{name}: no candidate meets the criteria: dropped with the cases that need it")
                continue
            cols = sorted(set(keys + ["donor_id", "sex", "development_stage", "tissue", "cell_type", "dataset_id"]))
            if not rehearsal:
                lines += ["", f"## {name} ({organism})", "", "| rank | " + " | ".join(keys) + " | title | donors >= "
                          f"{MIN_CELLS} cells | cells per donor |" + (" sexes |" if sexes else ""),
                          "|---" * (len(keys) + 4 + int(sexes)) + "|"]
                for i, row in cand.head(15).iterrows():
                    lines.append(f"| {i + 1} | " + " | ".join(str(row[k]) for k in keys)
                                 + f" | {titles.get(row['dataset_id'], '?')} | {row['donors']} | {row['cells_per_donor']:.0f} |"
                                 + (f" {row['sexes']} |" if sexes else ""))
            # the raw-counts criterion is checked on the extraction, in the order of the choice rule
            for rank in range(1, min(len(cand), MAX_RANKS) + 1):
                choice = cand.iloc[rank - 1].to_dict()
                group = _census_group(census, organism, choice, keys)
                if rehearsal:
                    ids = sample_cells(group, choice, keys, max_donors=2, max_cells=50)
                    target = out_dir / f"{name}-rehearsal.npz"
                else:
                    ids = sample_cells(group, choice, keys, max_donors, MAX_CELLS)
                    target = out_dir / f"{name}.npz"
                rec = extract(census, organism, ids, target, cols)
                if rec["raw"]:
                    break
                lines.append(f"{name}: rank {rank} fails the raw-counts criterion (values that are not "
                             "non-negative integers)" + ("" if rehearsal else f": {choice['dataset_id']}"))
            else:
                lines.append(f"{name}: none of the first {MAX_RANKS} candidates has raw counts: dropped with the "
                             "cases that need it")
                continue
            if rehearsal:
                target.unlink()
                digest = hashlib.sha256(json.dumps({k: str(choice[k]) for k in keys}).encode()).hexdigest()
                lines.append(f"{name}: {len(cand)} candidates; chosen rank {rank}, sha256 {digest[:16]}, with "
                             f"{int(choice['donors'])} qualifying donors; a {rec['cells']}-cell test extraction of "
                             f"{rec['genes']} genes with raw counts; {time.time() - t0:.0f} s")
                continue
            spec[name] = dict(file=f"{name}.npz", sha256=rec["sha256"], donor="donor_id", counts="X",
                              source=dict(census_release=release, organism=organism, rank=rank,
                                          **{k: str(choice[k]) for k in keys},
                                          title=titles.get(choice["dataset_id"]), qualifying_donors=int(choice["donors"]),
                                          extracted_cells=rec["cells"], genes=rec["genes"],
                                          rule=f"{'all' if max_donors is None else f'max {max_donors}'} donors, "
                                               f"{MAX_CELLS} cells each, seed {SEED}"))
            lines.append(f"\nChosen: rank {rank}; extracted {rec['cells']} cells x {rec['genes']} genes, "
                         f"sha256 {rec['sha256']}")
    return spec, lines


# --------------------------------------------------------------------------- #
# B3 and B4: the named candidates from their public sources
# --------------------------------------------------------------------------- #
B3_SOURCES = {  # E-MTAB-2805: one counts table per sorted phase (genes x cells, ERCC rows included)
    phase: [f"https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/805/E-MTAB-2805/Files/{phase}_singlecells_counts.txt",
            f"https://www.ebi.ac.uk/biostudies/files/E-MTAB-2805/{phase}_singlecells_counts.txt",
            f"https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2805/{phase}_singlecells_counts.txt"]
    for phase in ("G1", "S", "G2M")}


def _download(urls: list[str], timeout: float = 120.0) -> tuple[bytes, str]:
    last = None
    for url in urls:
        for attempt in range(3):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "metric-autopsy-panel"})
                with urllib.request.urlopen(req, timeout=timeout) as fh:
                    data = fh.read()
                if url.endswith(".gz"):
                    data = gzip.decompress(data)
                return data, url
            except Exception as exc:  # recorded if every source fails
                last = f"{url}: {exc!r}"
                time.sleep(2 ** attempt)
    raise RuntimeError(last)


def describe_table(text: str, name: str, sep: str = "\t") -> list[str]:
    """What a downloaded table looks like (for the rehearsal: the format of a named source is
    checked before the tag): its first lines, its shape, and per column the dtype and, for numeric
    columns, the missing, negative and non-integer values."""
    head = [ln[:400] for ln in text.splitlines()[:3]]
    out = [f"{name}: first lines:"] + [f"    {ln}" for ln in head]
    try:
        df = pd.read_csv(io.StringIO(text), sep=sep)
    except Exception as exc:
        return out + [f"{name}: not parsed ({exc!r})"]
    out.append(f"{name}: {df.shape[0]} rows x {df.shape[1]} columns")
    num = df.select_dtypes("number")
    v = num.to_numpy(dtype=float) if num.shape[1] else np.zeros((0, 0))
    out.append(f"{name}: {num.shape[1]} numeric columns; missing {int(np.isnan(v).sum())}, negative "
               f"{int((v < 0).sum())}, non-integer {int((np.abs(v - np.round(v)) > 1e-9).sum())}")
    for c in list(df.columns[:6]) + ([df.columns[-1]] if df.shape[1] > 6 else []):
        col = df[c]
        if np.issubdtype(col.dtype, np.number):
            x = col.to_numpy(dtype=float)
            out.append(f"    column {c!r}: {col.dtype}, missing {int(np.isnan(x).sum())}, non-integer "
                       f"{int((np.abs(x - np.round(x)) > 1e-9).sum())}, range {np.nanmin(x):g}..{np.nanmax(x):g}")
        else:
            out.append(f"    column {c!r}: {col.dtype}, e.g. {list(map(str, col.dropna().unique()[:3]))}")
    return out


def parse_buettner(tables: dict) -> tuple[np.ndarray, list, pd.DataFrame]:
    """Genes x cells tables (first columns: gene ID, symbol and optionally length; the rest cells)
    into cells x genes counts with each cell's phase."""
    blocks, obs, genes = [], [], None
    for phase, text in tables.items():
        df = pd.read_csv(io.StringIO(text), sep="\t")
        # the cells are the columns named as cells; the leading others describe the genes (IDs,
        # symbol, length)
        cells = [c for c in df.columns if "cell" in str(c).lower()]
        meta = [c for c in df.columns if c not in cells]
        if not cells or not meta:
            raise ValueError(f"{phase}: no cell or no gene columns among {list(df.columns[:6])}")
        sym = next((c for c in meta if "name" in c.lower() or "symbol" in c.lower()), meta[0])
        g = [str(s) if isinstance(s, str) and s else str(i) for s, i in zip(df[sym], df[meta[0]])]
        if genes is None:
            genes = g
        elif g != genes:
            raise ValueError("the phase tables list different genes")
        counts = df[cells].to_numpy(dtype=float).T
        blocks.append(counts)
        obs += [dict(phase=phase, batch=phase) for _ in range(counts.shape[0])]
    return np.vstack(blocks), genes, pd.DataFrame(obs)


def fetch_b3(out_dir: Path) -> tuple[dict | None, list[str]]:
    try:
        texts, used = {}, []
        for phase, urls in B3_SOURCES.items():
            data, url = _download(urls)
            texts[phase] = data.decode()
            used.append(url)
        X, genes, obs = parse_buettner(texts)
    except Exception as exc:
        diag = [ln for ph, t in texts.items() for ln in describe_table(t, f"B3 {ph}")] if texts else []
        return None, [f"B3: not obtained ({exc}); dropped with R2"] + diag
    ercc = sum(g.upper().startswith("ERCC") for g in genes)
    per_phase = obs["phase"].value_counts().to_dict()
    raw = bool(np.all(X >= 0) and np.allclose(X, np.round(X)))
    fails = ([] if ercc else ["no ERCC rows"]) + ([] if min(per_phase.values()) >= 50 else [f"cells per phase {per_phase}"]) \
        + ([] if raw else ["not raw counts"])
    if fails:
        return None, ([f"B3: fails the criteria ({'; '.join(fails)}); dropped with R2"]
                      + [ln for ph, t in texts.items() for ln in describe_table(t, f"B3 {ph}")])
    P.save_npz(out_dir / "B3.npz", np.round(X).astype(np.int64), obs, genes)
    spec = dict(file="B3.npz", sha256=P.sha256(out_dir / "B3.npz"), donor="batch", counts="X",
                source=dict(study="E-MTAB-2805", urls=used, cells_per_phase=per_phase, ercc_genes=ercc))
    return spec, [f"B3: E-MTAB-2805, cells per phase {per_phase}, {ercc} ERCC genes; sha256 {spec['sha256']}"]


B4_DEPOSIT = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146773/suppl/"  # Mahdessian et al. 2021
B4_FILES = dict(fucci="GSE146773_fucci_coords.csv.gz", counts="GSE146773_Counts.csv.gz")
B4_MIN_CELLS = 300


def check_b4(counts: pd.DataFrame, fucci: pd.DataFrame) -> tuple[list[str], dict]:
    """The B4 criteria on the deposit's tables (cells x genes counts; one FUCCI row per cell, named
    ``<well>_<plate>``): raw counts, FUCCI intensities for every cell, >= 300 cells, plate IDs."""
    cells = pd.Index(counts.index.astype(str))
    X = counts.to_numpy(dtype=float)
    nonint, neg = int((np.abs(X - np.round(X)) > 1e-9).sum()), int((X < 0).sum())
    f = fucci.assign(cell=fucci["cell"].astype(str)).set_index("cell")
    measured = f[["raw_green530", "raw_red585"]].notna().all(axis=1)
    with_fucci = int(cells.isin(f.index[measured]).sum())
    plates = cells.str.extract(r"_(\d+)$")[0]
    fails = []
    if nonint or neg:
        fails.append(f"not raw counts ({nonint:,} non-integer and {neg:,} negative values)")
    if with_fucci < len(cells):
        fails.append(f"FUCCI intensities for {with_fucci} of {len(cells)} cells")
    if len(cells) < B4_MIN_CELLS:
        fails.append(f"{len(cells)} cells")
    if plates.isna().any():
        fails.append(f"no plate ID for {int(plates.isna().sum())} cells")
    facts = dict(cells=len(cells), genes=int(X.shape[1]), non_integer=nonint, with_fucci=with_fucci,
                 plates=sorted(set(plates.dropna())))
    return fails, facts


def fetch_b4(out_dir: Path) -> tuple[dict | None, list[str]]:
    """The named B4 candidate's deposit, checked against the criteria. Its counts are RSEM expected
    counts in the rehearsal, so B4 and R3 are dropped unless the deposit changes; R3 has no runner."""
    try:
        tables = {k: pd.read_csv(io.StringIO(_download([B4_DEPOSIT + f], timeout=300)[0].decode()))
                  for k, f in B4_FILES.items()}
        fails, facts = check_b4(tables["counts"].set_index(tables["counts"].columns[0]), tables["fucci"])
    except Exception as exc:
        return None, [f"B4: not obtained ({exc!r}); dropped with R3"]
    where = f"GEO GSE146773 ({', '.join(B4_FILES.values())}: {facts['cells']} cells, {len(facts['plates'])} plates)"
    if fails:
        return None, [f"B4: {where} fails the criteria ({'; '.join(fails)}); dropped with R3"]
    return None, [f"B4: {where} meets the criteria, but no R3 runner was pre-registered (its counts were not "
                  "raw in the rehearsal); dropped with R3"]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--out-dir", type=Path, default=Path("data"))
    p.add_argument("--spec", type=Path, default=Path("backgrounds.json"))
    p.add_argument("--report", type=Path, default=Path("selection.md"))
    p.add_argument("--rehearsal", action="store_true",
                   help="run the Census rules and the B3/B4 downloads, print no identities, write nothing")
    args = p.parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    spec, lines = select_census(args.out_dir, args.rehearsal)
    for fetch in (fetch_b3, fetch_b4):
        s, more = fetch(args.out_dir)
        lines += more
        if s is not None and not args.rehearsal:
            spec[more[0].split(":")[0]] = s
    if args.rehearsal:
        print("\n".join(lines))
        return
    args.spec.write_text(json.dumps(spec, indent=1, sort_keys=True))
    args.report.write_text("# Background selection (validation/prereg/v1.md, section 3.1)\n\n" + "\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
