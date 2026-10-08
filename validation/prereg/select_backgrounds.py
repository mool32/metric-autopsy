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
* B4 — Mahdessian et al. 2021: U2OS FUCCI, full-length scRNA-seq with raw counts and FUCCI
  intensities for every cell; >= 300 cells; plate IDs.

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
MISSING_DONOR = {"", "na", "n/a", "nan", "none", "unknown", "not applicable"}
SEED = P.PLAN_SEED


# --------------------------------------------------------------------------- #
# the rules, on an obs table (pure functions: tested on synthetic tables)
# --------------------------------------------------------------------------- #
def _known(donor: pd.Series) -> pd.Series:
    return ~donor.astype(str).str.strip().str.lower().isin(MISSING_DONOR)


def candidates(obs: pd.DataFrame, keys: list[str], min_donors: int, both_sexes: bool = False) -> pd.DataFrame:
    """Every group of `keys` (e.g. dataset and cell type) with >= `min_donors` donors of >= MIN_CELLS
    cells, with the facts that decide the choice: qualifying donors, cells per qualifying donor,
    (for B2) the sexes among them. Sorted by the choice rule: most donors, then most cells per
    donor, then the keys in order."""
    obs = obs[_known(obs["donor_id"])]
    per = obs.groupby(keys + ["donor_id"], observed=True).size().rename("cells").reset_index()
    per = per[per["cells"] >= MIN_CELLS]
    if both_sexes:
        sex = (obs.groupby(keys + ["donor_id"], observed=True)["sex"].agg(lambda s: s.astype(str).mode().iat[0])
               .rename("sex").reset_index())
        per = per.merge(sex, on=keys + ["donor_id"])
    agg = dict(donors=("donor_id", "nunique"), cells=("cells", "sum"))
    if both_sexes:
        agg["sexes"] = ("sex", lambda s: ",".join(sorted(set(map(str, s)))))
    out = per.groupby(keys, observed=True).agg(**agg).reset_index()
    out["cells_per_donor"] = out["cells"] / out["donors"]
    out = out[out["donors"] >= min_donors]
    if both_sexes:
        out = out[out["sexes"].str.contains("female") & out["sexes"].str.split(",").map(lambda x: "male" in x)]
    return out.sort_values(["donors", "cells_per_donor", *keys], ascending=[False, False, *([True] * len(keys))],
                           kind="mergesort").reset_index(drop=True)


def sample_cells(obs: pd.DataFrame, choice: dict, keys: list[str], max_donors: int, max_cells: int,
                 seed: int = SEED) -> np.ndarray:
    """The soma_joinids to extract: the chosen group's qualifying donors (at most `max_donors`, a
    public-seed sample of the sorted IDs) and at most `max_cells` of each one's cells (a
    public-seed sample of its sorted joinids)."""
    rng = np.random.default_rng(seed)
    sel = obs[np.logical_and.reduce([obs[k].astype(str) == str(choice[k]) for k in keys]) & _known(obs["donor_id"])]
    counts = sel.groupby("donor_id", observed=True).size()
    donors = sorted(str(d) for d, n in counts.items() if n >= MIN_CELLS)
    if len(donors) > max_donors:
        donors = sorted(rng.choice(donors, size=max_donors, replace=False).tolist())
    out = []
    for d in donors:
        ids = np.sort(sel.loc[sel["donor_id"].astype(str) == d, "soma_joinid"].to_numpy())
        if len(ids) > max_cells:
            ids = np.sort(rng.choice(ids, size=max_cells, replace=False))
        out.append(ids)
    return np.concatenate(out)


# --------------------------------------------------------------------------- #
# the Census: candidates and extraction
# --------------------------------------------------------------------------- #
def _census_obs(census, organism: str, columns: list[str]) -> pd.DataFrame:
    import cellxgene_census
    flt = "is_primary_data == True and assay in [" + ", ".join(f'"{a}"' for a in DROPLET_3P) + "]"
    return cellxgene_census.get_obs(census, organism, value_filter=flt, column_names=["soma_joinid", *columns])


def extract(census, organism: str, joinids: np.ndarray, path: Path, obs_columns: list[str]) -> dict:
    """The cells' raw counts (all genes, sparse) with gene symbols and their obs, written by
    ``panel.save_npz``; returns the file's sha256 and shape."""
    import cellxgene_census
    ad = cellxgene_census.get_anndata(census, organism, X_name="raw", obs_coords=joinids,
                                      obs_column_names=obs_columns, var_column_names=["feature_id", "feature_name"])
    order = np.argsort(ad.obs["soma_joinid"].to_numpy(), kind="mergesort")
    ad = ad[order]
    obs = ad.obs[obs_columns].astype(str).reset_index(drop=True)
    P.save_npz(path, ad.X.tocsr(), obs, list(ad.var["feature_name"].astype(str)))
    return dict(sha256=P.sha256(path), cells=int(ad.n_obs), genes=int(ad.n_vars))


def _datasets(census) -> pd.DataFrame:
    return census["census_info"]["datasets"].read().concat().to_pandas()


def select_census(out_dir: Path, rehearsal: bool = False) -> tuple[dict, list[str]]:
    """B1 and B2 from the Census stable release. Returns their backgrounds.json entries and the
    report lines (in a rehearsal: only that each step ran, with a hash of the choice)."""
    import cellxgene_census
    version = cellxgene_census.get_census_version_description("stable")
    spec, lines = {}, [f"Census release: stable = {version.get('release_build')} "
                       f"(LTS {version.get('lts', False)})"]
    with cellxgene_census.open_soma(census_version="stable") as census:
        titles = _datasets(census).set_index("dataset_id")["dataset_title"].to_dict()
        rules = (("B1", "Homo sapiens", ["dataset_id", "cell_type"], B1_MIN_DONORS, False, MAX_DONORS),
                 ("B2", "Mus musculus", ["dataset_id", "tissue", "cell_type"], B2_MIN_MICE, True, MAX_DONORS))
        for name, organism, keys, k_min, sexes, max_donors in rules:
            t0 = time.time()
            cols = sorted(set(keys + ["donor_id", "sex", "development_stage", "tissue", "cell_type", "dataset_id"]))
            obs = _census_obs(census, organism, cols)
            cand = candidates(obs, keys, k_min, both_sexes=sexes)
            if cand.empty:
                lines.append(f"{name}: no candidate meets the criteria: dropped with the cases that need it")
                continue
            choice = cand.iloc[0].to_dict()
            digest = hashlib.sha256(json.dumps({k: str(choice[k]) for k in keys}).encode()).hexdigest()
            if rehearsal:
                lines.append(f"{name}: {len(cand)} candidates; choice sha256 {digest[:16]}; {time.time() - t0:.0f} s")
                continue
            lines += ["", f"## {name} ({organism})", "", "| rank | " + " | ".join(keys) + " | title | donors >= "
                      f"{MIN_CELLS} cells | cells per donor |" + (" sexes |" if sexes else ""),
                      "|---" * (len(keys) + 4 + int(sexes)) + "|"]
            for i, row in cand.head(15).iterrows():
                lines.append(f"| {i + 1} | " + " | ".join(str(row[k]) for k in keys)
                             + f" | {titles.get(row['dataset_id'], '?')} | {row['donors']} | {row['cells_per_donor']:.0f} |"
                             + (f" {row['sexes']} |" if sexes else ""))
            ids = sample_cells(obs, choice, keys, max_donors, MAX_CELLS)
            rec = extract(census, organism, ids, out_dir / f"{name}.npz", cols)
            spec[name] = dict(file=f"{name}.npz", sha256=rec["sha256"], donor="donor_id", counts="X",
                              source=dict(census_release=version.get("release_build"), organism=organism,
                                          **{k: str(choice[k]) for k in keys},
                                          title=titles.get(choice["dataset_id"]), qualifying_donors=int(choice["donors"]),
                                          extracted_cells=rec["cells"], genes=rec["genes"],
                                          rule=f"max {max_donors} donors, {MAX_CELLS} cells each, seed {SEED}"))
            lines.append(f"\nChosen: rank 1; extracted {rec['cells']} cells x {rec['genes']} genes, sha256 {rec['sha256']}")
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


def parse_buettner(tables: dict) -> tuple[np.ndarray, list, pd.DataFrame]:
    """Genes x cells tables (first columns: gene ID, symbol and optionally length; the rest cells)
    into cells x genes counts with each cell's phase."""
    blocks, obs, genes = [], [], None
    for phase, text in tables.items():
        df = pd.read_csv(io.StringIO(text), sep="\t")
        meta = [c for c in df.columns[:3] if not np.issubdtype(df[c].dtype, np.number) or "length" in c.lower()]
        sym = next((c for c in meta if "name" in c.lower() or "symbol" in c.lower()), meta[0])
        g = [str(s) if isinstance(s, str) and s else str(i) for s, i in zip(df[sym], df[meta[0]])]
        if genes is None:
            genes = g
        elif g != genes:
            raise ValueError("the phase tables list different genes")
        counts = df.drop(columns=meta).to_numpy(dtype=float).T
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
        return None, [f"B3: not obtained ({exc}); dropped with R2"]
    ercc = sum(g.upper().startswith("ERCC") for g in genes)
    per_phase = obs["phase"].value_counts().to_dict()
    raw = bool(np.all(X >= 0) and np.allclose(X, np.round(X)))
    fails = ([] if ercc else ["no ERCC rows"]) + ([] if min(per_phase.values()) >= 50 else [f"cells per phase {per_phase}"]) \
        + ([] if raw else ["not raw counts"])
    if fails:
        return None, [f"B3: fails the criteria ({'; '.join(fails)}); dropped with R2"]
    P.save_npz(out_dir / "B3.npz", np.round(X).astype(np.int64), obs, genes)
    spec = dict(file="B3.npz", sha256=P.sha256(out_dir / "B3.npz"), donor="batch", counts="X",
                source=dict(study="E-MTAB-2805", urls=used, cells_per_phase=per_phase, ercc_genes=ercc))
    return spec, [f"B3: E-MTAB-2805, cells per phase {per_phase}, {ercc} ERCC genes; sha256 {spec['sha256']}"]


B4_SOURCES: dict = {}  # Mahdessian et al. 2021: no machine-readable source with FUCCI intensities per cell is fixed
B4_LISTINGS = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146773/suppl/",
               "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146773&targ=self&form=text&view=brief")


def fetch_b4(out_dir: Path) -> tuple[dict | None, list[str]]:
    if not B4_SOURCES:
        return None, ["B4: no source with raw counts and per-cell FUCCI intensities is fixed for the named "
                      "candidate (Mahdessian et al. 2021); dropped with R3"]
    return None, ["B4: dropped with R3"]


def list_b4_sources() -> list[str]:
    """The public listings of the named B4 candidate's deposit (printed by the rehearsal, so that
    a source can be fixed before the tag)."""
    import re
    out = []
    for url in B4_LISTINGS:
        try:
            text = _download([url], timeout=60)[0].decode(errors="replace")
        except Exception as exc:
            out.append(f"B4 listing {url}: {exc}")
            continue
        names = sorted(set(re.findall(r'href="([^"?/][^"]*)"', text))) or text.splitlines()[:40]
        out.append(f"B4 listing {url}:")
        out += [f"  {n}" for n in names[:60]]
    return out


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
        print("\n".join(lines + list_b4_sources()))
        return
    args.spec.write_text(json.dumps(spec, indent=1, sort_keys=True))
    args.report.write_text("# Background selection (validation/prereg/v1.md, section 3.1)\n\n" + "\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
