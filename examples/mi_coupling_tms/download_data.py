#!/usr/bin/env python3
"""Fetch the datasets for the MI-coupling worked example.

Both datasets are public and pulled from the CELLxGENE Census when `cellxgene-census` is
installed; otherwise the script prints the manual portal instructions rather than pretending
to succeed. Nothing here is committed — see ../../DATASETS.md.

    pip install cellxgene-census
    python download_data.py --which all          # -> data/raw/*.h5ad
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

RAW = Path(__file__).resolve().parents[2] / "data" / "raw"

TMS_MANUAL = """\
Tabula Muris Senis (FACS / SmartSeq2), mouse, young 3mo vs old 24mo:
  - CELLxGENE portal: https://cellxgene.cziscience.com/  (search "Tabula Muris Senis")
  - or figshare collection 'Tabula Muris Senis' (per-tissue .h5ad)
  Save the FACS object to data/raw/tms_facs.h5ad with obs columns: age, sex, tissue, cell_type.
"""

SKIN_MANUAL = """\
Human skin fibroblasts (10x), 179 donors, ages 18-79:
  - CELLxGENE portal: https://cellxgene.cziscience.com/  (search "skin fibroblast")
  Save to data/raw/human_skin.h5ad with obs columns: age (or age group), sex, cell_type.
"""


def _resolve_tms_dataset_ids(census) -> list[str]:
    """Resolve the Tabula Muris Senis FACS (Smart-seq2) dataset_id(s) from the Census
    datasets table.

    dataset_id is a UUID, so it must be looked up by name — the previous filter used a
    hardcoded slug ('tabula-muris-senis') and 'collection_name' as if it were a per-cell obs
    column; neither exists in the obs schema, so the query returned nothing.
    """
    ds = census["census_info"]["datasets"].read().concat().to_pandas()
    is_tms = ds["collection_name"].str.contains("Tabula Muris Senis", case=False, na=False)
    is_facs = ds["dataset_title"].str.contains("facs|smart-seq", case=False, na=False)
    hit = ds[is_tms & is_facs]
    if hit.empty:                      # fall back to any TMS dataset; the obs assay filter narrows it
        hit = ds[is_tms]
    return hit["dataset_id"].tolist()


def _make_metric_autopsy_ready(adata):
    """Rename Census QC columns to what `metric_autopsy.qc.per_cell_qc` expects and derive a
    simple `age` token, so the downloaded object runs through the gates and the documented
    `--groups 3m 24m` CLI with no extra munging.
    """
    obs = adata.obs
    if "nnz" in obs:
        obs["n_genes_by_counts"] = obs["nnz"].astype(float)   # per-cell detected genes
    if "raw_sum" in obs:
        obs["total_counts"] = obs["raw_sum"].astype(float)
    if "development_stage" in obs and "age" not in obs:
        def _age_token(s):
            s = str(s)
            m = re.search(r"(\d+)\s*-?\s*(month|year)", s)
            return f"{m.group(1)}{'m' if m.group(2) == 'month' else 'y'}" if m else s
        obs["age"] = obs["development_stage"].map(_age_token)
    return adata


def _census_pull(which: str, out: Path):
    try:
        import cellxgene_census  # noqa
    except ImportError:
        print("cellxgene-census not installed — `pip install cellxgene-census`, "
              "or download manually:\n")
        print(TMS_MANUAL if which == "tms" else SKIN_MANUAL)
        return False
    import cellxgene_census

    # Only the columns the worked example needs (keeps the pull lighter). QC columns
    # (nnz, raw_sum) let the gates skip recomputing n_genes from the full matrix.
    obs_cols = ["assay", "dataset_id", "cell_type", "tissue", "tissue_general",
                "sex", "development_stage", "donor_id", "is_primary_data", "nnz", "raw_sum"]
    print(f"Opening CELLxGENE Census (stable) to fetch '{which}' → {out} …")
    with cellxgene_census.open_soma(census_version="stable") as census:
        if which == "tms":
            organism = "Mus musculus"
            ids = _resolve_tms_dataset_ids(census)
            if not ids:
                print("Could not locate Tabula Muris Senis in the Census datasets table.\n")
                print(TMS_MANUAL)
                return False
            id_list = ", ".join(f"'{i}'" for i in ids)
            value_filter = (f"dataset_id in [{id_list}] and assay == 'Smart-seq2' "
                            "and is_primary_data == True")
        else:
            organism = "Homo sapiens"
            value_filter = ("tissue_general == 'skin' and cell_type == 'fibroblast' "
                            "and is_primary_data == True")
        adata = cellxgene_census.get_anndata(
            census, organism=organism, obs_value_filter=value_filter,
            column_names={"obs": obs_cols},
        )

    adata = _make_metric_autopsy_ready(adata)
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    print(f"  wrote {adata.n_obs} cells × {adata.n_vars} genes → {out}")
    age_levels = sorted({str(a) for a in adata.obs.get("age", [])})
    print(f"  age levels present (pick two for --groups): {age_levels[:10]}")
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--which", choices=["tms", "skin", "all"], default="all")
    args = ap.parse_args(argv)

    targets = {
        "tms": RAW / "tms_facs.h5ad",
        "skin": RAW / "human_skin.h5ad",
    }
    picks = ["tms", "skin"] if args.which == "all" else [args.which]
    ok = all(_census_pull(w, targets[w]) for w in picks)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
