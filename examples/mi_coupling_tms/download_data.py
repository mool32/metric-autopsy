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


def _census_pull(which: str, out: Path):
    try:
        import cellxgene_census  # noqa
    except ImportError:
        print("cellxgene-census not installed — `pip install cellxgene-census`, "
              "or download manually:\n")
        print(TMS_MANUAL if which == "tms" else SKIN_MANUAL)
        return False
    import cellxgene_census

    print(f"Opening CELLxGENE Census (stable) to fetch '{which}' → {out} …")
    with cellxgene_census.open_soma(census_version="stable") as census:
        if which == "tms":
            organism = "Mus musculus"
            value_filter = (
                "dataset_id in ['tabula-muris-senis'] or "
                "collection_name == 'Tabula Muris Senis'"
            )
        else:
            organism = "Homo sapiens"
            value_filter = "tissue_general == 'skin' and cell_type == 'fibroblast'"
        adata = cellxgene_census.get_anndata(
            census, organism=organism, obs_value_filter=value_filter
        )
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    print(f"  wrote {adata.n_obs} cells × {adata.n_vars} genes → {out}")
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
