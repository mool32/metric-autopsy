"""The Census metadata of validation v1's backgrounds B1 and B2: the suspension type (cell or
nucleus) and the assay of the chosen group's cells, read with the selection's own filter. v1's
selection rule (v1.md section 3.1) did not use them, and the extraction did not store them; the
README states them from this script's output.

    python validation/exploratory/census_metadata.py   # needs the Census (requirements-select.txt)
"""
from __future__ import annotations

import sys
from pathlib import Path

import cellxgene_census

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "prereg"))
from select_backgrounds import CENSUS_FILTER  # noqa: E402  (the selection's filter: primary droplet-3' cells)

RELEASE = "2025-11-08"  # the release of the selection (selection.md)
GROUPS = {  # the chosen groups of selection.md / backgrounds.json
    "B1": ("homo_sapiens", 'dataset_id == "37a17b78-4864-4a42-b67b-31c00962795a" and cell_type == "oligodendrocyte"'),
    "B2": ("mus_musculus", 'dataset_id == "49e4ffcc-5444-406d-bdee-577127404ba8" and tissue == "islet of Langerhans"'
                           ' and cell_type == "type B pancreatic cell"'),
}


def main() -> int:
    with cellxgene_census.open_soma(census_version=RELEASE) as census:
        for name, (organism, flt) in GROUPS.items():
            obs = (census["census_data"][organism].obs
                   .read(value_filter=f"{CENSUS_FILTER} and {flt}", column_names=["suspension_type", "assay"])
                   .concat().to_pandas())
            cols = ["suspension_type", "assay"]  # the read also returns soma_joinid
            counts = {f"{s} / {a}": int(n) for (s, a), n in obs[cols].value_counts().sort_index().items()}
            print(f"{name} ({organism}): {len(obs)} primary droplet-3' cells of the chosen group; "
                  f"suspension type / assay: {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
