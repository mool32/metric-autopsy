"""The Census metadata of validation v1's backgrounds B1 and B2: the suspension type (cell or
nucleus) and the assay of the chosen group's cells, read with the selection's own filter. v1's
selection rule (v1.md section 3.1) did not use them, and the extraction did not store them; the
README states them from this script's output.

For B2, which mixes 10x 3' v2 and v3, it also prints each donor id's assay and its cells' median
depth (`raw_sum`, all genes) and genes detected (`nnz`), and the same per assay: the owner asked
(2026-10-11) whether depth differs by chemistry, since all of GATE 0's refusals fell on N7, B2's
null. `b2_depth_by_chemistry.py` joins these lines with the compact B2 and v1's reports.

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
        per_donor(census, *GROUPS["B2"])
    return 0


def per_donor(census, organism: str, flt: str) -> None:
    """B2's donor ids: assay, cells, median raw_sum and nnz; then the same per assay."""
    obs = (census["census_data"][organism].obs
           .read(value_filter=f"{CENSUS_FILTER} and {flt}",
                 column_names=["donor_id", "assay", "raw_sum", "nnz"])
           .concat().to_pandas())
    print("B2 per donor id (all cells of the chosen group; raw_sum over all genes):")
    for donor, g in obs.groupby("donor_id", observed=True):
        assays = {str(a): int(n) for a, n in g["assay"].value_counts().items() if n}
        print(f"  donor {donor} | cells {len(g)} | assay {assays} | median raw_sum "
              f"{g['raw_sum'].median():.0f} | median nnz {g['nnz'].median():.0f}")
    print("B2 per assay:")
    for assay, g in obs.groupby("assay", observed=True):
        print(f"  assay {assay} | donor ids {g['donor_id'].nunique()} | cells {len(g)} | median raw_sum "
              f"{g['raw_sum'].median():.0f} | median nnz {g['nnz'].median():.0f}")


if __name__ == "__main__":
    raise SystemExit(main())
