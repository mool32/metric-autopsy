# Datasets — manifest + how to re-fetch

> This file is what makes deleting `data/` safe: every dataset here can be restored from its
> source. Commit this file; do NOT commit the data itself (see .gitignore).

The engine and tests need **no external data** — they run on a synthetic planted confound
(`metric_autopsy.cli.demo_data`, `metric-autopsy --demo`). The datasets below are only for the
worked example in `examples/mi_coupling_tms/`, which reproduces the real MI-coupling autopsy.

| Dataset | Tier | Source / accession | Size | Re-fetch |
|---------|------|--------------------|------|----------|
| Tabula Muris Senis (FACS/SmartSeq2) | raw | figshare / CELLxGENE Census `tabula-muris-senis` | ~2–5 GB | `python examples/mi_coupling_tms/download_data.py --which tms` |
| Human skin fibroblasts (10x) | raw | CELLxGENE `cellxgene.cziscience.com` (skin, 179 donors) | ~1–2 GB | `python examples/mi_coupling_tms/download_data.py --which skin` |

## Tiers
- **raw** — public, re-downloadable; delete freely, restore from source.
- **derived** — computed by `src/` / example scripts; re-generate by re-running them.
- **irreplaceable** — none in this project (all inputs are public).

## Genes used by the worked example
`Smad3`, `Col1a1` (the coupling pair under test) · `Actb`, `Gapdh` (positive/housekeeping
control) · a random expression-matched pair (negative control). Factor columns expected in
`obs`: `age`, `sex`, and where available `tissue`, `cell_type`.

## Re-fetch all
```bash
python examples/mi_coupling_tms/download_data.py --which all
```
