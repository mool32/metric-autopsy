# Worked example — MI coupling on Tabula Muris Senis

The center of gravity of this repo: the real analysis that motivated it, re-run through the
gates so you can watch a plausible metric **fail every applicable gate and get killed**.

## The claim under autopsy

> SMAD→ECM mutual information (`MI(Smad3; Col1a1)`, 3-bin zero/low/high) declines with age —
> "coupling loss" — while housekeeping coupling is preserved.

It is wrong, and it is wrong in an instructive way: the decline follows **the detection rate**.
Old male cells detect fewer genes than young ones, and the MI falls with detection in the same
mice; the audit ([`REPORT.md`](../../validation/flagship_audit/REPORT.md)) has the numbers and
what they do not establish (whether the drop is technical, and the female arm).

## What the gates find

| Gate | Result on `mi_3bin` | Why |
|---|---|---|
| 0 Mathematical independence | **FAIL** | MI with a zero-bin moves ~60% under simulated dropout — it *is* a detection metric |
| 1 QC parity | **FAIL** | the male stratum is flagged (1.94× on the synthetic demo, 1.65× on the real data); the female stratum passes |
| 2 n_genes matching | **STOP** (in males) | young/old n_genes distributions don't overlap — incomparable |
| 3 Raw visibility | dropout-driven | the "shape change" in the Smad3×Col1a1 scatter is points collapsing onto the zero axes |
| 5 Controls | HK also declines in males | the negative/housekeeping control is confounded the same way |
| 6 Replication | **not run by the tool** | by-hand result from the internal checklist only: human skin (10x) showed the effect unmatched, but n_genes-matched it vanished even in the young. The automated GATE 6 run on 10x skin has not been done. |

*Where the numbers come from:* 1.94× is the synthetic demo (`scripts/run_gates.py --demo`);
1.65× is the real-data run of the preprint's §4, pooled across tissues, which the audit
reproduces. Its age groups and replicate counts are in the audit; see
[`validation/flagship_audit/REPORT.md`](../../validation/flagship_audit/REPORT.md).

**Verdict: the metric was not measuring biology.** What *did* survive: the sex dimorphism in
cell quality itself (male-old cells detect far fewer genes) — real, whether biological or
technical — and `Actb–Gapdh` co-expression, which survives matching in both sexes.

## Run it

```bash
pip install -e "../..[dev]" cellxgene-census
python download_data.py --which all          # -> ../../data/raw/*.h5ad  (public, not committed)

# autopsy on the real data
metric-autopsy --h5ad ../../data/raw/tms_facs.h5ad \
    --metric mi_3bin --gene-a Smad3 --gene-b Col1a1 \
    --group-col age --groups '3m' '24m' --within sex tissue \
    --pos-pair Actb Gapdh --neg-pair Malat1 Xist \
    --data2 ../../data/raw/human_skin.h5ad --no-stop
```

No data yet? The same failure runs on a synthetic planted confound with zero downloads:

```bash
metric-autopsy --demo --no-stop
```

## Notebook

[`notebook.ipynb`](notebook.ipynb) is the narrated end-to-end: load → QC table + per-stratum
n_genes histogram → the pre-registration → gates 0/1/2/3/5 → the full autopsy and verdict, with
figures ([`figures/`](figures/)). It runs on the synthetic confound out of the box (executed
outputs are committed); point the data cell at `../../data/raw/tms_facs.h5ad` to run it on the
real TMS FACS data — every downstream cell is unchanged because the gates take any AnnData.

```bash
pip install -e "../..[dev]"
jupyter nbconvert --to notebook --execute --inplace notebook.ipynb   # reproduce outputs
```

## Status

- [x] Synthetic reproduction (`--demo` and the notebook) — the 0/N failure, no data required.
- [x] `download_data.py` — CELLxGENE Census pull (falls back to manual instructions).
- [x] `notebook.ipynb` — executed, with committed outputs and figures; real-TMS-ready.
- [ ] Turnkey GATE 6 second-platform replication on human skin (v1.1 — needs both datasets local).
