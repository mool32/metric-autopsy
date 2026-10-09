"""Simulated backgrounds for the panel code's tests, the timing pilot and the workflow's dry run.

Not part of the panel: the confirmatory panel runs on the real backgrounds of section 3.1.
Negative-binomial counts with cell-size variation; genes at the three expression levels of
``panel.level_of`` (mean about 4-8, 1-1.6 and 0.15-0.45 counts per cell), each level with
coupled pairs (a shared per-cell factor), and 20 of the G2M genes.

The dry run's stand-in for the selection writes them as the files and backgrounds.json entries
that select_backgrounds.py writes, so that the pilot and everything after it read files, as in
the real run:

    python simulate.py write --out-dir data --spec backgrounds.json --report selection.md
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import panel as P

# public: the dry run's assignment is not the panel's (the panel's key is a drand round's
# randomness named in the run tag before it exists; beacon.py)
DRY_RUN_KEY = "0123456789abcdef" * 4


def simulated_background(name: str = "B1", donors: int = 20, cells: int = 220, n_genes: int = 420,
                         seed: int = 0, coupled: int = 10, plan: bool = True) -> P.Background:
    rng = np.random.default_rng(seed)
    third = n_genes // 3
    mu = np.concatenate([rng.uniform(4, 8, third), rng.uniform(1.0, 1.6, third),
                         rng.uniform(0.15, 0.45, n_genes - 2 * third)])
    genes = [f"GENE{i}" for i in range(n_genes)]
    for k, g in enumerate(P.G2M_GENES[:20]):  # module genes among the high and medium genes
        genes[third - 40 + 2 * k] = g
    pairs = [(lv * third + 2 * k, lv * third + 2 * k + 1) for lv in range(3) for k in range(coupled)]
    X, donor = [], []
    for d in range(donors):
        size = np.exp(rng.normal(0, 0.35, cells)) * np.exp(rng.normal(0, 0.1))
        m = np.outer(size, mu)
        for a, b in pairs:
            z = rng.normal(0, 1, cells)
            m[:, a] *= np.exp(0.5 * z)
            m[:, b] *= np.exp(0.5 * z)
        X.append(rng.negative_binomial(3.0, 3.0 / (3.0 + m)))
        donor += [f"{name}d{d}"] * cells
    bg = P.Background(np.vstack(X).astype(float), genes, np.asarray(donor), name)
    if plan:
        P.plan_background(bg)
    return bg


def dry_backgrounds(n_genes: int = 2500, plan: bool = True) -> dict:
    """B1 and B2 of the planned sizes (B1: 24 donors; B2: 12 mice), for timing and the dry run."""
    return {"B1": simulated_background("B1", donors=24, cells=220, n_genes=n_genes, plan=plan),
            "B2": simulated_background("B2", donors=12, cells=220, n_genes=n_genes, seed=1, plan=plan)}


def b3_stand_in(seed: int = 3, per_phase: int = 96, n_genes: int = 300, ercc: int = 100):
    """A B3-shaped table (cells sorted into G1, S and G2M, one capture batch each, ERCC rows): the
    counts, the genes and the cells' phase and batch."""
    rng = np.random.default_rng(seed)
    genes = [g.capitalize() for g in P.G2M_GENES[:20]] + [f"Gene{i}" for i in range(n_genes - 20)] \
        + [f"ERCC-{i:05d}" for i in range(ercc)]
    mu = rng.uniform(0.2, 6.0, len(genes))
    phases = ("G1", "S", "G2M")
    X = np.vstack([rng.negative_binomial(3.0, 3.0 / (3.0 + np.outer(np.exp(rng.normal(0, 0.3, per_phase)), mu)))
                   for _ in phases])
    obs = pd.DataFrame(dict(phase=np.repeat(phases, per_phase), batch=np.repeat(phases, per_phase)))
    return X, genes, obs


def write_backgrounds(out_dir: Path, spec_path: Path, report: Path | None = None,
                      backgrounds: dict | None = None) -> dict:
    """The dry run's selection: B1 and B2 (default: `dry_backgrounds`) as .npz files with the
    backgrounds.json entries select_backgrounds.py writes (file, sha256, donor column, counts),
    and B3 - the rehearsal's E-MTAB-2805 file where it left one in `out_dir`, else a simulated
    stand-in - which the panel's scripts must skip (only anchors.py reads B3)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    spec, lines = {}, ["Dry run: simulated backgrounds (simulate.py), not the rule of section 3.1."]
    for name, bg in (backgrounds or dry_backgrounds(plan=False)).items():
        path = out_dir / f"{name}.npz"
        P.save_npz(path, bg.X, pd.DataFrame({"donor_id": bg.donor}), bg.genes)
        spec[name] = dict(file=path.name, sha256=P.sha256(path), donor="donor_id", counts="X",
                          source=dict(simulated="simulate.dry_backgrounds", donors=len(set(bg.donor)),
                                      cells=int(bg.X.shape[0]), genes=len(bg.genes)))
        lines.append(f"{name}: simulated, {len(set(bg.donor))} donors, {bg.X.shape[0]} cells x {len(bg.genes)} genes")
    b3 = out_dir / "B3.npz"
    if b3.exists():
        origin = "the rehearsal's E-MTAB-2805 download"
    else:
        X, genes, obs = b3_stand_in()
        P.save_npz(b3, X, obs, genes)
        origin = "a simulated stand-in (simulate.b3_stand_in)"
    spec["B3"] = dict(file=b3.name, sha256=P.sha256(b3), donor="batch", counts="X", source=dict(dry_run=origin))
    lines.append(f"B3: {origin}")
    spec_path = Path(spec_path)
    spec_path.write_text(json.dumps(spec, indent=1, sort_keys=True))
    if report is not None:
        Path(report).write_text("# Background selection (dry run)\n\n" + "\n".join(lines) + "\n")
    return spec


def dry_pilot(bgs: dict, sesoi: float = 0.1, dose: float = 2.0) -> dict:
    """A stand-in pilot.json with every field the panel and the scoring read (not the oracle's
    pilot): SESOI `sesoi` and saturation dose `dose` everywhere, the truth 'valid' for the high
    and medium pairs and 'blind' for the low ones (as probe p14 found on these backgrounds), Δ*
    above the SESOI (E2's and E3's at their depth too); no case is marked establishable (the
    oracle has not run)."""
    out = dict(dry_run=True, dropped=[], pool={}, pool_size={}, sesoi={}, saturation_dose={}, truth={},
               e_dose={lv: dose for lv in P.LEVELS}, key_dose={lv: dose for lv in P.LEVELS},
               key_dose_found={lv: True for lv in P.LEVELS}, establishable={})
    for name, bg in bgs.items():
        pool = bg.plan["pool"]
        out["pool"][name] = [dict(index=pe["index"], level=pe["level"], pair=pe["pair"]) for pe in pool]
        out["pool_size"][name] = len(pool)
        out["sesoi"][name] = {lv: sesoi for lv in P.LEVELS}
        out["saturation_dose"][name] = {lv: dose for lv in P.LEVELS}
        dm = P.DELTA_MIN_FRACTION * sesoi
        out["truth"][name] = {str(k): dict(response=(0.0 if pe["level"] == "low" else 4 * dm), se=0.001,
                                           dose=dose, delta_min=dm,
                                           **{"class": "blind" if pe["level"] == "low" else "valid"})
                              for k, pe in enumerate(pool)}
    out["delta"] = {str(k): {**{f"{f:g}": dict(value=0.2 * f, se=0.01) for f in (0.25, 0.5, 1.0, 1.5)},
                             f"capture={P.E_CAPTURE:g}": dict(value=0.13, se=0.01)}  # E2, E3 at their depth
                    for k in range(len(bgs["B1"].plan["pool"]))}
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description="simulated backgrounds for the dry run")
    sub = p.add_subparsers(dest="cmd", required=True)
    w = sub.add_parser("write", help="the dry run's selection: simulated files and backgrounds.json")
    w.add_argument("--out-dir", type=Path, required=True)
    w.add_argument("--spec", type=Path, required=True)
    w.add_argument("--report", type=Path)
    args = p.parse_args(argv)
    spec = write_backgrounds(args.out_dir, args.spec, args.report)
    print(json.dumps({k: dict(file=v["file"], sha256=v["sha256"]) for k, v in spec.items()}, indent=1))


if __name__ == "__main__":
    main()
