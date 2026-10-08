"""Simulated backgrounds for the panel code's tests, the timing pilot and the workflow's dry run.

Not part of the panel: the confirmatory panel runs on the real backgrounds of section 3.1.
Negative-binomial counts with cell-size variation; genes at the three expression levels of
``panel.level_of`` (mean about 4-8, 1-1.6 and 0.15-0.45 counts per cell), each level with
coupled pairs (a shared per-cell factor), and 20 of the G2M genes.
"""
from __future__ import annotations

import numpy as np

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


def dry_backgrounds(n_genes: int = 2500) -> dict:
    """B1 and B2 of the planned sizes (B1: 24 donors; B2: 12 mice), for timing and the dry run."""
    return {"B1": simulated_background("B1", donors=24, cells=220, n_genes=n_genes),
            "B2": simulated_background("B2", donors=12, cells=220, n_genes=n_genes, seed=1)}


def dry_pilot(bgs: dict, sesoi: float = 0.1, dose: float = 2.0) -> dict:
    """A stand-in pilot.json with every field the panel and the scoring read (not the oracle's
    pilot): SESOI `sesoi` and saturation dose `dose` everywhere, the truth 'valid' for the high
    and medium pairs and 'blind' for the low ones (as probe p14 found on these backgrounds),
    every case establishable."""
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
    out["delta"] = {str(k): {f"{f:g}": dict(value=0.2 * f, se=0.01) for f in (0.25, 0.5, 1.0, 1.5)}
                    for k in range(len(bgs["B1"].plan["pool"]))}
    return out
