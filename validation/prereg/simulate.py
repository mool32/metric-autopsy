"""Simulated backgrounds for the panel code's tests, the timing pilot and the workflow's dry run.

Not part of the panel: the confirmatory panel runs on the real backgrounds of section 3.1.
Negative-binomial counts with cell-size variation; genes at the three expression levels of
``panel.level_of`` (mean about 4-8, 1-1.6 and 0.15-0.45 counts per cell), each level with
coupled pairs (a shared per-cell factor), and 20 of the G2M genes.
"""
from __future__ import annotations

import numpy as np

import panel as P

DRY_RUN_KEY = "0123456789abcdef0123456789abcdef"  # public: the dry run's assignment is not the panel's


def simulated_background(name: str = "B1", donors: int = 20, cells: int = 220, n_genes: int = 420,
                         seed: int = 0, coupled: int = 6, plan: bool = True) -> P.Background:
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
    """A stand-in pilot.json with every field the panel reads (not the oracle's pilot)."""
    pool = bgs["B1"].plan["pool"]
    return dict(dry_run=True, sesoi={lv: sesoi for lv in P.LEVELS}, key_dose={lv: dose for lv in P.LEVELS},
                delta={str(k): {f"{f:g}": dict(value=0.2 * f, se=0.01) for f in (0.25, 0.5, 1.0, 1.5)}
                       for k in range(len(pool))},
                establishable={}, dropped=[])
