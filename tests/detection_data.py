"""The synthetic data of v0.4's detection (the owner's decisions of 2026-10-19): mice nested in age whose old
group is deeper but less complex, as on preprint §4 (a block of five genes takes most of its counts, so the
other genes get fewer: more counts and fewer genes per cell). The generator of
`validation/exploratory/differential_sensitivity.py` (the rejected change 17) and of the acceptance of change
17' (`validation/exploratory/detection_equalization.py`); used by `tests/test_v04_detection.py`.

`capture_sd` (added for the acceptance's case (vi), its default the generator's 0.3): the sd of each cell's
log capture, per group; a wider spread in one group makes the dropout correlated within its cells. The draws
are the same as before at the default (numpy scales its standard normals)."""
from __future__ import annotations

from functools import partial

import numpy as np
import pandas as pd

from metric_autopsy import SimpleData, metrics

N_FILLER = 150
GENES = ["A", "B"] + [f"D{i}" for i in range(5)] + [f"F{i}" for i in range(N_FILLER)]
PREREG = {"estimand": "composition", "direction": "decrease", "sesoi": 0.1, "judgment_pending": False}


def deeper_less_complex(rho=(0.5, 0.5), dom=(0.1, 0.75), depth=(3000.0, 4200.0), mu_claim=-3.0, n_mice=8,
                        cells=300, var_bio=0.8, seed=0, capture_sd=(0.3, 0.3)):
    """Mice nested in age (young, old). The pair A, B: log expected expression mu_claim + a shared factor whose
    share of the variance `var_bio` is the coupling rho (a correlation of log expected expression, the
    coupling probe's construct; the genes' variance stays `var_bio`); 150 fillers of lower expression. Per
    group a depth, the share `dom` of the counts that five dominant genes take, and the sd of each cell's log
    capture: the old group deeper but less complex by default (more counts, fewer genes per cell, the pair
    detected less), as on §4."""
    rng = np.random.default_rng(seed)
    fill_mu = np.random.default_rng(123).normal(-3.0, 1.5, N_FILLER)
    blocks, rows = [], []
    for gi, age in enumerate(("young", "old")):
        for m in range(n_mice):
            load, noise = np.sqrt(rho[gi] * var_bio), np.sqrt((1 - rho[gi]) * var_bio)
            lat = rng.normal(size=cells)
            la = mu_claim + load * lat + noise * rng.normal(size=cells)
            lb = mu_claim + load * lat + noise * rng.normal(size=cells)
            lf = fill_mu[None, :] + rng.normal(0, 0.5, (cells, N_FILLER))
            rel = np.column_stack([np.exp(la), np.exp(lb), np.exp(lf)])
            rel /= rel.sum(axis=1, keepdims=True)
            D = depth[gi] * np.exp(rng.normal(0, capture_sd[gi], cells))
            rest = (1 - dom[gi]) * D[:, None] * rel
            lam = np.column_stack([rest[:, :2], dom[gi] * D[:, None] * np.full((cells, 5), 0.2), rest[:, 2:]])
            blocks.append(rng.poisson(lam).astype(float))
            rows += [dict(age=age, mouse=f"{age}{m}")] * cells
    return SimpleData(np.vstack(blocks), pd.DataFrame(rows), GENES)


SCENARIOS = {
    "gap, equal coupling": dict(rho=(0.5, 0.5)),
    "equal detection, equal coupling": dict(rho=(0.5, 0.5), dom=(0.1, 0.1), depth=(3000.0, 3000.0)),
    "equal detection, decrease": dict(rho=(0.6, 0.3), dom=(0.1, 0.1), depth=(3000.0, 3000.0)),
    "gap, decrease": dict(rho=(0.6, 0.3)),
    "smaller gap, decrease": dict(rho=(0.6, 0.3), dom=(0.1, 0.6)),
}
METRICS = {"mi_3bin": (partial(metrics.mi_3bin, gene_a="A", gene_b="B"), -3.0),
           "norm_pearson": (partial(metrics.norm_pearson, gene_a="A", gene_b="B"), -2.5)}
