"""A more realistic simulator than the 40-gene test generator: 10x-like depth, NB counts,
a cell-cycle module, and biology that can change RNA content (cycling cells are bigger)."""
import numpy as np, pandas as pd
from metric_autopsy import SimpleData

G = 1500
N_S, N_G2M = 40, 40

def base_means(rng, G=G):
    # log-normal gene means (relative abundances), scaled so a cell has ~3000 UMIs at size 1
    mu = np.exp(rng.normal(-1.0, 1.6, G))
    return mu / mu.sum() * 3000.0

def simulate_cells(rng, n, mu, *, phase, content=(1.0, 1.4, 1.9), cc_fold=4.0,
                   size_sd=0.35, theta=5.0, capture=1.0, coupled=None, coupling=0.0):
    """phase: array of 'G1'/'S'/'G2M' per cell (or a scalar). content: RNA content per phase.
    coupled: (i, j) gene indices sharing a latent log-scale factor with sd `coupling`."""
    if np.isscalar(phase):
        phase = np.array([phase] * n)
    lam = np.tile(mu, (n, 1))
    s_idx = np.arange(0, N_S); g_idx = np.arange(N_S, N_S + N_G2M)
    is_s = phase == "S"; is_g2m = phase == "G2M"
    lam[np.ix_(is_s, s_idx)] *= cc_fold
    lam[np.ix_(is_g2m, g_idx)] *= cc_fold
    # renormalize composition, then apply biological RNA content and technical size/capture
    lam = lam / lam.sum(axis=1, keepdims=True) * mu.sum()
    cont = np.where(is_s, content[1], np.where(is_g2m, content[2], content[0]))
    size = np.exp(rng.normal(0, size_sd, n))
    if coupled is not None and coupling > 0:
        z = rng.normal(0, 1, n)
        for g in coupled:
            lam[:, g] *= np.exp(coupling * z - coupling**2 / 2)
    lam = lam * (cont * size * capture)[:, None]
    # NB via gamma-poisson
    lam = rng.gamma(theta, lam / theta)
    return rng.poisson(lam).astype(float)

GENES = [f"S{i}" for i in range(N_S)] + [f"M{i}" for i in range(N_G2M)] + \
        [f"g{i}" for i in range(G - N_S - N_G2M)]

def make(blocks, seed=0, G_=G):
    rng = np.random.default_rng(seed)
    mu = base_means(rng)
    Xs, obs = [], []
    for spec in blocks:
        spec = dict(spec)
        n = spec.pop("n"); labels = spec.pop("obs")
        phase = spec.pop("phase")
        if callable(phase):
            phase = phase(rng, n)
        Xs.append(simulate_cells(rng, n, mu, phase=phase, **spec))
        o = pd.DataFrame({k: [v] * n for k, v in labels.items()})
        o["phase"] = phase if not np.isscalar(phase) else [phase] * n
        obs.append(o)
    return SimpleData(np.vstack(Xs), pd.concat(obs, ignore_index=True), GENES), mu

def lognorm(X):
    tot = X.sum(axis=1, keepdims=True); tot[tot == 0] = 1
    return np.log1p(X / tot * 1e4)

def g2m_score_cells(data, ctrl=None):
    """Tirosh/Seurat-style: mean lognorm of G2M genes minus mean of control genes."""
    from metric_autopsy.core import as_dense
    X = as_dense(data.X)
    names = list(data.var_names)
    idx = {g: i for i, g in enumerate(names)}
    L = lognorm(X)
    gm = [idx[g] for g in names if g.startswith("M")]
    if ctrl is None:
        ctrl = [g for g in names if g.startswith("g")][:200]
    ci = [idx[g] for g in ctrl if g in idx]
    return L[:, gm].mean(axis=1) - L[:, ci].mean(axis=1)

def mean_g2m_score(data):
    return float(np.mean(g2m_score_cells(data)))

def frac_cycling(data, t=0.25):
    return float(np.mean(g2m_score_cells(data) > t))
