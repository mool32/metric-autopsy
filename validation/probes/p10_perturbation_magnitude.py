"""How much does the headline GATE 0 verdict on mi_3bin depend on the (uncalibrated) choice
of perturbation and its magnitude?"""
from probe_common import *
import metric_autopsy.gates as G
from metric_autopsy.cli import demo_data

orig = G._perturb
def make_perturb(drop_frac, thin_p):
    def _p(data, kind, rng, protect=frozenset()):
        X = data.X.copy()
        if kind == "extra_dropout":
            nz = np.argwhere(X > 0); k = int(drop_frac * len(nz))
            if k:
                pick = nz[rng.choice(len(nz), size=k, replace=False)]; X[pick[:, 0], pick[:, 1]] = 0.0
            return SimpleData(X, data.obs, data.var_names), {}
        if kind == "depth_downsample":
            X = rng.binomial(X.astype(np.int64), thin_p).astype(float)
            return SimpleData(X, data.obs, data.var_names), {}
        return orig(data, kind, rng, protect)
    return _p

d = demo_data()
mi = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")
print("demo data, mi_3bin. rows: extra_dropout fraction; cols: binomial thinning keep-prob")
print("              " + "".join(f"thin={p:<10}" for p in [0.8, 0.5, 0.3]))
for f in [0.05, 0.10, 0.20]:
    cells = []
    for p in [0.8, 0.5, 0.3]:
        G._perturb = make_perturb(f, p)
        g = G.gate0_independence(mi, d)
        rs = g.detail["responses"]
        cells.append(f"{g.status.value}({rs['extra_dropout']['rel_change']:.0%}/{rs['depth_downsample']['rel_change']:.0%})")
    print(f"  drop={f:4.2f}:  " + "  ".join(f"{c:18s}" for c in cells))
G._perturb = orig
