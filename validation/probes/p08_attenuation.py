"""TRUE coupling, no confound anywhere. Does GATE 0 fail a sensible correlation metric just
because measurement noise attenuates it (reliability), at realistic 10x expression levels?"""
from probe_common import *
from probe_sim import *

rng0 = np.random.default_rng(5)
mu = base_means(rng0)
order = np.argsort(mu)
def pick(target):
    cand = [i for i in order if i >= N_S + N_G2M]
    return min(cand, key=lambda i: abs(mu[i] - target))

def lognorm_pearson(data, *, gene_a, gene_b):
    from metric_autopsy.core import as_dense, unique_col_index
    X = as_dense(data.X); v = list(data.var_names)
    L = lognorm(X)
    a, b = L[:, unique_col_index(v, gene_a)], L[:, unique_col_index(v, gene_b)]
    if a.std() == 0 or b.std() == 0: return 0.0
    return float(np.corrcoef(a, b)[0, 1])

print("mean UMI/cell of the coupled pair -> GATE 0 status, worst perturbation, rel shift")
for target in [0.3, 1.0, 3.0, 10.0]:
    i = pick(target); j = pick(target * 1.1)
    if j == i: j = order[list(order).index(i) + 1]
    rng = np.random.default_rng(11)
    X = simulate_cells(rng, 2000, mu, phase="G1", coupled=(i, j), coupling=0.8)
    d = SimpleData(X, pd.DataFrame({"age": ["young"] * 1000 + ["old"] * 1000}), GENES)
    ga, gb = GENES[i], GENES[j]
    for name, fn in [("norm_pearson", metrics.norm_pearson), ("lognorm_pearson(all cells)", lognorm_pearson)]:
        m = partial(fn, gene_a=ga, gene_b=gb)
        g = gate0_independence(m, d)
        rs = g.detail["responses"]
        w = max(rs, key=lambda k: rs[k]["rel_change"])
        print(f"  mean~{target:5.1f}  {name:27s} base={g.detail['baseline']:.3f} -> {g.status.value:4s}  "
              + "  ".join(f"{k[:9]}={rs[k]['rel_change']:.0%}" for k in rs))
