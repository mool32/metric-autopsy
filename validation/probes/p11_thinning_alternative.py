"""Prototype of an alternative to GATE 2's n_genes *cell selection*: equalize depth by
*binomial thinning* (quantile-matched library sizes) -- keeps every cell, removes only the
technical depth difference. Does it get all four ground-truth cases right?"""
from probe_common import *
from test_gates import _block, _assemble
from probe_sim import *
from metric_autopsy.cli import demo_data
from metric_autopsy.core import as_dense

def thin_to_match(data, group_col, groups, rng):
    X = as_dense(data.X).copy()
    g = np.asarray(data.obs[group_col])
    a, b = g == groups[0], g == groups[1]
    tot = X.sum(1)
    hi, lo = (a, b) if np.median(tot[a]) >= np.median(tot[b]) else (b, a)
    hi_idx = np.where(hi)[0]
    q = (np.argsort(np.argsort(tot[hi_idx])) + 0.5) / len(hi_idx)       # quantile of each high-depth cell
    target = np.quantile(tot[lo], q)
    p = np.clip(target / np.maximum(tot[hi_idx], 1), 0, 1)
    X[hi_idx] = rng.binomial(X[hi_idx].astype(np.int64), p[:, None]).astype(float)
    return SimpleData(X, data.obs, data.var_names)

def report(name, metric, data, group_col, groups, truth):
    rng = np.random.default_rng(0)
    eff = lambda d: metric(d[np.asarray(d.obs[group_col]) == groups[0]]) - metric(d[np.asarray(d.obs[group_col]) == groups[1]])
    e0 = eff(data)
    e_thin = np.mean([eff(thin_to_match(data, group_col, groups, rng)) for _ in range(5)])
    m = gate2_ngenes_matching(metric, data, group_col, groups)
    e_sel = m.detail.get("matched_effect", float("nan"))
    print(f"  {name:44s} truth={truth:9s} effect={e0:+.3f} | n_genes-selection: {e_sel:+.3f} ({abs(e_sel)/abs(e0):4.0%}) [{m.status.value}] "
          f"| thinning: {e_thin:+.3f} ({abs(e_thin)/abs(e0):4.0%})")

print("retained fraction of the between-group effect after depth equalization\n")
# 1. real biology that moves QC (proliferation decline; cycling cells carry 2.4x RNA)
def phases(p_s, p_g2m):
    def f(rng, n):
        u = rng.random(n); return np.where(u < p_s, "S", np.where(u < p_s + p_g2m, "G2M", "G1"))
    return f
blocks = [dict(n=600, obs=dict(age="young"), phase=phases(0.15, 0.20), content=(1.0, 1.6, 2.4)),
          dict(n=600, obs=dict(age="old"), phase=phases(0.02, 0.03), content=(1.0, 1.6, 2.4))]
d, _ = make(blocks, seed=2)
report("proliferation decline, mean G2M score", mean_g2m_score, d, "age", ("young", "old"), "REAL")
# 2. pure QC artifact: demo male stratum, mi_3bin (biology identical, old-male capture 30%)
demo = demo_data(); males = demo[np.asarray(demo.obs["sex"]) == "male"]
report("demo male stratum, mi_3bin (QC artifact)", partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1"),
       males, "age", ("young", "old"), "ARTIFACT")
# 3. Xist in the QC-degraded old stratum
rng = np.random.default_rng(7)
female = np.asarray(demo.obs["sex"]) == "female"; old = np.asarray(demo.obs["age"]) == "old"
eff = np.where((~female) & old, 0.30, 1.0)
xist = rng.poisson(np.where(female, 25.0, 0.0) * eff).astype(float)
dx = SimpleData(np.column_stack([demo.X, xist]), demo.obs, list(demo.var_names) + ["Xist"])
olds = dx[old]
def mean_lognorm_xist(data):
    X = as_dense(data.X); tot = X.sum(1); tot[tot == 0] = 1
    return float(np.mean(np.log1p(X[:, list(data.var_names).index("Xist")] / tot * 1e4)))
report("old stratum, Xist female vs male", mean_lognorm_xist, olds, "sex", ("female", "male"), "REAL")
# 4. real coupling difference with QC confound layered on top (old-male degraded)
rng = np.random.default_rng(3)
mix = _assemble([_block(rng, 400, 2.0, 1.0, "young", "male"), _block(rng, 400, 0.8, 0.3, "old", "male")])
report("coupling loss + QC confound, norm_pearson", partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1"),
       mix, "age", ("young", "old"), "REAL")
report("coupling loss + QC confound, mi_3bin", partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1"),
       mix, "age", ("young", "old"), "MIXED")
