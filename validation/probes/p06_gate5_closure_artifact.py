from probe_common import *
from test_gates import _block, _assemble
vals = []
for r in range(30):
    rng = np.random.default_rng(r)
    d = _assemble([_block(rng, 400, coupling=1.5, efficiency=1.0, age="young", sex="x")])
    pos = metrics.norm_pearson(d, gene_a="Actb", gene_b="Gapdh")
    neg = metrics.norm_pearson(d, gene_a="Gene0", gene_b="Gene1")
    # same pair, but without library normalisation (codetected spearman on raw counts)
    neg_raw = metrics.codetected_spearman(d, gene_a="Gene0", gene_b="Gene1")
    # all filler pairs: mean induced correlation
    vals.append((pos, neg, neg_raw))
v = np.array(vals)
print("norm_pearson pos (Actb-Gapdh): mean %.3f   band neg_max = 0.2*pos = %.3f" % (v[:,0].mean(), 0.2*v[:,0].mean()))
print("norm_pearson neg (Gene0-Gene1): mean %.3f sd %.3f  -> fraction |neg|>band: %.0f%%" % (v[:,1].mean(), v[:,1].std(), 100*np.mean(np.abs(v[:,1]) > 0.2*v[:,0])))
print("codetected_spearman neg on raw counts: mean %.3f sd %.3f" % (v[:,2].mean(), v[:,2].std()))
# Average over all filler pairs to show the systematic induced correlation from CP10k closure
rng = np.random.default_rng(0)
d = _assemble([_block(rng, 400, coupling=1.5, efficiency=1.0, age="young", sex="x")])
fill = [g for g in GENES if g.startswith("Gene")]
rs = [metrics.norm_pearson(d, gene_a=fill[i], gene_b=fill[j]) for i in range(len(fill)) for j in range(i+1, len(fill))]
rr = [metrics.codetected_spearman(d, gene_a=fill[i], gene_b=fill[j]) for i in range(len(fill)) for j in range(i+1, len(fill))]
print("all %d independent filler pairs: norm_pearson mean r = %.3f ; raw codetected spearman mean = %.3f" % (len(rs), np.mean(rs), np.mean(rr)))
