"""Same-data positive control: on the project's own demo data (male-old QC-degraded, as in
TMS), add Xist (expressed in female cells only). Pre-registered claim: 'Xist is higher in
female cells'. Ground truth: certain. What does the autopsy say?"""
from probe_common import *
from metric_autopsy.cli import demo_data
from metric_autopsy.core import as_dense

d = demo_data()
rng = np.random.default_rng(7)
female = np.asarray(d.obs["sex"]) == "female"
old_male = (~female) & (np.asarray(d.obs["age"]) == "old")
eff = np.where(old_male, 0.30, 1.0)                      # same capture efficiency as the demo's confound
xist = rng.poisson(np.where(female, 25.0, 0.0) * eff).astype(float)
X = np.column_stack([d.X, xist])
dx = SimpleData(X, d.obs, list(d.var_names) + ["Xist"])

def mean_lognorm_xist(data):
    X = as_dense(data.X); tot = X.sum(1); tot[tot == 0] = 1
    j = list(data.var_names).index("Xist")
    return float(np.mean(np.log1p(X[:, j] / tot * 1e4)))

print("metric = mean log1p(CP10k) of Xist; group_col=sex (female vs male); within=['age']")
a = run_autopsy(mean_lognorm_xist, dx, group_col="sex", groups=("female", "male"), within=["age"],
                stop_on_first_fail=False)
show(a)
print("\nper-stratum QC table from GATE 1:")
g1 = next(r for r in a.results if r.gate == 1)
for row in g1.detail["table"]:
    print("  ", row["stratum"], f"median n_genes female={row['median_n_genes_a']:.0f} male={row['median_n_genes_b']:.0f} ratio={row['n_genes_ratio']:.2f} overlap={row['overlap']:.2f} flagged={row['flagged']}")
print("\nArithmetic on the manuscript's real TMS table (§4): old stratum female/male = 3413/1701 = %.2fx > 1.5x" % (3413/1701))
