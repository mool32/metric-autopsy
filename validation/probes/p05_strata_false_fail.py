"""NULL data (no confound anywhere): how often do GATE 1 and GATE 5 fail anyway,
as a function of the number of strata and cells per stratum? (AGENTS.md: 'always stratify
by every factor'.)"""
from probe_common import *
from test_gates import _block, _assemble

def gate1_null_rate(K, n_per_group, reps=200, seed=0):
    rng = np.random.default_rng(seed)
    fails = 0
    for r in range(reps):
        rows = []
        for k in range(K):
            for grp in ("young", "old"):
                ng = np.round(np.exp(rng.normal(np.log(2000), 0.45, n_per_group)))
                rows += [dict(age=grp, stratum=f"s{k}", n_genes_by_counts=v) for v in ng]
        obs = pd.DataFrame(rows)
        d = SimpleData(np.ones((len(obs), 2)), obs, ["a", "b"])
        g = gate1_qc_parity(d, "age", ("young", "old"), within=["stratum"])
        fails += g.status == GateStatus.FAIL
    return fails / reps

print("GATE 1 false-FAIL rate on null data (identical QC distributions in both groups)")
print("  cells/group/stratum ->", [5, 10, 20, 50, 200])
for K in [1, 4, 16, 64]:
    rates = [gate1_null_rate(K, n, reps=100 if K < 64 else 40) for n in [5, 10, 20, 50, 200]]
    print(f"  K={K:3d} strata: " + "  ".join(f"{x:5.0%}" for x in rates))

def gate5_null_rate(K, n_per_stratum, reps=60, seed=0, pair_metric=metrics.norm_pearson):
    fails = 0
    for r in range(reps):
        rng = np.random.default_rng(seed + r)
        blocks = []
        for k in range(K):
            c, o = _block(rng, n_per_stratum, coupling=1.5, efficiency=1.0, age="young", sex="x")
            o["stratum"] = f"s{k}"
            blocks.append((c, o))
        d = _assemble(blocks)
        g = gate5_controls(pair_metric, d, ("Actb", "Gapdh"), ("Gene0", "Gene1"), within=["stratum"])
        fails += g.status == GateStatus.FAIL
    return fails / reps

print("\nGATE 5 false-FAIL rate on null data (controls behave identically in every stratum)")
for pm in [metrics.norm_pearson, metrics.mi_3bin]:
    print(f"  pair_metric={pm.__name__}; cells/stratum ->", [10, 30, 100, 400])
    for K in [1, 4, 16]:
        rates = [gate5_null_rate(K, n, pair_metric=pm) for n in [10, 30, 100, 400]]
        print(f"   K={K:3d} strata: " + "  ".join(f"{x:5.0%}" for x in rates))
