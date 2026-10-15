from probe_common import *
d = make_clean()
print("P1a: constant metric lambda d: 0.0 on clean data, via run_autopsy defaults")
a = run_autopsy(lambda data: 0.0, d, group_col="age", groups=("young","old"), within=["sex"])
show(a)

print("\nP1b: random-number metric (pure noise, no data dependence)")
rng = np.random.default_rng(0)
a = run_autopsy(lambda data: float(rng.normal()), d, group_col="age", groups=("young","old"), within=["sex"])
show(a)

print("\nP1c: GATE 0 location-invariance: mi_3bin vs mi_3bin + C on clean data")
mi = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")
for C in [0, 0.5, 1, 10]:
    g = gate0_independence(lambda data, C=C: mi(data) + C, d)
    w = max(g.detail["responses"].items(), key=lambda kv: kv[1]["rel_change"])
    print(f"  C={C:5}: GATE0 {g.status.value}  worst={w[0]} rel={w[1]['rel_change']:.3f} z={w[1]['z']:.1f}")
print("\nP1d: GATE 0 scale-invariance? mi_3bin * k")
for k in [1, 100, 0.01]:
    g = gate0_independence(lambda data, k=k: mi(data) * k, d)
    print(f"  k={k:6}: GATE0 {g.status.value}")
print("\nP1e: affine reparametrisation of norm_pearson: Fisher z / (1 - r)")
npr = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
for name, f in [("r", lambda x: npr(x)), ("1-r (a 'decoupling' score)", lambda x: 1 - npr(x)),
                ("atanh(r)", lambda x: float(np.arctanh(npr(x))))]:
    g = gate0_independence(f, d)
    w = max(g.detail["responses"].items(), key=lambda kv: kv[1]["rel_change"])
    print(f"  {name:28s}: GATE0 {g.status.value} worst={w[0]} rel={w[1]['rel_change']:.3f} base={g.detail['baseline']:.3f}")
