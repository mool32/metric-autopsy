"""NULL at the level that matters: no age effect, but mice differ from each other (biological
replicate variance). 3 young + 3 old mice. Does the tool certify an age effect?"""
from probe_common import *
from test_gates import _block, _assemble

def one(seed, n_mice=3, cells=300, mouse_sd=0.35):
    rng = np.random.default_rng(seed)
    blocks = []
    for age in ("young", "old"):
        for m in range(n_mice):
            c = rng.normal(1.5, mouse_sd)           # mouse-level coupling, SAME distribution for both ages
            X, o = _block(rng, cells, coupling=max(c, 0.05), efficiency=1.0, age=age, sex="female")
            o["mouse"] = f"{age}{m}"
            blocks.append((X, o))
    return _assemble(blocks)

npr = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
res = dict(survives=0, no_effect=0, other=0)
mouse_level_p = []
for s in range(40):
    d = one(s)
    g2 = gate2_ngenes_matching(npr, d, "age", ("young", "old"))
    if g2.status == GateStatus.PASS and "survives" in g2.message:
        res["survives"] += 1
    elif "nothing to preserve" in g2.message:
        res["no_effect"] += 1
    else:
        res["other"] += 1
    # what a replicate-aware test would say: per-mouse metric, exact permutation over mice
    import itertools
    mice = sorted(set(d.obs["mouse"]))
    per = {m: npr(d[np.asarray(d.obs["mouse"]) == m]) for m in mice}
    young = [per[m] for m in mice if m.startswith("young")]; old = [per[m] for m in mice if m.startswith("old")]
    obs_diff = np.mean(young) - np.mean(old)
    allv = young + old
    diffs = []
    for comb in itertools.combinations(range(6), 3):
        a = [allv[i] for i in comb]; b = [allv[i] for i in range(6) if i not in comb]
        diffs.append(np.mean(a) - np.mean(b))
    mouse_level_p.append(np.mean(np.abs(diffs) >= abs(obs_diff) - 1e-12))
print("40 null simulations (no age effect; 3 vs 3 mice; mouse-to-mouse sd in coupling):")
print("  GATE 2 'effect survives matching' (i.e. tool reports a real, QC-robust age effect):", res["survives"], "/ 40")
print("  GATE 2 'no meaningful effect':", res["no_effect"], "/ 40   other:", res["other"])
print("  replicate-level exact permutation p<0.05: %d / 40  (min attainable p with 3v3 = %.2f)"
      % (sum(p < 0.05 for p in mouse_level_p), 2/20))

d = one(0)
a = run_autopsy(npr, d, group_col="age", groups=("young", "old"), within=["mouse"] if False else [])
print("\nexample full autopsy on null seed 0 (Python API defaults):"); show(a)
