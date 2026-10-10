"""GATE 2 on TRUE, QC-free effects of modest size: how often is a real effect killed
(or waved through as 'no effect') purely by estimator noise? And how stable is the
20-permutation null floor?"""
from probe_common import *
from test_gates import _block, _assemble

npr = partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1")
for c_old in [1.5, 1.2, 0.9]:
    out = {"PASS-survives": 0, "PASS-no-effect": 0, "FAIL": 0, "STOP": 0}
    floors = []
    for s in range(40):
        rng = np.random.default_rng(100 + s)
        d = _assemble([_block(rng, 300, 1.5, 1.0, "young", "m"), _block(rng, 300, 1.5, 1.0, "young", "f"),
                       _block(rng, 300, c_old, 1.0, "old", "m"), _block(rng, 300, c_old, 1.0, "old", "f")])
        g = gate2_ngenes_matching(npr, d, "age", ("young", "old"), seed=s)
        floors.append(g.detail.get("min_effect", np.nan))
        if g.status == GateStatus.PASS:
            out["PASS-survives" if "survives" in g.message else "PASS-no-effect"] += 1
        else:
            out[g.status.value] += 1
    eff = g.detail["unmatched_effect"]
    print(f"young coupling 1.5 vs old {c_old}: (last effect {eff:+.3f})  {out}   "
          f"null floor across seeds: median {np.nanmedian(floors):.3f}, range {np.nanmin(floors):.3f}-{np.nanmax(floors):.3f}")
