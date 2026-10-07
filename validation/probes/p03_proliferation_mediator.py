from probe_common import *
from probe_sim import *

def phases(p_s, p_g2m):
    def f(rng, n):
        u = rng.random(n)
        return np.where(u < p_s, "S", np.where(u < p_s + p_g2m, "G2M", "G1"))
    return f

print("P3b: REAL BIOLOGY THAT MOVES QC — young tissue proliferates (35% cycling), old does not (5%).")
print("     Cycling cells carry more RNA (biology). metric = fraction of G2M-high cells / mean G2M score.")
for content in [(1.0, 1.4, 1.9), (1.0, 1.6, 2.4)]:
    blocks = []
    for sex in ["male", "female"]:
        blocks.append(dict(n=600, obs=dict(age="young", sex=sex), phase=phases(0.15, 0.20), content=content))
        blocks.append(dict(n=600, obs=dict(age="old", sex=sex), phase=phases(0.02, 0.03), content=content))
    d, mu = make(blocks, seed=2)
    for name, m in [("frac_cycling", frac_cycling), ("mean_g2m_score", mean_g2m_score)]:
        a = run_autopsy(m, d, group_col="age", groups=("young", "old"), within=["sex"],
                        stop_on_first_fail=False)
        print(f"\n  content(G1,S,G2M)={content}  metric={name}")
        show(a)
