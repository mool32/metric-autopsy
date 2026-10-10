from probe_common import *
from probe_sim import *

print("P3a: POSITIVE CONTROL — sorted G1 vs sorted G2M cells (FUCCI/Hoechst-style ground truth)")
print("     metric = mean G2M score (normalized, Tirosh-style). True effect is huge by construction.")
for g2m_content in [1.0, 1.4, 1.9, 2.4]:
    blocks = []
    for batch in ["b1", "b2"]:
        blocks.append(dict(n=500, obs=dict(sorted_phase="G1", batch=batch), phase="G1",
                           content=(1.0, 1.4, g2m_content)))
        blocks.append(dict(n=500, obs=dict(sorted_phase="G2M", batch=batch), phase="G2M",
                           content=(1.0, 1.4, g2m_content)))
    d, mu = make(blocks, seed=1)
    a = run_autopsy(mean_g2m_score, d, group_col="sorted_phase", groups=("G2M", "G1"),
                    within=["batch"], stop_on_first_fail=False)
    print(f"\n  G2M RNA content x{g2m_content}")
    show(a)
