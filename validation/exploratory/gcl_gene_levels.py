"""How the genes a top-expressed gene-set metric reads are spread over validation v1's expression
levels (exploratory; for the plan of real claims, validation/plans/real_claims.md).

The global coordination level (GCL; Levy et al. 2020, Nat Metab 2:1305-1315) is computed, in its
method paper (Amit et al. 2022, Sci Rep), on the top 2,000 expressed genes. v1's backgrounds keep
the most expressed genes of their cells, next to the genes v1 named (panel.N_GENES, 2,000 in all).
v1's pair metric was valid only at the high level (detected in a majority of cells, mean >= 2
counts per cell; panel.level_of); this script counts each background's 1,000 and 1,500 most expressed
genes by level. It says where a pair metric was blind, not whether GCL is.

    python validation/exploratory/gcl_gene_levels.py --compact COMPACT > gcl_gene_levels.log
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation" / "prereg"))
import numpy as np  # noqa: E402

import panel as P  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--compact", type=Path, required=True, help="the run's compact backgrounds (B1.npz, B2.npz)")
    args = ap.parse_args(argv)
    out = ["# The most expressed genes of v1's backgrounds by v1's expression level (exploratory)", ""]
    for name in P.BACKGROUNDS:
        path = args.compact / f"{name}.npz"
        bg = P.load_compact(path)
        mean, det = P._col_stats(bg.X)
        levels = np.array([P.level_of(m, d) or "below (detected in < 10% of cells)" for m, d in zip(mean, det)])
        order = np.argsort(-mean, kind="stable")
        out.append(f"{name}: {bg.X.shape[0]} cells, {bg.X.shape[1]} genes (file sha256 "
                   f"{hashlib.sha256(path.read_bytes()).hexdigest()[:16]}...)")
        named = set(bg.plan.get("g2m", [])) | set(bg.plan.get("random_genes", []))
        low_named = sum(g in named and d < 0.10 for g, d in zip(bg.genes, det))
        out.append(f"  the 2,000 kept genes include {len(named & set(bg.genes))} named ones (the G2M module, N6c's random "
                   f"genes), {low_named} of them detected in < 10% of cells, so the kept set is not the 2,000 most "
                   "expressed; counted here are the 1,000 and 1,500 most expressed")
        for label, idx in (("the 1,000 most expressed", order[:1000]), ("the 1,500 most expressed", order[:1500])):
            cnt = Counter(levels[idx])
            parts = [f"{lv} {cnt.get(lv, 0)} ({100 * cnt.get(lv, 0) / len(idx):.1f}%)"
                     for lv in ("high", "medium", "low", "below (detected in < 10% of cells)")]
            out.append(f"  {label:<26}: " + "; ".join(parts))
        out.append(f"  the 1,000th most expressed gene: mean {mean[order[999]]:.3f} counts per cell, detected in "
                   f"{100 * det[order[999]]:.1f}% of cells; the 1,500th: {mean[order[1499]]:.3f}, {100 * det[order[1499]]:.1f}%")
        out.append("")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
