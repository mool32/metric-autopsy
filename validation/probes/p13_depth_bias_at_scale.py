"""GATE 0 on log-normalized Pearson at a realistic scale (exploratory, found on v0.3.0.dev0).

A truly coupled pair (a shared latent factor), groups without biology, negative-binomial counts
with cell-size variation (log-SD 0.4) and 16 donors. On the 40-gene toy data `norm_pearson`
passes GATE 0. Here GATE 0's null (genes shuffled within depth bins) carries the CP10k ratio
correlation of p06, and it grows when depth is halved. GATE 0 sizes that shift in units of the
pair's signal above the null: once there are enough cells to resolve it (z > 4) and it exceeds
tol = 25% of the signal, the depth response is classified as bias and the metric is invalid.

Part 2 asks whether the same depth response would reach a verdict if it did not block: with a
capture loss in one group (identical biology, 8 vs 8 donors), the effect field alone (thinning to
equal depth, permutation over donors) should explain the raw difference by depth.

Decided 2026-10-08 (JOURNAL.md, D4): a depth bias that the declared correction removes between
groups is reported, not blocking; the regression is
test_probes.py::test_p13_norm_pearson_at_800_cells_per_donor_is_not_blocked. The log records the
finding under the earlier rule (git ac8f1a8); run directly, GATE 0 now reports these biases as
unsized (no SESOI is passed here) instead of failing.

    python p13_depth_bias_at_scale.py > p13_depth_bias_at_scale.log   # ~20 min
"""
import platform
import subprocess
from pathlib import Path

from probe_common import *  # noqa: F401,F403  (repo path, numpy, pandas, metrics, gates)

HERE = Path(__file__).resolve().parent
CONFIGS = ((2000, 200), (2000, 400), (2000, 800), (5000, 200))  # (genes, cells per donor)
SEEDS = range(4)


def simulate(n_genes, cells, donors=16, seed=0, capture_b=1.0):
    """Counts for `donors` donors in two label groups with identical biology.

    Genes A and B share a latent factor (log-scale loading 0.4); P1 and P2 another (0.6). Gene
    means are log-normal; each cell has a size factor (log-SD 0.4) and each donor a size offset
    (log-SD 0.15). NB dispersion r = 2. `capture_b` scales group B's means, which for NB counts
    equals binomial thinning (a pure capture loss).
    """
    rng = np.random.default_rng(seed)
    mu = np.exp(rng.normal(-1.0, 1.5, n_genes))
    names = [f"g{i}" for i in range(n_genes)]
    names[:4] = ["A", "B", "P1", "P2"]
    mu[:4] = [3.0, 3.0, 6.0, 6.0]
    X, rows = [], []
    for d in range(donors):
        size = np.exp(rng.normal(0, 0.4, cells)) * np.exp(rng.normal(0, 0.15))
        lat, lat2 = rng.normal(0, 1, cells), rng.normal(0, 1, cells)
        m = np.outer(size, mu) * (capture_b if d >= donors // 2 else 1.0)
        m[:, 0] *= np.exp(0.4 * lat)
        m[:, 1] *= np.exp(0.4 * lat)
        m[:, 2] *= np.exp(0.6 * lat2)
        m[:, 3] *= np.exp(0.6 * lat2)
        X.append(rng.negative_binomial(2.0, 2.0 / (2.0 + m)))
        rows += [dict(group="A" if d < donors // 2 else "B", donor=f"d{d}")] * cells
    X = np.vstack(X).astype(float)
    obs = pd.DataFrame(rows)
    obs["total_counts"] = X.sum(1)
    obs["n_genes_by_counts"] = (X > 0).sum(1)
    return SimpleData(X, obs, names)


def main():
    rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, capture_output=True,
                         text=True).stdout.strip() or "?"
    print(f"# p13 GATE 0 on norm_pearson at scale — git {rev}, python {platform.python_version()}, "
          f"numpy {np.__version__}")
    print("# per dataset: GATE 0 status / depth_downsample class (z of the metric's shift; null shift "
          "as a share of the signal above the null)")
    metric = partial(metrics.norm_pearson, gene_a="A", gene_b="B")
    for n_genes, cells in CONFIGS:
        out, n_fail = [], 0
        for seed in SEEDS:
            d = simulate(n_genes, cells, seed=seed)
            g0 = gate0_independence(metric, d, protect_genes=("A", "B"), seed=0)
            det = g0.detail
            r = det["responses"]["depth_downsample"]
            share = (f"{r['null_shift'] / det['signal']:+.0%}" if "null_shift" in r else "not tested")
            out.append(f"{g0.status.value}/{r['classification']} (z {r['z']:.1f}; null {share})")
            n_fail += g0.status == GateStatus.FAIL
            if seed == 0:
                print(f"  {n_genes} genes x {16 * cells} cells: null centre {det['null_center']:+.3f}, "
                      f"signal above null {det['signal']:.3f} (seed 0)")
        print(f"  {n_genes} genes, {cells} cells/donor: GATE 0 FAIL {n_fail}/{len(SEEDS)} — "
              + "; ".join(out), flush=True)

    from metric_autopsy.effect import estimate_effect
    from metric_autopsy.stats import fmt_rate
    print("# part 2: capture 0.5 in group B, identical biology, 2000 genes, 8 vs 8 donors x 400 cells; "
          "effect field only (composition estimand, SESOI 0.1)")
    n, explained, detected = 10, 0, 0
    for seed in range(n):
        d = simulate(2000, 400, seed=100 + seed, capture_b=0.5)
        eff, _ = estimate_effect(metric, d, group_col="group", groups=("A", "B"), replicate_col="donor",
                                 estimand="composition", sesoi=0.1, seed=seed)
        e = eff.detail
        explained += bool(e.get("explained_by_depth"))
        detected += eff.status == "DETECTED"
        print(f"  dataset {seed}: raw {e['raw_effect']:+.4f} (p={e.get('p_raw', float('nan')):.3f}) -> "
              f"at equal depth {e['effect']:+.4f}: {eff.status}"
              + ("; explained by depth" if e.get("explained_by_depth") else ""), flush=True)
    print(f"  explained by depth {fmt_rate(explained, n)}; corrected effect DETECTED {fmt_rate(detected, n)}")


if __name__ == "__main__":
    main()
