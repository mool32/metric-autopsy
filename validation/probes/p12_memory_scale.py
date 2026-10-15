"""Does the engine run at atlas scale? Sparse AnnData, realistic mouse var_names (with mt-/Rps genes)."""
import resource, sys, time
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
from functools import partial
from metric_autopsy import metrics, gate0_independence, gate1_qc_parity

n_cells, n_genes = int(sys.argv[1]), 20000
rng = np.random.default_rng(0)
X = sp.random(n_cells, n_genes, density=0.05, format="csr", random_state=0,
              data_rvs=lambda k: rng.poisson(3, k) + 1).astype(np.float32)
names = ["Smad3", "Col1a1", "Actb", "Gapdh", "mt-Co1", "Rps3"] + [f"g{i}" for i in range(n_genes - 6)]
obs = pd.DataFrame({"age": rng.choice(["young", "old"], n_cells), "sex": rng.choice(["m", "f"], n_cells)},
                   index=[f"c{i}" for i in range(n_cells)])
obs["n_genes_by_counts"] = np.diff(X.indptr); obs["total_counts"] = np.asarray(X.sum(1)).ravel()
a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=names))
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
print(f"{n_cells} cells x {n_genes} genes, sparse nnz={X.nnz/1e6:.1f}M ({X.data.nbytes/1e6+X.indices.nbytes/1e6:.0f} MB); peak RSS so far {rss():.2f} GB")
t = time.time(); gate1_qc_parity(a, "age", ("young", "old"), within=["sex"])
print(f"  after GATE 1 (QC cols precomputed in obs!): peak RSS {rss():.2f} GB  ({time.time()-t:.0f}s)")
t = time.time(); gate0_independence(partial(metrics.norm_pearson, gene_a="Smad3", gene_b="Col1a1"), a, n_baseline=5, n_perturb=3)
print(f"  after GATE 0 (5 bootstraps, 3 perturbations each): peak RSS {rss():.2f} GB  ({time.time()-t:.0f}s)")
dense_gb = n_cells * n_genes * 8 / 1e9
print(f"  dense float64 copy of X = {dense_gb:.2f} GB; TMS FACS (110,824 x ~22,900) would be {110824*22900*8/1e9:.0f} GB per copy")
