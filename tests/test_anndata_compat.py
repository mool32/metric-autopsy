"""AnnData interchangeability: the gates must give the SAME answer on a real AnnData
(with a sparse X and string obs index) as on the bundled SimpleData. Skipped if anndata
is not installed; the CI 'full' job runs it."""
from functools import partial

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")
sparse = pytest.importorskip("scipy.sparse")

from metric_autopsy import (
    SimpleData, GateStatus, metrics, run_autopsy,
    gate1_qc_parity, gate2_ngenes_matching,
)
from metric_autopsy.cli import demo_data


def _to_anndata(sd: SimpleData) -> "anndata.AnnData":
    import pandas as pd
    ad = anndata.AnnData(
        X=sparse.csr_matrix(sd.X),                      # sparse -> exercises as_dense()
        obs=sd.obs.copy(),
        var=pd.DataFrame(index=list(sd.var_names)),
    )
    ad.obs_names = [f"cell_{i}" for i in range(sd.n_obs)]  # string barcodes, like real data
    return ad


@pytest.fixture(scope="module")
def pair():
    sd = demo_data()
    return sd, _to_anndata(sd)


def test_pairwise_metric_matches(pair):
    sd, ad = pair
    for name in ("mi_3bin", "norm_pearson", "pearson"):
        fn = getattr(metrics, name)
        v_sd = fn(sd, gene_a="Smad3", gene_b="Col1a1")
        v_ad = fn(ad, gene_a="Smad3", gene_b="Col1a1")
        assert np.isclose(v_sd, v_ad, atol=1e-9), f"{name}: {v_sd} vs {v_ad}"


def test_gate1_agrees(pair):
    sd, ad = pair
    s = gate1_qc_parity(sd, "age", ("young", "old"), within=["sex"])
    a = gate1_qc_parity(ad, "age", ("young", "old"), within=["sex"])
    # the male-old QC confound is found on both (v0.3: a WARN diagnostic, not a FAIL)
    assert s.status == a.status == GateStatus.WARN
    assert ({r["stratum"]["sex"] for r in s.detail["flagged"]}
            == {r["stratum"]["sex"] for r in a.detail["flagged"]} == {"male"})


def test_gate2_agrees_on_male_stratum(pair):
    sd, ad = pair
    m = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")
    s = gate2_ngenes_matching(m, sd[np.asarray(sd.obs["sex"]) == "male"],
                              "age", ("young", "old")).status
    a = gate2_ngenes_matching(m, ad[np.asarray(ad.obs["sex"]) == "male"],
                              "age", ("young", "old")).status
    assert s == a == GateStatus.STOP


def test_full_autopsy_verdict_agrees(pair):
    """Same verdict, same four fields and same numbers on SimpleData and on AnnData. The demo
    has identical biology, so the verdict must not be SUPPORTED (v0.1.1 checked 'FAIL')."""
    sd, ad = pair
    m = partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")
    kw = dict(group_col="age", groups=("young", "old"), within=["sex"],
              gene_pair=("Smad3", "Col1a1"), replicate_col="mouse",
              pair_metric=metrics.mi_3bin, pos_pair=("Actb", "Gapdh"),
              neg_pair=("Gene0", "Gene1"), prereg={"estimand": "composition"},
              stop_on_first_fail=False)
    a_sd = run_autopsy(m, sd, **kw)
    a_ad = run_autopsy(m, ad, **kw)
    assert a_sd.verdict == a_ad.verdict
    assert not a_sd.verdict.startswith("SUPPORTED")
    assert ({k: v.status for k, v in a_sd.fields().items()}
            == {k: v.status for k, v in a_ad.fields().items()})
    assert np.isclose(a_sd.effect.detail["effect"], a_ad.effect.detail["effect"], atol=1e-12)
