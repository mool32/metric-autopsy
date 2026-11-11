"""v0.4 change 10 (validation/plans/v0.4_spec.md, section 10; coded on the owner's decision of 2026-11-05): sparse
input.

A scipy.sparse X is kept sparse, and the same data stored dense or sparse get the same report, number for number:
the data hash does not depend on the storage, every draw is keyed by it, and every whole-matrix draw runs on the
stored entries with the dense path's draws (a canonical CSR matrix stores its entries in the row-major order of
the dense matrix's nonzeros, and numpy's binomial draws nothing for a zero count). No path densifies the whole
matrix. Skipped without scipy (the CI's core-only job)."""
from __future__ import annotations

import json
from functools import partial

import numpy as np
import pandas as pd
import pytest

sp = pytest.importorskip("scipy.sparse")

from metric_autopsy import SimpleData, injected_signal, metrics, run_autopsy  # noqa: E402
from metric_autopsy import core, equalize, gates, qc, specificity  # noqa: E402
from metric_autopsy import injected_signal as inj  # noqa: E402
from metric_autopsy.cli import demo_data  # noqa: E402

from test_v03_verdict import log_total_endogenous, with_ercc  # noqa: E402
from test_v04_domain import typed_mice  # noqa: E402
from test_v04_set import SET, mean_log_level, mean_norm_pearson, set_data  # noqa: E402

PAIR = ("Smad3", "Col1a1")
PREREG = dict(estimand="composition", direction="decrease", sesoi=0.1)


def csr(d: SimpleData) -> SimpleData:
    return SimpleData(sp.csr_matrix(d.X), d.obs, d.var_names)


def report(a) -> dict:
    """The JSON report without the time it was written."""
    j = json.loads(json.dumps(a.to_dict(), default=str))
    j["provenance"].pop("timestamp_utc", None)
    return j



# --------------------------------------------------------------------------- #
# the container
# --------------------------------------------------------------------------- #
def test_a_sparse_x_is_held_as_a_canonical_csr_matrix_and_reads_as_the_dense_one():
    d = demo_data()
    for X in (sp.csr_matrix(d.X), sp.csc_matrix(d.X), sp.coo_matrix(d.X), sp.csr_array(d.X)):
        s = SimpleData(X, d.obs, d.var_names)
        assert isinstance(s.X, sp.csr_matrix) and s.X.has_canonical_format and s.X.dtype == np.float64
        assert np.array_equal(s.X.toarray(), d.X)
    s = csr(d)
    idx = np.random.default_rng(0).integers(0, d.n_obs, d.n_obs)  # a bootstrap: rows repeated
    assert core.issparse(s[idx].X) and np.array_equal(s[idx].X.toarray(), d[idx].X)
    assert np.array_equal(s.gene_vector("Col1a1"), d.gene_vector("Col1a1"))
    assert np.array_equal(core.row_sums(s.X), d.X.sum(axis=1))
    assert np.array_equal(core.n_detected(s.X), (d.X > 0).sum(axis=1))
    block = np.random.default_rng(1).poisson(3, (d.n_obs, 2)).astype(float)
    assert np.array_equal(core.with_columns(s.X, [5, 1], block).toarray(), core.with_columns(d.X, [5, 1], block))
    assert qc.per_cell_qc(s).equals(qc.per_cell_qc(d))


# --------------------------------------------------------------------------- #
# the whole-matrix draws: the dense path's, draw for draw
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("kind", ["extra_dropout", "depth_downsample", "library_scale", "gene_subsample"])
def test_gate0s_perturbations_draw_what_the_dense_path_draws(kind):
    d = demo_data()
    a, na = gates._perturb(d, kind, np.random.default_rng(3), frozenset(PAIR))
    b, nb = gates._perturb(csr(d), kind, np.random.default_rng(3), frozenset(PAIR))
    assert core.issparse(b.X) and na == nb and list(a.var_names) == list(b.var_names)
    assert np.array_equal(a.X, b.X.toarray())


def test_gate0s_null_and_the_claims_null_draw_what_the_dense_path_draws():
    d = demo_data()
    s = csr(d)
    assert gates._null_mode(s) == "depth-matched"
    a, b = gates._gene_null(d, np.random.default_rng(4)), gates._gene_null(s, np.random.default_rng(4))
    assert core.issparse(b.X) and np.array_equal(a.X, b.X.toarray())
    a, b = gates.claim_null(d, list(PAIR), np.random.default_rng(5)), gates.claim_null(s, list(PAIR), np.random.default_rng(5))
    assert np.array_equal(a.X, b.X.toarray())
    # chunks of 256 genes: the null of a wide matrix
    w = SimpleData(np.random.default_rng(6).poisson(0.2, (60, 600)).astype(float), d.obs.iloc[:60], [f"g{i}" for i in range(600)])
    assert np.array_equal(gates._gene_null(w, np.random.default_rng(7)).X, gates._gene_null(csr(w), np.random.default_rng(7)).X.toarray())


def test_non_count_input_is_perturbed_and_nulled_as_the_dense_path_does():
    d = demo_data()
    lg = SimpleData(np.log1p(d.X), d.obs, d.var_names)
    ls = csr(lg)
    assert gates._null_mode(ls) == "depth bins" and not equalize.looks_like_counts(ls.X)
    a, na = gates._perturb(lg, "depth_downsample", np.random.default_rng(1))
    b, nb = gates._perturb(ls, "depth_downsample", np.random.default_rng(1))
    assert na == nb and np.array_equal(a.X, b.X.toarray())
    bins = gates._depth_bins(lg)
    assert np.array_equal(bins, gates._depth_bins(ls))
    a, b = gates._gene_null(lg, np.random.default_rng(2), bins), gates._gene_null(ls, np.random.default_rng(2), bins)
    assert core.issparse(b.X) and np.array_equal(a.X, b.X.toarray())


def test_equal_depth_thinning_draws_what_the_dense_path_draws():
    d = demo_data()
    tot = d.X.sum(axis=1)
    kw = dict(totals=tot, within=["sex"])
    a, ia = equalize.thin_to_match(d, "age", ("young", "old"), rng=np.random.default_rng(6), **kw)
    b, ib = equalize.thin_to_match(csr(d), "age", ("young", "old"), rng=np.random.default_rng(6), **kw)
    assert core.issparse(b.X) and np.array_equal(a.X, b.X.toarray()) and ia == ib
    assert equalize.looks_like_counts(csr(d).X)
    big = SimpleData(np.random.default_rng(0).poisson(0.3, (500, 500)).astype(float), d.obs.iloc[:500], [f"g{i}" for i in range(500)])
    big.X[7, 9] = 0.5  # non-count entries the sampled check may or may not meet: the dense check's answer either way
    assert equalize.looks_like_counts(csr(big).X) == equalize.looks_like_counts(big.X)


PROBES = {
    "coupling": lambda: injected_signal.coupling(*PAIR, strength=1.3),
    "coupling_set": lambda: injected_signal.coupling_set(["Smad3", "Col1a1", "Actb", "Gapdh"], strength=1.1, couple="half"),
    "module": lambda: injected_signal.module(["Actb", "Gapdh"], fold=2.0, frac=0.6),
}


@pytest.mark.parametrize("name", list(PROBES))
def test_injections_draw_what_the_dense_path_draws(name):
    d = demo_data()
    s = csr(d)
    probe = PROBES[name]()
    m = partial(metrics.norm_pearson, gene_a="Actb", gene_b="Gapdh")
    a = gates.injected_deltas(m, d, probe, 6, np.random.default_rng(8))
    b = gates.injected_deltas(m, s, probe, 6, np.random.default_rng(8))
    assert len(a) == 6 and a == b
    x, y = probe(d, np.random.default_rng(9)), probe(s, np.random.default_rng(9))  # the injector's own copy
    assert core.issparse(y.X) and np.array_equal(x.X, y.X.toarray())
    if name == "module":
        assert probe.depth_loss(d) == probe.depth_loss(s)


def test_gate5s_null_and_positive_control_draw_what_the_dense_path_draws():
    d = demo_data()
    s = csr(d)
    a = gates._shuffled_null(metrics.norm_pearson, d, ("Actb", "Gapdh"), 12, np.random.default_rng(10))
    b = gates._shuffled_null(metrics.norm_pearson, s, ("Actb", "Gapdh"), 12, np.random.default_rng(10))
    assert np.array_equal(a, b)
    a = gates._control_response(metrics.norm_pearson, d, ("Actb", "Gapdh"), 2.0, 6, np.random.default_rng(11))
    b = gates._control_response(metrics.norm_pearson, s, ("Actb", "Gapdh"), 2.0, 6, np.random.default_rng(11))
    assert a == b


# --------------------------------------------------------------------------- #
# the same data stored sparse: the same report
# --------------------------------------------------------------------------- #
def _demo_pair(d):
    return run_autopsy(partial(metrics.mi_3bin, gene_a=PAIR[0], gene_b=PAIR[1]), d, group_col="age",
                       groups=("young", "old"), within=["sex"], gene_pair=PAIR, replicate_col="mouse",
                       pair_metric=metrics.mi_3bin, pos_pair=("Actb", "Gapdh"), neg_pair=("Gene0", "Gene1"),
                       signal_test=injected_signal.coupling(*PAIR), prereg=dict(PREREG), log_path="off")


def _typed_pair(d):
    return run_autopsy(partial(metrics.norm_pearson, gene_a=PAIR[0], gene_b=PAIR[1]), d, group_col="age",
                       groups=("young", "old"), gene_pair=PAIR, replicate_col="mouse", composition=["cell_type"],
                       pair_metric=metrics.norm_pearson, signal_test=injected_signal.coupling(*PAIR),
                       prereg=dict(PREREG, judgment_pending=False), log_path="off")


def _set(d):
    return run_autopsy(partial(mean_norm_pearson, genes=SET), d, group_col="age", groups=("young", "old"),
                       replicate_col="mouse", signal_test=injected_signal.coupling_set(SET), prereg=dict(PREREG),
                       log_path="off", specificity=False)


def _module(d):
    return run_autopsy(partial(mean_log_level, genes=SET), d, group_col="age", groups=("young", "old"),
                       replicate_col="mouse", signal_test=injected_signal.module(SET), prereg=dict(PREREG),
                       gene_metric=mean_log_level, log_path="off")


def _global(d):
    return run_autopsy(metrics.spectral_entropy, d, group_col="age", groups=("young", "old"), within=["sex"],
                       replicate_col="mouse", prereg=dict(PREREG), log_path="off")


def _content(d):
    return run_autopsy(log_total_endogenous, d, group_col="age", groups=("young", "old"), replicate_col="mouse",
                       prereg=dict(estimand="content", direction="decrease", sesoi=0.1, sesoi_scale="observed"),
                       log_path="off")


CLAIMS = {
    "a pair with controls and a coupling probe (the demo)": (demo_data, _demo_pair),
    "a pair in its domain of cell types, with the check of specificity": (lambda: typed_mice(n_mice=4), _typed_pair),
    "a set of genes with the set's coupling probe": (lambda: set_data(rho=(0.5, 0.15), seed=2), _set),
    "a set's level with the module probe and the check of specificity": (lambda: set_data(seed=3), _module),
    "a whole-matrix metric (no genes)": (demo_data, _global),
    "a content estimand (spike-in capture)": (lambda: with_ercc(content=(1.0, 0.6), capture=(1.0, 0.5)), _content),
    "non-count input": (lambda: (lambda d: SimpleData(np.log1p(d.X), d.obs, d.var_names))(demo_data()), _demo_pair),
}


@pytest.mark.parametrize("name", list(CLAIMS))
def test_the_same_data_stored_sparse_get_the_same_report(name):
    make, run = CLAIMS[name]
    d = make()
    assert report(run(d)) == report(run(csr(d)))


def test_an_anndata_with_a_sparse_x_gets_the_dense_simpledatas_report():
    anndata = pytest.importorskip("anndata")
    d = demo_data()
    ad = anndata.AnnData(X=sp.csr_matrix(d.X), obs=d.obs.copy(), var=pd.DataFrame(index=list(d.var_names)))
    ad.obs_names = [f"cell_{i}" for i in range(d.n_obs)]
    assert report(_demo_pair(ad)) == report(_demo_pair(d))


# --------------------------------------------------------------------------- #
# no path densifies the whole matrix
# --------------------------------------------------------------------------- #
WIDE = 2500  # wider than any block of columns the engine densifies (chunks of 256, specificity's 2 x 1,000 genes)
# Blocks of every gene the engine may densify: `core.ROW_CHUNK` rows at a time (the sums of non-integer rows), and
# the one documented dense step, GATE 0's null in depth bins (of non-count data, and of the library-size
# perturbation, whose counts are no longer integers), which permutes every gene within one depth bin at a time on
# that bin's dense block to draw the dense path's random numbers (a tenth of the cells by every gene).
ALLOWED = dict(on=False)


def _whole(shape) -> bool:
    return shape[1] > WIDE and shape[0] > core.ROW_CHUNK and not ALLOWED["on"]


class Guarded(sp.csr_matrix):
    """A CSR matrix that refuses to be densified whole (its slices and copies keep the class)."""

    def toarray(self, *a, **k):
        if _whole(self.shape):
            raise AssertionError(f"the whole matrix densified: {self.shape}")
        return super().toarray(*a, **k)


def wide(n_mice=4, cells=100, n_genes=3000, seed=0):
    """Mice nested in age: a coupled pair, a coupled control pair and `n_genes` sparse fillers."""
    rng = np.random.default_rng(seed)
    mu = np.exp(rng.normal(-2.5, 1.0, n_genes))
    rows, blocks = [], []
    for age in ("young", "old"):
        for m in range(n_mice):
            lat, hk = rng.normal(size=(cells, 1)), rng.normal(size=(cells, 1))
            pair = np.exp(1.2 + 0.9 * lat + rng.normal(0, 0.3, (cells, 2)))
            ctrl = np.exp(1.6 + 0.9 * hk + rng.normal(0, 0.2, (cells, 2)))
            fill = mu[None, :] * np.exp(rng.normal(0, 0.3, (cells, 1)))
            blocks.append(rng.poisson(np.hstack([pair, ctrl, fill])).astype(float))
            rows += [dict(age=age, mouse=f"{age}{m}", sex="female")] * cells
    names = list(PAIR) + ["Actb", "Gapdh"] + [f"G{i}" for i in range(n_genes)]
    return SimpleData(np.vstack(blocks), pd.DataFrame(rows), names)


def test_no_path_densifies_the_whole_matrix(monkeypatch):
    d = wide()
    orig = core.as_dense

    def guarded(x):
        if core.issparse(x) and _whole(x.shape):
            raise AssertionError(f"the whole matrix densified: {x.shape}")
        return orig(x)
    for mod in (core, gates, inj, qc, equalize, specificity):
        monkeypatch.setattr(mod, "as_dense", guarded)
    shuffle, bins_seen = gates._shuffle_null_sparse, []

    def one_bin_at_a_time(sd, bins, rng):
        bins_seen.append(int(np.bincount(bins).max()))
        ALLOWED["on"] = True
        try:
            return shuffle(sd, bins, rng)
        finally:
            ALLOWED["on"] = False
    monkeypatch.setattr(gates, "_shuffle_null_sparse", one_bin_at_a_time)
    s = SimpleData(Guarded(sp.csr_matrix(d.X)), d.obs, d.var_names)
    assert isinstance(s.X, Guarded)

    def run(x):
        return run_autopsy(partial(metrics.norm_pearson, gene_a=PAIR[0], gene_b=PAIR[1]), x, group_col="age",
                           groups=("young", "old"), gene_pair=PAIR, replicate_col="mouse",
                           pair_metric=metrics.norm_pearson, pos_pair=("Actb", "Gapdh"), neg_pair=("G0", "G1"),
                           signal_test=injected_signal.coupling(*PAIR), prereg=dict(PREREG), log_path="off")
    assert report(run(s)) == report(run(d))
    assert not bins_seen or max(bins_seen) <= d.n_obs // 5  # where the depth-bin null ran: a bin, not the matrix
