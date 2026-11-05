"""v0.4, the owner's second decisions of 2026-10-19 (validation/plans/v0.4_spec.md, changes 18 and 19):

- change 18: one direction for every difference of the dossier, groups[0] - groups[1], named in its header;
- change 19: the correction's line, with a warning where equal depth widened the groups' gap in genes per cell.

On the synthetic data of `tests/detection_data.py` (mice nested in age, the old group deeper but less complex,
as on preprint §4). These tests were in `tests/test_v04_detection.py` with change 17', reverted by the owner on
2026-11-05; the reference row is again on the uncorrected data."""
from __future__ import annotations

import functools

import pytest

from detection_data import METRICS, PREREG, SCENARIOS, deeper_less_complex
from metric_autopsy import injected_signal, run_autopsy
from metric_autopsy import report as R


@functools.lru_cache(maxsize=None)
def claim(scenario: str, metric: str):
    fn, mu = METRICS[metric]
    d = deeper_less_complex(mu_claim=mu, **SCENARIOS[scenario])
    return run_autopsy(fn, d, group_col="age", groups=("young", "old"), replicate_col="mouse", gene_pair=("A", "B"),
                       signal_test=injected_signal.coupling("A", "B"), prereg=dict(PREREG), log_path="off",
                       specificity=False)


# --------------------------------------------------------------------------- #
# change 18: one direction in every line of the dossier
# --------------------------------------------------------------------------- #
def test_every_difference_of_the_dossier_is_groups0_minus_groups1():
    a = claim("gap, equal coupling", "mi_3bin")
    e = a.effect.detail
    # the reference row's metric is the effect's own uncorrected difference, in the same direction
    assert a.reference["rows"][-1]["metric"]["difference"] == pytest.approx(e["raw_effect"])
    tech = {(r["cells"], r["covariate"]): r for r in a.technical["rows"]}
    gpc = tech[("all cells", "genes per cell")]
    assert gpc["difference"] == pytest.approx(gpc["first"] - gpc["second"]) and gpc["difference"] > 0
    md = a.to_markdown()
    assert "Every difference below is young − old (groups[0] − groups[1]" in md
    assert "young − old |" in md  # the reference table's columns
    assert a.to_dict()["differences"] == "young - old"


def test_the_specificity_line_gives_t_as_groups0_minus_groups1():
    from metric_autopsy.specificity import describe
    line = dict(applicable=True, statistic="the t", data="uncorrected", side="in the claimed direction (increase "
                "from young to old)", groups=["young", "old"], sign=-1.0, t_claim=-2.5, t_claim_signed=2.5,
                in_direction=dict(k=1, n=10, p=0.2, share_ci=(0, 0.4), outcome="undecided"),
                absolute=dict(k=1, n=10, p=0.2, outcome="undecided"),
                nearest=[dict(genes=["X", "Y"], t=-1.0, t_signed=1.0, at_least_as_far=False)],
                shift=dict(median_t=0.5, ci=(0.1, 0.9), ci_level=0.95, n_relabelings=70, p_one_sided=0.02))
    out = describe(line)
    assert "young − old: t = -2.5" in out[0] and "a negative t" in out[0]
    assert "X–Y -1" in out[2] and "median t -0.5" in out[3] and "[-0.9, -0.1]" in out[3]


# --------------------------------------------------------------------------- #
# change 19: the correction's line
# --------------------------------------------------------------------------- #
def test_the_correction_line_warns_where_equal_depth_widens_the_gap_in_detection():
    a = claim("gap, equal coupling", "mi_3bin")
    cap = a.capture
    row = cap["rows"][0]
    assert row["cells"] == "all cells" and row["genes_raw"] > 0
    assert row["widened"] == (abs(row["genes_equal_depth"]) > abs(row["genes_raw"])) and cap["widened"]
    md = a.to_markdown()
    assert "## The correction and the groups' detection" in md and "WARNING: the correction widened" in md
    lines = R.describe_capture(dict(groups=["young", "old"], rows=[
        dict(cells="all cells", genes_raw=0.24, genes_equal_depth=0.29, depth_raw=-0.54, depth_equal_depth=-0.06,
             widened=True),
        dict(cells="the domain's cells", genes_raw=0.17, genes_equal_depth=0.12, depth_raw=-0.6, depth_equal_depth=-0.1,
             widened=False)]))
    assert "+24.0% uncorrected, +29.0% at equal depth" in lines[0] and "WARNING" in lines[0]
    assert "WARNING" not in lines[1]


def test_no_line_of_change_17_prime_is_left_in_the_report():
    # change 17' reverted (the owner, 2026-11-05): no equalized detection, GATE 4 on the data as they are
    a = claim("gap, equal coupling", "mi_3bin")
    assert "detection" not in a.effect.detail and "equal_depth_effect" not in a.effect.detail
    g4 = next(r for r in a.results if r.gate == 4)
    assert "data" not in g4.detail and "equalization" not in g4.detail
    assert "## The pair's detection, equalized" not in a.to_markdown()
