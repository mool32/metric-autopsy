"""Tests of the fourth run's rules (real_pilot_6.py, real_pilot_6_spec.py, real_pilot_6_text.py), on fabricated
records, fabricated texts and random counts only. The workflow real-pilot-6.yml runs them in its dry run, with
the engine of the pull request's head:

    python -m pytest -q -p no:cacheprovider -c /dev/null validation/exploratory/test_real_pilot_6.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import real_pilot_5_spec as S5  # noqa: E402
import real_pilot_6 as R  # noqa: E402
import real_pilot_6_spec as S  # noqa: E402
import real_pilot_6_text as T  # noqa: E402


# --------------------------------------------------------------------------- #
# the claims and the criteria are the plan's
# --------------------------------------------------------------------------- #
def test_the_claims_are_the_third_runs_with_the_new_criteria():
    by5 = {c["id"]: c for c in S5.ACCEPT}
    assert [c["id"] for c in S.ACCEPT] == ["A1", "A2", "A3", "A4"]
    for c in S.ACCEPT:
        same = {k: v for k, v in c.items() if k not in ("criterion", "attempt")}
        assert same == {k: v for k, v in by5[c["id"]].items() if k != "criterion"}
    assert {c["id"]: c["attempt"] for c in S.ACCEPT} == {"A1": 1, "A2": 1, "A3": 2, "A4": 2}
    assert S.CRITERIA["A8"] == "no SUPPORTED; INCONCLUSIVE allowed"
    assert set(S.PARTS) == {"case2", "case1m", "case1", "case3a", "case3b", "case3c"}
    assert "nnz +1298, p 0,016; raw_sum p 0,008" in S.PREDICTION


def test_the_expected_domain_is_the_third_runs_table_by_the_rule_in_each_group():
    """A3_DOMAIN_EXPECTED is read off the third run's logged table (a1) with the rule fixed by the owner."""
    log = (HERE / "real_pilot_5_summary.log").read_text().splitlines()
    i = next(k for k, ln in enumerate(log) if ln.startswith("## (a1) the detections per type"))
    rows = []
    for ln in log[i + 4:]:
        if not ln.startswith("| "):
            break
        cells = [c.strip() for c in ln.strip("|").split("|")]
        rows.append((cells[0], *map(float, cells[2:6])))
    assert len(rows) == 35
    each = sorted(t for t, sy, cy, so, co in rows if min(sy, cy, so, co) >= 0.10)
    assert each == sorted(S.A3_DOMAIN_EXPECTED)


# --------------------------------------------------------------------------- #
# the scores
# --------------------------------------------------------------------------- #
def _tech(stratum, genes, p_genes, counts, p_counts):
    return dict(rows=[dict(cells="all cells", covariate="genes per cell", stratum=stratum, first=1.0, second=0.0,
                           difference=genes, p=p_genes, method="exact", n_relabelings=126, replicates={}),
                      dict(cells="all cells", covariate="counts per cell", stratum=stratum, first=0.0, second=1.0,
                           difference=counts, p=p_counts, method="exact", n_relabelings=126, replicates={})])


def _rec(label, cause="x", g0="PASS", equivalent=None, tech=None, domain=None, effect_status="DETECTED"):
    return dict(label=label, cause=cause, domain=domain if domain is not None else dict(types=["t"]),
                dossier=dict(effect_status=effect_status, dossier_only=False, technical=tech or _tech("all", 1, 0.01, -1, 0.01),
                             gate0=dict(status=g0, between_groups=None if equivalent is None else dict(equivalent=equivalent))))


def test_the_criteria_and_the_prediction_are_scored_as_the_plan_fixes_them():
    blocked = dict(label="NOT SUPPORTED", cause="metric_invalid_gate0", g0="FAIL", equivalent=False)
    recs = {"A1": _rec("SUPPORTED"), "A2": _rec("SUPPORTED"),
            "A3": _rec(**blocked, tech=_tech("all", 1298.0, 0.0159, -1.9e5, 0.0079)),
            "A4": _rec(**blocked, tech=_tech("male", 1298.0, 0.0159, -1.9e5, 0.0079)),
            **{f"A8-3a{i}": _rec("NO DETECTABLE EFFECT") for i in range(1, 9)},
            **{f"A8-3b{i}": _rec("NO DETECTABLE EFFECT") for i in range(1, 9)},
            **{f"A8-3c{i}": _rec("NO DETECTABLE EFFECT") for i in range(1, 5)}}
    recs["A8-3a8"] = _rec("INCONCLUSIVE")
    sc = R.score(recs)
    assert all(sc[k]["met"] for k in ("A1", "A2", "A3", "A4")) and sc["A8"]["met"]
    assert sc["prediction"]["score"] == "holds"
    assert sc["A8"]["changed_from_the_third_run"] == ["A8-3c3"]  # the third run's EFFECT BELOW SESOI
    # change 12 showed the groups equal: GATE 0 does not block, the first part fails
    recs["A4"] = _rec("INCONCLUSIVE", g0="WARN", equivalent=True, tech=_tech("male", 1298.0, 0.0159, -1.9e5, 0.0079))
    sc = R.score(recs)
    assert sc["prediction"]["score"] == "partly holds" and not sc["prediction"]["gate0_blocks"]
    # the covariates of the males do not differ: the second part fails too
    recs["A3"] = _rec(**blocked, tech=_tech("all", 1298.0, 0.2, -1.9e5, 0.0079))
    assert R.score(recs)["prediction"]["score"] == "does not hold"
    # A3's stratum is "all", A4's "male": a row of another stratum does not count
    recs["A3"] = _rec(**blocked, tech=_tech("female", 1298.0, 0.0159, -1.9e5, 0.0079))
    assert R.score(recs)["prediction"]["per_claim"]["A3"]["covariates_differ"] is False
    # A8: one SUPPORTED, or a claim that raised, fails the criterion; a missing claim too
    recs["A8-3b2"] = _rec("SUPPORTED")
    assert R.score(recs)["A8"]["supported"] == ["A8-3b2"] and R.score(recs)["A8"]["met"] is False
    recs["A8-3b2"] = dict(id="A8-3b2", error="ValueError()")
    assert R.score(recs)["A8"]["exceptions"] == ["A8-3b2"] and not R.score(recs)["A8"]["met"]
    del recs["A8-3b2"]
    assert R.score(recs)["A8"]["n"] == 19 and not R.score(recs)["A8"]["met"]
    # A3 without the dossier's covariates, or SUPPORTED, fails its criterion
    recs["A3"] = _rec(**blocked, tech=dict(rows=[]))
    assert R.score(recs)["A3"]["met"] is False
    recs["A3"] = _rec("SUPPORTED")
    assert R.score(recs)["A3"]["met"] is False and R.score(recs)["A3"]["not_supported"] is False


# --------------------------------------------------------------------------- #
# the dossier read from an engine's report
# --------------------------------------------------------------------------- #
def test_the_dossier_is_read_from_a_real_report_of_the_engine():
    from functools import partial

    from metric_autopsy import SimpleData, injected_signal, metrics, run_autopsy

    rng = np.random.default_rng(0)
    rows, blocks = [], []
    for age in ("young", "old"):
        for m in range(4):
            for t in ("A", "B"):
                lat = rng.normal(0, 1, 60)
                pair = np.exp(1.0 + (0.9 * lat[:, None] if t == "A" else 0) + rng.normal(0, 0.3, (60, 2)))
                if t == "B":
                    pair[:, 1] = 0.0
                filler = np.exp(0.4 + rng.normal(0, 0.5, (60, 20)))
                blocks.append(rng.poisson(np.column_stack([pair, filler])).astype(float))
                rows += [dict(age=age, mouse=f"{age}{m}", cell_type=t)] * 60
    d = SimpleData(np.vstack(blocks), pd.DataFrame(rows), ["Smad3", "Col1a1"] + [f"G{i}" for i in range(20)])
    a = run_autopsy(lambda x: -partial(metrics.mi_3bin, gene_a="Smad3", gene_b="Col1a1")(x), d, group_col="age",
                    groups=("young", "old"), replicate_col="mouse", gene_pair=("Smad3", "Col1a1"),
                    composition=["cell_type"], signal_test=injected_signal.coupling("Smad3", "Col1a1"),
                    prereg={"estimand": "composition", "direction": "decrease", "sesoi": 0.1, "judgment_pending": False},
                    log_path="off")
    rep = json.loads(json.dumps(a.to_dict()))
    dz = R.dossier(rep)
    assert set(dz) >= {"dossier_only", "effect_status", "gate0", "unmeasured", "one_group", "technical"}
    assert dz["technical"]["rows"] and dz["domain_rule"].startswith("types where, in each group")
    if a.metric_validity.status == "FAIL":
        assert dz["dossier_only"]
    rec = dict(id="X", claim=dict(attempt=2), verdict=rep["verdict"], cause=rep["cause"], gates=rep["gates"],
               effect={}, domain=None, gate4=None, specificity=None, dossier=dz)
    assert R.claim_lines(rec)[0].startswith("### X (attempt 2):")


# --------------------------------------------------------------------------- #
# what the TMS authors write: the parser, on fabricated texts
# --------------------------------------------------------------------------- #
def test_the_text_check_quotes_matching_sentences_with_their_section_and_writes_no_text(tmp_path):
    assert T.main(["--out", str(tmp_path), "--synthetic"]) == 0
    files = sorted(p.name for p in tmp_path.iterdir())
    assert files == ["real_pilot_6_text.json", "real_pilot_6_text.log"]
    res = json.loads((tmp_path / "real_pilot_6_text.json").read_text())
    epmc = next(r for r in res["sources"] if r["source"] == "tms-epmc")
    texts = [q["text"] for q in epmc["quotes"]]
    assert "SYNTHETIC: the 3-month data were previously published as the Tabula Muris." in texts
    assert not any("about nothing" in t for t in texts)  # a sentence that matches no pattern is not quoted
    assert {q["section"] for q in epmc["quotes"]} == {"Results", "Methods"}
    bio = next(r for r in res["sources"] if r["source"] == "tms-biorxiv")
    assert bio["version"] == "2" and bio["jats_url"].endswith("v2.source.xml")  # the newest version


def test_the_patterns_match_the_origin_and_not_the_atlas_name():
    pats = [re.compile(p, re.I) for p in S.TEXT_PATTERNS]

    def hit(s):
        return any(p.search(s) for p in pats)
    assert hit("data from 3-month-old mice") and hit("the 3m time point") and hit("three-month-old")
    assert hit("as previously published (Tabula Muris Consortium, 2018)")
    assert hit("the Tabula Muris data") and not hit("Tabula Muris Senis is an atlas of 23 tissues")
    assert not hit("at 13, 18, 24 and 30 months the cells")  # the 3-month pattern needs a 3 of its own
    assert not any(p.search("We sequenced 30 million reads") for p in pats)


def test_a_failed_fetch_is_recorded_and_never_raises(monkeypatch):
    def boom(*a, **k):
        raise OSError("no route")
    monkeypatch.setattr(T.urllib.request, "urlopen", boom)
    monkeypatch.setattr(T.time, "sleep", lambda s: None)
    rec = T.check("tms-epmc", "x", "https://example.org/x", synthetic=False)
    assert rec["outcome"].startswith("FAIL OSError") and "quotes" not in rec


@pytest.mark.parametrize("xml, n", [("<p>One. Two.</p><title>S</title><p>Three.</p>", 2),
                                    ("<collection><passage><text>A 3m one.</text></passage></collection>", 1)])
def test_the_blocks_of_jats_and_bioc(xml, n):
    assert len(T.blocks(xml)) == n
