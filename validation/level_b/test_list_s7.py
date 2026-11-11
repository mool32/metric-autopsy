"""list_s7.py: 12.3's counts read from a dataset's report (fabricated reports)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import list_s7 as L  # noqa: E402


def _rep(raw=True, fails=(), source="balanced", passes=True):
    return dict(counts=dict(raw=raw), sex=dict(fails=list(fails)), source=dict(verdict=source, sources=["a", "b"]),
                donors=dict(design="nested", donors={"A": 4, "B": 4}, rule4=dict(passes=passes)))


def test_the_check_of_origin_and_the_test_over_donors():
    assert not L.origin(_rep())["fails"]
    assert L.origin(_rep(raw=False))["why"] == ["counts not raw integers"]
    assert L.origin(_rep(fails=["d1"]))["fails"] and L.origin(_rep(source="fails"))["fails"]
    assert not L.origin(_rep(source="partly"))["fails"]  # partly is listed with its shares, not a failure
    assert L.donors_cannot_run(_rep()) is None and "rule 4 fails" in L.donors_cannot_run(_rep(passes=False))


def test_the_list_on_the_repository_s_files(tmp_path):
    res = L.assemble(tmp_path)  # no results branches: every paper's data not checked yet
    assert res["kept_counts"]["papers"] == len(res["papers"]) and res["controls"] is None
    assert len(res["data_not_checked"]["papers"]) == len(res["papers"])
    assert "## 7. 12.3's two counts so far" in L.markdown(res)


def test_the_verdicts_part_from_the_passes_comparison(tmp_path):
    import json
    root = tmp_path / "root"
    (root / "verdicts" / "merge").mkdir(parents=True)
    (root / "verdicts" / "batches").mkdir(parents=True)
    assert L.verdicts(root) is None
    claim = dict(claim="R00189:1", kind="pair", records=3, by_citations=2, by_terms=1, kept=2, coded=2,
                 verdicts_agreed=1, agreed_codes={"VERDICT/causal/holds": 1}, to_owner=1)
    (root / "verdicts" / "merge" / "merge.json").write_text(json.dumps(dict(
        lines=2, kappa_decision=0.0, kappa_code=0.0, verdicts_agreed=1, to_owner=1, claims=[claim])))
    (root / "verdicts" / "merge" / "verdicts.tsv").write_text(
        "claim\trecord\tkind\tcode\tquote_1\twhere_1\tquote_2\twhere_2\n"
        "R00189:1\tMED:1\tcausal\tholds\tknockdown lowered it\tFig 2\tthe knockdown lowered it\tFig 2B\n")
    (root / "verdicts" / "merge" / "to_owner.tsv").write_text("claim\trecord\tdecision_1\n" "R00189:1\tMED:2\tVERDICT\n")
    (root / "verdicts" / "batches" / "manifest.json").write_text(json.dumps(dict(batches=["batch_01"])))
    v = L.verdicts(root)
    assert v["merge"]["verdicts_agreed"] == 1 and v["to_owner"][0]["record"] == "MED:2"
    md = L.markdown(dict(L.assemble(tmp_path), verdicts=v))
    assert "R00189:1 (pair): found 3 (by citations 2, by the terms 1), kept by the filter 2, coded 2" in md
    assert '  - MED:1: causal, holds — "knockdown lowered it" (Fig 2); pass 2: "the knockdown lowered it" (Fig 2B)' in md
