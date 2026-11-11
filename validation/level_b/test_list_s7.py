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
