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
        lines=2, kappa_decision=0.0, kappa_code=0.0, agreement_decision=0.5, agreement_code=0.5, verdicts_agreed=1,
        to_owner=1, claims=[claim])))
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
    assert "the decisions: agreement 0.500, kappa 0.000" in md
    # every line NO_VERDICT in both passes (S4, J31): the kappa is undefined, said so
    v["merge"].update(kappa_decision=float("nan"), kappa_code=float("nan"), agreement_decision=1.0, agreement_code=1.0)
    md = L.markdown(dict(L.assemble(tmp_path), verdicts=v))
    assert "the decisions: agreement 1.000, kappa undefined (one category in both passes)" in md


def _report(raw=True, rule4=True, missing=(), source="balanced"):
    return dict(counts=dict(nonnegative=True, integers=raw, raw=raw, max=10.0), cells=10, cells_per_group={"A": 5, "B": 5},
                donors=dict(design="nested", donors={"A": 4, "B": 4}, composite=[],
                            rule4=dict(design="nested", relabelings=70, min_p=1 / 70 if rule4 else 0.2, passes=rule4)),
                sex=dict(both_sexes=False, fails=[]),
                genes={str(i): dict(present=["G1"], missing=list(missing), passes=not missing) for i in range(1, 7)},
                source=dict(verdict=source, sources=["s"]), files=[dict(file="f", bytes=1, sha256="0" * 64, url="u")])


def test_a_claim_on_several_datasets_and_the_branches_checked(tmp_path):
    """R02175's claims are checked on three datasets (its specs): one failing the check of origin makes them fail on
    some, not every, dataset; a second attempt supersedes its first; a branch's copy is checked against SHA256SUMS."""
    import hashlib
    import json
    res = tmp_path / "res"
    for branch, name, rep in (("level-b-data-3", "R02175", _report()), ("level-b-data-3", "R02175_d2", _report(missing=["X"])),
                              ("level-b-data-4", "R02175_d2", _report(raw=False))):
        (res / branch / name).mkdir(parents=True, exist_ok=True)
        (res / branch / name / "report.json").write_text(json.dumps(rep))
    for branch in ("level-b-data-3", "level-b-data-4"):
        files = sorted(p for p in (res / branch).rglob("report.json"))
        (res / branch / "SHA256SUMS").write_text("".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  "
                                                         f"{p.relative_to(res / branch)}\n" for p in files))
    heads = tmp_path / "heads.tsv"
    heads.write_text("results/level-b-data-3\tabc\nresults/level-b-data-4\tdef\n")
    L1 = L.assemble(res, heads=heads, commit="123")
    assert L1["origin_fails"]["every"] == [] and L1["origin_fails"]["some"] == [f"R02175:{i}" for i in range(1, 7)]
    assert L1["genes_missing"]["every"] == [] and L1["genes_missing"]["some"] == []  # the first attempt's misses superseded
    p = next(x for x in L1["papers"] if x["rid"] == "R02175")
    assert p["data"]["R02175_d2"]["report"]["_supersedes"] == "level-b-data-3"
    assert L1["provenance"]["branches"]["level-b-data-4"] == dict(L1["provenance"]["branches"]["level-b-data-4"], commit="def", bad=[])
    assert "R00765:1" in L1["data_not_checked"]["genes_not_listed"]
    md = L.markdown(L1)
    assert "a second attempt (the first on level-b-data-3)" in md and "origin FAILS: counts not raw integers" in md
    (res / "level-b-data-4" / "R02175_d2" / "report.json").write_text("{}")
    assert L.branch_check(res / "level-b-data-4")["bad"] == ["R02175_d2/report.json"]


def test_the_owner_s_part_from_the_repository_s_files():
    fo = L.ft_to_owner()
    assert len(fo["without_fulltext"]) == 35 and len(fo["decided_differently"]) == 7
    r = next(x for x in fo["decided_differently"] if x["rid"] == "R01215")
    assert r["passes"]["2"].startswith("KEEP") and r["claims"]["2"][0].startswith("claim 1: SET: Female-1 module")
    j = {e["entry"]: e for e in L.journal()}
    assert j["J3"]["deviation"] and j["J31"]["deviation"] and not j["J9"]["deviation"]  # J9: "No deviation."
    one = {c["claim"]: c for c in L.one_pass_claims(L._json(L.A1 / "kept_claims.json"))}
    assert one["R00653:1"]["same_genes_as"] == ["3 (passes 2)"] and one["R00765:2"]["same_genes_as"] == []
