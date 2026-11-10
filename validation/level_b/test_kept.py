"""kept.py: the kept papers' claims matched across the two passes (fabricated decisions, nothing read)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import kept as K  # noqa: E402
from screening import FT_HEADER  # noqa: E402


def test_pair_candidates_and_keys():
    assert K.pair_candidates("SRGN, RAB31") == {frozenset({"SRGN", "RAB31"})}
    assert K.pair_candidates("SRGN-RAB31") == {frozenset({"SRGN", "RAB31"})}
    assert K.pair_candidates("ACE2 and TMPRSS2") == {frozenset({"ACE2", "TMPRSS2"})}
    assert K.pair_candidates("SERPINE2; TFF3") == {frozenset({"SERPINE2", "TFF3"})}
    assert K.pair_candidates("CSNK2B, HLA-C") == {frozenset({"CSNK2B", "HLA-C"})}
    assert K.pair_candidates("Ace2 – Tmprss2") == {frozenset({"ACE2", "TMPRSS2"})}
    assert K.pair_candidates("HLA-C-B2M") == {frozenset({"HLA", "C-B2M"}), frozenset({"HLA-C", "B2M"})}
    assert K.keys("SET: M3 module, 87 genes, Table S4") == {("set", "m3 module")}
    assert K.set_parts("SET: Mitochondria, 1,204 genes, Supplementary Data 2") == dict(
        set_name="Mitochondria", set_size=1204, set_where="Supplementary Data 2")


def _line(rid, claim, genes):
    row = dict.fromkeys(FT_HEADER, "x")
    row.update(rid=rid, decision="CLAIM", code="-", claim=str(claim), genes=genes, rule4="yes", note="-")
    return "\t".join(row[h] for h in FT_HEADER)


def _decision(rid, decision, code="-"):
    row = dict.fromkeys(FT_HEADER, "-")
    row.update(rid=rid, decision=decision, code=code, claim="0", note="a note")
    return "\t".join(row[h] for h in FT_HEADER)


def _pass(d: Path, rows: dict):
    d.mkdir(parents=True)
    rids = list(rows)
    (d / "batch_01.list").write_text("".join(f"{r}\t/texts/{r}.txt\n" for r in rids))
    lines = ["\t".join(FT_HEADER)]
    for r, (decision, claims) in rows.items():
        lines.append(_decision(r, decision, "F1" if decision == "EXCLUDE" else "-"))
        lines += [_line(r, i, g) for i, g in enumerate(claims, start=1)]
    (d / "batch_01.tsv").write_text("\n".join(lines) + "\n")
    (d / "without_fulltext.txt").write_text("")


def test_kept_claims_matches_both_passes_and_keeps_the_old_papers(tmp_path):
    _pass(tmp_path / "p1", {"R00010": ("KEEP", ["SRGN, RAB31", "CXCL16, HEXB", "SET: M3 module, 87 genes, Table S4"]),
                            "R00020": ("KEEP", ["HLA-C-B2M"]),
                            "R00030": ("EXCLUDE", [])})
    _pass(tmp_path / "p2", {"R00010": ("KEEP", ["CXCL16-HEXB", "SET: M3 Module, 87 genes, Table S4", "A1, B1"]),
                            "R00020": ("KEEP", ["HLA-C, B2M"]),
                            "R00030": ("KEEP", ["X, Y"])})
    _pass(tmp_path / "o1", {"R00001": ("KEEP", ["ACE2, TMPRSS2"]), "R00002": ("EXCLUDE", [])})
    _pass(tmp_path / "o2", {"R00001": ("KEEP", ["ACE2-TMPRSS2"]), "R00002": ("KEEP", ["P, Q"])})
    (tmp_path / "ids.tsv").write_text("rid\tpmid\tpmcid\tdoi\tepmc_ids\nR00010\t1\tPMC1\t10.1/a\tMED:1,PPR:PPR9\n"
                                      "R00001\t2\t\t10.1/b\tMED:2\n")
    (tmp_path / "all.tsv").write_text("\t".join(["rid"] + ["c"] * 11) + "\n" + "\t".join(["R00010"] + ["-"] * 10 + ["A title"]) + "\n")
    res = K.kept_claims(tmp_path / "p1", tmp_path / "p2", tmp_path / "o1", tmp_path / "o2", ["R00001"],
                        [tmp_path / "ids.tsv"], [tmp_path / "all.tsv"])
    assert [p["rid"] for p in res["papers"]] == ["R00001", "R00010", "R00020"]
    old, a, b = res["papers"]
    assert old["read"] == "2026-10-19" and old["claims"][0]["genes"] == ["ACE2", "TMPRSS2"]
    assert old["claims"][0]["passes"] == [1, 2]
    assert a["pmcid"] == "PMC1" and a["epmc_ids"] == ["MED:1", "PPR:PPR9"] and a["title"] == "A title"
    kinds = [(c["kind"], c["genes"], c.get("set_name"), c["passes"]) for c in a["claims"]]
    assert kinds == [("pair", ["CXCL16", "HEXB"], None, [1, 2]), ("set", [], "M3 module", [1, 2]),
                     ("pair", ["SRGN", "RAB31"], None, [1]), ("pair", ["A1", "B1"], None, [2])]
    assert [c["id"] for c in a["claims"]] == ["1", "2", "3", "4"]
    assert b["claims"][0]["genes"] == ["HLA-C", "B2M"] and not b["claims"][0]["ambiguous"]
    c = res["counts"]
    assert c["papers"] == 3 and c["claims"] == 6 and c["sets"] == 1 and c["one_pass"] == 2 and c["ambiguous"] == 0
    assert not res["drawn"] and res["fewer_than_target"]


def test_an_unsettled_pair_is_marked_ambiguous(tmp_path):
    claims = K.paper_claims([dict(genes="HLA-C-B2M")], [])
    assert claims[0]["ambiguous"] and claims[0]["genes"] == [] and claims[0]["passes"] == [1]


def test_main_writes_the_json(tmp_path):
    for d in ("p1", "p2", "o1", "o2"):
        _pass(tmp_path / d, {"R00001": ("KEEP", ["G1, G2"])})
    (tmp_path / "ids.tsv").write_text("rid\tpmid\tpmcid\tdoi\tepmc_ids\n")
    (tmp_path / "all.tsv").write_text("rid\n")
    out = tmp_path / "kept.json"
    assert K.main(["--ft1", str(tmp_path / "p1"), "--ft2", str(tmp_path / "p2"), "--old-ft1", str(tmp_path / "o1"),
                   "--old-ft2", str(tmp_path / "o2"), "--old-kept", "R00001", "--ids", str(tmp_path / "ids.tsv"),
                   "--titles", str(tmp_path / "all.tsv"), "--out", str(out)]) == 0
    res = json.loads(out.read_text())
    assert res["counts"]["papers"] == 1 and res["papers"][0]["read"] == "amendment 1"
