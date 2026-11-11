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


def test_aliases_and_descriptions_in_parentheses_and_msigdb_names(tmp_path):
    """J21: a gene's alias and a set's description in parentheses are dropped; an MSigDB identifier keys a set."""
    assert K.pair_candidates("TMBIM6, HSPA5 (BIP)") == {frozenset({"TMBIM6", "HSPA5"})}
    assert K.keys("SET: Mitochondria NDCG (mitochondrial electron transport chain), gene count not given, Data 2") == \
        K.keys("SET: Mitochondria NDCG, gene count not given in text, Supplementary Data 2")
    assert K.keys("SET: HALLMARK_INTERFERON_GAMMA_RESPONSE weighted network, hub IFITM3, MSigDB") == \
        {("set", "hallmark_interferon_gamma_response")}
    assert K.keys("SET: IFITM3 with its partners in the IFN-gamma gene set (HALLMARK_INTERFERON_GAMMA_RESPONSE), 118 genes") == \
        {("set", "hallmark_interferon_gamma_response")}
    assert K.set_parts("SET: UPR NDCG, gene count not given in text (671 genes over the 11 NDCGs), Data 2")["set_size"] is None
    claims = K.paper_claims([dict(genes="TMBIM6, HSPA5")], [dict(genes="TMBIM6, HSPA5 (BIP)")])
    assert len(claims) == 1 and claims[0]["passes"] == [1, 2] and claims[0]["genes"] == ["TMBIM6", "HSPA5"]
    assert K.set_parts("SET: Proliferation NDCG (cell proliferation; includes PAX6, ISL1), gene count not given, Data 2")[
        "set_name"] == "Proliferation NDCG (cell proliferation; includes PAX6, ISL1)"
    assert K.keys("SET: Proliferation NDCG (cell proliferation; includes PAX6, ISL1), n, D") == K.keys("SET: Proliferation NDCG, n, D")


def test_ft_record_counts_each_session_from_its_file(tmp_path):
    import ft_record
    _pass(tmp_path / "ft" / "pass1", {"R00010": ("KEEP", ["G1, G2", "G3, G4"]), "R00020": ("EXCLUDE", [])})
    _pass(tmp_path / "ft" / "pass2", {"R00010": ("EXCLUDE", []), "R00020": ("EXCLUDE", [])})
    assert ft_record.main(["--ft", str(tmp_path / "ft"), "--out", str(tmp_path / "s.tsv")]) == 0
    lines = (tmp_path / "s.tsv").read_text().splitlines()
    assert lines[1] == "1\t01\t2\t1\t0\t1\tF1 1\t2" and lines[2] == "2\t01\t2\t0\t0\t2\tF1 2\t0"
