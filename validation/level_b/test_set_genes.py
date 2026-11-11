"""set_genes.py: the sets' genes from a supplement's tables (fabricated archives, nothing fetched)."""
from __future__ import annotations

import gzip
import io
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import set_genes as G  # noqa: E402
from test_data_run import _xlsx  # noqa: E402


def _archive(tmp_path: Path, rid: str, members: dict) -> Path:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        for name, data in members.items():
            z.writestr(name, data)
    raw = tmp_path / "raw" / "setgenes"
    raw.mkdir(parents=True, exist_ok=True)
    (raw / f"{rid}-supplementary.zip.gz").write_bytes(gzip.compress(buf.getvalue(), mtime=0))
    (tmp_path / "setgenes_index.tsv").write_text(
        "rid\tpmcid\tstatus\tfile\tbytes\tsha256\tnote\n"
        f"{rid}\tPMC1\t200\traw/setgenes/{rid}-supplementary.zip.gz\t1\t-\t-\n")
    return tmp_path


def test_the_dry_run(tmp_path):
    res = G.dry(tmp_path / "dry")
    assert res["sets"][0]["genes"] == ["GENE1", "GENE2", "GENE3"]


def test_an_excel_sheet_with_a_title_row_a_filter_and_cells_holding_several_genes(tmp_path):
    _xlsx(tmp_path / "s.xlsx", {"Notes": [["read me"]],
                               "Networks": [["Supplementary Data 2: the networks"], ["network", "genes"],
                                            ["Mitochondria", "ATP5F1A; NDUFA4"], ["Glycolysis", "PKM; ENO1"],
                                            ["Mitochondria", "COX7C; NDUFA4"]]})
    files = _archive(tmp_path, "R02175", {"41467_2024_1_MOESM4_ESM.xlsx": (tmp_path / "s.xlsx").read_bytes()})
    kept = dict(papers=[dict(rid="R02175", claims=[dict(id="1", kind="set", genes=[], set_name="Mitochondria",
                                                        set_size=3)])])
    spec = dict(rid="R02175", claim="1", member="41467_2024_1_MOESM4_ESM.xlsx", sheet="Networks", header_row=1,
                column="genes", rows=dict(network=["Mitochondria"]), split=";")
    res = G.extract([spec], files, kept)
    assert res[0]["genes"] == ["ATP5F1A", "NDUFA4", "COX7C"] and res[0]["n"] == 3 and res[0]["stated"] == 3
    G.fill(kept, res)
    assert kept["papers"][0]["claims"][0]["genes"] == ["ATP5F1A", "NDUFA4", "COX7C"]
    text = G.show(files, tmp_path / "show.md")
    assert "Networks: 5 rows" in text and "Notes: 1 rows" in text


def test_a_spec_for_a_claim_that_is_not_a_kept_set_stops(tmp_path):
    files = _archive(tmp_path, "R00001", {"t.tsv": "gene\nA\n"})
    kept = dict(papers=[dict(rid="R00001", claims=[dict(id="1", kind="pair", genes=["A", "B"])])])
    try:
        G.extract([dict(rid="R00001", claim="1", member="t.tsv", column="gene")], files, kept)
    except SystemExit as e:
        assert "not a kept set's claim" in str(e)
    else:
        raise AssertionError("no stop")


def test_the_extra_files_list_and_a_gmt_set(tmp_path):
    (tmp_path / "x.tsv").write_text("rid\turl\nR01795\thttps://example.org/h.all.v1.Hs.symbols.gmt\n")
    assert G.read_extra(tmp_path / "x.tsv") == [("R01795", "https://example.org/h.all.v1.Hs.symbols.gmt")]
    gmt = b"HALLMARK_A\thttp://x\tG1\tG2\nHALLMARK_B\thttp://x\tG3\n"
    assert G.gmt_set(gmt, "HALLMARK_B") == ["G3"]
    try:
        G.gmt_set(gmt, "HALLMARK_C")
    except SystemExit as e:
        assert "no set" in str(e)
    else:
        raise AssertionError("no stop")


def test_a_set_listed_in_the_paper_s_text_and_a_set_not_listed(tmp_path):
    import hashlib
    texts = tmp_path / "texts"
    texts.mkdir()
    body = "Results.\nWe identified a set of 3 genes (COX5A, COX5B,\nXBP1) that exhibited higher centrality.\n"
    (texts / "R00765.txt").write_text(body)
    sha = hashlib.sha256(body.encode()).hexdigest()
    files = _archive(tmp_path, "R00765", {"x.tsv": "gene\nA\n"})
    kept = dict(papers=[dict(rid="R00765", claims=[dict(id="1", kind="set", genes=[], set_name="Dozer-A3 module", set_size=57),
                                                   dict(id="2", kind="set", genes=[], set_name="neurodegeneration", set_size=3)])])
    specs = [dict(rid="R00765", claim="1", not_listed="shown only in figures"),
             dict(rid="R00765", claim="2", text="R00765.txt", sha256=sha, quote="a set of 3 genes (COX5A, COX5B, XBP1) that",
                  list_from="(", list_to=")", split=",")]
    res = G.extract(specs, files, kept, texts)
    assert res[0]["genes"] == [] and res[0]["note"] == "not listed: shown only in figures"
    assert res[1]["genes"] == ["COX5A", "COX5B", "XBP1"] and res[1]["n"] == res[1]["stated"] == 3
    G.fill(kept, res)
    assert kept["papers"][0]["claims"][0]["genes_from"] == "not listed: shown only in figures"
    assert kept["papers"][0]["claims"][1]["genes_from"] == "set_genes.json (the text R00765.txt)"
    for bad, why in ((dict(specs[1], sha256="0" * 64), "not the text S3 read"),
                     (dict(specs[1], quote="a set of 4 genes"), "the quote is not in")):
        try:
            G.extract([bad], files, kept, texts)
        except SystemExit as e:
            assert why in str(e)
        else:
            raise AssertionError("no stop")


def test_another_version_s_supplement_an_archive_inside_it_and_several_columns(tmp_path):
    """2026-11-11: a record's journal version fetched by a list; its genes in an archive inside the archive."""
    import search as S
    first = _archive(tmp_path / "one", "R00653", {"media-1.pdf": b"%PDF"})
    (tmp_path / "v.tsv").write_text("rid\tpmcid\nR00653\tPMC200\n")
    out = S.Out(tmp_path / "two")
    rows = G.fetch(G.Fabricated(out.log), out, dict(papers=[]), G.read_versions(tmp_path / "v.tsv"))
    assert rows[0]["file"] == "raw/setgenes/R00653-PMC200-supplementary.zip.gz"
    # an archive inside the journal version's archive, as Nature's Supplementary Data come
    _xlsx(tmp_path / "d.xlsx", {"Modules": [["module", "genes"], ["IFN", "ISG15"], ["IFN", "MX1"], ["AP", "HLA-A"]]})
    inner = io.BytesIO()
    with zipfile.ZipFile(inner, "w") as z:
        z.writestr("Supplementary Data 3.xlsx", (tmp_path / "d.xlsx").read_bytes())
    outer = io.BytesIO()
    with zipfile.ZipFile(outer, "w") as z:
        z.writestr("41467_MOESM4_ESM.zip", inner.getvalue())
    (tmp_path / "two" / rows[0]["file"]).write_bytes(gzip.compress(outer.getvalue(), mtime=0))
    arch = G.archives([first, tmp_path / "two"])
    assert set(arch) == {"R00653", "R00653:PMC1", "R00653:PMC200"} and arch["R00653"] == arch["R00653:PMC1"]
    text = G.show([first, tmp_path / "two"], tmp_path / "show.md")
    assert "`41467_MOESM4_ESM.zip!/Supplementary Data 3.xlsx`" in text and "Modules: 4 rows" in text
    kept = dict(papers=[dict(rid="R00653", claims=[dict(id="1", kind="set", genes=[], set_name="IFN module", set_size=None)])])
    spec = dict(rid="R00653", claim="1", archive="R00653:PMC200", member="41467_MOESM4_ESM.zip!/Supplementary Data 3.xlsx",
                sheet="Modules", column="genes", rows=dict(module=["IFN"]))
    assert G.extract([spec], [first, tmp_path / "two"], kept)[0]["genes"] == ["ISG15", "MX1"]
    # several columns: their union, in the spec's order
    _xlsx(tmp_path / "m.xlsx", {"S": [["title"], ["(5) Microfilaments", "(9) Microtubuli"], ["CAPN2", "AP1M2"], ["MARCKS", "CAPN2"]]})
    files = _archive(tmp_path / "three", "R02175", {"m.xlsx": (tmp_path / "m.xlsx").read_bytes()})
    kept2 = dict(papers=[dict(rid="R02175", claims=[dict(id="6", kind="set", genes=[], set_name="Microfilaments and Microtubuli")])])
    spec2 = dict(rid="R02175", claim="6", member="m.xlsx", sheet="S", header_row=1, column=["(5) Microfilaments", "(9) Microtubuli"])
    assert G.extract([spec2], files, kept2)[0]["genes"] == ["CAPN2", "MARCKS", "AP1M2"]


def test_show_says_which_members_it_cannot_read_and_goes_on(tmp_path):
    """2026-11-11 (J26): a member named .xlsx that is no workbook, and a .zip that is no archive, stopped `show`."""
    _xlsx(tmp_path / "ok.xlsx", {"S": [["gene"], ["A"]]})
    files = _archive(tmp_path, "R00653", {"bad.xlsx": b"not a workbook", "bad.zip": b"not an archive",
                                          "ok.xlsx": (tmp_path / "ok.xlsx").read_bytes()})
    text = G.show(files, tmp_path / "show.md")
    assert "### `bad.xlsx`: not read (BadZipFile" in text and "not read as an archive (BadZipFile" in text
    assert "### `ok.xlsx` / S: 2 rows" in text
