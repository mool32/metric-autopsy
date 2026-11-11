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
