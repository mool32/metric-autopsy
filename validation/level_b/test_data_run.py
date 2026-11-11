"""Level B, step S5 (`plan_s5.md`): data_run.py's loaders and checks on fabricated deposits (nothing fetched)."""
from __future__ import annotations

import gzip
import io
import json
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import data_run as R  # noqa: E402

GENES = ["XIST", "RPS4Y1", "DDX3Y", "ACE2", "TMPRSS2", "G1", "G2"]


def _cells(rng, n, female):
    X = rng.poisson(2.0, size=(n, len(GENES))).astype(float)
    X[:, 0] = rng.poisson(3.0, n) if female else 0
    X[:, 1:3] = 0 if female else rng.poisson(2.0, (n, 2))
    return X


def _mtx(M) -> bytes:  # genes x cells, as 10x writes it
    r, c = np.nonzero(M)
    s = "%%MatrixMarket matrix coordinate integer general\n%\n" + f"{M.shape[0]} {M.shape[1]} {len(r)}\n"
    s += "".join(f"{i + 1} {j + 1} {int(M[i, j])}\n" for i, j in zip(r, c))
    return gzip.compress(s.encode(), mtime=0)


def _series(path: Path, gsms, titles, chars: dict):
    lines = ["!Series_title\t\"FABRICATED\"", "!Sample_title\t" + "\t".join(f'"{t}"' for t in titles),
             "!Sample_geo_accession\t" + "\t".join(f'"{g}"' for g in gsms),
             "!Sample_source_name_ch1\t" + "\t".join('"tissue"' for _ in gsms)]
    for k, vals in chars.items():
        lines.append("!Sample_characteristics_ch1\t" + "\t".join(f'"{k}: {v}"' for v in vals))
    path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode(), mtime=0))


def test_csv_counts_nested_groups_and_rule4(tmp_path):
    rng = np.random.default_rng(0)
    rows, meta = [], []
    for k in range(8):  # 4 tumour patients, 4 normal ones, 30 cells each
        X = _cells(rng, 30, female=True)
        rows.append(X)
        meta += [dict(cell=f"P{k}_{i}", patient=f"P{k}", site="PT" if k < 4 else "FT", type="Epithelial") for i in range(30)]
    X = np.vstack(rows)
    meta = pd.DataFrame(meta)
    pd.DataFrame(X.T, index=GENES, columns=meta["cell"]).to_csv(tmp_path / "counts.csv.gz")
    meta.set_index("cell").to_csv(tmp_path / "meta.csv.gz")
    spec = dict(rid="R00001", organism="homo_sapiens", loader="csv_counts",
                files=dict(counts="https://example.org/counts.csv.gz", meta="https://example.org/meta.csv.gz"),
                groups=dict(values=["tumour", "normal"], rules=[dict(when=dict(site=["PT"]), value="tumour"),
                                                               dict(when=dict(site=["FT"]), value="normal")]),
                donor="patient", source=None, sex_label=None, claims=[dict(id="1", genes=["ACE2", "TMPRSS2"])])
    (tmp_path / "spec.json").write_text(json.dumps(spec))
    assert R.main(["--spec", str(tmp_path / "spec.json"), "--out", str(tmp_path / "out"), "--files", str(tmp_path)]) == 0
    res = json.loads((tmp_path / "out" / "report.json").read_text())
    assert res["counts"]["raw"] and res["cells"] == 240
    assert res["donors"]["design"] == "nested" and res["donors"]["donors"] == {"tumour": 4, "normal": 4}
    assert res["donors"]["rule4"]["passes"] and res["genes"]["1"]["passes"] and not res["sex"]["both_sexes"]
    assert res["source"]["verdict"] == "no source named"


def test_tenx_samples_with_the_series_matrix_paired_and_a_source(tmp_path):
    rng = np.random.default_rng(1)
    gsms, titles, status, patient, chem = [], [], [], [], []
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as tar:
        for k in range(6):
            for s, st in enumerate(("Involved", "Uninvolved")):
                gsm = f"GSM{100 + 2 * k + s}"
                gsms.append(gsm)
                titles.append(f"Ileal {st} {k}")
                status.append(st)
                patient.append(str(k))
                chem.append("V1" if k == 0 else "V2")
                X = _cells(rng, 20, female=k % 2 == 0)
                files = {f"{gsm}_s_barcodes.tsv.gz": gzip.compress("".join(f"B{i}-1\n" for i in range(20)).encode()),
                         f"{gsm}_s_genes.tsv.gz": gzip.compress("".join(f"E{i}\t{g}\n" for i, g in enumerate(GENES)).encode()),
                         f"{gsm}_s_matrix.mtx.gz": _mtx(X.T)}
                for name, data in files.items():
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    tar.addfile(info, io.BytesIO(data))
    (tmp_path / "raw.tar").write_bytes(buf.getvalue())
    _series(tmp_path / "series.txt.gz", gsms, titles, {"status": status, "patient id": patient, "10x chemistry": chem})
    spec = dict(rid="R00002", organism="homo_sapiens", loader="tenx_samples",
                files=dict(raw="https://example.org/raw.tar", series="https://example.org/series.txt.gz"),
                groups=dict(column="status", values=["Involved", "Uninvolved"]), donor="patient id",
                source="10x chemistry", sex_label=None, claims=[dict(id="1", genes=["ACE2", "XIAP"])])
    (tmp_path / "spec.json").write_text(json.dumps(spec))
    assert R.main(["--spec", str(tmp_path / "spec.json"), "--out", str(tmp_path / "out"), "--files", str(tmp_path)]) == 0
    res = json.loads((tmp_path / "out" / "report.json").read_text())
    assert res["donors"]["design"] == "paired" and res["donors"]["rule4"]["passes"]
    assert res["genes"]["1"]["missing"] == ["XIAP"] and res["source"]["verdict"] == "balanced"
    assert res["sex"]["both_sexes"] and res["sex"]["fails"] == []


def test_a_large_matrix_in_a_stream_with_libraries_from_the_series(tmp_path):
    rng = np.random.default_rng(2)
    libs = [("pPDCN100DAPIA", "Ctrl", "Female", "100"), ("pPDCN101DAPIA", "Ctrl", "Male", "101"),
            ("pPD102PosA", "Disease", "Female", "102"), ("pPD103PosA", "Disease", "Male", "103")]
    blocks, bcd = [], []
    for lib, _dis, sex, _d in libs:
        blocks.append(_cells(rng, 25, female=sex == "Female"))
        bcd += [f"{lib}_AAAC{i:04d}-1" for i in range(25)]
    X = np.vstack(blocks)
    (tmp_path / "matrix.mtx.gz").write_bytes(_mtx(X.T))
    (tmp_path / "features.tsv.gz").write_bytes(gzip.compress("".join(f"E{i}\t{g}\n" for i, g in enumerate(GENES)).encode()))
    (tmp_path / "barcodes.tsv.gz").write_bytes(gzip.compress(("\n".join(bcd) + "\n").encode()))
    _series(tmp_path / "s1.txt.gz", ["GSM1", "GSM2"], [libs[0][0], libs[2][0]],
            {"disease": ["Ctrl", "Disease"], "sex": ["Female", "Female"]})
    _series(tmp_path / "s2.txt.gz", ["GSM3", "GSM4"], [libs[1][0], libs[3][0]],
            {"disease": ["Ctrl", "Disease"], "sex": ["Male", "Male"]})
    spec = dict(rid="R00003", organism="homo_sapiens", loader="mtx_libraries",
                files=dict(matrix="https://example.org/matrix.mtx.gz", features="https://example.org/features.tsv.gz",
                           barcodes="https://example.org/barcodes.tsv.gz",
                           series=["https://example.org/s1.txt.gz", "https://example.org/s2.txt.gz"]),
                groups=dict(column="disease", values=["Disease", "Ctrl"]),
                donor=dict(column="library", regex=r"(\d+)"), source=None, sex_label="sex",
                claims=[dict(id="1", genes=["ACE2", "TMPRSS2"])])
    (tmp_path / "spec.json").write_text(json.dumps(spec))
    assert R.main(["--spec", str(tmp_path / "spec.json"), "--out", str(tmp_path / "out"), "--files", str(tmp_path)]) == 0
    res = json.loads((tmp_path / "out" / "report.json").read_text())
    assert res["counts"]["raw"] and res["cells"] == 100
    assert res["donors"]["design"] == "nested" and res["donors"]["donors"] == {"Disease": 2, "Ctrl": 2}
    assert not res["donors"]["rule4"]["passes"] and res["sex"]["both_sexes"] and res["sex"]["fails"] == []
    assert "2. Donors: nested" in (tmp_path / "out" / "report.md").read_text()


def _xlsx(path: Path, sheets: dict):
    """A minimal Excel workbook (numbers as numbers, text as shared strings), written with the standard library."""
    import zipfile
    from xml.sax.saxutils import escape
    M = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    Rn = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    P = "http://schemas.openxmlformats.org/package/2006/relationships"
    shared, idx, xml = [], {}, {}

    def col(j):
        s, j = "", j + 1
        while j:
            j, r = divmod(j - 1, 26)
            s = chr(65 + r) + s
        return s

    for k, rows in enumerate(sheets.values(), start=1):
        out = []
        for i, row in enumerate(rows, start=1):
            cells = []
            for j, v in enumerate(row):
                if v is None:
                    continue
                if isinstance(v, (int, float)):
                    cells.append(f'<c r="{col(j)}{i}"><v>{v}</v></c>')
                else:
                    idx.setdefault(v, len(shared))
                    if idx[v] == len(shared):
                        shared.append(v)
                    cells.append(f'<c r="{col(j)}{i}" t="s"><v>{idx[v]}</v></c>')
            out.append(f'<row r="{i}">{"".join(cells)}</row>')
        xml[k] = f'<worksheet xmlns="{M}"><sheetData>{"".join(out)}</sheetData></worksheet>'
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("xl/workbook.xml", f'<workbook xmlns="{M}" xmlns:r="{Rn}"><sheets>' + "".join(
            f'<sheet name="{escape(n)}" sheetId="{k}" r:id="rId{k}"/>' for k, n in enumerate(sheets, start=1)) + "</sheets></workbook>")
        z.writestr("xl/_rels/workbook.xml.rels", f'<Relationships xmlns="{P}">' + "".join(
            f'<Relationship Id="rId{k}" Type="{Rn}/worksheet" Target="worksheets/sheet{k}.xml"/>' for k in xml) + "</Relationships>")
        z.writestr("xl/sharedStrings.xml", f'<sst xmlns="{M}">' + "".join(f"<si><t>{escape(v)}</t></si>" for v in shared) + "</sst>")
        for k, x in xml.items():
            z.writestr(f"xl/worksheets/sheet{k}.xml", x)


def test_xlsx_sheet_reads_the_named_sheet_as_strings(tmp_path):
    _xlsx(tmp_path / "t.xlsx", {"Other": [["a"], [1]],
                                "Sheet1": [["Sample_ID", "status", "Patient.ID"], [69, "Involved", "pat. 5"],
                                           [68, "Uninvolved", "pat. 5"], [None], ["Column Legend:"],
                                           ["Sample_ID", "Sample Identifier"]]})
    t = R.xlsx_sheet(tmp_path / "t.xlsx", "Sheet1")
    assert list(t.columns) == ["Sample_ID", "status", "Patient.ID"]
    assert t["Sample_ID"].tolist() == ["69", "68", "Column Legend:", "Sample_ID"]
    assert t["Patient.ID"].tolist()[:2] == ["pat. 5", "pat. 5"]
    assert R.xlsx_sheet(tmp_path / "t.xlsx")["a"].tolist() == ["1"]  # no name: the first sheet


def test_tenx_raw_barcodes_with_a_threshold_and_the_patients_from_the_paper_s_table(tmp_path):
    """Cell Ranger's raw matrices (cells and near-empty barcodes), each sample's patient from the paper's table, one
    pair of samples the table does not name: its cells are listed, not checked."""
    rng = np.random.default_rng(3)
    gsms, titles, status, chem, table = [], [], [], [], [["Sample_ID", "tissue", "status", "Patient.ID"]]
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as tar:
        for k in range(7):
            for s, st in enumerate(("Involved", "Uninvolved")):
                num = 100 + 2 * k + s
                gsm = f"GSM{num}"
                gsms.append(gsm)
                titles.append(f"Ileal {st} {num}")
                status.append(st)
                chem.append("V1" if k == 0 else "V2")
                if k < 6:
                    table.append([num, "ILEUM", st, f"pat. {k}"])
                cells = _cells(rng, 20, female=k % 2 == 0) * 10
                empty = rng.poisson(0.1, size=(30, len(GENES))).astype(float)
                X = np.vstack([cells, empty])
                files = {f"{gsm}_{num}_barcodes.tsv.gz": gzip.compress("".join(f"B{i}-1\n" for i in range(50)).encode()),
                         f"{gsm}_{num}_genes.tsv.gz": gzip.compress("".join(f"E{i}\t{g}\n" for i, g in enumerate(GENES)).encode()),
                         f"{gsm}_{num}_matrix.mtx.gz": _mtx(X.T)}
                for name, data in files.items():
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    tar.addfile(info, io.BytesIO(data))
    (tmp_path / "raw.tar").write_bytes(buf.getvalue())
    _series(tmp_path / "series.txt.gz", gsms, titles, {"tissue": ["ileal"] * len(gsms), "status": status,
                                                       "10x chemistry": chem})
    _xlsx(tmp_path / "pairs.xlsx", {"Sheet1": table})
    spec = dict(rid="R00004", organism="homo_sapiens", loader="tenx_samples",
                files=dict(raw="https://example.org/raw.tar", series="https://example.org/series.txt.gz",
                           pairs="https://example.org/pairs.xlsx"),
                loader_args=dict(min_umis=20),
                tables=[dict({"as": "patient"}, file="pairs", sheet="Sheet1", key="Sample_ID", value="Patient.ID",
                             on=dict(column="title", regex=r"(\d+)$"))],
                cells=dict(tissue=["ileal"]), groups=dict(column="status", values=["Involved", "Uninvolved"]),
                donor="patient", source="10x chemistry", sex_label=None,
                claims=[dict(id="1", genes=["ACE2", "TMPRSS2"])])
    (tmp_path / "spec.json").write_text(json.dumps(spec))
    assert R.main(["--spec", str(tmp_path / "spec.json"), "--out", str(tmp_path / "out"), "--files", str(tmp_path)]) == 0
    res = json.loads((tmp_path / "out" / "report.json").read_text())
    assert res["cells"] == 6 * 2 * 20  # the near-empty barcodes and the unnamed pair left out
    assert res["without_donor"] == {"Ileal Involved 112 / Involved": 20, "Ileal Uninvolved 113 / Uninvolved": 20}
    assert res["donors"]["design"] == "paired" and res["donors"]["donors"] == {"Involved": 6, "Uninvolved": 6}
    assert res["donors"]["rule4"]["passes"] and res["source"]["verdict"] == "balanced"
    assert res["sex"]["both_sexes"] and res["sex"]["fails"] == []
    md = (tmp_path / "out" / "report.md").read_text()
    assert "cells no donor names, left out" in md and "Ileal Involved 112" in md
    log = (tmp_path / "out" / "data_run.log").read_text()
    assert "GSM100: 20 of 50 barcodes with at least 20 counts" in log and "table pairs: patient for 240 of 280 cells" in log


def test_an_r_object_s_export_with_the_paper_s_filters_and_its_minimum_of_cells(tmp_path):
    """r_export.R's files (counts genes by cells, genes, cells), the paper's cells (a tissue, treatments, responders)
    and its rule of at least 30 cells before and after: a patient with fewer is left out and listed."""
    rng = np.random.default_rng(4)
    cells, blocks = [], []
    for k in range(10):
        for grp in ("Pre-treatment", "Post-treatment"):
            n = 20 if (k == 9 and grp == "Post-treatment") else 40
            blocks.append(_cells(rng, n, female=True))
            cells += [dict(cell=f"C{len(cells) + i}", Patient=f"P{k:03d}", Group=grp,
                           Tissue="liver" if k == 8 else "breast", Treatment="PTX" if k == 7 else "Nab-PTX",
                           Response="R") for i in range(n)]
    X = np.vstack(blocks)
    exp = tmp_path / "export"
    exp.mkdir()
    (exp / "counts.mtx.gz").write_bytes(_mtx(X.T))
    (exp / "genes.tsv.gz").write_bytes(gzip.compress(("\n".join(GENES) + "\n").encode(), mtime=0))
    pd.DataFrame(cells).to_csv(exp / "cells.csv.gz", index=False)
    (tmp_path / "object.rds.gz").write_bytes(b"FABRICATED")
    spec = dict(rid="R00005", organism="homo_sapiens", loader="r_export",
                files=dict(object="https://example.org/object.rds.gz"), loader_args=dict(export_dir=str(exp)),
                cells=dict(Tissue=["breast"], Treatment=["Nab-PTX", "Nab-PTX+Anti-PD-L1"], Response=["R"]),
                groups=dict(column="Group", values=["Post-treatment", "Pre-treatment"]), min_cells_each_group=30,
                donor="Patient", source=None, sex_label=None, claims=[dict(id="1", genes=["ACE2", "TMPRSS2"])])
    (tmp_path / "spec.json").write_text(json.dumps(spec))
    assert R.main(["--spec", str(tmp_path / "spec.json"), "--out", str(tmp_path / "out"), "--files", str(tmp_path)]) == 0
    res = json.loads((tmp_path / "out" / "report.json").read_text())
    assert res["donors"]["design"] == "paired" and res["donors"]["donors"] == {"Post-treatment": 7, "Pre-treatment": 7}
    assert res["donors_left_out"] == {"P009": {"Post-treatment": 20, "Pre-treatment": 40}}
    assert res["cells"] == 7 * 2 * 40 and res["donors"]["rule4"]["passes"]
    assert "left out as the paper did" in (tmp_path / "out" / "report.md").read_text()
    # --fetch-only: the files and their record, nothing read
    assert R.main(["--spec", str(tmp_path / "spec.json"), "--out", str(tmp_path / "f"), "--files", str(tmp_path),
                   "--fetch-only"]) == 0
    assert (tmp_path / "f" / "files.tsv").exists() and not (tmp_path / "f" / "report.json").exists()


def test_an_r_export_named_by_its_symbol_column(tmp_path):
    exp = tmp_path / "export"
    exp.mkdir()
    X = _cells(np.random.default_rng(5), 6, female=True)
    ids = [f"ENSG{i:011d}" for i in range(len(GENES))]
    (exp / "counts.mtx.gz").write_bytes(_mtx(X.T))
    (exp / "genes.tsv.gz").write_bytes(gzip.compress(("\n".join(ids) + "\n").encode(), mtime=0))
    pd.DataFrame(dict(gene=ids, Symbol=GENES)).to_csv(exp / "rowdata.csv.gz", index=False)
    pd.DataFrame(dict(cell=[f"c{i}" for i in range(6)])).to_csv(exp / "cells.csv.gz", index=False)
    _, genes, obs = R.load_r_export({}, dict(export_dir=str(exp), symbol_column="Symbol"), lambda s: None)
    assert genes == GENES and len(obs) == 6
    _, genes, _ = R.load_r_export({}, dict(export_dir=str(exp)), lambda s: None)
    assert genes == ids


def test_a_10x_hdf5_matrix_with_the_deposit_s_table_of_cells_batches_and_a_mapped_sex_label(tmp_path):
    import h5py
    from scipy import sparse
    rng = np.random.default_rng(6)
    blocks, meta = [], []
    for k in range(8):
        female = k % 2 == 0
        blocks.append(_cells(rng, 30, female=female))
        meta += [dict(Barcode=f"B{len(meta) + i}-{k + 1}", SampleID=f"Sample-{k}", Diagnosis="AD" if k < 4 else "Control",
                      Batch=1 + k % 2, **{"Cell.Type": "ODC"}, Sex="F" if female else "M") for i in range(30)]
    X = np.vstack(blocks)
    M = sparse.csc_matrix(X.T)
    with h5py.File(tmp_path / "m.h5", "w") as f:
        g = f.create_group("matrix")
        g["data"], g["indices"], g["indptr"], g["shape"] = M.data, M.indices, M.indptr, np.array(M.shape)
        g["barcodes"] = np.array([m["Barcode"] for m in meta] + ["EMPTY-1"], dtype="S")[:M.shape[1]]
        ft = g.create_group("features")
        ft["name"] = np.array(GENES, dtype="S")
        ft["id"] = np.array([f"E{i}" for i in range(len(GENES))], dtype="S")
        ft["feature_type"] = np.array(["Gene Expression"] * len(GENES), dtype="S")
    pd.DataFrame(meta).to_csv(tmp_path / "meta.csv.gz", index=False)
    spec = dict(rid="R00006", organism="homo_sapiens", loader="tenx_h5",
                files=dict(matrix="https://example.org/m.h5", meta="https://example.org/meta.csv.gz"),
                cells={"Cell.Type": ["ODC"]}, groups=dict(column="Diagnosis", values=["AD", "Control"]),
                donor="SampleID", source="Batch", sex_label=dict(column="Sex", map=dict(M="male", F="female")),
                claims=[dict(id="1", genes=["ACE2", "TMPRSS2"])])
    (tmp_path / "spec.json").write_text(json.dumps(spec))
    assert R.main(["--spec", str(tmp_path / "spec.json"), "--out", str(tmp_path / "out"), "--files", str(tmp_path)]) == 0
    res = json.loads((tmp_path / "out" / "report.json").read_text())
    assert res["cells"] == 240 and res["donors"]["design"] == "nested" and res["donors"]["donors"] == {"AD": 4, "Control": 4}
    assert res["source"]["verdict"] == "balanced" and res["sex"]["both_sexes"] and res["sex"]["fails"] == []
    assert {d["label"] for d in res["sex"]["donors"]} == {"female", "male"}


def test_per_sample_count_tables_and_one_table_of_series_samples(tmp_path):
    rng = np.random.default_rng(7)
    gsms, titles, dis = [], [], []
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as tar:
        for k in range(8):
            gsm = f"GSM{500 + k}"
            gsms.append(gsm)
            titles.append(f"Donor{100 + k}_{'t2d' if k < 4 else 'control'}")
            dis.append("T2D" if k < 4 else "non-T2D")
            X = _cells(rng, 12, female=True)
            t = pd.DataFrame(X.T, index=GENES, columns=[f"A{i}" for i in range(12)])
            t.index.name = "gene"
            for name, data in {f"{gsm}_sample{100 + k}_counts.tab.gz": gzip.compress(t.to_csv(sep="\t").encode(), mtime=0),
                               f"{gsm}_sample{100 + k}_rpkms.tab.gz": gzip.compress(b"x", mtime=0)}.items():
                info = tarfile.TarInfo(name)
                info.size = len(data)
                tar.addfile(info, io.BytesIO(data))
    (tmp_path / "raw.tar").write_bytes(buf.getvalue())
    _series(tmp_path / "series.txt.gz", gsms, titles, {"disease state": dis})
    spec = dict(rid="R00007", organism="homo_sapiens", loader="tab_samples",
                files=dict(raw="https://example.org/raw.tar", series="https://example.org/series.txt.gz"),
                groups={"column": "disease state", "values": ["T2D", "non-T2D"]},
                donor=dict(column="title", regex=r"(Donor\d+)"), source=None, sex_label=None,
                claims=[dict(id="1", genes=["ACE2", "TMPRSS2"])])
    (tmp_path / "spec.json").write_text(json.dumps(spec))
    assert R.main(["--spec", str(tmp_path / "spec.json"), "--out", str(tmp_path / "out"), "--files", str(tmp_path)]) == 0
    res = json.loads((tmp_path / "out" / "report.json").read_text())
    assert res["cells"] == 96 and res["counts"]["raw"] and res["donors"]["donors"] == {"T2D": 4, "non-T2D": 4}
    # one table whose columns are the series' samples (a cell each), non-integer values: the counts check fails
    cells = [f"HP{100 + k // 3}{'T2D' if k < 12 else ''}_A{k}" for k in range(24)]
    tab = pd.DataFrame(rng.gamma(1.0, 2.0, size=(len(GENES), 24)), index=GENES, columns=cells)
    tab.to_csv(tmp_path / "rpkms.txt.gz", sep="\t")
    _series(tmp_path / "s2.txt.gz", [f"GSM{900 + k}" for k in range(24)], cells,
            {"disease state": ["T2D" if k < 12 else "Healthy" for k in range(24)]})
    spec2 = dict(rid="R00008", organism="homo_sapiens", loader="table_series",
                 files=dict(table="https://example.org/rpkms.txt.gz", series="https://example.org/s2.txt.gz"),
                 groups={"column": "disease state", "values": ["T2D", "Healthy"]},
                 donor=dict(column="cell", regex=r"^(HP\d+)"), source=None, sex_label=None,
                 claims=[dict(id="1", genes=["ACE2", "TMPRSS2"])])
    (tmp_path / "spec2.json").write_text(json.dumps(spec2))
    assert R.main(["--spec", str(tmp_path / "spec2.json"), "--out", str(tmp_path / "out2"), "--files", str(tmp_path)]) == 0
    res2 = json.loads((tmp_path / "out2" / "report.json").read_text())
    assert res2["cells"] == 24 and not res2["counts"]["raw"] and res2["donors"]["donors"] == {"T2D": 4, "Healthy": 4}
