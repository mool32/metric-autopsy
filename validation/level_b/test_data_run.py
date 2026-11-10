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
