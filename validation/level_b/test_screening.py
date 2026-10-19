"""Tests of level B's screening tools (`screening.py`) on fabricated records and decisions.

    python -m pytest -q -p no:cacheprovider validation/level_b/test_screening.py
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import screening as S  # noqa: E402


def _search(tmp: Path, n: int = 450) -> Path:
    d = tmp / "search"
    d.mkdir()
    recs = [dict(rid=f"R{i:05d}", title=f"Title {i}", abstract=f"Abstract {i}", journal="J", year=2020, pub_types=[],
                 pmid=str(1000 + i), doi=f"10.1/{i}", preprint=False, found_by=["pubmed"], pmids=[str(1000 + i)])
            for i in range(1, n + 1)]
    with gzip.open(d / "records.jsonl.gz", "wt") as f:
        for r in recs:
            f.write(json.dumps(r) + "\n")
    (d / "sample.txt").write_text("\n".join(r["rid"] for r in recs) + "\n")
    lines = [f"{hashlib.sha256((d / name).read_bytes()).hexdigest()}  {name}" for name in ("records.jsonl.gz", "sample.txt")]
    (d / "SHA256SUMS").write_text("\n".join(lines) + "\n")
    return d


def _write_pass(batches: Path, out: Path, decide):
    out.mkdir()
    man = json.loads((batches / "manifest.json").read_text())
    for b in man["batches"]:
        rids = S.read_batch(batches / b["name"])
        rows = ["\t".join(S.HEADER)] + ["\t".join(decide(r)) for r in rids]
        (out / b["name"].replace(".jsonl", ".tsv")).write_text("\n".join(rows) + "\n")


def test_the_batches_hold_the_sample_in_order_without_how_records_were_found(tmp_path):
    man = S.make_batches(_search(tmp_path), tmp_path / "b")
    assert [b["records"] for b in man["batches"]] == [200, 200, 50] and man["n_screened"] == 450
    first = json.loads((tmp_path / "b" / "batch_01.jsonl").read_text().splitlines()[0])
    assert list(first) == list(S.FIELDS) and first["rid"] == "R00001" and "found_by" not in first


def test_the_batches_refuse_a_list_that_is_not_the_searchs(tmp_path):
    d = _search(tmp_path)
    (d / "sample.txt").write_text("R00001\n")
    try:
        S.make_batches(d, tmp_path / "b")
    except SystemExit as e:
        assert "sample.txt" in str(e)
    else:
        raise AssertionError("a changed sample was accepted")


def test_cohen_kappa():
    assert S.cohen_kappa(list("KKEE"), list("KKEE"), ("K", "E")) == 1.0
    assert S.cohen_kappa(list("KE"), list("EK"), ("K", "E")) == -1.0
    # 20 items: 15 agree (K 5, E 10), 5 differ; pe = (7*8 + 13*12) / 400 = 0.53
    a = list("K" * 5 + "E" * 10 + "K" * 2 + "E" * 3)
    b = list("K" * 5 + "E" * 10 + "E" * 2 + "K" * 3)
    assert math.isclose(S.cohen_kappa(a, b, ("K", "E")), (0.75 - 0.53) / (1 - 0.53))


def test_check_finds_what_breaks_the_format(tmp_path):
    S.make_batches(_search(tmp_path), tmp_path / "b")
    _write_pass(tmp_path / "b", tmp_path / "p", lambda r: [r, "KEEP", "X4", "pair", "a note"] if r == "R00003" else
                [r, "EXCLUDE", "X4", "-", "no claim"])
    d, problems = S.read_pass(tmp_path / "b", tmp_path / "p")
    assert len(d) == 450 and problems == ["batch_01.tsv: R00003: code 'X4' with KEEP"]
    (tmp_path / "p" / "batch_03.tsv").unlink()
    _, problems = S.read_pass(tmp_path / "b", tmp_path / "p")
    assert "batch_03.tsv: missing or empty" in problems


def test_the_merge_counts_both_passes_and_sends_disagreements_and_unclear_to_the_owner(tmp_path):
    S.make_batches(_search(tmp_path), tmp_path / "b")

    def one(r):
        i = int(r[1:])
        if i <= 10:
            return [r, "KEEP", "-", "pair", "a claim"]
        if i <= 15:
            return [r, "UNCLEAR", "-", "-", "no abstract"]
        return [r, "EXCLUDE", "X6" if i <= 40 else "X4", "set" if i <= 40 else "-", "no pair"]

    def two(r):
        i = int(r[1:])
        if i <= 8 or 11 <= i <= 12:
            return [r, "KEEP", "-", "pair", "a claim"]
        return [r, "EXCLUDE", "X6" if i <= 40 else "X4", "set" if i <= 40 else "-", "no pair"]
    _write_pass(tmp_path / "b", tmp_path / "p1", one)
    _write_pass(tmp_path / "b", tmp_path / "p2", two)
    res = S.merge(tmp_path / "b", tmp_path / "p1", tmp_path / "p2", tmp_path / "m")
    assert res["both_keep"] == 8 and res["to_owner"] == 7 and res["fulltext"] == 15
    assert res["table"]["pass 1"]["KEEP"] == 10 and res["table"]["pass 1"]["X6"] == 25
    assert res["sets"]["both passes"]["code_X6"] == 25
    assert res["crosstab"] == {"EXCLUDE / EXCLUDE": 435, "KEEP / EXCLUDE": 2, "KEEP / KEEP": 8,
                               "UNCLEAR / EXCLUDE": 3, "UNCLEAR / KEEP": 2}
    owner = (tmp_path / "m" / "to_owner.tsv").read_text().splitlines()
    assert len(owner) == 8 and owner[1].startswith("R00009\tKEEP")
    assert (tmp_path / "m" / "fulltext.txt").read_text().split() == [f"R{i:05d}" for i in range(1, 16)]


def test_the_screeners_text_holds_each_record_in_short_lines_and_check_one_reads_it(tmp_path):
    S.make_batches(_search(tmp_path), tmp_path / "b")
    text = (tmp_path / "b" / "batch_01.txt").read_text()
    assert text.startswith("=== R00001 ===\nTITLE: Title 1\nJOURNAL: J | YEAR: 2020 | TYPES: - | PREPRINT: no\n"
                           "ABSTRACT: Abstract 1\n")
    assert text.count("=== R") == 200 and max(len(line) for line in text.splitlines()) <= S.WIDTH
    long = S.as_text(dict(rid="R1", title="t", abstract="word " * 300))
    assert all(len(line) <= S.WIDTH for line in long.splitlines()) and long.count("word") == 300
    dec = tmp_path / "d.tsv"
    dec.write_text("\t".join(S.HEADER) + "\n" + "\n".join(f"R{i:05d}\tEXCLUDE\tX1\t-\ta review" for i in range(1, 201)) + "\n")
    assert S.main(["check-one", "--batch", str(tmp_path / "b" / "batch_01.txt"), "--decisions", str(dec)]) == 0
    dec.write_text("\t".join(S.HEADER) + "\nR00001\tEXCLUDE\tX1\t-\ta review\n")
    assert S.main(["check-one", "--batch", str(tmp_path / "b" / "batch_01.txt"), "--decisions", str(dec)]) == 1


def test_the_fulltext_batches_and_their_check(tmp_path):
    m = tmp_path / "m"
    m.mkdir()
    head = "rid\tdecision 1\tcode 1\tkinds 1\tnote 1\tdecision 2\tcode 2\tkinds 2\tnote 2\tpmid\tdoi\ttitle"
    rows = [f"R{i:05d}\t{'KEEP' if i <= 10 else 'EXCLUDE'}\t{'-' if i <= 10 else 'X4'}\t-\tn\t"
            f"{'KEEP' if i % 2 else 'EXCLUDE'}\t{'-' if i % 2 else 'X4'}\t-\tn\t\t\tt" for i in range(1, 21)]
    (m / "all.tsv").write_text(head + "\n" + "\n".join(rows) + "\n")
    texts = tmp_path / "texts"
    texts.mkdir()
    for i in range(1, 10):  # R00010 has no full text
        (texts / f"R{i:05d}.txt").write_text(f"=== R{i:05d} ===\n")
    res = S.fulltext_batches(m, 1, texts, tmp_path / "b1")
    assert res == dict(kept=10, with_fulltext=9, without=1, batches=["batch_01.list", "batch_02.list"])
    assert (tmp_path / "b1" / "without_fulltext.txt").read_text() == "R00010\n"
    assert S.fulltext_batches(m, 2, texts, tmp_path / "b2")["kept"] == 10
    batch = tmp_path / "b1" / "batch_02.list"
    assert batch.read_text().split("\t")[0] == "R00009"
    claim = ["A, B", "mouse", "T cells", "young; old", "lower in old", "Pearson", "results, Fig. 2", "a quote", "GSE1",
             "4 v 4", "yes", "-"]
    good = ["\t".join(S.FT_HEADER), "\t".join(["R00009", "KEEP", "-", "0"] + ["-"] * 12 + ["a pair falls with age"]),
            "\t".join(["R00009", "CLAIM", "-", "1"] + claim + ["-"])]
    dec = tmp_path / "d.tsv"
    dec.write_text("\n".join(good) + "\n")
    assert S.check_fulltext(batch, dec) == []
    dec.write_text("\n".join(good[:2]) + "\n")
    assert S.check_fulltext(batch, dec) == ["R00009: KEEP without a claim line"]
    dec.write_text("\n".join([good[0], good[1].replace("KEEP\t-", "EXCLUDE\t-")]) + "\n")
    assert S.check_fulltext(batch, dec) == ["R00009: code '-' with EXCLUDE"]
