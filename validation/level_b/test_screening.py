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

import pytest

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


def test_the_fulltext_merge_keeps_what_both_passes_keep_and_sends_the_rest_to_the_owner(tmp_path):
    m = tmp_path / "m"
    m.mkdir()
    head = "rid\tdecision 1\tcode 1\tkinds 1\tnote 1\tdecision 2\tcode 2\tkinds 2\tnote 2\tpmid\tdoi\ttitle"
    ta = {"R00001": ("KEEP", "KEEP"), "R00002": ("KEEP", "KEEP"), "R00003": ("KEEP", "EXCLUDE"), "R00004": ("KEEP", "KEEP")}
    (m / "all.tsv").write_text(head + "\n" + "\n".join(
        f"{r}\t{x}\t{'-' if x == 'KEEP' else 'X4'}\t-\tn\t{y}\t{'-' if y == 'KEEP' else 'X4'}\t-\tn\t\t\tTitle {r}"
        for r, (x, y) in ta.items()) + "\n")
    claim = ["A, B", "mouse", "T cells", "young; old", "lower in old", "Pearson", "results", "q", "GSE1", "4 v 4", "yes", "-", "-"]

    def write(ft, rows, without):
        ft.mkdir()
        rids = [r for r, *_ in rows]
        (ft / "batch_01.list").write_text("".join(f"{r}\t/x/{r}.txt\n" for r in rids))
        lines = ["\t".join(S.FT_HEADER)]
        for r, dec, code in rows:
            lines.append("\t".join([r, dec, code, "0"] + ["-"] * 12 + ["note"]))
            if dec == "KEEP":
                lines.append("\t".join([r, "CLAIM", "-", "1"] + claim))
        (ft / "batch_01.tsv").write_text("\n".join(lines) + "\n")
        (ft / "without_fulltext.txt").write_text("".join(f"{r}\n" for r in without))
    write(tmp_path / "f1", [("R00001", "KEEP", "-"), ("R00002", "EXCLUDE", "F7"), ("R00003", "KEEP", "-")], ["R00004"])
    write(tmp_path / "f2", [("R00001", "KEEP", "-"), ("R00002", "KEEP", "-")], ["R00004"])
    res = S.merge_fulltext(m, tmp_path / "f1", tmp_path / "f2", tmp_path / "out")
    assert res["read_by_both"] == 2 and res["kept_by_both"] == 1 and res["to_owner"] == 3
    assert res["crosstab"] == {"EXCLUDE / KEEP": 1, "KEEP / KEEP": 1}
    assert res["counts"]["pass 1"]["F7"] == 1 and res["counts"]["pass 2"]["without a full text"] == 1
    kept = (tmp_path / "out" / "ft_kept.tsv").read_text().splitlines()
    assert len(kept) == 3 and kept[1].startswith("1\tR00001\t1\tA, B")
    owner = (tmp_path / "out" / "ft_to_owner.txt").read_text()
    assert "R00003" in owner and "pass 2: not read (titles and abstracts: EXCLUDE)" in owner
    assert "R00004" in owner and "not read (no full text)" in owner


def _search_of(tmp: Path, name: str, texts: list) -> Path:
    """A search's list of records with these (title, abstract), numbered in this order."""
    d = tmp / name
    d.mkdir()
    with gzip.open(d / "records.jsonl.gz", "wt") as f:
        for i, (t, a) in enumerate(texts, 1):
            f.write(json.dumps(dict(rid=f"R{i:05d}", title=t, abstract=a, journal="J", year=2020, pub_types=[],
                                    pmid=str(1000 + i), doi="", preprint=False)) + "\n")
    (d / "sample.txt").write_text("\n".join(f"R{i:05d}" for i in range(1, len(texts) + 1)) + "\n")
    (d / "SHA256SUMS").write_text("".join(f"{hashlib.sha256((d / n).read_bytes()).hexdigest()}  {n}\n"
                                          for n in ("records.jsonl.gz", "sample.txt")))
    return d


def test_the_carry_keeps_the_decisions_of_unchanged_texts_and_screens_the_rest(tmp_path):
    first = [(f"Title {i}", f"Abstract {i}") for i in range(1, 451)]
    hidden = [(f"Hidden {i}", f"A paper a false merge hid {i}") for i in range(1, 11)]
    corrected = hidden[:5] + first[:100] + hidden[5:] + first[100:]  # the numbers shift
    S.make_batches(_search_of(tmp_path, "s1", first), tmp_path / "b1")
    S.make_batches(_search_of(tmp_path, "s2", corrected), tmp_path / "b2")

    def one(r):
        return [r, "EXCLUDE", "X1", "-", "Corrigendum to R00004" if r == "R00005" else "a review"]
    _write_pass(tmp_path / "b1", tmp_path / "p1", one)
    _write_pass(tmp_path / "b1", tmp_path / "p2", lambda r: [r, "KEEP", "-", "pair", "a claim"] if r == "R00450" else one(r))
    res = S.carry(tmp_path / "b1", tmp_path / "p1", tmp_path / "p2", tmp_path / "b2", tmp_path / "c")
    assert res["carried"] == 450 and res["to_screen"] == 10 and res["first_list_not_carried"] == []
    assert res["note_numbers_translated"] == 2 and res["note_numbers_untranslated"] == []
    carried = dict(line.split("\t") for line in (tmp_path / "c" / "carried.tsv").read_text().splitlines()[1:])
    assert carried["R00001"] == "-" and carried["R00006"] == "R00001" and carried["R00460"] == "R00450"
    p1 = (tmp_path / "c" / "carried_pass1.tsv").read_text()
    assert "R00010\tEXCLUDE\tX1\t-\tCorrigendum to R00009\n" in p1  # R00005 is now R00010, R00004 R00009
    rescreen = [json.loads(line)["title"] for line in (tmp_path / "c" / "rescreen" / "batch_01.jsonl").read_text().splitlines()]
    assert rescreen == [t for t, _ in hidden]
    _write_pass(tmp_path / "c" / "rescreen", tmp_path / "q1", lambda r: [r, "EXCLUDE", "X4", "-", "no claim"])
    _write_pass(tmp_path / "c" / "rescreen", tmp_path / "q2", lambda r: [r, "UNCLEAR", "-", "-", "no abstract"])
    out = S.combine(tmp_path / "b2", tmp_path / "c", tmp_path / "q1", tmp_path / "q2", tmp_path / "all")
    assert out["pass 1"] == dict(carried=450, screened_again=10, decisions=460)
    d2, problems = S.read_pass(tmp_path / "b2", tmp_path / "all" / "pass2")
    assert not problems and d2["R00460"]["decision"] == "KEEP" and d2["R00001"]["decision"] == "UNCLEAR"
    m = S.merge(tmp_path / "b2", tmp_path / "all" / "pass1", tmp_path / "all" / "pass2", tmp_path / "m")
    assert m["n_screened"] == 460 and m["to_owner"] == 11 and m["fulltext"] == 11


def test_the_claims_list_draws_30_lists_small_designs_and_the_papers_verdicts(tmp_path):
    import numpy as np
    m = tmp_path / "m"
    m.mkdir()
    n = 35
    rids = [f"R{i:05d}" for i in range(1, n + 1)]
    head = "rid\tdecision 1\tcode 1\tkinds 1\tnote 1\tdecision 2\tcode 2\tkinds 2\tnote 2\tpmid\tdoi\ttitle"
    (m / "all.tsv").write_text(head + "\n" + "\n".join(f"{r}\tKEEP\t-\tpair\tn\tKEEP\t-\tpair\tn\t{i}\t10.1/{i}\tTitle {r}"
                                                       for i, r in enumerate(rids, 1)) + "\n")
    table = {**{d: 0 for d in S.DECISIONS}, **{c: 0 for c in S.CODES}, "KEEP": n}
    (m / "merge.json").write_text(json.dumps(dict(n_screened=n, kappa=1.0, agreement=1.0, both_keep=n, to_owner=0,
                                                  fulltext=n, table={"pass 1": table, "pass 2": table},
                                                  sets={k: dict(code_X6=0, code_X7=0, kind_set=0, kind_global=0)
                                                        for k in ("pass 1", "pass 2")} |
                                                       {"both passes": dict(code_X6=0, code_X7=0, either_X6_or_X7=0)})))
    donors = {"R00001": "2 against 8", "R00002": "eight per group", "R00003": "6 paired donors"}

    def write(ft):
        ft.mkdir()
        (ft / "batch_01.list").write_text("".join(f"{r}\t/x/{r}.txt\n" for r in rids))
        lines = ["\t".join(S.FT_HEADER)]
        for r in rids:
            lines.append("\t".join([r, "KEEP", "-", "0"] + ["-"] * 12 + ["note"]))
            verdict = "smFISH: the correlation falls too" if r == "R00004" else "-"
            lines.append("\t".join([r, "CLAIM", "-", "1", "A, B", "mouse", "T cells", "young; old", "lower in old", "Pearson",
                                    "results", "q", "GSE1", donors.get(r, "4 vs 4"), "yes", verdict, "-"]))
        (ft / "batch_01.tsv").write_text("\n".join(lines) + "\n")
        (ft / "without_fulltext.txt").write_text("")
    write(tmp_path / "f1")
    write(tmp_path / "f2")
    S.merge_fulltext(m, tmp_path / "f1", tmp_path / "f2", tmp_path / "ftm")
    res = S.claims_list(m, tmp_path / "ftm", tmp_path / "f1", tmp_path / "f2", tmp_path / "out")
    assert res["papers_kept"] == 35 and res["papers_listed"] == 30 and res["drawn"] and not res["fewer_than_target"]
    assert res["small_designs"] == 1 and res["donors_unread"] == 1 and res["verdicts_found"] == 1
    listed = json.loads((tmp_path / "out" / "claims.json").read_text())["listed"]
    assert listed == sorted(rids[i] for i in np.random.default_rng(20261020).choice(35, 30, replace=False))
    md = (tmp_path / "out" / "claims.md").read_text()
    assert "Kept, not drawn: " in md and "- R00001, pass 1, claim 1: **A, B**" in md
    assert "R00002, pass 2, claim 1: eight per group" in md and "smFISH: the correlation falls too" in md


def test_the_owners_sheet_is_the_screening_files_own(capsys):
    """`decision_sheet.md` is what `decision_sheet.py` prints from the screening's files (the owner, 2026-11-10):
    the 12 and the 19 records sent to the owner, 23 in all, in four groups that hold every one of them once."""
    import decision_sheet as D
    assert D.main() == 0
    printed = capsys.readouterr().out
    assert printed == (HERE / "decision_sheet.md").read_text()
    ta, ft, index, allrows = D.read_lists()
    union = sorted(set(ta) | set(ft))
    groups = {g: [r for r in union if D.group_of(r, ta, ft) == g] for g, _, _ in D.GROUPS}
    assert {g: len(v) for g, v in groups.items()} == {"A": 9, "B": 8, "C": 2, "D": 4}
    assert sorted(r for v in groups.values() for r in v) == union and len(union) == 23


# --------------------------------------------------------------------------- #
# amendment 1 (plan_s1.md, plan_s2.md): the preprints merged with the list, and batches of some records
# --------------------------------------------------------------------------- #
def test_the_preprints_merge_by_id_doi_and_title_and_number_the_new_ones():
    import merge_preprints as M
    listing = [dict(rid="R00001", title="A coexpression study of single cells in disease", pmids=["111"], dois=["10.1/a"],
                    epmc_ids=["MED:111", "PPR:PPR500"]),
               dict(rid="R00002", title="Another record of the list with a long title", pmids=["222"], dois=["10.1/b"],
                    epmc_ids=["MED:222"])]
    pre = [dict(epmc_id="PPR:PPR500", title="whatever", doi="10.1101/x"),                         # in the list by id
           dict(epmc_id="PPR:PPR900", title="Another record of the list, with a long title!", doi="10.1101/y"),  # by title
           dict(epmc_id="PPR:PPR700", title="A new preprint about network rewiring in nuclei", doi="10.1101/z"),
           dict(epmc_id="PPR:PPR800", title="A new preprint about network rewiring in nuclei", doi="10.1101/w"),  # its twin
           dict(epmc_id="PPR:PPR600", title="Short", doi="10.1/B")]                                  # by DOI (case)
    fresh, rows = M.merge_preprints(listing, pre, first_new=3507)
    st = {r["epmc_id"]: (r["status"], r["rid"]) for r in rows}
    assert st["PPR:PPR500"] == ("in list", "R00001") and st["PPR:PPR900"] == ("in list", "R00002")
    assert st["PPR:PPR600"] == ("in list", "R00002")
    assert st["PPR:PPR700"] == ("new", "R03507") and st["PPR:PPR800"] == ("duplicate", "R03507")
    assert [r["rid"] for r in fresh] == ["R03507"] and fresh[0]["preprint"] and fresh[0]["found_by"] == ["preprints"]


def test_batches_of_some_records_and_the_new_ones(tmp_path):
    recs = [dict(rid=f"R{i:05d}", title=f"T{i}", abstract=f"A{i}", journal="J", year=2020, pub_types=[], pmid=str(i),
                 doi=f"10.1/{i}", preprint=False) for i in range(1, 11)]
    search = tmp_path / "list"
    search.mkdir()
    with gzip.open(search / "records.jsonl.gz", "wt") as f:
        for r in recs:
            f.write(json.dumps(r) + "\n")
    sums = tmp_path / "sums"
    sums.write_text(f"{hashlib.sha256((search / 'records.jsonl.gz').read_bytes()).hexdigest()}  records.jsonl.gz\n")
    (tmp_path / "rids.txt").write_text("R00007\nR00002\nR00009\n")
    extra = tmp_path / "new.jsonl.gz"
    with gzip.open(extra, "wt") as f:
        f.write(json.dumps(dict(rid="R03507", title="P", abstract="Q", journal="bioRxiv", year=2024, pub_types=["Preprint"],
                                pmid="", doi="10.1101/p", preprint=True)) + "\n")
    man = S.make_subset_batches(search, tmp_path / "rids.txt", tmp_path / "b", extra, size=2, prefix="sets", sums_file=sums)
    assert man["n_screened"] == 4 and [b["name"] for b in man["batches"]] == ["sets_01.jsonl", "sets_02.jsonl"]
    assert S.read_batch(tmp_path / "b" / "sets_01.jsonl") == ["R00002", "R00007"]
    assert S.read_batch(tmp_path / "b" / "sets_02.jsonl") == ["R00009", "R03507"]
    assert "PREPRINT: yes" in (tmp_path / "b" / "sets_02.txt").read_text()
    sums.write_text("0" * 64 + "  records.jsonl.gz\n")
    with pytest.raises(SystemExit):
        S.make_subset_batches(search, tmp_path / "rids.txt", tmp_path / "c", sums_file=sums)


def test_the_set_records_are_the_screenings_x6_or_x7_in_either_pass():
    rids = S.set_records(HERE / "screening" / "all.tsv")
    assert len(rids) == 210 and rids == sorted(rids)
    assert {"R00220", "R02780", "R02858"} <= set(rids)  # one pass kept each (the owner's group B)
    assert (HERE / "amendment1" / "set_records.txt").read_text().split() == rids


def test_amendment_1s_full_texts_list_the_new_records_send_the_rest_to_the_publisher_and_batch_all(tmp_path):
    recs = [dict(rid="R00001", pmid="1", pmids=["1"], pmcid="PMC1", doi="10.1/a", epmc_ids=["MED:1"]),
            dict(rid="R00002", pmid="2", pmids=["2"], pmcid="", doi="10.1/b", epmc_ids=["MED:2"])]
    search = tmp_path / "list"
    search.mkdir()
    with gzip.open(search / "records.jsonl.gz", "wt") as f:
        for r in recs:
            f.write(json.dumps(r) + "\n")
    extra = tmp_path / "new.jsonl.gz"
    with gzip.open(extra, "wt") as f:
        f.write(json.dumps(dict(rid="R03507", pmid="", pmids=[], pmcid="", doi="10.1101/p", epmc_ids=["PPR:PPR9"])) + "\n")
    (tmp_path / "rids.txt").write_text("R00001\nR00002\nR03507\n")
    assert S.fulltext_list(search, tmp_path / "rids.txt", tmp_path / "list.tsv", extra) == 3
    lines = (tmp_path / "list.tsv").read_text().splitlines()
    assert lines[3].split("\t") == ["R03507", "", "", "10.1101/p", "PPR:PPR9"]
    (tmp_path / "index.tsv").write_text("rid\ttried\tfulltext\tfile\tnote\nR00001\tPMC1\tPMC1\tf\t-\n"
                                        "R00002\t-\t-\t-\tno core record\nR03507\tPPR9\t-\t-\tPPR9: HTTP 404\n")
    assert S.publisher_list(tmp_path / "index.tsv", tmp_path / "list.tsv", tmp_path / "pub.tsv") == 2
    pub = [line.split("\t") for line in (tmp_path / "pub.tsv").read_text().splitlines()]
    assert pub[0] == ["rid", "pmid", "pmcid", "doi", "action", "repository"]
    assert [r[0] for r in pub[1:]] == ["R00002", "R03507"] and {r[4] for r in pub[1:]} == {"publisher"}
    texts = tmp_path / "texts"
    texts.mkdir()
    for rid in ("R00002", "R00001", "R03507"):
        (texts / f"{rid}.txt").write_text("x")
    res = S.fulltext_batches_all(["R03507", "R00001", "R00009", "R00002"], texts, tmp_path / "b", size=2)
    assert res["with_fulltext"] == 3 and res["without"] == 1 and res["batches"] == ["batch_01.list", "batch_02.list"]
    assert [line.split("\t")[0] for line in (tmp_path / "b" / "batch_01.list").read_text().splitlines()] == ["R00001", "R00002"]
    assert (tmp_path / "b" / "without_fulltext.txt").read_text() == "R00009\n"
