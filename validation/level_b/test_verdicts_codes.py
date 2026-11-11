"""verdicts_codes.py: S4's batches, the check of a coder's file and the comparison of the passes (fabricated data:
the search's dry run, nothing fetched)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import verdicts as V  # noqa: E402
import verdicts_codes as C  # noqa: E402

LINE = {"1": dict(species="human", cells="enterocytes", conditions="inflamed vs not", direction="higher in inflamed",
                  metric="Pearson", where="Fig 6G", quote="co-expression of ACE2 and TMPRSS2 is elevated")}
PAIR = dict(id="1", kind="pair", genes=["ACE2", "TMPRSS2"], lines=LINE)
SET = dict(id="1", kind="set", genes=["IFI6", "MX1", "OAS1"], set_name="interferon module", lines=LINE)


def _searched(tmp_path: Path) -> tuple[Path, Path]:
    kept = dict(papers=[dict(rid="R00189", title="t", pmid="12345612", pmcid="", doi="10.1/x", epmc_ids=[],
                             claims=[PAIR, dict(PAIR, id="2")]),
                        dict(rid="R09999", title="u", pmid="", pmcid="", doi="", epmc_ids=["PPR:PPR13"], claims=[SET])])
    claims = tmp_path / "kept.json"
    claims.write_text(json.dumps(kept))
    assert V.main(["--claims", str(claims), "--out", str(tmp_path / "res"), "--fabricated"]) == 0
    return tmp_path / "res", claims


def _codes(batch: Path, decide) -> str:
    out = ["\t".join(C.HEADER)]
    for claim, record in C.read_list(batch):
        out.append("\t".join([claim, record] + decide(claim, record)))
    return "\n".join(out) + "\n"


NO = lambda c, r: ["NO_VERDICT", "-", "-", "-", "-", "same method"]  # noqa: E731


def test_the_batches_from_the_search_s_dry_run(tmp_path):
    res, claims = _searched(tmp_path)
    m = C.batches(res, claims, tmp_path / "texts", tmp_path / "b", size=3)
    by = {c["claim"]: c for c in m["claims"]}
    assert by["R00189:1"]["kept"] == 2 and by["R00189:2"]["kept"] == 2 and by["R09999:1"]["kept"] >= 1
    # R00189's four lines exceed the size: batches of their own; R09999's lines open the next batch
    assert m["batches"][:2] == ["batch_01", "batch_02"] and m["lines"] == sum(c["kept"] for c in m["claims"])
    first = C.read_list(tmp_path / "b" / "batch_01.list")
    assert len(first) == 3 and first[0][0] == "R00189:1"
    assert all(k[0].startswith("R09999") for k in C.read_list(tmp_path / "b" / f"{m['batches'][-1]}.list"))
    md = (tmp_path / "b" / "batch_01.md").read_text()
    assert "## Claim R00189:1 — a pair: **ACE2, TMPRSS2**" in md and "co-expression of ACE2 and TMPRSS2" in md
    # a record with a full text points to it, written as text and summed
    texts = sorted((tmp_path / "texts").glob("V*.txt"))
    assert texts and "FABRICATED" in texts[0].read_text()
    assert (tmp_path / "b" / "texts.sha256").read_text().count(".txt") == len(texts) == m["texts"]
    assert str(texts[0].resolve()) in "".join((tmp_path / "b" / f"{n}.list").read_text() for n in m["batches"])
    last = (tmp_path / "b" / f"{m['batches'][-1]}.md").read_text()
    assert "a set: **interferon module** (3 genes: IFI6, MX1, OAS1)" in last


def test_the_check_of_a_coder_s_file(tmp_path):
    res, claims = _searched(tmp_path)
    C.batches(res, claims, tmp_path / "texts", tmp_path / "b", size=20)
    batch = tmp_path / "b" / "batch_01.list"
    good = tmp_path / "good.tsv"
    good.write_text(_codes(batch, NO))
    assert C.check_codes(batch, good) == []
    verdict = lambda c, r: ["VERDICT", "causal", "holds", "knockdown of ACE2 lowered the pair's correlation",  # noqa: E731
                            "Results, Fig 2", "ACE2 knocked down; correlation measured"]
    good.write_text(_codes(batch, verdict))
    assert C.check_codes(batch, good) == []
    bad = tmp_path / "bad.tsv"
    bad.write_text(_codes(batch, lambda c, r: ["VERDICT", "causal", "maybe", "-", "-", "-"]))
    p = "\n".join(C.check_codes(batch, bad))
    assert "code 'maybe'" in p and "a verdict needs its quote" in p and "a note is needed" in p
    bad.write_text(_codes(batch, lambda c, r: ["NO_VERDICT", "causal", "-", "x " * 41, "abstract", "other cells"]))
    p = "\n".join(C.check_codes(batch, bad))
    assert "NO_VERDICT takes kind - and code -" in p and "41 words" in p
    lines = good.read_text().splitlines()
    bad.write_text("\n".join([lines[0]] + lines[2:] + lines[1:2]) + "\n")
    assert any("not the batch's in order" in x for x in C.check_codes(batch, bad))
    assert C.main(["check-codes", "--batch", str(batch), "--codes", str(good)]) == 0


def test_the_comparison_of_the_passes(tmp_path):
    res, claims = _searched(tmp_path)
    m = C.batches(res, claims, tmp_path / "texts", tmp_path / "b", size=20)
    rows = [k for n in m["batches"] for k in C.read_list(tmp_path / "b" / f"{n}.list")]
    first = rows[0]
    yes = lambda q: ["VERDICT", "replication", "holds", q, "Results", "smFISH of both genes"]  # noqa: E731
    for p, decide in (("p1", lambda c, r: yes("seen by smFISH") if (c, r) == first else NO(c, r)),
                      ("p2", lambda c, r: yes("smFISH again") if (c, r) in (first, rows[1]) else NO(c, r))):
        (tmp_path / p).mkdir()
        for n in m["batches"]:
            (tmp_path / p / f"{n}.tsv").write_text(_codes(tmp_path / "b" / f"{n}.list", decide))
    out = C.merge(tmp_path / "b", tmp_path / "p1", tmp_path / "p2", tmp_path / "merge")
    assert out["lines"] == len(rows) and out["verdicts_agreed"] == 1 and out["to_owner"] == 1
    assert out["pairs"]["VERDICT/replication/holds / VERDICT/replication/holds"] == 1
    assert 0 < out["kappa_decision"] < 1
    agreed = (tmp_path / "merge" / "verdicts.tsv").read_text().splitlines()
    assert len(agreed) == 2 and "seen by smFISH\tResults\tsmFISH again" in agreed[1]
    owner = (tmp_path / "merge" / "to_owner.tsv").read_text().splitlines()
    assert len(owner) == 2 and owner[1].startswith(f"{rows[1][0]}\t{rows[1][1]}\tNO_VERDICT")
    per = {c["claim"]: c for c in out["claims"]}
    assert per[first[0]]["verdicts_agreed"] == 1 and "claims_with_agreed_verdict" not in out
    # a file failing its check stops the comparison
    (tmp_path / "p2" / f"{m['batches'][0]}.tsv").write_text("claim\n")
    try:
        C.merge(tmp_path / "b", tmp_path / "p1", tmp_path / "p2", tmp_path / "merge2")
    except SystemExit as e:
        assert "fail their check" in str(e)
    else:
        raise AssertionError("no stop")
