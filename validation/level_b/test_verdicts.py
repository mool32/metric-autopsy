"""Level B, step S4 (`plan_s4.md`): the verdicts' search, its query, its filter and its dry run."""
from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import verdicts as V  # noqa: E402

PAIR = dict(id="c1", kind="pair", genes=["ACE2", "TMPRSS2"])
SET = dict(id="c2", kind="set", genes=["IFI6", "MX1", "OAS1"], set_name="interferon module")


def test_the_query_names_both_genes_and_the_terms():
    q = V.claim_query(PAIR)
    assert q.startswith('(TITLE:"ACE2" OR ABSTRACT:"ACE2") AND (TITLE:"TMPRSS2" OR ABSTRACT:"TMPRSS2") AND (')
    for t in ("coexpression", "differential correlation", "knockdown", "smFISH", "in situ"):
        assert f'TITLE:"{t}" OR ABSTRACT:"{t}"' in q
    assert V.claim_query(SET).startswith('(TITLE:"interferon module" OR ABSTRACT:"interferon module") AND (')
    assert V.claim_query(dict(SET, set_name="")) is None


def test_the_filter_wants_both_genes_as_words_or_the_sets_name_or_two_genes():
    assert V.passes(dict(title="ACE2 and TMPRSS2 in the gut", abstract=""), PAIR)
    assert not V.passes(dict(title="ACE2 alone", abstract="TMPRSS2X is another gene"), PAIR)
    assert V.passes(dict(title="", abstract="tmprss2 rises with ace2"), PAIR)
    assert V.passes(dict(title="The Interferon Module of monocytes", abstract=""), SET)
    assert V.passes(dict(title="IFI6 and MX1", abstract=""), SET)
    assert not V.passes(dict(title="IFI6 only", abstract=""), SET)


def test_a_papers_citations_key():
    assert V.paper_key(dict(pmid="1", pmcid="PMC2")) == dict(source="MED", id="1")
    assert V.paper_key(dict(pmid="", pmcid="PMC2")) == dict(source="PMC", id="PMC2")
    assert V.paper_key(dict(pmid="", pmcid="", epmc_ids=["PPR:PPR77"])) == dict(source="PPR", id="PPR77")
    assert V.paper_key(dict(pmid="", pmcid="", epmc_ids=[])) is None


def test_the_dry_run(tmp_path):
    kept = dict(papers=[dict(rid="R00189", pmid="12345612", pmcid="", doi="10.1/x", epmc_ids=[], claims=[PAIR]),
                        dict(rid="R09999", pmid="", pmcid="", doi="", epmc_ids=["PPR:PPR13"], claims=[SET])])
    (tmp_path / "kept.json").write_text(json.dumps(kept))
    assert V.main(["--claims", str(tmp_path / "kept.json"), "--out", str(tmp_path / "out"), "--fabricated"]) == 0
    res = json.loads((tmp_path / "out" / "verdicts.json").read_text())
    by = {(c["rid"], c["claim"]): c for c in res["claims"]}
    assert by[("R00189", "c1")]["citing"] == 2 and by[("R00189", "c1")]["by_terms"] == 1
    assert by[("R00189", "c1")]["kept"] == 2  # the citing record naming both genes, and the query's knockdown study
    assert by[("R09999", "c2")]["kept"] >= 1
    found = [json.loads(x) for x in gzip.open(tmp_path / "out" / "found.jsonl.gz", "rt")]
    assert [f["claim"] for f in found] == ["c1", "c2"]
    assert (tmp_path / "out" / "fulltext_index.tsv").read_text().startswith("epmc_id\ttried\tfulltext")
    assert "FABRICATED" in (tmp_path / "out" / "verdicts.log").read_text()


def test_every_version_a_record_joins_has_its_citations_asked():
    """2026-11-11: a record that joins a preprint and its journal article is cited under both PMIDs."""
    keys = V.paper_keys(dict(pmid="32676577", pmcid="PMC7728249", epmc_ids=["MED:33156843"]))
    assert keys == [dict(source="MED", id="32676577"), dict(source="MED", id="33156843")]
    assert V.paper_keys(dict(pmid="1", epmc_ids=["MED:1", "PPR:PPR9"])) == [dict(source="MED", id="1")]
    assert V.paper_keys(dict(pmid="", pmcid="PMC2", epmc_ids=[])) == [dict(source="PMC", id="PMC2")]
    assert V.paper_keys(dict(pmid="", pmcid="", epmc_ids=[])) == []



def test_the_full_texts_with_short_attempts_a_time_budget_and_the_outage_rule(tmp_path):
    """2026-11-11 (J30): a full text that fails is recorded and the next asked; the time spent stops the asking."""
    import search as S

    class Net:
        def __init__(self, fail):
            self.fail, self.calls = fail, []

        def get(self, url, params, post=False, attempts=6, timeout=180):
            self.calls.append((url, attempts, timeout))
            fid = url.rsplit("/", 2)[-2]
            if fid in self.fail:
                raise S.Unavailable("no answer", None)
            return b"<article/>"

    out = S.Out(tmp_path / "o")
    recs = [dict(epmc_id=f"MED:{i}", pmcid=f"PMC{i}") for i in range(8)]
    rows = V.fulltexts(Net({"PMC1"}), out, recs)
    assert [r["fulltext"] for r in rows].count("-") == 1 and "no answer on both attempts" in rows[1]["note"]
    net = Net(set())
    V.fulltexts(net, out, recs[:1])
    assert net.calls[0][1:] == (V.FT_ATTEMPTS, V.FT_TIMEOUT)
    ticks = iter(range(0, 1000, 10))
    rows = V.fulltexts(Net(set()), out, recs, budget=25, clock=lambda: next(ticks))
    assert [r["note"] for r in rows][-1] == "not fetched: the run's time budget" and rows[0]["fulltext"] == "PMC0"
    try:
        V.fulltexts(Net({f"PMC{i}" for i in range(7)}), out, recs)
    except S.Unavailable as e:
        assert "an outage" in str(e)
    else:
        raise AssertionError("no stop")
