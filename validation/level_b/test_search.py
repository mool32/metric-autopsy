"""Tests of level B's search (`search.py`) on fabricated responses: nothing is fetched.

    python -m pytest -q -p no:cacheprovider validation/level_b/test_search.py
"""
from __future__ import annotations

import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import search as S  # noqa: E402
from fabricated import MISSING_PMID, N_PREPRINT_TWIN  # noqa: E402

PROTOCOL = HERE.parent / "plans" / "level_b_protocol.md"


def test_the_pubmed_query_is_the_protocols():
    text = PROTOCOL.read_text()
    section = text[text.index("## 3. Terms"):text.index("## 4. Dates")]
    block = re.search(r"```\n(.*?)\n```", section, re.S).group(1)
    assert " ".join(block.split()) == S.PUBMED_TERM
    assert "condition*" not in S.PUBMED_TERM and '"differential correlation"[tiab]' in S.PUBMED_TERM


def test_the_europepmc_query_holds_every_term_in_the_title_or_the_abstract():
    q = S.EUROPEPMC_QUERY
    for t in S.SINGLE_CELL + S.COORDINATION:
        assert f'TITLE:"{t}" OR ABSTRACT:"{t}"' in q
    assert "*" not in q.split("FIRST_PDATE")[0]  # no truncation inside a phrase
    assert 'LANG:"eng"' in q and "FIRST_PDATE:[2017-01-01 TO 2026-09-30]" in q
    # every PubMed term of coordination is there (truncation spelled out)
    for t in re.findall(r'"?([^"()]+?)"?\[tiab\]', S.PUBMED_TERM.split(") AND (")[1].split(") AND English")[0]):
        t = t.strip().replace("OR ", "").strip()
        assert t.rstrip("*") in q


def test_the_sample_is_the_protocols_draw():
    assert S.draw_sample(4000) == list(range(4000))
    s = S.draw_sample(6123)
    assert len(s) == 5000 and s == sorted(np.random.default_rng(20261019).choice(6123, 5000, replace=False).tolist())


def test_merging_by_pmid_doi_and_title():
    base = dict(abstract="", journal="J", year=2020, pub_types=[], languages=[], source="MED")
    pm = [dict(base, pmid="1", doi="10.1/a", title="A long enough title about co-expression"),
          dict(base, pmid="2", doi="", title="Another title about single cell coupling")]
    ep = [dict(base, pmid="1", doi="10.1/a", title="A long enough title about co-expression", epmc_id="MED:1"),
          dict(base, pmid="", doi="10.1/b", title="Another title about single-cell coupling!", source="PPR",
               epmc_id="PPR:9"),  # the same title once normalized: the preprint of record 2
          dict(base, pmid="", doi="10.1/c", title="short", source="PPR", epmc_id="PPR:8")]
    ci = [dict(base, pmid="", doi="10.1/A", title="x", cites="levy2020", epmc_id="MED:1")]
    ci[0]["doi"] = ci[0]["doi"].lower()
    m = S.merge(dict(pubmed=pm, europepmc=ep, citations=ci))
    assert len(m) == 3
    one = next(r for r in m if r["pmids"] == ["1"])
    assert one["found_by"] == ["cites:levy2020", "europepmc", "pubmed"]
    two = next(r for r in m if r["pmids"] == ["2"])
    assert two["dois"] == ["10.1/b"] and not two["preprint"]
    assert [r["rid"] for r in m] == ["R00001", "R00002", "R00003"] and m[2]["preprint"]


def test_pubmed_xml_is_parsed():
    xml = b"""<?xml version="1.0"?><PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>123</PMID><Article>
    <Journal><JournalIssue><PubDate><MedlineDate>2019 Jan-Feb</MedlineDate></PubDate></JournalIssue>
    <Title>Nature Methods</Title><ISOAbbreviation>Nat Methods</ISOAbbreviation></Journal>
    <ArticleTitle>Measures of <i>association</i></ArticleTitle><ELocationID EIdType="doi">10.1/XY</ELocationID>
    <Abstract><AbstractText Label="AIM">To test.</AbstractText><AbstractText>More.</AbstractText></Abstract>
    <Language>eng</Language><PublicationTypeList><PublicationType>Journal Article</PublicationType></PublicationTypeList>
    </Article></MedlineCitation><PubmedData><ArticleIdList><ArticleId IdType="pmc">PMC9</ArticleId></ArticleIdList>
    </PubmedData></PubmedArticle></PubmedArticleSet>"""
    (r,) = S.parse_pubmed_xml(xml)
    assert r["pmid"] == "123" and r["year"] == 2019 and r["doi"] == "10.1/xy" and r["pmcid"] == "PMC9"
    assert r["title"] == "Measures of association" and r["abstract"] == "AIM: To test.\nMore."
    assert r["journal"] == "Nat Methods" and r["pub_types"] == ["Journal Article"]


@pytest.fixture(scope="module")
def dry(tmp_path_factory):
    out = tmp_path_factory.mktemp("dry")
    assert S.main(["lookup", "--out", str(out / "lookup"), "--fabricated"]) == 0
    papers = out / "lookup" / "method_papers.verified.json"
    rc = S.main(["search", "--out", str(out / "search"), "--fabricated", "--papers", str(papers),
                 "--window-limit", "300"])
    return out, rc


def test_the_lookup_verifies_one_journal_article_per_paper(dry):
    out, _ = dry
    res = json.loads((out / "lookup" / "lookup.json").read_text())
    assert res["ok"] and all(p["n_matching"] == 1 and p["verified"]["source"] == "MED" for p in res["papers"])
    assert all(p["n_records"] == 3 for p in res["papers"])  # a preprint and an erratum left aside
    assert (out / "lookup" / "raw" / "lookup" / "levy2020.json.gz").exists()


def test_the_lookup_fails_where_a_paper_does_not_match(tmp_path):
    papers = json.loads(S.PAPERS.read_text())
    papers["papers"][3]["expect"]["year"] = 2021  # CS-CORE in the wrong year
    f = tmp_path / "papers.json"
    f.write_text(json.dumps(papers))
    assert S.main(["lookup", "--out", str(tmp_path / "out"), "--fabricated", "--papers", str(f)]) == 1
    assert "NOT VERIFIED" in (tmp_path / "out" / "lookup.log").read_text()


def test_the_search_refuses_unverified_papers(tmp_path):
    papers = json.loads(S.PAPERS.read_text())
    papers["papers"][4]["verified"] = None
    f = tmp_path / "papers.json"
    f.write_text(json.dumps(papers))
    assert S.main(["search", "--out", str(tmp_path / "out"), "--fabricated", "--papers", str(f)]) == 1
    assert "not all verified" in (tmp_path / "out" / "search.log").read_text()


def test_the_repositorys_method_papers_are_verified_by_the_lookup():
    papers = json.loads(S.PAPERS.read_text())["papers"]
    assert [p["key"] for p in papers] == ["levy2020", "gcl2021", "amit2022", "cscore2023", "skinnider2019"]
    assert all(p["verified"]["source"] == "MED" and p["verified"]["year"] == p["expect"]["year"] for p in papers)
    assert {p["key"]: p["verified"]["id"] for p in papers} == dict(
        levy2020="33139959", gcl2021="34040065", amit2022="35534606", cscore2023="37563115", skinnider2019="30962620")


def test_the_dry_run_of_the_search(dry):
    out, rc = dry
    assert rc == 0
    d = out / "search"
    summary = json.loads((d / "search.json").read_text())
    log = (d / "search.log").read_text()
    with gzip.open(d / "records.jsonl.gz", "rt") as f:
        recs = [json.loads(line) for line in f]
    # PubMed by windows, split where a window held more than the limit
    assert "split" in log and summary["pubmed"]["unique"] == 5600 and summary["pubmed"]["fetched"] == 5599
    assert f"['{MISSING_PMID}']" in log
    assert summary["europepmc"]["fetched"] == summary["europepmc"]["hit_count"] == 5800
    assert all(v["hit_count"] == 40 and v["unique"] == 40 for v in summary["citations"].values())
    assert sum(v["without_core"] for v in summary["citations"].values()) >= 5
    # merged: 5,600 PubMed records, 500 only in Europe PMC, the preprints without a twin, the new citing ones
    assert summary["merged"] == len(recs) == 5600 + 500 + (100 - N_PREPRINT_TWIN) + 5 * (8 + 2)
    assert [r["rid"] for r in recs] == [f"R{i:05d}" for i in range(1, len(recs) + 1)]
    twins = [r for r in recs if r["pmids"] and any(e.startswith("PPR:") for e in r["epmc_ids"])]
    assert len(twins) == N_PREPRINT_TWIN
    missing = next(r for r in recs if MISSING_PMID in r["pmids"])
    assert missing["found_by"] == ["europepmc"]  # efetch returned nothing, Europe PMC had it
    # the sample
    sample = (d / "sample.txt").read_text().split()
    assert len(sample) == 5000 == summary["sample"] and summary["sample_seed"] == 20261019
    idx = sorted(np.random.default_rng(20261019).choice(len(recs), 5000, replace=False).tolist())
    assert sample == [recs[i]["rid"] for i in idx]
    assert sum(bool(r.get("in_sample")) for r in recs) == 5000
    # every response written as it came, every request recorded
    raws = list((d / "raw").rglob("*.gz"))
    assert len(raws) > 50 and (d / "requests.tsv").read_text().count("\n") > 50
    assert summary["numpy"] == np.__version__


def test_the_dry_run_of_the_full_texts(dry, tmp_path):
    import screening
    out, _ = dry
    with gzip.open(out / "search" / "records.jsonl.gz", "rt") as f:
        recs = [json.loads(line) for line in f]
    pick = recs[:60] + [r for r in recs if r["preprint"]][:10]
    (tmp_path / "rids.txt").write_text("\n".join(r["rid"] for r in pick) + "\n")
    assert screening.fulltext_list(out / "search", tmp_path / "rids.txt", tmp_path / "list.tsv") == 70
    assert S.main(["fulltext", "--out", str(tmp_path / "ft"), "--fabricated", "--list", str(tmp_path / "list.tsv")]) == 0
    summary = json.loads((tmp_path / "ft" / "fulltext.json").read_text())
    index = [line.split("\t") for line in (tmp_path / "ft" / "fulltext_index.tsv").read_text().splitlines()[1:]]
    assert summary["n"] == len(index) == 70 and summary["fetched"] + summary["none"] == 70
    assert 0 < summary["fetched"] < 70 and summary["core_asked"] > 0
    with_pmc = [r for r in pick if r.get("pmcid")]
    assert with_pmc and all(row[1] != "-" for row in index if row[0] in {r["rid"] for r in with_pmc})
    assert any(row[2].startswith("PPR") for row in index) and any("HTTP 404" in row[4] for row in index)
    res = screening.fulltext_text(tmp_path / "ft", tmp_path / "text")
    assert res["converted"] == summary["fetched"] and not res["failed"]
    one = next((tmp_path / "text").glob("*.txt")).read_text()
    assert "### Results" in one and "[Figure 1] The pair by age." in one and "A reference" not in one
    assert "## Abstract" in one and "Eight mice per group" in one
