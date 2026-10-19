"""Fabricated responses for the dry run of `search.py` (nothing is fetched): PubMed's esearch and efetch,
Europe PMC's search and citations endpoint, as `search.Net.get` would return them.

The fabricated world, fixed: 5,600 PubMed records over 2017–2026 (5,200 of them in Europe PMC too); 500
records only Europe PMC finds and 100 preprints, 30 of them with the exact title of a PubMed record (the
preprint of a paper); for each method paper 40 citing records, most of them already found, one without a
core record; one PMID whose efetch returns nothing. With every source merged the list holds more than
5,000 records, so the sample is drawn. Every request is checked: PubMed's term must be the protocol's
query, Europe PMC's search the protocol's or a lookup's or a batch of ids."""
from __future__ import annotations

import json
import re
import urllib.parse
import xml.etree.ElementTree as ET
from datetime import date, timedelta

from search import EUROPEPMC_QUERY, PAPERS, PUBMED_TERM

N_PUBMED, N_BOTH, N_EPMC_ONLY, N_PREPRINT, N_PREPRINT_TWIN, N_CITING = 5600, 5200, 500, 100, 30, 40
MISSING_PMID = "31000077"


def _day(i: int) -> date:
    return date(2017, 1, 1) + timedelta(days=(i * 37) % 3560)


class Fabricated:
    def __init__(self, log, papers: dict):
        self.log, self.rows, self.papers = log, [], papers
        self.pm = [dict(pmid=str(31000000 + 7 * i), doi=f"10.5555/fab.{i}", day=_day(i),
                        title=f"Fabricated record {i}: gene co-expression in single cells, condition {i % 9}",
                        abstract=f"Background: fabricated abstract {i}. Results: the coexpression of GENE{i}A and "
                                 f"GENE{i}B changes between conditions.", journal=["Nat Commun", "Cell Rep", "Elife"][i % 3])
                   for i in range(N_PUBMED)]
        self.epmc = [self._core(r, "MED") for r in self.pm[:N_BOTH]]
        self.epmc += [self._core(dict(pmid=str(39000000 + i), doi=f"10.5555/epmc.{i}", day=_day(i + 11),
                                      title=f"Europe PMC only record {i}: decoupling of transcriptional coordination",
                                      abstract=f"Fabricated abstract of Europe PMC record {i}.", journal="Sci Rep"), "MED")
                      for i in range(N_EPMC_ONLY)]
        for j in range(N_PREPRINT):
            twin = self.pm[j * 3] if j < N_PREPRINT_TWIN else None
            self.epmc.append(dict(id=f"PPR{800000 + j}", source="PPR", doi=f"10.1101/fab.{j}",
                                  title=twin["title"] if twin else f"Preprint {j}: network rewiring in single nuclei",
                                  abstractText=f"Preprint abstract {j}.", pubYear=str(_day(j).year),
                                  firstPublicationDate=_day(j).isoformat(), pubTypeList=dict(pubType=["Preprint"]),
                                  language="eng", journalInfo=None, bookOrReportDetails=dict(publisher="bioRxiv")))
        self.method = {}
        # the fabricated truth is the repository's method_papers.json, whatever the run is told to expect
        truth = {p["key"]: p["expect"] for p in json.loads(PAPERS.read_text())["papers"]}
        for k, p in enumerate(papers["papers"]):
            e = truth.get(p["key"], p["expect"])
            pmid = str(30000000 + k)
            title = " ".join(e.get("title_words", [])) or f"Method paper {p['key']}"
            self.method[p["key"]] = self._core(dict(pmid=pmid, doi=e.get("doi", f"10.5555/method.{k}"), day=date(e["year"], 6, 1),
                                                    title=f"Fabricated {title}", abstract="Method.",
                                                    journal=e["journal"]), "MED", cited_by=N_CITING)
        self.citing = {}
        for k, p in enumerate(papers["papers"]):
            rows = [dict(id=self.pm[(k * 400 + 13 * i) % N_PUBMED]["pmid"], source="MED") for i in range(30)]
            rows += [dict(id=str(38000000 + 100 * k + i), source="MED") for i in range(8)]
            rows += [dict(id=f"PPR{900000 + 10 * k + i}", source="PPR") for i in range(2)]
            self.citing[p["key"]] = rows
        self.citing_only = {}
        for k, p in enumerate(papers["papers"]):
            for i in range(8):
                pmid = str(38000000 + 100 * k + i)
                self.citing_only[("MED", pmid)] = self._core(dict(pmid=pmid, doi=f"10.5555/cite.{k}.{i}", day=_day(k + i),
                                                                  title=f"Citing record {k}.{i}: differential correlation",
                                                                  abstract="Citing.", journal="Bioinformatics"), "MED")
            self.citing_only[("PPR", f"PPR{900000 + 10 * k}")] = dict(
                id=f"PPR{900000 + 10 * k}", source="PPR", doi=f"10.1101/cite.{k}", title=f"Citing preprint {k}",
                abstractText="Citing preprint.", pubYear="2024", pubTypeList=dict(pubType=["Preprint"]), language="eng")
            # PPR{900000 + 10k + 1} has no core record: the citation's own fields are kept
        self.by_key = {("MED", r["pmid"]): r for r in self.epmc if r.get("source") == "MED"}
        self.by_key.update({("PPR", r["id"]): r for r in self.epmc if r.get("source") == "PPR"})
        self.by_key.update(self.citing_only)
        self.pm_by_id = {r["pmid"]: r for r in self.pm}

    @staticmethod
    def _core(r, source, cited_by=0):
        return dict(id=r["pmid"], source=source, pmid=r["pmid"], doi=r["doi"], title=r["title"], abstractText=r["abstract"],
                    pubYear=str(r["day"].year), firstPublicationDate=r["day"].isoformat(),
                    journalInfo=dict(journal=dict(isoabbreviation=r["journal"], title=r["journal"])),
                    pubTypeList=dict(pubType=["research-article", "Journal Article"]), language="eng", citedByCount=cited_by)

    # ------------------------------------------------------------------ the one entry point
    def get(self, url: str, params: dict, post: bool = False) -> bytes:
        self.rows.append(dict(time="fabricated", method="POST" if post else "GET", url=url,
                              params=json.dumps({k: str(v)[:80] for k, v in params.items()}), status=200, bytes=0,
                              sha256="-", attempt=1))
        path = urllib.parse.urlsplit(url).path
        if path.endswith("/esearch.fcgi"):
            return self._esearch(params)
        if path.endswith("/efetch.fcgi"):
            assert post, "efetch goes by POST"
            return self._efetch(params)
        if path.endswith("/search"):
            return self._epmc_search(params)
        m = re.search(r"/(MED|PPR|PMC|AGR)/([^/]+)/citations$", path)
        if m:
            return self._citations(m.group(1), m.group(2), params)
        raise AssertionError(f"no fabricated response for {url}")

    def _esearch(self, params):
        assert params["term"] == PUBMED_TERM, "PubMed's term must be the protocol's query"
        assert params.get("tool") and "email" not in params
        if "mindate" in params:
            a = date(*map(int, params["mindate"].split("/")))
            b = date(*map(int, params["maxdate"].split("/")))
            ids = [r["pmid"] for r in self.pm if a <= r["day"] <= b]
        else:
            ids = [r["pmid"] for r in self.pm]
        n = int(params["retmax"])
        res = dict(count=str(len(ids)), retmax=str(min(n, len(ids))), retstart="0", idlist=ids[:n],
                   querytranslation="(fabricated translation)")
        return json.dumps(dict(header=dict(type="esearch"), esearchresult=res)).encode()

    def _efetch(self, params):
        root = ET.Element("PubmedArticleSet")
        for pmid in params["id"].split(","):
            r = self.pm_by_id.get(pmid)
            if r is None or pmid == MISSING_PMID:
                continue
            art = ET.SubElement(root, "PubmedArticle")
            mc = ET.SubElement(art, "MedlineCitation")
            ET.SubElement(mc, "PMID").text = pmid
            a = ET.SubElement(mc, "Article")
            j = ET.SubElement(a, "Journal")
            ji = ET.SubElement(j, "JournalIssue")
            pd = ET.SubElement(ji, "PubDate")
            ET.SubElement(pd, "Year").text = str(r["day"].year)
            ET.SubElement(j, "Title").text = r["journal"] + " (full title)"
            ET.SubElement(j, "ISOAbbreviation").text = r["journal"]
            t = ET.SubElement(a, "ArticleTitle")
            t.text = r["title"].split(":")[0] + ":"
            ET.SubElement(t, "i").text = " gene"
            t[0].tail = r["title"].split(":", 1)[1].replace(" gene", "", 1)
            ab = ET.SubElement(a, "Abstract")
            for part in r["abstract"].split(". "):
                lab, _, txt = part.partition(": ")
                at = ET.SubElement(ab, "AbstractText", Label=lab.upper()) if txt else ET.SubElement(ab, "AbstractText")
                at.text = txt or part
            ET.SubElement(a, "Language").text = "eng"
            ptl = ET.SubElement(a, "PublicationTypeList")
            ET.SubElement(ptl, "PublicationType").text = "Journal Article"
            pdat = ET.SubElement(art, "PubmedData")
            idl = ET.SubElement(pdat, "ArticleIdList")
            ET.SubElement(idl, "ArticleId", IdType="pubmed").text = pmid
            ET.SubElement(idl, "ArticleId", IdType="doi").text = r["doi"].upper()
        return ET.tostring(root, encoding="utf-8", xml_declaration=True)

    def _epmc_search(self, params):
        q = params["query"]
        size = int(params.get("pageSize", 25))
        if q == EUROPEPMC_QUERY:
            start = 0 if params.get("cursorMark") == "*" else int(params["cursorMark"])
            page = self.epmc[start:start + size]
            nxt = str(start + size) if start + size < len(self.epmc) else params.get("cursorMark")
            return json.dumps(dict(hitCount=len(self.epmc), nextCursorMark=nxt, resultList=dict(result=page))).encode()
        for p in self.papers["papers"]:
            if q == p["query"]:
                rec = self.method[p["key"]]
                decoy = dict(rec, id=str(int(rec["id"]) + 1), pmid=str(int(rec["id"]) + 1), pubTypeList=dict(pubType=["Erratum"]))
                preprint = dict(rec, id=f"PPR7{rec['id'][-3:]}", source="PPR", pmid="", pubTypeList=dict(pubType=["Preprint"]))
                return json.dumps(dict(hitCount=3, resultList=dict(result=[preprint, rec, decoy]))).encode()
        pairs = re.findall(r"EXT_ID:(\S+) AND SRC:(\w+)", q)
        assert pairs, f"no fabricated search for {q[:100]}"
        res = [self.by_key[(src, cid)] for cid, src in pairs if (src, cid) in self.by_key]
        return json.dumps(dict(hitCount=len(res), resultList=dict(result=res[:size]))).encode()

    def _citations(self, src, cid, params):
        key = next(k for k, r in self.method.items() if r["id"] == cid and r["source"] == src)
        rows = self.citing[key]
        page, size = int(params["page"]), int(params["pageSize"])
        size = min(size, 16)  # pages smaller than asked, so that the paging runs
        chunk = rows[(page - 1) * size:page * size]
        cits = [dict(c, title=f"Citation {c['id']}", journalAbbreviation="J Fab", pubYear=2023) for c in chunk]
        return json.dumps(dict(hitCount=len(rows), citationList=dict(citation=cits))).encode()
