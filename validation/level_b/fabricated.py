"""Fabricated responses for the dry run of `search.py` (nothing is fetched): PubMed's esearch and efetch,
Europe PMC's search and citations endpoint, as `search.Net.get` would return them.

The fabricated world, fixed: 5,600 PubMed records over 2017–2026 (5,200 of them in Europe PMC too); 500
records only Europe PMC finds and 100 preprints, 30 of them with the exact title of a PubMed record (the
preprint of a paper); for each method paper 40 citing records, most of them already found, one without a
core record; one PMID whose efetch returns nothing. With every source merged the list holds more than
5,000 records, so the sample is drawn. Every request is checked: PubMed's term must be the protocol's
query, Europe PMC's search the protocol's or a lookup's or a batch of ids. For the full texts: every third
PubMed record that Europe PMC holds has a PMCID in its core record (but not in PubMed's), and the full text
of every other such PMCID is open (PubMed gives these PMCIDs too), but one answers a server error on every
attempt; the even preprints' full texts are open.
Every fifth PubMed record lists its references with their ids, as PubMed does: the next record's DOI and PMID
and a PMCID of their own (a parser that took them for the record's own would merge the two). Each response
is recorded as `search.Net` records it, with its SHA256, so that a dry run can be replayed (`remerge`).
Amendment 1's second fetch and preprints (`plan_s1.md`): a DOI whose last character is an even digit opens
on a page that names its PDF (open), any other DOI is refused (HTTP 403); NCBI's links from a paper to GEO
are none for an odd PMID and one data set for an even one; the supplementary files of a PMCID are a small
zip; NGDC's BioProject page names its accession; the preprints' query answers with the fabricated
preprints, without a language (as Europe PMC gives them)."""
from __future__ import annotations

import hashlib
import io
import json
import re
import zipfile
import urllib.error
import urllib.parse
import xml.etree.ElementTree as ET
from datetime import date, timedelta

from search import EUROPEPMC_QUERY, PAPERS, PREPRINT_QUERY, PUBMED_TERM, Unavailable, shown_params

N_PUBMED, N_BOTH, N_EPMC_ONLY, N_PREPRINT, N_PREPRINT_TWIN, N_CITING = 5600, 5200, 500, 100, 30, 40
MISSING_PMID = "31000077"
SERVER_ERROR_FULLTEXT = "PMC7000006"  # open, but its full text answers 500 on every attempt


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
        self.epmc = [self._core(dict(r, pmcid=f"PMC{7000000 + i}" if i % 3 == 0 else ""), "MED")
                     for i, r in enumerate(self.pm[:N_BOTH])]
        self.open_fulltext = {f"PMC{7000000 + i}" for i in range(0, N_BOTH, 6)}
        self.open_fulltext |= {f"PPR{800000 + j}" for j in range(0, N_PREPRINT, 2)}
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
        core = dict(id=r["pmid"], source=source, pmid=r["pmid"], doi=r["doi"], title=r["title"], abstractText=r["abstract"],
                    pubYear=str(r["day"].year), firstPublicationDate=r["day"].isoformat(),
                    journalInfo=dict(journal=dict(isoabbreviation=r["journal"], title=r["journal"])),
                    pubTypeList=dict(pubType=["research-article", "Journal Article"]), language="eng", citedByCount=cited_by)
        if r.get("pmcid"):
            core.update(pmcid=r["pmcid"], fullTextIdList=dict(fullTextId=[r["pmcid"]]))
        return core

    # ------------------------------------------------------------------ the one entry point
    def get(self, url: str, params: dict, post: bool = False) -> bytes:
        self.rows.append(dict(time="fabricated", method="POST" if post else "GET", url=url, params=shown_params(params),
                              status=200, bytes=0, sha256=hashlib.sha256(b"").hexdigest(), attempt=1))
        body = self._respond(url, params, post)
        self.rows[-1].update(bytes=len(body), sha256=hashlib.sha256(body).hexdigest())
        return body

    def get_once(self, url: str) -> tuple:
        """`search.Net.get_once`: one GET, (status, final URL, body), an error answered, not raised."""
        status, final, body = self._once(url)
        self.rows.append(dict(time="fabricated", method="GET", url=url, params=shown_params({}), status=status,
                              bytes=len(body), sha256=hashlib.sha256(body).hexdigest(), attempt=1))
        return status, final, body

    @staticmethod
    def _once(url):
        parts = urllib.parse.urlsplit(url)
        if parts.hostname == "doi.org":
            doi = parts.path.lstrip("/")
            if doi[-1:].isdigit() and int(doi[-1]) % 2 == 0:
                final = f"https://publisher.example/article/{doi}"
                return 200, final, (f'<html><head><meta name="citation_pdf_url" content="/pdf/{doi}.pdf">'
                                    f"</head><body><h1>Fabricated article {doi}</h1><h2>Results</h2></body></html>").encode()
            return 403, f"https://publisher.example/article/{doi}", b"<html>Access denied</html>"
        if parts.hostname == "publisher.example" and parts.path.startswith("/pdf/"):
            return 200, url, b"%PDF-1.4 fabricated"
        raise AssertionError(f"no fabricated one-attempt response for {url}")

    def _respond(self, url, params, post):
        path = urllib.parse.urlsplit(url).path
        if path.endswith("/elink.fcgi"):
            return self._elink(params)
        if path.endswith("/supplementaryFiles"):
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w") as z:
                z.writestr("Table1.docx", "fabricated supplementary table")
            return buf.getvalue()
        if urllib.parse.urlsplit(url).hostname == "ngdc.cncb.ac.cn":
            acc = path.rsplit("/", 1)[-1]
            return f"<html><body>BioProject {acc}: fabricated page, no data released</body></html>".encode()
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
        m = re.search(r"/((?:PMC|PPR)\d+)/fullTextXML$", path)
        if m:
            return self._fulltext(url, m.group(1))
        raise AssertionError(f"no fabricated response for {url}")

    def _fulltext(self, url, fid):
        if fid == SERVER_ERROR_FULLTEXT:  # as search.Net raises after six attempts answered 500
            self.rows[-1]["status"] = 500
            raise Unavailable(f"the request failed six times (HTTP 500): {url}", 500)
        if fid not in self.open_fulltext:
            self.rows[-1]["status"] = 404
            raise urllib.error.HTTPError(url, 404, "Not Found", None, None)
        return (f'<?xml version="1.0"?><article><front><article-meta><article-id pub-id-type="pmcid">{fid}</article-id>'
                f'<title-group><article-title>Fabricated full text {fid}</article-title></title-group>'
                '<abstract><p>An abstract.</p></abstract></article-meta></front><body>'
                '<sec><title>Results</title><p>The correlation of GENE1A and GENE1B fell with age '
                '(<xref ref-type="fig" rid="f1">Fig. 1</xref>).</p>'
                '<fig id="f1"><label>Figure 1</label><caption><p>The pair by age.</p></caption></fig></sec>'
                '<sec><title>Methods</title><p>Eight mice per group; counts on GEO.</p></sec></body>'
                '<back><ref-list><ref><mixed-citation>A reference.</mixed-citation></ref></ref-list></back>'
                '</article>').encode()

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
            i = (int(pmid) - 31000000) // 7
            if i % 6 == 0 and i < N_BOTH:  # PubMed gives the PMCID of every other record whose core record has one
                ET.SubElement(idl, "ArticleId", IdType="pmc").text = f"PMC{7000000 + i}"
            if i % 5 == 1:  # the references, each with its own ids, after the record's
                refs = ET.SubElement(pdat, "ReferenceList")
                nxt = self.pm[(i + 1) % N_PUBMED]
                for k, (rdoi, rpmid) in enumerate([("10.9999/unrelated", "20000001"), (nxt["doi"], nxt["pmid"])]):
                    ref = ET.SubElement(refs, "Reference")
                    ET.SubElement(ref, "Citation").text = f"Reference {k + 1}."
                    ril = ET.SubElement(ref, "ArticleIdList")
                    ET.SubElement(ril, "ArticleId", IdType="doi").text = rdoi
                    ET.SubElement(ril, "ArticleId", IdType="pubmed").text = rpmid
                    ET.SubElement(ril, "ArticleId", IdType="pmc").text = f"PMC{5000000 + 10 * i + k}"
        return ET.tostring(root, encoding="utf-8", xml_declaration=True)

    @staticmethod
    def _elink(params):
        assert params["dbfrom"] == "pubmed" and params["db"] == "gds" and "email" not in params
        pmid = params["id"]
        dbs = [dict(dbto="gds", linkname="pubmed_gds", links=[str(200000000 + int(pmid) % 1000)])] if int(pmid) % 2 == 0 else []
        return json.dumps(dict(header=dict(type="elink"), linksets=[dict(dbfrom="pubmed", ids=[pmid], linksetdbs=dbs)])).encode()

    def _epmc_search(self, params):
        q = params["query"]
        size = int(params.get("pageSize", 25))
        if q == PREPRINT_QUERY:
            pre = [dict(r, language=None) for r in self.epmc if r.get("source") == "PPR"]
            start = 0 if params.get("cursorMark") == "*" else int(params["cursorMark"])
            size = min(size, 40)  # pages smaller than asked, so that the paging runs
            page = pre[start:start + size]
            nxt = str(start + size) if start + size < len(pre) else params.get("cursorMark")
            return json.dumps(dict(hitCount=len(pre), nextCursorMark=nxt, resultList=dict(result=page))).encode()
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
