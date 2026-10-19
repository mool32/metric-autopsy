#!/usr/bin/env python3
"""Level B's search for claims (`validation/plans/level_b_protocol.md`, approved by the owner on 2026-10-19).

Two runs, each once and by its own label in GitHub Actions (this environment reaches neither database):

- ``lookup``: the method papers whose forward citations are searched (`method_papers.json`), each found in
  Europe PMC by its own query and checked against what the protocol expects (its journal's ISO
  abbreviation, its year, its title's words or its DOI), exactly once among the journal articles. The
  records are written; their Europe PMC source and id are then committed into `method_papers.json`
  ("verified") before the search runs.
- ``search``: PubMed through E-utilities (``esearch`` by date windows of at most 9,999 ids, a window split
  in two while it holds more, then ``efetch`` by PMID, 200 at a time), Europe PMC's REST search (cursor
  pages of 1,000, every source, preprints included), and each verified method paper's forward citations
  (Europe PMC's citations endpoint, then each citing record's core record, abstract included, 25 at a
  time). Every response is written as it came (gzip) with every request's URL, parameters, status, size
  and SHA256 (``requests.tsv``). The records are then merged (the same PMID, else DOI, else title once
  lower-cased and stripped to letters and digits, is one record), numbered in a fixed order (PMID, then
  DOI, then title) and, if more than 5,000, a sample of 5,000 is drawn for screening with the protocol's
  seed: ``numpy.random.default_rng(20261019).choice(n, 5000, replace=False)``.

- ``fulltext``: after the screening of titles and abstracts (`screening.md`), the full texts of the records
  either pass kept or marked unclear (the committed list `fulltext_list.tsv`): Europe PMC's open-access
  full text (JATS XML, ``/{PMCID}/fullTextXML``), the PMCID taken from the record or, where it has none,
  from its core record in Europe PMC (by PMID, else by its Europe PMC ids; a preprint's own id is tried
  too). A record without an accessible full text is listed so (the protocol: it goes to the owner).

- ``remerge``: a finished search's records rebuilt from its own raw responses, nothing fetched: the search's
  code runs again with each request answered by the response the run recorded for it (``requests.tsv`` and
  ``raw/``, matched by URL, parameters and SHA256), so a fix of the parsing or the merging takes effect on
  the same data. Written to a new directory; the run's own files are left as they are.

``--fabricated`` replaces the network with fabricated responses (the dry run): every path of the code runs
and its outputs are checked; nothing is fetched. NCBI's ``tool`` parameter names this script; no e-mail is
sent.

    python -I search.py lookup --out OUT [--fabricated]
    python -I search.py search --out OUT [--papers method_papers.json] [--fabricated]
    python -I search.py fulltext --out OUT [--list fulltext_list.tsv] [--fabricated]
    python -I search.py remerge --search RESULTS --out OUT [--papers method_papers.json]
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import platform
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPERS = HERE / "method_papers.json"
FULLTEXT_LIST = HERE / "fulltext_list.tsv"
FULLTEXT_COLUMNS = ["rid", "pmid", "pmcid", "doi", "epmc_ids"]

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
TOOL = "metric-autopsy-level-b"
HEADERS = {"User-Agent": "metric-autopsy level B search (https://github.com/mool32/metric-autopsy)"}
PAUSE = {"eutils.ncbi.nlm.nih.gov": 0.4, "www.ebi.ac.uk": 0.2}  # NCBI: 3 requests a second without a key

FIRST, LAST = date(2017, 1, 1), date(2026, 9, 30)
WINDOW_LIMIT = 9999  # the most ids PubMed's esearch returns for one query
EFETCH_BATCH = 200
CITING_BATCH = 25
SAMPLE_SEED, SAMPLE_SIZE = 20261019, 5000

# the protocol's section 3, the PubMed block as approved (a test checks it against the protocol)
PUBMED_QUERY = """("single-cell"[tiab] OR "single cell"[tiab] OR "single-nucleus"[tiab] OR "single nucleus"[tiab]
 OR scRNA-seq[tiab] OR snRNA-seq[tiab])
AND
(coexpression[tiab] OR co-expression[tiab] OR "gene-gene correlation*"[tiab] OR "gene coordination"[tiab]
 OR "transcriptional coordination"[tiab] OR "coordinated expression"[tiab] OR decoupling[tiab]
 OR "loss of coordination"[tiab] OR "network rewiring"[tiab] OR "differential coexpression"[tiab]
 OR "co-regulation"[tiab] OR "transcriptional coherence"[tiab] OR "differential correlation"[tiab])
AND English[lang] AND ("2017/01/01"[dp] : "2026/09/30"[dp])"""
PUBMED_TERM = " ".join(PUBMED_QUERY.split())

# Europe PMC: the same two blocks, each term in the title or the abstract; no truncation inside a phrase
SINGLE_CELL = ["single-cell", "single cell", "single-nucleus", "single nucleus", "scRNA-seq", "snRNA-seq"]
COORDINATION = ["coexpression", "co-expression", "gene-gene correlation", "gene-gene correlations",
                "gene coordination", "transcriptional coordination", "coordinated expression", "decoupling",
                "loss of coordination", "network rewiring", "differential coexpression", "co-regulation",
                "transcriptional coherence", "differential correlation"]


def _epmc_block(terms) -> str:
    return "(" + " OR ".join(f'TITLE:"{t}" OR ABSTRACT:"{t}"' for t in terms) + ")"


EUROPEPMC_QUERY = (f'{_epmc_block(SINGLE_CELL)} AND {_epmc_block(COORDINATION)} AND LANG:"eng" '
                   f"AND FIRST_PDATE:[{FIRST.isoformat()} TO {LAST.isoformat()}]")


# --------------------------------------------------------------------------- #
# the network, recorded
# --------------------------------------------------------------------------- #
def shown_params(params: dict) -> str:
    """A request's parameters as `requests.tsv` records them (a long value cut, with its length)."""
    shown = {k: (v if len(str(v)) <= 120 else f"{str(v)[:60]}… ({len(str(v))} characters)") for k, v in params.items()}
    return json.dumps(shown, ensure_ascii=False)


class Net:
    """GET and POST with retries on network errors and on 429 and 5xx; every request recorded."""

    def __init__(self, log):
        self.log, self.rows, self.last = log, [], {}

    def get(self, url: str, params: dict, post: bool = False) -> bytes:
        host = urllib.parse.urlsplit(url).hostname or ""
        query = urllib.parse.urlencode(params, doseq=True)
        last_error = None
        for attempt in range(6):
            wait = PAUSE.get(host, 0.2) - (time.monotonic() - self.last.get(host, 0.0))
            if wait > 0:
                time.sleep(wait)
            self.last[host] = time.monotonic()
            req = (urllib.request.Request(url, data=query.encode(), headers=HEADERS) if post
                   else urllib.request.Request(f"{url}?{query}", headers=HEADERS))
            try:
                with urllib.request.urlopen(req, timeout=180) as r:
                    body, status = r.read(), r.status
                self._record(url, params, post, status, body, attempt)
                return body
            except urllib.error.HTTPError as e:
                last_error = f"HTTP {e.code}"
                if e.code not in (429, 500, 502, 503, 504):
                    self._record(url, params, post, e.code, b"", attempt)
                    raise
            except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
                last_error = f"{type(e).__name__}: {e}"
            self.log(f"  retry {attempt + 1} after {last_error}: {url}")
            time.sleep(min(60, 2 ** (attempt + 1)))
        raise RuntimeError(f"the request failed six times ({last_error}): {url}")

    def _record(self, url, params, post, status, body, attempt):
        self.rows.append(dict(time=datetime.now(timezone.utc).isoformat(timespec="seconds"), method="POST" if post else "GET",
                              url=url, params=shown_params(params), status=status, bytes=len(body),
                              sha256=hashlib.sha256(body).hexdigest(), attempt=attempt + 1))


class Replay:
    """The responses a finished run recorded, in its order: each request must be the next one of its
    requests.tsv (method, URL and parameters as recorded), and gets the raw response whose SHA256 the run
    recorded for it. Nothing is fetched."""

    def __init__(self, log, root: Path):
        self.log, self.rows, self.last = log, [], {}
        self.bodies = {}
        for p in sorted((root / "raw").rglob("*.gz")):
            with gzip.open(p) as g:
                body = g.read()
            self.bodies[hashlib.sha256(body).hexdigest()] = body
        lines = (root / "requests.tsv").read_text().splitlines()
        cols = lines[0].split("\t")
        self.recorded = [dict(zip(cols, line.split("\t"))) for line in lines[1:] if line.strip()]
        self.i = 0

    def get(self, url: str, params: dict, post: bool = False) -> bytes:
        if self.i >= len(self.recorded):
            raise RuntimeError(f"the run recorded no more requests; asked for {url}")
        rec = self.recorded[self.i]
        self.i += 1
        shown = shown_params(params)
        if (rec["method"], rec["url"], rec["params"]) != ("POST" if post else "GET", url, shown):
            raise RuntimeError(f"request {self.i} is not the run's: {url} {shown}, recorded {rec['url']} {rec['params']}")
        if rec["status"] != "200":
            raise urllib.error.HTTPError(url, int(rec["status"]), "as recorded", None, None)
        body = self.bodies.get(rec["sha256"])
        if body is None:
            raise RuntimeError(f"request {self.i}: no raw response with the recorded SHA256 {rec['sha256']}")
        self.rows.append(dict(rec))
        return body


# --------------------------------------------------------------------------- #
# outputs
# --------------------------------------------------------------------------- #
class Out:
    def __init__(self, root: Path, write_raw: bool = True):
        self.root, self.write_raw = root, write_raw
        root.mkdir(parents=True, exist_ok=True)
        self.lines = []

    def log(self, line: str = ""):
        self.lines.append(line)
        print(line, flush=True)

    def raw(self, rel: str, body: bytes) -> str:
        p = self.root / "raw" / f"{rel}.gz"
        if not self.write_raw:  # a replay: the run's raw responses are where they were
            return str(p.relative_to(self.root))
        p.parent.mkdir(parents=True, exist_ok=True)
        with gzip.GzipFile(filename="", mode="wb", fileobj=p.open("wb"), mtime=0) as g:
            g.write(body)
        return str(p.relative_to(self.root))

    def finish(self, name: str, net: Net, summary: dict):
        with (self.root / "requests.tsv").open("w") as f:
            cols = ["time", "method", "url", "params", "status", "bytes", "sha256", "attempt"]
            f.write("\t".join(cols) + "\n")
            for r in net.rows:
                f.write("\t".join(str(r[c]) for c in cols) + "\n")
        (self.root / f"{name}.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False) + "\n")
        (self.root / f"{name}.log").write_text("\n".join(self.lines) + "\n")


def _header(out: Out, mode: str, fabricated: bool):
    out.log(f"# Level B, the {mode} — {datetime.now(timezone.utc).isoformat(timespec='seconds')}"
            f"{' (FABRICATED responses: the dry run, nothing fetched)' if fabricated else ''}")
    out.log(f"python {platform.python_version()}; the protocol validation/plans/level_b_protocol.md (approved "
            "2026-10-19); this script validation/level_b/search.py")
    out.log()


# --------------------------------------------------------------------------- #
# parsing
# --------------------------------------------------------------------------- #
def _text(el) -> str:
    return " ".join("".join(el.itertext()).split()) if el is not None else ""


def parse_pubmed_xml(body: bytes) -> list:
    """PubMed's efetch XML: one record per article (PMID, title, abstract, journal, year, DOI, PMCID,
    publication types, languages)."""
    root = ET.fromstring(body)
    out = []
    for art in root.iter("PubmedArticle"):
        mc, pd = art.find("MedlineCitation"), art.find("PubmedData")
        a = mc.find("Article") if mc is not None else None
        if mc is None or a is None:
            continue
        # the article's own ids only: PubmedData/ReferenceList holds the references' ArticleIds too
        ids = {i.get("IdType"): (i.text or "").strip() for i in (pd.findall("ArticleIdList/ArticleId") if pd is not None else [])}
        doi = ids.get("doi") or next((e.text for e in a.iter("ELocationID") if e.get("EIdType") == "doi"), "")
        parts = []
        for t in a.iter("AbstractText"):
            label = t.get("Label")
            parts.append(f"{label}: {_text(t)}" if label else _text(t))
        year = (a.findtext("Journal/JournalIssue/PubDate/Year") or (a.findtext("Journal/JournalIssue/PubDate/MedlineDate")
                or "")[:4] or a.findtext("ArticleDate/Year") or "")
        out.append(dict(pmid=(mc.findtext("PMID") or "").strip(), doi=(doi or "").strip().lower(),
                        pmcid=ids.get("pmc", ""), title=_text(a.find("ArticleTitle")), abstract="\n".join(parts),
                        journal=a.findtext("Journal/ISOAbbreviation") or a.findtext("Journal/Title") or "",
                        year=int(year) if year.isdigit() else None,
                        pub_types=[_text(p) for p in a.iter("PublicationType")],
                        languages=[(x.text or "") for x in a.iter("Language")], source="MED"))
    for book in root.iter("PubmedBookArticle"):
        doc = book.find("BookDocument")
        if doc is None:
            continue
        year = doc.findtext("Book/PubDate/Year") or ""
        out.append(dict(pmid=(doc.findtext("PMID") or "").strip(), doi="", pmcid="",
                        title=_text(doc.find("ArticleTitle")) or _text(doc.find("Book/BookTitle")),
                        abstract="\n".join(_text(t) for t in doc.iter("AbstractText")), journal="(book)",
                        year=int(year) if year.isdigit() else None, pub_types=["Book"], languages=[], source="MED"))
    return out


_TAG = re.compile(r"<[^>]+>")


def parse_epmc_record(r: dict) -> dict:
    """A Europe PMC core record, its fields as the merged list keeps them."""
    j = (r.get("journalInfo") or {}).get("journal") or {}
    year = r.get("pubYear") or ""
    types = (r.get("pubTypeList") or {}).get("pubType") or []
    return dict(pmid=str(r.get("pmid") or ""), doi=str(r.get("doi") or "").lower(), pmcid=str(r.get("pmcid") or ""),
                title=" ".join(_TAG.sub(" ", r.get("title") or "").split()),
                abstract=" ".join(_TAG.sub(" ", r.get("abstractText") or "").split()),
                journal=(j.get("isoabbreviation") or j.get("medlineAbbreviation") or j.get("title")
                         or (r.get("bookOrReportDetails") or {}).get("publisher") or ""),
                year=int(year) if str(year).isdigit() else None,
                pub_types=types if isinstance(types, list) else [types], languages=[r.get("language") or ""],
                source=str(r.get("source") or ""), epmc_id=f"{r.get('source')}:{r.get('id')}",
                first_date=r.get("firstPublicationDate") or "", cited_by=r.get("citedByCount"))


def norm_title(t: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (t or "").lower())


def _norm_journal(j: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (j or "").lower())


# --------------------------------------------------------------------------- #
# the lookup of the method papers
# --------------------------------------------------------------------------- #
NOT_ARTICLES = ("erratum", "correction", "comment", "retraction", "editorial", "news", "published erratum")


def lookup(net, out: Out, papers: dict) -> dict:
    """Each method paper, found by its query and checked: exactly one journal article (source MED, not an
    erratum or a comment) with the expected journal, year and title words (or DOI)."""
    results, ok = [], True
    for p in papers["papers"]:
        body = net.get(f"{EPMC}/search", dict(query=p["query"], resultType="core", format="json", pageSize=25))
        rel = out.raw(f"lookup/{p['key']}.json", body)
        recs = [parse_epmc_record(r) for r in (json.loads(body).get("resultList") or {}).get("result", [])]
        e = p["expect"]

        def matches(r):
            if r["source"] != "MED" or any(t.lower() in NOT_ARTICLES for t in r["pub_types"]):
                return False
            if _norm_journal(r["journal"]) != _norm_journal(e["journal"]) or r["year"] != e["year"]:
                return False
            if e.get("doi") and r["doi"] != e["doi"].lower():
                return False
            return all(w.lower() in r["title"].lower() for w in e.get("title_words", []))
        hits = [r for r in recs if matches(r)]
        found = len(hits) == 1
        ok &= found
        rec = hits[0] if found else None
        out.log(f"- {p['key']} ({p['label']}): query {p['query']!r} -> {len(recs)} records, {len(hits)} matching "
                f"{e['journal']} {e['year']}" + (f": {rec['epmc_id']} PMID {rec['pmid']} doi {rec['doi']} — "
                                                 f"{rec['title']} ({rec['journal']} {rec['year']}; cited by "
                                                 f"{rec['cited_by']})" if rec else " — NOT VERIFIED"))
        for r in recs[:10]:
            out.log(f"    {r['epmc_id']} {r['journal']} {r['year']} {r['pub_types']}: {r['title'][:140]}")
        results.append(dict(key=p["key"], label=p["label"], query=p["query"], expect=e, raw=rel, n_records=len(recs),
                            n_matching=len(hits), verified=dict(source=rec["source"], id=rec["epmc_id"].split(":", 1)[1],
                                                                pmid=rec["pmid"], doi=rec["doi"], title=rec["title"],
                                                                journal=rec["journal"], year=rec["year"],
                                                                cited_by=rec["cited_by"]) if rec else None))
    out.log()
    out.log("All the method papers verified." if ok else "NOT ALL the method papers verified: the search must not run "
            "until method_papers.json is corrected.")
    return dict(ok=bool(ok), papers=results)


# --------------------------------------------------------------------------- #
# the search
# --------------------------------------------------------------------------- #
def _windows(first: date, last: date):
    for y in range(first.year, last.year + 1):
        yield max(first, date(y, 1, 1)), min(last, date(y, 12, 31))


def pubmed(net, out: Out, window_limit: int = WINDOW_LIMIT) -> tuple[list, dict]:
    """The query over the whole window (its count, PubMed's translation and warnings), then by date windows of
    at most `window_limit` ids, each split in two while it holds more, then the records by PMID."""
    common = dict(db="pubmed", term=PUBMED_TERM, retmode="json", tool=TOOL)
    body = net.get(f"{EUTILS}/esearch.fcgi", dict(common, retmax=0))
    out.raw("pubmed/esearch-total.json", body)
    res = json.loads(body).get("esearchresult", {})
    total = int(res.get("count", 0))
    out.log(f"PubMed: {total} records for the query; translation: {res.get('querytranslation')}")
    for k in ("errorlist", "warninglist"):
        if res.get(k):
            out.log(f"  PubMed's {k}: {json.dumps(res[k])}")
    if res.get("errorlist") and any(res["errorlist"].get(x) for x in ("phrasesnotfound", "fieldsnotfound")):
        out.log("  (a phrase or field PubMed does not know; the search goes on, the log keeps it)")
    pmids, windows, stack = [], [], list(_windows(FIRST, LAST))[::-1]
    while stack:
        a, b = stack.pop()
        body = net.get(f"{EUTILS}/esearch.fcgi", dict(common, datetype="pdat", mindate=a.strftime("%Y/%m/%d"),
                                                     maxdate=b.strftime("%Y/%m/%d"), retmax=window_limit))
        r = json.loads(body).get("esearchresult", {})
        n = int(r.get("count", 0))
        if n > window_limit:
            out.raw(f"pubmed/esearch-{a:%Y%m%d}-{b:%Y%m%d}-split.json", body)
            if a == b:
                raise RuntimeError(f"one day, {a}, holds {n} records, more than {window_limit}")
            mid = a + (b - a) // 2
            stack += [(mid + timedelta(days=1), b), (a, mid)]
            out.log(f"  window {a}..{b}: {n} records, split")
            continue
        out.raw(f"pubmed/esearch-{a:%Y%m%d}-{b:%Y%m%d}.json", body)
        ids = r.get("idlist", [])
        if len(ids) != n:
            raise RuntimeError(f"window {a}..{b}: {n} counted, {len(ids)} returned")
        windows.append(dict(first=a.isoformat(), last=b.isoformat(), count=n))
        pmids += ids
        out.log(f"  window {a}..{b}: {n} records")
    unique = list(dict.fromkeys(pmids))
    out.log(f"  {len(pmids)} ids over {len(windows)} windows, {len(unique)} unique (the whole query counted {total})")
    records = []
    for i in range(0, len(unique), EFETCH_BATCH):
        batch = unique[i:i + EFETCH_BATCH]
        body = net.get(f"{EUTILS}/efetch.fcgi", dict(db="pubmed", id=",".join(batch), retmode="xml", tool=TOOL), post=True)
        out.raw(f"pubmed/efetch-{i // EFETCH_BATCH:04d}.xml", body)
        got = parse_pubmed_xml(body)
        records += got
        missing = set(batch) - {g["pmid"] for g in got}
        if missing:
            out.log(f"  efetch batch {i // EFETCH_BATCH}: {len(missing)} ids without a record: {sorted(missing)[:10]}")
    out.log(f"  {len(records)} records fetched")
    return records, dict(total=total, translation=res.get("querytranslation"), errorlist=res.get("errorlist"),
                         warninglist=res.get("warninglist"), windows=windows, ids=len(pmids), unique=len(unique),
                         fetched=len(records))


def europepmc(net, out: Out) -> tuple[list, dict]:
    """The query in Europe PMC's syntax, every source, cursor pages of 1,000 core records."""
    cursor, page, records, hit = "*", 0, [], None
    while True:
        body = net.get(f"{EPMC}/search", dict(query=EUROPEPMC_QUERY, resultType="core", format="json", pageSize=1000,
                                              cursorMark=cursor))
        out.raw(f"europepmc/search-{page:04d}.json", body)
        d = json.loads(body)
        hit = int(d.get("hitCount", 0)) if hit is None else hit
        res = (d.get("resultList") or {}).get("result", [])
        records += [parse_epmc_record(r) for r in res]
        nxt = d.get("nextCursorMark")
        page += 1
        if not res or not nxt or nxt == cursor:
            break
        cursor = nxt
    by_src = {}
    for r in records:
        by_src[r["source"]] = by_src.get(r["source"], 0) + 1
    out.log(f"Europe PMC: {hit} records for the query; {len(records)} fetched over {page} pages; by source "
            + ", ".join(f"{k} {v}" for k, v in sorted(by_src.items())))
    if len(records) != hit:
        out.log(f"  the pages hold {len(records)} records, the count {hit}")
    return records, dict(hit_count=hit, fetched=len(records), pages=page, by_source=by_src, query=EUROPEPMC_QUERY)


def citations(net, out: Out, papers: dict) -> tuple[list, dict]:
    """Each verified method paper's citing records (the citations endpoint, pages of 1,000), then their core
    records from Europe PMC's search, `CITING_BATCH` at a time; a record the batch misses is asked alone."""
    records, summary = [], {}
    for p in papers["papers"]:
        v = p.get("verified")
        if not v:
            raise RuntimeError(f"{p['key']} is not verified: run the lookup and record its result in method_papers.json")
        listed, page, hit = [], 1, None
        while True:
            body = net.get(f"{EPMC}/{v['source']}/{v['id']}/citations", dict(page=page, pageSize=1000, format="json"))
            out.raw(f"europepmc/citations-{p['key']}-{page:03d}.json", body)
            d = json.loads(body)
            hit = int(d.get("hitCount", 0)) if hit is None else hit
            got = (d.get("citationList") or {}).get("citation", [])
            listed += got
            if not got or len(listed) >= hit:
                break
            page += 1
        keys = list(dict.fromkeys((str(c.get("source")), str(c.get("id"))) for c in listed if c.get("id")))
        cores = {}
        for i in range(0, len(keys), CITING_BATCH):
            batch = keys[i:i + CITING_BATCH]
            q = " OR ".join(f"(EXT_ID:{cid} AND SRC:{src})" for src, cid in batch)
            body = net.get(f"{EPMC}/search", dict(query=q, resultType="core", format="json", pageSize=100))
            out.raw(f"europepmc/citing-{p['key']}-{i // CITING_BATCH:04d}.json", body)
            for r in (json.loads(body).get("resultList") or {}).get("result", []):
                cores[(str(r.get("source")), str(r.get("id")))] = r
        alone = [k for k in keys if k not in cores]
        for src, cid in alone:
            body = net.get(f"{EPMC}/search", dict(query=f"EXT_ID:{cid} AND SRC:{src}", resultType="core", format="json",
                                                  pageSize=10))
            out.raw(f"europepmc/citing-{p['key']}-{src}-{cid}.json", body)
            for r in (json.loads(body).get("resultList") or {}).get("result", []):
                cores[(str(r.get("source")), str(r.get("id")))] = r
        missing = [k for k in keys if k not in cores]
        for k in keys:
            if k in cores:
                rec = parse_epmc_record(cores[k])
            else:  # no core record: the citation's own fields, without an abstract
                c = next(c for c in listed if (str(c.get("source")), str(c.get("id"))) == k)
                rec = dict(pmid=str(c.get("id")) if c.get("source") == "MED" else "", doi="", pmcid="",
                           title=c.get("title") or "", abstract="", journal=c.get("journalAbbreviation") or "",
                           year=int(c["pubYear"]) if str(c.get("pubYear", "")).isdigit() else None, pub_types=[],
                           languages=[], source=k[0], epmc_id=f"{k[0]}:{k[1]}", first_date="", cited_by=None,
                           no_core_record=True)
            rec["cites"] = p["key"]
            records.append(rec)
        summary[p["key"]] = dict(hit_count=hit, listed=len(listed), unique=len(keys), core=len(keys) - len(missing),
                                 without_core=len(missing), verified_cited_by=v.get("cited_by"))
        out.log(f"Citations of {p['key']} ({v['source']}:{v['id']}): {hit} counted, {len(keys)} unique, "
                f"{len(keys) - len(missing)} with a core record, {len(missing)} without (cited by "
                f"{v.get('cited_by')} at the lookup)")
    return records, summary


def merge(sources: dict) -> list:
    """One record per paper: the same PMID, else the same DOI, else the same normalized title (of at least 20
    characters). The fields are the PubMed record's where there is one, else Europe PMC's MED record's,
    else the first; every id and every source is kept. Numbered by PMID, then DOI, then title."""
    items = []
    for name, recs in sources.items():
        for r in recs:
            items.append((name, r))
    parent = list(range(len(items)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    seen = {}
    for i, (_, r) in enumerate(items):
        keys = [("pmid", r["pmid"])] if r.get("pmid") else []
        keys += [("doi", r["doi"])] if r.get("doi") else []
        t = norm_title(r.get("title"))
        keys += [("title", t)] if len(t) >= 20 else []
        for k in keys:
            if k in seen:
                parent[find(i)] = find(seen[k])
            else:
                seen[k] = i
    groups = {}
    for i in range(len(items)):
        groups.setdefault(find(i), []).append(items[i])
    merged = []
    rank = {"pubmed": 0, "europepmc": 1}
    for members in groups.values():
        members.sort(key=lambda m: (rank.get(m[0], 2), 0 if m[1].get("source") == "MED" else 1))
        best = dict(members[0][1])
        pmids = sorted({m[1]["pmid"] for m in members if m[1].get("pmid")}, key=lambda x: (len(x), x))
        dois = sorted({m[1]["doi"] for m in members if m[1].get("doi")})
        best.update(pmids=pmids, dois=dois, epmc_ids=sorted({m[1]["epmc_id"] for m in members if m[1].get("epmc_id")}),
                    found_by=sorted({m[0] if m[0] != "citations" else f"cites:{m[1]['cites']}" for m in members}),
                    preprint=all(m[1].get("source") == "PPR" for m in members))
        if not best.get("abstract"):
            best["abstract"] = next((m[1]["abstract"] for m in members if m[1].get("abstract")), "")
        for k in ("cites", "cited_by", "no_core_record"):
            best.pop(k, None)
        merged.append(best)

    def order(r):
        return (int(r["pmids"][0]) if r["pmids"] and r["pmids"][0].isdigit() else float("inf"),
                r["dois"][0] if r["dois"] else "~", norm_title(r["title"]))
    merged.sort(key=order)
    for i, r in enumerate(merged, 1):
        r["rid"] = f"R{i:05d}"
    return merged


def draw_sample(n: int, size: int = SAMPLE_SIZE, seed: int = SAMPLE_SEED) -> list:
    """The records' numbers (0-based, sorted) the screening reads: all, or `size` of them drawn without
    replacement by numpy's default generator with the protocol's seed."""
    import numpy as np
    if n <= size:
        return list(range(n))
    return sorted(int(i) for i in np.random.default_rng(seed).choice(n, size, replace=False))


def search(net, out: Out, papers: dict, window_limit: int = WINDOW_LIMIT) -> dict:
    out.log("## The queries")
    out.log(f"PubMed: {PUBMED_TERM}")
    out.log(f"Europe PMC: {EUROPEPMC_QUERY}")
    out.log()
    out.log("## The sources")
    pm, pm_s = pubmed(net, out, window_limit)
    ep, ep_s = europepmc(net, out)
    ci, ci_s = citations(net, out, papers)
    merged = merge(dict(pubmed=pm, europepmc=ep, citations=ci))
    import numpy as np
    sample = draw_sample(len(merged))
    for i in sample:
        merged[i]["in_sample"] = True
    content = "".join(json.dumps(r, ensure_ascii=False, sort_keys=True) + "\n" for r in merged).encode()
    with gzip.GzipFile(filename="", mode="wb", fileobj=(out.root / "records.jsonl.gz").open("wb"), mtime=0) as g:
        g.write(content)  # no name or time in the header: the same records give the same file
    (out.root / "sample.txt").write_text("".join(f"{merged[i]['rid']}\n" for i in sample))
    found = {}
    for r in merged:
        key = "+".join(sorted({"cites" if s.startswith("cites:") else s for s in r["found_by"]}))
        found[key] = found.get(key, 0) + 1
    out.log()
    out.log("## The merged list")
    out.log(f"{len(pm) + len(ep) + len(ci)} records from the sources -> {len(merged)} unique records "
            f"({sum(r['preprint'] for r in merged)} preprints only); found by: "
            + ", ".join(f"{k} {v}" for k, v in sorted(found.items())))
    out.log(f"The sample for screening: {len(sample)} of {len(merged)}"
            + (f" (numpy {np.__version__}, default_rng({SAMPLE_SEED}).choice({len(merged)}, {SAMPLE_SIZE}, "
               "replace=False))" if len(merged) > SAMPLE_SIZE else " (all the records: no more than 5,000)"))
    return dict(pubmed=pm_s, europepmc=ep_s, citations=ci_s, merged=len(merged), found_by=found,
                sample=len(sample), sample_seed=SAMPLE_SEED if len(merged) > SAMPLE_SIZE else None,
                numpy=np.__version__, records_content_sha256=hashlib.sha256(content).hexdigest())


# --------------------------------------------------------------------------- #
# the full texts (after the screening of titles and abstracts)
# --------------------------------------------------------------------------- #
def read_fulltext_list(path: Path) -> list:
    lines = [line.split("\t") for line in path.read_text().splitlines() if line.strip()]
    if not lines or lines[0] != FULLTEXT_COLUMNS:
        raise SystemExit(f"{path}: the header must be {FULLTEXT_COLUMNS}")
    items = [dict(zip(FULLTEXT_COLUMNS, row + [""] * (len(FULLTEXT_COLUMNS) - len(row)))) for row in lines[1:]]
    for it in items:
        it["epmc_ids"] = [e for e in it["epmc_ids"].split(",") if e]
    return items


def _pmcid_of(core: dict) -> str:
    ids = [str(i) for i in ((core.get("fullTextIdList") or {}).get("fullTextId") or []) if str(i).startswith("PMC")]
    return str(core.get("pmcid") or (ids[0] if ids else ""))


def fulltext(net, out: Out, items: list) -> dict:
    """Each listed record's open-access full text from Europe PMC, as it came; one row per record of what was
    tried and found (`fulltext_index.tsv`)."""
    out.log(f"## The full texts of {len(items)} records (Europe PMC, open access, JATS XML)")
    keys = []  # the core records of those without a PMCID: by PMID, else by their Europe PMC ids
    for it in items:
        if not it["pmcid"]:
            keys += [("MED", it["pmid"])] if it["pmid"] else [tuple(e.split(":", 1)) for e in it["epmc_ids"]]
    keys = list(dict.fromkeys(keys))
    cores = {}
    for i in range(0, len(keys), CITING_BATCH):
        batch = keys[i:i + CITING_BATCH]
        q = " OR ".join(f"(EXT_ID:{cid} AND SRC:{src})" for src, cid in batch)
        body = net.get(f"{EPMC}/search", dict(query=q, resultType="core", format="json", pageSize=100))
        out.raw(f"fulltext/core-{i // CITING_BATCH:04d}.json", body)
        for r in (json.loads(body).get("resultList") or {}).get("result", []):
            cores[(str(r.get("source")), str(r.get("id")))] = r
    out.log(f"Core records asked for {len(keys)} records without a PMCID: {sum(k in cores for k in keys)} found")
    rows, counts = [], dict(fetched=0, none=0)
    for it in items:
        tried, got, note = [], None, ""
        if it["pmcid"]:
            tried.append(it["pmcid"])
        else:
            mine = [("MED", it["pmid"])] if it["pmid"] else [tuple(e.split(":", 1)) for e in it["epmc_ids"]]
            for k in mine:
                c = cores.get(k)
                if c is not None and _pmcid_of(c):
                    tried.append(_pmcid_of(c))
            tried += [e.split(":", 1)[1] for e in it["epmc_ids"] if e.startswith("PPR:")]  # "PPR:PPR123" -> "PPR123"
            if not any(k in cores for k in mine):
                note = "no core record"
        tried = list(dict.fromkeys(tried))
        for fid in tried:
            try:
                body = net.get(f"{EPMC}/{fid}/fullTextXML", {})
            except urllib.error.HTTPError as e:
                if e.code in (404, 400):
                    note = (note + "; " if note else "") + f"{fid}: HTTP {e.code}"
                    continue
                raise
            if body.strip():
                got = fid
                path = out.raw(f"fulltext/{it['rid']}.xml", body)
                break
        counts["fetched" if got else "none"] += 1
        rows.append(dict(rid=it["rid"], tried=",".join(tried) or "-", fulltext=got or "-",
                         file=path if got else "-", note=note or "-"))
    with (out.root / "fulltext_index.tsv").open("w") as f:
        f.write("rid\ttried\tfulltext\tfile\tnote\n")
        for r in rows:
            f.write("\t".join(r[c] for c in ("rid", "tried", "fulltext", "file", "note")) + "\n")
    out.log(f"Full texts: {counts['fetched']} of {len(items)} fetched; {counts['none']} without an accessible full text")
    return dict(n=len(items), **counts, core_asked=len(keys), core_found=sum(k in cores for k in keys))


# --------------------------------------------------------------------------- #
# the command line
# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("mode", choices=("lookup", "search", "fulltext", "remerge"))
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--papers", type=Path, default=PAPERS)
    ap.add_argument("--fabricated", action="store_true", help="fabricated responses: the dry run, nothing fetched")
    ap.add_argument("--window-limit", type=int, default=WINDOW_LIMIT, help=argparse.SUPPRESS)
    ap.add_argument("--list", type=Path, default=FULLTEXT_LIST, help="fulltext: the records whose full texts are fetched")
    ap.add_argument("--search", type=Path, help="remerge: the finished search's directory (raw/, requests.tsv)")
    a = ap.parse_args(argv)
    if a.mode == "remerge":
        if a.search is None:
            ap.error("remerge needs --search")
        out = Out(a.out, write_raw=False)
        _header(out, a.mode, False)
        out.log(f"The search's raw responses: {a.search} (each request answered by the response it recorded)")
        out.log()
        own = a.search / "method_papers.json"  # the results branch keeps the papers the run used
        papers = json.loads((own if own.exists() else a.papers).read_text())
        net = Replay(out.log, a.search)
        res = search(net, out, papers, a.window_limit)
        if net.i != len(net.recorded):
            raise RuntimeError(f"{len(net.recorded) - net.i} recorded requests were not replayed")
        out.log(f"\nAll {net.i} recorded requests replayed, each answered by its recorded response.")
        out.log(f"records.jsonl.gz holds {res['merged']} records; SHA256 of its content (unzipped) "
                f"{res['records_content_sha256']}")
        out.finish("remerge", net, dict(res, ok=True, replayed=net.i, papers=papers["papers"]))
        names = ["records.jsonl.gz", "sample.txt", "requests.tsv", "remerge.log", "remerge.json"]
        (out.root / "SHA256SUMS").write_text("".join(f"{hashlib.sha256((out.root / n).read_bytes()).hexdigest()}  {n}\n"
                                                     for n in names))
        return 0
    out = Out(a.out)
    _header(out, a.mode, a.fabricated)
    papers = json.loads(a.papers.read_text())
    if a.fabricated:
        from fabricated import Fabricated
        net = Fabricated(out.log, papers)
    else:
        net = Net(out.log)
    if a.mode == "fulltext":
        items = read_fulltext_list(a.list)
        res = fulltext(net, out, items)
        out.finish("fulltext", net, dict(res, ok=True, list=str(a.list)))
        return 0
    if a.mode == "lookup":
        res = lookup(net, out, papers)
        verified = dict(papers, papers=[dict(p, verified=r["verified"]) for p, r in zip(papers["papers"], res["papers"])])
        (a.out / "method_papers.verified.json").write_text(json.dumps(verified, indent=2, ensure_ascii=False) + "\n")
        out.finish("lookup", net, res)
        return 0 if res["ok"] else 1
    if not all(p.get("verified") for p in papers["papers"]):
        out.log("The method papers are not all verified in " + str(a.papers) + ": the search does not run.")
        out.finish("search", net, dict(ok=False, error="method papers not verified"))
        return 1
    res = search(net, out, papers, a.window_limit)
    out.finish("search", net, dict(res, ok=True, papers=papers["papers"]))
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(HERE))
    raise SystemExit(main())
