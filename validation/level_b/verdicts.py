#!/usr/bin/env python3
"""Level B, step S4 (`plan_s4.md`): the search for each kept claim's independent verdict (the protocol's section 7).

Run once in GitHub Actions (`.github/workflows/level-b-search.yml`, the label `level-b-verdicts`); this
environment does not reach Europe PMC. For every kept paper of `amendment1/kept_claims.json` (written from S3's
list) and each of its claims, Europe PMC, every source, no window, every response kept as it came:

* **the forward citations of the claim's paper** (its citations endpoint, by its PMID, else its PMCID, else its
  preprint's id; every page; then each citing record's core record, as `search.py`'s citations);
* **the claim's genes with the terms**: a pair's two genes each in the title or the abstract, and one of section 3's
  terms of coordination or of perturbation (`PERTURBATION`) in the title or the abstract; a set by its name as a
  phrase where the paper names it (else by its forward citations only);
* **the filter fixed before any reading** (`passes`): a record is read for a claim where its title or abstract names
  both genes of the pair (each as a word, case ignored), or the set's name or two of its genes; the claim's own
  paper is left out;
* **the full texts** of the records the filter keeps (Europe PMC's open access, by PMCID or preprint id).

    python -I validation/level_b/verdicts.py --claims validation/level_b/amendment1/kept_claims.json --out OUT
    python -I validation/level_b/verdicts.py --claims ... --out OUT --fabricated     # the dry run, nothing fetched

OUT gets verdicts.log, verdicts.json (the counts per claim), found.jsonl.gz (per claim the records found, by source,
and whether the filter keeps each), fulltext_index.tsv, raw/ (every response) and requests.tsv.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import urllib.error
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import search as S  # noqa: E402  the search's network, raw store, parser and citations, unchanged

PERTURBATION = ["knockdown", "knockout", "overexpression", "smFISH", "in situ"]
TERMS = S._epmc_block(S.COORDINATION + PERTURBATION)


def claim_query(claim: dict) -> str | None:
    """The claim's Europe PMC query, or None for a set without a name (its forward citations only)."""
    if claim["kind"] == "set":
        name = (claim.get("set_name") or "").strip()
        return f'(TITLE:"{name}" OR ABSTRACT:"{name}") AND {TERMS}' if name else None
    a, b = claim["genes"][:2]
    return f'(TITLE:"{a}" OR ABSTRACT:"{a}") AND (TITLE:"{b}" OR ABSTRACT:"{b}") AND {TERMS}'


def _word(text: str, word: str) -> bool:
    return re.search(rf"(?<![A-Za-z0-9]){re.escape(word)}(?![A-Za-z0-9])", text, flags=re.IGNORECASE) is not None


def passes(rec: dict, claim: dict) -> bool:
    """The filter fixed before any reading (`plan_s4.md`, section 2)."""
    text = f"{rec.get('title') or ''} {rec.get('abstract') or ''}"
    if claim["kind"] == "set":
        name = (claim.get("set_name") or "").strip()
        return bool(name and _word(text, name)) or sum(_word(text, g) for g in claim["genes"]) >= 2
    return all(_word(text, g) for g in claim["genes"][:2])


def _own(rec: dict, paper: dict) -> bool:
    ids = {("pmid", paper.get("pmid") or ""), ("doi", (paper.get("doi") or "").lower()), ("pmcid", paper.get("pmcid") or "")}
    return any(v and (k, v) in ids for k, v in (("pmid", rec.get("pmid")), ("doi", rec.get("doi")), ("pmcid", rec.get("pmcid"))))


def paper_key(paper: dict) -> dict | None:
    """The citations endpoint's source and id for a paper: MED and its PMID, else PMC and its PMCID, else its
    preprint's PPR id."""
    if paper.get("pmid"):
        return dict(source="MED", id=paper["pmid"])
    if paper.get("pmcid"):
        return dict(source="PMC", id=paper["pmcid"])
    ppr = [e.split(":", 1)[1] for e in paper.get("epmc_ids") or [] if e.startswith("PPR:")]
    return dict(source="PPR", id=ppr[0]) if ppr else None


def search_records(net, out, query: str, tag: str) -> list:
    cursor, page, records = "*", 0, []
    while True:
        body = net.get(f"{S.EPMC}/search", dict(query=query, resultType="core", format="json", pageSize=1000,
                                                cursorMark=cursor))
        out.raw(f"verdicts/search-{tag}-{page:04d}.json", body)
        d = json.loads(body)
        res = (d.get("resultList") or {}).get("result", [])
        records += [S.parse_epmc_record(r) for r in res]
        nxt, page = d.get("nextCursorMark"), page + 1
        if not res or not nxt or nxt == cursor:
            return records
        cursor = nxt


def fulltexts(net, out, recs: list) -> list:
    rows = []
    for r in recs:
        ids = [r["pmcid"]] if r.get("pmcid") else []
        ids += [r["epmc_id"].split(":", 1)[1]] if r.get("epmc_id", "").startswith("PPR:") else []
        got, note = None, ""
        for fid in ids:
            try:
                body = net.get(f"{S.EPMC}/{fid}/fullTextXML", {})
            except urllib.error.HTTPError as e:
                if e.code not in (400, 404):
                    raise
                note += f"{fid}: HTTP {e.code}; "
                continue
            except S.Unavailable as e:
                if e.status is None:
                    raise
                note += f"{fid}: HTTP {e.status} on all 6 attempts; "
                continue
            if body.strip():
                got = (fid, out.raw(f"verdicts/fulltext-{fid}.xml", body))
                break
        rows.append(dict(epmc_id=r["epmc_id"], tried=",".join(ids) or "-", fulltext=got[0] if got else "-",
                         file=got[1] if got else "-", note=note.strip("; ") or "-"))
    return rows


def run(net, out, kept: dict) -> dict:
    found, summary, to_read = [], [], {}
    for paper in kept["papers"]:
        key = paper_key(paper)
        citing = []
        if key:
            recs, _ = S.citations(net, out, dict(papers=[dict(key=paper["rid"], verified=key)]))
            citing = recs
        for claim in paper["claims"]:
            q = claim_query(claim)
            by_terms = search_records(net, out, q, f"{paper['rid']}-{claim['id']}") if q else []
            seen, rows = set(), []
            for source, recs in (("citations", citing), ("terms", by_terms)):
                for r in recs:
                    if _own(r, paper):
                        continue
                    k = r["epmc_id"]
                    if k in seen:
                        next(x for x in rows if x["epmc_id"] == k)["found_by"].append(source)
                        continue
                    seen.add(k)
                    keep = passes(r, claim)
                    rows.append(dict(epmc_id=k, found_by=[source], kept=keep, pmid=r.get("pmid"), pmcid=r.get("pmcid"),
                                     doi=r.get("doi"), title=r.get("title"), abstract=r.get("abstract"),
                                     year=r.get("year"), journal=r.get("journal")))
                    if keep:
                        to_read[k] = r
            found.append(dict(rid=paper["rid"], claim=claim["id"], kind=claim["kind"], query=q, records=rows))
            summary.append(dict(rid=paper["rid"], claim=claim["id"], citing=len(citing), by_terms=len(by_terms),
                                records=len(rows), kept=sum(x["kept"] for x in rows)))
            out.log(f"{paper['rid']} {claim['id']} ({claim['kind']}): citing {len(citing)}, by the terms {len(by_terms)}, "
                    f"records {len(rows)}, kept by the filter {summary[-1]['kept']}")
    ft = fulltexts(net, out, list(to_read.values()))
    with gzip.GzipFile(filename="", mode="wb", fileobj=(out.root / "found.jsonl.gz").open("wb"), mtime=0) as g:
        for f in found:
            g.write((json.dumps(f, ensure_ascii=False, sort_keys=True) + "\n").encode())
    with (out.root / "fulltext_index.tsv").open("w") as f:
        f.write("epmc_id\ttried\tfulltext\tfile\tnote\n")
        for r in ft:
            f.write("\t".join(r[c] for c in ("epmc_id", "tried", "fulltext", "file", "note")) + "\n")
    out.log(f"Records the filter keeps, over the claims: {len(to_read)}; their full texts fetched: "
            f"{sum(r['fulltext'] != '-' for r in ft)}")
    return dict(claims=summary, to_read=len(to_read), fulltexts=sum(r["fulltext"] != "-" for r in ft))


class Fabricated:
    """Europe PMC's answers for the dry run (FABRICATED): two citing records per paper, one of which names a pair's
    genes; one record per claim's query; a full text for every second PMCID."""

    def __init__(self, log):
        self.log, self.rows = log, []

    def _record(self, i: int, title: str, abstract: str) -> dict:
        return dict(id=str(900000 + i), source="MED", pmid=str(900000 + i), pmcid=f"PMC{800000 + i}" if i % 2 == 0 else "",
                    doi=f"10.9999/fab.{i}", title=title, abstractText=abstract, pubYear="2024",
                    journalInfo=dict(journal=dict(isoabbreviation="Fab J")))

    def get(self, url: str, params: dict, post: bool = False) -> bytes:
        body = self._answer(url, params)
        self.rows.append(dict(time="FABRICATED", method="GET", url=url, params=S.shown_params(params), status=200,
                              bytes=len(body), sha256="-", attempt=1))
        return body

    def _answer(self, url: str, params: dict) -> bytes:
        if url.endswith("/citations"):
            pid = url.rsplit("/", 2)[-2]
            return json.dumps(dict(hitCount=2, citationList=dict(citation=[
                dict(id=str(900100 + int(pid[-2:]) * 2 + k), source="MED", title=f"citing {k}") for k in (0, 1)]))).encode()
        if url.endswith("/fullTextXML"):
            return b"<article><body><sec><title>Results</title><p>FABRICATED</p></sec></body></article>"
        q = params.get("query", "")
        if "EXT_ID:" in q:
            ids = re.findall(r"EXT_ID:(\d+)", q)
            return json.dumps(dict(resultList=dict(result=[
                self._record(int(i) - 900000, "A citing paper about ACE2 and TMPRSS2 coexpression" if int(i) % 2 == 0
                             else "A citing paper about something else", "FABRICATED") for i in ids]))).encode()
        if params.get("cursorMark") != "*":
            return json.dumps(dict(resultList=dict(result=[]), nextCursorMark=params.get("cursorMark"))).encode()
        genes = re.findall(r'TITLE:"([^"]+)" OR ABSTRACT:"\1"\) AND', q)
        title = "A knockdown study of " + " and ".join(genes[:2]) if genes else "nothing"
        return json.dumps(dict(hitCount=1, nextCursorMark="AoE", resultList=dict(result=[
            self._record(500 + len(q) % 97, title, "FABRICATED")]))).encode()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--claims", type=Path, required=True, help="amendment1/kept_claims.json")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--fabricated", action="store_true", help="fabricated answers: the dry run, nothing fetched")
    a = ap.parse_args(argv)
    out = S.Out(a.out)
    out.log(f"# Level B, S4: the independent verdicts' search{' — FABRICATED' if a.fabricated else ''} (verdicts.py)")
    out.log("")
    kept = json.loads(a.claims.read_text())
    net = Fabricated(out.log) if a.fabricated else S.Net(out.log)
    res = run(net, out, kept)
    out.finish("verdicts", net, dict(res, ok=True, fabricated=a.fabricated, claims_file=str(a.claims)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
