#!/usr/bin/env python3
"""The first list of level B's search against the corrected one (`search.py remerge`): what the PubMed
parser's bug changed. Nothing is fetched; both lists and the raw responses are read as they are.

1. PubMed's records parsed as the first list parsed them (every `ArticleId` under `PubmedData`, the
   references' included) and as now (`PubmedData/ArticleIdList` only): how many got another DOI or PMCID.
2. The first list's records that held several PMIDs, and how many of them the corrected list splits.
3. The corrected list's records with several PMIDs: whether their PubMed records share a title and one is
   a preprint.
4. The records' texts as the screeners read them (`screening.as_text` without the number): how many of the
   first list's are in the corrected list unchanged, and how many of the corrected list's are new.

    python -I validation/level_b/remerge_compare.py <results/level-b-search> <remerge dir> > validation/level_b/runs/remerge_compare.log
"""
from __future__ import annotations

import gzip
import json
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import screening  # noqa: E402
from search import norm_title, parse_pubmed_xml  # noqa: E402

PREPRINT_DOI = ("10.1101/", "10.21203/", "10.26434/", "10.20944/", "10.2139/", "10.31219/", "10.36227/")


def _first_list_ids(art) -> dict:
    """The ids as the first list's parser read them: every ArticleId under PubmedData, the last one winning."""
    pd = art.find("PubmedData")
    return {i.get("IdType"): (i.text or "").strip() for i in (pd.iter("ArticleId") if pd is not None else [])}


def _records(path: Path) -> list:
    with gzip.open(path, "rt") as f:
        return [json.loads(line) for line in f]


def main(results: Path, remerge: Path) -> int:
    print(f"# Level B: the first list against the corrected one — {results.name}, {remerge.name}\n")
    # 1. the parser, both ways
    n = doi_diff = pmc_diff = pmc_from_ref = 0
    pubmed = {}
    for f in sorted((results / "raw" / "pubmed").glob("efetch-*.xml.gz")):
        body = gzip.open(f).read()
        now = {r["pmid"]: r for r in parse_pubmed_xml(body)}
        pubmed.update(now)
        for art in ET.fromstring(body).iter("PubmedArticle"):
            pmid = (art.findtext("MedlineCitation/PMID") or "").strip()
            if pmid not in now:
                continue
            n += 1
            old = _first_list_ids(art)
            a = art.find("MedlineCitation/Article")
            eloc = next((e.text for e in a.iter("ELocationID") if e.get("EIdType") == "doi"), "") if a is not None else ""
            old_doi = (old.get("doi") or eloc or "").strip().lower()
            doi_diff += old_doi != now[pmid]["doi"]
            pmc_diff += old.get("pmc", "") != now[pmid]["pmcid"]
            pmc_from_ref += bool(old.get("pmc")) and not now[pmid]["pmcid"]
    print("## 1. PubMed's records, parsed as the first list did and as now")
    print(f"  {n} records: another DOI in the first list {doi_diff}; another PMCID {pmc_diff} (of which a PMCID where "
          f"the record has none: {pmc_from_ref})\n")
    first, corrected = _records(results / "records.jsonl.gz"), _records(remerge / "records.jsonl.gz")
    # 2. the first list's records with several PMIDs
    where = {p: r["rid"] for r in corrected for p in r["pmids"]}
    multi = [r for r in first if len(r["pmids"]) > 1]
    split = [r for r in multi if len({where.get(p) for p in r["pmids"]}) > 1]
    print("## 2. The first list's records with several PMIDs")
    print(f"  {len(multi)} records; the corrected list splits {len(split)} of them, which held "
          f"{sum(len({where.get(p) for p in r['pmids']}) for r in split)} papers (the corrected list's records); PMIDs "
          f"per split record: " + ", ".join(f"{k}: {v}" for k, v in sorted(Counter(len(r['pmids']) for r in split).items())))
    print(f"  the others ({len(multi) - len(split)}) are records of the corrected list as they were\n")
    # 3. the corrected list's records with several PMIDs
    cm = [r for r in corrected if len(r["pmids"]) > 1]
    same_title = [r for r in cm if len({norm_title(pubmed[p]["title"]) for p in r["pmids"] if p in pubmed}) == 1]
    with_preprint = [r for r in cm if any("Preprint" in (pubmed[p].get("pub_types") or []) or
                                          pubmed[p]["doi"].startswith(PREPRINT_DOI) for p in r["pmids"] if p in pubmed)]
    print("## 3. The corrected list's records with several PMIDs")
    print(f"  {len(cm)} records ({', '.join(f'{k} PMIDs: {v}' for k, v in sorted(Counter(len(r['pmids']) for r in cm).items()))}); "
          f"their PubMed records share one title (normalized) in {len(same_title)}; one of them is a preprint (type "
          f"'Preprint' or a preprint server's DOI) in {len(with_preprint)}")
    cores = {}
    for f in sorted((results / "raw" / "europepmc").glob("*.json.gz")):
        for c in (json.loads(gzip.open(f).read()).get("resultList") or {}).get("result", []):
            cores[f"{c.get('source')}:{c.get('id')}"] = c
    for r in cm:
        if r not in same_title or r not in with_preprint:
            print(f"    {r['rid']} {r['pmids']} {r['dois']}:")
            for p in r["pmids"]:
                print(f"      PubMed {p}: DOI {pubmed[p]['doi']}; {pubmed[p]['title']}")
            for e in r["epmc_ids"]:
                c = cores.get(e, {})
                print(f"      Europe PMC {e}: DOI {(c.get('doi') or '').lower()}; {c.get('title')}")
    print()
    # 4. the texts the screeners read
    t_first = Counter(screening.screened_text(r) for r in first)
    t_corr = Counter(screening.screened_text(r) for r in corrected)
    print("## 4. The records' texts as the screeners read them")
    print(f"  the first list: {len(first)} records, {len(t_first)} distinct texts; the corrected list: {len(corrected)} "
          f"records, {len(t_corr)} distinct texts")
    print(f"  the first list's texts in the corrected list unchanged: {sum(1 for t in t_first if t in t_corr)}; "
          f"the corrected list's texts not in the first list: {sum(1 for t in t_corr if t not in t_first)}")
    by_text = {screening.screened_text(r): r for r in corrected}
    other_doi = sum(1 for r in first if r["doi"] != by_text[screening.screened_text(r)]["doi"])
    other_pmc = sum(1 for r in first if (r.get("pmcid") or "") != (by_text[screening.screened_text(r)].get("pmcid") or ""))
    print(f"  of the first list's records, {other_doi} have another DOI in the corrected list and {other_pmc} another PMCID")
    return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]), Path(sys.argv[2])))
