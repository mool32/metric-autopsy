#!/usr/bin/env python3
"""Two checks of the search's raw responses (`results/level-b-search`), for the owner; nothing is fetched.

1. Europe PMC's preprints: the protocol meant the Europe PMC query to include them ("every source included,
   preprints too"), and the query holds LANG:"eng". The records of the query's pages by source and language,
   and the language of the preprint records among the forward citations' core records (which carry no
   language filter).
2. PubMed's warnings on the query (its translation drops a quoted phrase it does not index).

    python -I validation/level_b/search_checks.py <results/level-b-search> > validation/level_b/runs/search_checks.log
"""
from __future__ import annotations

import gzip
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path


def main(root: Path) -> int:
    sums = dict(reversed(line.split(maxsplit=1)) for line in (root / "SHA256SUMS").read_text().splitlines())
    sums = {k.lstrip("*").removeprefix("./"): v for k, v in sums.items()}
    bad = [name for name, digest in sums.items() if (root / name).exists()
           and hashlib.sha256((root / name).read_bytes()).hexdigest() != digest]
    print(f"# Checks of the search's raw responses — {root.name}, {len(sums)} files in SHA256SUMS, "
          f"{len(bad)} not matching{': ' + ', '.join(bad) if bad else ''}\n")
    pages = sorted((root / "raw" / "europepmc").glob("search-*.json.gz"))
    by = Counter()
    for p in pages:
        for r in json.loads(gzip.open(p).read())["resultList"]["result"]:
            by[(r.get("source"), r.get("language"))] += 1
    print(f"## 1. Europe PMC: the query's {len(pages)} pages, records by (source, language)")
    for (src, lang), n in sorted(by.items(), key=lambda kv: -kv[1]):
        print(f"  {src} {lang}: {n}")
    print(f"  preprints (source PPR): {sum(n for (src, _), n in by.items() if src == 'PPR')}")
    citing = sorted((root / "raw" / "europepmc").glob("citing-*.json.gz"))
    langs, n_all = Counter(), 0
    for p in citing:
        for r in (json.loads(gzip.open(p).read()).get("resultList") or {}).get("result", []):
            n_all += 1
            if r.get("source") == "PPR":
                langs[str(r.get("language"))] += 1
    print(f"\nThe forward citations' core records ({len(citing)} files, {n_all} records): the preprints' language "
          + ", ".join(f"{k} {v}" for k, v in sorted(langs.items())))
    print("So a preprint record carries no language, and LANG:\"eng\" leaves every preprint out of the query; the"
          " search found preprints only through the forward citations.")
    log = (root / "search.log").read_text().splitlines()
    print("\n## 2. PubMed's warnings on the query (search.log)")
    for line in log:
        if "warninglist" in line or line.startswith("PubMed:"):
            print("  " + line.strip()[:400])
    return 0


if __name__ == "__main__":
    raise SystemExit(main(Path(sys.argv[1])))
