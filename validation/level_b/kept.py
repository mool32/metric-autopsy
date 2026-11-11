#!/usr/bin/env python3
"""Level B, the end of S3 (`plan_s3.md`, section 4, "The claims of a kept paper"): the kept papers and their claims,
for S4's search of verdicts, S5's data checks and S7's list (`amendment1/kept_claims.json`).

A paper is kept where both passes keep it at its full text: amendment 1's passes (S3) for the records they read,
the passes of 2026-10-19 for the four papers kept then and not read again. A kept paper's claims are its claim
lines in both passes, matched across the passes by their genes:

* **a pair** by its two symbols, order and case ignored: the field split at a comma, a semicolon, a slash, "&",
  "+", the word "and", an en or em dash, or a hyphen with a space on each side; a field with no such separator
  is split at a space, or at one of its hyphens (each hyphen a candidate: "SRGN-RAB31"), and two lines match where
  their candidates share a pair;
* **a set** (`SET: <name>, <n> genes, <where listed>`) by its name, case and spacing ignored;
* a field read neither way is matched only by the same text (letters and digits, case ignored).

What parentheses hold is dropped first, from a gene ("HSPA5 (BIP)") and from a set's name ("Mitochondria NDCG
(the electron transport chain)"); a set whose name holds an MSigDB-style identifier (capitals and underscores, as
HALLMARK_INTERFERON_GAMMA_RESPONSE) is matched by it (J21: found on the merge's output, before any use of it).

The claims both passes list come first, in pass 1's order, then those only pass 1 lists, then those only pass 2
lists; each carries the passes that list it and both passes' lines. A pair whose two symbols cannot be told (its
lines give several candidates and no other line settles them) is marked `ambiguous`, its genes empty. A set's
genes are filled later from where the paper lists them (not here). More than 30 kept papers: the protocol's draw
(section 9) names the 30 listed.

    python -I validation/level_b/kept.py --ft1 validation/level_b/amendment1/ft/pass1 \
        --ft2 validation/level_b/amendment1/ft/pass2 --old-ft1 validation/level_b/fulltext_screening/pass1 \
        --old-ft2 validation/level_b/fulltext_screening/pass2 --old-kept R00189 R01225 R01962 R02693 \
        --ids validation/level_b/fulltext_list.tsv validation/level_b/amendment1/fulltext_list.tsv \
        --titles validation/level_b/screening/all.tsv validation/level_b/amendment1/ta/all.tsv \
        --out validation/level_b/amendment1/kept_claims.json
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from screening import PAPERS_MAX, PAPERS_MIN, PAPERS_SEED, read_fulltext_pass  # noqa: E402

_STRONG = re.compile(r"\s*(?:[,;/&+]|\band\b|\s[-–—]\s|[–—])\s*", flags=re.IGNORECASE)


def _clean(field: str) -> str:
    return field.strip().strip("*`").strip()


def _no_parens(text: str) -> str:
    """The text without what parentheses hold (an alias, a description): 'HSPA5 (BIP)' -> 'HSPA5' (J21)."""
    prev = None
    while prev != text:
        prev, text = text, re.sub(r"\s*\([^()]*\)", "", text)
    return " ".join(text.split())


_MSIGDB = re.compile(r"\b[A-Z][A-Z0-9]*(?:_[A-Z0-9]+)+\b")


def is_set(field: str) -> bool:
    return _clean(field).upper().startswith("SET:")


def set_parts(field: str) -> dict:
    """A set's line: its name (before the first comma), its number of genes and where they are listed."""
    body = _clean(field)[4:].strip()
    depth, cut = 0, len(body)
    for i, ch in enumerate(body):  # the first comma outside parentheses ends the name
        depth += (ch == "(") - (ch == ")")
        if ch == "," and depth == 0:
            cut = i
            break
    name, rest = body[:cut], body[cut + 1:]
    m = re.search(r"(\d[\d,]*)\s+genes?", _no_parens(rest))
    where = _no_parens(rest)
    return dict(set_name=name.strip(), set_size=int(m.group(1).replace(",", "")) if m else None,
                set_where=(where[m.end():] if m else where).strip(" ,;"))


def pair_candidates(field: str) -> set:
    """The unordered pairs of symbols a pair's `genes` field can name (upper case)."""
    toks = [t.strip().upper() for t in _STRONG.split(_no_parens(_clean(field))) if t and t.strip()]
    if len(toks) == 2:
        return {frozenset(toks)} if toks[0] != toks[1] else set()
    if len(toks) != 1:
        return set()
    t, out = toks[0], set()
    ws = t.split()
    if len(ws) == 2 and ws[0] != ws[1]:
        out.add(frozenset(ws))
    if len(ws) == 1:
        for m in re.finditer("-", t):
            a, b = t[:m.start()], t[m.end():]
            if a and b and a != b:
                out.add(frozenset((a, b)))
    return out


def keys(field: str) -> set:
    if is_set(field):
        name = set_parts(field)["set_name"]
        ids = _MSIGDB.findall(_no_parens(name)) or _MSIGDB.findall(name)
        return {("set", ids[0].lower() if ids else _no_parens(name).lower())}
    cands = pair_candidates(field)
    if cands:
        return {("pair", c) for c in cands}
    return {("raw", re.sub(r"[^A-Za-z0-9]", "", field).upper())}


def _ordered(field: str, pair: frozenset) -> list:
    """The pair's symbols in the order the line names them."""
    up = _no_parens(_clean(field)).upper()
    return sorted(pair, key=lambda g: up.find(g) if up.find(g) >= 0 else len(up))


def paper_claims(lines1: list, lines2: list) -> list:
    """Both passes' claim lines of one paper, matched by their genes (module docstring)."""
    k1 = [keys(c["genes"]) for c in lines1]
    k2 = [keys(c["genes"]) for c in lines2]
    used2, matched = set(), []
    for i in range(len(lines1)):
        j = next((j for j in range(len(lines2)) if j not in used2 and k1[i] & k2[j]), None)
        if j is not None:
            used2.add(j)
        matched.append((i, j))
    groups = ([(i, j) for i, j in matched if j is not None] + [(i, None) for i, j in matched if j is None]
              + [(None, j) for j in range(len(lines2)) if j not in used2])
    claims = []
    for n, (i, j) in enumerate(groups, start=1):
        lines = {}
        if i is not None:
            lines["1"] = lines1[i]
        if j is not None:
            lines["2"] = lines2[j]
        common = (k1[i] if i is not None else k2[j]) & (k2[j] if j is not None else k1[i])
        first = lines.get("1", lines.get("2"))
        c = dict(id=str(n), passes=[int(p) for p in lines], lines=lines)
        if any(k[0] == "set" for k in common):
            c.update(kind="set", genes=[], **set_parts(first["genes"]))
        else:
            pairs = [k[1] for k in common if k[0] == "pair"]
            c.update(kind="pair", genes=_ordered(first["genes"], pairs[0]) if len(pairs) == 1 else [],
                     ambiguous=len(pairs) != 1)
        claims.append(c)
    return claims


def _table(paths: list, key_col: int, cols: dict) -> dict:
    out = {}
    for p in paths:
        lines = Path(p).read_text().splitlines()
        for line in lines[1:]:
            r = line.split("\t")
            if r and r[0] and r[0] not in out:
                out[r[0]] = {k: (r[i] if i < len(r) else "") for k, i in cols.items()}
    return out


def kept_claims(ft1: Path, ft2: Path, old1: Path, old2: Path, old_kept: list, ids: list, titles: list) -> dict:
    import numpy as np
    d1, c1, p1 = read_fulltext_pass(ft1)
    d2, c2, p2 = read_fulltext_pass(ft2)
    o1, oc1, q1 = read_fulltext_pass(old1)
    o2, oc2, q2 = read_fulltext_pass(old2)
    if p1 or p2:
        raise SystemExit("amendment 1's full-text files do not check:\n" + "\n".join(p1 + p2))
    new = sorted(r for r in set(d1) & set(d2) if d1[r]["decision"] == d2[r]["decision"] == "KEEP")
    old = sorted(r for r in old_kept if r not in d1 and r not in d2
                 and o1.get(r, {}).get("decision") == o2.get(r, {}).get("decision") == "KEEP")
    missing = sorted(set(old_kept) - set(old) - set(d1) - set(d2))
    if missing:
        raise SystemExit(f"not kept by both passes of 2026-10-19: {missing}")
    id_of = _table(ids, 0, dict(pmid=1, pmcid=2, doi=3, epmc_ids=4))
    title_of = _table(titles, 0, dict(title=11))
    kept = sorted(set(new) | set(old))
    papers = []
    for rid in kept:
        src = "2026-10-19" if rid in old else "amendment 1"
        a, b = (oc1, oc2) if rid in old else (c1, c2)
        i = id_of.get(rid, {})
        papers.append(dict(rid=rid, read=src, title=title_of.get(rid, {}).get("title", ""), pmid=i.get("pmid", ""),
                           pmcid=i.get("pmcid", ""), doi=i.get("doi", ""),
                           epmc_ids=[e for e in i.get("epmc_ids", "").split(",") if e],
                           claims=paper_claims(a.get(rid, []), b.get(rid, []))))
    listed = ([kept[k] for k in sorted(np.random.default_rng(PAPERS_SEED).choice(len(kept), PAPERS_MAX, replace=False))]
              if len(kept) > PAPERS_MAX else list(kept))
    cl = [c for p in papers for c in p["claims"]]
    with_pairs = sum(any(c["kind"] == "pair" for c in p["claims"]) for p in papers)
    counts = dict(papers=len(papers), read_in_amendment_1=len(new), kept_on_2026_10_19=len(old),
                  papers_with_pairs=with_pairs, papers_with_sets=sum(any(c["kind"] == "set" for c in p["claims"]) for p in papers),
                  claims=len(cl), pairs=sum(c["kind"] == "pair" for c in cl), sets=sum(c["kind"] == "set" for c in cl),
                  both_passes=sum(len(c["passes"]) == 2 for c in cl), one_pass=sum(len(c["passes"]) == 1 for c in cl),
                  ambiguous=sum(bool(c.get("ambiguous")) for c in cl))
    return dict(counts=counts, drawn=len(kept) > PAPERS_MAX, fewer_than_target=with_pairs < PAPERS_MIN,
                draw=f"numpy.random.default_rng({PAPERS_SEED}).choice({len(kept)}, {PAPERS_MAX}, replace=False)"
                if len(kept) > PAPERS_MAX else None, numpy=np.__version__, listed=listed, papers=papers)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ft1", type=Path, required=True)
    ap.add_argument("--ft2", type=Path, required=True)
    ap.add_argument("--old-ft1", type=Path, required=True)
    ap.add_argument("--old-ft2", type=Path, required=True)
    ap.add_argument("--old-kept", nargs="+", required=True)
    ap.add_argument("--ids", nargs="+", required=True)
    ap.add_argument("--titles", nargs="+", required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    res = kept_claims(a.ft1, a.ft2, a.old_ft1, a.old_ft2, a.old_kept, a.ids, a.titles)
    a.out.write_text(json.dumps(res, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps(dict(counts=res["counts"], drawn=res["drawn"], fewer_than_target=res["fewer_than_target"]), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
