#!/usr/bin/env python3
"""Level B's screening of titles and abstracts (`validation/plans/level_b_protocol.md`, section 5): the
batches the screeners read, the check of what each pass wrote, and the comparison of the two passes.

    python -I validation/level_b/screening.py batches --search <results/level-b-search> --out <dir>
    python -I validation/level_b/screening.py check --batches <dir> --decisions <pass dir>
    python -I validation/level_b/screening.py merge --batches <dir> --pass1 <dir> --pass2 <dir> --out <dir>
    python -I validation/level_b/screening.py carry --old-batches <dir> --pass1 <dir> --pass2 <dir> --batches <dir> --out <dir>
    python -I validation/level_b/screening.py combine --batches <dir> --carried <dir> --pass1 <dir> --pass2 <dir> --out <dir>
    python -I validation/level_b/screening.py claims --merge <dir> --ft-merge <dir> --ft1 <dir> --ft2 <dir> --out <dir>
    python -I validation/level_b/screening.py fulltext-list --search <dir> --rids <fulltext.txt> --out <tsv>
    python -I validation/level_b/screening.py fulltext-text --fulltext <results/level-b-fulltext> --out <dir>
    python -I validation/level_b/screening.py fulltext-batches --merge <dir> --pass 1 --texts <dir> --out <dir>
    python -I validation/level_b/screening.py check-fulltext --batch <batch_NN.list> --decisions <tsv>
    python -I validation/level_b/screening.py merge-fulltext --merge <dir> --ft1 <dir> --ft2 <dir> --out <dir>

`batches` reads the search's merged list (`records.jsonl.gz`, checked against the search's SHA256SUMS) and its
sample (`sample.txt`, every record where the list holds at most 5,000), and writes the screened records in
the order of their numbers, `BATCH` to a batch, with the fields the screeners read (`screening.md`) and
nothing of how a record was found: `batch_NN.jsonl` for the tools, `batch_NN.txt` for the screener (the
same records as text in short lines). `check` checks a pass's decision files against `screening.md`'s
format; `check-one` one batch's file (the screener runs it before finishing).
`carry` takes a first list's decisions to a corrected list of the same search: a record of the corrected list
whose text, as the screener reads it (`as_text` without its number), a record of the first list had gets that
record's decisions from each pass, the records' numbers named in its notes translated; the rest are written
as batches that both passes screen as any batch (`rescreen/`). `combine` writes each pass's decisions on the
corrected list, the carried ones and the rescreen's, per batch of the corrected list, for `check` and `merge`.
`merge` compares the passes: Cohen's kappa over the three decisions, the counts of every step (a
PRISMA-style table) per pass, the records for the owner (decided differently by the two passes, or unclear
in either), the counts of claims about sets and global measures (section 1), and the records whose full text
is fetched (kept or unclear in either pass). `fulltext-list` writes those records' ids for `search.py fulltext`
(`fulltext_list.tsv`); `fulltext-text` turns each fetched full text (JATS XML) into the text a screener reads:
the title, the abstract, every section with its heading, the figures' and tables' captions and the tables'
cells, the back matter (data availability, notes); the references are left out. `fulltext-batches` gives
each pass the full texts of the records it kept (`FT_BATCH` to a batch; a record without a full text is
listed apart, for the owner), `check-fulltext` checks a full-text decision file against
`screening_fulltext.md`'s format, and `merge-fulltext` compares the passes' full-text decisions: Cohen's
kappa over the records both passes read, the counts per pass, the papers both keep with each pass's claims,
and what goes to the owner (decided differently, unclear, read by one pass only, or without a full text).
`claims` writes the list for the owner (`claims.md`, `claims.json`): both steps' counts per pass, the sets and
global measures, the papers both passes keep with each pass's claims (more than 30: 30 drawn by the protocol's
seed), the designs with fewer than `MIN_REPLICATES` donors in a group listed apart, the verdicts the papers
give, and the records for the owner.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
import textwrap
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

BATCH = 200
FIELDS = ("rid", "title", "abstract", "journal", "year", "pub_types", "pmid", "doi", "preprint")
WIDTH = 110
HEADER = ["rid", "decision", "code", "kinds", "note"]
DECISIONS = ("KEEP", "EXCLUDE", "UNCLEAR")
CODES = tuple(f"X{i}" for i in range(1, 11))
KINDS = ("pair", "set", "global")
FT_BATCH = 8
FT_HEADER = ["rid", "decision", "code", "claim", "genes", "species", "cells", "conditions", "direction", "metric",
             "where", "quote", "data", "donors", "rule4", "verdict", "note"]
FT_CODES = tuple(f"F{i}" for i in range(1, 11))
CODE_NAMES = {
    "X1": "not a paper with findings of its own", "X2": "no single-cell RNA counts",
    "X3": "development data or a control", "X4": "no changed coordination between conditions",
    "X5": "cell states within donors", "X6": "SET: a set's coordination", "X7": "GLOBAL: a global measure",
    "X8": "too few donors (the abstract shows it)", "X9": "data not public or not counts",
    "X10": "duplicate",
}
FT_CODE_NAMES = {
    "F1": "no claim about a named pair with a direction", "F2": "SET or GLOBAL only", "F3": "a metric v0.4 cannot check",
    "F4": "not single-cell RNA counts", "F5": "cell states within donors", "F6": "development data or a control",
    "F7": "data not public as counts", "F8": "no donor ids", "F9": "too few donors", "F10": "duplicate",
}
PAPERS_MIN, PAPERS_MAX, PAPERS_SEED = 20, 30, 20261020  # the protocol's section 9
MIN_REPLICATES = 3  # the engine's default minimum of replicates per group (report.DEFAULT_PREREG's min_replicates)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def as_text(r: dict) -> str:
    """One record as the screener reads it: its number, title, journal, year, types and abstract."""
    wrap = textwrap.TextWrapper(width=WIDTH, break_long_words=True, break_on_hyphens=False)
    meta = (f"JOURNAL: {r.get('journal') or '-'} | YEAR: {r.get('year') or '-'} | TYPES: "
            f"{'; '.join(r.get('pub_types') or []) or '-'} | PREPRINT: {'yes' if r.get('preprint') else 'no'}")
    parts = [f"=== {r['rid']} ===", *wrap.wrap("TITLE: " + " ".join(str(r.get("title") or "-").split())), meta]
    parts += wrap.wrap("ABSTRACT: " + (" ".join(str(r.get("abstract") or "").split()) or "(none)"))
    return "\n".join(parts) + "\n"


def screened_text(r: dict) -> str:
    """What the screener read of a record, without its number."""
    return as_text({k: r.get(k) for k in FIELDS}).split("\n", 1)[1]


def _write_batch(out: Path, name: str, records: list) -> dict:
    with open(out / name, "w") as f, open(out / name.replace(".jsonl", ".txt"), "w") as t:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
            t.write(as_text(r) + "\n")
    return dict(name=name, records=len(records), first=records[0]["rid"], last=records[-1]["rid"], sha256=sha256(out / name),
                text_sha256=sha256(out / name.replace(".jsonl", ".txt")))


def _batch_records(batches: Path) -> list:
    manifest = json.loads((batches / "manifest.json").read_text())
    return [json.loads(line) for b in manifest["batches"] for line in (batches / b["name"]).read_text().splitlines()
            if line.strip()]


def make_batches(search: Path, out: Path, size: int = BATCH) -> dict:
    sums = {}
    for line in (search / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        sums[name.lstrip("*").removeprefix("./")] = digest
    for name in ("records.jsonl.gz", "sample.txt"):
        got = sha256(search / name)
        if sums.get(name) != got:
            raise SystemExit(f"{name}: SHA256 {got} is not the search's {sums.get(name)}")
    with gzip.open(search / "records.jsonl.gz", "rt") as f:
        recs = {r["rid"]: r for r in map(json.loads, f)}
    sample = (search / "sample.txt").read_text().split()
    missing = [r for r in sample if r not in recs]
    if missing:
        raise SystemExit(f"{len(missing)} records of the sample are not in the list: {missing[:5]}")
    rids = sorted(sample)
    out.mkdir(parents=True, exist_ok=True)
    files = [_write_batch(out, f"batch_{k // size + 1:02d}.jsonl", [{key: recs[rid].get(key) for key in FIELDS}
                                                                    for rid in rids[k:k + size]])
             for k in range(0, len(rids), size)]
    manifest = dict(records_sha256=sums["records.jsonl.gz"], sample_sha256=sums["sample.txt"], n_list=len(recs),
                    n_screened=len(rids), batch_size=size, batches=files)
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    return manifest


def read_batch(path: Path) -> list:
    return [json.loads(line)["rid"] for line in path.read_text().splitlines() if line.strip()]


def read_decisions(path: Path, rids: list) -> tuple[dict, list]:
    """A pass's decisions for one batch, and the problems with them (screening.md's format)."""
    problems, out = [], {}
    lines = path.read_text().splitlines() if path.exists() else []
    if not lines:
        return out, [f"{path.name}: missing or empty"]
    if lines[0].split("\t") != HEADER:
        problems.append(f"{path.name}: the header is {lines[0]!r}")
    rows = [line.split("\t") for line in lines[1:] if line.strip()]
    if [r[0] for r in rows] != rids:
        seen = [r[0] for r in rows]
        problems.append(f"{path.name}: the rids are not the batch's in order ({len(seen)} rows for {len(rids)} records; "
                        f"missing {sorted(set(rids) - set(seen))[:5]}, extra {sorted(set(seen) - set(rids))[:5]})")
    for r in rows:
        if len(r) != 5:
            problems.append(f"{path.name}: {r[0]}: {len(r)} fields")
            continue
        rid, decision, code, kinds, note = r
        if decision not in DECISIONS:
            problems.append(f"{path.name}: {rid}: decision {decision!r}")
        if (decision == "EXCLUDE") != (code in CODES) or (decision != "EXCLUDE" and code != "-"):
            problems.append(f"{path.name}: {rid}: code {code!r} with {decision}")
        ks = kinds.split(",") if kinds != "-" else []
        if any(k not in KINDS for k in ks) or len(set(ks)) != len(ks):
            problems.append(f"{path.name}: {rid}: kinds {kinds!r}")
        if not note.strip():
            problems.append(f"{path.name}: {rid}: no note")
        out[rid] = dict(decision=decision, code=code, kinds=ks, note=note.strip())
    return out, problems


def read_pass(batches: Path, decisions: Path) -> tuple[dict, list]:
    manifest = json.loads((batches / "manifest.json").read_text())
    allp, alld = [], {}
    for b in manifest["batches"]:
        rids = read_batch(batches / b["name"])
        d, p = read_decisions(decisions / b["name"].replace(".jsonl", ".tsv"), rids)
        alld.update(d)
        allp += p
    return alld, allp


def cohen_kappa(a: list, b: list, categories=DECISIONS) -> float:
    """Cohen's kappa of two raters' labels over the same items."""
    n = len(a)
    if n == 0 or n != len(b):
        return float("nan")
    po = sum(x == y for x, y in zip(a, b)) / n
    ca, cb = Counter(a), Counter(b)
    pe = sum(ca[c] * cb[c] for c in categories) / n ** 2
    return (po - pe) / (1 - pe) if pe < 1 else float("nan")


def merge(batches: Path, pass1: Path, pass2: Path, out: Path) -> dict:
    manifest = json.loads((batches / "manifest.json").read_text())
    d1, p1 = read_pass(batches, pass1)
    d2, p2 = read_pass(batches, pass2)
    if p1 or p2:
        raise SystemExit("the passes' files do not check:\n" + "\n".join(p1 + p2))
    rids = [rid for b in manifest["batches"] for rid in read_batch(batches / b["name"])]
    recs = {}
    for b in manifest["batches"]:
        for line in (batches / b["name"]).read_text().splitlines():
            r = json.loads(line)
            recs[r["rid"]] = r
    a, b = [d1[r]["decision"] for r in rids], [d2[r]["decision"] for r in rids]
    both_ex = [r for r in rids if d1[r]["decision"] == d2[r]["decision"] == "EXCLUDE"]
    res = dict(n_screened=len(rids), kappa=cohen_kappa(a, b),
               agreement=sum(x == y for x, y in zip(a, b)) / len(rids),
               kappa_keep_or_unclear_vs_exclude=cohen_kappa(["EXCLUDE" if x == "EXCLUDE" else "KEEP" for x in a],
                                                            ["EXCLUDE" if x == "EXCLUDE" else "KEEP" for x in b],
                                                            ("KEEP", "EXCLUDE")),
               code_agreement_where_both_exclude=(sum(d1[r]["code"] == d2[r]["code"] for r in both_ex) / len(both_ex)
                                                  if both_ex else float("nan")),
               table={f"pass {k}": {**Counter(d[r]["decision"] for r in rids),
                                   **{c: sum(d[r]["code"] == c for r in rids) for c in CODES}}
                      for k, d in ((1, d1), (2, d2))},
               crosstab={f"{x} / {y}": n for (x, y), n in sorted(Counter(zip(a, b)).items())})
    both_keep = [r for r in rids if d1[r]["decision"] == d2[r]["decision"] == "KEEP"]
    owner = [r for r in rids if d1[r]["decision"] != d2[r]["decision"] or "UNCLEAR" in (d1[r]["decision"], d2[r]["decision"])]
    fulltext = [r for r in rids if {d1[r]["decision"], d2[r]["decision"]} & {"KEEP", "UNCLEAR"}]
    sets = {f"pass {k}": dict(code_X6=sum(d[r]["code"] == "X6" for r in rids), code_X7=sum(d[r]["code"] == "X7" for r in rids),
                              kind_set=sum("set" in d[r]["kinds"] for r in rids),
                              kind_global=sum("global" in d[r]["kinds"] for r in rids))
            for k, d in ((1, d1), (2, d2))}
    sets["both passes"] = dict(code_X6=sum(d1[r]["code"] == d2[r]["code"] == "X6" for r in rids),
                               code_X7=sum(d1[r]["code"] == d2[r]["code"] == "X7" for r in rids),
                               either_X6_or_X7=sum(bool({d1[r]["code"], d2[r]["code"]} & {"X6", "X7"}) for r in rids))
    res.update(both_keep=len(both_keep), to_owner=len(owner), fulltext=len(fulltext), sets=sets)
    out.mkdir(parents=True, exist_ok=True)
    (out / "merge.json").write_text(json.dumps(res, indent=1) + "\n")

    def row(r):
        x, y = d1[r], d2[r]
        return "\t".join([r, x["decision"], x["code"], ",".join(x["kinds"]) or "-", x["note"], y["decision"], y["code"],
                          ",".join(y["kinds"]) or "-", y["note"], str(recs[r].get("pmid") or ""),
                          str(recs[r].get("doi") or ""), (recs[r].get("title") or "").replace("\t", " ")])
    head = "rid\tdecision 1\tcode 1\tkinds 1\tnote 1\tdecision 2\tcode 2\tkinds 2\tnote 2\tpmid\tdoi\ttitle"
    (out / "all.tsv").write_text(head + "\n" + "\n".join(row(r) for r in rids) + "\n")
    (out / "to_owner.tsv").write_text(head + "\n" + "\n".join(row(r) for r in owner) + "\n")
    (out / "fulltext.txt").write_text("\n".join(fulltext) + "\n")
    return res


def carry(old_batches: Path, pass1: Path, pass2: Path, batches: Path, out: Path, size: int = BATCH) -> dict:
    """The first list's decisions on the corrected list where a record's screened text is unchanged; the
    others as batches to screen (`rescreen/`)."""
    old = _batch_records(old_batches)
    new = _batch_records(batches)
    by_text = {}
    for r in old:
        by_text.setdefault(screened_text(r), []).append(r["rid"])
    if any(len(v) > 1 for v in by_text.values()):
        raise SystemExit("records of the first list share a text: the carry would be ambiguous")
    old_of = {r["rid"]: (by_text.get(screened_text(r)) or [None])[0] for r in new}
    new_of = {o: n for n, o in old_of.items() if o}
    if len(new_of) != sum(1 for o in old_of.values() if o):
        raise SystemExit("one record of the first list matches several of the corrected list")
    decided = {}
    for k, d in ((1, pass1), (2, pass2)):
        dec, problems = read_pass(old_batches, d)
        if problems:
            raise SystemExit(f"pass {k}'s files do not check:\n" + "\n".join(problems))
        decided[k] = dec
    out.mkdir(parents=True, exist_ok=True)
    untranslated, translated = set(), 0

    def translate(note):
        def one(m):
            nonlocal translated
            if m.group(0) in new_of:
                translated += 1
                return new_of[m.group(0)]
            untranslated.add(m.group(0))
            return m.group(0)
        return re.sub(r"\bR\d{5}\b", one, note)
    carried = [r["rid"] for r in new if old_of[r["rid"]]]
    for k in (1, 2):
        rows = ["\t".join(HEADER)]
        for rid in carried:
            x = decided[k][old_of[rid]]
            rows.append("\t".join([rid, x["decision"], x["code"], ",".join(x["kinds"]) or "-", translate(x["note"])]))
        (out / f"carried_pass{k}.tsv").write_text("\n".join(rows) + "\n")
    (out / "carried.tsv").write_text("rid\tfirst_list_rid\n" + "".join(f"{r['rid']}\t{old_of[r['rid']] or '-'}\n" for r in new))
    left = [r for r in new if not old_of[r["rid"]]]
    re_dir = out / "rescreen"
    re_dir.mkdir(exist_ok=True)
    files = [_write_batch(re_dir, f"batch_{k // size + 1:02d}.jsonl", left[k:k + size]) for k in range(0, len(left), size)]
    man_new = json.loads((batches / "manifest.json").read_text())
    (re_dir / "manifest.json").write_text(json.dumps(dict(
        records_sha256=man_new["records_sha256"], sample_sha256=man_new["sample_sha256"], n_list=man_new["n_list"],
        n_screened=len(left), batch_size=size, batches=files), indent=1) + "\n")
    res = dict(first_list=len(old), corrected_list=len(new), carried=len(carried), to_screen=len(left),
               first_list_not_carried=sorted(set(r["rid"] for r in old) - set(new_of)),
               note_numbers_translated=translated, note_numbers_untranslated=sorted(untranslated),
               rescreen_batches=len(files))
    (out / "carry.json").write_text(json.dumps(res, indent=1) + "\n")
    return res


def combine(batches: Path, carried: Path, pass1: Path, pass2: Path, out: Path) -> dict:
    """Each pass's decisions on the corrected list (the carried and the rescreen's), per batch of the list."""
    manifest = json.loads((batches / "manifest.json").read_text())
    rows = [line.split("\t") for line in (carried / "carried.tsv").read_text().splitlines()[1:] if line.strip()]
    kept = [rid for rid, o in rows if o != "-"]
    res = {}
    for k, d in ((1, pass1), (2, pass2)):
        dec, problems = read_decisions(carried / f"carried_pass{k}.tsv", kept)
        more, p2 = read_pass(carried / "rescreen", d)
        if problems or p2:
            raise SystemExit(f"pass {k}'s files do not check:\n" + "\n".join(problems + p2))
        if set(dec) & set(more):
            raise SystemExit(f"pass {k}: records both carried and screened again")
        dec.update(more)
        (out / f"pass{k}").mkdir(parents=True, exist_ok=True)
        n = 0
        for b in manifest["batches"]:
            rids = read_batch(batches / b["name"])
            missing = [r for r in rids if r not in dec]
            if missing:
                raise SystemExit(f"pass {k}: no decision for {len(missing)} records of {b['name']}: {missing[:5]}")
            lines = ["\t".join(HEADER)] + ["\t".join([r, dec[r]["decision"], dec[r]["code"], ",".join(dec[r]["kinds"]) or "-",
                                                     dec[r]["note"]]) for r in rids]
            (out / f"pass{k}" / b["name"].replace(".jsonl", ".tsv")).write_text("\n".join(lines) + "\n")
            n += len(rids)
        res[f"pass {k}"] = dict(carried=len(kept), screened_again=len(more), decisions=n)
    return res


def fulltext_list(search: Path, rids_file: Path, out: Path) -> int:
    """The ids `search.py fulltext` needs for the records named in `rids_file` (one rid a line)."""
    with gzip.open(search / "records.jsonl.gz", "rt") as f:
        recs = {r["rid"]: r for r in map(json.loads, f)}
    rids = [r for r in rids_file.read_text().split() if r]
    rows = ["rid\tpmid\tpmcid\tdoi\tepmc_ids"]
    for rid in rids:
        r = recs[rid]
        rows.append("\t".join([rid, (r.get("pmids") or [r.get("pmid") or ""])[0] or "", r.get("pmcid") or "",
                               r.get("doi") or "", ",".join(r.get("epmc_ids") or [])]))
    out.write_text("\n".join(rows) + "\n")
    return len(rids)


def _flat(el) -> str:
    return " ".join("".join(el.itertext()).split())


def jats_text(xml: bytes) -> str:
    """A JATS full text as a screener reads it, in short lines; the references left out."""
    root = ET.fromstring(xml)
    for parent in list(root.iter()):
        for child in list(parent):
            if child.tag in ("ref-list", "ref"):
                parent.remove(child)
    wrap = textwrap.TextWrapper(width=WIDTH, break_long_words=True, break_on_hyphens=False)
    out = []

    def para(text, prefix=""):
        text = " ".join(text.split())
        if text:
            out.extend(wrap.wrap(prefix + text))

    def walk(el, depth):
        for c in el:
            if c.tag == "title":
                continue
            if c.tag in ("sec", "app", "ack", "notes", "fn-group", "glossary"):
                t = c.find("title")
                out.append("")
                out.append("#" * min(depth + 2, 5) + " " + (_flat(t) if t is not None else c.get("sec-type") or c.tag))
                walk(c, depth + 1)
            elif c.tag in ("fig", "table-wrap", "supplementary-material"):
                label = c.find("label")
                cap = c.find("caption")
                para((_flat(cap) if cap is not None else ""), f"[{_flat(label) if label is not None else c.tag}] ")
                for row in c.iter("tr"):
                    para(" | ".join(_flat(cell) for cell in row if cell.tag in ("td", "th")), "  ")
            elif c.tag in ("p", "list", "disp-quote", "statement", "def-list", "boxed-text"):
                para(_flat(c))
            else:
                walk(c, depth)

    meta = root.find(".//article-meta")
    title = meta.find(".//article-title") if meta is not None else None
    out.append("# " + (_flat(title) if title is not None else "(no title)"))
    for ab in (meta.findall("abstract") if meta is not None else []):
        out.append("")
        out.append("## Abstract" + (f" ({ab.get('abstract-type')})" if ab.get("abstract-type") else ""))
        walk(ab, 1)
        if not any(c.tag in ("p", "sec") for c in ab):
            para(_flat(ab))
    for part in ("body", "back"):
        el = root.find(part)
        if el is not None:
            out.append("")
            out.append(f"## [{part}]")
            walk(el, 1)
    return "\n".join(out) + "\n"


def fulltext_text(fulltext: Path, out: Path) -> dict:
    index = [line.split("\t") for line in (fulltext / "fulltext_index.tsv").read_text().splitlines()[1:] if line.strip()]
    out.mkdir(parents=True, exist_ok=True)
    done, failed = 0, []
    for rid, _tried, fid, file, _note in index:
        if file == "-":
            continue
        with gzip.open(fulltext / file, "rb") as f:
            xml = f.read()
        try:
            text = jats_text(xml)
        except ET.ParseError as e:
            failed.append(f"{rid}: {e}")
            continue
        (out / f"{rid}.txt").write_text(f"=== {rid} ({fid}) ===\n" + text)
        done += 1
    return dict(converted=done, failed=failed, listed=len(index))


def fulltext_batches(merge_dir: Path, which: int, texts: Path, out: Path, size: int = FT_BATCH) -> dict:
    """Pass `which`'s kept records with a full text, `size` to a batch (`batch_NN.list`: rid, text file), and
    those without one (`without_fulltext.txt`)."""
    rows = [line.split("\t") for line in (merge_dir / "all.tsv").read_text().splitlines()[1:] if line.strip()]
    col = 1 if which == 1 else 5
    kept = [r[0] for r in rows if r[col] == "KEEP"]
    have = [rid for rid in kept if (texts / f"{rid}.txt").exists()]
    without = [rid for rid in kept if rid not in set(have)]
    out.mkdir(parents=True, exist_ok=True)
    names = []
    for k in range(0, len(have), size):
        name = f"batch_{k // size + 1:02d}.list"
        (out / name).write_text("".join(f"{rid}\t{(texts / f'{rid}.txt').resolve()}\n" for rid in have[k:k + size]))
        names.append(name)
    (out / "without_fulltext.txt").write_text("".join(f"{rid}\n" for rid in without))
    return dict(kept=len(kept), with_fulltext=len(have), without=len(without), batches=names)


def check_fulltext(batch: Path, decisions: Path) -> list:
    rids = [line.split("\t")[0] for line in batch.read_text().splitlines() if line.strip()]
    lines = decisions.read_text().splitlines() if decisions.exists() else []
    if not lines:
        return [f"{decisions.name}: missing or empty"]
    problems = [] if lines[0].split("\t") == FT_HEADER else [f"{decisions.name}: the header is {lines[0]!r}"]
    rows = [line.split("\t") for line in lines[1:] if line.strip()]
    order, last, claims = [], None, 0
    for r in rows:
        if len(r) != len(FT_HEADER):
            problems.append(f"{r[0]}: {len(r)} fields, not {len(FT_HEADER)}")
            continue
        rid, decision, code, claim = r[:4]
        if decision == "CLAIM":
            claims += 1
            if last is None or last[0] != rid or last[1] != "KEEP":
                problems.append(f"{rid}: a claim line without its KEEP line just before")
            if claim != str(claims) or claims > 5:
                problems.append(f"{rid}: claim {claim!r} (claims are numbered 1 to 5 in order)")
            if code != "-" or r[14] not in ("yes", "no") or any(not x.strip() for x in r[4:16]):
                problems.append(f"{rid}: claim {claim}: code, rule4 or an empty field")
            continue
        if last is not None and last[1] == "KEEP" and claims == 0:
            problems.append(f"{last[0]}: KEEP without a claim line")
        last, claims = (rid, decision), 0
        order.append(rid)
        if decision not in DECISIONS:
            problems.append(f"{rid}: decision {decision!r}")
        elif (decision == "EXCLUDE") != (code in FT_CODES) or (decision != "EXCLUDE" and code != "-"):
            problems.append(f"{rid}: code {code!r} with {decision}")
        if claim != "0" or not r[16].strip():
            problems.append(f"{rid}: a decision line needs claim 0 and a note")
    if last is not None and last[1] == "KEEP" and claims == 0:
        problems.append(f"{last[0]}: KEEP without a claim line")
    if order != rids:
        problems.append(f"the decisions' rids {order[:3]}... are not the batch's {rids[:3]}... in order "
                        f"({len(order)} for {len(rids)})")
    return problems


def read_fulltext_pass(ft_dir: Path) -> tuple[dict, dict, list]:
    """A pass's full-text decisions (`ft_dir`: its batch lists and, beside each, `batch_NN.tsv`): per rid the
    decision line and the claim lines, and the problems found."""
    decisions, claims, problems = {}, {}, []
    for lst in sorted(ft_dir.glob("batch_*.list")):
        tsv = lst.with_suffix(".tsv")
        problems += [f"{tsv.name}: {p}" for p in check_fulltext(lst, tsv)]
        if not tsv.exists():
            continue
        for line in tsv.read_text().splitlines()[1:]:
            r = line.split("\t")
            if len(r) != len(FT_HEADER):
                continue
            row = dict(zip(FT_HEADER, r))
            if row["decision"] == "CLAIM":
                claims.setdefault(row["rid"], []).append(row)
            else:
                decisions[row["rid"]] = row
    return decisions, claims, problems


def merge_fulltext(merge_dir: Path, ft1: Path, ft2: Path, out: Path) -> dict:
    d1, c1, p1 = read_fulltext_pass(ft1)
    d2, c2, p2 = read_fulltext_pass(ft2)
    if p1 or p2:
        raise SystemExit("the passes' full-text files do not check:\n" + "\n".join(p1 + p2))
    rows = [line.split("\t") for line in (merge_dir / "all.tsv").read_text().splitlines()[1:] if line.strip()]
    ta = {r[0]: (r[1], r[5]) for r in rows}
    titles = {r[0]: r[11] for r in rows}
    without = {k: [x for x in (ft / "without_fulltext.txt").read_text().split() if x]
               for k, ft in ((1, ft1), (2, ft2))}
    both = sorted(set(d1) & set(d2))
    a, b = [d1[r]["decision"] for r in both], [d2[r]["decision"] for r in both]
    counts = {f"pass {k}": {**Counter(x["decision"] for x in d.values()),
                           **{c: sum(x["code"] == c for x in d.values()) for c in FT_CODES},
                           "kept at titles and abstracts": sum(ta[r][k - 1] == "KEEP" for r in ta),
                           "without a full text": len(without[k])}
              for k, d in ((1, d1), (2, d2))}
    kept = [r for r in both if d1[r]["decision"] == d2[r]["decision"] == "KEEP"]
    owner = sorted(set(r for r in both if d1[r]["decision"] != d2[r]["decision"] or "UNCLEAR" in (d1[r]["decision"],
                                                                                                    d2[r]["decision"]))
                   | (set(d1) ^ set(d2)) | set(without[1]) | set(without[2]))
    res = dict(read_by_both=len(both), kappa=cohen_kappa(a, b),
               agreement=(sum(x == y for x, y in zip(a, b)) / len(both)) if both else float("nan"),
               crosstab={f"{x} / {y}": n for (x, y), n in sorted(Counter(zip(a, b)).items())},
               counts=counts, kept_by_both=len(kept), to_owner=len(owner),
               sets_or_global={f"pass {k}": sum(x["code"] == "F2" for x in d.values()) for k, d in ((1, d1), (2, d2))})
    out.mkdir(parents=True, exist_ok=True)
    (out / "ft_merge.json").write_text(json.dumps(res, indent=1) + "\n")

    def dline(k, d, rid):
        x = d.get(rid)
        if x is None:
            return f"pass {k}: not read" + (" (no full text)" if rid in without[k] else
                                             f" (titles and abstracts: {ta.get(rid, ('?', '?'))[k - 1]})")
        return f"pass {k}: {x['decision']} {x['code']} — {x['note']}"
    cols = [c for c in FT_HEADER if c not in ("decision", "code")]
    with (out / "ft_kept.tsv").open("w") as f:
        f.write("pass\t" + "\t".join(cols) + "\ttitle\n")
        for rid in kept:
            for k, c in ((1, c1), (2, c2)):
                for cl in c.get(rid, []):
                    f.write(f"{k}\t" + "\t".join(cl[x] for x in cols) + f"\t{titles.get(rid, '')}\n")
    with (out / "ft_to_owner.txt").open("w") as f:
        for rid in owner:
            f.write(f"{rid}  {titles.get(rid, '')}\n  {dline(1, d1, rid)}\n  {dline(2, d2, rid)}\n")
            for k, c in ((1, c1), (2, c2)):
                for cl in c.get(rid, []):
                    f.write(f"    pass {k} claim {cl['claim']}: {cl['genes']} | {cl['conditions']} | {cl['direction']} | "
                            f"{cl['metric']} | {cl['donors']} (rule 4 {cl['rule4']}) | {cl['data']}\n")
            f.write("\n")
    return res


def donor_counts(text: str):
    """The numbers a claim's `donors` field names, and whether it says the same donors are in both conditions;
    None where it names no number (the owner reads it). A hint for the owner, not a parse of the design."""
    t = re.sub(r"\bn\s*[12]\b", " ", text.lower())
    nums = [int(x) for x in re.findall(r"(?<![\w.])(\d+)(?![\w.])", t)]
    paired = any(w in t for w in ("paired", "same donors", "same animals", "same patients", "both conditions"))
    return (nums, paired) if nums else None


def claims_list(merge_dir: Path, ft_merge: Path, ft1: Path, ft2: Path, out: Path) -> dict:
    """The list for the owner (protocol, section 11.5): the counts of both steps, the sets and global measures
    (section 1), the papers both passes keep with each pass's claims (more than 30: 30 drawn, section 9), the
    designs the engine cannot run, the verdicts the papers give, and what the owner decides."""
    import numpy as np
    ta = json.loads((merge_dir / "merge.json").read_text())
    ft = json.loads((ft_merge / "ft_merge.json").read_text())
    rows = {r[0]: r for r in (line.split("\t") for line in (merge_dir / "all.tsv").read_text().splitlines()[1:]) if r[0]}
    d1, c1, _ = read_fulltext_pass(ft1)
    d2, c2, _ = read_fulltext_pass(ft2)
    kept = sorted(r for r in set(d1) & set(d2) if d1[r]["decision"] == d2[r]["decision"] == "KEEP")
    drawn = (sorted(kept[i] for i in np.random.default_rng(PAPERS_SEED).choice(len(kept), PAPERS_MAX, replace=False))
             if len(kept) > PAPERS_MAX else list(kept))
    small, unread, verdicts = [], [], []
    for rid in kept:
        for k, c in ((1, c1), (2, c2)):
            for cl in c.get(rid, []):
                n = donor_counts(cl["donors"])
                if n is None:
                    unread.append((rid, k, cl))
                elif not n[1] and min(n[0]) < MIN_REPLICATES:
                    small.append((rid, k, cl))
                if cl["verdict"].strip() not in ("-", ""):
                    verdicts.append((rid, k, cl))
    res = dict(titles_and_abstracts=dict(n=ta["n_screened"], kappa=ta["kappa"], agreement=ta["agreement"],
                                         both_keep=ta["both_keep"], to_owner=ta["to_owner"], fulltext=ta["fulltext"]),
               full_texts=dict(read_by_both=ft["read_by_both"], kappa=ft["kappa"], kept_by_both=ft["kept_by_both"],
                               to_owner=ft["to_owner"]),
               sets=ta["sets"], sets_at_full_text=ft["sets_or_global"], papers_kept=len(kept), papers_listed=len(drawn),
               drawn=len(kept) > PAPERS_MAX, fewer_than_target=len(kept) < PAPERS_MIN,
               claims={f"pass {k}": sum(len(c.get(r, [])) for r in kept) for k, c in ((1, c1), (2, c2))},
               small_designs=len({r for r, *_ in small}), donors_unread=len({r for r, *_ in unread}),
               verdicts_found=len({r for r, *_ in verdicts}))
    out.mkdir(parents=True, exist_ok=True)
    (out / "claims.json").write_text(json.dumps(dict(res, kept=kept, listed=drawn), indent=1) + "\n")

    def cl_line(k, cl):
        return (f"  - pass {k}, claim {cl['claim']}: **{cl['genes']}** ({cl['species']}; {cl['cells']}); {cl['conditions']}; "
                f"{cl['direction']}; {cl['metric']}; {cl['where']}; data {cl['data']}; donors {cl['donors']} (rule 4 "
                f"{cl['rule4']}); verdict {cl['verdict']}\n    > {cl['quote']}")
    L = ["# Level B: the claims for the owner", "",
         f"Titles and abstracts: {ta['n_screened']} records, kappa {ta['kappa']:.3f}, agreement {ta['agreement']:.4f}; "
         f"kept by both {ta['both_keep']}, to the owner {ta['to_owner']}, full texts asked {ta['fulltext']}.",
         f"Full texts: read by both passes {ft['read_by_both']}, kappa {ft['kappa']:.3f}; kept by both {ft['kept_by_both']}, "
         f"to the owner {ft['to_owner']}.", "",
         "## The counts per pass", "", "| step | pass 1 | pass 2 |", "|---|---|---|"]
    for key in ("KEEP", "UNCLEAR", "EXCLUDE"):
        L.append(f"| titles and abstracts: {key} | {ta['table']['pass 1'].get(key, 0)} | {ta['table']['pass 2'].get(key, 0)} |")
    for c in CODES:
        L.append(f"| titles and abstracts: {c} ({CODE_NAMES[c]}) | {ta['table']['pass 1'].get(c, 0)} | "
                 f"{ta['table']['pass 2'].get(c, 0)} |")
    for key in ("KEEP", "UNCLEAR", "EXCLUDE", "without a full text"):
        L.append(f"| full texts: {key} | {ft['counts']['pass 1'].get(key, 0)} | {ft['counts']['pass 2'].get(key, 0)} |")
    for c in FT_CODES:
        L.append(f"| full texts: {c} ({FT_CODE_NAMES[c]}) | {ft['counts']['pass 1'].get(c, 0)} | {ft['counts']['pass 2'].get(c, 0)} |")
    st, s1, s2 = ta["sets"], ta["sets"]["pass 1"], ta["sets"]["pass 2"]
    L += ["", "## Sets and global measures (protocol, section 1: counted, not kept)", "",
          f"Titles and abstracts: X6 (a set) {s1['code_X6']} and {s2['code_X6']}, both passes "
          f"{st['both passes']['code_X6']}; X7 (a global measure) {s1['code_X7']} and {s2['code_X7']}, both "
          f"{st['both passes']['code_X7']}; X6 or X7 in either pass {st['both passes']['either_X6_or_X7']} records. "
          f"Records whose kinds name a set: {s1['kind_set']} and {s2['kind_set']}; a global measure: {s1['kind_global']} and "
          f"{s2['kind_global']}. Full texts, F2 (sets or global measures only): {ft['sets_or_global']['pass 1']} and "
          f"{ft['sets_or_global']['pass 2']}.", "",
          f"## The papers both passes keep: {len(kept)}", ""]
    if len(kept) > PAPERS_MAX:
        L.append(f"More than {PAPERS_MAX}: {PAPERS_MAX} drawn by numpy.random.default_rng({PAPERS_SEED}).choice({len(kept)}, "
                 f"{PAPERS_MAX}, replace=False) over the kept papers in the order of their records; the others listed below.")
    if len(kept) < PAPERS_MIN:
        L.append(f"**Fewer than {PAPERS_MIN} papers with pairs** (section 9): the owner decides about change 9 from the counts "
                 "of sets above.")
    L.append("")
    for rid in drawn:
        r = rows.get(rid, [rid] + [""] * 11)
        L.append(f"### {rid}: {r[11]}")
        L.append(f"PMID {r[9] or '-'}, DOI {r[10] or '-'}")
        for k, c in ((1, c1), (2, c2)):
            L += [cl_line(k, cl) for cl in c.get(rid, [])]
        L.append("")
    if len(kept) > PAPERS_MAX:
        L += ["Kept, not drawn: " + ", ".join(r for r in kept if r not in set(drawn)), ""]
    L += [f"## Two groups whose donors field names a number below {MIN_REPLICATES} (the engine's minimum of replicates per "
          "group; a 2-against-8 design passes rule 4 but not the engine): listed apart, the owner checks each", ""]
    L += [f"- {rid}, " + cl_line(k, cl).strip()[2:] for rid, k, cl in small] or ["- none"]
    L += ["", "## Claims whose donors field gives no numbers (the owner reads them)", ""]
    L += [f"- {rid}, pass {k}, claim {cl['claim']}: {cl['donors']}" for rid, k, cl in unread] or ["- none"]
    L += ["", "## The verdicts the papers themselves give (section 7; the search of section 7 is not yet run)", ""]
    L += [f"- {rid}, pass {k}, claim {cl['claim']} ({cl['genes']}): {cl['verdict']}" for rid, k, cl in verdicts] or ["- none"]
    L += ["", "## To the owner", "", f"At titles and abstracts: `{merge_dir.name}/to_owner.tsv` (both passes' notes). At the full "
          "texts:", ""]
    L += ["```", (ft_merge / "ft_to_owner.txt").read_text().rstrip(), "```", ""]
    (out / "claims.md").write_text("\n".join(L) + "\n")
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)
    b = sub.add_parser("batches")
    b.add_argument("--search", type=Path, required=True)
    b.add_argument("--out", type=Path, required=True)
    b.add_argument("--size", type=int, default=BATCH)
    c = sub.add_parser("check")
    c.add_argument("--batches", type=Path, required=True)
    c.add_argument("--decisions", type=Path, required=True)
    o = sub.add_parser("check-one")
    o.add_argument("--batch", type=Path, required=True, help="batch_NN.jsonl or batch_NN.txt")
    o.add_argument("--decisions", type=Path, required=True)
    fl = sub.add_parser("fulltext-list")
    fl.add_argument("--search", type=Path, required=True)
    fl.add_argument("--rids", type=Path, required=True)
    fl.add_argument("--out", type=Path, required=True)
    ft = sub.add_parser("fulltext-text")
    ft.add_argument("--fulltext", type=Path, required=True)
    ft.add_argument("--out", type=Path, required=True)
    fb = sub.add_parser("fulltext-batches")
    fb.add_argument("--merge", type=Path, required=True)
    fb.add_argument("--pass", dest="which", type=int, choices=(1, 2), required=True)
    fb.add_argument("--texts", type=Path, required=True)
    fb.add_argument("--out", type=Path, required=True)
    fb.add_argument("--size", type=int, default=FT_BATCH)
    cf = sub.add_parser("check-fulltext")
    cf.add_argument("--batch", type=Path, required=True)
    cf.add_argument("--decisions", type=Path, required=True)
    mf = sub.add_parser("merge-fulltext")
    mf.add_argument("--merge", type=Path, required=True)
    mf.add_argument("--ft1", type=Path, required=True)
    mf.add_argument("--ft2", type=Path, required=True)
    mf.add_argument("--out", type=Path, required=True)
    m = sub.add_parser("merge")
    m.add_argument("--batches", type=Path, required=True)
    m.add_argument("--pass1", type=Path, required=True)
    m.add_argument("--pass2", type=Path, required=True)
    m.add_argument("--out", type=Path, required=True)
    cy = sub.add_parser("carry")
    cy.add_argument("--old-batches", type=Path, required=True)
    cy.add_argument("--pass1", type=Path, required=True)
    cy.add_argument("--pass2", type=Path, required=True)
    cy.add_argument("--batches", type=Path, required=True)
    cy.add_argument("--out", type=Path, required=True)
    cl = sub.add_parser("claims")
    cl.add_argument("--merge", type=Path, required=True)
    cl.add_argument("--ft-merge", type=Path, required=True)
    cl.add_argument("--ft1", type=Path, required=True)
    cl.add_argument("--ft2", type=Path, required=True)
    cl.add_argument("--out", type=Path, required=True)
    cb = sub.add_parser("combine")
    cb.add_argument("--batches", type=Path, required=True)
    cb.add_argument("--carried", type=Path, required=True)
    cb.add_argument("--pass1", type=Path, required=True)
    cb.add_argument("--pass2", type=Path, required=True)
    cb.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.mode == "carry":
        print(json.dumps(carry(a.old_batches, a.pass1, a.pass2, a.batches, a.out), indent=1))
        return 0
    if a.mode == "claims":
        print(json.dumps(claims_list(a.merge, a.ft_merge, a.ft1, a.ft2, a.out), indent=1))
        return 0
    if a.mode == "combine":
        print(json.dumps(combine(a.batches, a.carried, a.pass1, a.pass2, a.out), indent=1))
        return 0
    if a.mode == "batches":
        man = make_batches(a.search, a.out, a.size)
        print(f"{man['n_screened']} of {man['n_list']} records in {len(man['batches'])} batches; records.jsonl.gz "
              f"{man['records_sha256']}")
        return 0
    if a.mode == "fulltext-list":
        print(f"{fulltext_list(a.search, a.rids, a.out)} records listed in {a.out}")
        return 0
    if a.mode == "fulltext-text":
        res = fulltext_text(a.fulltext, a.out)
        print(json.dumps(res, indent=1))
        return 1 if res["failed"] else 0
    if a.mode == "fulltext-batches":
        print(json.dumps(fulltext_batches(a.merge, a.which, a.texts, a.out, a.size), indent=1))
        return 0
    if a.mode == "check-fulltext":
        problems = check_fulltext(a.batch, a.decisions)
        for p in problems:
            print(p)
        print(f"{len(problems)} problems")
        return 1 if problems else 0
    if a.mode == "merge-fulltext":
        print(json.dumps(merge_fulltext(a.merge, a.ft1, a.ft2, a.out), indent=1))
        return 0
    if a.mode == "check-one":
        batch = a.batch.with_suffix(".jsonl")
        d, problems = read_decisions(a.decisions, read_batch(batch))
        for p in problems:
            print(p)
        print(f"{len(d)} decisions for {len(read_batch(batch))} records; {len(problems)} problems")
        return 1 if problems else 0
    if a.mode == "check":
        d, problems = read_pass(a.batches, a.decisions)
        for p in problems:
            print(p)
        print(f"{len(d)} decisions; {len(problems)} problems; " + ", ".join(f"{k} {v}" for k, v in
                                                                         sorted(Counter(x['decision'] for x in d.values()).items())))
        return 1 if problems else 0
    res = merge(a.batches, a.pass1, a.pass2, a.out)
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
