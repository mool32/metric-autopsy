#!/usr/bin/env python3
"""Level B's screening of titles and abstracts (`validation/plans/level_b_protocol.md`, section 5): the
batches the screeners read, the check of what each pass wrote, and the comparison of the two passes.

    python -I validation/level_b/screening.py batches --search <results/level-b-search> --out <dir>
    python -I validation/level_b/screening.py check --batches <dir> --decisions <pass dir>
    python -I validation/level_b/screening.py merge --batches <dir> --pass1 <dir> --pass2 <dir> --out <dir>

`batches` reads the search's merged list (`records.jsonl.gz`, checked against the search's SHA256SUMS) and its
sample (`sample.txt`, every record where the list holds at most 5,000), and writes the screened records in
the order of their numbers, `BATCH` to a batch, with the fields the screeners read (`screening.md`) and
nothing of how a record was found: `batch_NN.jsonl` for the tools, `batch_NN.txt` for the screener (the
same records as text in short lines). `check` checks a pass's decision files against `screening.md`'s
format; `check-one` one batch's file (the screener runs it before finishing).
`merge` compares the passes: Cohen's kappa over the three decisions, the counts of every step (a
PRISMA-style table) per pass, the records for the owner (decided differently by the two passes, or unclear
in either), the counts of claims about sets and global measures (section 1), and the records whose full text
is fetched (kept or unclear in either pass).
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
import textwrap
from collections import Counter
from pathlib import Path

BATCH = 200
FIELDS = ("rid", "title", "abstract", "journal", "year", "pub_types", "pmid", "doi", "preprint")
WIDTH = 110
HEADER = ["rid", "decision", "code", "kinds", "note"]
DECISIONS = ("KEEP", "EXCLUDE", "UNCLEAR")
CODES = tuple(f"X{i}" for i in range(1, 11))
KINDS = ("pair", "set", "global")
CODE_NAMES = {
    "X1": "not a paper with findings of its own", "X2": "no single-cell RNA counts",
    "X3": "development data or a control", "X4": "no changed coordination between conditions",
    "X5": "cell states within donors", "X6": "SET: a set's coordination", "X7": "GLOBAL: a global measure",
    "X8": "too few donors (the abstract shows it)", "X9": "data not public or not counts",
    "X10": "duplicate",
}


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
    files = []
    for k in range(0, len(rids), size):
        name = f"batch_{k // size + 1:02d}.jsonl"
        with open(out / name, "w") as f, open(out / name.replace(".jsonl", ".txt"), "w") as t:
            for rid in rids[k:k + size]:
                r = {key: recs[rid].get(key) for key in FIELDS}
                f.write(json.dumps(r, ensure_ascii=False) + "\n")
                t.write(as_text(r) + "\n")
        files.append(dict(name=name, records=len(rids[k:k + size]), first=rids[k], last=rids[min(k + size, len(rids)) - 1],
                          sha256=sha256(out / name), text_sha256=sha256(out / name.replace(".jsonl", ".txt"))))
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
    m = sub.add_parser("merge")
    m.add_argument("--batches", type=Path, required=True)
    m.add_argument("--pass1", type=Path, required=True)
    m.add_argument("--pass2", type=Path, required=True)
    m.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.mode == "batches":
        man = make_batches(a.search, a.out, a.size)
        print(f"{man['n_screened']} of {man['n_list']} records in {len(man['batches'])} batches; records.jsonl.gz "
              f"{man['records_sha256']}")
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
