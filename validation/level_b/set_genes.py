#!/usr/bin/env python3
"""Level B, `plan_s5.md` section 5: the kept sets' genes, taken from where each paper lists them.

* `fetch` (in Actions, once, by the label `level-b-setgenes`): for every kept paper with a set's claim
  (`amendment1/kept_claims.json`), Europe PMC's supplementary files of its open-access paper, by its PMCID, as they
  came (`raw/setgenes/<rid>-supplementary.zip.gz`), with `setgenes_index.tsv`; a paper without a PMCID, or whose
  supplement Europe PMC does not give, is listed so.
* `show`: each archive's members and, of each table, its sheets and first rows — what is looked at to write the
  specs, nothing more of the files.
* `fetch --extra extra.tsv` also fetches the files a list names beforehand (rid, URL), each as it came with its
  SHA256 (`raw/setgenes/<rid>-<file>.gz`, `setgenes_extra.tsv`): a set the paper names from a public collection
  without listing it (MSigDB's Hallmark sets: the collection's file of the release the paper used).
* `fetch --versions versions.tsv` fetches, instead, the supplement of another version of a kept record (its
  journal article where S3 read its preprint; `raw/setgenes/<rid>-<pmcid>-supplementary.zip.gz`); `show` and
  `extract` then take both runs' directories, the first run's first, and a spec names that archive `rid:pmcid`.
  A member of an archive inside the archive is named `outer.zip!/inner.xlsx`.
* `extract`: each set's genes by its spec (`set_specs.json`: the archive's member, its sheet, the row of its header,
  the column of the genes, the rows that are the set's; or an extra file and the set's name in it, a GMT; or the
  sentence of the paper's text that lists them, `--texts`; or "not listed" and why), their number against the
  number the paper states, and the kept claims with the sets' genes filled.
* `dry`: the three on fabricated answers (FABRICATED, nothing fetched).

    python -I validation/level_b/set_genes.py fetch --claims validation/level_b/amendment1/kept_claims.json --out OUT
    python -I validation/level_b/set_genes.py show --files OUT --out OUT/show.md
    python -I validation/level_b/set_genes.py extract --specs validation/level_b/set_specs.json --files OUT \\
        --claims validation/level_b/amendment1/kept_claims.json --out validation/level_b/amendment1/set_genes.json \\
        --claims-out validation/level_b/amendment1/kept_claims.json
    python -I validation/level_b/set_genes.py dry --out DIR
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import sys
import urllib.error
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import search as S  # noqa: E402  the search's network, raw store and request log, unchanged

TABLES = (".xlsx", ".csv", ".tsv", ".txt")
SHOW_ROWS = 5
_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


# --------------------------------------------------------------------------- #
# fetch
# --------------------------------------------------------------------------- #
def _error_message(body: bytes) -> str:
    m = re.search(rb"<errMsg>(.*?)</errMsg>", body)
    return m.group(1).decode("utf-8", "replace") if m else f"not an archive ({len(body)} bytes)"


def fetch_extra(net, out, rows: list) -> list:
    """The files named beforehand (`extra.tsv`: rid, url), each fetched once with the search's retries."""
    got = []
    for rid, url in rows:
        name = url.rsplit("/", 1)[-1]
        body = net.get(url, {})
        got.append(dict(rid=rid, url=url, file=out.raw(f"setgenes/{rid}-{name}", body), bytes=len(body),
                        sha256=hashlib.sha256(body).hexdigest()))
        out.log(f"{rid}: {name}, {len(body):,} bytes")
    with (out.root / "setgenes_extra.tsv").open("w") as f:
        cols = ["rid", "url", "file", "bytes", "sha256"]
        f.write("\t".join(cols) + "\n")
        f.writelines("\t".join(str(r[c]) for c in cols) + "\n" for r in got)
    return got


def read_extra(path: Path) -> list:
    rows = [line.split("\t") for line in path.read_text().splitlines() if line.strip() and not line.startswith("#")]
    if not rows or rows[0][:2] != ["rid", "url"]:
        raise SystemExit(f"{path}: the header must be rid<TAB>url")
    return [(r[0], r[1]) for r in rows[1:]]


def read_versions(path: Path) -> list:
    """Another version's PMCID per kept record (`rid<TAB>pmcid`, a header first): its supplement fetched instead."""
    return [tuple(line.split("\t")[:2]) for line in path.read_text().splitlines()[1:] if line.strip()]


def fetch(net, out, kept: dict, versions: list | None = None) -> list:
    rows = []
    wanted = ([dict(rid=r, pmcid=v, claims=[dict(kind="set")], version=True) for r, v in versions] if versions
              else kept["papers"])
    for p in wanted:
        if not any(c["kind"] == "set" for c in p["claims"]):
            continue
        row = dict(rid=p["rid"], pmcid=p.get("pmcid") or "-", status="-", file="-", bytes=0, sha256="-", note="-")
        rows.append(row)
        if not p.get("pmcid"):
            row["note"] = "no PMCID"
            continue
        stem = f"{p['rid']}-{p['pmcid']}" if p.get("version") else p["rid"]
        try:
            body = net.get(f"{S.EPMC}/{p['pmcid']}/supplementaryFiles", {})
        except urllib.error.HTTPError as e:
            row.update(status=e.code, note=f"HTTP {e.code}")
            continue
        except S.Unavailable as e:
            if e.status is None:  # the network itself: a technical failure, not an answer
                raise
            row.update(status=e.status, note=f"HTTP {e.status} on all 6 attempts")
            continue
        is_zip = body[:2] == b"PK"
        row.update(status=200, bytes=len(body), sha256=hashlib.sha256(body).hexdigest(),
                   file=out.raw(f"setgenes/{stem}-supplementary.{'zip' if is_zip else 'xml'}", body),
                   note="-" if is_zip else _error_message(body))
        out.log(f"{p['rid']} {p['pmcid']}: {'an archive' if is_zip else row['note']}, {len(body):,} bytes")
    with (out.root / "setgenes_index.tsv").open("w") as f:
        cols = ["rid", "pmcid", "status", "file", "bytes", "sha256", "note"]
        f.write("\t".join(cols) + "\n")
        f.writelines("\t".join(str(r[c]) for c in cols) + "\n" for r in rows)
    return rows


# --------------------------------------------------------------------------- #
# reading a supplement's tables
# --------------------------------------------------------------------------- #
def archives(files) -> dict:
    """The fetched archives (the gzip of each, as the fetch stored it), by rid (the first run's, the version S3
    read) and by `rid:pmcid` (any run's); `files` one results directory or several."""
    out = {}
    for d in ([files] if isinstance(files, Path) else files):
        for line in (d / "setgenes_index.tsv").read_text().splitlines()[1:]:
            r = line.split("\t")
            if len(r) >= 4 and r[3].endswith(".zip.gz"):
                out.setdefault(r[0], d / r[3])
                out[f"{r[0]}:{r[1]}"] = d / r[3]
    return out


def member(z: zipfile.ZipFile, name: str) -> bytes:
    """A member's bytes; `outer.zip!/inner.xlsx` names a member of an archive inside the archive."""
    if "!/" in name:
        outer, inner = name.split("!/", 1)
        with zipfile.ZipFile(io.BytesIO(z.read(outer))) as zz:
            return zz.read(inner)
    return z.read(name)


def _zip(path: Path) -> zipfile.ZipFile:
    with gzip.open(path, "rb") as g:
        return zipfile.ZipFile(io.BytesIO(g.read()))


def sheet_names(data: bytes) -> list:
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
    return [e.get("name") for e in wb.iter(_M + "sheet")]


def table_rows(name: str, data: bytes, sheet: str | None = None) -> list:
    """A table's rows as lists of strings: an Excel sheet (by data_run.py's reader, the standard library's), or a
    CSV, TSV or text file (a tab in its first line makes it tab-separated)."""
    if name.lower().endswith(".xlsx"):
        import tempfile
        from data_run import xlsx_rows
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "t.xlsx"
            p.write_bytes(data)
            return [["" if v is None else v for v in r] for r in xlsx_rows(p, sheet)]
    text = data.decode("utf-8-sig", "replace")
    first = text.split("\n", 1)[0]
    return [r for r in csv.reader(io.StringIO(text), delimiter="\t" if "\t" in first else ",")]


def _show_tables(L: list, name: str, data: bytes):
    sheets = sheet_names(data) if name.lower().endswith(".xlsx") else [None]
    for sh in sheets:
        try:
            rows = table_rows(name, data, sh)
        except Exception as e:  # a table this reader cannot read is said so
            L += ["", f"### `{name}`{' / ' + sh if sh else ''}: not read ({type(e).__name__}: {e})"]
            continue
        L += ["", f"### `{name}`{' / ' + sh if sh else ''}: {len(rows)} rows", "", "```"]
        L += ["\t".join(r)[:300] for r in rows[:SHOW_ROWS]] + ["```"]


def show(files, out: Path) -> str:
    """Each archive's members (and those of an archive inside it, one level), each table's sheets and first rows."""
    L = ["# The supplements' members, and each table's sheets and first rows (set_genes.py show)", ""]
    seen = set()
    for key, path in sorted(archives(files).items()):
        if path in seen:
            continue
        seen.add(path)
        z = _zip(path)
        L += [f"## {key} ({path.name})", ""]
        inner = {}
        for info in z.infolist():
            L.append(f"- `{info.filename}` ({info.file_size:,} bytes)")
            if info.filename.lower().endswith(".zip"):
                with zipfile.ZipFile(io.BytesIO(z.read(info.filename))) as zz:
                    inner[info.filename] = [(i.filename, i.file_size) for i in zz.infolist()]
                L += [f"  - `{info.filename}!/{n}` ({b:,} bytes)" for n, b in inner[info.filename]]
        names = [i.filename for i in z.infolist()] + [f"{o}!/{n}" for o, ms in inner.items() for n, _b in ms]
        for name in names:
            if name.lower().endswith(TABLES):
                _show_tables(L, name, member(z, name))
        L.append("")
    text = "\n".join(L) + "\n"
    out.write_text(text)
    return text


# --------------------------------------------------------------------------- #
# extract
# --------------------------------------------------------------------------- #
def extract_one(z: zipfile.ZipFile, spec: dict) -> list:
    rows = table_rows(spec["member"], member(z, spec["member"]), spec.get("sheet"))
    h = int(spec.get("header_row", 0))
    header, body = rows[h], rows[h + 1:]
    names = spec.get("column", 0)
    cols = [header.index(c) for c in (names if isinstance(names, list) else [names])] if names != 0 else [0]
    filt = {header.index(k): {str(x) for x in v} for k, v in (spec.get("rows") or {}).items()}
    genes = []
    for col in cols:  # several columns: their union, in the spec's order
        for r in body:
            if any((r[j] if j < len(r) else "") not in vals for j, vals in filt.items()):
                continue
            v = r[col].strip() if col < len(r) else ""
            for g in (v.split(spec["split"]) if spec.get("split") else [v]):
                g = g.strip()
                if g and g.lower() not in ("nan", "none", "na") and g not in genes:
                    genes.append(g)
    return genes


def gmt_set(data: bytes, name: str) -> list:
    """A set's genes from a GMT file (name, description, genes, tab-separated)."""
    for line in data.decode("utf-8", "replace").splitlines():
        f = line.split("\t")
        if f and f[0] == name:
            return [g for g in f[2:] if g.strip()]
    raise SystemExit(f"no set {name!r} in the GMT file")


def _norm(text: str) -> str:
    return " ".join(text.split())


def text_set(texts: Path, spec: dict) -> list:
    """A set the paper lists in its text (`plan_s5.md`, section 5, as precised on 2026-11-11): the sentence the spec
    quotes, found in the paper's text as S3 read it (its SHA256 checked), and the list it holds between `list_from`
    and `list_to`, cut at `split`."""
    path = texts / spec["text"]
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != spec["sha256"]:
        raise SystemExit(f"{spec['rid']} claim {spec['claim']}: {spec['text']} is not the text S3 read (SHA256)")
    quote = _norm(spec["quote"])
    if quote not in _norm(data.decode("utf-8")):
        raise SystemExit(f"{spec['rid']} claim {spec['claim']}: the quote is not in {spec['text']}")
    start = quote.index(spec["list_from"]) + len(spec["list_from"])
    inner = quote[start:quote.index(spec["list_to"], start)]
    genes = []
    for g in inner.split(spec["split"]):
        g = g.strip()
        if g and g not in genes:
            genes.append(g)
    return genes


def extract(specs: list, files, kept: dict, texts: Path | None = None) -> list:
    arch = archives(files)
    stated = {(p["rid"], c["id"]): c for p in kept["papers"] for c in p["claims"] if c["kind"] == "set"}
    res = []
    for s in specs:
        c = stated.get((s["rid"], str(s["claim"])))
        if c is None:
            raise SystemExit(f"{s['rid']} claim {s['claim']}: not a kept set's claim")
        note = None
        if s.get("not_listed"):
            genes, note = [], f"not listed: {s['not_listed']}"
        elif s.get("text"):
            if texts is None:
                raise SystemExit(f"{s['rid']} claim {s['claim']}: a text's list needs --texts")
            genes = text_set(texts, s)
        elif s.get("extra"):
            first = files if isinstance(files, Path) else files[0]
            with gzip.open(first / "raw" / "setgenes" / f"{s['rid']}-{s['extra']}.gz", "rb") as g:
                genes = gmt_set(g.read(), s["gmt_set"])
        else:
            genes = extract_one(_zip(arch[s.get("archive") or s["rid"]]), s)
        res.append(dict(rid=s["rid"], claim=str(s["claim"]), set_name=c.get("set_name"), stated=c.get("set_size"),
                        n=len(genes), genes=genes, note=note, spec=s))
    return res


def fill(kept: dict, res: list) -> dict:
    by = {(r["rid"], r["claim"]): r for r in res}
    for p in kept["papers"]:
        for c in p["claims"]:
            r = by.get((p["rid"], c["id"]))
            if r is not None:
                sp = r["spec"]
                c["genes"] = r["genes"]
                c["genes_from"] = (r["note"] if r.get("note") else
                                   f"set_genes.json ({sp.get('member') or sp.get('extra') or 'the text ' + str(sp.get('text'))})")
    return kept


# --------------------------------------------------------------------------- #
# the dry run
# --------------------------------------------------------------------------- #
class Fabricated:
    """Europe PMC's answers (FABRICATED): an archive with a CSV list for an even PMCID, "not open access" else."""

    def __init__(self, log):
        self.log, self.rows = log, []

    def get(self, url: str, params: dict, post: bool = False) -> bytes:
        if url.endswith(".gmt"):
            body = b"HALLMARK_FAB\thttp://fab\tGENE1\tGENE5\tGENE7\nHALLMARK_OTHER\thttp://fab\tGENE2\n"
            self.rows.append(dict(time="FABRICATED", method="GET", url=url, params="{}", status=200, bytes=len(body),
                                  sha256="-", attempt=1))
            return body
        pmcid = url.rsplit("/", 2)[-2]
        if int(pmcid[3:]) % 2 == 0:
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w") as z:
                z.writestr("supp/Table_S4.csv", "module,gene\nM3,GENE1\nM3,GENE2\nM4,GENE9\nM3,GENE3\n")
                z.writestr("supp/readme.txt", "FABRICATED\n")
            body = buf.getvalue()
        else:
            body = (b'<?xml version="1.0"?><errorBean><errCode>0</errCode><errMsg>Article with id ' + pmcid.encode()
                    + b' is not open access one</errMsg></errorBean>')
        self.rows.append(dict(time="FABRICATED", method="GET", url=url, params=S.shown_params(params), status=200,
                              bytes=len(body), sha256="-", attempt=1))
        return body


def dry(out_dir: Path) -> dict:
    out = S.Out(out_dir)
    out.log("# Level B, the sets' genes — FABRICATED answers (the dry run, nothing fetched)")
    kept = dict(papers=[
        dict(rid="R09001", pmcid="PMC100", claims=[dict(id="1", kind="set", genes=[], set_name="M3 module", set_size=3)]),
        dict(rid="R09002", pmcid="PMC101", claims=[dict(id="1", kind="set", genes=[], set_name="M7", set_size=10)]),
        dict(rid="R09003", pmcid="", claims=[dict(id="1", kind="set", genes=[], set_name="X", set_size=5)]),
        dict(rid="R09004", pmcid="PMC102", claims=[dict(id="1", kind="pair", genes=["A", "B"])])])
    kept["papers"][1]["claims"].append(dict(id="2", kind="set", genes=[], set_name="HALLMARK_FAB", set_size=3))
    net = Fabricated(out.log)
    rows = fetch(net, out, kept)
    extra = fetch_extra(net, out, [("R09002", "https://example.org/h.all.v0.Hs.symbols.gmt")])
    text = show(out_dir, out_dir / "show.md")
    res = extract([dict(rid="R09001", claim="1", member="supp/Table_S4.csv", column="gene", rows=dict(module=["M3"])),
                   dict(rid="R09002", claim="2", extra="h.all.v0.Hs.symbols.gmt", gmt_set="HALLMARK_FAB")],
                  out_dir, kept)
    fill(kept, res)
    assert [r["note"] for r in rows] == ["-", "Article with id PMC101 is not open access one", "no PMCID"], rows
    assert "supp/Table_S4.csv" in text and res[0]["genes"] == ["GENE1", "GENE2", "GENE3"] and res[0]["n"] == 3
    assert kept["papers"][0]["claims"][0]["genes"] == ["GENE1", "GENE2", "GENE3"]
    assert extra[0]["bytes"] > 0 and res[1]["genes"] == ["GENE1", "GENE5", "GENE7"] and res[1]["n"] == res[1]["stated"] == 3
    out.log("the dry run's checks pass")
    out.finish("setgenes", net, dict(fabricated=True, rows=rows, sets=res))
    return dict(rows=rows, sets=res)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="mode", required=True)
    f = sub.add_parser("fetch")
    f.add_argument("--claims", type=Path, required=True)
    f.add_argument("--out", type=Path, required=True)
    f.add_argument("--extra", type=Path, help="files named beforehand: rid<TAB>url")
    f.add_argument("--versions", type=Path, help="another version's PMCID per record (rid<TAB>pmcid): those only")
    s = sub.add_parser("show")
    s.add_argument("--files", type=Path, nargs="+", required=True)
    s.add_argument("--out", type=Path, required=True)
    e = sub.add_parser("extract")
    e.add_argument("--specs", type=Path, required=True)
    e.add_argument("--files", type=Path, nargs="+", required=True, help="the results directories, the first run's first")
    e.add_argument("--claims", type=Path, required=True)
    e.add_argument("--out", type=Path, required=True)
    e.add_argument("--claims-out", type=Path)
    e.add_argument("--texts", type=Path, help="S3's texts (a set the paper lists in its text)")
    d = sub.add_parser("dry")
    d.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.mode == "fetch":
        out = S.Out(a.out)
        out.log("# Level B, the sets' genes: Europe PMC's supplementary files of the kept papers with a set (set_genes.py)")
        net = S.Net(out.log)
        rows = fetch(net, out, json.loads(a.claims.read_text()), read_versions(a.versions) if a.versions else None)
        extra = fetch_extra(net, out, read_extra(a.extra)) if a.extra else []
        out.finish("setgenes", net, dict(fabricated=False, claims_file=str(a.claims), rows=rows, extra=extra))
        return 0
    if a.mode == "show":
        show(a.files if len(a.files) > 1 else a.files[0], a.out)
        return 0
    if a.mode == "extract":
        kept = json.loads(a.claims.read_text())
        res = extract(json.loads(a.specs.read_text()), a.files if len(a.files) > 1 else a.files[0], kept, a.texts)
        a.out.write_text(json.dumps(res, indent=1, ensure_ascii=False) + "\n")
        for r in res:
            print(f"{r['rid']} claim {r['claim']} ({r['set_name']}): {r['n']} genes, the paper states {r['stated']}"
                  + (f"; {r['note']}" if r.get("note") else ""))
        if a.claims_out:
            a.claims_out.write_text(json.dumps(fill(kept, res), indent=1, ensure_ascii=False) + "\n")
        return 0
    dry(a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
