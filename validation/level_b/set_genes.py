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
* `extract`: each set's genes by its spec (`set_specs.json`: the archive's member, its sheet, the row of its header,
  the column of the genes, the rows that are the set's; or an extra file and the set's name in it, a GMT), their
  number against the number the paper states, and the kept claims with the sets' genes filled.
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


def fetch(net, out, kept: dict) -> list:
    rows = []
    for p in kept["papers"]:
        if not any(c["kind"] == "set" for c in p["claims"]):
            continue
        row = dict(rid=p["rid"], pmcid=p.get("pmcid") or "-", status="-", file="-", bytes=0, sha256="-", note="-")
        rows.append(row)
        if not p.get("pmcid"):
            row["note"] = "no PMCID"
            continue
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
                   file=out.raw(f"setgenes/{p['rid']}-supplementary.{'zip' if is_zip else 'xml'}", body),
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
def archives(files: Path) -> dict:
    """The fetched archives, by rid (the gzip of each, as the fetch stored it)."""
    out = {}
    for line in (files / "setgenes_index.tsv").read_text().splitlines()[1:]:
        r = line.split("\t")
        if len(r) >= 4 and r[3].endswith(".zip.gz"):
            out[r[0]] = files / r[3]
    return out


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


def show(files: Path, out: Path) -> str:
    L = ["# The supplements' members, and each table's sheets and first rows (set_genes.py show)", ""]
    for rid, path in sorted(archives(files).items()):
        z = _zip(path)
        L += [f"## {rid} ({path.name})", ""]
        for info in z.infolist():
            L.append(f"- `{info.filename}` ({info.file_size:,} bytes)")
        for info in z.infolist():
            if not info.filename.lower().endswith(TABLES):
                continue
            data = z.read(info.filename)
            sheets = sheet_names(data) if info.filename.lower().endswith(".xlsx") else [None]
            for sh in sheets:
                try:
                    rows = table_rows(info.filename, data, sh)
                except Exception as e:  # a table this reader cannot read is said so
                    L += ["", f"### `{info.filename}`{' / ' + sh if sh else ''}: not read ({type(e).__name__}: {e})"]
                    continue
                L += ["", f"### `{info.filename}`{' / ' + sh if sh else ''}: {len(rows)} rows", "", "```"]
                L += ["\t".join(r)[:300] for r in rows[:SHOW_ROWS]] + ["```"]
        L.append("")
    text = "\n".join(L) + "\n"
    out.write_text(text)
    return text


# --------------------------------------------------------------------------- #
# extract
# --------------------------------------------------------------------------- #
def extract_one(z: zipfile.ZipFile, spec: dict) -> list:
    rows = table_rows(spec["member"], z.read(spec["member"]), spec.get("sheet"))
    h = int(spec.get("header_row", 0))
    header, body = rows[h], rows[h + 1:]
    col = header.index(spec["column"]) if spec.get("column") is not None else 0
    filt = {header.index(k): {str(x) for x in v} for k, v in (spec.get("rows") or {}).items()}
    genes = []
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


def extract(specs: list, files: Path, kept: dict) -> list:
    arch = archives(files)
    stated = {(p["rid"], c["id"]): c for p in kept["papers"] for c in p["claims"] if c["kind"] == "set"}
    res = []
    for s in specs:
        c = stated.get((s["rid"], str(s["claim"])))
        if c is None:
            raise SystemExit(f"{s['rid']} claim {s['claim']}: not a kept set's claim")
        if s.get("extra"):
            with gzip.open(files / "raw" / "setgenes" / f"{s['rid']}-{s['extra']}.gz", "rb") as g:
                genes = gmt_set(g.read(), s["gmt_set"])
        else:
            genes = extract_one(_zip(arch[s["rid"]]), s)
        res.append(dict(rid=s["rid"], claim=str(s["claim"]), set_name=c.get("set_name"), stated=c.get("set_size"),
                        n=len(genes), genes=genes, spec=s))
    return res


def fill(kept: dict, res: list) -> dict:
    by = {(r["rid"], r["claim"]): r for r in res}
    for p in kept["papers"]:
        for c in p["claims"]:
            r = by.get((p["rid"], c["id"]))
            if r is not None:
                c["genes"] = r["genes"]
                c["genes_from"] = f"set_genes.json ({r['spec'].get('member') or r['spec'].get('extra')})"
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
    s = sub.add_parser("show")
    s.add_argument("--files", type=Path, required=True)
    s.add_argument("--out", type=Path, required=True)
    e = sub.add_parser("extract")
    e.add_argument("--specs", type=Path, required=True)
    e.add_argument("--files", type=Path, required=True)
    e.add_argument("--claims", type=Path, required=True)
    e.add_argument("--out", type=Path, required=True)
    e.add_argument("--claims-out", type=Path)
    d = sub.add_parser("dry")
    d.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.mode == "fetch":
        out = S.Out(a.out)
        out.log("# Level B, the sets' genes: Europe PMC's supplementary files of the kept papers with a set (set_genes.py)")
        net = S.Net(out.log)
        rows = fetch(net, out, json.loads(a.claims.read_text()))
        extra = fetch_extra(net, out, read_extra(a.extra)) if a.extra else []
        out.finish("setgenes", net, dict(fabricated=False, claims_file=str(a.claims), rows=rows, extra=extra))
        return 0
    if a.mode == "show":
        show(a.files, a.out)
        return 0
    if a.mode == "extract":
        kept = json.loads(a.claims.read_text())
        res = extract(json.loads(a.specs.read_text()), a.files, kept)
        a.out.write_text(json.dumps(res, indent=1, ensure_ascii=False) + "\n")
        for r in res:
            print(f"{r['rid']} claim {r['claim']} ({r['set_name']}): {r['n']} genes, the paper states {r['stated']}")
        if a.claims_out:
            a.claims_out.write_text(json.dumps(fill(kept, res), indent=1, ensure_ascii=False) + "\n")
        return 0
    dry(a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
