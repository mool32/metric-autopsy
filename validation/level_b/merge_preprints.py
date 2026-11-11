"""Amendment 1, 12.2 (`plan_s1.md`, section 2): Europe PMC's preprints merged with level B's corrected list.

A preprint is the list's paper where it shares the list record's Europe PMC id, PMID or DOI, or its title
once lower-cased and stripped to letters and digits (at least 20 characters; the protocol's section 5 and
`search.merge`); preprints that share one of these among themselves are one paper (the smallest Europe PMC
id kept). The others are new: numbered from R03507 in the order of their Europe PMC ids, with the fields the
screening reads. Nothing is fetched; both inputs are checked against their SHA256 first.

    python -I validation/level_b/merge_preprints.py --list <the corrected list's directory> \
        --preprints <results/level-b-preprints> --out <dir>
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
LIST_SUMS = HERE / "runs" / "remerge.SHA256SUMS"
FIRST_NEW = 3507  # the corrected list holds R00001-R03506


def norm_title(t) -> str:
    return re.sub(r"[^a-z0-9]", "", str(t or "").lower())


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _sums(path: Path) -> dict:
    out = {}
    for line in path.read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        out[name.lstrip("*").removeprefix("./")] = digest
    return out


def _read(path: Path) -> list:
    with gzip.open(path, "rt") as f:
        return [json.loads(line) for line in f if line.strip()]


def _ppr_number(epmc_id: str) -> int:
    m = re.search(r"(\d+)$", epmc_id or "")
    return int(m.group(1)) if m else 10 ** 12


def merge_preprints(listing: list, preprints: list, first_new: int = FIRST_NEW) -> tuple[list, list]:
    """(the new records, one row per fetched preprint: its id, `in list` / `duplicate` / `new`, and the rid)."""
    keys = {}
    for r in listing:
        for e in r.get("epmc_ids") or []:
            keys.setdefault(("epmc", e), r["rid"])
        for p in r.get("pmids") or ([r["pmid"]] if r.get("pmid") else []):
            keys.setdefault(("pmid", p), r["rid"])
        for d in r.get("dois") or ([r["doi"]] if r.get("doi") else []):
            keys.setdefault(("doi", d.lower()), r["rid"])
        t = norm_title(r.get("title"))
        if len(t) >= 20:
            keys.setdefault(("title", t), r["rid"])

    def keys_of(p):
        ks = [("epmc", p["epmc_id"])] if p.get("epmc_id") else []
        ks += [("pmid", p["pmid"])] if p.get("pmid") else []
        ks += [("doi", p["doi"].lower())] if p.get("doi") else []
        t = norm_title(p.get("title"))
        ks += [("title", t)] if len(t) >= 20 else []
        return ks

    rows, fresh = [], []
    for p in sorted(preprints, key=lambda p: (_ppr_number(p.get("epmc_id")), p.get("epmc_id") or "")):
        hit = next((keys[k] for k in keys_of(p) if k in keys), None)
        if hit is not None:
            rows.append(dict(epmc_id=p.get("epmc_id"), status="in list" if not hit.startswith("new:") else "duplicate",
                             rid=hit.removeprefix("new:")))
            continue
        rid = f"R{first_new + len(fresh):05d}"
        rec = dict(rid=rid, title=p.get("title") or "", abstract=p.get("abstract") or "", journal=p.get("journal") or "",
                   year=p.get("year"), pub_types=p.get("pub_types") or [], pmid=p.get("pmid") or "",
                   doi=(p.get("doi") or "").lower(), preprint=True, pmids=[p["pmid"]] if p.get("pmid") else [],
                   dois=[p["doi"].lower()] if p.get("doi") else [], epmc_ids=[p["epmc_id"]] if p.get("epmc_id") else [],
                   found_by=["preprints"], languages=p.get("languages") or [""])
        fresh.append(rec)
        for k in keys_of(p):
            keys.setdefault(k, f"new:{rid}")
        rows.append(dict(epmc_id=p.get("epmc_id"), status="new", rid=rid))
    return fresh, rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--list", required=True, type=Path, help="the corrected list's directory (records.jsonl.gz)")
    ap.add_argument("--preprints", required=True, type=Path, help="results/level-b-preprints (preprints.jsonl.gz)")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args(argv)
    want = _sums(LIST_SUMS)["records.jsonl.gz"]
    got = _sha256(a.list / "records.jsonl.gz")
    if got != want:
        raise SystemExit(f"the list's records.jsonl.gz has SHA256 {got}, not the corrected list's {want}")
    pre_sums = _sums(a.preprints / "SHA256SUMS")
    got = _sha256(a.preprints / "preprints.jsonl.gz")
    if pre_sums.get("preprints.jsonl.gz") != got:
        raise SystemExit(f"preprints.jsonl.gz has SHA256 {got}, not its branch's {pre_sums.get('preprints.jsonl.gz')}")
    listing, preprints = _read(a.list / "records.jsonl.gz"), _read(a.preprints / "preprints.jsonl.gz")
    fresh, rows = merge_preprints(listing, preprints)
    a.out.mkdir(parents=True, exist_ok=True)
    with gzip.GzipFile(filename="", mode="wb", fileobj=(a.out / "preprints_new.jsonl.gz").open("wb"), mtime=0) as g:
        for r in fresh:
            g.write((json.dumps(r, ensure_ascii=False, sort_keys=True) + "\n").encode())
    with (a.out / "preprints_merge.tsv").open("w") as f:
        f.write("epmc_id\tstatus\trid\n")
        for r in rows:
            f.write(f"{r['epmc_id']}\t{r['status']}\t{r['rid']}\n")
    counts = {s: sum(r["status"] == s for r in rows) for s in ("in list", "duplicate", "new")}
    langs = {}
    for r in fresh:
        lang = (r["languages"] or [""])[0] or "none"
        langs[lang] = langs.get(lang, 0) + 1
    lines = [
        "# Level B, amendment 1: the preprints merged with the corrected list (merge_preprints.py)",
        "",
        f"The corrected list: {len(listing)} records (records.jsonl.gz, SHA256 {want}).",
        f"The preprints: {len(preprints)} records (preprints.jsonl.gz, SHA256 {pre_sums['preprints.jsonl.gz']}).",
        f"Already in the list: {counts['in list']}; the same paper as an earlier preprint: {counts['duplicate']}; "
        f"new: {counts['new']}"
        + (f", numbered {fresh[0]['rid']} to {fresh[-1]['rid']}" if fresh else "") + ".",
        "The new preprints by language: " + (", ".join(f"{k} {v}" for k, v in sorted(langs.items())) or "-") + ".",
        f"preprints_new.jsonl.gz SHA256 {_sha256(a.out / 'preprints_new.jsonl.gz')}.",
    ]
    (a.out / "preprints_merge.log").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
