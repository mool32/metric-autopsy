#!/usr/bin/env python3
"""What the authors of Tabula Muris Senis write about the origin of the 3-month FACS data, the batches and
the time points (the owner's second decisions of 2026-10-18, item 8; validation/exploratory/real_pilot_6_plan.md).

The run fetches each source of real_pilot_6_spec.TEXT_SOURCES once (the Nature paper's NIH author manuscript
from Europe PMC and from NCBI's BioC, and the bioRxiv preprint's JATS text), and for each records its URL,
its outcome, its size and the SHA256 of what was fetched, and quotes every sentence that matches one of
real_pilot_6_spec.TEXT_PATTERNS, with the section it stands in (at most MAX_QUOTES, each cut at
MAX_QUOTE_CHARS). The texts themselves are not written anywhere: the author manuscript has no Creative
Commons licence. The standard library only; it imports neither the engine nor numpy. What is fetched is
data: it is parsed with regular expressions, never executed or rendered.

    python real_pilot_6_text.py --out OUT
    python real_pilot_6_text.py --out OUT --synthetic   # fabricated responses, no network (the dry run)
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import platform
import re
import sys
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import real_pilot_6_spec as S  # noqa: E402

UA = "metric-autopsy validation (https://github.com/mool32/metric-autopsy)"
_BLOCK = re.compile(r"<(p|text|title|abstract)\b[^>]*>(.*?)</\1>", re.S | re.I)
_TAG = re.compile(r"<[^>]+>")
_SENT = re.compile(r"(?<=[.!?])\s+(?=[A-Z(\[0-9])")


def fetch(url: str, synthetic: bool, tries: int = 4) -> tuple[bytes | None, str]:
    """(the body, the outcome), never raising."""
    if synthetic:
        return fabricated(url), "ok (SYNTHETIC)"
    last = ""
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=120) as r:  # noqa: S310 (fixed https URLs of the spec)
                return r.read(), f"ok {r.status}"
        except Exception as exc:  # recorded, the next try after a pause
            last = f"{type(exc).__name__}: {exc}"
            time.sleep(5 * (i + 1))
    return None, f"FAIL {last}"


def fabricated(url: str) -> bytes:
    """The dry run's responses: a JATS text, a BioC text and bioRxiv's details, all SYNTHETIC."""
    if "api.biorxiv.org" in url:
        return json.dumps({"collection": [
            {"version": "1", "license": "cc_by", "jatsxml": "https://example.org/v1.source.xml"},
            {"version": "2", "license": "cc_by", "jatsxml": "https://example.org/v2.source.xml"}]}).encode()
    body = ("<article><front><title>SYNTHETIC</title><abstract><p>SYNTHETIC abstract.</p></abstract></front><body>"
            "<sec><title>Results</title><p>SYNTHETIC: the 3-month data were previously published as the Tabula Muris. "
            "A sentence about nothing. The time points were 1, 3, 18, 21, 24 and 30 months.</p></sec>"
            "<sec><title>Methods</title><p>SYNTHETIC: each batch was sorted on one day.</p></sec></body></article>")
    if "bionlp" in url:
        body = ("<collection><document><passage><infon key=\"section_type\">METHODS</infon><text>SYNTHETIC: the 3m "
                "mice came from an earlier batch.</text></passage></document></collection>")
    return body.encode()


def blocks(xml: str) -> list[tuple[str, str]]:
    """(section, paragraph text) for every <p>, BioC <text> and <abstract> block, the section being the last
    <title> before it (or BioC's section_type)."""
    out, section = [], ""
    for m in _BLOCK.finditer(xml):
        tag, inner = m.group(1).lower(), m.group(2)
        text = re.sub(r"\s+", " ", html.unescape(_TAG.sub(" ", inner))).strip()
        if tag == "title":
            section = text[:120]
            continue
        if tag == "text":
            st = re.findall(r'<infon key="section_type">([^<]*)</infon>', xml[max(0, m.start() - 400):m.start()])
            section = st[-1] if st else section
        if text:
            out.append((section, text))
    return out


def quotes(xml: str) -> tuple[list[dict], int]:
    pats = [re.compile(p, re.I) for p in S.TEXT_PATTERNS]
    seen, out, n_sent = set(), [], 0
    for section, text in blocks(xml):
        for sent in _SENT.split(text):
            n_sent += 1
            hits = [p.pattern for p in pats if p.search(sent)]
            key = sent.strip().lower()
            if hits and key not in seen:
                seen.add(key)
                out.append(dict(section=section, matched=hits, text=sent.strip()[:S.MAX_QUOTE_CHARS],
                                cut=len(sent.strip()) > S.MAX_QUOTE_CHARS))
    return out, n_sent


def check(name: str, what: str, url: str, synthetic: bool) -> dict:
    rec = dict(source=name, what=what, url=url)
    body, outcome = fetch(url, synthetic)
    rec.update(outcome=outcome)
    if body is None:
        return rec
    if name == "tms-biorxiv":  # the details first, then the newest version's JATS text
        try:
            versions = json.loads(body.decode("utf-8", "replace")).get("collection") or []
        except ValueError as exc:
            rec.update(outcome=f"FAIL the details are not JSON ({exc})")
            return rec
        if not versions:
            rec.update(outcome="FAIL bioRxiv lists no version")
            return rec
        last = versions[-1]
        rec.update(details_sha256=hashlib.sha256(body).hexdigest(), version=last.get("version"),
                   license=last.get("license"), jats_url=last.get("jatsxml"))
        if not last.get("jatsxml"):
            rec.update(outcome="FAIL the newest version has no JATS text")
            return rec
        body, outcome = fetch(last["jatsxml"], synthetic)
        rec.update(outcome=outcome)
        if body is None:
            return rec
    text = body.decode("utf-8", "replace")
    q, n_sent = quotes(text)
    rec.update(bytes=len(body), sha256=hashlib.sha256(body).hexdigest(), paragraphs=len(blocks(text)),
               sentences=n_sent, matched=len(q), quotes=q[:S.MAX_QUOTES], quotes_cut=max(0, len(q) - S.MAX_QUOTES))
    if not rec["paragraphs"]:
        rec["outcome"] += "; no paragraph found (not a full text?)"
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--synthetic", action="store_true")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    tag = "SYNTHETIC " if args.synthetic else ""
    recs = [check(n, w, u, args.synthetic) for n, w, u in S.TEXT_SOURCES]
    lines = [f"# {tag}What the TMS authors write (validation/exploratory/real_pilot_6_plan.md; real_pilot_6_text.py)",
             "", f"python {platform.python_version()}; the paper doi:{S.TMS_DOI} ({S.TMS_PMCID}), the preprint "
             f"doi:{S.TMS_PREPRINT_DOI}; the patterns {list(S.TEXT_PATTERNS)}; the texts are not archived", ""]
    for r in recs:
        lines.append(f"## {r['source']}: {r['what']}")
        lines.append(f"- {r['url']}: {r['outcome']}")
        for k in ("version", "license", "jats_url"):
            if r.get(k):
                lines.append(f"- {k}: {r[k]}")
        if "sha256" in r:
            lines.append(f"- {r['bytes']} bytes, sha256 {r['sha256']}; {r['paragraphs']} paragraphs, {r['sentences']} "
                         f"sentences, {r['matched']} matched" + (f" ({r['quotes_cut']} beyond the first "
                                                                  f"{S.MAX_QUOTES} not quoted)" if r["quotes_cut"] else ""))
            for q in r["quotes"]:
                lines.append(f"  - [{q['section'] or 'no section'}] \"{q['text']}\"" + (" [cut]" if q["cut"] else ""))
        lines.append("")
    ok = [r for r in recs if r.get("paragraphs")]
    lines.append(f"{tag}sources with a text: {len(ok)} of {len(recs)}")
    (args.out / "real_pilot_6_text.json").write_text(json.dumps(dict(synthetic=args.synthetic, sources=recs), indent=1))
    (args.out / "real_pilot_6_text.log").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
