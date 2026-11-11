#!/usr/bin/env python3
"""The independence of the screening's sessions, from their transcripts: for each session (its pass, its
batch, its transcript), how many times the transcript names the other pass's directory and each other
forbidden path, and which batches it names. Counts only; a session is "ok" where it names no forbidden
path and no batch but its own.

    python -I validation/level_b/independence.py --transcripts DIR --agents agents.tsv \
        --pass1 DIR1 --pass2 DIR2 [--forbid PATH ...] > independence.log

`agents.tsv`: one line per session, "pass<TAB>batch<TAB>session id"; the transcript is DIR/<id>.output.
The transcripts are the agent sessions' own records (not in the repository).

`--calls-only PATH ...`: paths the instructions themselves name (every session reads them, so its transcript
holds them): a session fails only where its own tool calls (the inputs of its requests, from the JSONL
transcript) name one. Added 2026-11-10 for S3, where the instructions cite the protocol (`JOURNAL.md`, J14).

`--text-pattern`: the full texts' file names (S3: `R<nnnnn>.txt`, the default; S4: `V<nnnnn>.txt`, its texts' keys,
the first column of its batch lists). Added 2026-11-11 for S4, before its passes.

`--exclude-attachment TYPE ...`: attachments of these types (JSONL lines `{"type": "attachment", "attachment":
{"type": TYPE, ...}}`) are counted apart: what they name is reported, not checked, and `--excluded-out` writes their
distinct contents for the record. Added 2026-11-11 for S4, after its passes (`JOURNAL.md`, J31): the harness put the
parent session's task list (`task_reminder`) in every S4 session's context; it names paths, no decision and no code.
The strict check (without the option) stays the first line of the record.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


def calls_text(s: str) -> str:
    """The session's own requests: the inputs of its tool calls, one per line, from a JSONL transcript."""
    out = []
    for line in s.splitlines():
        try:
            d = json.loads(line)
        except ValueError:
            continue
        m = d.get("message") if isinstance(d, dict) else None
        if not isinstance(m, dict) or m.get("role") != "assistant" or not isinstance(m.get("content"), list):
            continue
        for b in m["content"]:
            if isinstance(b, dict) and b.get("type") == "tool_use":
                out.append(json.dumps(b.get("input"), ensure_ascii=False))
    return "\n".join(out)


def split_attachments(s: str, types: set[str]) -> tuple[str, list[str]]:
    """The transcript without its attachments of `types`, and those attachments' contents (one JSON text each)."""
    keep, out = [], []
    for line in s.splitlines():
        try:
            d = json.loads(line)
        except ValueError:
            keep.append(line)
            continue
        a = d.get("attachment") if isinstance(d, dict) and d.get("type") == "attachment" else None
        if isinstance(a, dict) and a.get("type") in types:
            out.append(json.dumps(a.get("content"), ensure_ascii=False, sort_keys=True))
        else:
            keep.append(line)
    return "\n".join(keep), out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--transcripts", type=Path, required=True)
    ap.add_argument("--agents", type=Path, required=True)
    ap.add_argument("--pass1", required=True)
    ap.add_argument("--pass2", required=True)
    ap.add_argument("--forbid", nargs="*", default=[], help="paths no session may name (the other decisions)")
    ap.add_argument("--calls-only", nargs="*", default=[],
                    help="paths the instructions name: a session may not name them in its own tool calls")
    ap.add_argument("--lists", nargs=2, metavar=("DIR1", "DIR2"),
                    help="full texts: each pass's batch lists; a session may name no text file outside its own list")
    ap.add_argument("--batch-pattern", default=r"batches/batch_(\d+)",
                    help=r"the regular expression of a batch's path, its number the group (full texts: /batch_(\d+)\.list)")
    ap.add_argument("--text-pattern", default=r"\b(R\d{5})\.txt",
                    help=r"the regular expression of a full text's file name, its key the group (S4: \b(V\d{5})\.txt)")
    ap.add_argument("--exclude-attachment", nargs="*", default=[],
                    help="attachment types the harness put in the sessions' context: counted apart, not checked")
    ap.add_argument("--excluded-out", type=Path, help="the excluded attachments' distinct contents (JSON lines)")
    a = ap.parse_args(argv)
    dirs = {"1": a.pass1, "2": a.pass2}
    bad = n = 0
    distinct: dict[str, None] = {}
    for line in a.agents.read_text().splitlines():
        if not line.strip():
            continue
        p, b, sid = line.split("\t")
        n += 1
        f = a.transcripts / f"{sid.strip()}.output"
        if not f.exists():
            print(f"pass {p} batch {b}: no transcript")
            bad += 1
            continue
        s = f.read_text(errors="replace")
        excluded = ""
        if a.exclude_attachment:
            s, parts = split_attachments(s, set(a.exclude_attachment))
            distinct.update(dict.fromkeys(parts))
            excluded = "\n".join(parts)
        other = s.count(dirs["2" if p == "1" else "1"])
        forbidden = {x: s.count(x) for x in a.forbid}
        batches = sorted(set(re.findall(a.batch_pattern, s)))
        ok = other == 0 and not any(forbidden.values()) and batches == [b]
        extra = "".join(f"; {x} named {k} times" for x, k in forbidden.items())
        if a.calls_only:
            calls = calls_text(s)
            called = {x: calls.count(x) for x in a.calls_only}
            ok = ok and not any(called.values())
            extra += "".join(f"; {x} named {s.count(x)} times, in its own calls {k}" for x, k in called.items())
        if a.lists:
            own = {line.split("\t")[0] for line in (Path(a.lists[int(p) - 1]) / f"batch_{b}.list").read_text().splitlines()
                   if line.strip()}
            named = set(re.findall(a.text_pattern, s))
            ok = ok and named <= own
            extra += f"; texts named {len(named)}, outside its list {sorted(named - own)}"
        if a.exclude_attachment:
            names = [dirs["1"], dirs["2"]] + a.forbid + a.calls_only
            named = ", ".join(f"{x} {excluded.count(x)}" for x in names if excluded.count(x))
            extra += f"; excluded attachments {len(parts)} ({len(set(parts))} distinct), naming: {named or 'none'}"
        bad += not ok
        print(f"pass {p} batch {b}: other pass's directory named {other} times; own {s.count(dirs[p])}{extra}; "
              f"batches named {batches}; {'ok' if ok else 'CHECK'}")
    print(f"{n} sessions: " + ("all independent" if bad == 0 else f"{bad} to check"))
    if a.excluded_out:
        a.excluded_out.write_text("".join(x + "\n" for x in distinct))
        print(f"excluded attachments' distinct contents: {len(distinct)}, written to {a.excluded_out.name}")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
