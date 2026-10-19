#!/usr/bin/env python3
"""The independence of the screening's sessions, from their transcripts: for each session (its pass, its
batch, its transcript), how many times the transcript names the other pass's directory and each other
forbidden path, and which batches it names. Counts only; a session is "ok" where it names no forbidden
path and no batch but its own.

    python -I validation/level_b/independence.py --transcripts DIR --agents agents.tsv \
        --pass1 DIR1 --pass2 DIR2 [--forbid PATH ...] > independence.log

`agents.tsv`: one line per session, "pass<TAB>batch<TAB>session id"; the transcript is DIR/<id>.output.
The transcripts are the agent sessions' own records (not in the repository).
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--transcripts", type=Path, required=True)
    ap.add_argument("--agents", type=Path, required=True)
    ap.add_argument("--pass1", required=True)
    ap.add_argument("--pass2", required=True)
    ap.add_argument("--forbid", nargs="*", default=[], help="paths no session may name (the other decisions)")
    ap.add_argument("--batch-pattern", default=r"batches/batch_(\d+)",
                    help=r"the regular expression of a batch's path, its number the group (full texts: /batch_(\d+)\.list)")
    a = ap.parse_args(argv)
    dirs = {"1": a.pass1, "2": a.pass2}
    bad = n = 0
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
        other = s.count(dirs["2" if p == "1" else "1"])
        forbidden = {x: s.count(x) for x in a.forbid}
        batches = sorted(set(re.findall(a.batch_pattern, s)))
        ok = other == 0 and not any(forbidden.values()) and batches == [b]
        bad += not ok
        extra = "".join(f"; {x} named {k} times" for x, k in forbidden.items())
        print(f"pass {p} batch {b}: other pass's directory named {other} times; own {s.count(dirs[p])}{extra}; "
              f"batches named {batches}; {'ok' if ok else 'CHECK'}")
    print(f"{n} sessions: " + ("all independent" if bad == 0 else f"{bad} to check"))
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
