#!/usr/bin/env python3
"""Level B, S3 (`plan_s3.md`): each full-text session's counts, from the decision files themselves (not from what
the sessions said), for the record (`amendment1/ft/README.md`).

    python -I validation/level_b/ft_record.py --ft validation/level_b/amendment1/ft --out validation/level_b/amendment1/ft/sessions.tsv
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from screening import FT_CODES, FT_HEADER  # noqa: E402


def batch_counts(tsv: Path) -> dict:
    rows = [dict(zip(FT_HEADER, line.split("\t"))) for line in tsv.read_text().splitlines()[1:] if line.strip()]
    dec = [r for r in rows if r["decision"] != "CLAIM"]
    c = Counter(r["decision"] for r in dec)
    codes = Counter(r["code"] for r in dec if r["decision"] == "EXCLUDE")
    return dict(records=len(dec), KEEP=c.get("KEEP", 0), UNCLEAR=c.get("UNCLEAR", 0), EXCLUDE=c.get("EXCLUDE", 0),
                codes=", ".join(f"{k} {codes[k]}" for k in FT_CODES if codes.get(k)) or "-",
                claim_lines=sum(r["decision"] == "CLAIM" for r in rows))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ft", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    cols = ["pass", "batch", "records", "KEEP", "UNCLEAR", "EXCLUDE", "codes", "claim_lines"]
    lines = ["\t".join(cols)]
    for p in (1, 2):
        for tsv in sorted((a.ft / f"pass{p}").glob("batch_*.tsv")):
            r = dict(batch_counts(tsv), **{"pass": p, "batch": tsv.stem.split("_")[1]})
            lines.append("\t".join(str(r[c]) for c in cols))
    a.out.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
