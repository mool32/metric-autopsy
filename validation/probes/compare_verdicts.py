#!/usr/bin/env python3
"""The dev set's verdicts at two engines, line by line: two logs of `verdicts_v03.py` (the same script, so the
same cases in the same order), and every line that differs, under its probe and case.

    python validation/probes/compare_verdicts.py BEFORE.log AFTER.log
"""
from __future__ import annotations

import difflib
import sys
from pathlib import Path


def main(before: Path, after: Path) -> int:
    a, b = before.read_text().splitlines(), after.read_text().splitlines()
    print(f"# The dev set's verdicts: {before.name} ({a[0][2:]}) against {after.name} ({b[0][2:]})\n")
    a, b = a[1:], b[1:]
    changed = 0
    section, case = "", ""
    sm = difflib.SequenceMatcher(a=a, b=b, autojunk=False)
    context = {}
    for i, line in enumerate(a):  # the probe and case each line of the earlier log sits under
        if line.startswith("==="):
            section, case = line, ""
        elif line.startswith("  ") and not line.startswith("    "):
            case = line.strip()
        context[i] = (section, case)
    shown = set()
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal":
            continue
        sec, cas = context.get(i1, context.get(max(i1 - 1, 0), ("", "")))
        if (sec, cas) not in shown:
            print(f"{sec}\n  [{cas}]" if cas else sec)
            shown.add((sec, cas))
        for line in a[i1:i2]:
            print(f"  - {line.strip()}")
        for line in b[j1:j2]:
            print(f"  + {line.strip()}")
        changed += max(i2 - i1, j2 - j1)
        if any("VERDICT:" in line for line in a[i1:i2] + b[j1:j2]):
            print("    (a verdict changed)")
    n_verdicts = sum("VERDICT:" in line for line in a)
    n_changed = sum(1 for tag, i1, i2, j1, j2 in sm.get_opcodes() if tag != "equal"
                    for line in a[i1:i2] if "VERDICT:" in line)
    print(f"\n{changed} lines differ; VERDICT lines: {n_verdicts} before, {n_changed} of them changed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(Path(sys.argv[1]), Path(sys.argv[2])))
