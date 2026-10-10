"""Validation v1's false SUPPORTED on the pure nulls (N1, N5, N7) and on N8, by the metric's truth
and level. Exploratory, not a criterion, decided after the results.

v1_by_level.py reads only scores.json, which tabulates the outcomes by condition and level but not
by condition and truth. This script rebuilds every card's row as score.py does (`score.card_rows`:
the card assignment from the key, the oracle's truth from pilot.json, the outcome from the
report), so a rate can be read on the cards whose metric is valid. It changes nothing in v1's
result.

    python validation/exploratory/v1_nulls_by_truth.py --results RESULTS --pilot pilot.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation" / "prereg"))
sys.path.insert(0, str(ROOT / "validation" / "exploratory"))
import panel as P  # noqa: E402
import score  # noqa: E402
from v1_by_level import fmt  # noqa: E402  (k/n = rate (95% Clopper-Pearson CI), rounded half up)

PURE = ("N1", "N5", "N7")
LEVELS = ("high", "medium", "low")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--results", type=Path, required=True, help="results/panel-v1, its reports unpacked")
    ap.add_argument("--pilot", type=Path, required=True)
    args = ap.parse_args(argv)
    manifest = json.loads((args.results / "manifest.json").read_text())
    pilot_raw = args.pilot.read_bytes()
    pilot = json.loads(pilot_raw)
    key = manifest["key"]
    reports, _ = score.read_results(args.results, P.check_key(key))
    entries = P.assign(key, pilot.get("dropped", ()), pilot.get("pool_size"))
    cond_of = {e["id"]: e["condition"] for e in entries}
    rows = score.card_rows(entries, reports, pilot)
    for r in rows:
        r["condition"] = cond_of[r["dataset"]]
    out = ["# Validation v1: false SUPPORTED on the pure nulls and on N8, by the metric's truth and level",
           "# Exploratory, not a criterion; rows rebuilt as score.py builds them (score.card_rows)",
           f"# key {key}; pilot.json sha256 {hashlib.sha256(pilot_raw).hexdigest()}; {len(rows)} cards", ""]
    for name, conds in (("the pure nulls N1, N5 and N7", PURE), ("N8", ("N8",))):
        sel = [r for r in rows if r["condition"] in conds]
        out.append(f"{name}: {len(sel)} cards")
        for truth in ("valid", "ambiguous", "blind"):
            for level in LEVELS:
                cell = [r for r in sel if r["truth"] == truth and r["level"] == level]
                if cell:
                    k = sum(r["outcome"] == P.SUPPORTED for r in cell)
                    out.append(f"  {truth:9s} {level:6s}: false SUPPORTED {fmt(k, len(cell))}")
        valid = [r for r in sel if r["truth"] == "valid"]
        k = sum(r["outcome"] == P.SUPPORTED for r in valid)
        out.append(f"  every card with a valid metric: false SUPPORTED {fmt(k, len(valid))}")
        out.append("")
    useless = [r for r in rows if r["truth"] in ("useless", "constant")]
    k = sum(r["outcome"] == P.SUPPORTED for r in useless)
    out.append(f"useless metrics (N6a-c, the constant included): SUPPORTED on {k} of {len(useless)} cards")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
