"""Validation v1's pilot, as the cards used it, and the real effects' denominators. Read after the
results from pilot.json (fixed in the run tag before the key) and scores.json; it changes nothing.

* Per background and expression level: the SESOI and the saturation dose, which every card passed to
  the engine (the SESOI on the metric's observed scale, delta_min = 0.5 x SESOI; GATE 4's coupling
  at the saturation dose).
* B1's key dose per level (v1.md section 3.2: the smallest grid dose at which the oracle detects E1
  with power >= 0.9 and the level's mean Δ* is >= 1.25 x SESOI); where none is found, the real
  effects use the saturation dose. Δ* at E1's four doses (0.25, 0.5, 1 and 1.5 times that base).
* N2's raw power at c = 0.5 and the levels where it reaches 0.9 (`informative_levels`: there N2
  tests the engine's correction of depth).
* The real effects' cards: by condition, level and the metric's truth; the cards with a valid
  metric (S5), the establishable ones (S3) and those at or above the SESOI (S6).

    python validation/exploratory/v1_pilot_table.py --results RESULTS --pilot pilot.json --scores scores.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation" / "prereg"))
import panel as P  # noqa: E402
import score  # noqa: E402

LEVELS = ("high", "medium", "low")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--results", type=Path, required=True, help="results/panel-v1, its reports unpacked")
    ap.add_argument("--pilot", type=Path, required=True)
    ap.add_argument("--scores", type=Path, required=True)
    args = ap.parse_args(argv)
    pilot_raw, scores_raw = args.pilot.read_bytes(), args.scores.read_bytes()
    pilot, scores = json.loads(pilot_raw), json.loads(scores_raw)
    out = ["# Validation v1's pilot as the cards used it, and the real effects' denominators (read after the results)",
           f"# pilot.json sha256 {hashlib.sha256(pilot_raw).hexdigest()}; scores.json sha256 "
           f"{hashlib.sha256(scores_raw).hexdigest()}", "",
           "## Per background and level: the SESOI and GATE 4's dose (the saturation dose) every card passed",
           f"  {'':4s} {'level':7s} {'SESOI':>6s} {'found':>6s} {'delta_min':>10s} {'saturation dose':>16s}"]
    for bg in ("B1", "B2"):
        for lv in LEVELS:
            s = pilot["sesoi"][bg][lv]
            out.append(f"  {bg:4s} {lv:7s} {s:6.3f} {str(pilot['sesoi_found'][bg][lv]):>6s} {0.5 * s:10.4f} "
                       f"{pilot['saturation_dose'][bg][lv]:16.2f}")
    out += ["", "## B1's real effects: the key dose and Δ* (the level's mean over its pairs)"]
    for lv in LEVELS:
        s = pilot["sesoi"]["B1"][lv]
        scan = pilot["dose_scan"][lv]
        best_p = max(scan, key=lambda d: d["power"])
        best_d = max(scan, key=lambda d: d["delta"])
        meets = [d["dose"] for d in scan if d["power"] >= 0.9 and d["delta"] >= 1.25 * s]
        dl = pilot["delta_level"][lv]
        out.append(f"  {lv}: key dose found {pilot['key_dose_found'][lv]} (doses meeting power >= 0.9 and "
                   f"Δ* >= 1.25 x SESOI = {1.25 * s:.3f}: {meets or 'none'}); the scan's largest power "
                   f"{best_p['power']:.2f} (dose {best_p['dose']:g}), largest Δ* {best_d['delta']:.3f} (dose "
                   f"{best_d['dose']:g}); base dose of the real effects {pilot['e_dose'][lv]:g} (the saturation dose)")
        out.append("      Δ* at E1's doses (x the base): " + ", ".join(
            f"x{f} {dl[f]['value']:.3f}" for f in ("0.25", "0.5", "1", "1.5")) +
            f"; at the base dose {dl['1']['value']:.3f} against the SESOI {s:.3f} (ratio {dl['1']['value'] / s:.2f})")
    out += ["", "## N2's raw power at c = 0.5 (the oracle's power to detect N2's raw difference)",
            "  " + ", ".join(f"{lv} {pilot['n2_raw_power'][lv]:.2f}" for lv in LEVELS) +
            f"; levels where it reaches 0.9 (`informative_levels`): {pilot['informative_levels'] or 'none'}",
            "  so at no level did N2 test the engine's correction of depth", ""]

    manifest = json.loads((args.results / "manifest.json").read_text())
    reports, _ = score.read_results(args.results, P.check_key(manifest["key"]))
    entries = P.assign(manifest["key"], pilot.get("dropped", ()), pilot.get("pool_size"))
    rows = [r for r in score.card_rows(entries, reports, pilot) if r["condition"] in ("E1", "E2", "E3")]
    out.append(f"## The real effects' cards: {len(rows)}")
    for c in ("E1", "E2", "E3"):
        rs = [r for r in rows if r["condition"] == c]
        out.append(f"  {c}: {len(rs)} cards; " + "; ".join(
            f"{lv} " + "/".join(str(sum(r['level'] == lv and r['truth'] == t for r in rs)) for t in ("valid", "ambiguous", "blind"))
            for lv in LEVELS) + " (valid/ambiguous/blind)")
    tr = Counter(r["truth"] for r in rows)
    cells = scores["criteria"]
    out += [f"  by truth: valid {tr['valid']}, ambiguous {tr['ambiguous']}, blind {tr['blind']}",
            f"  S5's effect cell (a valid metric): {cells['S5']['cells']['S5:effect']['n']}; "
            f"S3's effect stratum (establishable): {cells['S3']['strata']['effect']['n']}; "
            f"S6's effect cell (at or above the SESOI, where NO DETECTABLE EFFECT is an error): "
            f"{cells['S6']['cells']['S6:effect']['n']}; S2's real-effect cell (SUPPORTED an error): "
            f"{cells['S2']['cells']['S2:E']['n']}"]
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
