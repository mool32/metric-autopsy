#!/usr/bin/env python3
"""Criterion (ii) of change 17''s acceptance failed: at an equal detection the verdicts differ from 1337ead's in
15 of 80 seeds for mi_3bin and 2 of 80 for norm_pearson (detection_equalization.log). This reads the run's JSON
(no new run) and sets each differing seed's numbers beside each other: the verdict and its cause, the effect,
the effect at equal depth alone, s_m and GATE 4's status, in both engines, and the same numbers over the seeds
whose verdicts agree.

    python validation/exploratory/detection_equalization_flips.py > validation/exploratory/detection_equalization_flips.log
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
EQUAL = ("eq0", "eq3", "eq5", "eq7")


def f(x, n=4):
    return "—" if x is None else f"{x:+.{n}g}"


def main() -> int:
    d = json.loads((HERE / "detection_equalization.json").read_text())
    base = {(r["scenario"], r["metric"], r["seed"]): r for r in d["base"]}
    print("# Criterion (ii) of change 17''s acceptance: the seeds whose verdicts differ, from detection_equalization.json\n")
    for metric in ("mi_3bin", "norm_pearson"):
        rows = [(r, base[(r["scenario"], r["metric"], r["seed"])]) for r in d["new"]
                if r["metric"] == metric and r["scenario"].split(":")[0] in EQUAL]
        diff = [(r, b) for r, b in rows if (r["cls"], r["cause"]) != (b["cls"], b["cause"])]
        same = [(r, b) for r, b in rows if (r["cls"], r["cause"]) == (b["cls"], b["cause"])]
        print(f"## {metric}: {len(diff)} of {len(rows)} seeds differ\n")
        print("| scenario | seed | 1337ead: verdict [cause] | effect | s_m | GATE 4 | new: verdict [cause] | effect | "
              "effect at equal depth alone | s_m | GATE 4 |")
        print("|---|---|---|---|---|---|---|---|---|---|---|")
        for r, b in diff:
            print(f"| {r['scenario'].split(':')[0]} | {r['seed']} | {b['cls']} [{b['cause']}] | {f(b['effect'])} | "
                  f"{f(b['s_m'])} | {b['g4']} | {r['cls']} [{r['cause']}] | {f(r['effect'])} | {f(r['equal_depth_effect'])} | "
                  f"{f(r['s_m'])} | {r['g4']} |")
        for label, group in (("differ", diff), ("agree", same)):
            ds = [abs(r["s_m"] - b["s_m"]) / b["s_m"] for r, b in group if r["s_m"] and b["s_m"]]
            de = [abs(r["effect"] - b["effect"]) for r, b in group if r["effect"] is not None and b["effect"] is not None]
            same_depth = sum(1 for r, b in group if r["equal_depth_effect"] is not None and b["effect"] is not None
                             and abs(r["equal_depth_effect"] - b["effect"]) <= 1e-12 * max(1.0, abs(b["effect"])))
            g0 = sum(1 for r, b in group if "metric_invalid_gate0" in (r["cause"], b["cause"]))
            print(f"\nThe {len(group)} seeds that {label}: s_m's relative change, mean {sum(ds) / len(ds):.3f} "
                  f"(max {max(ds):.3f}); |effect's change| mean {sum(de) / len(de):.4g} (max {max(de):.4g}); the effect at "
                  f"equal depth alone equal to 1337ead's effect in {same_depth}; a GATE 0 FAIL in either engine in {g0}"
                  if ds and de else f"\nThe {len(group)} seeds that {label}: —")
        flips = {}
        for r, b in diff:
            k = f"{b['cls']} [{b['cause']}] -> {r['cls']} [{r['cause']}]"
            flips[k] = flips.get(k, 0) + 1
        print("\nThe changes: " + "; ".join(f"{k}: {v}" for k, v in sorted(flips.items())) + "\n")
    # the log's line "equals 1337ead's effect in 281 of 280 runs" counted scenario v above and not below
    pairs = [(r, base[(r["scenario"], r["metric"], r["seed"])]) for r in d["new"]
             if (r["scenario"], r["metric"], r["seed"]) in base]
    same = sum(1 for r, b in pairs if r["equal_depth_effect"] is not None and b["effect"] is not None
               and abs(r["equal_depth_effect"] - b["effect"]) <= 1e-12 * max(1.0, abs(b["effect"])))
    print(f"The new engine's effect at equal depth alone equals 1337ead's effect in {same} of the {len(pairs)} runs both "
          "engines made (the log's line \"281 of 280\" counted the clean data's run in the first number only).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
