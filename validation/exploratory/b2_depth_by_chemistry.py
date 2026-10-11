"""B2's depth by chemistry (10x 3' v2 or v3) and GATE 0's refusals on N7. Read after the results
(the owner's question of 2026-10-11: all of GATE 0's refusals fell on N7, B2's null, and B2 mixes
the two chemistries); it changes nothing in v1's result.

* Each donor id's assay comes from the Census (census_metadata.log, workflow run 38122733079).
* Depth is a cell's total count over the 2,000 genes the panel kept (the compact B2 the shards
  read), as the engine sees it; the Census's raw_sum over all genes is in census_metadata.log.
* Every N7 dataset is rebuilt from the key exactly as the shards built it (panel.build), and its
  card's outcome read from its report: the share of 10x 3' v2 cells in each group, and the ratio of
  the groups' median depths, for the refused cards and the others.
* GATE 0's flagged nuisance and its bias against the tolerance, from each refused report.

    python validation/exploratory/b2_depth_by_chemistry.py --results RESULTS --pilot pilot.json \
        --compact compact/B2.npz [--census validation/exploratory/census_metadata.log]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation" / "prereg"))
import numpy as np  # noqa: E402

import panel as P  # noqa: E402
import score  # noqa: E402

V2, V3 = "10x 3' v2", "10x 3' v3"


def census_assays(path: Path) -> dict:
    out = {}
    for line in path.read_text().splitlines():
        m = re.match(r"\s+donor (\S+) \| cells \d+ \| assay \{\"(10x 3' v[23])\": \d+\}", line)
        if m:
            out[m.group(1)] = m.group(2)
    return out


def describe(vals) -> str:
    v = np.asarray(vals, dtype=float)
    return f"median {np.median(v):.3f}, range {v.min():.3f} to {v.max():.3f}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--results", type=Path, required=True, help="results/panel-v1, its reports unpacked")
    ap.add_argument("--pilot", type=Path, required=True)
    ap.add_argument("--compact", type=Path, required=True, help="compact/B2.npz of results/panel-v1")
    ap.add_argument("--census", type=Path, default=ROOT / "validation" / "exploratory" / "census_metadata.log")
    args = ap.parse_args(argv)
    from scipy import stats
    assay = census_assays(args.census)
    bg = P.load_compact(args.compact)
    pilot_raw = args.pilot.read_bytes()
    pilot = json.loads(pilot_raw)
    manifest = json.loads((args.results / "manifest.json").read_text())
    key = manifest["key"]
    out = ["# B2's depth by chemistry and GATE 0's refusals on N7 (read after the results; v1 unchanged)",
           f"# compact B2 sha256 {hashlib.sha256(args.compact.read_bytes()).hexdigest()}; pilot.json sha256 "
           f"{hashlib.sha256(pilot_raw).hexdigest()}; key {key}; assays from {args.census.name}", ""]

    depth = np.asarray(bg.X.sum(axis=1), dtype=float)
    genes_detected = np.asarray((bg.X > 0).sum(axis=1), dtype=float)
    donors = list(bg.donors)
    missing = [d for d in donors if d not in assay]
    if missing:
        raise SystemExit(f"no assay for {missing}")
    out.append("## B2 by chemistry (depth: counts over the 2,000 kept genes)")
    per_donor = {d: float(np.median(depth[bg.donor == d])) for d in donors}
    for a in (V2, V3):
        ds = [d for d in donors if assay[d] == a]
        cells = np.isin(bg.donor, ds)
        out.append(f"  {a}: {len(ds)} donor ids, {int(cells.sum())} cells; per cell: median depth "
                   f"{np.median(depth[cells]):.0f}, median genes detected {np.median(genes_detected[cells]):.0f}; "
                   f"donor medians of depth {min(per_donor[d] for d in ds):.0f} to {max(per_donor[d] for d in ds):.0f}")
    studies = Counter((d.split("__")[1], assay[d]) for d in donors)
    out.append("  studies: " + ", ".join(f"{s} {n} ({a[-2:]})" for (s, a), n in sorted(studies.items())))
    v2d = [per_donor[d] for d in donors if assay[d] == V2]
    v3d = [per_donor[d] for d in donors if assay[d] == V3]
    mw = stats.mannwhitneyu(v3d, v2d, alternative="two-sided")
    out += [f"  v3 against v2, the donor ids' median depths: ratio of medians {np.median(v3d) / np.median(v2d):.2f}, "
            f"Mann-Whitney two-sided p {mw.pvalue:.4f}", ""]

    reports, _ = score.read_results(args.results, P.check_key(key))
    entries = [e for e in P.assign(key, pilot.get("dropped", ()), pilot.get("pool_size")) if e["condition"] == "N7"]
    rows = []
    for e in entries:
        X, obs, genes, cards = P.build(e, {"B2": bg}, pilot)
        tot = np.asarray(obs["total_counts"], dtype=float)
        v2 = np.array([assay[d] == V2 for d in obs["donor"]])
        grp = np.asarray(obs["group"])
        share = {g: float(v2[grp == g].mean()) for g in ("A", "B")}
        med = {g: float(np.median(tot[grp == g])) for g in ("A", "B")}
        rec = reports.get(cards[0]["id"])
        verdict, cause = (rec[0], rec[1]) if rec and not isinstance(rec, dict) else (None, None)
        outcome = P.outcome(verdict, cause)
        level = P.pool_level(P.conditions()["N7"], int(e["pair"]), pilot)
        msg = ""
        if outcome == P.REFUSAL:
            rep = json.loads((args.results / "reports" / f"{cards[0]['id']}.json").read_text())
            msg = next((g.get("message", "") for g in rep.get("gates", []) if g.get("gate") == 0), "")
        rows.append(dict(id=e["id"], pair=int(e["pair"]), genes=tuple(cards[0]["gene_pair"]), level=level,
                         refused=outcome == P.REFUSAL, v2_all=float(v2.mean()), v2_gap=abs(share["A"] - share["B"]),
                         log2_depth_ratio=abs(float(np.log2(med["B"] / med["A"]))), depth=float(np.median(tot)),
                         msg=msg))
    ref = [r for r in rows if r["refused"]]
    out.append(f"## GATE 0's refusals on N7: {len(ref)} of {len(rows)} cards "
               f"({', '.join(f'{lv} {sum(r['level'] == lv for r in ref)} of {sum(r['level'] == lv for r in rows)}' for lv in ('high', 'medium', 'low'))})")
    nuis = Counter(re.match(r"'(\w+)'", r["msg"]).group(1) if re.match(r"'(\w+)'", r["msg"]) else "other" for r in ref)
    out.append("  the flagged nuisance: " + ", ".join(f"{k} {v}" for k, v in nuis.items()))
    for lv in ("high", "medium"):
        sub = [r for r in rows if r["level"] == lv]
        a, b = [r for r in sub if r["refused"]], [r for r in sub if not r["refused"]]
        if not a:
            continue
        out.append(f"  {lv} level, refused ({len(a)}) against the others ({len(b)}):")
        for f, what in (("v2_all", "share of 10x 3' v2 cells in the dataset"),
                        ("v2_gap", "|difference of the groups' v2 shares|"),
                        ("log2_depth_ratio", "|log2 ratio of the groups' median depths|"),
                        ("depth", "the dataset's median depth")):
            p = stats.mannwhitneyu([r[f] for r in a], [r[f] for r in b], alternative="two-sided").pvalue
            out.append(f"    {what}: refused {describe([r[f] for r in a])}; others {describe([r[f] for r in b])}; "
                       f"Mann-Whitney p {p:.3f}")
    out.append("  refusals by gene pair (refused / N7 cards with the pair):")
    by_pair = Counter(r["pair"] for r in rows)
    for pair, n in sorted(Counter(r["pair"] for r in ref).items(), key=lambda kv: -kv[1]):
        g = next(r["genes"] for r in rows if r["pair"] == pair)
        lv = next(r["level"] for r in rows if r["pair"] == pair)
        out.append(f"    pair {pair:2d} {g[0]}-{g[1]} ({lv}): {n} / {by_pair[pair]}")
    biases = [re.search(r": ([\d.]+) \(95% lower bound ([\d.]+)\) against a tolerance of ([\d.]+)", r["msg"]) for r in ref]
    ratios = [float(m.group(1)) / float(m.group(3)) for m in biases if m]
    if ratios:
        out.append(f"  bias over its tolerance (0.5 x SESOI): {describe(ratios)} ({len(ratios)} refusals)")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
