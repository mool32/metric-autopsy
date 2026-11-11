#!/usr/bin/env python3
"""How many of validation v1's donor-level values of norm_pearson were its fallback 0.0 (the owner's
decision of 2026-10-18, item 4; exploratory: v1's datasets are development data now, v1.md section 6).

v0.3.0's `norm_pearson` returns 0.0 where fewer than 3 cells carry both genes, or where one gene does
not vary among them: a value that was not measured, reported as no correlation. Nothing here changes
v1's result or the engine.

`run` rebuilds every dataset of v1 from its key (panel.build), checks it against the manifest's
data_sha256, and for every claim card on norm_pearson counts, per donor (in N5's paired design per
donor and group), the cells that carry both genes of the card's pair and whether the fallback applies:
- on the uncorrected cells (the raw donor values);
- where v1's engine estimated the card's effect (its report holds the donor values), on each of the
  five thinnings to equal depth that the effect averaged, redrawn with the tag's engine and the card's
  seed exactly as `effect.estimate_effect` drew them; each donor's mean over its thinnings is checked
  against the corrected value in v1's report, and its raw value against the report's raw value.
It also counts the donor values below 20 cells with both genes, the per-replicate minimum of v0.4.

    python validation/exploratory/v1_fallbacks.py run --results RESULTS --compact COMPACT \
        --pilot pilot.json --engine-src TAG_SRC --out fallbacks.jsonl [--workers 4] [--datasets N]
    python validation/exploratory/v1_fallbacks.py report --results RESULTS --pilot pilot.json \
        --rows fallbacks.jsonl > v1_fallbacks.log

RESULTS is results/panel-v1 with its reports unpacked (manifest.json, reports/); COMPACT the run's
compact backgrounds (compact.tar unpacked); pilot.json the run tag's; TAG_SRC the tag's src/
(`git archive v0.3.0-prereg src`), put first on the path so that the thinning is the tag's.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation" / "prereg"))
sys.path.insert(0, str(ROOT / "validation" / "exploratory"))
import frozen  # noqa: E402  (standard library only)

frozen.pin_numerics()  # one BLAS thread per process, set before numpy loads, as the run did

import numpy as np  # noqa: E402

MIN_CODETECTED = 20  # v0.4's per-replicate minimum of cells with both genes (metrics.MIN_CODETECTED)
N_EQUALIZE = 5       # the tag's estimate_effect default: the corrected value averages five thinnings


def fallback_kind(a: np.ndarray, b: np.ndarray, both: np.ndarray) -> str | None:
    """The tag's norm_pearson on these cells: 'fewer than 3' or 'no variance' where it returns its
    fallback 0.0, None where it returns a correlation."""
    if both.sum() < 3:
        return "fewer than 3"
    if a[both].std() == 0 or b[both].std() == 0:
        return "no variance"
    return None


def _pair_arrays(X: np.ndarray, ia: int, ib: int):
    tot = X.sum(axis=1)
    tot = np.where(tot == 0, 1.0, tot)
    a_raw, b_raw = X[:, ia], X[:, ib]
    return (a_raw > 0) & (b_raw > 0), np.log1p(a_raw / tot * 1e4), np.log1p(b_raw / tot * 1e4)


def _pearson(a, b, both) -> float:
    return 0.0 if fallback_kind(a, b, both) else float(np.corrcoef(a[both], b[both])[0, 1])


def unit_masks(group: np.ndarray, donor: np.ndarray) -> tuple[str, list]:
    """The donor-level units: nested (each donor in one group) or paired (N5: every donor in both
    groups, one unit per donor and group), as (label, donor, group, mask)."""
    donors = list(dict.fromkeys(donor))
    member = {d: set(group[donor == d]) for d in donors}
    if all(len(m) == 2 for m in member.values()):
        return "paired", [(d, g, (donor == d) & (group == g)) for d in donors for g in ("A", "B")]
    if all(len(m) == 1 for m in member.values()):
        return "nested", [(d, next(iter(member[d])), donor == d) for d in donors]
    return "mixed", []


def card_units(X, genes, card, group, donor) -> tuple[str, list]:
    ia, ib = genes.index(card["gene_pair"][0]), genes.index(card["gene_pair"][1])
    both, a, b = _pair_arrays(X, ia, ib)
    design, units = unit_masks(group, donor)
    out = []
    for d, g, m in units:
        mb = m & both
        out.append(dict(donor=d, group=g, cells=int(m.sum()), both=int(mb.sum()),
                        fallback=fallback_kind(a, b, mb), value=_pearson(a, b, mb)))
    return design, out


def thinned_units(X, obs, genes, card, logged_values: list) -> dict:
    """The five thinnings the tag's estimate_effect averaged, redrawn with its own code and seed."""
    from metric_autopsy import effect as E
    from metric_autopsy.core import SimpleData
    from metric_autopsy.equalize import thin_to_match
    from metric_autopsy.report import normalize_prereg

    import run_panel as R
    _, kw = R.make_run_args(card)
    prereg = normalize_prereg(card["prereg"])
    groups, rep_col = tuple(card["groups"]), card["replicate_col"]
    data = SimpleData(X, obs, genes)
    sd = E._restrict(data, card["group_col"], groups)
    design = E.assess_design(sd, group_col=card["group_col"], groups=groups, within=[], replicate_col=rep_col,
                             estimand=prereg.get("estimand"), min_replicates=int(prereg["min_replicates"]),
                             spikein_prefix=prereg.get("spikein_prefix", "ERCC-"), _sd=sd)
    if not (design["identifiable"] and design["correction"] in E.THINNING):
        return dict(thinned=False, reason=design.get("reason") or design["correction"])
    rng = np.random.default_rng(kw["seed"])
    ia, ib = genes.index(card["gene_pair"][0]), genes.index(card["gene_pair"][1])
    group = np.asarray(sd.obs[card["group_col"]]).astype(str)
    donor = np.asarray(sd.obs[rep_col]).astype(str)
    _, units = unit_masks(group, donor)
    sets = []
    for _ in range(N_EQUALIZE):
        eq, _info = thin_to_match(sd, card["group_col"], groups, totals=design["totals"], within=[], rng=rng)
        both, a, b = _pair_arrays(np.asarray(eq.X, float), ia, ib)
        sets.append([(int((m & both).sum()), fallback_kind(a, b, m & both), _pearson(a, b, m & both))
                     for _, _, m in units])
    logged = {r["replicate"]: r for r in logged_values}
    out, match = [], True
    for j, (d, g, _) in enumerate(units):
        per = [s[j] for s in sets]
        if design["rep_design"]["kind"] == "paired":
            # a paired donor's value is its first group's minus its second's, per thinning
            out.append(dict(donor=d, group=g, both=[p[0] for p in per], fallbacks=sum(p[1] is not None for p in per)))
            continue
        mean = float(np.mean([p[2] for p in per]))
        lv = logged.get(d, {}).get("corrected")
        ok = lv is not None and abs(mean - float(lv)) <= 1e-9 * max(1.0, abs(mean))
        match &= ok
        out.append(dict(donor=d, group=g, both=[p[0] for p in per], fallbacks=sum(p[1] is not None for p in per),
                        corrected=mean, logged=lv, match=bool(ok)))
    if design["rep_design"]["kind"] == "paired":
        by_donor = defaultdict(list)
        for j, (d, g, _) in enumerate(units):
            by_donor[d].append(j)
        for d, (ja, jb) in by_donor.items():
            first = groups[0]
            sign = 1.0 if units[ja][1] == first else -1.0
            mean = float(np.mean([sign * (s[ja][2] - s[jb][2]) for s in sets]))
            lv = logged.get(d, {}).get("corrected")
            ok = lv is not None and abs(mean - float(lv)) <= 1e-9 * max(1.0, abs(mean))
            match &= ok
            for j in (ja, jb):
                out[j].update(corrected=mean, logged=lv, match=bool(ok))
    return dict(thinned=True, units=out, match=bool(match))


STATE: dict = {}


def analyze(entry: dict) -> list[dict]:
    import panel as P
    X, obs, genes, cards = P.build(entry, STATE["bgs"], STATE["pilot"])
    ok = P.dataset_sha256(X, obs, genes) == STATE["want"][entry["id"]]["data_sha256"]
    group, donor = np.asarray(obs["group"]).astype(str), np.asarray(obs["donor"]).astype(str)
    out = []
    for card in cards:
        if card["metric"] != "norm_pearson":
            continue
        rec = dict(card=card["id"], dataset=entry["id"], data_ok=bool(ok), pair=list(card["gene_pair"]))
        if ok:
            design, units = card_units(np.asarray(X, float), genes, card, group, donor)
            rec.update(design=design, raw=units)
            values = STATE["logged"].get(card["id"])
            if values:
                logged = {r["replicate"]: r for r in values}
                if design == "nested":
                    rec["raw_match"] = all(
                        logged.get(u["donor"]) is not None
                        and abs(float(logged[u["donor"]]["raw"]) - u["value"]) <= 1e-9 * max(1.0, abs(u["value"]))
                        for u in units)
                rec["eq"] = thinned_units(X, obs, genes, card, values)
        out.append(rec)
    return out


def run(args) -> int:
    import multiprocessing as mp

    sys.path.insert(0, str(args.engine_src.resolve()))  # the tag's engine, first on the path
    import blind
    import metric_autopsy
    import panel as P
    src = Path(metric_autopsy.__file__).resolve()
    if not str(src).startswith(str(args.engine_src.resolve())):
        raise SystemExit(f"the engine is not the tag's: {src}")
    manifest = json.loads((args.results / "manifest.json").read_text())
    pilot = json.loads(args.pilot.read_text())
    if manifest["pilot_sha256"] != blind.pilot_sha256(pilot):
        raise SystemExit("the pilot differs from the run's")
    bgs, info = blind.load_prepared(args.compact, None)
    if json.dumps(info, sort_keys=True) != json.dumps(manifest["backgrounds"], sort_keys=True):
        raise SystemExit("the compact backgrounds differ from the run's")
    key = P.check_key(manifest["key"])
    # the donor values each report holds (the cards whose effect the engine estimated)
    logged = {}
    for d in manifest["datasets"]:
        for c in d["cards"]:
            path = args.results / "reports" / f"{c['id']}.json"
            if path.exists():
                rep = json.loads(path.read_text())
                vals = (((rep.get("fields") or {}).get("effect") or {}).get("detail") or {}).get("replicate_values")
                if vals and rep.get("metric") == "norm_pearson":
                    logged[c["id"]] = vals
    want = {d["id"]: d for d in manifest["datasets"]}
    entries = [e for e in P.assign(key, pilot.get("dropped", ()), pilot.get("pool_size")) if e["id"] in want]
    if args.datasets:
        entries = entries[:args.datasets]
    done = set()
    if args.out.exists():
        done = {json.loads(line)["dataset"] for line in args.out.read_text().splitlines() if line.strip()}
    todo = [e for e in entries if e["id"] not in done]
    STATE.update(bgs=bgs, pilot=pilot, want=want, logged=logged)
    print(f"engine {metric_autopsy.__version__} from {src.parent}; {len(logged)} cards hold donor values; "
          f"{len(entries)} datasets, {len(done)} done "
          f"before, {len(todo)} to run on {args.workers} workers", file=sys.stderr)
    ctx = mp.get_context("fork")
    with open(args.out, "a") as fh, ctx.Pool(args.workers) as pool:
        for i, recs in enumerate(pool.imap_unordered(analyze, todo, chunksize=4), 1):
            for r in recs:
                fh.write(json.dumps(r) + "\n")
            if i % 250 == 0:
                fh.flush()
                print(f"{i}/{len(todo)}", file=sys.stderr)
    return 0


# --------------------------------------------------------------------------- #
# report
# --------------------------------------------------------------------------- #
def pct(k: int, n: int) -> str:
    return f"{k}/{n} ({100.0 * k / n:.2f}%)" if n else f"{k}/0"


def report(args) -> int:
    import panel as P
    import score
    manifest = json.loads((args.results / "manifest.json").read_text())
    pilot_raw = args.pilot.read_bytes()
    pilot = json.loads(pilot_raw)
    key = P.check_key(manifest["key"])
    reports, _ = score.read_results(args.results, key)
    entries = P.assign(key, pilot.get("dropped", ()), pilot.get("pool_size"))
    published = {r["id"]: r for r in score.card_rows(entries, reports, pilot)}
    by_entry = {e["id"]: e for e in entries}
    raw = args.rows.read_bytes()
    if args.rows.suffix == ".gz":
        raw = gzip.decompress(raw)
    rows = [json.loads(line) for line in raw.decode().splitlines() if line.strip()]
    levels = ("high", "medium", "low")
    out = ["# norm_pearson's fallback 0.0 in validation v1's donor-level values (v1_fallbacks.py; exploratory)",
           f"# key {key}; pilot.json sha256 {hashlib.sha256(pilot_raw).hexdigest()}; rows sha256 "
           f"{hashlib.sha256(raw).hexdigest()}",
           f"# claim cards on norm_pearson: {len(rows)}; datasets rebuilt identical to the manifest's: "
           f"{sum(r['data_ok'] for r in rows)} cards' datasets",
           "", "The fallback: fewer than 3 cells with both genes, or no variance among them (v0.3.0's metrics.norm_pearson)."]
    ok = [r for r in rows if r["data_ok"]]
    for r in ok:
        pub = published[r["card"]]
        r["level"], r["outcome"] = pub["level"], pub["outcome"]
        r["condition"] = by_entry[r["dataset"]]["condition"]

    def tally(sel, field="raw"):
        n = fb3 = fbv = below = 0
        for r in sel:
            for u in r[field]:
                n += 1
                fb3 += u["fallback"] == "fewer than 3"
                fbv += u["fallback"] == "no variance"
                below += u["both"] < MIN_CODETECTED
        return n, fb3, fbv, below

    out += ["", "## 1. The uncorrected donor values, every claim card on norm_pearson (rebuilt from the key)", "",
            "| level | cards | donor values | fallback: fewer than 3 | fallback: no variance | below 20 cells with both genes |",
            "|---|---|---|---|---|---|"]
    for lv in levels + ("all",):
        sel = [r for r in ok if lv == "all" or r["level"] == lv]
        n, fb3, fbv, below = tally(sel)
        out.append(f"| {lv} | {len(sel)} | {n} | {pct(fb3, n)} | {pct(fbv, n)} | {pct(below, n)} |")
    cards_with = Counter(r["level"] for r in ok if any(u["fallback"] for u in r["raw"]))
    out.append("\nCards with at least one uncorrected donor value at the fallback: "
               + ", ".join(f"{lv} {pct(cards_with[lv], sum(r['level'] == lv for r in ok))}" for lv in levels))

    est = [r for r in ok if r.get("eq", {}).get("thinned")]
    out += ["", "## 2. The cards whose effect v1's engine estimated (the donor values its verdicts used)", "",
            f"Cards: {len(est)}; the raw values recomputed here equal the reports' in "
            f"{sum(bool(r.get('raw_match')) for r in est if r['design'] == 'nested')} of "
            f"{sum(r['design'] == 'nested' for r in est)} nested cards; the corrected values (the mean over five "
            f"thinnings) equal the reports' in {sum(r['eq']['match'] for r in est)} of {len(est)} cards.", "",
            "| level | cards | donor values | raw at the fallback | corrected with >= 1 thinning at the fallback | "
            "corrected with all 5 at the fallback | cards with >= 1 corrected donor value touched |",
            "|---|---|---|---|---|---|---|"]
    for lv in levels + ("all",):
        sel = [r for r in est if lv == "all" or r["level"] == lv]
        n = sum(len(r["raw"]) for r in sel)
        raw_fb = sum(u["fallback"] is not None for r in sel for u in r["raw"])
        touched = sum(u["fallbacks"] >= 1 for r in sel for u in r["eq"]["units"])
        all5 = sum(u["fallbacks"] == N_EQUALIZE for r in sel for u in r["eq"]["units"])
        cards_t = sum(any(u["fallbacks"] >= 1 for u in r["eq"]["units"]) for r in sel)
        out.append(f"| {lv} | {len(sel)} | {n} | {pct(raw_fb, n)} | {pct(touched, n)} | {pct(all5, n)} | "
                   f"{pct(cards_t, len(sel))} |")
    out += ["", "By v1's outcome (the cards with an estimated effect): cards with at least one corrected donor value "
            "touched by the fallback", "", "| outcome | cards | touched |", "|---|---|---|"]
    for o in sorted({r["outcome"] for r in est}):
        sel = [r for r in est if r["outcome"] == o]
        t = sum(any(u["fallbacks"] >= 1 for u in r["eq"]["units"]) for r in sel)
        out.append(f"| {o} | {len(sel)} | {pct(t, len(sel))} |")
    out += ["", f"## 3. v0.4's minimum on the same cards: a donor value counts with at least {MIN_CODETECTED} cells with "
            "both genes", "",
            "v0.4's engine leaves out a replicate whose uncorrected value, or whose value at equal depth on every one of "
            "the five thinnings, has fewer cells with both genes, and makes the comparison UNIDENTIFIABLE where fewer "
            "than v1's minimum of 3 replicates per group (or pairs) are left; v1's draws, not v0.4's.", "",
            "| level | cards | donors left out | cards with >= 1 left out | cards left with fewer than 3 |",
            "|---|---|---|---|---|"]
    for lv in levels + ("all",):
        sel = [r for r in est if lv == "all" or r["level"] == lv]
        n_units = n_out = cards_out = cards_short = 0
        for r in sel:
            eq = {(u["donor"], u.get("group")): u for u in r["eq"]["units"]}
            short = {}
            for u in r["raw"]:
                e = eq.get((u["donor"], u.get("group"))) or eq.get((u["donor"], None))
                thin = e["both"] if e else []
                out_ = u["both"] < MIN_CODETECTED or (bool(thin) and all(b < MIN_CODETECTED for b in thin))
                short[(u["donor"], u.get("group"))] = out_
            n_units += len(short)
            n_out += sum(short.values())
            cards_out += any(short.values())
            if r["design"] == "paired":
                donors = {d for d, _ in short}
                left = sum(1 for d in donors if not any(v for (dd, _), v in short.items() if dd == d))
            else:
                left = min(sum(1 for (d, g), v in short.items() if g == grp and not v) for grp in ("A", "B"))
            cards_short += left < 3
        out.append(f"| {lv} | {len(sel)} | {pct(n_out, n_units)} | {pct(cards_out, len(sel))} | "
                   f"{pct(cards_short, len(sel))} |")
    out += ["", "By v1's outcome (the cards with an estimated effect): cards left with fewer than 3 replicates per "
            "group by v0.4's minimum", "", "| outcome | cards | fewer than 3 |", "|---|---|---|"]
    for o in sorted({r["outcome"] for r in est}):
        sel = [r for r in est if r["outcome"] == o]
        k = 0
        for r in sel:
            eq = {(u["donor"], u.get("group")): u for u in r["eq"]["units"]}
            short = {}
            for u in r["raw"]:
                e = eq.get((u["donor"], u.get("group"))) or eq.get((u["donor"], None))
                thin = e["both"] if e else []
                short[(u["donor"], u.get("group"))] = (u["both"] < MIN_CODETECTED
                                                       or (bool(thin) and all(b < MIN_CODETECTED for b in thin)))
            if r["design"] == "paired":
                donors = {d for d, _ in short}
                left = sum(1 for d in donors if not any(v for (dd, _), v in short.items() if dd == d))
            else:
                left = min(sum(1 for (d, g), v in short.items() if g == grp and not v) for grp in ("A", "B"))
            k += left < 3
        out.append(f"| {o} | {len(sel)} | {pct(k, len(sel))} |")
    out += ["", "By condition (all levels): donor values at the fallback, uncorrected, every claim card", "",
            "| condition | cards | donor values | at the fallback |", "|---|---|---|---|"]
    for c in sorted({r["condition"] for r in ok}):
        sel = [r for r in ok if r["condition"] == c]
        n, fb3, fbv, _ = tally(sel)
        out.append(f"| {c} | {len(sel)} | {n} | {pct(fb3 + fbv, n)} |")
    print("\n".join(out))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--results", type=Path, required=True)
    r.add_argument("--compact", type=Path, required=True)
    r.add_argument("--pilot", type=Path, required=True)
    r.add_argument("--engine-src", type=Path, required=True)
    r.add_argument("--out", type=Path, required=True)
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--datasets", type=int, default=0)
    p = sub.add_parser("report")
    p.add_argument("--results", type=Path, required=True)
    p.add_argument("--pilot", type=Path, required=True)
    p.add_argument("--rows", type=Path, required=True)
    args = ap.parse_args(argv)
    return run(args) if args.cmd == "run" else report(args)


if __name__ == "__main__":
    raise SystemExit(main())
