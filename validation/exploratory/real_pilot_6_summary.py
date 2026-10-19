#!/usr/bin/env python3
"""Preprint §4 on the fourth run of v0.4 (exploratory, written after the run; the owner's decisions of
2026-10-19, item 1): what the engine's reports of A3 (§4 in males) and A4 (§4 within sex) say about the
claim that `mi_3bin(Smad3, Col1a1)` falls from the 3-month to the 24-month mice of Tabula Muris Senis FACS,
and what the TMS authors write about the 3-month data.

Reads the run's files only (`results/real-pilot-6`), each checked first against the run's SHA256SUMS: the
engine's reports of A3 and A4, the authors' quotes, and case 1's source metadata (for the cells of the
composite unit). Computes from them, with the engine's own estimate and test (`metric_autopsy.stats`, as
the run's whole-mixture line), the mixture's estimate without the composite unit, after checking that
the same code gives the run's estimate with it; every other number is the run's, as its JSON holds it.
Nothing is re-run on the counts.

    git archive origin/results/real-pilot-6 | tar -x -C RUN_DIR
    python validation/exploratory/real_pilot_6_summary.py RUN_DIR > validation/exploratory/real_pilot_6_summary.md
"""
from __future__ import annotations

import csv
import difflib
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np

from metric_autopsy import __version__ as ENGINE
from metric_autopsy.stats import nested_t, permutation_test_nested

FILES = {"A3": "run/parts/part-case1m/reports/A3.json", "A4": "run/parts/part-case1/reports/A4.json",
         "text": "text/real_pilot_6_text.json", "obs": "extract/source_obs_case1.csv.gz",
         "readme": "README.md"}
RESULTS_COMMIT = "1e2f82402f8ced3f3a7e919503674ebc98f40ea9"  # results/real-pilot-6, the run's one commit
COMPOSITE = "3_10_M/3_11_M"
FIBROBLASTS = ("fibroblast", "fibroblast of cardiac tissue")
# the plan's rule for the authors' text (real_pilot_6_plan.md, "What the TMS authors write"): a sentence
# names the 3-month data or time point (TEXT_PATTERNS 0 or 1) together with the Tabula Muris without
# "Senis" (3) or with "previously published" (2)
P3M = (r"\b3[- ]?(?:m|mo|month)s?\b", r"\bthree[- ]month")
PORIGIN = (r"previously (?:published|described|generated|reported)", r"\bTabula Muris\b(?! Senis)")
# the sentences quoted in full below: the origin and its use as the annotation's reference, and the
# batch correction
SHOWN = (r"previously published 3m", r"as a reference for the automated pipeline", r"\bbbknn\b")


def f(x, n=4, sign=True):
    if x is None:
        return "—"
    try:
        x = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(x):
        return "—"
    return f"{x:+.{n}g}" if sign else f"{x:.{n}g}"


def ci(c, n=3):
    return f"[{f(c[0], n)}, {f(c[1], n)}]" if c else "—"


def pct(x):
    return "—" if x is None else f"{100 * float(x):.1f}%"


def units(n: dict) -> str:
    return " against ".join(f"{v} {k}" for k, v in n.items()) + " units"


def verify(run: Path) -> dict:
    sums = {}
    for line in (run / "SHA256SUMS").read_text().splitlines():
        h, name = line.split(None, 1)
        sums[name.strip().removeprefix("./")] = h
    out = {}
    for key, rel in FILES.items():
        h = hashlib.sha256((run / rel).read_bytes()).hexdigest()
        if sums.get(rel) != h:
            raise SystemExit(f"{rel}: sha256 {h} is not the run's ({sums.get(rel)})")
        out[key] = h
    return out


def mixture_estimate(values: list, groups) -> dict:
    """The whole mixture's estimate and test as the run's (`effect.py`: replicates nested in groups, one
    stratum: Welch's t for the interval, the exact permutation over the replicates for p)."""
    y = np.array([v["corrected"] for v in values], float)
    g = np.array([1 if v["group"] == groups[0] else 0 for v in values], int)
    s = np.zeros(len(y), int)
    t = nested_t(y, g, s)
    perm = permutation_test_nested(y, g, s)
    return dict(effect=t.est, ci=t.ci(0.95), se=t.se, df=t.df, p=perm.p, method=perm.method,
                n_relabelings=perm.n, min_p=perm.min_attainable_p, n={groups[0]: int(g.sum()), groups[1]: int((g == 0).sum())})


def composite_cells(path: Path) -> dict:
    """The composite unit's cells in case 1's source obs: their tissues, and the plate and mouse each
    cell's name gives (Tabula Muris 2018's names: well.plate.mouse.1.1)."""
    tissues, plates, mice, n = {}, {}, {}, 0
    with gzip.open(path, "rt", newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("donor_id") != COMPOSITE:
                continue
            n += 1
            tissues[row.get("tissue")] = tissues.get(row.get("tissue"), 0) + 1
            parts = (row.get("cell") or "").split(".")
            if len(parts) >= 3:
                plates[parts[1]] = plates.get(parts[1], 0) + 1
                mice[parts[2]] = mice.get(parts[2], 0) + 1
    return dict(n=n, tissues=tissues, plates=plates, mice=mice)


def type_rows(per_type: list, genes, groups, in_domain: bool) -> list:
    rows = []
    for t in per_type:
        if bool(t.get("in_domain")) != in_domain:
            continue
        c, r, d = t.get("corrected") or {}, t.get("raw") or {}, t.get("detection") or {}
        rows.append(dict(type=t["type"], units=t.get("units"), effect=c.get("effect"), p=c.get("p"),
                         tier=c.get("tier"), raw=r.get("effect"), ci=c.get("ci"),
                         det={(gr, g): (d.get(gr) or {}).get(g) for gr in groups for g in genes},
                         index={gr: (t.get("index") or {}).get(gr) for gr in groups}))
    return rows


def table(rows: list, genes, groups) -> list:
    head = ("| type | units | effect at equal depth | p | raw | "
            + " | ".join(f"{g} {gr}" for g in genes for gr in groups) + " |")
    lines = [head, "|" + "---|" * (5 + len(genes) * len(groups))]
    for r in rows:
        dets = " | ".join(pct(r["det"][(gr, g)]) for g in genes for gr in groups)
        lines.append(f"| {r['type']} | {r['units']} | {f(r['effect'])} | {f(r['p'], 3, sign=False)} ({r['tier']}) | "
                     f"{f(r['raw'])} | {dets} |")
    return lines


def count(rows, pred) -> int:
    return sum(1 for r in rows if pred(r))


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    run = Path(argv[0])
    hashes = verify(run)
    a3 = json.loads((run / FILES["A3"]).read_text())
    a4 = json.loads((run / FILES["A4"]).read_text())
    text = json.loads((run / FILES["text"]).read_text())
    groups = a3["params"]["groups"]
    genes = a3["params"]["gene_pair"]
    eff = a3["fields"]["effect"]["detail"]
    des = a3["fields"]["design_adequacy"]["detail"]
    mix = eff["mixture"]
    dom = eff["domain"]
    out = []
    w = out.append

    w("# Preprint §4 on the fourth run of v0.4 (exploratory)")
    w("")
    w(f"Written after the run by `real_pilot_6_summary.py` from the run's JSON (`results/real-pilot-6` at "
      f"{RESULTS_COMMIT[:7]}; workflow run 39163054704 on 57980d5, the engine v0.4.0.dev0; the plan "
      f"`real_pilot_6_plan.md`). The claim (preprint §4, the setup of the v0.1.1 audit): "
      f"`{a3['metric']}({', '.join(genes)})` falls from the 3-month to the 24-month mice of Tabula Muris "
      f"Senis FACS. The differences below are {groups[0]} − {groups[1]}, the engine's convention "
      f"(groups[0] − groups[1]): the claim's direction is a positive difference. A3 is the claim in males, "
      f"A4 within sex. Only the mixture without the composite unit is computed here (the engine "
      f"{ENGINE}'s `nested_t` and `permutation_test_nested`, as the run's mixture line); every other number "
      f"is the run's. Exploratory: not in the preprint, not in the manuscript.")
    w("")
    w("The files read, as the run's SHA256SUMS: " + "; ".join(f"`{FILES[k]}` {hashes[k][:12]}…" for k in FILES) + ".")
    w("")
    for cid, rep in (("A3", a3), ("A4", a4)):
        e = rep["fields"]["effect"]["detail"]
        w(f"- **{cid}:** {rep['verdict'].split(' — ')[0]} [cause `{rep['cause']}`]: within the domain "
          f"{f(e['effect'])} {ci(e['ci'])}, exact p {f(e['p_perm'], 3, sign=False)}.")
    w("")

    # ------------------------------------------------------------------ the mixture
    w("## 1. The whole mixture")
    w("")
    c, r = mix["corrected"], mix["raw"]
    w(f"At equal depth {f(c['effect'])} {ci(c['ci'])} (raw {f(r['effect'])}), exact p {f(c['p'], 3, sign=False)} "
      f"over {c['test'].split('n=')[-1].rstrip(')')} relabelings, {units(c['n_replicates'])}"
      + (": the smallest two-sided p these units can reach (2/126)" if abs(c["p"] - 2 / 126) < 1e-12 else "")
      + ". The units, at equal depth:")
    w("")
    group_of = {v["replicate"]: v["group"] for v in eff["replicate_values"]}
    comp = {x["unit"]: x for x in des.get("composite_units") or []}
    for v in mix["replicate_values"]:
        v["group"] = group_of.get(v["replicate"]) or (groups[0] if v["replicate"] in comp else None)
    vals = sorted(mix["replicate_values"], key=lambda v: (v["group"] != groups[0], -v["corrected"]))
    w("| unit | group | at equal depth | raw | note |")
    w("|---|---|---|---|---|")
    top = max(vals, key=lambda v: v["corrected"])["replicate"]
    for v in vals:
        note = []
        if v["replicate"] in comp:
            note.append("composite: its id names " + " and ".join(comp[v["replicate"]]["names"])
                        + ", each also a unit of its own")
        if v["replicate"] == top:
            note.append("the largest value")
        w(f"| {v['replicate']} | {v['group']} | {f(v['corrected'])} | {f(v['raw'])} | {'; '.join(note)} |")
    w("")
    full = mixture_estimate(vals, groups)
    same = (abs(full["effect"] - c["effect"]) < 1e-12 and abs(full["p"] - c["p"]) < 1e-12)
    w(f"The same code on the nine units gives {f(full['effect'])} {ci(full['ci'])}, exact p "
      f"{f(full['p'], 3, sign=False)} — {'the run' + chr(39) + 's estimate' if same else 'NOT the run' + chr(39) + 's estimate'}.")
    w("")
    cells = composite_cells(run / FILES["obs"])
    w(f"The composite unit `{COMPOSITE}` (the engine's flag `COMPOSITE_UNITS`) holds {cells['n']} cells, "
      f"all in {', '.join(f'{k} ({v})' for k, v in cells['tissues'].items())}, on plate(s) "
      f"{', '.join(f'{k} ({v})' for k, v in cells['plates'].items())}; each cell's name gives one mouse: "
      f"{', '.join(f'{k} {v}' for k, v in sorted(cells['mice'].items()))}. It enters only the mixture: no "
      f"type has 20 of its cells, so the domain and the types outside it are estimated over 4 against 4 units.")
    w("")
    rest = [v for v in vals if v["replicate"] not in comp]
    wo = mixture_estimate(rest, groups)
    w(f"**Without the composite unit** ({units(wo['n'])}): {f(wo['effect'])} {ci(wo['ci'])}, exact p "
      f"{f(wo['p'], 3, sign=False)} over {wo['n_relabelings']} relabelings"
      + (f" — again the smallest p these units can reach ({f(wo['min_p'], 3, sign=False)})"
         if abs(wo["p"] - wo["min_p"]) < 1e-12 else "")
      + (": every young unit lies above every old one" if min(v["corrected"] for v in rest if v["group"] == groups[0])
         > max(v["corrected"] for v in rest if v["group"] == groups[1]) else "")
      + f". The mixture's difference is smaller by {f(c['effect'] - wo['effect'], 3, sign=False)}.")
    w("")

    # ------------------------------------------------------------------ the domain
    w("## 2. Within the domain")
    w("")
    p = des["power"]
    w(f"The domain by the rule \"in each group\" ({dom['rule']['text']}): {', '.join(dom['types'])}. Combined, "
      f"the types weighted equally: {f(eff['effect'])} {ci(eff['ci'])} (raw {f(eff['raw_effect'])}), exact p "
      f"{f(eff['p_perm'], 3, sign=False)} over {eff['perm_n']} relabelings, {units(eff['n_replicates'])}: "
      f"{a3['fields']['effect']['status']}. The MDE at 80% power is {f(p['mde'], 3, sign=False)}, larger than "
      f"the mixture's effect ({f(c['effect'], 3)}, without the composite unit {f(wo['effect'], 3)}): within the "
      f"domain this design could not have detected a drop the size of the mixture's. The SESOI in the metric's "
      f"units is s_m {f(p['sesoi'], 3)}, below the design's floor {f(p['floor'], 3, sign=False)}.")
    w("")
    rows_in = type_rows(dom["per_type"], genes, groups, True)
    for ln in table(rows_in, genes, groups):
        w(ln)
    w("")
    fib = [r for r in rows_in if r["type"] in FIBROBLASTS]
    neg = all((r["effect"] or 0) < 0 for r in fib)
    w(("In both fibroblast types the difference is negative, against the claim: " if neg else
       "The fibroblast types: ")
      + "; ".join(f"{r['type']} {f(r['effect'])} (p {f(r['p'], 3, sign=False)})" for r in fib)
      + f". Positive in {count(rows_in, lambda r: (r['effect'] or 0) > 0)} of {len(rows_in)} types of the domain.")
    w("")

    # ------------------------------------------------------------------ outside
    w("## 3. Outside the domain")
    w("")
    o = dom["outside"]
    w(f"Combined over the {len(o['types'])} types outside the domain, the types weighted equally: "
      f"{f(o['corrected']['effect'])} {ci(o['corrected']['ci'])} (raw {f(o['raw']['effect'])}), exact p "
      f"{f(o['corrected']['p'], 3, sign=False)} over {o['corrected']['test'].split('n=')[-1].rstrip(')')} "
      f"relabelings, {units(o['corrected']['n_replicates'])}"
      + (" (the smallest p they can reach)" if abs(o["corrected"]["p"] - 2 / 70) < 1e-12 and
         o["corrected"]["n_replicates"] == {groups[0]: 4, groups[1]: 4} else "")
      + ". Per type, the claim genes' detection (the share of the type's cells with a count) in each group:")
    w("")
    rows_out = type_rows(dom["per_type"], genes, groups, False)
    for ln in table(rows_out, genes, groups):
        w(ln)
    w("")
    n = len(rows_out)
    gy, go = groups

    def more(gene):
        return count(rows_out, lambda r: (r["det"][(gy, gene)] or 0) > (r["det"][(go, gene)] or 0))
    w(f"- The sign: positive (the claim's direction) in {count(rows_out, lambda r: (r['effect'] or 0) > 0)} of {n} "
      f"types at equal depth ({count(rows_out, lambda r: (r['raw'] or 0) > 0)} of {n} raw); nominally p < 0.05 "
      f"in {count(rows_out, lambda r: r['p'] is not None and r['p'] < 0.05)}.")
    w(f"- {genes[1]} is detected more often in the {gy} in {more(genes[1])} of {n} types, {genes[0]} in "
      f"{more(genes[0])} of {n}; the index of the claim genes' counts per cell is higher in the {gy} in "
      f"{count(rows_out, lambda r: (r['index'][gy] or 0) > (r['index'][go] or 0))} of {n}.")
    hi = [r for r in rows_out if (r["det"][(gy, genes[1])] or 0) >= 0.05 and (r["det"][(go, genes[1])] or 0) < 0.02]
    w(f"- Outside the domain {genes[1]}, a gene of fibroblasts, is detected in more than 5% of the {gy} cells "
      f"but under 2% of the {go} in {len(hi)} types, among them "
      + ", ".join(f"{r['type']} ({pct(r['det'][(gy, genes[1])])} against {pct(r['det'][(go, genes[1])])})"
                  for r in hi[:6]) + ".")
    w("")

    # ------------------------------------------------------------------ technical covariates
    w("## 4. The technical covariates per unit (the dossier's line)")
    w("")
    w(f"The median cell of each unit, the {groups[0]} units' mean against the {groups[1]}'s, with the exact "
      f"permutation p (`Autopsy.technical`; genes and counts over the pool's {a3['specificity']['pool']['genes']:,} genes).")
    w("")
    w("| claim | cells | covariate | stratum | " + f"{groups[0]} | {groups[1]} | difference | p (relabelings) | "
      f"units, {groups[0]} / {groups[1]} |")
    w("|---|---|---|---|---|---|---|---|---|")
    for cid, rep in (("A3", a3), ("A4", a4)):
        for t in rep["technical"]["rows"]:
            if t.get("difference") is None:
                continue
            w(f"| {cid} | {t['cells']} | {t['covariate']} | {t['stratum']} | {t['first']:,.1f} | {t['second']:,.1f} | "
              f"{t['difference']:+,.1f} | {f(t['p'], 3, sign=False)} ({t['n_relabelings']}) | "
              f"{' / '.join(str(v) for v in t['replicates'].values())} |")
    w("")
    sig = [(cid, t) for cid, rep in (("A3", a3), ("A4", a4)) for t in rep["technical"]["rows"]
           if t.get("p") is not None and t["p"] < 0.05]
    w("With p < 0.05: " + "; ".join(f"{cid} {t['cells']}, {t['covariate']}, {t['stratum']} ({t['difference']:+,.0f})"
                                   for cid, t in sig)
      + f"; no female row. GATE 1's ratio of median genes passed the males at "
      f"{f(a3['gates'][1]['detail']['worst_ratio'], 3, sign=False)}.")
    w("")

    # ------------------------------------------------------------------ GATE 0, GATE 4, the correction
    w("## 5. Attenuation, the probe and the correction on §4")
    w("")
    g0 = a3["gates"][0]["detail"]
    xd = g0["responses"]["extra_dropout"]
    g4 = a3["gates"][4]["detail"]
    eq = mix["equalization"]["capture_check"]
    eqd = dom["capture_check"]
    w(f"- GATE 0: extra dropout removes {pct(xd['signal_loss'])} of the signal above the null (an attenuation, not "
      f"a bias: the null moves by {f(xd['null_shift'], 3)}, z {f(xd['null_shift_z'], 3, sign=False)}, against the "
      f"tolerance {f(g0['bias_tolerance_abs'], 3, sign=False)} = {g0['bias_tolerance']} × s_m). GATE 0 passes an "
      f"attenuation on as a loss of power, but where the groups differ in detection (section 4) it can create a "
      f"difference at equal coupling (the differential sensitivity, decided 2026-10-19).")
    w(f"- GATE 4 injected the coupling into the {g4['injected_group']} only, the group with more probed counts "
      f"per cell ({f(g4['probed_counts_per_cell'][groups[0]], 4, sign=False)} against "
      f"{f(g4['probed_counts_per_cell'][groups[1]], 4, sign=False)}): s_m {f(g4['statistic_change'], 3)} "
      f"{ci(g4['statistic_ci'])} is the {g4['injected_group']}'s response; the {groups[1]}'s was not measured.")
    w(f"- The correction thinned the {mix['equalization']['strata'][0]['thinned_group']} cells to the "
      f"{groups[0]}'s depth; the gap in genes per cell ({groups[0]} − {groups[1]}, relative) grew from "
      f"{pct(eq['raw']['relative_difference']['genes'])} to {pct(eq['equalized']['relative_difference']['genes'])} on "
      f"all cells, and from {pct(eqd['raw']['relative_difference']['genes'])} to "
      f"{pct(eqd['equalized']['relative_difference']['genes'])} on the domain's: equal depth makes the groups' "
      f"detection less equal.")
    w("")

    # ------------------------------------------------------------------ A4
    e4 = a4["fields"]["effect"]["detail"]
    m4, o4 = e4["mixture"]["corrected"], e4["domain"]["outside"]["corrected"]
    w("## 6. Within sex (A4)")
    w("")
    strata4 = {st["stratum"]: st for st in a4["composition"]["domain"]["strata"]}
    dom4 = a4["composition"]["domain"]["domain"]
    fem = {t["type"]: t["cells"] for t in strata4["female"]["types"] if t["type"] in dom4["male"]}
    w(f"The domain by stratum: males {', '.join(dom4['male'])}; females none — their {groups[1]} units have no "
      f"cell of these types ({'; '.join(f'{k} {v[groups[0]]} against {v[groups[1]]}' for k, v in fem.items())}"
      f" cells), so the domain's estimate is the males'. Within the domain {f(e4['effect'])} {ci(e4['ci'])} "
      f"(exact p {f(e4['p_perm'], 3, sign=False)}); the mixture {f(m4['effect'])} (p {f(m4['p'], 3, sign=False)}, "
      f"{m4['test']}); outside the domain {f(o4['effect'])} (p {f(o4['p'], 3, sign=False)}, {o4['test']}). GATE "
      f"4's response, injected into the {a4['gates'][4]['detail']['injected_group']}, was not established "
      f"({a4['gates'][4]['status']}).")
    w("")

    # ------------------------------------------------------------------ the authors' text
    w("## 7. What the TMS authors write")
    w("")
    for s in text["sources"]:
        extra = f"; {s['bytes']:,} bytes, sha256 {s['sha256'][:12]}…" if s.get("sha256") else ""
        ver = f"; version {s['version']}, licence {s['license']}" if s.get("version") else ""
        w(f"- {s['source']} ({s['what']}): {s['outcome']}{ver}{extra}.")
    w("")
    states = []
    for s in text["sources"]:
        for q in s.get("quotes") or []:
            m = q.get("matched") or []
            if any(x in m for x in P3M) and any(x in m for x in PORIGIN):
                states.append((s["source"], q))
    w(f"By the plan's rule (a sentence that names the 3-month data together with the Tabula Muris or with "
      f"\"previously published\") the authors state that the 3-month data are the original Tabula Muris's: "
      f"{len(states)} such sentences, in {len({x[0] for x in states})} sources. Quoted, with the batch "
      f"correction:")
    w("")
    shown = {s["source"]: [q for q in s.get("quotes") or [] if any(re.search(pat, q["text"]) for pat in SHOWN)]
             for s in text["sources"]}
    paper = shown.get("tms-bioc") or []
    for q in paper:
        w(f"> {q['text']}")
        w(">")
        w(f"> — the Nature paper (NCBI's BioC), {q['section']}")
        w("")
    pre = shown.get("tms-biorxiv") or []
    if pre:
        same = []
        for q in pre:
            best = max((difflib.SequenceMatcher(None, q["text"], x["text"]).ratio() for x in paper), default=0.0)
            same.append(f"{q['section']} ({best:.2f})")
        w(f"The preprint (bioRxiv, version {next(x for x in text['sources'] if x['source'] == 'tms-biorxiv')['version']}) "
          f"has {len(pre)} such sentences, each matched to the paper's by its similarity (difflib's ratio, 1 = the "
          f"same text; the figure and reference numbers differ): " + "; ".join(same) + ".")
        w("")
    w("bbknn corrects the neighbourhood graph between the two technologies (FACS and droplet); it does not "
      "change the counts, which the engine reads, and no quoted sentence names a correction between the "
      "3-month data and the later time points. The later time points' cells were annotated with the 3-month "
      "cells as the reference: a comparison within a cell type compares labels given by two procedures.")
    w("")

    # ------------------------------------------------------------------ conclusion
    w("## Conclusion")
    w("")
    w("Most likely a difference in co-detection between two generations of data (Tabula Muris 2018 against "
      "Tabula Muris Senis), not a coupling; the engine does not establish this — the decisive steps are the "
      "biology (in which cells Col1a1 is expressed) and the authors' text. (The owner's wording, 2026-10-19: "
      "«скорее всего, различие в со-детекции между двумя поколениями данных, а не связь; движок этого не "
      "устанавливает — решающие шаги: биология (в каких клетках экспрессируется Col1a1) и текст авторов».)")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
