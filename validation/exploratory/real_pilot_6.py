#!/usr/bin/env python3
"""The fourth run of v0.4 on the development data (validation/exploratory/real_pilot_6_plan.md): the engine
after the owner's second decisions of 2026-10-18 (the full dossier, change 12, the domain of a pair in each
group, norm_pearson's NaN) on the third run's claims. Exploratory, for development data only.

Parts (the workflow runs them as parallel jobs), each ``run_autopsy`` of the engine installed (this pull
request's src/) on the third run's claims, loaded and called as the third run did (real_pilot_5.py):

- ``case2``: A1 (Xist) and A2 (the Y genes), female against male at 3 months;
- ``case1m``: A3, §4 in males, the recorded second attempt; ``case1``: A4, §4 within sex, the same;
- ``case3a``, ``case3b``, ``case3c``: A8, the 20 claims of SEA-AD, a technical negative control;
- ``--merge``: the parts' outputs into real_pilot_6.json and real_pilot_6.log, every criterion and the
  owner's prediction scored as the plan fixes them, with the dossier's lines of each claim.

    python validation/exploratory/real_pilot_6.py --inputs INPUTS --full FULL --part PART --out OUT
    python validation/exploratory/real_pilot_6.py --merge PARTS_DIR --out OUT [--synthetic]

INPUTS and FULL are the third run's (real_pilot_5_inputs.sh; real_pilot_6_extract.py). Every exception is
caught per claim and reported with its traceback; a part that raises still writes what it has.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import real_pilot_5 as R5  # noqa: E402  (pins v1's numerics before numpy; the loaders and the engine's call)
import real_pilot_6_spec as S  # noqa: E402

np, pd, MA = R5.np, R5.pd, R5.MA


# --------------------------------------------------------------------------- #
# the dossier's lines of a report
# --------------------------------------------------------------------------- #
def dossier(rep: dict) -> dict:
    """What the engine reports beside the verdict since the second decisions: the effect after a FAIL
    (DOSSIER_ONLY), GATE 0's test between the groups (change 12), the unmeasured replicates and units,
    the types where a pair reaches the threshold in one group only, and the technical covariates."""
    ef = rep["fields"].get("effect") or {}
    det = ef.get("detail") or {}
    g0 = next((g for g in rep["gates"] if g["gate"] == 0), {})
    g0d = g0.get("detail") or {}
    b = g0d.get("between_groups")
    dom = (rep.get("composition") or {}).get("domain") or {}
    return dict(
        dossier_only="DOSSIER_ONLY" in (ef.get("flags") or []), effect_status=ef.get("status"),
        effect_reason=ef.get("reason"),
        gate0=dict(status=g0.get("status"), message=g0.get("message"), handling=g0d.get("bias_handling"),
                   tolerance=g0d.get("bias_tolerance_abs"),
                   between_groups=b and {k: b.get(k) for k in ("equivalent", "margin", "effect", "ci_2alpha", "status",
                                                              "test", "estimate_over", "n_replicates", "error")}),
        unmeasured=det.get("unmeasured"), domain_unmeasured=(det.get("domain") or {}).get("unmeasured"),
        domain_rule=(dom.get("rule") or {}).get("text"), one_group=dom.get("one_group"),
        technical=rep.get("technical") or {})


def tech_row(rec: dict, cells: str, covariate: str, stratum: str) -> dict | None:
    for r in (rec.get("dossier") or {}).get("technical", {}).get("rows") or []:
        if (r["cells"], r["covariate"], r["stratum"]) == (cells, covariate, stratum):
            return r
    return None


# --------------------------------------------------------------------------- #
# the parts
# --------------------------------------------------------------------------- #
def claims_of(part: str, inputs: Path) -> list:
    if part.startswith("case3"):
        return [dict(c, criterion=S.CRITERIA["A8"], attempt=1) for c in R5.case3_claims(inputs, part[4:])]
    return [c for c in S.ACCEPT if c["part"] == part]


def run_part(args, part: str, log) -> dict:
    out = dict(claims=[])
    name = {"case2": "case2", "case1m": "case1", "case1": "case1"}.get(part, "case3")
    males = (lambda o: o["sex"].astype(str).to_numpy() == "male") if part == "case1m" else None
    data, pool, sha = R5.load_case(args.inputs, args.full, name, args.synthetic, rows=males)
    out["case_sha256"] = sha
    for claim in claims_of(part, args.inputs):
        try:
            rec = R5.accept_run(claim, data, pool, args.out, log)
            rep = json.loads((args.out / rec["report"]).read_text())
            rec["dossier"] = dossier(rep)
            d = rec["dossier"]
            log(f"    dossier: effect {d['effect_status']}{' (DOSSIER_ONLY)' if d['dossier_only'] else ''}; GATE 0 "
                f"{d['gate0']['status']}, between the groups {d['gate0']['between_groups']}; unmeasured "
                f"{len(d['unmeasured'] or [])} replicates, {len(d['domain_unmeasured'] or [])} units of the domain; "
                f"one group only {[(x['stratum'], x['type']) for x in d['one_group'] or []]}")
            out["claims"].append(rec)
        except Exception as exc:
            out["claims"].append(dict(id=claim["id"], claim=claim, error=repr(exc),
                                      traceback=traceback.format_exc().splitlines()))
            log(f"- {claim['id']}: EXCEPTION {exc!r}")
    return out


# --------------------------------------------------------------------------- #
# the scores
# --------------------------------------------------------------------------- #
def score(recs: dict) -> dict:
    def ok(i):
        r = recs.get(i)
        return r is not None and "error" not in r
    out = {}
    for i in ("A1", "A2"):
        lab = recs[i]["label"] if ok(i) else None
        out[i] = dict(criterion=S.CRITERIA[i], label=lab, met=(lab == "SUPPORTED") if lab else None)
    for i in ("A3", "A4"):
        if not ok(i):
            out[i] = dict(criterion=S.CRITERIA[i], label=None, met=None)
            continue
        r, d = recs[i], recs[i]["dossier"]
        has = bool(d["effect_status"] is not None and d["technical"].get("rows") and r.get("domain"))
        out[i] = dict(criterion=S.CRITERIA[i], label=r["label"], not_supported=r["label"] != "SUPPORTED",
                      dossier=has, met=bool(r["label"] != "SUPPORTED" and has))
    pred = {}
    for i, stratum in (("A3", "all"), ("A4", "male")):
        if not ok(i):
            pred[i] = None
            continue
        r = recs[i]
        g0 = r["dossier"]["gate0"]
        blocks = bool(g0["status"] == "FAIL" and r["cause"] == "metric_invalid_gate0"
                      and not (g0.get("between_groups") or {}).get("equivalent"))
        genes = tech_row(r, "all cells", "genes per cell", stratum)
        counts = tech_row(r, "all cells", "counts per cell", stratum)
        differ = bool(genes and genes.get("difference") is not None and genes["difference"] > 0
                      and genes["p"] < S.ALPHA and counts and counts.get("difference") is not None
                      and counts["difference"] < 0 and counts["p"] < S.ALPHA)
        pred[i] = dict(gate0_blocks=blocks, covariates_differ=differ, genes=genes, counts=counts)
    both = [v for v in pred.values() if v]
    p1 = len(both) == 2 and all(v["gate0_blocks"] for v in both)
    p2 = len(both) == 2 and all(v["covariates_differ"] for v in both)
    out["prediction"] = dict(text=S.PREDICTION, per_claim=pred, gate0_blocks=p1, covariates_differ=p2,
                             score=("holds" if p1 and p2 else "partly holds" if (p1 or p2) else "does not hold")
                             if len(both) == 2 else "not scored: a claim raised")
    a8 = {k: v for k, v in recs.items() if k.startswith("A8-")}
    labels = {k: (v.get("label") if "error" not in v else "EXCEPTION") for k, v in a8.items()}
    out["A8"] = dict(criterion=S.CRITERIA["A8"], n=len(a8), labels=labels,
                     supported=sorted(k for k, lb in labels.items() if lb == "SUPPORTED"),
                     exceptions=sorted(k for k, lb in labels.items() if lb == "EXCEPTION"),
                     met=bool(len(a8) == 20 and all(lb not in ("SUPPORTED", "EXCEPTION") for lb in labels.values())),
                     changed_from_the_third_run=sorted(
                         k for k, lb in labels.items()
                         if lb != S.A8_THIRD_RUN.get(k[3:], "NO DETECTABLE EFFECT")))
    dom = (recs.get("A3") or {}).get("domain") or {}
    out["A3_domain"] = dict(expected=list(S.A3_DOMAIN_EXPECTED), types=dom.get("types"),
                            as_expected=sorted(dom.get("types") or []) == sorted(S.A3_DOMAIN_EXPECTED))
    return out


# --------------------------------------------------------------------------- #
# the merge
# --------------------------------------------------------------------------- #
def fmt(x, n=4):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:+.{n}g}"


def fx(s):
    if not s or s.get("effect") is None:
        return "—"
    return f"{fmt(s['effect'])} (p {s.get('p'):.3g}, {s.get('tier')})" if s.get("p") is not None else fmt(s["effect"])


def claim_lines(r: dict) -> list:
    if "error" in r:
        return [f"### {r['id']}: EXCEPTION {r['error']}"]
    d, ef = r["dossier"], r["effect"]
    g0 = d["gate0"]
    b = g0.get("between_groups") or {}
    out = [f"### {r['id']} (attempt {r['claim'].get('attempt', 1)}): {r['verdict']} [cause {r['cause']}]", "",
           "- gates: " + "; ".join(f"{g['gate']} {g['status']}" for g in r["gates"]),
           f"- GATE 0: {g0['status']}, handling {g0['handling']}, tolerance {fmt(g0['tolerance'])}"
           + (f"; the groups' technical difference {fmt(b.get('effect'))}, 90% interval "
              f"[{fmt((b.get('ci_2alpha') or [None, None])[0])}, {fmt((b.get('ci_2alpha') or [None, None])[1])}] "
              f"against ±{fmt(b.get('margin'))}: {'within' if b.get('equivalent') else 'not shown within'}"
              f" ({b.get('test')}; {b.get('estimate_over')})" if b else ""),
           f"- s_m: {fmt((r.get('gate4') or {}).get('sesoi_metric'))} {(r.get('gate4') or {}).get('sesoi_metric_ci')}",
           f"- the effect ({d['effect_status']}{', the dossier' if d['dossier_only'] else ''}): {fmt(ef.get('effect'))} "
           f"CI {ef.get('ci')}, raw {fmt(ef.get('raw_effect'))}, over {ef.get('estimate_over')}"]
    dom = r.get("domain") or {}
    if dom:
        out.append(f"- the domain: {dom.get('types')} ({d['domain_rule']})")
        mix = r.get("mixture") or {}
        out.append(f"- the whole mixture: corrected {fx(mix.get('corrected'))}, raw {fx(mix.get('raw'))}; outside the "
                   f"domain: corrected {fx((dom.get('outside') or {}).get('corrected'))}, raw "
                   f"{fx((dom.get('outside') or {}).get('raw'))}")
        for t in dom.get("per_type") or []:
            out.append(f"  - {t['type']} ({'in' if t['in_domain'] else 'outside'}): corrected {fx(t.get('corrected'))}, "
                       f"raw {fx(t.get('raw'))}; detection {t.get('detection')}; counts index {t.get('index')}")
    for x in d.get("one_group") or []:
        out.append(f"- one group only: {x['stratum']}: {x['type']} detected in {x['detected_in']} ({x['detection']}; "
                   f"index {x['index']})")
    if d.get("unmeasured") or d.get("domain_unmeasured"):
        out.append(f"- unmeasured: replicates {d.get('unmeasured')}; units of the domain "
                   f"{[(u['replicate'], u['type']) for u in d.get('domain_unmeasured') or [] if u.get('in_domain')]}")
    tech = d.get("technical") or {}
    for t in tech.get("rows") or []:
        if t.get("difference") is None:
            out.append(f"- technical, {t['cells']}, {t['covariate']}, {t['stratum']}: not compared ({t.get('reason')})")
        else:
            out.append(f"- technical, {t['cells']}, {t['covariate']}, {t['stratum']}: {fmt(t['first'], 6)} against "
                       f"{fmt(t['second'], 6)}, difference {fmt(t['difference'], 6)} ({t['method']} p {t['p']:.3g}, "
                       f"{t['n_relabelings']} relabelings, {t['replicates']})")
    if tech.get("error"):
        out.append(f"- technical: {tech['error']}")
    sp = r.get("specificity") or {}
    if sp.get("in_direction"):
        out.append(f"- specificity: {sp['in_direction'].get('k')} of {sp.get('n_controls')} controls at least as far "
                   f"in the direction; the claim's t {fmt(sp.get('t_claim'))}")
    return out + [""]


def merge(parts_dir: Path, out_dir: Path, synthetic: bool) -> int:
    parts, logs = {}, []
    for part in S.PARTS:
        f = parts_dir / f"part-{part}" / f"{part}.json"
        if f.exists():
            parts[part] = json.loads(f.read_text())
            lg = parts_dir / f"part-{part}" / f"{part}.log"
            logs += lg.read_text().splitlines() if lg.exists() else []
        else:
            parts[part] = dict(missing=True)
    recs = {r["id"]: r for p in S.PARTS for r in (parts[p].get("result") or {}).get("claims", [])}
    sc = score(recs)
    tag = "SYNTHETIC " if synthetic else ""
    res = dict(synthetic=synthetic, parts={k: {kk: vv for kk, vv in v.items() if kk != "result"} for k, v in parts.items()},
               results={k: v.get("result") for k, v in parts.items()}, scores=sc)
    (out_dir / "real_pilot_6.json").write_text(json.dumps(res, indent=1, default=str))
    pr = sc["prediction"]
    lines = [f"# {tag}The fourth run of v0.4 on the development data (validation/exploratory/real_pilot_6_plan.md)", "",
             "## The scores", "", "| test | criterion | verdict | met |", "|---|---|---|---|"]
    for k in ("A1", "A2", "A3", "A4"):
        v = sc[k]
        lines.append(f"| {k} | {v['criterion']} | {v.get('label')} | {v.get('met')} |")
    a8 = sc["A8"]
    lines.append(f"| A8 | {a8['criterion']} | {a8['n']} claims; SUPPORTED: {a8['supported'] or 'none'}; exceptions: "
                 f"{a8['exceptions'] or 'none'} | {a8['met']} |")
    lines += ["", f"The owner's prediction (A3, A4): «{pr['text']}» — {pr['score']} (GATE 0 blocks on both: "
              f"{pr['gate0_blocks']}; genes per cell higher and counts per cell lower in the young males' units, each "
              f"p < {S.ALPHA}: {pr['covariates_differ']})"]
    for i, v in (pr.get("per_claim") or {}).items():
        if v:
            lines.append(f"- {i}: GATE 0 blocks {v['gate0_blocks']}; genes per cell {fmt((v['genes'] or {}).get('difference'), 6)} "
                         f"(p {(v['genes'] or {}).get('p')}); counts per cell {fmt((v['counts'] or {}).get('difference'), 6)} "
                         f"(p {(v['counts'] or {}).get('p')})")
    dm = sc["A3_domain"]
    lines += ["", f"A3's domain by the rule in each group: {dm['types']}; expected from the third run's table "
              f"{dm['expected']}: {'as expected' if dm['as_expected'] else 'NOT as expected'}",
              f"A8's verdicts that differ from the third run's: {a8['changed_from_the_third_run'] or 'none'}",
              "", "## The claims", ""]
    for k in sorted(recs):
        lines += claim_lines(recs[k])
    lines += ["## A8's verdicts", ""] + [f"- {k}: {lb}" for k, lb in sorted(a8["labels"].items())]
    lines += ["", "## The parts' logs", ""] + logs
    (out_dir / "real_pilot_6.log").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:30]))
    return 0


# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--inputs", type=Path)
    p.add_argument("--full", type=Path)
    p.add_argument("--part", choices=S.PARTS)
    p.add_argument("--merge", type=Path, help="the directory of the parts' outputs (part-<name>/)")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--synthetic", action="store_true")
    args = p.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    if args.merge:
        return merge(args.merge, args.out, args.synthetic)
    lines = []

    def log(*msg):
        for m in msg:
            print(m, flush=True)
            lines.append(str(m))
    tag = "SYNTHETIC " if args.synthetic else ""
    log(f"# {tag}part {args.part} (validation/exploratory/real_pilot_6_plan.md)", "",
        f"engine {MA.__version__} from {Path(MA.__file__).resolve().parent}; python {platform.python_version()}; numpy "
        f"{np.__version__}; pandas {pd.__version__}; numerics {R5.NUMERICS}", "")
    t0, status = time.time(), 0
    try:
        result = run_part(args, args.part, log)
    except Exception as exc:
        result = dict(error=repr(exc), traceback=traceback.format_exc().splitlines())
        log(f"EXCEPTION {exc!r}", *result["traceback"])
        status = 1
    out = dict(part=args.part, synthetic=bool(args.synthetic), engine=MA.__version__, seconds=round(time.time() - t0, 1),
               result=result)
    (args.out / f"{args.part}.json").write_text(json.dumps(out, indent=1, default=str))
    log("", f"{tag}part {args.part}: {out['seconds']} s")
    (args.out / f"{args.part}.log").write_text("\n".join(lines) + "\n")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
