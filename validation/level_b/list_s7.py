#!/usr/bin/env python3
"""Level B, S7 (`plan_s7.md`): the list for the owner — the one stop — assembled from the repository's files and the
results branches (fetched beforehand into one directory, one sub-directory per branch), every number read from them.
Nothing here reads a claim's metric.

    python -I validation/level_b/list_s7.py --results DIR --out validation/level_b/amendment1/list

DIR holds the results branches it reads by their names (`level-b-data-<n>/`, `level-b-controls/`, and, when they
exist, `level-b-verdicts/` with the passes' codes). OUT gets list.md (the owner's page) and list.json.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
A1 = HERE / "amendment1"
SET_CAVEATS = ("v1 does not cover a set's path (set_norm_pearson with the coupling_set probe)",
               "the set's definition was not checked against Levy et al. 2020")


def _json(path: Path):
    return json.loads(path.read_text())


def counts(root: Path = HERE) -> dict:
    """Part 1: the screening's counts, both screenings."""
    ta0, ta1 = _json(root / "screening" / "merge.json"), _json(root / "amendment1" / "ta" / "merge.json")
    ft0, ft1 = _json(root / "fulltext_screening" / "ft_merge.json"), _json(root / "amendment1" / "ft" / "merge" / "ft_merge.json")
    pick = lambda d, ks: {k: d[k] for k in ks}  # noqa: E731
    return dict(titles_abstracts={"2026-10-19": pick(ta0, ["n_screened", "kappa", "agreement", "both_keep", "to_owner", "fulltext"]),
                                  "amendment 1": pick(ta1, ["n_screened", "kappa", "agreement", "both_keep", "to_owner", "fulltext"])},
                full_texts={"2026-10-19": pick(ft0, ["read_by_both", "kappa", "agreement", "kept_by_both", "to_owner"]),
                            "amendment 1": pick(ft1, ["read_by_both", "kappa", "agreement", "kept_by_both", "to_owner"])})


def data_reports(results: Path) -> dict:
    """Part 4: each dataset's checks, by spec name (`<rid>` or `<rid>_<x>`), from every `level-b-data-<n>` branch."""
    out = {}
    for d in sorted(results.glob("level-b-data-*")):
        for rep in sorted(d.glob("*/report.json")):
            r = _json(rep)
            r["_branch"] = d.name
            out[rep.parent.name] = r
    return out


def origin(rep: dict) -> dict:
    """12.3's check of origin on one dataset's report: fails where the counts are not raw, an id names two
    individuals, or a source one group's cells do not share with the other's."""
    why = []
    if not rep.get("counts", {}).get("raw"):
        why.append("counts not raw integers")
    if rep.get("sex", {}).get("fails"):
        why.append(f"ids naming two individuals or against their label: {rep['sex']['fails']}")
    if rep.get("source", {}).get("verdict") == "fails":
        why.append(f"the source confounds the groups: {rep['source'].get('sources')}")
    return dict(fails=bool(why), why=why)


def donors_cannot_run(rep: dict) -> str | None:
    """12.3's test over donors cannot run where the design gives fewer donors than the graded rule needs."""
    r4 = rep.get("donors", {}).get("rule4", {})
    return None if r4.get("passes") else f"rule 4 fails ({rep.get('donors', {}).get('design')}, {rep.get('donors', {}).get('donors')})"


def controls(results: Path) -> dict | None:
    p = results / "level-b-controls" / "claims.json"
    return _json(p) if p.exists() else None


def assemble(results: Path, root: Path = HERE) -> dict:
    kept = _json(root / "amendment1" / "kept_claims.json")
    reps = data_reports(results)
    papers = []
    for p in kept["papers"]:
        mine = {k: v for k, v in reps.items() if k == p["rid"] or k.startswith(p["rid"] + "_")}
        papers.append(dict(rid=p["rid"], title=p["title"], read=p["read"], pmid=p["pmid"], doi=p["doi"],
                           claims=p["claims"], data={k: dict(report=v, origin=origin(v), cannot_run=donors_cannot_run(v))
                                                     for k, v in mine.items()}))
    fail_origin = [p["rid"] for p in papers if any(d["origin"]["fails"] for d in p["data"].values())]
    cannot = [p["rid"] for p in papers if any(d["cannot_run"] for d in p["data"].values())]
    unchecked = [p["rid"] for p in papers if not p["data"]]
    n_claims = lambda rids: sum(len(p["claims"]) for p in papers if p["rid"] in rids)  # noqa: E731
    return dict(counts=counts(root), kept_counts=kept["counts"], fewer_than_target=kept["fewer_than_target"],
                drawn=kept["drawn"], papers=papers, controls=controls(results),
                origin_fails=dict(papers=fail_origin, claims=n_claims(fail_origin)),
                donors_cannot_run=dict(papers=cannot, claims=n_claims(cannot)),
                data_not_checked=dict(papers=unchecked, claims=n_claims(unchecked)))


def markdown(L: dict) -> str:
    c = L["counts"]
    out = ["# Level B: the list for the owner (S7; `plan_s7.md`)", "",
           "**Nothing here is validated; no claim's metric was read.** Every number below is read by "
           "`list_s7.py` from the repository's files and the results branches.", "", "## 1. The search and the screening", ""]
    for step, rows in (("titles and abstracts", c["titles_abstracts"]), ("full texts", c["full_texts"])):
        for when, r in rows.items():
            out.append(f"- {step}, {when}: " + "; ".join(f"{k} {v:.3f}" if isinstance(v, float) else f"{k} {v}"
                                                        for k, v in r.items()))
    k = L["kept_counts"]
    out += ["", f"## 2. The kept papers: {k['papers']} ({k['papers_with_pairs']} with pairs, {k['papers_with_sets']} with sets)", ""]
    if L["fewer_than_target"]:
        out.append(f"**Fewer than 20 papers with pairs** (the protocol's section 9): {k['papers_with_pairs']}.")
    out += ["", "## 3–4. Each paper's claims and data", ""]
    for p in L["papers"]:
        out += [f"### {p['rid']}: {p['title']}", f"PMID {p['pmid'] or '-'}, DOI {p['doi'] or '-'}; read {p['read']}", ""]
        for kind in ("pair", "set"):
            for cl in [x for x in p["claims"] if x["kind"] == kind]:
                what = ", ".join(cl["genes"]) if kind == "pair" else f"{cl.get('set_name')} ({len(cl['genes'])} genes taken)"
                out.append(f"- {kind} {cl['id']}: **{what}** — passes {cl['passes']}")
                if kind == "set":
                    out.append("  - caveats: " + "; ".join(SET_CAVEATS))
        for name, d in p["data"].items():
            r = d["report"]
            out.append(f"- data `{name}` ({r['_branch']}): counts raw {r['counts'].get('raw')}; donors {r['donors']['design']} "
                       f"{r['donors']['donors']}, rule 4 {'passes' if r['donors']['rule4'].get('passes') else 'fails'}; "
                       f"source {r['source']['verdict']}; origin {'FAILS: ' + '; '.join(d['origin']['why']) if d['origin']['fails'] else 'passes'}")
        if not p["data"]:
            out.append("- data: not checked yet")
        out.append("")
    out += ["## 6. The controls", ""]
    ctl = L["controls"]
    if ctl:
        for n in ctl["negative"]:
            out.append(f"- {n['id']}: {n['collection_name']} — {n['tissue']}, {n['cell_type']}; {n['condition_a']} against "
                       f"{n['condition_b']}; {n['donors']} donors, {n['cells']} cells; {len(n['claims'])} claims")
        out.append(f"- positive: {len(ctl['positive'])} dataset(s); failed: {ctl['positive_failed'] or 'none'}")
    else:
        out.append("- the controls' branch not given")
    out += ["", "## 7. 12.3's two counts so far", "",
            f"- fail the check of origin: {len(L['origin_fails']['papers'])} papers, {L['origin_fails']['claims']} claims "
            f"({', '.join(L['origin_fails']['papers']) or 'none'})",
            f"- the test over donors cannot run: {len(L['donors_cannot_run']['papers'])} papers, "
            f"{L['donors_cannot_run']['claims']} claims ({', '.join(L['donors_cannot_run']['papers']) or 'none'})",
            f"- data not checked: {len(L['data_not_checked']['papers'])} papers ({', '.join(L['data_not_checked']['papers']) or 'none'})",
            "", "Whether the test over donors supports a claim is known only after the runs."]
    return "\n".join(out) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    L = assemble(a.results)
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "list.json").write_text(json.dumps(L, indent=1, ensure_ascii=False, default=str) + "\n")
    (a.out / "list.md").write_text(markdown(L))
    print(json.dumps(dict(papers=len(L["papers"]), origin_fails=L["origin_fails"], donors_cannot_run=L["donors_cannot_run"],
                          data_not_checked=L["data_not_checked"]), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
