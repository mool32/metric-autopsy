#!/usr/bin/env python3
"""Level B, S7 (`plan_s7.md`): the list for the owner — the one stop — assembled from the repository's files and the
results branches (fetched beforehand into one directory, one sub-directory per branch), every number read from them.
Nothing here reads a claim's metric.

    python -I validation/level_b/list_s7.py --results DIR --out validation/level_b/amendment1/list \
        [--heads heads.tsv] [--commit SHA]

DIR holds the results branches it reads by their names (`level-b-data-<n>/`, `level-b-controls/`), each with its
SHA256SUMS, checked here. `heads.tsv`: "results/<branch><TAB><commit>", the commits they were taken from; `--commit`:
the repository's commit the list is assembled from. OUT gets list.md (the owner's page) and list.json.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
A1 = HERE / "amendment1"
SET_CAVEATS = ("v1 does not cover a set's path (set_norm_pearson with the coupling_set probe)",
               "the set's definition was not checked against Levy et al. 2020")
# `metrics.set_norm_pearson` (the engine of 5a89f98): NaN where any of the set's pairs is not measured
SET_METRIC_RULE = ("the set's metric is the mean of norm_pearson over all its pairs, not measured (NaN) where any pair "
                   "is not (fewer than min_codetected cells carry both genes) or a gene is not in the data; S5 checked "
                   "the genes' presence, not each pair's co-detection: the runs show it")
# plan_s7.md, section 1, part 8: the journal's open points, as the plan names them (the entries say more)
OPEN_POINTS = (
    ("J7", "the three records Europe PMC refuses (R00014, R00173, R00394), which NCBI's services may hold: none "
           "was asked there, the owner's decision having named Europe PMC"),
    ("J8", "N2's cell type is the Census's \"unknown\" (unannotated cells whose make-up may differ between the two "
           "chemistries; the ranked table holds the next collections), and the positive control gives one pair"),
    ("J15", "R00189's samples 189 and 190, in no row of Martin et al.'s Table S2: not paired by their numbers "
            "(an inference, not the authors' statement); their cells listed, not checked"),
    ("J20", "the bioRxiv preprints refused twice for the requests' rate: another way to their texts (bioRxiv's "
            "own API, another network) is the owner's to decide, not taken"),
)


def _json(path: Path):
    return json.loads(path.read_text())


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------------------------------------------- provenance

def branch_check(d: Path) -> dict:
    """A results branch's copy against its own SHA256SUMS: every file listed present and equal."""
    sums = d / "SHA256SUMS"
    if not sums.exists():
        return dict(sha256sums=None, files=0, bad=["no SHA256SUMS"])
    bad, n = [], 0
    for line in sums.read_text().splitlines():
        if not line.strip():
            continue
        want, name = line.split(maxsplit=1)
        name = name.lstrip("*")
        n += 1
        p = d / name
        if not p.exists() or _sha256(p) != want:
            bad.append(name)
    return dict(sha256sums=_sha256(sums), files=n, bad=bad)


def provenance(results: Path, heads: Path | None, commit: str | None) -> dict:
    h = {}
    if heads and heads.exists():
        for line in heads.read_text().splitlines():
            if line.strip():
                b, c = line.split("\t")
                h[b.removeprefix("results/")] = c.strip()
    branches = {d.name: dict(branch_check(d), commit=h.get(d.name)) for d in sorted(results.glob("level-b-*")) if d.is_dir()}
    return dict(commit=commit, branches=branches)


# ------------------------------------------------------------------------------------------- part 1: the counts

def sources(root: Path = HERE) -> dict:
    """Part 1: the records per source — the search of 2026-10-19 (its corrected list) and amendment 1's additions."""
    s = _json(root / "runs" / "remerge.json")
    pre = [line.split("\t") for line in (root / "amendment1" / "preprints" / "preprints_merge.tsv").read_text().splitlines()[1:]
           if line.strip()]
    status = {}
    for row in pre:
        status[row[1]] = status.get(row[1], 0) + 1
    sets = [x for x in (root / "amendment1" / "set_records.txt").read_text().split() if x]
    with gzip.open(root / "amendment1" / "preprints" / "preprints_new.jsonl.gz", "rt") as f:
        new = sum(1 for line in f if line.strip())
    return dict(search_2026_10_19=dict(pubmed=s["pubmed"]["total"], europepmc=s["europepmc"]["hit_count"],
                                       europepmc_by_source=s["europepmc"]["by_source"],
                                       citations={k: v["listed"] for k, v in s["citations"].items()},
                                       merged=s["merged"], found_by=s["found_by"]),
                amendment_1=dict(sets_records=len(sets), preprints_found=len(pre), preprints_status=status,
                                 preprints_new=new))


def counts(root: Path = HERE) -> dict:
    """Part 1: the screening's counts, both screenings."""
    ta0, ta1 = _json(root / "screening" / "merge.json"), _json(root / "amendment1" / "ta" / "merge.json")
    ft0, ft1 = _json(root / "fulltext_screening" / "ft_merge.json"), _json(root / "amendment1" / "ft" / "merge" / "ft_merge.json")
    pick = lambda d, ks: {k: d[k] for k in ks}  # noqa: E731
    return dict(titles_abstracts={"2026-10-19": pick(ta0, ["n_screened", "kappa", "agreement", "both_keep", "to_owner", "fulltext"]),
                                  "amendment 1": pick(ta1, ["n_screened", "kappa", "agreement", "both_keep", "to_owner", "fulltext"])},
                full_texts={"2026-10-19": pick(ft0, ["read_by_both", "kappa", "agreement", "kept_by_both", "to_owner"]),
                            "amendment 1": pick(ft1, ["read_by_both", "kappa", "agreement", "kept_by_both", "to_owner"])})


# ----------------------------------------------------------------------------------------- part 4: the datasets

def data_reports(results: Path) -> dict:
    """Part 4: each dataset's checks, by spec name (`<rid>` or `<rid>_<x>`), from every `level-b-data-<n>` branch; a
    later branch's report of the same spec (its second attempt) supersedes the earlier one, which is named."""
    out = {}
    for d in sorted(results.glob("level-b-data-*"), key=lambda p: int(p.name.rsplit("-", 1)[1])):
        for rep in sorted(d.glob("*/report.json")):
            r = _json(rep)
            r["_branch"] = d.name
            if rep.parent.name in out:
                r["_supersedes"] = out[rep.parent.name]["_branch"]
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


def dataset_status(rep: dict, claim_id: str) -> dict:
    """One claim on one dataset: the check of origin, the source, rule 4, and its genes (missing: the metric cannot
    be measured on the claim's genes as listed)."""
    g = rep.get("genes", {}).get(claim_id, {})
    return dict(origin=origin(rep), source=rep.get("source", {}).get("verdict"), rule4=not donors_cannot_run(rep),
                missing=g.get("missing", []), present=len(g.get("present", [])))


def specs(root: Path = HERE) -> dict:
    return {p.stem: _json(p) for p in sorted((root / "data_specs").glob("*.json"))}


# ------------------------------------------------------------------------------------------ part 5: the verdicts

def verdicts(root: Path = HERE) -> dict | None:
    """Part 5 (S4): per claim the records found, by source, kept by the filter and coded; the verdicts both passes
    code alike, with their quotes; the lines coded differently, for the owner. None before the passes' comparison."""
    m, b = root / "verdicts" / "merge", root / "verdicts" / "batches"
    if not (m / "merge.json").exists():
        return None
    rows = lambda p: [dict(zip(p.read_text().splitlines()[0].split("\t"), line.split("\t")))  # noqa: E731
                      for line in p.read_text().splitlines()[1:] if line.strip()]
    res = _json(m / "merge.json")
    manifest = _json(b / "manifest.json")
    return dict(merge={k: v for k, v in res.items() if k != "claims"}, claims=res["claims"],
                agreed=rows(m / "verdicts.tsv"), to_owner=rows(m / "to_owner.tsv"),
                batches=manifest["batches"], cap=manifest.get("cap"))


# ------------------------------------------------------------------------------------------ part 6: the controls

def controls(results: Path) -> dict | None:
    p = results / "level-b-controls" / "claims.json"
    return _json(p) if p.exists() else None


def _tsv(path: Path) -> list:
    lines = path.read_text().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"))) for line in lines[1:] if line.strip()]


def controls_detail(results: Path, kept: dict) -> dict | None:
    """Part 6's detail: each negative control's donors and cells per condition (from the selection's units), what the
    selection left out and why, and plan_s6's check that no control's collection is a kept paper's (by DOI and
    title)."""
    d = results / "level-b-controls"
    if not (d / "claims.json").exists():
        return None
    sel = _json(d / "selection.json")
    cand = _tsv(d / "candidates.tsv")
    tried = _tsv(d / "tried.tsv")
    doi = {c["collection_id"]: c.get("collection_doi", "") for c in cand}
    excluded = {}
    for c in cand:
        if c.get("excluded"):
            excluded[c["excluded"]] = excluded.get(c["excluded"], 0) + 1
    units = {}
    for ch in sel["negative"]["chosen"]:
        per = {}
        for u in ch["units"]:
            donor, cond = [x.strip() for x in u["unit"].split("|")]
            per.setdefault(cond, []).append(u["cells"])
        units[ch["id"]] = {c: dict(donors=len(v), cells=sum(v), cells_each=sorted(set(v))) for c, v in per.items()}
    norm = lambda t: re.sub(r"[^a-z0-9]+", " ", (t or "").lower()).strip()  # noqa: E731
    kept_dois = {p["doi"].lower(): p["rid"] for p in kept["papers"] if p.get("doi")}
    kept_titles = {norm(p["title"]): p["rid"] for p in kept["papers"]}
    same = []
    for ch in sel["negative"]["chosen"]:
        x = doi.get(ch["collection_id"], "").lower()
        if x and x in kept_dois:
            same.append(f"{ch['id']}: DOI {x} is {kept_dois[x]}'s")
        if norm(ch["collection_name"]) in kept_titles:
            same.append(f"{ch['id']}: its collection's title is {kept_titles[norm(ch['collection_name'])]}'s")
    pos = sel.get("positive", {})
    for k, v in sel["rules"].get("excluded_dois", {}).items():
        if "positive control" in v and k.lower() in kept_dois:
            same.append(f"the positive control: DOI {k} is {kept_dois[k.lower()]}'s")
    return dict(release=sel.get("release"), rules=sel["rules"], candidates=sel["negative"]["candidates"],
                excluded=sel["negative"]["excluded"], excluded_by=excluded, collections=sel["negative"]["collections"],
                tried=tried, units=units, dois={ch["id"]: doi.get(ch["collection_id"], "") for ch in sel["negative"]["chosen"]},
                same_as_kept=same, positive_source=pos.get("source"))


# --------------------------------------------------------------------------------------- part 8: for the owner

def ft_to_owner(root: Path = HERE) -> dict:
    """Amendment 1's records sent to the owner at the full texts (`ft_to_owner.txt`): those both passes could not read
    (no full text), and those they decide differently or mark unclear, with each pass's words."""
    blocks = [b for b in (root / "amendment1" / "ft" / "merge" / "ft_to_owner.txt").read_text().split("\n\n") if b.strip()]
    without, differ = [], []
    for b in blocks:
        lines = [x for x in b.splitlines() if x.strip()]
        rid, _, title = lines[0].partition("  ")
        passes, claims = {}, {}
        for x in lines[1:]:
            m = re.match(r"pass (\d)( claim \d+)?: (.*)$", x.strip())
            if m and m.group(2):
                claims.setdefault(m.group(1), []).append(f"claim{m.group(2)[6:]}: {m.group(3)}")
            elif m:
                passes[m.group(1)] = m.group(3)
        rec = dict(rid=rid.strip(), title=title.strip(), passes=passes, claims=claims)
        (without if all(v.startswith("not read (no full text)") for v in passes.values()) else differ).append(rec)
    return dict(without_fulltext=without, decided_differently=differ)


def one_pass_claims(kept: dict) -> list:
    """The claims one pass lists alone (`plan_s3.md`, section 4), each with the other pass's claims of the same paper
    that hold the same genes."""
    out = []
    for p in kept["papers"]:
        for c in p["claims"]:
            if len(c["passes"]) != 1:
                continue
            twins = [f"{o['id']} (passes {', '.join(map(str, o['passes']))})" for o in p["claims"]
                     if o is not c and c["genes"] and set(o["genes"]) == set(c["genes"])]
            out.append(dict(claim=f"{p['rid']}:{c['id']}", kind=c["kind"], pass_=c["passes"][0], genes=len(c["genes"]),
                            name=c.get("set_name") or ", ".join(c["genes"]), same_genes_as=twins))
    return out


def journal(root: Path = HERE) -> list:
    """Every entry of the journal (since the owner's decisions of 2026-11-10, J1), with whether its text records a
    deviation in that word."""
    text = (root / "JOURNAL.md").read_text()
    out = []
    for m in re.finditer(r"^## (J\d+) — (\S+): (.*)$", text, flags=re.M):
        end = text.find("\n## J", m.end())
        body = text[m.end(): end if end > 0 else len(text)]
        out.append(dict(entry=m.group(1), date=m.group(2), title=m.group(3).strip(),
                        deviation=bool(re.search(r"(?<!no )deviation", (m.group(3) + body).lower()))))
    return out


# ------------------------------------------------------------------------------------------------- the assembly

def assemble(results: Path, root: Path = HERE, heads: Path | None = None, commit: str | None = None) -> dict:
    kept = _json(root / "amendment1" / "kept_claims.json")
    reps = data_reports(results)
    sp = specs(root)
    papers, claims = [], []
    for p in kept["papers"]:
        mine = {k: v for k, v in reps.items() if k == p["rid"] or k.startswith(p["rid"] + "_")}
        data = {k: dict(report=v, spec=sp.get(k, {}), origin=origin(v), cannot_run=donors_cannot_run(v),
                        donors_rule4=not donors_cannot_run(v),
                        claims=[str(c["id"]) for c in sp.get(k, {}).get("claims", [])]) for k, v in mine.items()}
        for c in p["claims"]:
            on = {k: dataset_status(d["report"], str(c["id"])) for k, d in data.items() if str(c["id"]) in d["claims"]}
            claims.append(dict(claim=f"{p['rid']}:{c['id']}", rid=p["rid"], kind=c["kind"], datasets=on))
        papers.append(dict(rid=p["rid"], title=p["title"], read=p["read"], pmid=p["pmid"], doi=p["doi"],
                           claims=p["claims"], data=data))

    def tally(test) -> dict:
        every = [c["claim"] for c in claims if c["datasets"] and all(test(s) for s in c["datasets"].values())]
        some = [c["claim"] for c in claims if c["datasets"] and any(test(s) for s in c["datasets"].values())
                and c["claim"] not in every]
        return dict(every=every, some=some, papers_every=sorted({x.split(":")[0] for x in every}))

    def counted(test) -> dict:
        t = tally(test)
        return dict(t, papers=t["papers_every"], claims=len(t["every"]))

    unchecked = [c["claim"] for c in claims if not c["datasets"]]
    no_data = sorted({c["rid"] for c in claims if all(not x["datasets"] for x in claims if x["rid"] == c["rid"])})
    not_listed = [f"{p['rid']}:{c['id']}" for p in kept["papers"] for c in p["claims"] if not c["genes"]]
    structure_only = [f"{name} ({d['spec'].get('accession', '?')})" for p in papers for name, d in p["data"].items()
                      if not d["claims"]]
    return dict(provenance=provenance(results, heads, commit), sources=sources(root), counts=counts(root),
                kept_counts=kept["counts"], fewer_than_target=kept["fewer_than_target"], drawn=kept["drawn"],
                papers=papers, claims=claims, verdicts=verdicts(root), controls=controls(results),
                controls_detail=controls_detail(results, kept),
                origin_fails=counted(lambda s: s["origin"]["fails"]),
                donors_cannot_run=counted(lambda s: not s["rule4"]),
                genes_missing=counted(lambda s: bool(s["missing"])),
                source_not_named=counted(lambda s: s["source"] == "no source named"),
                source_partly=counted(lambda s: s["source"] == "partly"),
                data_not_checked=dict(claims_list=unchecked, claims=len(unchecked), papers=no_data,
                                      genes_not_listed=not_listed, datasets_without_claim=structure_only),
                ft_to_owner=ft_to_owner(root), one_pass=one_pass_claims(kept), journal=journal(root))


# ------------------------------------------------------------------------------------------------- the page

def _kappa(x) -> str:
    """A kappa, or why it is undefined: both passes code every line in one and the same category (chance agreement 1)."""
    return "undefined (one category in both passes)" if x is None or x != x else f"{x:.3f}"


def _n(x) -> str:
    return f"{x:,}" if isinstance(x, int) else str(x)


def _line(ln: dict) -> str:
    return (f"{ln.get('species')}; {ln.get('cells')}; {ln.get('conditions')}; {ln.get('direction')}; metric "
            f"{ln.get('metric')}; data {ln.get('data')}; donors {ln.get('donors')}; rule 4 {ln.get('rule4')} — "
            f"\"{ln.get('quote')}\" ({ln.get('where')})")


def _dataset_md(name: str, d: dict) -> list:
    r, s = d["report"], d["spec"]
    dn = r.get("donors", {})
    r4 = dn.get("rule4", {})
    src = r.get("source", {})
    sex = r.get("sex", {})
    out = [f"- **`{name}`** ({s.get('accession', '?')}; `results/{r['_branch']}`"
           + (f", superseding `results/{r['_supersedes']}`'s first attempt" if r.get("_supersedes") else "")
           + f"; claims {', '.join(d['claims']) or 'none'})"]
    for f in r.get("files", []):
        out.append(f"  - file `{f.get('file')}`: {_n(f.get('bytes'))} bytes, SHA256 `{f.get('sha256')}` ({f.get('url')})")
    cnt = r.get("counts", {})
    out.append(f"  - 1 counts: non-negative {cnt.get('nonnegative')}, integers {cnt.get('integers')}, raw "
               f"**{cnt.get('raw')}**, max {cnt.get('max')}; cells {_n(r.get('cells'))} ("
               + ", ".join(f"{k} {_n(v)}" for k, v in (r.get('cells_per_group') or {}).items()) + ")")
    p_min = f"{r4['min_p']:.4g}" if isinstance(r4.get("min_p"), (int, float)) else "-"
    out.append(f"  - 2 donors: {dn.get('design')}, " + ", ".join(f"{k} {v}" for k, v in (dn.get("donors") or {}).items())
               + f"; ids naming several: {dn.get('composite') or 'none'}; rule 4 **{'passes' if r4.get('passes') else 'fails'}** "
               f"({r4.get('relabelings')} relabelings, smallest p {p_min})")
    out.append(f"  - 3 one individual per id: both sexes {sex.get('both_sexes')}; fails {sex.get('fails') or 'none'}")
    genes = r.get("genes", {})
    out.append("  - 4 genes: " + ("; ".join(f"claim {k} {len(v.get('present', []))} present"
                                           + (f", **missing {', '.join(v['missing'])}**" if v.get("missing") else "")
                                           for k, v in genes.items()) or "no claim's genes checked"))
    shares = src.get("shares")
    out.append(f"  - 5 source: **{src.get('verdict')}**" + (f" ({', '.join(map(str, src.get('sources') or []))})" if src.get("sources") else "")
               + (f"; shares {json.dumps(shares, ensure_ascii=False)}" if shares and src.get("verdict") == "partly" else ""))
    out.append(f"  - 6 the design as the paper ran it: {s.get('design', '-')}")
    if r.get("without_donor"):
        out.append("  - cells left out, no donor named: " + "; ".join(f"{k} {_n(v)}" for k, v in r["without_donor"].items()))
    if r.get("donors_left_out"):
        out.append("  - donors left out (cells in one group only, or too few): "
                   + "; ".join(f"{k} {json.dumps(v)}" for k, v in r["donors_left_out"].items()))
    o = d["origin"]
    out.append(f"  - check of origin (12.3): **{'FAILS: ' + '; '.join(o['why']) if o['fails'] else 'passes'}**")
    if s.get("notes"):
        out.append(f"  - the spec's notes: {s['notes']}")
    return out


def flags(d: dict) -> list:
    """What one dataset's checks leave for the owner, in a few words each."""
    r, out = d["report"], []
    if d["origin"]["fails"]:
        out.append("origin FAILS: " + "; ".join(d["origin"]["why"]))
    src = r.get("source", {}).get("verdict")
    if src in ("no source named", "partly"):
        out.append(f"source {src}")
    miss = {k: v["missing"] for k, v in r.get("genes", {}).items() if v.get("missing")}
    if miss:
        out.append("genes missing: " + "; ".join(f"claim {k} {', '.join(v)}" for k, v in miss.items()))
    if not d["donors_rule4"]:
        out.append("rule 4 fails")
    if r.get("without_donor"):
        out.append(f"{sum(r['without_donor'].values()):,} cells without a donor, left out")
    if r.get("donors_left_out"):
        out.append(f"donors left out: {', '.join(r['donors_left_out'])}")
    if r.get("_supersedes"):
        out.append(f"a second attempt (the first on {r['_supersedes']})")
    return out or ["none"]


def markdown(L: dict) -> str:
    out = ["# Level B: the list for the owner (S7; `plan_s7.md`)", "",
           "**Nothing here is validated; no claim's metric was read.** Every number below is read by "
           "`list_s7.py` from the repository's files and the results branches. The runs on the engine of 5a89f98 start "
           "only after the owner approves this list; what the owner decides is gathered in part 8.", ""]
    pv = L.get("provenance") or {}
    if pv:
        out.append(f"Assembled from the repository at `{pv.get('commit') or '?'}` and:")
        for b, x in (pv.get("branches") or {}).items():
            out.append(f"- `results/{b}` at `{x.get('commit') or '?'}`: SHA256SUMS `{x.get('sha256sums')}`, "
                       f"{x.get('files')} files, " + ("all as listed" if not x.get("bad") else f"**not as listed: {x['bad']}**"))
        out.append("")
    out += ["## 1. The search and the screening", ""]
    src = L.get("sources")
    if src:
        s, a = src["search_2026_10_19"], src["amendment_1"]
        out.append(f"- the search of 2026-10-19: PubMed {_n(s['pubmed'])}; Europe PMC {_n(s['europepmc'])} "
                   f"({', '.join(f'{k} {v}' for k, v in s['europepmc_by_source'].items())}); the method papers' citations "
                   + ", ".join(f"{k} {v}" for k, v in s["citations"].items())
                   + f"; merged {_n(s['merged'])} records (the corrected list), found by "
                   + ", ".join(f"{k} {v}" for k, v in s["found_by"].items()))
        out.append(f"- amendment 1: the sets' records {a['sets_records']}; preprints found {a['preprints_found']} ("
                   + ", ".join(f"{k} {v}" for k, v in a["preprints_status"].items()) + f"), new {a['preprints_new']}")
    c = L["counts"]
    for step, rows in (("titles and abstracts", c["titles_abstracts"]), ("full texts", c["full_texts"])):
        for when, r in rows.items():
            out.append(f"- {step}, {when}: " + "; ".join(f"{k} {v:.3f}" if isinstance(v, float) else f"{k} {v}"
                                                        for k, v in r.items()))
    k = L["kept_counts"]
    out += ["", f"## 2. The kept papers: {k['papers']} ({k['papers_with_pairs']} with pairs, {k['papers_with_sets']} with sets)", "",
            f"{k['claims']} claims ({k['pairs']} pairs, {k['sets']} sets; {k['both_passes']} listed by both passes, "
            f"{k['one_pass']} by one); papers read in amendment 1 {k['read_in_amendment_1']}, kept on 2026-10-19 "
            f"{k['kept_on_2026_10_19']}; the draw of 30 {'made' if L['drawn'] else 'not needed'}."]
    if L["fewer_than_target"]:
        out.append(f"**Fewer than 20 papers with pairs** (the protocol's section 9): {k['papers_with_pairs']}.")
    out += ["", "## 3–4. Each paper's claims and data", ""]
    for p in L["papers"]:
        out += [f"### {p['rid']}: {p['title']}", f"PMID {p['pmid'] or '-'}, DOI {p['doi'] or '-'}; read {p['read']}", ""]
        for kind in ("pair", "set"):
            for cl in [x for x in p["claims"] if x["kind"] == kind]:
                if kind == "pair":
                    what = ", ".join(cl["genes"])
                else:
                    what = (f"{cl.get('set_name')}: {len(cl['genes'])} genes taken"
                            + (f" (the paper's number {cl['set_size']})" if cl.get("set_size") else "")
                            + f"; where listed: {cl.get('set_where')}; genes from {cl.get('genes_from')}")
                out.append(f"- **{kind} {cl['id']}: {what}** — passes {', '.join(map(str, cl['passes']))}")
                for ps, ln in sorted(cl.get("lines", {}).items()):
                    out.append(f"  - pass {ps}: {_line(ln)}")
                if kind == "set":
                    out.append("  - caveats: " + "; ".join(SET_CAVEATS))
        out.append("")
        for name, d in p["data"].items():
            out += _dataset_md(name, d)
        if not p["data"]:
            out.append("- data: not checked")
        out.append("")
    out += ["## 5. The independent verdicts (S4)", ""]
    v = L.get("verdicts")
    if v:
        mg = v["merge"]
        out.append(f"- coded by both passes: {mg['lines']} lines; the decisions: agreement {mg['agreement_decision']:.3f}, "
                   f"kappa {_kappa(mg['kappa_decision'])}; the full codes: agreement {mg['agreement_code']:.3f}, kappa "
                   f"{_kappa(mg['kappa_code'])}; agreed verdicts {mg['verdicts_agreed']}; to the owner {mg['to_owner']}")
        for c in v["claims"]:
            capped = f", read by its paper's terms {c.get('read')} (over the cap)" if c.get("terms_applied") else ""
            out.append(f"- {c['claim']} ({c['kind']}): found {c['records']} (by citations {c['by_citations']}, by the "
                       f"terms {c['by_terms']}), kept by the filter {c['kept']}{capped}, coded {c['coded']}; agreed "
                       f"verdicts {c['verdicts_agreed']} {c['agreed_codes'] or ''}; to the owner {c['to_owner']}")
            for a in [x for x in v["agreed"] if x["claim"] == c["claim"]]:
                out.append(f"  - {a['record']}: {a['kind']}, {a['code']} — \"{a['quote_1']}\" ({a['where_1']}); "
                           f"pass 2: \"{a['quote_2']}\" ({a['where_2']})")
        for t in v["to_owner"]:
            out.append(f"- coded differently: {t.get('claim')} {t.get('record')}: pass 1 {t.get('decision_1')} "
                       f"{t.get('kind_1', '')} {t.get('code_1', '')}; pass 2 {t.get('decision_2', '')} {t.get('kind_2', '')} "
                       f"{t.get('code_2', '')}")
        out.append("")
    else:
        out += ["- the verdicts not coded yet", ""]
    out += ["## 6. The controls", ""]
    ctl, cd = L["controls"], L.get("controls_detail")
    if ctl:
        out.append(f"The Census release {ctl.get('release')}.")
        for n in ctl["negative"]:
            u = (cd or {}).get("units", {}).get(n["id"], {})
            per = "; ".join(f"{cond}: {x['donors']} donors, {_n(x['cells'])} cells ({'/'.join(map(str, x['cells_each']))} each)"
                            for cond, x in u.items())
            out.append(f"- **{n['id']}**: \"{n['collection_name']}\" (collection {n['collection_id']}, DOI "
                       f"{(cd or {}).get('dois', {}).get(n['id']) or '-'}); {n['tissue']}, {n['cell_type']}; A {n['condition_a']} "
                       f"against B {n['condition_b']}; {n['donors']} donors, {_n(n['cells'])} cells ({per}); sex calls "
                       f"against labels: {n.get('sex_conflicts') or 'none'}; extraction SHA256 `{n.get('npz_sha256')}`")
            levels = {}
            for c in n["claims"]:
                levels.setdefault((c["kind"], c.get("level")), []).append(c)
            for (kind, lev), cs in levels.items():
                if kind == "pair":
                    out.append(f"  - pairs, {lev}: " + "; ".join(f"{'–'.join(c['genes'])} ({c['higher_in']})" for c in cs))
                else:
                    out.append(f"  - sets, {lev}: " + "; ".join(f"{c['id']} {len(c['genes'])} genes ({c['higher_in']})" for c in cs))
            pc = n.get("pool", {}).get("positive_control", {})
            out.append(f"  - the pool: {n.get('pool', {}).get('pairs_per_level')} pairs per level; GATE 5's positive pair "
                       f"{'–'.join(pc.get('pair', []))} (r {pc.get('r', float('nan')):.3f})")
        for p in ctl["positive"]:
            out.append(f"- **positive**: {p['source']}; {p['type']}; {_n(p['cells'])} cells, {_n(p['singlets'])} singlets, "
                       f"integers {p['integers']}; donors in both conditions {len(p['donors_in_both'])}; label problems "
                       f"{p['label_problems'] or 'none'}; sex conflicts {p['sex_conflicts'] or 'none'}; cells kept "
                       f"{_n(p['cells_kept'])}; extraction SHA256 `{p.get('npz_sha256')}`")
            out.append("  - genes: " + "; ".join(f"{g['gene']} {g['detection_ctrl']:.1%} → {g['detection_stim']:.1%}, induced in "
                                                 f"{g['induced_in']}/{g['donors']}{' (eligible)' if g['eligible'] else ''}"
                                                 for g in p["genes"] if g.get("present")))
            out.append("  - claims: " + "; ".join(f"{c['id']} {'–'.join(c['genes'])} higher in {c['higher_in']}" for c in p["claims"]))
        out.append(f"- positive controls failed: {ctl['positive_failed'] or 'none'}")
        if cd:
            out.append(f"- left out by the selection: {cd['candidates']} qualifying candidates, {cd['excluded']} excluded ("
                       + ", ".join(f"{k} {v}" for k, v in cd["excluded_by"].items()) + f"), {cd['collections']} other "
                       f"collections, " + "; ".join(f"rank {t['rank']} {t['outcome']}" for t in cd["tried"]))
            out.append("- **no control's collection is a kept paper's** (DOIs and titles compared)" if not cd["same_as_kept"]
                       else f"- **a control's collection is a kept paper's: {cd['same_as_kept']}**")
    else:
        out.append("- the controls' branch not given")
    tl = lambda t: (f"{len(t['every'])} claims on every dataset that checks them ({', '.join(t['every']) or 'none'}); "  # noqa: E731
                    f"{len(t['some'])} on some of them ({', '.join(t['some']) or 'none'})")
    dn = L["data_not_checked"]
    out += ["", "## 7. 12.3's two counts so far", "",
            f"- **(i) fail the check of origin**: {tl(L['origin_fails'])}",
            "- **(ii) the test over donors cannot run:**",
            f"  - by rule 4: {tl(L['donors_cannot_run'])}",
            f"  - the metric not measured on the genes as listed (a gene not in the data): {tl(L['genes_missing'])}",
            f"  - the metric not measured, the paper listing no genes: {len(dn.get('genes_not_listed', []))} claims "
            f"({', '.join(dn.get('genes_not_listed', [])) or 'none'})",
            f"  - no dataset checks the claim's genes: {dn['claims']} claims ({', '.join(dn['claims_list']) or 'none'}); "
            f"datasets checked for their structure without a claim's genes: "
            f"{', '.join(dn.get('datasets_without_claim', [])) or 'none'} (the spec's notes say why)",
            f"- the check of the source cannot be made (no source named; not a pass): {tl(L['source_not_named'])}",
            f"- the source partly shared (listed with its shares; not a failure): {tl(L['source_partly'])}",
            "", "Whether the test over donors supports a claim is known only after the runs.", ""]
    out += ["## 8. For the owner (nothing here is decided)", ""]
    fo = L.get("ft_to_owner")
    if fo:
        out.append(f"**8.1 The full texts' records of amendment 1** (`amendment1/ft/merge/ft_to_owner.txt`): "
                   f"{len(fo['decided_differently'])} decided differently or marked unclear, {len(fo['without_fulltext'])} "
                   f"without a full text (excluded by default, J1's rule; why each: J18, J20).")
        for r in fo["decided_differently"]:
            out.append(f"- {r['rid']} {r['title']}")
            for ps, t in sorted(r["passes"].items()):
                out.append(f"  - pass {ps}: {t}")
                out += [f"    - {c}" for c in r.get("claims", {}).get(ps, [])]
        out.append("- without a full text:")
        out += [f"  - {r['rid']} {r['title']}" for r in fo["without_fulltext"]]
        out.append("")
    op = L.get("one_pass") or []
    out.append(f"**8.2 The claims one pass lists alone:** {len(op)}; pairs whose two symbols cannot be told: "
               f"{L['kept_counts'].get('ambiguous', 0)}.")
    for c in op:
        out.append(f"- {c['claim']} ({c['kind']}, pass {c['pass_']}): {c['name']}"
                   + (f" ({c['genes']} genes)" if c["kind"] == "set" else "")
                   + (f" — the same genes as claim {', '.join(c['same_genes_as'])}" if c["same_genes_as"] else ""))
    out.append("")
    if v:
        mg = v["merge"]
        none = [c["claim"] for c in v["claims"] if not c["kept"]]
        over = [c for c in v["claims"] if c.get("terms_applied")]
        out.append(f"**8.3 The verdicts:** coded differently {mg['to_owner']}; agreed {mg['verdicts_agreed']}. No record "
                   f"kept for {len(none)} claims ({', '.join(none)}). Over the cap of {v.get('cap')}: "
                   + ("; ".join(f"{c['claim']} {c['kept']} kept, {c['read']} read by its paper's terms, "
                                f"{c['kept'] - c['read']} not read" for c in over) or "none")
                   + ". The independence check as the plan states it is not met: the harness's task list was in every "
                   "coder's context (J31, a deviation).")
        out.append("")
    out.append("**8.4 The data**, per dataset: the paper's donors (pass 1's words for the dataset's first claim) beside "
               "the deposit's, and what parts 4 and 7 flag; where the authors' cells are not identifiable from the deposit, "
               "the spec's notes (part 4) and the journal (J19, J27, J29) say so. For the sets: " + SET_METRIC_RULE + ".")
    out.append("")
    out.append("| Dataset | Claims | The paper's donors (pass 1) | The deposit's donors | Flags |")
    out.append("|---|---|---|---|---|")
    for p in L["papers"]:
        for name, d in p["data"].items():
            r = d["report"]
            first = next((c for c in p["claims"] if str(c["id"]) in d["claims"]), p["claims"][0] if p["claims"] else None)
            said = (first or {}).get("lines", {}).get("1", (first or {}).get("lines", {}).get("2", {})).get("donors", "-")
            dep = ", ".join(f"{k} {v}" for k, v in (r.get("donors", {}).get("donors") or {}).items())
            cell = lambda t: str(t).replace("|", "\\|")  # noqa: E731
            out.append(f"| `{name}` | {', '.join(d['claims']) or 'none'} | {cell(said)} | {r.get('donors', {}).get('design')}: "
                       f"{cell(dep)} | " + cell("; ".join(flags(d))) + " |")
    out.append("")
    if ctl:
        out.append(f"**8.5 The controls:** {len(ctl['negative'])} negative controls (cell types: "
                   + ", ".join(f"{n['id']} {n['cell_type']}" for n in ctl["negative"])
                   + f"); positive claims {sum(len(p['claims']) for p in ctl['positive'])} (J8).")
        out.append("")
    out.append("**8.6 The journal's open points** (`plan_s7.md`, part 8):")
    for j, t in OPEN_POINTS:
        out.append(f"- {j}: {t}")
    out.append("")
    out.append("**8.7 Every journal entry since the owner's decisions of 2026-11-10** (\"deviation\": its text records one "
               "in that word; the others record runs, technical failures and their second attempts, and precisions):")
    for e in L.get("journal") or []:
        out.append(f"- {e['entry']} ({e['date']}){' **deviation**' if e['deviation'] else ''}: {e['title']}")
    return "\n".join(out) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--heads", type=Path, help="results/<branch><TAB><commit> per line")
    ap.add_argument("--commit", help="the repository's commit the list is assembled from")
    a = ap.parse_args(argv)
    L = assemble(a.results, heads=a.heads, commit=a.commit)
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "list.json").write_text(json.dumps(L, indent=1, ensure_ascii=False, default=str) + "\n")
    (a.out / "list.md").write_text(markdown(L))
    bad = {b: x["bad"] for b, x in L["provenance"]["branches"].items() if x["bad"]}
    print(json.dumps(dict(papers=len(L["papers"]), claims=len(L["claims"]), branches_not_as_listed=bad,
                          origin_fails=L["origin_fails"], donors_cannot_run=L["donors_cannot_run"],
                          genes_missing=L["genes_missing"], source_not_named=L["source_not_named"],
                          data_not_checked=L["data_not_checked"]), indent=1))
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
