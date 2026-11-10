"""The owner's sheet of level B's undecided records (the owner, 2026-11-10: «Сведи 12 и 19 записей в лист
решений на одну страницу»).

The 12 records sent to the owner at titles and abstracts (`screening/to_owner.tsv`) and the 19 at the full
texts (`fulltext_screening/ft_to_owner.txt`) are 23 records, 8 of them in both lists. They fall into four
groups by what keeps them undecided: no full text; read at the full text by one pass only; unclear in both
at the full text; unclear at titles and abstracts and so read by no pass at the full text (a pass read the
full texts of the records it kept). For each: the passes' decisions as those files give them, whether a
full text was fetched (`runs/fulltext_index.tsv`) and a one-line summary of the passes' notes, written by
hand (`SUMMARY`; the notes themselves are in the files named). Nothing is decided here.

    python validation/level_b/decision_sheet.py > validation/level_b/decision_sheet.md
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

# One line per record: what the passes' notes say, shortened by hand (the notes are in the input files).
SUMMARY = {
    # no full text (both passes kept it at titles and abstracts)
    "R00014": "gene-gene correlation across macrophages changes between cytokine environments; pairs may be in the results",
    "R00173": "HSPC coexpression patterns in GATA2 deficiency against healthy donors; pairs may be in the results",
    "R00394": "the same donors before and after BCG: modules with CCL3 and CCL4 lose coordination (a set; pairs possible)",
    "R00438": "differentially co-expressed ligand-receptor pairs, HCC against normal tissue; bulk and scRNA; the measure open",
    "R00582": "MYBL2-ATAD2 co-expression in ovarian cancer, bulk and single-cell; a change between conditions open",
    "R01062": "FCRL3 correlates with TIGIT in Tregs of SLE patients; a comparison with controls possible",
    "R01338": "differential correlation of cuproptosis genes across MAFLD stages; whether single-cell, open",
    "R02068": "ENO1-GPX4 co-expression stronger in CML with a poor response to TKIs",
    "R02134": "CD58-CD147 (and CD63) co-expression up in pSS T cells; may rest on flow cytometry",
    # read at the full text by one pass only, which excluded it
    "R00220": "Drosha knockout against wild-type mESC: claims of sets or global measures only; one cell line, no donors",
    "R01910": "Spp1-Vegfa 'co-expression' in relapsed tumours is double-positive cells, not a correlation",
    "R02061": "the pair shifts are correlations across bulk samples; single-cell data used for expression only",
    "R02262": "LMNA-PRKDC 'co-expression' is the share of double-positive cells; one PDX line",
    "R02780": "the pairs are GRN edges (pySCENIC); the rest regulon activity and network rewiring (sets), 687 donors",
    "R02787": "lncRNA-neighbour correlations across pseudobulk samples; no change between conditions claimed",
    "R02831": "METTL3-PRMT5 correlation across bulk samples, IHC or cell lines; in scRNA-seq only UMAP plots",
    "R02858": "Oxt-Avp correlation higher in subordinate mice, 3 against 3 per sex (smallest p 0.1)",
    # unclear in both at the full text (rule 3: the data public as counts)
    "R03285": "LGALS9-HAVCR2 and LGALS9-TIGIT correlations higher in TB non-responders, 4 against 4; data promised at NGDC",
    "R03433": "APOE-PPARG correlation in astrocytes, 10 against 7 donors; the text names no repository",
    # unclear at titles and abstracts, so read by no pass at the full text
    "R00222": "no abstract; by its title a commentary on Levy et al. 2020 (X1 if so)",
    "R00599": "a letter without an abstract: COVID-19 multi-organ co-expression; new data and a pair's change unknown",
    "R01382": "no abstract: MAPK-pathway co-expression in lung adenocarcinoma; a pair's change unknown",
    "R01395": "a retracted paper; pass 2: no change of coordination (X4)",
}

GROUPS = [
    ("A", "No full text (both passes kept it at titles and abstracts)",
     "ask Europe PMC again for the 3 with a PMC id (it answered HTTP 500 to all 6 attempts) and look for the 6 "
     "others elsewhere (the publisher, the authors), then both passes read them; or exclude them, as without an "
     "accessible full text."),
    ("B", "Read at the full text by one pass only, which excluded it (the passes decided it differently at titles and abstracts)",
     "accept the reading pass's exclusion; or the other pass reads the full text. Under the amendment (section 12.1) "
     "R02780's network claims are read again by the passes over sets anyway; R00220 and R02858 fail rule 4 as sets too."),
    ("C", "Unclear in both at the full text (rule 3: the data public as counts)",
     "keep it if its data are public as counts before the list is approved, else exclude; or exclude it now."),
    ("D", "Unclear at titles and abstracts, so read by no pass at the full text",
     "both passes read the 3 fetched full texts (R00222 has none; excluded if it is a commentary); or exclude the 4."),
]


def read_lists():
    ta = {r["rid"]: r for r in csv.DictReader(open(HERE / "screening" / "to_owner.tsv"), delimiter="\t")}
    text = (HERE / "fulltext_screening" / "ft_to_owner.txt").read_text()
    ft = {}
    for block in re.split(r"\n\s*\n", text.strip()):
        lines = block.strip().splitlines()
        rid, title = lines[0].split(None, 1)
        passes = {}
        for line in lines[1:]:
            m = re.match(r"\s*pass (\d): (.*)", line)
            if m:
                passes[int(m.group(1))] = m.group(2).strip()
        ft[rid] = dict(title=title.strip(), passes=passes)
    index = {r["rid"]: r for r in csv.DictReader(open(HERE / "runs" / "fulltext_index.tsv"), delimiter="\t")}
    allrows = {r["rid"]: r for r in csv.DictReader(open(HERE / "screening" / "all.tsv"), delimiter="\t")}
    return ta, ft, index, allrows


def group_of(rid, ta, ft) -> str:
    if rid not in ft:
        return "D"
    p = ft[rid]["passes"]
    texts = [p.get(1, ""), p.get(2, "")]
    if all(t.startswith("not read (no full text)") for t in texts):
        return "A"
    if any(t.startswith("not read (titles and abstracts") for t in texts):
        return "B"
    if all(t.startswith("UNCLEAR") for t in texts):
        return "C"
    raise ValueError(f"{rid}: a record outside the four groups: {texts}")


def short(title: str, n: int = 58) -> str:
    title = title.rstrip(".")
    return title if len(title) <= n else title[: n - 1].rstrip() + "…"


def code_of(rid, ta, ft, allrows) -> str:
    """The passes' decisions, compact: at titles and abstracts, then at the full text where read."""
    r = ta.get(rid) or allrows[rid]
    at = f"{r['decision 1'][:4]}{'/' + r['code 1'] if r['code 1'] not in ('-', '') else ''}" \
         f" · {r['decision 2'][:4]}{'/' + r['code 2'] if r['code 2'] not in ('-', '') else ''}"
    out = f"TA {at}"
    if rid in ft:
        reads = []
        for k in (1, 2):
            t = ft[rid]["passes"].get(k, "")
            m = re.match(r"(EXCLUDE|KEEP|UNCLEAR) (F\d+|-)", t)
            if m:
                reads.append(f"p{k} {m.group(1)[:4]}{'/' + m.group(2) if m.group(2) != '-' else ''}")
        if reads:
            out += "; FT " + ", ".join(reads)
    return out


def fulltext_of(rid, index) -> str:
    row = index.get(rid)
    if row is None:
        return "not asked"
    if row["fulltext"] not in ("-", ""):
        return "fetched"
    if row["tried"] in ("-", ""):
        return "no PMC id"
    if "HTTP 500 on all 6 attempts" in row["note"]:
        return "HTTP 500 ×6"
    return f"not fetched: {row['note']}"


def main() -> int:
    ta, ft, index, allrows = read_lists()
    union = sorted(set(ta) | set(ft))
    both = sorted(set(ta) & set(ft))
    assert len(ta) == 12 and len(ft) == 19, (len(ta), len(ft))
    assert set(SUMMARY) == set(union), set(SUMMARY) ^ set(union)
    groups = {g: [r for r in union if group_of(r, ta, ft) == g] for g, _, _ in GROUPS}
    title = {r: (ft[r]["title"] if r in ft else ta[r]["title"]) for r in union}
    pmid = {r: (ta[r]["pmid"] if r in ta else allrows[r]["pmid"]) for r in union}
    print("# Level B: the owner's sheet of undecided records")
    print()
    print(f"{len(ta)} records at titles and abstracts (`screening/to_owner.tsv`) and {len(ft)} at the full texts "
          f"(`fulltext_screening/ft_to_owner.txt`): {len(union)} records, {len(both)} in both lists. Generated by "
          "`decision_sheet.py` from those files and `runs/fulltext_index.tsv`; the summaries are written by hand "
          "from the passes' notes. TA: titles and abstracts (pass 1 · pass 2); FT: the full text (p1, p2 where read); "
          "codes as in `screening.md` and `screening_fulltext.md`.")
    for g, head, options in GROUPS:
        rows = groups[g]
        print()
        print(f"**{g}. {head}: {len(rows)}.** Options: {options}")
        print()
        print("| Record | Paper | Passes | Full text | What the notes say |")
        print("|---|---|---|---|---|")
        for r in rows:
            print(f"| {r} (PMID {pmid[r]}) | {short(title[r])} | {code_of(r, ta, ft, allrows)} | "
                  f"{fulltext_of(r, index)} | {SUMMARY[r]} |")
    print()
    print("A decision can name a group (\"A: exclude\") or a record (\"R02780: the other pass reads it\").")
    return 0


if __name__ == "__main__":
    sys.exit(main())
