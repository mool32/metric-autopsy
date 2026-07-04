#!/usr/bin/env python3
"""Regenerate the arXiv LaTeX (main.tex) from manuscript.md, and the plain-text abstract.

    cd paper && python3 build_arxiv.py        # -> main.tex, _abstract_plain.txt
    tectonic main.tex                          # -> main.pdf  (or submit main.tex to arXiv)

Needs: pandoc. (tectonic optional, only to compile locally.) Run from the paper/ directory.
"""
import re
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
md = (HERE / "manuscript.md").read_text(encoding="utf-8")


def sanitize(s: str) -> str:
    """Map the few non-ASCII chars to pure-ASCII forms that render in text AND \\texttt
    (the × instances sit inside code spans, so dollar-math is unusable there)."""
    return (s.replace("×", "x").replace("·", " - ")
             .replace("§", "Section ").replace("—", "---").replace("–", "--"))


title = re.search(r"^# (.+)$", md, re.M).group(1).strip()
abstract = re.search(r"## Abstract\s*\n(.+?)\n\s*---", md, re.S).group(1).strip()
body = md[md.index("## 1. Introduction"):]
body = re.sub(r"^(#{2,3}) \d+(\.\d+)?\.?\s+", r"\1 ", body, flags=re.M)  # drop manual numbering

FIG1 = ("\n\n![Per-cell gene detection (n\\_genes) by sex and age in the synthetic worked "
        "example. Only male-old is QC-degraded (median 16 vs 31 genes), reproducing the "
        "Tabula Muris Senis interaction confound; pooling the sexes hides it."
        "](figures/fig1_qc_strata.png){width=78%}\n")
FIG2 = ("\n\n![Raw Smad3 vs Col1a1 counts, young vs old. The apparent between-group shape "
        "change is dropout collapsing points onto the zero axes (zero-fraction 0.21 vs 0.37), "
        "not a change in coupling.](figures/fig2_raw_scatter.png){width=88%}\n")


def insert_after(text, lead, fig):
    lines = text.split("\n")
    for i, ln in enumerate(lines):
        if ln.startswith(lead):
            lines.insert(i + 1, fig)
            return "\n".join(lines)
    raise SystemExit(f"anchor not found: {lead!r}")


body = insert_after(body, "**GATE 1 (QC parity).**", FIG1)
body = insert_after(body, "**GATE 3 (raw visibility).**", FIG2)

author = ("Theodor Spiro\\thanks{Vaika Inc.\\quad ORCID 0009-0004-5382-9346\\quad "
          "\\texttt{tspiro@vaika.org}\\quad Software: \\url{https://github.com/mool32/metric-autopsy}}")
meta = (f"---\ntitle: |\n  {title}\nauthor: '{author}'\ndate: \"\"\nabstract: |\n"
        + "".join("  " + ln + "\n" for ln in sanitize(abstract).split("\n")) + "---\n")

(HERE / "_body.md").write_text(sanitize(body), encoding="utf-8")
(HERE / "_meta.yaml").write_text(meta, encoding="utf-8")
plain = abstract.replace("`", "").replace("*", "").replace("–", "-").replace("—", "--")
(HERE / "_abstract_plain.txt").write_text(plain + "\n", encoding="utf-8")

subprocess.run([
    "pandoc", str(HERE / "_body.md"), "--metadata-file", str(HERE / "_meta.yaml"),
    "--standalone", "--shift-heading-level-by=-1", "--number-sections",
    "-V", "documentclass=article", "-V", "fontsize=11pt", "-V", "geometry:margin=1in",
    "-V", "colorlinks=true", "-V", "linkcolor=RoyalBlue", "-V", "urlcolor=RoyalBlue",
    "-o", str(HERE / "main.tex"),
], check=True)
(HERE / "_body.md").unlink(); (HERE / "_meta.yaml").unlink()  # keep only main.tex + abstract
print("wrote main.tex and _abstract_plain.txt")
