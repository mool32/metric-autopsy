# Submitting `metric-autopsy` to arXiv

Everything needed to post the preprint. The paper is a single-cell methods/tool paper; the
best-fit arXiv archive is **q-bio.QM** (Quantitative Methods).

## Files in this folder

| File | What | Upload to arXiv? |
|---|---|---|
| `arxiv-submission.tar.gz` | ready-to-upload source bundle: `main.tex` + `figures/` | **yes — upload this** |
| `main.tex` | the arXiv LaTeX source (pure ASCII, pdflatex-safe, numbered sections) | (inside the tarball) |
| `figures/fig1_qc_strata.png`, `fig2_raw_scatter.png` | the two figures `main.tex` references | (inside the tarball) |
| `main.pdf` | compiled preview — **15 pages, 4 figures** | **no** — arXiv builds the PDF from source |
| `build_arxiv.py` | regenerates `main.tex` from `manuscript.md` (needs `pandoc`) | no |

> Do **not** upload `main.pdf`. arXiv wants the LaTeX *source* and compiles it itself; uploading
> a PDF alongside the source confuses its AutoTeX. Upload only `arxiv-submission.tar.gz`.

## Before you start

- Log in at [arxiv.org](https://arxiv.org). You already have an arXiv paper, so an account exists.
- **Endorsement:** arXiv may require an endorsement for a *new* archive. If you have not posted
  to **q-bio** before, you might be prompted to request endorsement (a one-line request to an
  existing q-bio author, or arXiv auto-endorses established submitters). If prompted, follow the
  on-screen request flow; it is usually granted quickly. If you already posted to q-bio, skip this.
- **License:** we recommend **CC BY 4.0** to match the manuscript license in this repo.

## Step by step

1. **Start New Submission.**
2. **License** → choose **Creative Commons Attribution 4.0 (CC BY 4.0)**.
3. **Upload** `arxiv-submission.tar.gz`. arXiv detects LaTeX and compiles with pdfLaTeX. Wait for
   processing to succeed, then **view the generated PDF** and confirm it matches `main.pdf`
   (title page, abstract, 4 figures, ~15 pages). The source is standard LaTeX (a few UTF-8
   author names in the references compile fine under arXiv's pdfLaTeX via inputenc utf8).
4. **Add/verify metadata** — paste the blocks from the bottom of this file into the matching
   fields:
   - **Title** → the Title field.
   - **Authors** → `Spiro, Theodor` (arXiv "Surname, First" form). Add your **ORCID
     0009-0004-5382-9346** in your arXiv author profile so it links.
   - **Abstract** → the Abstract field (plain text below).
   - **Comments** → the Comments field (page count + software/Zenodo links).
   - **Primary category** → `q-bio.QM` (Quantitative Methods).
   - **Cross-list** (optional but recommended) → `q-bio.GN` (Genomics); optionally `stat.ME`
     (Methodology). Keep cross-lists minimal — 1–2 is plenty.
   - **Report-no / Journal-ref / DOI** → **leave blank.** (The Zenodo DOI goes in *Comments*, not
     the DOI field — that field is reserved for a journal-published version.)
   - **MSC / ACM class** → leave blank.
5. **Review → Submit.** New submissions go on a moderation hold and are announced at the next
   scheduled mailing (typically the next business day; the daily cutoff is ~14:00 US Eastern).

## After it is announced

Once arXiv assigns an identifier (e.g. `arXiv:26XX.XXXXX`), tell me and I will:
- point the README **Preprint** badge/link at the arXiv abstract page,
- add the arXiv ID to `CITATION.cff` (`preferred-citation`) and to `main.pdf`'s citation,
- update `PROJECT.md` and the portfolio master-index row,
- optionally add the entry to `mool32.github.io/_data/publications.yml` + `papers.bib`.

---

## Copy-paste blocks

### Title
```
Metric autopsy: a metric-agnostic gate system for separating biological signal from QC and technical artifacts in single-cell metrics
```

### Authors
```
Spiro, Theodor
```

### Abstract
```
Single-cell analyses routinely follow a "compute then believe" default: a scalar metric is calculated, it differs between conditions, and the difference is reported as biology. But a metric can move because of dropout, library size, a batch effect, a factorial-interaction confound, or plain mathematics -- with no change in the biological quantity of interest. This tendency is not hypothetical: three of my own analyses each survived weeks of work before a short quality-control (QC) check invalidated them -- in one case a 45-second QC check invalidated three weeks of downstream work. I present metric-autopsy, a metric-agnostic gate system that inverts the default from "compute then believe" to "state your commitments, red-team, then believe." The user supplies a metric as a black-box callable and the names of factorial obs columns; the harness probes the data and the metric's response to controlled perturbations of the data through eight gates (GATE 0-7). Six gates run automatically from the data -- mathematical independence, factorial QC parity, n_genes matching, raw-data visibility, stratified controls, and cross-dataset replication (the last only when a second dataset is supplied) -- and two are judgment gates the accompanying skill elicits rather than scripts. Because the automatic perturbations are single-cell-specific (dropout, depth, library size), the tool is metric-agnostic but domain-locked to scRNA-seq QC. I describe each gate faithful to the implementation, and demonstrate the workflow on a synthetic confound patterned on a mutual-information "coupling loss" driven by a sex-by-age detection artifact that was recorded, by hand, in an internal validation checklist; the real-data autopsy is not run here. The tool ships as a Claude Code skill and a pip package, and was itself subjected to an adversarial code audit (24 findings fixed). The gates are necessary, not sufficient: no correlation metric is fully depth-invariant under dropout.
```

### Comments
```
15 pages, 4 figures. Software (MIT): https://github.com/mool32/metric-autopsy ; archived at Zenodo https://doi.org/10.5281/zenodo.21195679
```

### Categories
```
Primary:   q-bio.QM   (Quantitative Methods)
Cross-list: q-bio.GN  (Genomics)   [optional: stat.ME]
```
