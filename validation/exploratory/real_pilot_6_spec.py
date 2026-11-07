"""What the fourth run's plan fixes before the run (validation/exploratory/real_pilot_6_plan.md): the claims,
the criteria, the owner's prediction, and the sources and patterns of the check of what the TMS authors
write. The standard library only, as real_pilot_5_spec: the extraction, the text check and the analyses all
import it, and none reads a count to fill it. The inputs, the genes and the claims A1-A4 are the third
run's (real_pilot_5_spec), unchanged.
"""
from __future__ import annotations

from real_pilot_5_spec import (ACCEPT as ACCEPT5, ALPHA, CASE_SHA256, PAIR, RELEASE, RUN1_COMMIT,
                               RUN2_COMMIT, S4_NEG, S4_POS, SEX_GENES, SESOI_MEAN, SESOI_PAIR, TMS_DATASET, Y_GENES)

__all__ = ["ALPHA", "CASE_SHA256", "PAIR", "RELEASE", "RUN1_COMMIT", "RUN2_COMMIT", "S4_NEG", "S4_POS", "SEX_GENES",
           "SESOI_MEAN", "SESOI_PAIR", "TMS_DATASET", "Y_GENES", "RESULTS_BRANCH", "LABEL", "CRITERIA", "ACCEPT", "PARTS",
           "A3_DOMAIN_EXPECTED", "A8_THIRD_RUN", "PREDICTION", "PREDICTION_EN", "TMS_DOI", "TMS_PMCID",
           "TMS_PREPRINT_DOI", "TEXT_SOURCES", "TEXT_PATTERNS", "MAX_QUOTES", "MAX_QUOTE_CHARS"]

RESULTS_BRANCH = "results/real-pilot-6"
LABEL = "pilot-run-6"

# --------------------------------------------------------------------------- #
# the claims: the third run's A1-A4 (A3 and A4 their recorded second attempt) and its A8
# --------------------------------------------------------------------------- #
CRITERIA = {
    "A1": "SUPPORTED",
    "A2": "SUPPORTED",
    "A3": "not SUPPORTED, with the dossier",
    "A4": "not SUPPORTED, with the dossier",
    "A8": "no SUPPORTED; INCONCLUSIVE allowed",
}
ACCEPT = tuple(dict(c, criterion=CRITERIA[c["id"]], attempt=2 if c["id"] in ("A3", "A4") else 1)
               for c in ACCEPT5 if c["id"] in ("A1", "A2", "A3", "A4"))
PARTS = ("case2", "case1m", "case1", "case3a", "case3b", "case3c")

# the domain of A3 by the rule "in each group", read off the third run's logged table (a1) before this run
# (real_pilot_5_summary.log): the types whose Smad3 and Col1a1 are both detected in >= 10% of the males'
# cells of the type in each age group
A3_DOMAIN_EXPECTED = ("fibroblast", "fibroblast of cardiac tissue", "granulocytopoietic cell", "mesenchymal stem cell",
                      "mesenchymal stem cell of adipose tissue")

# the third run's verdicts on A8 (real_pilot_5.log), against which this run's are set, not scored
A8_THIRD_RUN = {"3a8": "INCONCLUSIVE", "3c3": "EFFECT BELOW SESOI"}  # every other of the 20: NO DETECTABLE EFFECT

# the owner's prediction for A3 and A4 (the second decisions of 2026-10-18, item 2), verbatim
PREDICTION = ("GATE 0 по-прежнему блокирует, потому что у молодых и старых самцов различаются число генов и счётов "
              "на клетку (nnz +1298, p 0,016; raw_sum p 0,008).")
PREDICTION_EN = ("GATE 0 still blocks, because the young and old males differ in genes and counts per cell (nnz "
                 "+1,298, p 0.016; raw_sum p 0.008).")

# --------------------------------------------------------------------------- #
# what the TMS authors write (the owner's item 8): fetched by the run, quoted, the texts not archived
# --------------------------------------------------------------------------- #
TMS_DOI, TMS_PMCID, TMS_PREPRINT_DOI = "10.1038/s41586-020-2496-1", "PMC8240505", "10.1101/661728"
TEXT_SOURCES = (
    ("tms-epmc", "the Nature paper, Europe PMC's full text (an NIH author manuscript)",
     f"https://www.ebi.ac.uk/europepmc/webservices/rest/{TMS_PMCID}/fullTextXML"),
    ("tms-bioc", "the Nature paper, NCBI's BioC (PMC open access and author manuscripts)",
     f"https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/{TMS_PMCID}/unicode"),
    ("tms-biorxiv", "the bioRxiv preprint, its JATS full text (the newest version, by bioRxiv's API)",
     f"https://api.biorxiv.org/details/biorxiv/{TMS_PREPRINT_DOI}"),
)
# a sentence is quoted where it matches one of these (case-insensitive regular expressions): the 3-month
# data's origin, the batches and the time points; short quotes only, the texts are not archived
TEXT_PATTERNS = (
    r"\b3[- ]?(?:m|mo|month)s?\b", r"\bthree[- ]month", r"previously (?:published|described|generated|reported)",
    r"\bTabula Muris\b(?! Senis)", r"\bbatch(?:es)?\b", r"\btime ?points?\b",
)
MAX_QUOTES, MAX_QUOTE_CHARS = 40, 600
