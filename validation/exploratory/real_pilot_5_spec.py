"""What the third diagnostics' plan fixes before the run (validation/exploratory/real_pilot_diagnostics_3_plan.md):
the registered tables of Marsh et al. 2022, the lists of genes, the rules' constants and the acceptance
claims. The standard library only: the extraction (v1's selection environment) and the analyses (the
engine's environment) both import it, and neither reads a count to fill it.
"""
from __future__ import annotations

RELEASE = "2025-11-08"
SEED = 20261018
ALPHA = 0.05
N_PERM = 9999                # the permutations of (b3)
MIN_UNIT_CELLS, MIN_UNITS = 20, 3   # a unit counts for a type with >= 20 of its cells; a type needs >= 3 units per group

# --------------------------------------------------------------------------- #
# the inputs, by commit and hash
# --------------------------------------------------------------------------- #
RUN1_COMMIT = "6133dfdf4c05f083f012aaec4136e376715af05a"      # results/real-pilot: case3.npz, select.json
RUN2_COMMIT = "bad785a3c56e9de1ee24737953c1b4b63112aa7c"      # results/real-pilot-2: case1.npz, case2.npz
MARSH_META_COMMIT = "757c2ab48af792f63b1992c8eb98b05a17d4272f"  # results/marsh-metadata
CASE_SHA256 = {"case1.npz": "13e7a563d4ed8fa482f976ead77b492189e8c729eb3ab103a203ccfdf8f24a44",
               "case2.npz": "7191ad731063a45110991a7fa7d3b8c15d43325bf7e9b6d95072165dc1655beb",
               "case3.npz": "31c297fe83adac9d93ff905695301bf69f09d583f9306ef980ae4210bedc8568"}
TMS_DATASET = "98e5ea9f-16d6-47ec-a529-686e76515e39"   # Tabula Muris Senis, FACS: case 1's and case 2's cells
MARSH_TAR = "GSE152183_RAW.tar"
MARSH_TAR_BYTES = 76861440    # GEO's file list (results/marsh-metadata)
GEO_TMS, GEO_TM = "GSE132042", "GSE109774"

# --------------------------------------------------------------------------- #
# (b): §4's technical counts
# --------------------------------------------------------------------------- #
PAIR = ("Smad3", "Col1a1")
FIB_WORDS = ("fibroblast", "mesenchymal", "stromal")   # the second diagnostics' label rule
# Muhl et al. 2020 (Nat Commun 11:3953): the matrix genes of their short-list of common fibroblast
# markers, named in the text; and the owner's examples not among them (Col3a1, Dcn, Sparc)
MUHL_ECM = ("Col1a1", "Col1a2", "Col5a1", "Loxl1", "Lum", "Fbln1", "Fbln2")
OWNER_EXAMPLES = ("Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Sparc")
FIB_GENES = ("Col1a1", "Col1a2", "Col3a1", "Col5a1", "Dcn", "Fbln1", "Fbln2", "Loxl1", "Lum", "Sparc")
# the index of foreign counts: the list less the claim's gene and less Sparc (in Marsh et al.'s
# microglial identity score, Supplementary Table 4: expressed outside the fibroblasts by name)
INDEX_GENES = tuple(g for g in FIB_GENES if g not in ("Col1a1", "Sparc"))
INDEX_PLUS = tuple(g for g in FIB_GENES if g != "Sparc")
YOUNG_UNITS = ("3_8_M", "3_9_M", "3_10_M", "3_11_M", "3_38_F", "3_39_F", "3_56_F", "3_57_F")  # non-composite
OLD_UNITS = ("24_58_M", "24_59_M", "24_60_M", "24_61_M", "21_48_F", "21_54_F", "21_55_F")

# --------------------------------------------------------------------------- #
# Marsh et al. 2022: the registration (GSE152183's SOFT records and Supplementary Table 1)
# --------------------------------------------------------------------------- #
ARMS = ("ENZ-NONE", "ENZ-INHIB", "DNC-NONE", "DNC-INHIB")
# gsm: (GEO title, arm, batch)
MARSH_SAMPLES = {
    "GSM4605687": ("Dounce_NO_sample_01", "DNC-NONE", "Batch 1"),
    "GSM4605688": ("Miltenyi_NO_sample_04", "ENZ-NONE", "Batch 1"),
    "GSM4605689": ("Dounce_IN_sample_04", "DNC-INHIB", "Batch 1"),
    "GSM4605690": ("Miltenyi_IN_sample_01", "ENZ-INHIB", "Batch 1"),
    "GSM4605691": ("Dounce_NO_sample_02", "DNC-NONE", "Batch 1"),
    "GSM4605692": ("Miltenyi_NO_sample_05", "ENZ-NONE", "Batch 1"),
    "GSM4605693": ("Dounce_IN_sample_05", "DNC-INHIB", "Batch 2"),
    "GSM4605694": ("Miltenyi_IN_sample_02", "ENZ-INHIB", "Batch 2"),
    "GSM4605695": ("Dounce_NO_sample_03", "DNC-NONE", "Batch 2"),
    "GSM4605696": ("Miltenyi_NO_sample_06", "ENZ-NONE", "Batch 2"),
    "GSM4605697": ("Dounce_IN_sample_06", "DNC-INHIB", "Batch 2"),
    "GSM4605698": ("Miltenyi_IN_sample_03", "ENZ-INHIB", "Batch 2"),
}
# the SOFT characteristics each arm must carry
ARM_CHARACTERISTICS = {"ENZ-NONE": ("Enzymatic", "No"), "ENZ-INHIB": ("Enzymatic", "Yes"),
                       "DNC-NONE": ("Cold Dounce", "No"), "DNC-INHIB": ("Cold Dounce", "Yes")}
MARSH_POST_QC = {  # Supplementary Table 1: cells after the authors' QC, by arm and batch, for the report
    ("DNC-NONE", "Batch 1"): (1916, 991), ("ENZ-NONE", "Batch 1"): (1631, 1176), ("DNC-INHIB", "Batch 1"): (1632,),
    ("ENZ-INHIB", "Batch 1"): (1434,), ("DNC-INHIB", "Batch 2"): (2040, 2209), ("ENZ-INHIB", "Batch 2"): (1675, 1591),
    ("DNC-NONE", "Batch 2"): (1528,), ("ENZ-NONE", "Batch 2"): (1740,)}
MARSH_SET = ("Rgs1", "Hist2h2aa1", "Hist1h4i", "Nfkbiz", "Klf2", "Junb", "Dusp1", "Ccl3", "Hspa1a", "Hsp90aa1",
             "Fos", "Hspa1b", "Jun", "Jund", "Nfkbid", "Gem", "Ccl4", "Ier5", "Txnip", "Hist1h2bc", "Zfp36",
             "Hist1h1c", "Egr1", "Atf3", "Rhob")   # Supplementary Table 4, "Micro/Myeloid Shared Act. Score"
MARSH_PAIR = ("Fos", "Jun")
MARSH_TYPE = "CD45+/CD11b+ sorted cells"
CONTRASTS = {"M1": (("ENZ-INHIB", "ENZ-NONE"), "increase", "primary: the same exposure, the induction blocked"),
             "M2": (("DNC-NONE", "ENZ-NONE"), "increase", "secondary: against cold Dounce"),
             "M3": (("DNC-NONE", "DNC-INHIB"), "two-sided", "null: Dounce with the inhibitors against without")}
MARSH_CAVEAT = ("the genes' mean and detection move together with their coupling: this control checks that the "
                "engine does not refuse a real change of coupling, not that it tells a change of coupling from a "
                "change of detection")

# --------------------------------------------------------------------------- #
# the acceptance tests
# --------------------------------------------------------------------------- #
SESOI_PAIR, SESOI_MEAN = 0.1, 0.5      # on the construct's scale
Y_GENES = ("Ddx3y", "Eif2s3y", "Kdm5d", "Uty")
SEX_GENES = ("Xist", "Tsix", "Jpx", "Ftx",                                    # the first diagnostics' list
             "Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Zfy1", "Zfy2", "Sry", "Usp9y", "Ube1y1",
             "Ddx3x", "Eif2s3x", "Kdm5c", "Kdm6a", "Zfx", "Usp9x", "Uba1",
             "Mid1", "Erdr1", "Sts")
S4_POS, S4_NEG = ("Actb", "Gapdh"), ("Tst", "Lrrc42")
ACCEPT = (
    dict(id="A1", part="case2", case="case2", metric="log_cp10k_mean", genes=["Xist"], groups=("female", "male"),
         group_col="sex", replicate="mouse", within=[], direction="decrease", sesoi=SESOI_MEAN, probe="module",
         composition=["cell_type"], known_to_differ=True, criterion="SUPPORTED"),
    dict(id="A2", part="case2", case="case2", metric="log_cp10k_mean", genes=list(Y_GENES), groups=("female", "male"),
         group_col="sex", replicate="mouse", within=[], direction="increase", sesoi=SESOI_MEAN, probe="module",
         composition=["cell_type"], known_to_differ=True, criterion="SUPPORTED"),
    dict(id="A3", part="case1m", case="case1", males=True, metric="mi_3bin", genes=list(PAIR), groups=("young", "old"),
         group_col="age", replicate="mouse", within=[], direction="decrease", sesoi=SESOI_PAIR, probe="coupling",
         controls=True, composition=["cell_type", "tissue"], criterion="not SUPPORTED, no effect in the domain, the drop outside"),
    dict(id="A4", part="case1", case="case1", metric="mi_3bin", genes=list(PAIR), groups=("young", "old"),
         group_col="age", replicate="mouse", within=["sex"], direction="decrease", sesoi=SESOI_PAIR, probe="coupling",
         controls=True, composition=["cell_type", "tissue"], criterion="not SUPPORTED"),
    *[dict(id=f"A{4 + i}{m}", part="marsh", case="marsh", contrast=c, metric=metric, genes=list(MARSH_PAIR),
           groups=CONTRASTS[c][0], group_col="arm", replicate="mouse", within=[], direction=CONTRASTS[c][1],
           sesoi=SESOI_PAIR, probe="coupling", composition=["cell_type"],
           criterion=({"M1": "SUPPORTED", "M3": "not SUPPORTED"}.get(c, "reported") if not m else "reported"))
      for m, metric in (("", "norm_pearson"), ("m", "mi_3bin")) for i, c in enumerate(("M1", "M2", "M3"), 1)],
)
CASE3_PARTS = ("3a", "3b", "3c")   # A8: the first run's 20 claims, in three parts by prefix
PARTS = ("b", "a", "case2", "case1m", "case1", "marsh", *[f"case{p}" for p in CASE3_PARTS])
