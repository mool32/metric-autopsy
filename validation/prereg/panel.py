"""The blinded panel of the confirmatory validation (validation/prereg/v1.md, section 3).

This module holds the design (conditions, variants, numbers of datasets), the key assignment,
the backgrounds' fixed genes and pairs, the truth generators, the claim cards and the allowed
(label, cause) outcomes. The panel is built on the fly by ``blind.py`` (in GitHub Actions or one
local run): no dataset is stored, every dataset and claim card is identified by a canonical
sha256 (``dataset_sha256``, ``card_sha256``), and the datasets are reproducible from the key.

The key is the randomness of a public drand round named in the run tag before it exists
(``beacon.py``): 256 bits written as 64 lowercase hex characters. It decides every dataset's
condition, variant (dose or step), side, gene pair and seed, and the order of the dataset IDs;
the pair is drawn independently of the condition.

Independence: the truth is generated here, with this module's own binomial thinning (as in
seqgendiff, Gerard 2020), never with the engine's modules; this file does not import
``metric_autopsy``.
"""
from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

KEY_N = 790            # datasets per key null condition (oc.log: S1 as a whole, 5 conditions)
DONORS_PER_GROUP = 8   # N1-N3, N6, N8, E1-E3: 2 x 8 donors from B1
N7_MAX_MICE = 24       # N7: at most 24 of B2's qualifying mice per dataset, so its memory is bounded
CELLS_PER_DONOR = 200  # cells drawn per donor and dataset; donors with fewer are not used
E_CAPTURE = 0.5        # E2, E3: the capture of the side that loses it (the correction of depth brings the other to it)
N_GENES = 2000         # genes kept per background: the named genes, then the most expressed
RANDOM_SCORE_GENES = 50
LEVELS = ("high", "medium", "low")
MAX_PAIRS_PER_LEVEL = 8  # the pool takes 8 pairs per level where every level of the background has
MIN_PAIRS_PER_LEVEL = 4  # them, else the largest common number >= 4 (decided 2026-10-08)
LEVEL_CANDIDATES = 1000  # the most expressed eligible genes of a level enter its pair search
PLAN_CELLS_PER_DONOR = 200  # cells per donor (public seed) on which pair correlations are computed
MIN_DONOR_DETECTION = 0.05  # a pair gene is detected in >= 5% of the cells of every donor
NEG_MAX_R = 0.02       # a negative control's |partial correlation| is below this
NEG_CANDIDATES = 10    # of the 10 such pairs closest in mean expression, the most typical of GATE 5's null
# GATE 5's matched null for a negative control (copies of the engine's rule, checked against it in
# test_the_panels_copies_of_the_engines_rules_agree_with_it): pairs from the two genes' expression
# neighbourhoods, each the closest genes in mean count, at least 20 or 5% of the genes; 200 pairs
MATCHED_NEIGHBOURS_MIN = 20
MATCHED_NEIGHBOURS_FRAC = 0.05
NEG_NULL_PAIRS = 200
N8_BETA = (2.0, 2.0)   # N8: per-cell capture ~ Beta(2, 2) in one group (mean 0.5)
PLAN_SEED = 20261008   # public seed of the plan's cell subsample and N6c's random genes
DELTA_MIN_FRACTION = 0.5  # delta_min = 0.5 x SESOI: the smallest response to the injection that matters
TRUTH_BAND = 0.2       # the metric is valid at >= 1.2 delta_min, blind at <= 0.8 delta_min, else ambiguous

# Tirosh et al. 2016 G2/M genes (human symbols); the module the random-gene score claims to
# measure (N6c). Genes absent from a background are dropped.
G2M_GENES = (
    "HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 "
    "PIMREG SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 JPT1 "
    "CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN "
    "LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA").split()

# --------------------------------------------------------------------------- #
# the outcomes a verdict is scored by: its label and its cause (v1.md, section 3.2)
# --------------------------------------------------------------------------- #
SUPPORTED = "SUPPORTED"
NDE = "NO DETECTABLE EFFECT"
INCONCLUSIVE = "INCONCLUSIVE"
NS_INVALID = "NOT SUPPORTED: metric invalid (GATE 4/5)"
NS_DEPTH = "NOT SUPPORTED: explained by depth"
NS_OPPOSITE = "NOT SUPPORTED: opposite direction"
DEGENERATE = "DEGENERATE METRIC"
REFUSAL = "REFUSAL: nuisance bias (GATE 0)"
UNIDENTIFIABLE = "UNIDENTIFIABLE"
OTHER = "OTHER"  # a label and cause the panel does not expect (e.g. not replicated: no GATE 6 here)
ERROR = "ERROR"  # no report, an engine error, or a label that does not match its cause
OUTCOMES = (SUPPORTED, NDE, INCONCLUSIVE, NS_INVALID, NS_DEPTH, NS_OPPOSITE, DEGENERATE, REFUSAL,
            UNIDENTIFIABLE, OTHER, ERROR)
DEFINITE = frozenset({SUPPORTED, NDE, NS_INVALID, NS_DEPTH, NS_OPPOSITE, DEGENERATE, UNIDENTIFIABLE})
LABELS = ("SUPPORTED", "NOT SUPPORTED", "NO DETECTABLE EFFECT", "INCONCLUSIVE", "UNIDENTIFIABLE",
          "DEGENERATE METRIC")

# The engine's cause codes (metric_autopsy.report.CAUSES; test_the_panels_copies_of_the_engines_rules_agree_with_it)
# and the outcome each one is scored as. GATE 0's block for a nuisance bias is a refusal: allowed
# everywhere, never definite (decided 2026-10-08).
_INCONCLUSIVE_CAUSES = ("effect_not_evaluated", "insufficient_replication", "absence_untested_metric",
                        "effect_inconclusive", "detected_untested_metric", "bias_unsized", "no_direction",
                        "judgment_pending")
CAUSE_OUTCOME = {"replicated": SUPPORTED, "provisional": SUPPORTED, "no_detectable_effect": NDE,
                 "metric_invalid_gate4": NS_INVALID, "metric_invalid_gate5": NS_INVALID,
                 "metric_invalid_gate0": REFUSAL, "explained_by_depth": NS_DEPTH,
                 "opposite_direction": NS_OPPOSITE, "degenerate_metric": DEGENERATE,
                 "unidentifiable": UNIDENTIFIABLE, "not_replicated": OTHER, "metric_invalid": OTHER,
                 **{c: INCONCLUSIVE for c in _INCONCLUSIVE_CAUSES}}
_OUTCOME_LABEL = {SUPPORTED: "SUPPORTED", NDE: "NO DETECTABLE EFFECT", INCONCLUSIVE: "INCONCLUSIVE",
                  NS_INVALID: "NOT SUPPORTED", NS_DEPTH: "NOT SUPPORTED", NS_OPPOSITE: "NOT SUPPORTED",
                  REFUSAL: "NOT SUPPORTED", OTHER: "NOT SUPPORTED", DEGENERATE: "DEGENERATE METRIC",
                  UNIDENTIFIABLE: "UNIDENTIFIABLE"}


def label(verdict: str | None) -> str | None:
    """The verdict's label: the text before ' — ', without bracketed qualifiers.
    'SUPPORTED (provisional until replicated)' and 'SUPPORTED — replicated' are SUPPORTED."""
    if not verdict:
        return None
    head = verdict.split(" — ")[0].split(" [")[0].split(" (")[0].strip()
    return head if head in LABELS else None


def outcome(verdict: str | None, cause: str | None) -> str:
    """The outcome a report is scored as: its cause's outcome, if its label is that outcome's."""
    out = CAUSE_OUTCOME.get(cause or "")
    if out is None or label(verdict) != _OUTCOME_LABEL[out]:
        return ERROR
    return out


# --------------------------------------------------------------------------- #
# the design: every condition with its variants and number of datasets
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Condition:
    name: str
    background: str
    variants: tuple            # ((variant label, number of datasets), ...)
    data: str                  # the truth about the data: "null" (no biology) or "effect" (a real one)
    artifact: bool = False     # a technical artifact is planted (N2, N3, N8; E2, E3 with the effect)
    metric: str = "norm_pearson"  # or a useless metric: "random", "constant", "score" (N6)
    key: bool = False          # key null condition (criterion S1), at its key variant
    key_variant: str | None = None
    cards: int = 1             # claim cards per dataset (N4: without and with a replicate column)
    oracle: bool = True        # whether the oracle can establish a definite verdict (N4: no)


CONDITIONS = (
    Condition("N1", "B1", (("null", KEY_N),), "null", key=True, key_variant="null"),
    Condition("N2", "B1", (("c=0.5", KEY_N), ("c=0.9", 100), ("c=0.7", 100), ("c=0.3", 100)), "null",
              artifact=True, key=True, key_variant="c=0.5"),
    Condition("N3", "B1", (("f=0.1", 100), ("f=0.2", 100), ("f=0.4", 100)), "null", artifact=True),
    Condition("N4", "B1", (("3v3", 300),), "null", cards=2, oracle=False),
    Condition("N5", "B1", (("sham", KEY_N),), "null", key=True, key_variant="sham"),
    Condition("N6a", "B1", (("random", 300),), "null", metric="random"),
    Condition("N6b", "B1", (("constant", 50),), "null", metric="constant"),
    Condition("N6c", "B1", (("random-genes", KEY_N),), "null", metric="score", key=True,
              key_variant="random-genes"),
    Condition("N7", "B2", (("mice", 300),), "null"),
    Condition("N8", "B1", (("beta(2,2)", KEY_N),), "null", artifact=True, key=True, key_variant="beta(2,2)"),
    Condition("E1", "B1", (("dose=key", 200), ("dose=0.25", 100), ("dose=0.5", 100), ("dose=1.5", 100)),
              "effect"),
    Condition("E2", "B1", (("against", 200),), "effect", artifact=True),
    Condition("E3", "B1", (("with", 200),), "effect", artifact=True),
)
BACKGROUNDS = ("B1", "B2")
# Conditions whose data change the pair's counts (capture loss, dropout, variable capture, an
# injected coupling) or the dataset's size (N4: 2 x 3 donors; N5: 8 donors whose cells are split at
# random, so GATE 4's precision differs: the third review): the truth about the metric and GATE 4's odds are
# measured on their own datasets (pilot.json "truth_case"); the others share their background's
# null (N1 on B1, N7 on B2).
TRUTH_CASE_CONDITIONS = ("N2", "N3", "N4", "N5", "N8", "E1", "E2", "E3")

# The drop order of v1.md section 4 (compute budget): whole variants, never replicates of the
# rest, and never a key null condition or E1-E3 at the key dose.
DROP_ORDER = (("N3", "f=0.1"), ("N3", "f=0.2"), ("N3", "f=0.4"), ("E1", "dose=0.25"), ("E1", "dose=0.5"),
              ("E1", "dose=1.5"), ("N2", "c=0.9"), ("N2", "c=0.7"), ("N2", "c=0.3"), ("N4", "3v3"))


def conditions() -> dict:
    return {c.name: c for c in CONDITIONS}


def _dropped(name: str, variant: str, dropped) -> bool:
    background = conditions()[name].background if name in conditions() else None
    return any(d in (name, f"{name}:{variant}", f"{background}:*") for d in (dropped or ()))


def background_drops(present) -> list:
    """The drop of every background without a qualifying candidate (v1.md 3.1: dropped with the
    cases that need it): 'B2:*' drops N7 (and the anchors drop R1). B1 carries every other case, so
    without it there is no validation to run."""
    if "B1" not in present:
        raise SystemExit("B1 has no qualifying candidate: the cases N1-N6, N8 and E1-E3 cannot be built "
                         "and the validation cannot run (v1.md 3.1)")
    return [f"{b}:*" for b in BACKGROUNDS if b not in present]


def check_dropped(dropped) -> list:
    """Drops are allowed only as v1.md section 4 and 3.1 fix them: the backgrounds without a
    qualifying candidate (`background_drops`; only B2), then a prefix of DROP_ORDER."""
    dropped = list(dropped or ())
    backgrounds = [d for d in dropped if d.endswith(":*")]
    if backgrounds and (backgrounds != dropped[:len(backgrounds)] or set(backgrounds) - {"B2:*"}):
        raise ValueError(f"only B2 can be dropped as a background, before the order's drops; got {dropped}")
    rest = dropped[len(backgrounds):]
    allowed = [f"{n}:{v}" for n, v in DROP_ORDER]
    if rest != allowed[:len(rest)]:
        raise ValueError(f"drops must be a prefix of the pre-registered order {allowed}, got {rest}")
    return dropped


def n_datasets(dropped=()) -> int:
    return sum(n for c in CONDITIONS for v, n in c.variants if not _dropped(c.name, v, dropped))


def n_cards(dropped=(), data: str | None = None) -> int:
    return sum(n * c.cards for c in CONDITIONS for v, n in c.variants
               if not _dropped(c.name, v, dropped) and (data is None or c.data == data))


def pairs_per_level(pool: list) -> int:
    """The pool's pairs per level (the same for every level of a background)."""
    return len(pool) // len(LEVELS)


def level_of_pair(pair_index: int, pool_size: int) -> str:
    """Pool entries are ordered by level: the first pool_size / 3 are high, then medium, then low."""
    return LEVELS[int(pair_index) // (int(pool_size) // len(LEVELS))]


# --------------------------------------------------------------------------- #
# the key: condition, variant, side, gene pair and seed of every dataset ID
# --------------------------------------------------------------------------- #
def check_key(key: str) -> str:
    """The key: a drand round's randomness, 256 bits as 64 lowercase hex characters."""
    if not isinstance(key, str) or not re.fullmatch(r"[0-9a-f]{64}", key):
        raise ValueError("the key is a drand randomness: 64 lowercase hex characters (beacon.py)")
    return key


DEFAULT_POOL_SIZES = {b: len(LEVELS) * MAX_PAIRS_PER_LEVEL for b in BACKGROUNDS}


def assign(key: str, dropped=(), pool_sizes: dict | None = None) -> list[dict]:
    """Every dataset of the design, in an order and with sides, gene pairs and seeds decided by
    the key. Deterministic in the key; nothing about an entry is visible in its ID. The gene pair
    is drawn independently of the condition: a 62-bit integer per dataset, reduced modulo the
    size of its background's pool (fixed before the key, `pool_sizes`; pilot.json records it),
    so that every pair of the pool is equally likely."""
    dropped = check_dropped(dropped)
    sizes = dict(DEFAULT_POOL_SIZES, **(pool_sizes or {}))
    entries = [dict(condition=c.name, variant=v, index=i)
               for c in CONDITIONS for v, n in c.variants if not _dropped(c.name, v, dropped)
               for i in range(n)]
    ss = np.random.SeedSequence(int(check_key(key), 16))
    order_seq, *data_seqs = ss.spawn(len(entries) + 1)
    perm = np.random.default_rng(order_seq).permutation(len(entries))
    width = len(str(len(entries)))
    out = []
    conds = conditions()
    for rank, k in enumerate(perm):
        e = dict(entries[k])
        draw = np.random.default_rng(data_seqs[k].spawn(1)[0])
        side = ("A", "B")[int(draw.integers(2))]
        pair_draw = int(draw.integers(2 ** 62))
        e.update(id=f"D{rank + 1:0{width}d}", side=side, pair_draw=pair_draw,
                 pair=pair_draw % int(sizes[conds[e["condition"]].background]),
                 seed=int(data_seqs[k].generate_state(1, np.uint64)[0]))
        out.append(e)
    return out


def card_ids(entry: dict) -> list[str]:
    return [entry["id"]] if conditions()[entry["condition"]].cards == 1 else [f"{entry['id']}a", f"{entry['id']}b"]


# --------------------------------------------------------------------------- #
# the truth and the allowed outcomes (pilot.json fixes the truth before the key)
# --------------------------------------------------------------------------- #
def e1_factor(variant: str) -> float:
    return 1.0 if variant == "dose=key" else float(variant.split("=")[1])


def pool_level(cond: Condition, pair_index: int, pilot: dict) -> str:
    return pilot["pool"][cond.background][int(pair_index)]["level"]


def sesoi_of(cond: Condition, pair_index: int, pilot: dict) -> float:
    """The SESOI of the pair's level on the condition's background (one rule for every
    background, computed on that background: v1.md 3.2)."""
    return float(pilot["sesoi"][cond.background][pool_level(cond, pair_index, pilot)])


def delta_min_of(cond: Condition, pair_index: int, pilot: dict) -> float:
    return DELTA_MIN_FRACTION * sesoi_of(cond, pair_index, pilot)


def delta_of(cond_name: str, variant: str, pair_index: int, pilot: dict) -> dict:
    """Δ*, the population difference of the metric (signal side minus the other) in the version
    of the condition without its artifact, for the dataset's pair: E1 at its dose; E2 and E3 at
    the key dose and at the depth their analysis runs at, both sides at the capture E_CAPTURE (the
    correction of depth thins the side that kept its depth: the fourth review found Δ* at full
    depth overstating the effect their data carry, 0.12-0.69 of it on simulated data) (pilot.json,
    written before the key)."""
    if cond_name in ("E2", "E3"):
        return pilot["delta"][str(pair_index)][f"capture={E_CAPTURE:g}"]
    factor = e1_factor(variant) if cond_name == "E1" else 1.0
    return pilot["delta"][str(pair_index)][f"{factor:g}"]


def classify_response(response: float, delta_min: float, band: float = TRUTH_BAND) -> str:
    """The truth about the metric on a pair, from its population response to the GATE 4
    injection: valid at >= (1 + band) delta_min, blind at <= (1 - band) delta_min, else ambiguous."""
    if response >= (1 + band) * delta_min:
        return "valid"
    if response <= (1 - band) * delta_min:
        return "blind"
    return "ambiguous"


def truth_record(cond: Condition, variant: str, pair_index: int, pilot: dict) -> dict:
    """The oracle's record of the metric's population response on this case: the condition's own
    datasets where its data change the pair's counts (pilot.json "truth_case"), else the
    background's null (pilot.json "truth"). A pilot that measured its background's cases but not
    this one is refused (the fourth review: the background's null stood in for it silently); a
    stand-in without case records (oc.py's scenarios) takes the background's null."""
    cases = (pilot.get("truth_case") or {}).get(cond.background)
    if cond.name in TRUTH_CASE_CONDITIONS and cases:
        rec = (cases.get(f"{cond.name}:{variant}") or {}).get(str(int(pair_index)))
        if rec is None:
            raise ValueError(f"pilot.json has no truth for {cond.name}:{variant} on pair {pair_index} "
                             f"of {cond.background}, whose other cases it measured")
        return rec
    return pilot["truth"][cond.background][str(int(pair_index))]


def metric_truth(cond: Condition, variant: str, pair_index: int, pilot: dict) -> str:
    """'valid', 'blind' or 'ambiguous' for log-normalized Pearson on the pair, on the data of this
    condition and variant (the oracle's population response, pilot.json); a useless metric (N6)
    is 'useless', the constant 'constant'."""
    if cond.metric == "constant":
        return "constant"
    if cond.metric != "norm_pearson":
        return "useless"
    return truth_record(cond, variant, pair_index, pilot)["class"]


INVALID_ALLOWED = frozenset({NS_INVALID, INCONCLUSIVE, REFUSAL})


def data_allowed(cond: Condition, variant: str, pair_index: int, pilot: dict) -> frozenset:
    """The outcomes correct for a valid metric, given the data (v1.md 3.2): on null data NO
    DETECTABLE EFFECT, INCONCLUSIVE and NOT SUPPORTED for the opposite direction, plus NOT
    SUPPORTED explained by depth where an artifact is planted; on a real effect SUPPORTED and
    INCONCLUSIVE, plus NO DETECTABLE EFFECT where |Δ*| < SESOI. A refusal by GATE 0 always."""
    if cond.data == "null":
        ok = {NDE, INCONCLUSIVE, NS_OPPOSITE, REFUSAL} | ({NS_DEPTH} if cond.artifact else set())
    else:
        d = float(delta_of(cond.name, variant, pair_index, pilot)["value"])
        small = abs(d) < sesoi_of(cond, pair_index, pilot)
        # below the SESOI, "explained by depth" where an artifact is planted says what NDE says: with
        # a SESOI the engine calls depth only when the corrected effect is shown smaller (journal D7)
        ok = ({SUPPORTED, INCONCLUSIVE, REFUSAL} | ({NDE} if small else set())
              | ({NS_DEPTH} if small and cond.artifact else set()))
    return frozenset(ok)


def allowed(cond: Condition, variant: str, pair_index: int, pilot: dict) -> frozenset:
    """The outcomes correct for this card (the same rule for every condition, decided
    2026-10-08): a blind or useless metric allows NOT SUPPORTED (metric invalid), INCONCLUSIVE
    and a refusal, and the constant also DEGENERATE METRIC; a valid metric allows
    `data_allowed`; an ambiguous one the union of both sets."""
    truth = metric_truth(cond, variant, pair_index, pilot)
    if truth == "useless":
        return INVALID_ALLOWED
    if truth == "constant":
        return INVALID_ALLOWED | {DEGENERATE}
    if truth == "blind":
        return INVALID_ALLOWED
    ok = data_allowed(cond, variant, pair_index, pilot)
    return ok if truth == "valid" else ok | INVALID_ALLOWED


def definite(cond: Condition, variant: str, pair_index: int, pilot: dict) -> frozenset:
    """The correct definite outcomes: allowed, and neither INCONCLUSIVE nor a refusal."""
    return allowed(cond, variant, pair_index, pilot) & DEFINITE


# The effect verdicts, which the engine's rules exclude on some cards whatever the data (the fourth
# review): without a replicate unit (N4's first card: effect.py gives no effect verdict without one)
# and on the constant metric (DEGENERATE METRIC is decided first: report.decide_cause). There a
# sound validator never gives them: one is an unexpected verdict, counted in S7 (none allowed).
EFFECT_VERDICTS = frozenset({SUPPORTED, NDE, NS_OPPOSITE, NS_DEPTH})


def excluded(cond: Condition, card: int) -> frozenset:
    """The outcomes the engine's rules exclude on a dataset's `card`-th claim card (``card_ids``;
    N4's first card is the one without the replicate unit, ``build``)."""
    if (cond.cards == 2 and card == 0) or cond.metric == "constant":
        return EFFECT_VERDICTS
    return frozenset()


def card_allowed(cond: Condition, variant: str, pair_index: int, pilot: dict, card: int) -> frozenset:
    """The outcomes correct for a dataset's `card`-th claim card: ``allowed`` without those the
    engine's rules exclude on it (``excluded``)."""
    return allowed(cond, variant, pair_index, pilot) - excluded(cond, card)


# --------------------------------------------------------------------------- #
# backgrounds: genes, the pair pool and its controls, fixed before the key
# --------------------------------------------------------------------------- #
@dataclass
class Background:
    """Raw counts of one background (cells x genes; numpy or scipy.sparse) with the donor of
    every cell. `plan_background` fixes the plan and compacts it to the qualifying donors' cells
    and the kept genes (dense int32)."""
    X: object
    genes: list
    donor: np.ndarray
    name: str = "B?"
    extra_obs: pd.DataFrame | None = None
    plan: dict = field(default_factory=dict)

    @property
    def donors(self) -> list:
        vals, counts = np.unique(self.donor, return_counts=True)
        return [str(v) for v, c in zip(vals, counts) if c >= CELLS_PER_DONOR]


def _is_sparse(X) -> bool:
    return hasattr(X, "tocsr") and not isinstance(X, np.ndarray)


def _counts_matrix(ad, spec: dict):
    """The raw counts of an AnnData and gene symbols: spec['counts'] is 'raw' (raw.X), 'X', or
    'layer:<name>'; 'auto' (default) takes raw.X if present, else layers['counts'], else X.
    Symbols come from spec['gene_symbols'] (a var column; 'feature_name' if present), else var_names."""
    how = spec.get("counts", "auto")
    if how == "auto":
        how = "raw" if ad.raw is not None else ("layer:counts" if "counts" in ad.layers else "X")
    if how == "raw":
        X, var = ad.raw.X, ad.raw.var
    elif how.startswith("layer:"):
        X, var = ad.layers[how.split(":", 1)[1]], ad.var
    elif how == "X":
        X, var = ad.X, ad.var
    else:
        raise ValueError(f"unknown counts source {how!r}")
    col = spec.get("gene_symbols") or ("feature_name" if "feature_name" in var.columns
                                       or "feature_name" in ad.var.columns else None)
    if col is None:
        genes = [str(g) for g in var.index]
    elif col in var.columns:
        genes = [str(g) for g in var[col]]
    else:
        genes = [str(g) for g in ad.var.loc[var.index, col]]
    return X, genes


def unique_names(names) -> list[str]:
    """Gene symbols made unique as anndata's var_names_make_unique does: the second occurrence of a
    symbol becomes 'symbol-1', the third 'symbol-2', skipping a name already taken; first
    occurrences keep their symbol. Several Ensembl genes can share one symbol, and the engine
    refuses duplicate var_names (the second review's K2)."""
    names = [str(n) for n in names]
    taken, seen, out = set(names), {}, []
    for n in names:
        if n not in seen:
            seen[n] = 0
            out.append(n)
            continue
        k = seen[n]
        while f"{n}-{k + 1}" in taken:
            k += 1
        k += 1
        seen[n] = k
        taken.add(f"{n}-{k}")
        out.append(f"{n}-{k}")
    return out


def load_background(name: str, spec: dict) -> Background:
    """An .h5ad (raw counts as `_counts_matrix` finds them; cells filtered by spec['filter'] before
    anything is loaded into memory) or an .npz written by `save_npz`. Counts stay sparse; gene
    symbols are made unique (`unique_names`)."""
    path = Path(spec["path"])
    filters = spec.get("filter") or {}
    if path.suffix == ".npz":
        z = np.load(path, allow_pickle=False)
        obs = pd.DataFrame({k[4:]: z[k] for k in z.files if k.startswith("obs_")})
        if "X_data" in z.files:
            from scipy import sparse
            X = sparse.csr_matrix((z["X_data"], z["X_indices"], z["X_indptr"]), shape=tuple(z["X_shape"]))
        else:
            X = z["X"]
        genes = unique_names(z["genes"])
        keep = np.ones(len(obs), bool)
        for col, val in filters.items():
            keep &= np.asarray(obs[col]).astype(str) == str(val)
        X, obs = X[keep], obs[keep].reset_index(drop=True)
    else:
        import anndata  # lazily: only the panel run needs it
        ad = anndata.read_h5ad(path, backed="r")
        keep = np.ones(ad.n_obs, bool)
        for col, val in filters.items():
            keep &= np.asarray(ad.obs[col]).astype(str) == str(val)
        sub = ad[np.where(keep)[0]].to_memory()
        X, genes = _counts_matrix(sub, spec)
        genes = unique_names(genes)
        obs = sub.obs.reset_index(drop=True)
    X = X.tocsr() if _is_sparse(X) else np.asarray(X)
    vals = X.data if _is_sparse(X) else X
    if not (np.all(vals >= 0) and np.all(vals == np.round(vals))):
        raise ValueError(f"{name}: the panel needs raw counts")
    return Background(X, genes, np.asarray(obs[spec["donor"]]).astype(str), name, obs)


def resolve_background(spec_path: Path, name: str, spec: dict, data_dir=None) -> dict:
    """A backgrounds.json entry with its local path: "path" (relative to the JSON) or "file" in
    `data_dir` (default: next to the JSON); the file must match the entry's "sha256" if given."""
    spec = dict(spec)
    if "path" in spec:
        p = Path(spec["path"])
        spec["path"] = str(p if p.is_absolute() else Path(spec_path).parent / p)
    else:
        spec["path"] = str(Path(data_dir or Path(spec_path).parent) / spec["file"])
    if spec.get("sha256") and sha256(Path(spec["path"])) != spec["sha256"]:
        raise ValueError(f"{name}: {spec['path']} does not match the sha256 in {Path(spec_path).name}")
    return spec


def load_backgrounds(spec_path, data_dir=None) -> dict:
    """The panel's backgrounds (B1, B2) of a backgrounds JSON (written by the rule of v1.md section
    3.1: per background the url, file, sha256, cell filter, donor column, counts source and
    gene-symbol column), loaded from `data_dir`, checked against its sha256 and planned
    (`plan_background`). The anchors' backgrounds (B3, B4) in the same file are not the panel's:
    anchors.py loads them itself."""
    spec = json.loads(Path(spec_path).read_text())
    bgs = {}
    for k, v in spec.items():
        if k not in BACKGROUNDS:
            continue
        bgs[k] = load_background(k, resolve_background(spec_path, k, v, data_dir))
        plan_background(bgs[k])
    return bgs


def level_of(mean: float, detection: float) -> str | None:
    """A gene's expression level: low - detected in a minority of cells (10% to 50%); medium -
    in a majority, with a mean below 2 counts per cell; high - in a majority, mean >= 2."""
    if 0.10 <= detection < 0.50:
        return "low"
    if detection >= 0.50:
        return "high" if mean >= 2.0 else "medium"
    return None


def _col_stats(X):
    """Mean count and detected share per gene, from exact integer sums (the same floats for a
    dense and a sparse matrix of the same counts)."""
    n = X.shape[0]
    if _is_sparse(X):
        X = X.tocsr()
        total = np.asarray(X.sum(axis=0, dtype=np.int64)).ravel()
        nnz = np.bincount(X.indices[X.data > 0], minlength=X.shape[1])
    else:
        total = np.asarray(X).sum(axis=0).astype(np.int64)
        nnz = (np.asarray(X) > 0).sum(axis=0)
    return total / n, nnz / n


def _dense(X) -> np.ndarray:
    return X.toarray() if _is_sparse(X) else np.asarray(X)


def norm_pearson_cols(a: np.ndarray, b: np.ndarray, tot: np.ndarray) -> float:
    """The engine's log-normalized Pearson (metrics.norm_pearson) of two count columns with the
    cells' totals: CP10k and log1p, over the cells where both are detected."""
    both = (a > 0) & (b > 0)
    if both.sum() < 3:
        return 0.0
    t = np.where(tot == 0, 1.0, tot)
    x, y = np.log1p(a / t * 1e4)[both], np.log1p(b / t * 1e4)[both]
    if x.std() == 0 or y.std() == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def _neighbourhood(means: np.ndarray, genes: list, target: int, exclude: set) -> list:
    """The genes closest to `target` in mean count among `genes` (GATE 5's rule, gates._neighbourhood)."""
    k = max(MATCHED_NEIGHBOURS_MIN, int(np.ceil(MATCHED_NEIGHBOURS_FRAC * len(genes))))
    pool = [g for g in genes if g != target and g not in exclude]
    dist = np.array([abs(means[g] - means[target]) for g in pool])
    return [pool[j] for j in np.argsort(dist, kind="mergesort")[:k]]


def _gate5_typicality(value, means, kept, pair, exclude, rng) -> float:
    """|z| of a candidate negative control's log-normalized Pearson against GATE 5's matched null
    on the plan's cells (pairs from the two genes' expression neighbourhoods among the kept genes,
    NEG_NULL_PAIRS of them drawn): the third review found the controls chosen by partial
    correlation alone failing GATE 5 systematically on some pairs (|z| 2 to 3.4, at 48 mice in 5 of
    5 datasets). `value(x, y)` is the metric of genes x and y on the plan's cells."""
    na = _neighbourhood(means, kept, pair[0], exclude)
    nb = _neighbourhood(means, kept, pair[1], exclude)
    pool = sorted({tuple(sorted((x, y))) for x in na for y in nb if x != y})
    pick = [pool[i] for i in rng.permutation(len(pool))[:NEG_NULL_PAIRS]]
    null = np.array([value(x, y) for x, y in pick])
    sd = float(null.std(ddof=1)) if len(null) > 1 else 0.0
    if not sd > 0:
        return float("inf")
    return float(round(abs(value(*pair) - float(np.median(null))) / sd, 10))


def _partial_corr(blocks, libs) -> np.ndarray:
    """Pooled within-donor correlation of log-normalized expression over all cells, with the log
    library size partialled out per donor: the co-variation of two genes that depth does not
    explain. `blocks` holds one donor's log-normalized values each, `libs` its log totals."""
    G = blocks[0].shape[1]
    S = np.zeros((G, G))
    for Y, L in zip(blocks, libs):
        Yc = Y - Y.mean(axis=0)
        Lc = L - L.mean()
        ll = float(Lc @ Lc)
        R = Yc - np.outer(Lc, (Lc @ Yc) / ll) if ll > 0 else Yc
        S += R.T @ R
    sd = np.sqrt(np.clip(np.diag(S), 0, None))
    den = np.outer(sd, sd)
    return np.divide(S, den, out=np.full_like(S, np.nan), where=den > 0)


def _level_search(X, names, mean, level, order, tot, plan_rows):
    """Per level: the candidate genes, their partial-correlation matrix and the pairs ranked by it."""
    out = {}
    for lvl in LEVELS:
        cand = [int(j) for j in order if level[j] == lvl][:LEVEL_CANDIDATES]
        lib = [np.where(tot[r] > 0, tot[r], 1.0) for r in plan_rows]
        blocks = [np.log1p(_dense(X[r][:, cand]) / t[:, None] * 1e4) for r, t in zip(plan_rows, lib)]
        # rounded so that the summation order of the matrix products cannot reorder near-ties
        C = (np.round(_partial_corr(blocks, [np.log(t) for t in lib]), 10) if len(cand) > 1
             else np.full((len(cand), len(cand)), np.nan))
        iu = np.triu_indices(len(cand), 1)
        r = C[iu]
        top = [t for t in np.argsort(-np.nan_to_num(r, nan=-np.inf), kind="mergesort") if np.isfinite(r[t])]
        out[lvl] = dict(cand=cand, C=C, iu=iu, r=r, top=top)
    return out


def _choose_pairs(search: dict, mean, names, k: int, typical=None):
    """The pool with k pairs per level (and the positive control) by the rule of
    `plan_background`; ValueError if some level cannot provide them. `typical(c, d, a, b,
    positive, li, j)` scores a candidate negative control (lower is more typical of GATE 5's
    null); without it the closest pair in mean expression is taken."""
    used: set = set()
    pool, positive, levels = [], None, {}
    for li, lvl in enumerate(LEVELS):
        cand, C, iu, r, top = (search[lvl][x] for x in ("cand", "C", "iu", "r", "top"))
        need = k + (1 if lvl == "high" else 0)
        if len(cand) < 4 * need:
            raise ValueError(f"level {lvl!r} has {len(cand)} eligible genes, too few for {k} pairs and "
                             "their negative controls")
        chosen = []
        for t in top:
            a, b = cand[iu[0][t]], cand[iu[1][t]]
            if a in used or b in used:
                continue
            chosen.append((a, b, float(r[t])))
            used.update((a, b))
            if len(chosen) == need:
                break
        if len(chosen) < need:
            raise ValueError(f"level {lvl!r}: fewer than {need} disjoint pairs")
        if lvl == "high":
            positive = chosen.pop()
        avg = (mean[cand][:, None] + mean[cand][None, :]) / 2
        for j, (a, b, rab) in enumerate(chosen):
            free = np.array([g not in used for g in cand])
            ok = np.triu((np.abs(np.nan_to_num(C, nan=np.inf)) < NEG_MAX_R) & np.outer(free, free), 1)
            gap = np.where(ok, np.abs(avg - (mean[a] + mean[b]) / 2), np.inf)
            if not np.isfinite(gap).any():
                raise ValueError(f"no negative control with |r| < {NEG_MAX_R} for {names[a]}-{names[b]}")
            flat = np.argsort(gap, axis=None, kind="mergesort")[:NEG_CANDIDATES]
            closest = [np.unravel_index(int(t), gap.shape) for t in flat if np.isfinite(gap.flat[int(t)])]
            if typical is None:
                (i, jj), z = closest[0], None
            else:
                scored = [(typical(cand[i], cand[jj], a, b, positive, li, j), rank, (i, jj))
                          for rank, (i, jj) in enumerate(closest)]
                z, _, (i, jj) = min(scored)
            c, d = cand[i], cand[jj]
            used.update((c, d))
            pool.append(dict(level=lvl, index=li * k + j, pair=[names[a], names[b]], r=rab,
                             neg_pair=[names[c], names[d]], neg_r=float(C[i, jj]),
                             neg_gate5_z=z, _genes=(a, b, c, d)))
        levels[lvl] = dict(candidates=len(cand), median_r=float(np.nanmedian(r)) if len(r) else float("nan"))
    return pool, positive, used, levels


def plan_background(bg: Background, seed: int = PLAN_SEED) -> dict:
    """Fixed before the key: the pair pool, its controls, the kept genes and N6c's random genes.

    Rule (v1.md section 3.1): only donors with >= CELLS_PER_DONOR cells; a gene is eligible for
    a pair if it is not mitochondrial or ribosomal and is detected in >= 5% of the cells of every
    donor; eligible genes are split into levels (`level_of`); in each level the 1,000 most
    expressed enter the search. A pair's correlation is the pooled within-donor correlation of
    log-normalized expression over all cells with the log library size partialled out
    (`_partial_corr`), on 200 cells per donor (public seed): what the two genes share beyond
    depth. Per level the k most correlated disjoint pairs form the pool, k = 8 where every level
    has them (each with its negative control) and otherwise the largest k >= 4 that every level
    has; the next disjoint pair of the high level is the background's positive control (GATE 5),
    as a housekeeping pair would be; each pool pair's negative control is, of the 10 unused pairs of
    its level closest to it in mean expression whose |partial correlation| is below 0.02, the one
    whose log-normalized Pearson is most typical of GATE 5's matched null on the plan's cells
    (`_gate5_typicality`, the third review)."""
    rng = np.random.default_rng(seed)
    donor_all = np.asarray(bg.donor).astype(str)
    vals, counts = np.unique(donor_all, return_counts=True)
    donors = vals[counts >= CELLS_PER_DONOR]
    rows = np.where(np.isin(donor_all, donors))[0]
    X = bg.X[rows]
    donor = donor_all[rows]
    names = [str(g) for g in bg.genes]
    upper = [g.upper() for g in names]
    mean, det = _col_stats(X)
    donor_det = np.full(len(names), np.inf)
    for d in donors:
        donor_det = np.minimum(donor_det, _col_stats(X[donor == d])[1])
    skip = np.array([u.startswith(("MT-", "RPS", "RPL")) for u in upper])
    eligible = ~skip & (donor_det >= MIN_DONOR_DETECTION)
    level = np.array([level_of(m, f) if ok else None for m, f, ok in zip(mean, det, eligible)], dtype=object)
    order = np.argsort(-mean, kind="mergesort")
    tot = (np.asarray(X.sum(axis=1, dtype=np.int64)).ravel() if _is_sparse(X)
           else np.asarray(X).sum(axis=1).astype(np.int64)).astype(float)
    plan_rows = [np.sort(rng.choice(np.where(donor == d)[0], size=min(PLAN_CELLS_PER_DONOR, int((donor == d).sum())),
                                    replace=False)) for d in donors]
    search = _level_search(X, names, mean, level, order, tot, plan_rows)
    # GATE 5's view of a candidate negative control: on the plan's cells, with the totals and the
    # expression neighbourhoods of the genes a dataset keeps (the most expressed; the named genes
    # join them)
    cells = np.concatenate(plan_rows)
    Xc = X[cells]
    top = [int(j) for j in order[:N_GENES]]
    if _is_sparse(Xc):
        Xc = Xc.tocsc()
        Xc.eliminate_zeros()
        tot_kept = np.asarray(Xc[:, top].sum(axis=1, dtype=np.int64)).ravel().astype(float)
    else:
        Xc = np.asarray(Xc)
        tot_kept = Xc[:, top].sum(axis=1).astype(np.int64).astype(float)
    tot_kept = np.where(tot_kept == 0, 1.0, tot_kept)

    def detected(j):  # the cells where gene j is detected (sorted) and its counts there
        if _is_sparse(Xc):
            lo, hi = Xc.indptr[j], Xc.indptr[j + 1]
            return Xc.indices[lo:hi], Xc.data[lo:hi].astype(float)
        idx = np.nonzero(Xc[:, j] > 0)[0]
        return idx, Xc[idx, j].astype(float)

    def value(x, y):  # norm_pearson_cols on the plan's cells, from the co-detected entries only
        ix, vx = detected(x)
        iy, vy = detected(y)
        common, px, py = np.intersect1d(ix, iy, assume_unique=True, return_indices=True)
        return norm_pearson_cols(vx[px], vy[py], tot_kept[common])

    def typical(c, d, a, b, positive, li, j):
        exclude = {a, b, c, d} | ({positive[0], positive[1]} if positive else set())
        # a dataset's N_GENES kept genes as far as they are known here: the pair, the candidate and
        # the most expressed (the fourth review: the most expressed N_GENES and the four came to up to
        # 2,004, so GATE 5's neighbourhoods, 5% of the genes, to 101 genes, where the engine's are 100)
        kept = sorted(list(dict.fromkeys([a, b, c, d, *top]))[:N_GENES])
        z_rng = np.random.default_rng([seed, li, j, c, d])
        return _gate5_typicality(value, mean, kept, (c, d), exclude, z_rng)
    tried = []
    for k in range(MAX_PAIRS_PER_LEVEL, MIN_PAIRS_PER_LEVEL - 1, -1):
        try:
            pool, positive, used, levels = _choose_pairs(search, mean, names, k, typical)
            break
        except ValueError as exc:
            tried.append(f"{k}: {exc}")
    else:
        raise ValueError(f"{bg.name}: no pool of >= {MIN_PAIRS_PER_LEVEL} pairs per level ({'; '.join(tried)})")
    for pe in pool:
        a, b, c, d = pe.pop("_genes")
        pe.update(mean=[float(mean[a]), float(mean[b])], detection=[float(det[a]), float(det[b])],
                  neg_mean=[float(mean[c]), float(mean[d])])
    for lvl in LEVELS:
        levels[lvl]["eligible_genes"] = int(sum(level == lvl))
    pos_pair = [names[positive[0]], names[positive[1]]]
    for pe in pool:
        pe["pos_pair"] = pos_pair
    g2m = [g for g in G2M_GENES if g in upper]
    named = set(used) | {upper.index(g) for g in g2m}
    keep = list(dict.fromkeys([*sorted(named), *[int(j) for j in order if j not in named]]))
    keep = sorted(keep[:max(N_GENES, len(named))])
    free = [j for j in keep if j not in named]
    random_genes = sorted(rng.choice(free, size=min(RANDOM_SCORE_GENES, len(free)), replace=False).tolist())
    plan = dict(rule=dict(cells_per_donor=CELLS_PER_DONOR, min_donor_detection=MIN_DONOR_DETECTION,
                          level_rule="low: detected in [10%, 50%); medium: >= 50%, mean < 2; high: >= 50%, mean >= 2",
                          level_candidates=LEVEL_CANDIDATES, plan_cells_per_donor=PLAN_CELLS_PER_DONOR,
                          neg_max_r=NEG_MAX_R, max_pairs_per_level=MAX_PAIRS_PER_LEVEL,
                          min_pairs_per_level=MIN_PAIRS_PER_LEVEL, seed=seed),
                pairs_per_level=len(pool) // len(LEVELS), fewer_pairs_because=tried,
                donors=[str(d) for d in donors], cells=int(len(rows)), genes=[names[j] for j in keep],
                levels=levels, pool=pool,
                positive_control=dict(pair=pos_pair, r=positive[2],
                                      mean=[float(mean[positive[0]]), float(mean[positive[1]])]),
                g2m=[names[upper.index(g)] for g in g2m], random_genes=[names[j] for j in random_genes])
    kept = X[:, keep]
    bg.X = (kept.astype(np.int32).toarray() if _is_sparse(kept) else np.asarray(kept, dtype=np.int32))
    bg.genes = [names[j] for j in keep]
    bg.donor = donor
    if bg.extra_obs is not None and len(bg.extra_obs) == len(donor_all):
        bg.extra_obs = bg.extra_obs.iloc[rows].reset_index(drop=True)
    bg.plan = plan
    return plan


def background_sha256(bg: Background) -> str:
    """Canonical sha256 of a planned background: its counts, genes, donors and plan."""
    h = hashlib.sha256(b"metric-autopsy panel background v1\n")
    h.update(json.dumps(dict(name=bg.name, shape=list(bg.X.shape), genes=bg.genes, plan=bg.plan),
                        sort_keys=True).encode())
    h.update(np.ascontiguousarray(bg.X, dtype="<i4").tobytes())
    h.update("\x1f".join(map(str, bg.donor)).encode())
    return h.hexdigest()


def save_compact(bg: Background, path: Path) -> None:
    """A planned background as one .npz (counts of the kept genes, donors, plan)."""
    np.savez_compressed(path, X=np.asarray(bg.X, dtype=np.int32), genes=np.asarray(bg.genes),
                        donor=np.asarray(bg.donor).astype(str), plan=np.asarray(json.dumps(bg.plan)),
                        name=np.asarray(bg.name))


def load_compact(path: Path) -> Background:
    z = np.load(path, allow_pickle=False)
    return Background(z["X"], [str(g) for g in z["genes"]], z["donor"].astype(str), str(z["name"]),
                      plan=json.loads(str(z["plan"])))


# --------------------------------------------------------------------------- #
# truth generators: binomial thinning (seqgendiff-style), this module's own code
# --------------------------------------------------------------------------- #
def thin(X: np.ndarray, p, rng) -> np.ndarray:
    """Binomial thinning: every molecule kept with probability p (scalar or broadcastable)."""
    return rng.binomial(np.round(X).astype(np.int64), np.clip(p, 0.0, 1.0)).astype(np.float64)


def drop_out(X: np.ndarray, frac: float, rng) -> np.ndarray:
    """Zero a fraction `frac` of the detected entries (extra dropout)."""
    out = X.copy()
    nz = np.argwhere(out > 0)
    k = int(round(frac * len(nz)))
    if k:
        pick = nz[rng.choice(len(nz), size=k, replace=False)]
        out[pick[:, 0], pick[:, 1]] = 0.0
    return out


def variable_capture(X: np.ndarray, rng, beta=N8_BETA) -> np.ndarray:
    """N8: every cell keeps its molecules with its own probability, drawn from Beta(a, b)."""
    return thin(X, rng.beta(beta[0], beta[1], size=X.shape[0])[:, None], rng)


def _keep(z: np.ndarray, dose: float) -> np.ndarray:
    """Keep probability 1 / (1 + exp(-dose z)): mean 0.5 at every dose, spread growing with it."""
    return 1.0 / (1.0 + np.exp(-dose * z))


def inject_coupling(X: np.ndarray, ia: int, ib: int, dose: float, rng) -> np.ndarray:
    """Couple genes a and b by thinning both with one per-cell keep probability driven by a
    shared z_i ~ N(0, 1): cells with a high z keep more of both."""
    out = X.copy()
    p = _keep(rng.standard_normal(X.shape[0]), dose)
    out[:, ia] = thin(X[:, ia], p, rng)
    out[:, ib] = thin(X[:, ib], p, rng)
    return out


def sham_coupling(X: np.ndarray, ia: int, ib: int, dose: float, rng) -> np.ndarray:
    """The same marginal thinning of a and b with an independent z per gene: no coupling."""
    out = X.copy()
    for j in (ia, ib):
        out[:, j] = thin(X[:, j], _keep(rng.standard_normal(X.shape[0]), dose), rng)
    return out


# --------------------------------------------------------------------------- #
# one dataset
# --------------------------------------------------------------------------- #
def _draw_cells(bg: Background, donors, rng):
    rows, donor_of = [], []
    for d in donors:
        idx = np.where(bg.donor == d)[0]
        pick = rng.choice(idx, size=min(CELLS_PER_DONOR, len(idx)), replace=False)
        rows.append(pick)
        donor_of += [d] * len(pick)
    rows = np.concatenate(rows)
    return np.asarray(bg.X[rows], dtype=np.float64), np.asarray(donor_of)


def _two_groups(bg: Background, rng, per_group: int, total: int | None = None):
    """`total` donors (default 2 x per_group) drawn without replacement, the first `per_group` of a
    permutation in A and the rest in B."""
    donors = rng.choice(bg.donors, size=total or 2 * per_group, replace=False)
    X, donor = _draw_cells(bg, donors, rng)
    group_of = {d: ("A" if k < per_group else "B") for k, d in enumerate(rng.permutation(donors))}
    return X, donor, np.array([group_of[d] for d in donor])


def pool_entry(bg: Background, entry: dict) -> dict:
    pool = bg.plan["pool"]
    if not 0 <= int(entry["pair"]) < len(pool):
        raise ValueError(f"{entry.get('id')}: pair {entry['pair']} outside {bg.name}'s pool of {len(pool)}")
    return pool[int(entry["pair"])]


def e_dose(level: str, pilot: dict) -> float:
    """The base dose of the real effects at a level of B1: the key dose where the pilot found one,
    else the level's saturation dose on B1 (decided 2026-10-08)."""
    return float(pilot["e_dose"][level])


def build(entry: dict, bgs: dict, pilot: dict) -> tuple[np.ndarray, pd.DataFrame, list, list[dict]]:
    """Counts, obs, genes and claim card(s) of one dataset. `pilot` holds, per background and
    level, the SESOI and the saturation dose (the GATE 4 injection of the card), the base dose of
    the real effects at each level of B1, and Δ* per pair (whose sign is the true direction)."""
    cond = conditions()[entry["condition"]]
    bg = bgs[cond.background]
    rng = np.random.default_rng(entry["seed"])
    pe = pool_entry(bg, entry)
    level = pe["level"]
    genes = list(bg.genes)
    ia, ib = genes.index(pe["pair"][0]), genes.index(pe["pair"][1])
    side, other = entry["side"], ("B" if entry["side"] == "A" else "A")
    variant = entry["variant"]
    if cond.name == "N5":
        donors = rng.choice(bg.donors, size=DONORS_PER_GROUP, replace=False)
        X, donor = _draw_cells(bg, donors, rng)
        group = rng.choice(np.array(["A", "B"]), size=len(donor))
    elif cond.name == "N4":
        X, donor, group = _two_groups(bg, rng, 3)
    elif cond.name == "N7":  # N7_MAX_MICE mice (all where B2 has fewer); with an odd number B has one more
        mice = min(len(bg.donors), N7_MAX_MICE)
        X, donor, group = _two_groups(bg, rng, mice // 2, total=mice)
    else:
        X, donor, group = _two_groups(bg, rng, DONORS_PER_GROUP)
    in_side = group == side
    if cond.name == "N2":
        X[in_side] = thin(X[in_side], float(variant.split("=")[1]), rng)
    elif cond.name == "N3":
        X[in_side] = drop_out(X[in_side], float(variant.split("=")[1]), rng)
    elif cond.name == "N8":
        X[in_side] = variable_capture(X[in_side], rng)
    elif cond.data == "effect":
        factor = e1_factor(variant) if cond.name == "E1" else 1.0
        dose = factor * e_dose(level, pilot)
        X[in_side] = inject_coupling(X[in_side], ia, ib, dose, rng)
        X[~in_side] = sham_coupling(X[~in_side], ia, ib, dose, rng)
        if cond.name == "E2":   # artifact against the effect: capture loss where the signal is
            X[in_side] = thin(X[in_side], E_CAPTURE, rng)
        elif cond.name == "E3":  # artifact with the effect: capture loss on the other side
            X[~in_side] = thin(X[~in_side], E_CAPTURE, rng)
    obs = pd.DataFrame({"group": group, "donor": donor})
    obs["total_counts"] = X.sum(axis=1)
    obs["n_genes_by_counts"] = (X > 0).sum(axis=1)
    # the claim's direction: for a real effect the true one (the signal side is higher where
    # Δ* > 0; an effect dataset whose Δ* is unknown is refused, not given the signal side's
    # direction: the third review); for a null a coin, so that the card does not reveal the condition
    if cond.data == "effect":
        try:
            positive = float(delta_of(cond.name, variant, entry["pair"], pilot)["value"]) >= 0
        except (KeyError, TypeError) as exc:
            raise ValueError(f"{entry.get('id')}: Δ* of pair {entry['pair']} ({cond.name} {variant}) is unknown; "
                             "the claim's direction is its sign") from exc
        higher = side if positive else other
    else:
        higher = ("A", "B")[int(rng.integers(2))]
    direction = "decrease" if higher == "A" else "increase"  # change from A to B
    sesoi = float(pilot["sesoi"][cond.background][level])
    strength = float(pilot["saturation_dose"][cond.background][level])
    base = dict(id=entry["id"], background=cond.background, metric=cond.metric, gene_pair=pe["pair"],
                pos_pair=pe["pos_pair"], neg_pair=pe["neg_pair"],
                score_genes=bg.plan["random_genes"] if cond.metric == "score" else None,
                group_col="group", groups=["A", "B"], replicate_col="donor",
                signal_test=(dict(kind="module", genes=bg.plan["g2m"]) if cond.metric == "score"
                             else dict(kind="coupling", genes=pe["pair"], strength=strength)),
                prereg=dict(estimand="composition", direction=direction, sesoi=sesoi, sesoi_scale="observed",
                            delta_min=DELTA_MIN_FRACTION * sesoi, judgment_pending=False))
    cards = [base]
    if cond.cards == 2:  # N4: the same cells analysed without and with the replicate unit
        cards = [dict(base, id=f"{entry['id']}a", replicate_col=None), dict(base, id=f"{entry['id']}b")]
    return X, obs, genes, cards


# --------------------------------------------------------------------------- #
# canonical hashes and files
# --------------------------------------------------------------------------- #
def dataset_sha256(X: np.ndarray, obs: pd.DataFrame, genes: list) -> str:
    """Canonical sha256 of a dataset (counts as little-endian int32, genes, obs columns); the
    same bytes on any machine, so a revealed key reproduces it."""
    h = hashlib.sha256(b"metric-autopsy panel dataset v1\n")
    Xi = np.ascontiguousarray(np.round(np.asarray(X)), dtype="<i4")
    h.update(json.dumps(dict(shape=list(Xi.shape), genes=list(genes), obs=list(obs.columns))).encode())
    h.update(Xi.tobytes())
    for c in obs.columns:
        v = np.asarray(obs[c])
        if np.issubdtype(v.dtype, np.number):
            h.update(np.ascontiguousarray(v, dtype="<f8").tobytes())
        else:
            h.update("\x1f".join(map(str, v)).encode())
    return h.hexdigest()


def card_sha256(card: dict) -> str:
    return hashlib.sha256(json.dumps(card, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def save_npz(path: Path, X, obs: pd.DataFrame, genes: list) -> None:
    """Plain arrays only (strings as fixed-width unicode), loadable without pickle; a sparse X is
    stored as CSR arrays."""
    arrays = {}
    for c in obs.columns:
        v = np.asarray(obs[c])
        arrays[f"obs_{c}"] = v.astype(str) if v.dtype == object else v
    if _is_sparse(X):
        X = X.tocsr()
        arrays.update(X_data=X.data.astype(np.int32), X_indices=X.indices, X_indptr=X.indptr,
                      X_shape=np.asarray(X.shape))
        np.savez_compressed(path, genes=np.asarray(genes), **arrays)
    else:
        np.savez_compressed(path, X=np.asarray(X, dtype=np.int32), genes=np.asarray(genes), **arrays)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
