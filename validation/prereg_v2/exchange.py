"""Validation v2's exchange: what the gates add over a comparator, and the intervals of every rate
(v2.md sections 5 and 6; `validation/plans/v2_requirements.md` section 1, the owner's decisions of
2026-10-11, 2026-10-12 and 2026-10-17).

Every analyzer's outcome on a claim card — the engine's report, (b) the test over donors, (c) the engine
with GATE 0, 2, 4 and 5 off — is first read as one of five calls about the claim (`CALLS`):

- ``effect``  an effect claimed in the claim's direction (SUPPORTED; (b)'s detected effect);
- ``against`` an effect claimed against it (NOT SUPPORTED with ``opposite_direction``);
- ``absent``  an absence claimed (NO DETECTABLE EFFECT);
- ``below``   a detected effect shown smaller than the SESOI (EFFECT BELOW SESOI, v0.4 change 14);
- ``none``    no claim about the data (INCONCLUSIVE, metric invalid, explained by depth, ...).

Claims in v2 are directional, and on an effect card the claim names the true direction (the mirrored
wording restates it), so "against" is against the truth there.

The card is then scored by its truth (`Truth`): the data (null, or a real effect), whether the metric is
at risk (valid or ambiguous by v2.md section 3.3's truth), and on an effect card whether |Δ*| >= SESOI.
A **false conclusion**: an effect claimed in either direction on null data, or against the true direction
on a real effect; an absence claimed on a real effect with |Δ*| >= SESOI or with a metric that is not at
risk; EFFECT BELOW SESOI on a real effect with |Δ*| >= SESOI (change 14: on null data it is no false
effect, below the SESOI it is correct). An absence claimed on null data with a metric not at risk is true
of the data but not shown by the metric: an **unwarranted absence**, counted apart. A **detection**: an
effect claimed in the true direction on a real effect; EFFECT BELOW SESOI is not one.

The exchange (the owner, 2026-10-12): one false conclusion is worth four lost detections. With a share pi
of true effects among the claims, the gates pay off against a comparator where
(1 - pi) r + pi r1 >= pi l / 4, r and r1 the false conclusions they remove per null and per effect card,
l the detections they lose per effect card; pi* = r / (r + l/4 - r1), the largest share at which they pay
off (`pi_star`, the formula of `validation/exploratory/v1_correction_rule.py`).

The intervals: π* by percentile bootstrap over cards (null and effect cards resampled apart) and over
clusters (gene pairs; all of a cluster's cards together); a rate by Clopper-Pearson at its effective size
under the design effect by donor, by pair or by both (`overlap_deff`, `joint_deff`: the estimators of
`validation/prereg/score.py` and `validation/exploratory/v1_pair_deff.py`, requirement 4).

Imports numpy and this directory's `num` only, never the engine.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from num import clopper_pearson

CALLS = ("effect", "against", "absent", "below", "none")

# the engine's causes (report.CAUSES of v0.4) that make a call, and the label each must carry
_CALL_OF_CAUSE = {
    "replicated": ("effect", "SUPPORTED"),
    "provisional": ("effect", "SUPPORTED"),
    "opposite_direction": ("against", "NOT SUPPORTED"),
    "no_detectable_effect": ("absent", "NO DETECTABLE EFFECT"),
    "effect_below_sesoi": ("below", "EFFECT BELOW SESOI"),
}


def label_of(verdict: str) -> str:
    """A verdict's label: its text before " — ", without a parenthesis ("SUPPORTED (provisional ...)")."""
    head = str(verdict).split(" — ", 1)[0].strip()
    return head.split(" (", 1)[0].strip()


def engine_call(verdict: str, cause: str) -> str:
    """The call of an engine report's verdict and cause (`report.decide_cause`). A cause that makes a call
    must come with its label; a mismatch raises (it is a rule violation, scored apart, not a call)."""
    if cause in _CALL_OF_CAUSE:
        call, label = _CALL_OF_CAUSE[cause]
        if label_of(verdict) != label:
            raise ValueError(f"the cause {cause!r} with the label {label_of(verdict)!r}")
        return call
    return "none"


@dataclass(frozen=True)
class Truth:
    """A card's truth: ``data`` "null" or "effect"; ``at_risk``, the metric valid or ambiguous on the
    card's pool item; ``big``, on an effect card |Δ*| >= SESOI."""

    data: str
    at_risk: bool
    big: bool = False

    def __post_init__(self):
        if self.data not in ("null", "effect"):
            raise ValueError(f"data must be 'null' or 'effect', not {self.data!r}")
        if self.data == "null" and self.big:
            raise ValueError("a null card has no effect to be big")


def kinds(call: str, truth: Truth) -> tuple[bool, bool, bool]:
    """(a false effect, a false absence, an unwarranted absence) for one call on one card."""
    if call not in CALLS:
        raise ValueError(f"unknown call {call!r}")
    if truth.data == "null":
        return call in ("effect", "against"), False, call == "absent" and not truth.at_risk
    false_absence = ((call == "absent" and (truth.big or not truth.at_risk))
                     or (call == "below" and truth.big))
    return call == "against", false_absence, False


def wrong(call: str, truth: Truth) -> bool:
    """A false conclusion (a false effect or a false absence)."""
    fe, fa, _ = kinds(call, truth)
    return fe or fa


def detection(call: str, truth: Truth) -> bool:
    """An effect claimed in the true direction on a real effect."""
    return truth.data == "effect" and call == "effect"


def pi_star(r0: float, r1: float, lost: float) -> float:
    """The largest share of true effects at which the gates pay off, (1 - pi) r0 + pi r1 >= pi lost / 4:
    1.0 where they pay off at every share, 0.0 where at none, nan where only above a share."""
    a, b = r0, lost / 4 - r1
    if a + b > 0:
        return min(max(a / (a + b), 0.0), 1.0)
    return 1.0 if a >= 0 else float("nan")


def _ci(values) -> dict:
    v = np.asarray(values, float)
    undefined = int(np.isnan(v).sum())
    if undefined == len(v):
        return dict(lo=float("nan"), hi=float("nan"), undefined=undefined)
    lo, hi = np.nanpercentile(v, [2.5, 97.5])
    return dict(lo=float(lo), hi=float(hi), undefined=undefined)


def exchange(ref_calls, calls, truths, clusters, *, rng, n_boot: int = 2000) -> dict:
    """π* of an analyzer (`calls`) against a comparator (`ref_calls`) on the same cards, with its 95%
    percentile intervals over cards and over clusters (`clusters`: a card's gene pair, or any unit whose
    cards go together). The bootstrap draws as `v1_correction_rule.py` section 4 drew: first `n_boot`
    resamples of the null and the effect cards, then `n_boot` of the clusters, from `rng`."""
    truths = list(truths)
    n = len(truths)
    if not (len(ref_calls) == len(calls) == len(clusters) == n):
        raise ValueError("ref_calls, calls, truths and clusters must have one entry per card")
    is_null = np.array([t.data == "null" for t in truths])
    n0, n1 = int(is_null.sum()), int((~is_null).sum())
    if n0 == 0 or n1 == 0:
        raise ValueError("the exchange needs null cards and effect cards")
    dE = np.array([float(wrong(a, t)) - float(wrong(b, t)) for a, b, t in zip(ref_calls, calls, truths)])
    dS = np.array([float(detection(a, t)) - float(detection(b, t)) for a, b, t in zip(ref_calls, calls, truths)])
    i0, i1 = np.flatnonzero(is_null), np.flatnonzero(~is_null)
    e0, e1, ls = int(dE[i0].sum()), int(dE[i1].sum()), int(dS[i1].sum())
    est = pi_star(e0 / n0, e1 / n1, ls / n1)
    card = []
    for _ in range(n_boot):
        j0, j1 = rng.choice(i0, n0), rng.choice(i1, n1)
        card.append(pi_star(dE[j0].sum() / n0, dE[j1].sum() / n1, dS[j1].sum() / n1))
    units = sorted({str(c) for c in clusters})
    at = {u: i for i, u in enumerate(units)}
    uidx = np.array([at[str(c)] for c in clusters])
    agg = np.zeros((len(units), 5))
    np.add.at(agg, uidx, np.column_stack([is_null, ~is_null, dE * is_null, dE * ~is_null, dS]))
    clus = []
    for _ in range(n_boot):
        t = np.bincount(rng.integers(0, len(units), len(units)), minlength=len(units)) @ agg
        clus.append(pi_star(t[2] / t[0], t[3] / t[1], t[4] / t[1]) if t[0] and t[1] else float("nan"))
    share = n1 / n
    return dict(null_cards=n0, effect_cards=n1, clusters=len(units),
                removed_null=e0, removed_effect=e1, lost=ls,
                r=e0 / n0, r1=e1 / n1, l=ls / n1, pi_star=est,
                ci_cards=_ci(card), ci_clusters=_ci(clus),
                share=share, pays_off_at_share=bool(e0 + e1 >= ls / 4))


# --------------------------------------------------------------------------- #
# design effects and the intervals of rates (requirement 4)
# --------------------------------------------------------------------------- #
def _units(sets: list) -> np.ndarray:
    """The cards x units incidence matrix (float64) of each card's set of units."""
    names = sorted({u for st in sets for u in st})
    at = {u: i for i, u in enumerate(names)}
    M = np.zeros((len(sets), len(names)))
    for j, st in enumerate(sets):
        M[j, [at[u] for u in st]] = 1.0
    return M


def _moments(M: np.ndarray, e: np.ndarray) -> tuple[float, float, float]:
    """For O = M M^T without its diagonal: (sum of O, sum of O^2, e' O e), through M'M."""
    S = M.sum(axis=1)
    G = M.T @ M
    tot = M.sum(axis=0)
    return (float(tot @ tot - S.sum()), float((G * G).sum() - (S * S).sum()),
            float((M.T @ e) @ (M.T @ e) - (e * e * S).sum()))


def overlap_deff(y, sets: list) -> dict:
    """The design effect of a rate whose cards share units (donors; or a gene pair, one unit a card):
    Cov(y_j, y_k) = beta x |S_j & S_k|, beta the moment estimate over all pairs of cards, truncated at 0
    (`score.overlap_interval`'s estimator, computed through M'M as `v1_pair_deff.overlap_deff`)."""
    y = np.asarray(y, float)
    n = len(y)
    p = float(y.mean())
    vb = p * (1 - p) / n
    if vb == 0:
        return dict(deff=1.0, beta=0.0)
    e = y - p
    s_o, s_o2, eoe = _moments(_units(sets), e)
    beta = max(0.0, eoe / s_o2) if s_o2 > 0 else 0.0
    var = (float(e @ e) + beta * s_o) / n ** 2
    return dict(deff=max(1.0, var / vb), beta=beta)


def joint_deff(y, donor_sets: list, units: list) -> dict:
    """The design effect by donor and by pair together: Cov(y_j, y_k) = b_d |S_j & S_k| + b_p [unit_j ==
    unit_k] for j != k, least squares on the residuals' cross-products, each b >= 0
    (`v1_pair_deff.joint_deff`)."""
    y = np.asarray(y, float)
    n = len(y)
    p = float(y.mean())
    vb = p * (1 - p) / n
    if vb == 0:
        return dict(deff=1.0, beta_donor=0.0, beta_pair=0.0)
    e = y - p
    Md, Mp = _units(donor_sets), _units([[u] for u in units])
    sd, a11, r1 = _moments(Md, e)
    sp, a22, r2 = _moments(Mp, e)
    a12 = 0.0
    S = Md.sum(axis=1)
    for c in range(Mp.shape[1]):
        rows = Mp[:, c] > 0
        t = Md[rows].sum(axis=0)
        a12 += float(t @ t - S[rows].sum())
    det = a11 * a22 - a12 * a12
    bd, bp = ((r1 * a22 - r2 * a12) / det, (r2 * a11 - r1 * a12) / det) if det > 0 else (0.0, 0.0)
    if bd < 0 or bp < 0:
        cand = [(max(0.0, r1 / a11) if a11 > 0 else 0.0, 0.0), (0.0, max(0.0, r2 / a22) if a22 > 0 else 0.0)]

        def sse(b):
            return -2 * (b[0] * r1 + b[1] * r2) + b[0] ** 2 * a11 + b[1] ** 2 * a22 + 2 * b[0] * b[1] * a12
        bd, bp = min(cand, key=sse)
    var = (float(e @ e) + bd * sd + bp * sp) / n ** 2
    return dict(deff=max(1.0, var / vb), beta_donor=float(bd), beta_pair=float(bp))


def rate(y, *, donor_sets: list | None = None, units: list | None = None, alpha: float = 0.05) -> dict:
    """A rate k/n with its Clopper-Pearson interval, and its design effects by donor, by pair and by both
    where the cards' donors and pairs are given; the interval at the largest of them is the
    Korn-Graubard interval at the effective size n / deff."""
    y = np.asarray(y, float)
    n = len(y)
    if n == 0:
        return dict(k=0, n=0, rate=float("nan"), ci95=[float("nan")] * 2)
    k = float(y.sum())
    out = dict(k=int(k), n=n, rate=k / n, ci95=list(clopper_pearson(k, n, alpha)))
    deffs = {}
    if donor_sets is not None:
        deffs["donor"] = overlap_deff(y, donor_sets)
    if units is not None:
        deffs["pair"] = overlap_deff(y, [[u] for u in units])
    if donor_sets is not None and units is not None:
        deffs["both"] = joint_deff(y, donor_sets, units)
    if deffs:
        worst = max(d["deff"] for d in deffs.values())
        n_eff = n / worst
        out.update(deff=deffs, deff_max=worst, n_eff=n_eff,
                   ci95_deff=list(clopper_pearson(out["rate"] * n_eff, n_eff, alpha)))
    return out
