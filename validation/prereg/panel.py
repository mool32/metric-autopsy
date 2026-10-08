"""The blinded panel of the confirmatory validation (validation/prereg/v1.md, section 3).

    python validation/prereg/panel.py --key-seed SEED --backgrounds backgrounds.json \
        --pilot pilot.json --out DIR

The key seed (held by the project owner, never in the repository) decides every dataset's
condition, variant (dose or step), side and seed, and the order of the dataset IDs. For each ID
the script writes ``DIR/<ID>.npz`` (counts, obs, genes) and ``DIR/<ID>.json`` (the claim card,
which does not depend on the condition); ``DIR/manifest.json`` lists the IDs with the sha256 of
both files. ``score.py`` re-derives the assignment from the same seed after the engine has run.

Independence: the truth is generated here, with this module's own binomial thinning (as in
seqgendiff, Gerard 2020), never with the engine's modules; this file does not import
``metric_autopsy``. ``backgrounds.json`` maps "B1" and "B2" to {"path": .h5ad or .npz, "donor":
obs column, "filter": {obs column: value}}; ``pilot.json`` holds the SESOI, the dose ladder and
the establishable conditions (written by ``oracle.py pilot``).
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

KEY_N = 600            # datasets per key null condition (oc.log: S1 needs >= 506)
DONORS_PER_GROUP = 8   # N1-N3, N6, E1-E3: 2 x 8 donors from B1
CELLS_PER_DONOR = 200  # cells drawn per donor and dataset
N_GENES = 2000         # genes kept per background: the most expressed, plus the named genes
RANDOM_SCORE_GENES = 50

# Tirosh et al. 2016 G2/M genes (human symbols); the module the random-gene score claims to
# measure (N6c). Genes absent from a background are dropped.
G2M_GENES = (
    "HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 "
    "PIMREG SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 JPT1 "
    "CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN "
    "LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA").split()

# --------------------------------------------------------------------------- #
# the design: every condition with its variants and number of datasets
# --------------------------------------------------------------------------- #
NULL_ALLOWED = ("NO DETECTABLE EFFECT", "INCONCLUSIVE", "NOT SUPPORTED")
EFFECT_ALLOWED = ("SUPPORTED", "INCONCLUSIVE")


@dataclass(frozen=True)
class Condition:
    name: str
    background: str
    variants: tuple            # ((variant label, number of datasets), ...)
    allowed: tuple
    definite: str | None       # expected verdict where the oracle finds it establishable
    null: bool                 # scored for false SUPPORTED
    key: bool = False          # key null condition (criterion S1), at its key variant
    key_variant: str | None = None
    cards: int = 1             # claim cards per dataset (N4: without and with a replicate column)


CONDITIONS = (
    Condition("N1", "B1", (("null", KEY_N),), NULL_ALLOWED, "NO DETECTABLE EFFECT", True, True, "null"),
    Condition("N2", "B1", (("c=0.5", KEY_N), ("c=0.9", 100), ("c=0.7", 100), ("c=0.3", 100)),
              NULL_ALLOWED, "NOT SUPPORTED", True, True, "c=0.5"),
    Condition("N3", "B1", (("f=0.1", 100), ("f=0.2", 100), ("f=0.4", 100)),
              NULL_ALLOWED, "NOT SUPPORTED", True),
    Condition("N4", "B1", (("3v3", 100),), ("INCONCLUSIVE", "NO DETECTABLE EFFECT", "NOT SUPPORTED"),
              None, True, cards=2),
    Condition("N5", "B1", (("sham", KEY_N),), ("NO DETECTABLE EFFECT", "INCONCLUSIVE"),
              "NO DETECTABLE EFFECT", True, True, "sham"),
    Condition("N6a", "B1", (("random", 100),), ("NOT SUPPORTED", "INCONCLUSIVE"), "NOT SUPPORTED", True),
    Condition("N6b", "B1", (("constant", 50),), ("DEGENERATE METRIC",), "DEGENERATE METRIC", True),
    Condition("N6c", "B1", (("random-genes", KEY_N),), ("NOT SUPPORTED", "INCONCLUSIVE"),
              "NOT SUPPORTED", True, True, "random-genes"),
    Condition("N7", "B2", (("mice", 200),), NULL_ALLOWED, "NO DETECTABLE EFFECT", True),
    Condition("E1", "B1", (("dose=key", 200), ("dose=0.25", 100), ("dose=0.5", 100), ("dose=1.5", 100)),
              EFFECT_ALLOWED, "SUPPORTED", False),
    Condition("E2", "B1", (("against", 200),), EFFECT_ALLOWED, "SUPPORTED", False),
    Condition("E3", "B1", (("with", 200),), EFFECT_ALLOWED, "SUPPORTED", False),
)


def conditions() -> dict:
    return {c.name: c for c in CONDITIONS}


def n_datasets() -> int:
    return sum(n for c in CONDITIONS for _, n in c.variants)


# --------------------------------------------------------------------------- #
# the key: condition, variant, side and seed of every dataset ID
# --------------------------------------------------------------------------- #
def assign(key_seed: int) -> list[dict]:
    """Every dataset of the design, in an order and with seeds decided by the key seed.

    Deterministic in the seed; nothing about an entry is visible in its ID."""
    entries = [dict(condition=c.name, variant=v, index=i)
               for c in CONDITIONS for v, n in c.variants for i in range(n)]
    ss = np.random.SeedSequence(int(key_seed))
    order_seq, *data_seqs = ss.spawn(len(entries) + 1)
    rng = np.random.default_rng(order_seq)
    perm = rng.permutation(len(entries))
    out = []
    width = len(str(len(entries)))
    for rank, k in enumerate(perm):
        e = dict(entries[k])
        side_rng = np.random.default_rng(data_seqs[k].spawn(1)[0])
        e.update(id=f"D{rank + 1:0{width}d}", side=("A", "B")[int(side_rng.integers(2))],
                 seed=int(data_seqs[k].generate_state(1)[0]))
        out.append(e)
    return out


# --------------------------------------------------------------------------- #
# backgrounds: genes, analysed pair and controls, fixed before the key
# --------------------------------------------------------------------------- #
@dataclass
class Background:
    """Raw counts of one background (cells x genes) with the donor of every cell."""
    X: np.ndarray
    genes: list
    donor: np.ndarray
    name: str = "B?"
    extra_obs: pd.DataFrame | None = None
    plan: dict = field(default_factory=dict)
    plan_cols: list = field(default_factory=list)

    @property
    def donors(self) -> list:
        vals, counts = np.unique(self.donor, return_counts=True)
        return [str(v) for v, c in zip(vals, counts) if c >= CELLS_PER_DONOR]


def load_background(name: str, spec: dict) -> Background:
    """An .h5ad (raw counts in X or layers['counts']) or an .npz written by `save_npz`."""
    path = Path(spec["path"])
    if path.suffix == ".npz":
        z = np.load(path, allow_pickle=False)
        X, genes, obs = z["X"], list(z["genes"]), pd.DataFrame({k[4:]: z[k] for k in z.files if k.startswith("obs_")})
    else:
        import anndata  # lazily: only the panel run needs it
        ad = anndata.read_h5ad(path)
        X = ad.layers["counts"] if "counts" in ad.layers else ad.X
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        genes, obs = [str(g) for g in ad.var_names], ad.obs.reset_index(drop=True)
    for col, val in (spec.get("filter") or {}).items():
        keep = np.asarray(obs[col]).astype(str) == str(val)
        X, obs = X[keep], obs[keep].reset_index(drop=True)
    X = np.asarray(X, dtype=np.float64)
    if not (np.all(X >= 0) and np.allclose(X, np.round(X))):
        raise ValueError(f"{name}: the panel needs raw counts")
    return Background(X, genes, np.asarray(obs[spec["donor"]]).astype(str), name, obs)



def load_backgrounds(spec_path) -> dict:
    """Every background of a backgrounds JSON ({"B1": {"path", "donor", "filter"}, ...}, written by
    the rule of v1.md section 3.1), loaded and planned (`plan_background`)."""
    spec = json.loads(Path(spec_path).read_text())
    bgs = {k: load_background(k, v) for k, v in spec.items()}
    for bg in bgs.values():
        plan_background(bg)
    return bgs

def _norm_log(X: np.ndarray) -> np.ndarray:
    tot = X.sum(axis=1, keepdims=True)
    return np.log1p(X / np.where(tot > 0, tot, 1.0) * 1e4)


def _corr_cols(Y: np.ndarray) -> np.ndarray:
    Y = Y - Y.mean(axis=0)
    sd = np.sqrt((Y ** 2).sum(axis=0))
    C = (Y.T @ Y) / np.outer(np.where(sd > 0, sd, 1.0), np.where(sd > 0, sd, 1.0))
    return C


def plan_background(bg: Background, seed: int = 20261008) -> dict:
    """Fixed before the key: kept genes, analysed pair, positive and negative control pairs and
    the random genes of N6c. Rule: among the 300 most expressed genes other than mitochondrial
    and ribosomal ones, the analysed pair is the most correlated pair of log-normalized
    expression; the positive control the most correlated pair disjoint from it; the negative
    control the pair closest in mean expression to the analysed pair whose correlation is below
    0.02 in absolute value."""
    rng = np.random.default_rng(seed)
    means = bg.X.mean(axis=0)
    names = [str(g) for g in bg.genes]
    upper = [g.upper() for g in names]
    skip = np.array([u.startswith(("MT-", "RPS", "RPL")) for u in upper])
    order = np.argsort(-means, kind="mergesort")
    cand = [j for j in order if not skip[j]][:300]
    sub = bg.X if bg.X.shape[0] <= 5000 else bg.X[rng.choice(bg.X.shape[0], 5000, replace=False)]
    C = _corr_cols(_norm_log(sub)[:, cand])
    np.fill_diagonal(C, -np.inf)
    flat = np.argsort(-C, axis=None, kind="mergesort")
    pairs = [(cand[i // len(cand)], cand[i % len(cand)]) for i in flat[: 4 * len(cand)]]
    pairs = [(a, b) for a, b in pairs if a < b]
    a1, b1 = pairs[0]
    pos = next((a, b) for a, b in pairs if len({a, b, a1, b1}) == 4)
    np.fill_diagonal(C, 1.0)
    target = (means[a1] + means[b1]) / 2
    neg, best = None, np.inf
    at = {g: k for k, g in enumerate(cand)}
    for i, a in enumerate(cand):
        for b in cand[i + 1:]:
            if len({a, b, a1, b1, *pos}) < 6 or abs(C[at[a], at[b]]) >= 0.02:
                continue
            gap = abs((means[a] + means[b]) / 2 - target)
            if gap < best:
                neg, best = (a, b), gap
    if neg is None:
        raise ValueError(f"{bg.name}: no expression-matched pair with |r| < 0.02 for the negative control")
    g2m = [g for g in G2M_GENES if g in upper]
    named = {a1, b1, *pos, *neg} | {upper.index(g) for g in g2m}
    keep = list(dict.fromkeys([*named, *[j for j in order if j not in named]]))[:max(N_GENES, len(named))]
    pool = [j for j in keep if j not in named]
    random_genes = sorted(rng.choice(pool, size=min(RANDOM_SCORE_GENES, len(pool)), replace=False).tolist())
    plan = dict(genes=[names[j] for j in sorted(keep)],
                pair=[names[a1], names[b1]], pos_pair=[names[pos[0]], names[pos[1]]],
                neg_pair=[names[neg[0]], names[neg[1]]],
                g2m=[names[upper.index(g)] for g in g2m], random_genes=[names[j] for j in random_genes])
    bg.plan = plan
    bg.plan_cols = [names.index(g) for g in plan["genes"]]
    return plan


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
    return bg.X[np.ix_(rows, bg.plan_cols)].copy(), np.asarray(donor_of)


def _two_groups(bg: Background, rng, per_group: int):
    donors = rng.choice(bg.donors, size=2 * per_group, replace=False)
    X, donor = _draw_cells(bg, donors, rng)
    group_of = {d: ("A" if k < per_group else "B") for k, d in enumerate(rng.permutation(donors))}
    return X, donor, np.array([group_of[d] for d in donor])


def build(entry: dict, bgs: dict, pilot: dict) -> tuple[np.ndarray, pd.DataFrame, list, list[dict]]:
    """Counts, obs, genes and claim card(s) of one dataset."""
    cond = conditions()[entry["condition"]]
    bg = bgs[cond.background]
    rng = np.random.default_rng(entry["seed"])
    genes = list(bg.plan["genes"])
    ia, ib = genes.index(bg.plan["pair"][0]), genes.index(bg.plan["pair"][1])
    side, other = entry["side"], ("B" if entry["side"] == "A" else "A")
    variant = entry["variant"]
    replicate_col = "donor"
    if cond.name == "N5":
        donors = rng.choice(bg.donors, size=DONORS_PER_GROUP, replace=False)
        X, donor = _draw_cells(bg, donors, rng)
        group = rng.choice(np.array(["A", "B"]), size=len(donor))
    elif cond.name == "N4":
        X, donor, group = _two_groups(bg, rng, 3)
    elif cond.name == "N7":
        per = len(bg.donors) // 2
        X, donor, group = _two_groups(bg, rng, per)
    else:
        X, donor, group = _two_groups(bg, rng, DONORS_PER_GROUP)
    in_side = group == side
    if cond.name == "N2":
        c = float(variant.split("=")[1])
        X[in_side] = thin(X[in_side], c, rng)
    elif cond.name == "N3":
        X[in_side] = drop_out(X[in_side], float(variant.split("=")[1]), rng)
    elif cond.name.startswith("E"):
        factor = 1.0 if cond.name != "E1" or variant == "dose=key" else float(variant.split("=")[1])
        dose = factor * float(pilot["key_dose"])
        X[in_side] = inject_coupling(X[in_side], ia, ib, dose, rng)
        X[~in_side] = sham_coupling(X[~in_side], ia, ib, dose, rng)
        if cond.name == "E2":   # artifact against the effect: capture loss where the signal is
            X[in_side] = thin(X[in_side], 0.5, rng)
        elif cond.name == "E3":  # artifact with the effect: capture loss on the other side
            X[~in_side] = thin(X[~in_side], 0.5, rng)
    obs = pd.DataFrame({"group": group, "donor": donor})
    obs["total_counts"] = X.sum(axis=1)
    obs["n_genes_by_counts"] = (X > 0).sum(axis=1)
    # the claim's direction: the signal side is the higher one (E); a coin for the nulls, so the
    # card does not reveal the condition
    higher = side if cond.name.startswith("E") else ("A", "B")[int(rng.integers(2))]
    direction = "decrease" if higher == "A" else "increase"  # change from A to B
    metric = {"N6a": "random", "N6b": "constant", "N6c": "score"}.get(cond.name, "norm_pearson")
    base = dict(id=entry["id"], background=cond.background, metric=metric,
                gene_pair=bg.plan["pair"], pos_pair=bg.plan["pos_pair"], neg_pair=bg.plan["neg_pair"],
                score_genes=bg.plan["random_genes"] if metric == "score" else None,
                group_col="group", groups=["A", "B"], replicate_col=replicate_col,
                signal_test=(dict(kind="module", genes=bg.plan["g2m"]) if metric == "score"
                             else dict(kind="coupling", genes=bg.plan["pair"])),
                prereg=dict(estimand="composition", direction=direction, sesoi=float(pilot["sesoi"]),
                            judgment_pending=False))
    cards = [base]
    if cond.cards == 2:  # N4: the same cells analysed without and with the replicate unit
        cards = [dict(base, id=f"{entry['id']}a", replicate_col=None),
                 dict(base, id=f"{entry['id']}b")]
    return X, obs, genes, cards


# --------------------------------------------------------------------------- #
# files
# --------------------------------------------------------------------------- #
def save_npz(path: Path, X: np.ndarray, obs: pd.DataFrame, genes: list) -> None:
    """Plain arrays only (strings as fixed-width unicode), loadable without pickle."""
    arrays = {}
    for c in obs.columns:
        v = np.asarray(obs[c])
        arrays[f"obs_{c}"] = v.astype(str) if v.dtype == object else v
    np.savez_compressed(path, X=np.asarray(X, dtype=np.int32), genes=np.asarray(genes), **arrays)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_panel(key_seed: int, bgs: dict, pilot: dict, out: Path, only: list | None = None) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    entries = assign(key_seed)
    if only is not None:
        entries = [e for e in entries if e["id"] in set(only)]
    manifest = dict(n_datasets=n_datasets(), pilot=pilot, backgrounds={k: bg.plan for k, bg in bgs.items()},
                    datasets=[])
    for e in entries:
        X, obs, genes, cards = build(e, bgs, pilot)
        data_path = out / f"{e['id']}.npz"
        save_npz(data_path, X, obs, genes)
        for card in cards:
            card_path = out / f"{card['id']}.json"
            card_path.write_text(json.dumps(dict(card, data=data_path.name), indent=1))
            manifest["datasets"].append(dict(id=card["id"], data=data_path.name, data_sha256=sha256(data_path),
                                             card=card_path.name, card_sha256=sha256(card_path)))
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1))
    return manifest


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--key-seed", type=int, required=True)
    p.add_argument("--backgrounds", required=True, help="JSON: {B1: {path, donor, filter}, B2: {...}}")
    p.add_argument("--pilot", required=True, help="pilot.json from oracle.py")
    p.add_argument("--out", required=True)
    args = p.parse_args(argv)
    bgs = load_backgrounds(args.backgrounds)
    pilot = json.loads(Path(args.pilot).read_text())
    m = write_panel(args.key_seed, bgs, pilot, Path(args.out))
    print(f"{len(m['datasets'])} claim cards from {m['n_datasets']} datasets written to {args.out}")


if __name__ == "__main__":
    main()
