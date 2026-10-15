"""Background B2 by sex: which donor ids hold one mouse of one sex, where both sexes share a
development stage, and how GATE 4's probe for anchor R1's Xist claim responds over all cells and
over the cells with Xist. Exploratory, after R1's result (DEVIATIONS.md D1): nothing here re-runs
R1 or changes a verdict.

R1 (validation/prereg/anchors.py) compares female with male mice within development stage, with
the donor id as the mouse. Two inputs:

* the original B2 (the backgrounds-data artifact of v1's workflow run 37990113304, checked by the
  run tag's backgrounds.json), loaded exactly as R1 loads it (anchors._load: Xist and the four Y
  genes, then the most expressed genes, 2,000 in all), with each cell's annotated sex and
  development stage;
* the compact B2 (compact/B2.npz of results/panel-v1: the panel's 2,000 genes, with Xist but
  without the Y genes, and no sex or stage), which classifies the ids by Xist alone.

An id's class is set by the share of its cells with Xist > 0: one male (at most 5%), one female
(at least 95%), or both sexes (in between). The gaps between the classes are printed, so another
threshold can be checked against them. GATE 4's probe is the engine's own
(gates.gate4_signal_response) with R1's settings, as run_autopsy calls it: module(Xist), fold 2,
frac 1.0, against its sham, 200 injections, seed 0, delta_min 0.25 (0.5 x SESOI), the declared
direction an increase. It runs on all cells (R1's Xist claim: FAIL), on the cells with Xist > 0
and, with the original B2, on the cells annotated female: the data of R1's sham claim (female
against female), whose GATE 4 passed (metric validity PASS in anchors_r1.json), so this line
repeats it. Next to each the response without the random thinning: the metric's difference at
the expected counts (the injection keeps Xist and halves the other genes, the sham halves every
gene). With the original B2, GATE 4 runs on it only (R1's data).

    python validation/exploratory/b2_sex_structure.py --backgrounds backgrounds.json --data-dir DATA \
        [--compact compact/B2.npz]
    python validation/exploratory/b2_sex_structure.py --compact compact/B2.npz
"""
from __future__ import annotations

import argparse
import hashlib
import json
from functools import partial
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation" / "prereg"))
import anchors  # noqa: E402  (sets the numerical environment before numpy loads, as anchors.main)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from metric_autopsy import SimpleData, __version__, gates, injected_signal  # noqa: E402
from metric_autopsy.effect import _replicate_design  # noqa: E402  (run_autopsy's design rule)

PREFIX = "mouse_pancreatic_islet_atlas_Hrovatin__"  # every B2 donor id starts so
MALE_MAX, FEMALE_MIN = 0.05, 0.95  # an id's class by the share of its cells with Xist > 0
CLASSES = ("one male", "one female", "both sexes")
DELTA_MIN = 0.5 * anchors.SESOI  # report.delta_min_of: 0.5 x SESOI, R1's delta_min 0.25
ALPHA, SEED = 0.05, 0  # run_autopsy's defaults, as R1 ran


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def src_digest() -> str:
    """sha256 of the list of src/'s files' sha256, as DEVIATIONS.md D1 records it
    (find src -type f | LC_ALL=C sort | xargs sha256sum | sha256sum), byte-code caches left out."""
    files = sorted(p.relative_to(ROOT).as_posix() for p in (ROOT / "src").rglob("*")
                   if p.is_file() and "__pycache__" not in p.parts and not any(s.endswith(".egg-info") for s in p.parts))
    return hashlib.sha256("".join(f"{sha256(ROOT / f)}  {f}\n" for f in files).encode()).hexdigest()


def pct(x: float) -> str:
    return f"{100 * x:.1f}%"


def load_original(spec: Path, data_dir: Path) -> SimpleData:
    """B2 as anchors.r1 builds R1's data: sex, age (the development stage), mouse (the donor id)."""
    X, genes, obs = anchors._load(spec, "B2", data_dir, named=("Xist", *anchors.Y_GENES))
    obs = obs.rename(columns={"donor_id": "mouse", "development_stage": "age"})
    return SimpleData(X, obs[["sex", "age", "mouse"]].astype(str), genes)


def load_compact(path: Path) -> SimpleData:
    z = np.load(path, allow_pickle=False)
    return SimpleData(z["X"], pd.DataFrame({"mouse": z["donor"].astype(str)}), [str(g) for g in z["genes"]])


def per_id(data: SimpleData) -> pd.DataFrame:
    """One row per donor id: its cells, the shares with Xist > 0 and (if kept) with any Y gene > 0,
    the cells of each annotated sex and the stages (original B2), and its class."""
    X, obs = data.X, data.obs
    xist = X[:, anchors._cols(data.var_names, ["Xist"])].sum(axis=1)
    ycols = anchors._cols(data.var_names, anchors.Y_GENES)
    y = X[:, ycols].sum(axis=1) if ycols else None
    rows = []
    for mouse in sorted(set(obs["mouse"])):
        m = (obs["mouse"] == mouse).to_numpy()
        group, _, sample = mouse.removeprefix(PREFIX).partition("__")
        r = dict(id=mouse, group=group, sample=sample, cells=int(m.sum()), xist=float((xist[m] > 0).mean()))
        if y is not None:
            r["y"] = float((y[m] > 0).mean())
        if "sex" in obs:
            r["female_cells"] = int((obs["sex"][m] == "female").sum())
            r["male_cells"] = int((obs["sex"][m] == "male").sum())
        if "age" in obs:
            r["stage"] = "; ".join(sorted(set(obs["age"][m])))
        r["class"] = CLASSES[0] if r["xist"] <= MALE_MAX else CLASSES[1] if r["xist"] >= FEMALE_MIN else CLASSES[2]
        rows.append(r)
    return pd.DataFrame(rows)


def structure(ids: pd.DataFrame, data: SimpleData) -> list[str]:
    out = [f"{data.n_obs} cells, {len(ids)} donor ids, {data.n_vars} genes "
           f"(Y genes kept: {', '.join(g for g in anchors.Y_GENES if anchors._cols(data.var_names, [g])) or 'none'})",
           "", "per donor id (the id's study group and sample, as the id names them):"]
    head = f"  {'group':16s} {'sample':14s} {'cells':>5s} {'Xist>0':>7s}"
    head += f" {'Y>0':>6s}" if "y" in ids else ""
    head += f" {'annotated F/M cells':>20s}" if "female_cells" in ids else ""
    head += f"  {'development stage':30s}" if "stage" in ids else ""
    out.append(head + "  class")
    for _, r in ids.sort_values(["class", "group", "sample"]).iterrows():
        line = f"  {r.group:16s} {r['sample']:14s} {r.cells:5d} {pct(r.xist):>7s}"
        line += f" {pct(r.y):>6s}" if "y" in ids else ""
        line += f" {f'{r.female_cells}/{r.male_cells}':>20s}" if "female_cells" in ids else ""
        line += f"  {r.stage:30s}" if "stage" in ids else ""
        out.append(line + f"  {r['class']}")
    xist = data.X[:, anchors._cols(data.var_names, ["Xist"])].sum(axis=1)
    out += ["", f"cells with Xist > 0: {int((xist > 0).sum())} of {data.n_obs} ({100 * (xist > 0).mean():.2f}%)", "",
            f"classes by the share of an id's cells with Xist > 0 (one male <= {pct(MALE_MAX)}, one female >= "
            f"{pct(FEMALE_MIN)}, both sexes in between):"]
    for c in CLASSES:
        sub = ids[ids["class"] == c]
        if sub.empty:
            out.append(f"  {c}: no id")
            continue
        groups = ", ".join(f"{g} {n}" for g, n in sub["group"].value_counts().sort_index().items())
        line = f"  {c}: {len(sub)} ids, Xist > 0 in {pct(sub.xist.min())}-{pct(sub.xist.max())} of their cells"
        line += f", a Y gene > 0 in {pct(sub.y.min())}-{pct(sub.y.max())}" if "y" in sub else ""
        out.append(line + f"; groups: {groups}")
    shares = {c: ids.loc[ids["class"] == c, "xist"] for c in CLASSES}
    if all(len(s) for s in shares.values()):
        out.append(f"  the same classes for a male bound in [{pct(shares[CLASSES[0]].max())}, {pct(shares[CLASSES[2]].min())}) "
                   f"and a female bound in ({pct(shares[CLASSES[2]].max())}, {pct(shares[CLASSES[1]].min())}]")
    return out


def stage_days(stage: str) -> tuple:
    """A development stage's age in days, for ordering ('2-week-old stage', '20-month-old stage and
    over'); a stage it cannot read sorts last."""
    num, _, rest = stage.partition("-")
    unit = 7.0 if rest.startswith("week") else 30.4 if rest.startswith("month") else None
    if num.isdigit() and unit is not None:
        return (0, float(num) * unit, stage)
    return (1, 0.0, stage)


def stages(ids: pd.DataFrame, data: SimpleData) -> list[str]:
    """Where the sexes meet: by development stage, the ids of each class and of each annotated sex,
    R1's allowed set as anchors.r1_allowed computes it, and run_autopsy's design rule on R1's obs."""
    obs = data.obs
    out = ["development stages (each id lies in one stage):",
           f"  {'stage':30s} {'one male':>8s} {'one female':>10s} {'both sexes':>10s}   ids with cells annotated female / male"]
    for stage in sorted(set(obs["age"]), key=stage_days):
        sub = ids[ids["stage"] == stage]
        n = {c: int((sub["class"] == c).sum()) for c in CLASSES}
        out.append(f"  {stage:30s} {n[CLASSES[0]]:8d} {n[CLASSES[1]]:10d} {n[CLASSES[2]]:10d}   "
                   f"{int((sub.female_cells > 0).sum())} / {int((sub.male_cells > 0).sum())}")
    single = ids[ids["class"] != CLASSES[2]].groupby("stage")["class"].nunique()
    meet = sorted(single[single > 1].index)
    both = ids[(ids.female_cells > 0) & (ids.male_cells > 0)]
    out += [f"  stages where an id of one male and an id of one female meet: {', '.join(meet) or 'none'}",
            "  stages with cells of both annotated sexes: " + ("; ".join(
                f"{s} ({', '.join(f'{r.group}__{r['sample']} {r['class']}' for _, r in g.iterrows())})"
                for s, g in both.groupby("stage")) or "none")]
    counts = obs.groupby(["age", "sex"])["mouse"].nunique()
    per_sex = obs.groupby("sex")["mouse"].nunique()
    out += ["", "R1's allowed set (anchors.r1_allowed): ids per annotated sex over all stages: "
            + ", ".join(f"{s} {int(n)}" for s, n in per_sex.items())
            + f" -> {', '.join(anchors.r1_allowed(obs))}; an id with cells of both annotated sexes counts for each",
            f"  ids per (stage, annotated sex), as R1's note counts them: {int(counts.sum())} over {obs['mouse'].nunique()} ids"]
    design = _replicate_design(obs, "sex", ("female", "male"), "mouse", ["age"])
    out.append(f"run_autopsy's design rule (effect._replicate_design) on R1's obs, sex within stage: "
               f"{design['kind']}" + (f": {design['reason']}" if "reason" in design else ""))
    single_ids = set(ids.loc[ids["class"] != CLASSES[2], "id"])
    keep = obs["mouse"].isin(single_ids).to_numpy()
    d2 = _replicate_design(obs[keep], "sex", ("female", "male"), "mouse", ["age"])
    if "strata" in d2:
        sexes = {}
        for r in d2["reps"]:
            sexes.setdefault(d2["strata"][r], set()).update(d2["membership"][r])
        shared = sorted((s for s, g in sexes.items() if len(g) == 2), key=stage_days)
        out.append(f"  the same rule without the {len(ids) - len(single_ids)} ids of both sexes: {d2['kind']}; "
                   f"stages holding both sexes: {', '.join(shared) or 'none'}, so no stratum gives a contrast "
                   "(stats.stratum_weights counts only strata holding both groups)")
    else:
        out.append(f"  the same rule without the ids of both sexes: {d2['kind']}: {d2.get('reason', '')}")
    return out


def expected_response(data: SimpleData, cols: list) -> float:
    """The probe's response without the random thinning: the metric (log1p CP10k of Xist, averaged
    over cells) at the expected counts of the injection (Xist kept, the other genes halved) minus at
    those of the sham (every gene halved)."""
    X = data.X
    x = X[:, cols].sum(axis=1)
    tot = X.sum(axis=1)
    fold = anchors.MODULE_FOLD
    inj_tot, sham_tot = x + (tot - x) / fold, tot / fold
    inj = x / np.where(inj_tot > 0, inj_tot, 1.0) * 1e4
    sham = (x / fold) / np.where(sham_tot > 0, sham_tot, 1.0) * 1e4
    return float(np.mean(np.log1p(inj)) - np.mean(np.log1p(sham)))


def gate4(data: SimpleData, sets: list[tuple[str, np.ndarray]]) -> list[str]:
    cols = anchors._cols(data.var_names, ["Xist"])
    metric = partial(anchors.log_cp10k_mean, cols=cols)
    inject = injected_signal.module([data.var_names[j] for j in cols], fold=anchors.MODULE_FOLD, frac=anchors.MODULE_FRAC)
    out = [f"GATE 4's probe for R1's Xist claim: {inject.description} against its sham, {gates.GATE4_N_REP} injections, "
           f"seed {SEED}, delta_min {DELTA_MIN:g}, declared direction an increase "
           "(gates.gate4_signal_response, as run_autopsy calls it)",
           f"  {'cells':34s} {'n':>6s} {'Xist>0':>7s} {'response (95% CI)':>34s}  outcome  {'without thinning':>16s}"]
    x = data.X[:, cols].sum(axis=1)
    for label, mask in sets:
        sub = data if mask.all() else data[mask]
        g = gates.gate4_signal_response(metric, sub, inject, direction="increase", delta_min=DELTA_MIN, alpha=ALPHA, seed=SEED)
        d = g.detail
        lo, hi = d["ci"]
        resp = f"{d['mean_response']:+.4f} ({lo:+.4f} to {hi:+.4f})"
        out.append(f"  {label:34s} {int(mask.sum()):6d} {pct(float((x[mask] > 0).mean())):>7s} {resp:>34s}  "
                   f"{g.status.value:7s}  {expected_response(sub, cols):+16.4f}")
        out.append(f"    {g.status.value}: {g.message}")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--backgrounds", type=Path, help="the run tag's backgrounds.json (the original B2)")
    ap.add_argument("--data-dir", type=Path, help="the folder holding B2.npz (the backgrounds-data artifact)")
    ap.add_argument("--compact", type=Path, help="compact/B2.npz of results/panel-v1")
    args = ap.parse_args(argv)
    if not (args.backgrounds or args.compact) or bool(args.backgrounds) != bool(args.data_dir):
        ap.error("give --backgrounds with --data-dir (the original B2), --compact, or both")
    out = ["# B2 by sex: donor ids, development stages, and GATE 4's probe for R1's Xist claim",
           "# Exploratory, after R1's result (DEVIATIONS.md D1): R1 is not re-run, and no verdict changes.",
           f"# metric_autopsy {__version__}; src/ sha256 of its files' sha256 list {src_digest()} "
           "(DEVIATIONS.md D1 records the tag's)",
           f"# script sha256 {sha256(Path(__file__))}; anchors.py sha256 {sha256(Path(anchors.__file__))}"]
    original = None
    if args.backgrounds:
        spec = json.loads(args.backgrounds.read_text())["B2"]
        original = load_original(args.backgrounds, args.data_dir)  # resolve_background checks the sha256
        ids = per_id(original)
        out += [f"# the original B2: {spec['file']} sha256 {spec['sha256']} (the run tag's backgrounds.json; "
                f"{spec['source']['title']}, Census {spec['source']['census_release']}, dataset {spec['source']['dataset_id']})",
                "", "## The original B2, loaded as R1 loads it", ""]
        out += structure(ids, original) + [""] + stages(ids, original)
    if args.compact:
        compact = load_compact(args.compact)
        cids = per_id(compact)
        out += ["", f"## The compact B2 ({args.compact.name}, sha256 {sha256(args.compact)}): Xist only", ""]
        out += structure(cids, compact)
        if original is not None:
            same = ids.set_index("id")[["cells", "xist"]].equals(cids.set_index("id")[["cells", "xist"]])
            out.append(f"  the same ids, cells and Xist > 0 shares as the original B2: {'yes' if same else 'no'}")
    data = original if original is not None else compact
    x = data.X[:, anchors._cols(data.var_names, ["Xist"])].sum(axis=1)
    sets = [("all cells", np.ones(data.n_obs, bool)), ("cells with Xist > 0", x > 0)]
    if "sex" in data.obs:
        sets.append(("female cells (R1's sham claim)", (data.obs["sex"] == "female").to_numpy()))
    out += ["", f"## GATE 4 ({'the original' if original is not None else 'the compact'} B2)", ""] + gate4(data, sets)
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
