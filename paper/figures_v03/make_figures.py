"""The five figures of paper/manuscript_v03.md, each drawn from a published file of the repository.

    python paper/figures_v03/make_figures.py --scores scores.json --audit audit_results.json \
        --b2-log validation/exploratory/b2_sex_structure.log \
        --p16-log validation/probes/p16_variable_capture_n8.log --out paper/figures_v03

* scores.json: branch results/panel-v1 (commit adc7d65), the scores of validation v1.
* audit_results.json: branch results/flagship-audit (commit d2a5420), the flagship audit on TMS.
* b2_sex_structure.log: B2 by sex and GATE 4's probe for R1's Xist claim (workflow run 38045140802).
* p16_variable_capture_n8.log: the development probe p16 on simulated data.

Fig. 1 the design of v1 (claim cards per condition); Fig. 2 the criteria against their limits;
Fig. 3 false SUPPORTED by condition at the high expression level; Fig. 4 TMS per mouse; Fig. 5 the
real positive controls of v0.1.1 and v0.3. Every number drawn is also in the manuscript's text or
tables. Needs matplotlib (an optional dependency of the package).
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "validation" / "flagship_audit"))
from metric_autopsy.stats import clopper_pearson  # noqa: E402
from audit_tms import exact_permutation_p  # noqa: E402  (the audit's own test, section E)

# the reference palette (dataviz skill, light mode), validated for these three slots all-pairs
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 9,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "grid.linestyle": "-",
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "legend.labelcolor": INK2,
})


def _limit_line(ax, x, label):
    """A solid reference line with its label just above the x-axis, clear of the title."""
    ax.axvline(x, color=INK2, linewidth=1.0, zorder=1)
    ax.text(x, 0.015, f" {label}", transform=ax.get_xaxis_transform(), color=INK2, fontsize=7.5,
            va="bottom", ha="left", zorder=4, bbox=dict(boxstyle="square,pad=0.1", fc=SURFACE, ec="none"))


def fig1_design(scores: dict, out: Path):
    """Claim cards per condition, grouped by what the data are."""
    pc = scores["per_condition"]
    order = [("Nulls and artifacts", ["N1", "N2", "N3", "N4", "N5", "N7", "N8"]),
             ("Useless metrics", ["N6a", "N6b", "N6c"]),
             ("Real effects", ["E1", "E2", "E3"])]
    rows = []
    for group, conds in order:
        for c in conds:
            n = sum(int(v["n"]) for k, v in pc.items() if k.split(":")[0] == c)
            rows.append((group, c, n))
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    y = list(range(len(rows)))[::-1]
    ax.barh(y, [r[2] for r in rows], height=0.55, color=BLUE, zorder=2)
    for yi, (_, c, n) in zip(y, rows):
        ax.text(n + 20, yi, f"{n:,}", va="center", color=INK2, fontsize=8)
    ax.set_yticks(y, [f"{r[1]}" for r in rows])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("claim cards")
    total = sum(r[2] for r in rows)
    ax.set_title(f"Validation v1: {total:,} claim cards in 9,750 datasets")
    # group labels at the left margin
    start = 0
    for group, conds in order:
        ys = y[start:start + len(conds)]
        ax.text(-0.17, (max(ys) + min(ys)) / 2, group, transform=ax.get_yaxis_transform(), ha="right",
                va="center", color=INK, fontsize=8.5, fontweight="bold")
        start += len(conds)
    ax.set_xlim(0, max(r[2] for r in rows) * 1.12)
    fig.subplots_adjust(left=0.30, right=0.97, top=0.90, bottom=0.13)
    fig.savefig(out / "fig1_design.png", dpi=200)
    plt.close(fig)


def fig2_criteria(scores: dict, out: Path):
    """Every cell of the error criteria as a share of its allowance, and S3 as a share of its minimum."""
    cr = scores["criteria"]
    err = [(f"S1 {c}", v["k"], v["max_allowed"]) for c, v in cr["S1"]["conditions"].items()]
    for s in ("S2", "S4", "S5", "S6"):
        for cell, v in cr[s]["cells"].items():
            err.append((cell.replace(":", " "), v["k"], v["max_allowed"]))
    err.append(("S7b crashes", cr["S7b"]["k"], cr["S7b"]["max_allowed"]))
    s3 = [(f"S3 {st}", v["k"], v["min_required"]) for st, v in cr["S3"]["strata"].items()]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.4, 5.6), gridspec_kw=dict(width_ratios=[3, 2]))
    y = list(range(len(err)))[::-1]
    a.scatter([k / m for _, k, m in err], y, s=36, color=BLUE, edgecolor=SURFACE, linewidth=1.5, zorder=3)
    a.set_yticks(y, [f"{lab}: {k} (≤ {m})" for lab, k, m in err], fontsize=7.5)
    a.set_xlim(0, 1.08)
    a.grid(axis="y", visible=False)
    _limit_line(a, 1.0, "limit")
    a.set_xlabel("errors as a share of the allowance")
    a.set_title("(a) error criteria: every cell below its limit")
    y3 = list(range(len(s3)))[::-1]
    b.scatter([k / m for _, k, m in s3], y3, s=36, color=BLUE, edgecolor=SURFACE, linewidth=1.5, zorder=3)
    b.set_yticks(y3, [f"{lab}: {k} (≥ {m})" for lab, k, m in s3], fontsize=7.5)
    b.set_xlim(0.9, 1.5)
    b.set_ylim(-0.8, len(s3) - 0.2)
    b.grid(axis="y", visible=False)
    _limit_line(b, 1.0, "minimum")
    b.set_xlabel("correct definite outcomes\nas a share of the minimum")
    b.set_title("(b) S3: decisiveness")
    fig.text(0.01, 0.01, "S7a (rule violations): 0 of 10,540, where none is allowed.", color=INK2, fontsize=8)
    fig.subplots_adjust(left=0.27, right=0.97, top=0.93, bottom=0.12, wspace=0.95)
    fig.savefig(out / "fig2_criteria.png", dpi=200)
    plt.close(fig)


def fig3_by_level(scores: dict, p16_log: Path, out: Path):
    """False SUPPORTED at the high expression level, with 95% Clopper-Pearson intervals."""
    pl = scores["per_level"]

    def sup(key):
        e = pl[key]
        return int(e["outcomes"]["SUPPORTED"]["k"]), int(e["n"])
    rows = [("N1 (null)", *sup("N1:null:high")), ("N5 (sham)", *sup("N5:sham:high")),
            ("N7 (null over donor ids)", *sup("N7:mice:high"))]
    k = sum(r[1] for r in rows)
    n = sum(r[2] for r in rows)
    rows.append(("pure nulls together", k, n))
    rows += [("N3, dropout f = 0.4", *sup("N3:f=0.4:high")), ("N8, per-cell capture", *sup("N8:beta(2,2):high"))]
    m = re.search(r"N8 beta\(2,2\).*?\n\s*false SUPPORTED by the engine (\d+)/(\d+)", p16_log.read_text())
    dev = (int(m.group(1)), int(m.group(2)))
    fig, ax = plt.subplots(figsize=(7.0, 3.9))
    y = list(range(len(rows) + 1))[::-1]
    for yi, (lab, k, n) in zip(y[1:], rows):
        lo, hi = clopper_pearson(k, n)
        ax.plot([100 * lo, 100 * hi], [yi, yi], color=BLUE, linewidth=2, solid_capstyle="round", zorder=2)
        ax.scatter([100 * k / n], [yi], s=48, color=BLUE, edgecolor=SURFACE, linewidth=1.5, zorder=3)
        ax.text(100 * hi + 0.4, yi, f"{k}/{n}", va="center", color=INK2, fontsize=8)
    lo, hi = clopper_pearson(*dev)
    ax.plot([100 * lo, 100 * hi], [y[0], y[0]], color=ORANGE, linewidth=2, solid_capstyle="round", zorder=2)
    ax.scatter([100 * dev[0] / dev[1]], [y[0]], s=48, color=ORANGE, edgecolor=SURFACE, linewidth=1.5, zorder=3)
    ax.text(100 * hi + 0.4, y[0], f"{dev[0]}/{dev[1]}", va="center", color=INK2, fontsize=8)
    ax.set_yticks(y, ["N8, development probe p16\n(simulated data, all levels)"] + [r[0] for r in rows])
    ax.grid(axis="y", visible=False)
    _limit_line(ax, 2.5, "nominal 2.5%")
    ax.set_xlim(0, 20)
    ax.set_xlabel("false SUPPORTED, % of cards (95% CI)")
    ax.set_title("False SUPPORTED at high expression")
    ax.set_ylim(-1.1, len(rows) + 0.5)
    fig.text(0.33, 0.015, "At medium and low expression no card was SUPPORTED (0 of 6,988).", color=INK2, fontsize=8)
    ax.scatter([], [], s=48, color=BLUE, label="validation v1, high level (real backgrounds)")
    ax.scatter([], [], s=48, color=ORANGE, label="development probe (simulated)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.18), ncol=2, fontsize=8)
    fig.subplots_adjust(left=0.33, right=0.97, top=0.90, bottom=0.27)
    fig.savefig(out / "fig3_by_level.png", dpi=200)
    plt.close(fig)


def fig4_tms(audit: dict, out: Path):
    """TMS: per-mouse median genes detected against mi_3bin, young and old, by sex. The Census's stages are "3m"
    and "20m" ("20-month-old stage and over"); the legend gives each group's age from its mouse ids (24_60_M is 24
    months)."""
    mice = [m for m in audit["E"]["mice"] if "/" not in m["donor_id"] and m["age"] in ("3m", "20m")]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)
    for ax, sex in zip(axes, ("male", "female")):
        sub = [m for m in mice if m["sex"] == sex]
        young = [m for m in sub if m["age"] == "3m"]
        old = [m for m in sub if m["age"] == "20m"]
        for grp, col in ((young, BLUE), (old, ORANGE)):
            ages = "/".join(str(a) for a in sorted({int(m["donor_id"].split("_")[0]) for m in grp}))
            ax.scatter([m["median_nnz"] for m in grp], [m["mi_3bin"] for m in grp], s=56, color=col,
                       edgecolor=SURFACE, linewidth=1.5, zorder=3, label=f"{ages} months ({len(grp)} mice)")
        p_nnz, _ = exact_permutation_p([m["median_nnz"] for m in young], [m["median_nnz"] for m in old])
        p_mi, _ = exact_permutation_p([m["mi_3bin"] for m in young], [m["mi_3bin"] for m in old])
        ax.set_title(f"{sex}s: exact p = {p_nnz:.3f} (genes), {p_mi:.3f} (MI)", fontsize=9)
        ax.set_xlabel("median genes detected per cell")
        ax.legend(loc="upper left", fontsize=7.5)
    axes[0].set_ylabel("mi_3bin (Smad3–Col1a1)")
    fig.subplots_adjust(left=0.10, right=0.98, top=0.88, bottom=0.16, wspace=0.08)
    fig.savefig(out / "fig4_tms_mice.png", dpi=200)
    plt.close(fig)


def _b2(log: Path) -> tuple[dict, list]:
    """GATE 4's three rows and the per-id rows (stage, class) from b2_sex_structure.log."""
    text = log.read_text()
    g4 = {}
    for lab, key in (("all cells", "all"), ("cells with Xist > 0", "xist"), ("female cells (R1's sham claim)", "female")):
        m = re.search(re.escape(lab) + r"\s+(\d+)\s+([\d.]+)%\s+([+-][\d.]+) \(([+-][\d.]+) to ([+-][\d.]+)\)\s+(\w+)", text)
        g4[key] = dict(n=int(m.group(1)), share=float(m.group(2)), mean=float(m.group(3)), lo=float(m.group(4)),
                       hi=float(m.group(5)), status=m.group(6))
    original = text.split("## The compact B2")[0]
    ids = re.findall(r"^\s{2}(\S+)\s+(\S+)\s+\d+\s+[\d.]+%\s+[\d.]+%\s+\d+/\d+\s+(.+?stage(?: and over)?)\s+"
                     r"(one male|one female|both sexes)$", original, flags=re.M)
    return g4, ids


def fig5_controls(audit: dict, b2_log: Path, out: Path):
    """v0.1.1 on TMS (GATE 1) and v0.3 on B2 (GATE 4 and the design), side by side."""
    g1 = {r["stratum"]["age"]: r["n_genes_ratio"] for r in audit["F"]["Xist"]["gate1_table"]}
    kept = (audit["F"]["Xist"]["gate2"]["retained_frac"], audit["F"]["Y_score"]["gate2"]["retained_frac"])
    g4, ids = _b2(b2_log)
    fig = plt.figure(figsize=(7.6, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.25], hspace=0.7, wspace=1.05)
    a = fig.add_subplot(gs[0, 0])
    strata = ["3m", "18m", "20m"]
    ya = list(range(len(strata)))[::-1]
    a.scatter([g1[s] for s in strata], ya, s=56, color=BLUE, edgecolor=SURFACE, linewidth=1.5, zorder=3)
    for yi, s in zip(ya, strata):
        a.text(g1[s], yi + 0.22, f"{g1[s]:.2f}×", ha="center", va="bottom", color=INK2, fontsize=8)
    a.set_yticks(ya, [f"{s.replace('m', ' months')}" for s in strata])
    a.set_ylim(-0.9, len(strata) - 0.5)
    a.set_xlim(0.9, 2.5)
    a.grid(axis="y", visible=False)
    _limit_line(a, 1.5, "GATE 1 limit 1.5×")
    a.set_xlabel("genes detected, female ÷ male (median)")
    a.set_title("(a) v0.1.1 on TMS: GATE 1 kills\nXist and the Y score", fontsize=9)
    a.text(0.0, -0.42, f"GATE 2 keeps {100 * kept[0]:.0f}% (Xist) and {100 * kept[1]:.0f}% (Y score)",
           transform=a.transAxes, color=INK2, fontsize=8)
    b = fig.add_subplot(gs[0, 1])
    rows = [("all cells (R1)", g4["all"]), ("cells with Xist", g4["xist"]), ("female cells (sham)", g4["female"])]
    yb = list(range(len(rows)))[::-1]
    for yi, (lab, r) in zip(yb, rows):
        b.plot([r["lo"], r["hi"]], [yi, yi], color=BLUE, linewidth=2, zorder=2)
        b.scatter([r["mean"]], [yi], s=56, color=BLUE, edgecolor=SURFACE, linewidth=1.5, zorder=3)
    b.set_yticks(yb, [f"{lab}: {r['mean']:+.3f}, {r['status']}\n({r['share']:.1f}% with Xist)" for lab, r in rows],
                 fontsize=7.5)
    b.set_xlim(0, 0.8)
    b.set_ylim(-0.9, len(rows) - 0.5)
    b.grid(axis="y", visible=False)
    _limit_line(b, 0.25, "$\\delta_{min}$ 0.25")
    b.set_xlabel("GATE 4 response to 2-fold Xist (95% CI)")
    b.set_title("(b) v0.3 on B2: GATE 4 over all\ncells calls Xist invalid", fontsize=9)
    c = fig.add_subplot(gs[1, :])
    stage_order = ["2-week-old stage", "5-week-old stage", "8-week-old stage", "2-month-old stage",
                   "3-month-old stage", "4-month-old stage", "6-month-old stage", "20-month-old stage and over"]
    classes = (("one male", BLUE, "single male"), ("one female", ORANGE, "single female"),
               ("both sexes", AQUA, "pool of both sexes"))
    counts = {s: {cl: 0 for cl, _, _ in classes} for s in stage_order}
    for _, _, stage, cl in ids:
        counts[stage.strip()][cl] += 1
    yc = list(range(len(stage_order)))[::-1]
    left = [0] * len(stage_order)
    for cl, col, lab in classes:
        vals = [counts[s][cl] for s in stage_order]
        c.barh(yc, vals, left=left, height=0.55, color=col, edgecolor=SURFACE, linewidth=2, zorder=2, label=lab)
        for yi, v, l0 in zip(yc, vals, left):
            if v:
                c.text(l0 + v + 0.2, yi, f"{v}", va="center", color=INK2, fontsize=8)
        left = [l0 + v for l0, v in zip(left, vals)]
    c.set_yticks(yc, [s.replace("-old stage", "").replace(" and over", "+").replace("-", " ") for s in stage_order])
    c.set_xlim(0, max(left) + 2)
    c.grid(axis="y", visible=False)
    c.set_xlabel("B2 donor ids (samples)")
    c.set_title("(c) B2 by development stage: single males and single females never share a stage", fontsize=9)
    c.legend(loc="lower right", fontsize=8)
    fig.subplots_adjust(left=0.13, right=0.97, top=0.92, bottom=0.08)
    fig.savefig(out / "fig5_positive_controls.png", dpi=200)
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--scores", type=Path, required=True)
    ap.add_argument("--audit", type=Path, required=True)
    ap.add_argument("--b2-log", type=Path, required=True)
    ap.add_argument("--p16-log", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    args = ap.parse_args(argv)
    scores = json.loads(args.scores.read_text())
    audit = json.loads(args.audit.read_text())
    args.out.mkdir(parents=True, exist_ok=True)
    fig1_design(scores, args.out)
    fig2_criteria(scores, args.out)
    fig3_by_level(scores, args.p16_log, args.out)
    fig4_tms(audit, args.out)
    fig5_controls(audit, args.b2_log, args.out)
    print("\n".join(sorted(str(p.relative_to(args.out)) for p in args.out.glob("fig*.png"))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
