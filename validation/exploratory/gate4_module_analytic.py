"""GATE 4's module probe on a mean-expression metric, by arithmetic: why the default fails at a SESOI
of 0.5, why R1's probe failed on B2, and why markers of single counts pass at a low depth. Read
after R1's result; the engine is not changed. gate4_module_dilution.log's simulations confirm the
numbers; this script does not simulate.

The metric is R1's: the mean over cells of f(x) = log1p(x), x the marker's CP10k in the cell.

1. An identity, not a finding. A cell without the marker has x = 0 under the injection and under
   the sham, so f = 0 in both; the response over all cells is therefore exactly the share of cells
   with the marker times the mean response in those cells, in every draw.
2. A bound. The module raises the marker's share 2-fold in a random `frac` of cells by halving
   every other gene there, so x becomes at most 2x, and f(2x) - f(x) = ln((1 + 2x) / (1 + x)) is
   below ln 2 for every x. Where the marker has many counts the sham (every gene halved in another
   random `frac` of cells) leaves f's mean unchanged, so the response is at most
   frac x share x ln 2. With fewer counts the sham's thinning of the marker lowers f's mean too, so
   the response exceeds the formula, and at a low depth it exceeds ln 2 (section 3): the bound holds
   for a well-detected marker, such as Xist in B2.
3. A single count at a low depth. With one count of the marker in a cell of depth d, the sham
   keeps it with probability 1/2 and otherwise sets x to 0: the sham's loss adds to the response
   about as much as the 2-fold change.

    python validation/exploratory/gate4_module_analytic.py [--dilution-log gate4_module_dilution.log]
        [--b2-log b2_sex_structure.log]
"""
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
DELTA_MIN = 0.25  # R1's: 0.5 x SESOI 0.5
DEFAULT_FRAC, R1_FRAC = 0.3, 1.0  # injected_signal.module's default; anchors.py's MODULE_FRAC


def f(x: float) -> float:
    return math.log1p(x)


def cp10k(count: float, depth: float) -> float:
    return 1e4 * count / depth


def simulated(log: Path, depth: str, label: str, share: str = "100%") -> tuple[float, float] | None:
    """(frac 0.3, frac 1.0) mean responses of gate4_module_dilution.log at a depth and marker level."""
    text = log.read_text()
    block = text.split(f"## Depth about {depth} counts per cell")[1].split("##")[0]
    m = re.search(re.escape(label) + r"\s+" + re.escape(share) + r"\s+([+-][\d.]+) .*?\s([+-][\d.]+) \(", block)
    return (float(m.group(1)), float(m.group(2))) if m else None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dilution-log", type=Path, default=HERE / "gate4_module_dilution.log")
    ap.add_argument("--b2-log", type=Path, default=HERE / "b2_sex_structure.log")
    args = ap.parse_args(argv)
    ln2 = math.log(2.0)
    out = ["# GATE 4's module probe on a mean-expression metric, by arithmetic (read after R1's result)",
           f"# metric: mean over cells of log1p(CP10k of the marker); delta_min {DELTA_MIN:g} (0.5 x SESOI 0.5)", ""]

    out.append("## 1. The identity: response over all cells = share of cells with the marker x response in them")
    b2 = args.b2_log.read_text()
    m_share = re.search(r"cells with Xist > 0: (\d+) of (\d+)", b2)
    m_all = re.search(r"\n\s+all cells\s+\d+\s+[\d.]+%\s+([+-][\d.]+)", b2)
    m_pos = re.search(r"cells with Xist > 0\s+\d+\s+[\d.]+%\s+([+-][\d.]+)", b2)
    share = int(m_share.group(1)) / int(m_share.group(2)) if m_share else float("nan")
    if m_share and m_all and m_pos:
        out.append(f"  B2 (b2_sex_structure.log): {share:.4f} x {float(m_pos.group(1)):+.4f} = "
                   f"{share * float(m_pos.group(1)):+.4f}; measured over all cells {float(m_all.group(1)):+.4f} "
                   "(separate runs, so equal up to their Monte Carlo error)")
    out.append("")

    out.append("## 2. The bound for a well-detected marker: frac x share x ln 2")
    out.append(f"  ln 2 = {ln2:.4f}; f(2x) - f(x) = ln((1 + 2x)/(1 + x)) < ln 2 for every x, e.g. x = 6: "
               f"{math.log(13 / 7):.4f}, x = 50: {math.log(101 / 51):.4f}")
    out.append(f"  the default (frac {DEFAULT_FRAC:g}) with the marker in every cell: at most {DEFAULT_FRAC * ln2:.4f}, "
               f"below delta_min {DELTA_MIN:g}: it cannot pass")
    out.append(f"  R1's probe (frac {R1_FRAC:g}) needs a share of at least {DELTA_MIN:g} / ln 2 = {DELTA_MIN / ln2:.4f} of "
               f"cells with the marker; B2 has {share:.4f}: the failure follows from B2's composition and "
               "anchors.py's constants" if m_share else "")
    if m_pos:
        pos = float(m_pos.group(1))
        out.append(f"  with B2's measured response in its Xist cells ({pos:+.4f}, below ln 2: Xist is well detected "
                   f"there), the share needed is {DELTA_MIN:g} / {pos:.4f} = {DELTA_MIN / pos:.4f}")
    for depth, label, x in (("20,000", "high (50)", 50.0), ("20,000", "as B2's Xist (6)", 6.0)):
        sim = simulated(args.dilution_log, depth, label)
        if sim:
            out.append(f"  at depth {depth}, marker at CP10k {x:g}: frac 0.3 x ln((1 + 2x)/(1 + x)) = "
                       f"{DEFAULT_FRAC * math.log((1 + 2 * x) / (1 + x)):.4f}; simulated {sim[0]:+.4f} "
                       f"(frac 1.0: {math.log((1 + 2 * x) / (1 + x)):.4f}; simulated {sim[1]:+.4f})")
    out.append("  where the simulation exceeds the formula, the excess is the sham's thinning of the marker's counts, "
               "which the formula leaves out")
    out.append("")

    out.append("## 3. One count of the marker in a cell of depth 2,000")
    d = 2000.0
    base = f(cp10k(1, d + 1))
    inj = f(cp10k(1, d / 2 + 1))  # every other gene halved, the marker kept
    sham = 0.5 * f(cp10k(1, d / 2 + 1))  # every gene halved: the count kept with probability 1/2
    out += [f"  untouched: CP10k {cp10k(1, d + 1):.3f}, f {base:.4f}",
            f"  injected (other genes halved): CP10k {cp10k(1, d / 2 + 1):.3f}, f {inj:.4f}; the 2-fold change "
            f"adds {inj - base:.4f}",
            f"  sham (every gene halved): f {sham:.4f} on average; the sham's loss of the single count adds "
            f"{base - sham:.4f}"]
    for frac in (R1_FRAC, DEFAULT_FRAC):
        resp = frac * (inj - sham)
        out.append(f"  frac {frac:g}: response {frac:g} x ({inj - base:.4f} + {base - sham:.4f}) = {resp:.4f}, of which "
                   f"the 2-fold change {frac * (inj - base):.4f} and the sham's loss {frac * (base - sham):.4f}")
    sim = simulated(args.dilution_log, "2,000", "low (1)")
    if sim:
        out.append(f"  simulated (gate4_module_dilution.log, depth about 2,000, a marker of one count): frac 0.3 "
                   f"{sim[0]:+.4f}, frac 1.0 {sim[1]:+.4f}")
    out.append(f"  so the default passes there only through the sham: the 2-fold change alone gives "
               f"{DEFAULT_FRAC * (inj - base):.4f}, below delta_min {DELTA_MIN:g}")
    out.append(f"  and the response in the marker's cells (frac 1) is {inj - sham:.4f}, above ln 2 = {ln2:.4f}: "
               "the bound of section 2 does not hold for a marker of single counts")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
