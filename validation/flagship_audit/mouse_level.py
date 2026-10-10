"""Section E of the flagship audit without the composite donor ids: the per-mouse exact permutation
test within sex, young (3m) against old (20m), recomputed from the audit's published per-mouse
table.

Two donor ids name two mice each (3_10_M/3_11_M, 3_38_F/3_39_F), and both mice are also counted on
their own. The audit counted these ids as units (REPORT.md); here they are dropped, as the owner
decided on 2026-10-10. The test is the audit's own (`audit_tms.exact_permutation_p`), on its
per-mouse values (`E.mice` of audit_results.json, branch results/flagship-audit).

    python validation/flagship_audit/mouse_level.py --results audit_results.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_tms import exact_permutation_p  # noqa: E402  (the audit's test, section E)

GROUPS = ("3m", "20m")  # young, old: the audit's groups


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--results", type=Path, required=True, help="the audit's audit_results.json")
    args = ap.parse_args(argv)
    raw = args.results.read_bytes()
    d = json.loads(raw)
    prov = d["meta"]["provenance"]
    mice = d["E"]["mice"]
    composite = sorted(m["donor_id"] for m in mice if "/" in m["donor_id"])
    kept = [m for m in mice if "/" not in m["donor_id"]]
    out = ["# Flagship audit, section E without the composite donor ids",
           f"# audit_results.json sha256 {hashlib.sha256(raw).hexdigest()}; the audit's data_sha256 "
           f"{prov['data_sha256']}, script_sha256 {prov['script_sha256']}",
           f"dropped: {', '.join(composite)} (each names two mice that are also counted on their own)", ""]
    for sex in sorted({m["sex"] for m in kept}):
        young = [m for m in kept if m["sex"] == sex and m["age"] == GROUPS[0]]
        old = [m for m in kept if m["sex"] == sex and m["age"] == GROUPS[1]]
        out.append(f"{sex}: young {', '.join(m['donor_id'] for m in young)}; old {', '.join(m['donor_id'] for m in old)}")
        for col in ("median_nnz", "mi_3bin"):
            a, b = [float(m[col]) for m in young], [float(m[col]) for m in old]
            p, pmin = exact_permutation_p(a, b)
            splits = math.comb(len(a) + len(b), len(a))
            separated = max(b) < min(a) or max(a) < min(b)
            out.append(f"  {col}: {len(a)} young against {len(b)} old; means {sum(a) / len(a):.4g} against "
                       f"{sum(b) / len(b):.4g}; exact two-sided p = {p:.4f} ({round(p * splits)}/{splits}), "
                       f"smallest attainable {pmin:.4f}; complete separation: {'yes' if separated else 'no'}")
    print("\n".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
