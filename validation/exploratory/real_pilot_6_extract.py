#!/usr/bin/env python3
"""The fourth run's extraction (validation/exploratory/real_pilot_6_plan.md): the full set of genes for the
cells of cases 1, 2 and 3 and case 1's metadata, read again from the CELLxGENE Census (release 2025-11-08)
by the third run's code, unchanged (real_pilot_5_extract.run_census); the Marsh matrices are not read (Marsh
et al. 2022 is no longer a control for coupling). v1's selection environment; it does not import the engine.

    python real_pilot_6_extract.py --data INPUTS/data --out OUT [--synthetic] [--chunk N]
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import real_pilot_5_extract as E5  # noqa: E402  (the third run's extraction: run_census, its checks)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--data", type=Path, required=True, help="the case files, select.json (run 2), select1.json (run 1)")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--synthetic", action="store_true")
    p.add_argument("--chunk", type=int, default=E5.CHUNK)
    args = p.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    lines = []

    def log(*msg):
        for m in msg:
            print(m, flush=True)
            lines.append(str(m))
    tag = "SYNTHETIC " if args.synthetic else ""
    log(f"# {tag}The fourth run's extraction (validation/exploratory/real_pilot_6_plan.md; the third run's "
        "real_pilot_5_extract.run_census)", "",
        f"python {platform.python_version()}; numpy {E5.np.__version__}; pandas {E5.pd.__version__}; Census release "
        f"{'SYNTHETIC' if args.synthetic else E5.S.RELEASE}", "")
    out = dict(synthetic=bool(args.synthetic), release="SYNTHETIC" if args.synthetic else E5.S.RELEASE,
               python=platform.python_version(), numpy=E5.np.__version__, pandas=E5.pd.__version__, cases={})
    t0, status = time.time(), 0
    try:
        E5.run_census(args, log, out)
    except Exception as exc:  # reported as it is; the analyses do not start
        out.update(error=repr(exc), traceback=traceback.format_exc().splitlines())
        log(f"EXCEPTION {exc!r}", *out["traceback"])
        status = 1
    out["seconds"] = round(time.time() - t0, 1)
    log("", f"{tag}extraction: {out['seconds']} s")
    (args.out / "extract.json").write_text(json.dumps(out, indent=1, default=str))
    (args.out / "extract.log").write_text("\n".join(lines) + "\n")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
