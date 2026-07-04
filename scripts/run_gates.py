#!/usr/bin/env python3
"""Thin wrapper the `metric-autopsy` skill runs from a source checkout.

All logic lives in `metric_autopsy.cli`; this shim just puts `src/` on the path so it works
without `pip install`. See `metric-autopsy --help` (installed) or `python scripts/run_gates.py
--demo`.
"""
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
if _SRC.exists():
    sys.path.insert(0, str(_SRC))

from metric_autopsy.cli import main  # noqa: E402

if __name__ == "__main__":
    main()
