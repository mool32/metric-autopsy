"""The frozen files of validation v1 and the guard every job of the real run runs first
(validation/prereg/v1.md, section 2): at the commit a job runs on, the files of the tag
v0.3.0-prereg must be unchanged; only the data the workflow commits may differ. Standard library
only, so that a job can run it before it installs anything; the workflow runs the tag's own copy
(``git show v0.3.0-prereg:validation/prereg/frozen.py``) on the checkout (``--repo``), so that a
changed copy in the checkout cannot pass itself.

    python validation/prereg/frozen.py --frozen-tag v0.3.0-prereg [--repo .]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
FROZEN_PATHS = ("src", "pyproject.toml", "validation/prereg/v1.md", "validation/prereg/*.py",
                "validation/prereg/requirements-*.txt", ".github/workflows/validation.yml")
# One numerical path on every machine (the workflow sets the same): one BLAS thread, OpenBLAS's
# Haswell kernels, numpy's dispatch up to AVX2 (x86-64-v3). The scripts set these before numpy
# loads, unless the caller set them (`pin_numerics`), so that anyone who re-runs a card gets the
# same bytes (blind.py verify). numpy only warns about features a machine lacks.
NUMERIC_ENV = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
               "OPENBLAS_CORETYPE": "Haswell", "NPY_DISABLE_CPU_FEATURES": "X86_V4 AVX512_ICL AVX512_SPR"}


def pin_numerics() -> dict:
    """Set NUMERIC_ENV where the environment does not already set it (before numpy loads)."""
    for k, v in NUMERIC_ENV.items():
        os.environ.setdefault(k, v)
    return {k: os.environ[k] for k in NUMERIC_ENV}


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True, check=True).stdout


def guard(frozen_tag: str, repo: Path | str | None = None) -> dict:
    """The checked-out commit of `repo` (default: this file's repository) descends from the frozen
    tag and its engine, protocol, panel code, pinned requirements and workflow are byte-identical
    to the tag's; only data and documents (backgrounds.json, pilot.json, logs) may differ."""
    repo = Path(repo) if repo is not None else HERE
    head = _git(repo, "rev-parse", "HEAD").strip()
    tag_commit = _git(repo, "rev-list", "-n", "1", frozen_tag).strip()
    if subprocess.run(["git", "merge-base", "--is-ancestor", tag_commit, head], cwd=repo).returncode != 0:
        raise SystemExit(f"guard: HEAD {head[:12]} does not descend from {frozen_tag}")
    changed = [f for f in _git(repo, "diff", "--name-only", tag_commit, head, "--", *FROZEN_PATHS).split() if f]
    if changed:
        raise SystemExit(f"guard: frozen files differ from {frozen_tag}: {changed}")
    return dict(head=head, frozen_tag=frozen_tag, frozen_commit=tag_commit)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--frozen-tag", required=True)
    p.add_argument("--repo", type=Path, help="the checkout to check (default: this file's repository)")
    args = p.parse_args(argv)
    print(json.dumps(guard(args.frozen_tag, args.repo), indent=1))


if __name__ == "__main__":
    main()
