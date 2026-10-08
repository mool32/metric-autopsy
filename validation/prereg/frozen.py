"""The frozen files of validation v1 and the guard every job of the real run runs first
(validation/prereg/v1.md, section 2): at the commit a job runs on, the files of the tag
v0.3.0-prereg must be unchanged; only the data the workflow commits may differ. Standard library
only, so that a job can run it before it installs anything.

    python validation/prereg/frozen.py --frozen-tag v0.3.0-prereg
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
FROZEN_PATHS = ("src", "pyproject.toml", "validation/prereg/v1.md", "validation/prereg/*.py",
                "validation/prereg/requirements-*.txt", ".github/workflows/validation.yml")


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=HERE, capture_output=True, text=True, check=True).stdout


def guard(frozen_tag: str) -> dict:
    """The checked-out commit descends from the frozen tag and its engine, protocol, panel code,
    pinned requirements and workflow are byte-identical to the tag's; only data and documents
    (backgrounds.json, pilot.json, logs) may differ."""
    head = _git("rev-parse", "HEAD").strip()
    tag_commit = _git("rev-list", "-n", "1", frozen_tag).strip()
    if subprocess.run(["git", "merge-base", "--is-ancestor", tag_commit, head], cwd=HERE).returncode != 0:
        raise SystemExit(f"guard: HEAD {head[:12]} does not descend from {frozen_tag}")
    changed = [f for f in _git("diff", "--name-only", tag_commit, head, "--", *FROZEN_PATHS).split() if f]
    if changed:
        raise SystemExit(f"guard: frozen files differ from {frozen_tag}: {changed}")
    return dict(head=head, frozen_tag=frozen_tag, frozen_commit=tag_commit)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--frozen-tag", required=True)
    print(json.dumps(guard(p.parse_args(argv).frozen_tag), indent=1))


if __name__ == "__main__":
    main()
