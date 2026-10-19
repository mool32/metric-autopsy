"""Provenance for verdicts: data and pre-registration hashes, versions, and a run log.

Every number that leaves the engine should be traceable to the exact data, the exact
pre-registration and the exact code. A run log (JSON lines) also records how many times a
claim has been run, so that "re-run until it passes" is visible rather than silent.

Every run with a pre-registration is logged, in every interface (Python API, CLI, MCP):
agents work through the API, which is where "re-run until it passes" is most tempting. The
path is ``run_autopsy(log_path=...)`` (``off`` disables), else ``$METRIC_AUTOPSY_LOG``, else
``metric_autopsy_runs.jsonl`` in the working directory. Runs without a pre-registration (which
cannot reach an effect verdict: no estimand) and the demo are not logged unless a path is given.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import platform
from pathlib import Path
from typing import Iterable

import numpy as np

SCHEMA = "metric-autopsy/report/v1"
DEFAULT_LOG = "metric_autopsy_runs.jsonl"
DEFAULT_CLI_LOG = DEFAULT_LOG  # backward-compatible name


def _canonical(obj) -> str:
    return json.dumps(jsonable(obj), sort_keys=True, separators=(",", ":"), default=str)


def sha256_json(obj) -> str:
    return hashlib.sha256(_canonical(obj).encode()).hexdigest()


_MIRROR = {"increase": "decrease", "decrease": "increase"}


def canonical_claim(metric: str, genes, groups, direction) -> dict:
    """A claim's canonical form: the metric, its genes, the groups' labels sorted and the direction
    stated for that order. "Higher in females" and "lower in males" are one claim, so the mirrored
    wording (the groups and the direction swapped) has the same form (v0.4 change 1, decided
    2026-10-12)."""
    given = [str(g) for g in groups]
    order = sorted(given)
    if order != given and direction in _MIRROR:
        direction = _MIRROR[direction]
    return dict(metric=str(metric), genes=[str(g) for g in (genes or ())], groups=order,
                direction=direction)


def claim_seed(data_sha256: str, claim: dict, seed: int, stream: str) -> int:
    """The seed of one random stream of the engine (a gate's draws), keyed by the data, the claim's
    canonical form and the user's seed: a mirrored wording gets the same draws, and two claims on
    the same data get different ones, so that their Monte Carlo errors are independent."""
    key = sha256_json(dict(data=data_sha256, claim=claim, seed=int(seed), stream=str(stream)))
    return int(key[:16], 16)


def sha256_data(data, obs_columns: Iterable[str] = ()) -> str:
    """Hash of X's values, var_names and the named obs columns.

    The values are hashed whatever their storage (dense or sparse, any numeric dtype): row by row,
    the column index and the float64 value of every nonzero entry, as a canonical CSR. The same data
    held as a SimpleData or as an AnnData with a sparse X therefore get the same hash, which since
    v0.4 also keys the engine's random draws (``claim_seed``).
    """
    h = hashlib.sha256()
    X = data.X
    n_obs, n_vars = (int(v) for v in X.shape)
    h.update(f"matrix({n_obs},{n_vars})".encode())
    csr = X.tocsr() if hasattr(X, "tocsr") else None
    step = max(1, 2_000_000 // max(n_vars, 1))
    for i in range(0, n_obs, step):
        if csr is not None:
            block = csr[i:i + step].copy()
            block.sum_duplicates()
            block.eliminate_zeros()
            block.sort_indices()
            indptr, indices, values = block.indptr, block.indices, block.data
        else:
            dense = np.asarray(X[i:i + step])
            nz = dense != 0
            rows, indices = np.nonzero(nz)  # row-major: by row, then by column
            indptr = np.concatenate([[0], np.cumsum(nz.sum(axis=1))])
            values = dense[rows, indices]
        for arr, dtype in ((indptr, np.int64), (indices, np.int64), (values, np.float64)):
            h.update(np.ascontiguousarray(arr, dtype=dtype).tobytes())
    h.update("\x1f".join(map(str, data.var_names)).encode())
    obs = data.obs
    for col in sorted({c for c in obs_columns if c in obs.columns}):
        h.update(col.encode() + b"\x1e")
        h.update("\x1f".join(map(str, np.asarray(obs[col]))).encode())
    return h.hexdigest()


def environment() -> dict:
    from . import __version__
    env = dict(metric_autopsy=__version__, python=platform.python_version(),
               platform=platform.platform(), numpy=np.__version__)
    from importlib import metadata
    for mod in ("pandas", "scipy", "anndata"):
        try:
            env[mod] = metadata.version(mod)
        except Exception:
            env[mod] = None
    return env


def utc_now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


def resolve_log_path(log_path=None, default=None):
    """Explicit path > METRIC_AUTOPSY_LOG > `default` (None = no log). 'off' disables."""
    candidate = log_path if log_path is not None else os.environ.get("METRIC_AUTOPSY_LOG", default)
    if candidate is None or str(candidate).strip().lower() in ("", "off", "none", "0"):
        return None
    return Path(candidate)


def read_log(path) -> list[dict]:
    path = Path(path)
    if not path.exists():
        return []
    out = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line:
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return out


def append_log(path, record: dict) -> None:
    path = Path(path)
    if path.parent and not path.parent.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(jsonable(record), sort_keys=True, default=str) + "\n")


def jsonable(o):
    """Recursively convert to strict-JSON-safe Python (NaN/inf -> None, numpy -> python)."""
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, set)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return [jsonable(v) for v in o.tolist()]
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        f = float(o)
        return f if np.isfinite(f) else None
    if hasattr(o, "value") and hasattr(o, "name") and not isinstance(o, str):  # enums
        return o.value
    if hasattr(o, "to_dict") and not isinstance(o, type):
        try:
            return jsonable(o.to_dict())
        except Exception:
            return str(o)
    if o is None or isinstance(o, (str, int)):
        return o
    return str(o)
