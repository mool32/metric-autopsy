#!/usr/bin/env python3
"""The extraction of the second diagnostics (validation/exploratory/real_pilot_diagnostics_2_plan.md),
check (ii′): the counts of the full set of genes for the claims' cells, read again from the CELLxGENE
Census (release 2025-11-08), and the controls drawn by the plan's rule. Exploratory, development data;
it does not import the engine (it runs in the selection's environment, requirements-select.txt).

The cells are those of the second run's files (results/real-pilot-2 at bad785a, checked by
real_pilot_3_inputs.sh): the males of case1.npz (claim 1a-males) and the cells of case2.npz (claims
2a and 2b), by their soma_joinid. The genes are every gene the Census marks as measured in every
dataset of the case (real_pilot_select.present_genes), read 10,000 cells at a time and checked to be
raw counts (real_pilot_select.read_blocks). The counts of the claim genes must equal the case files'.

    python validation/exploratory/real_pilot_4_extract.py --data RUN2_DATA --string STRING --out OUT
    python validation/exploratory/real_pilot_4_extract.py --data DATA --string STRING --out OUT --synthetic

In OUT: full_case1.npz and full_case2.npz (the claim genes and the controls' genes as columns, dense
int32, with the cells' obs; case 2 also the cells' totals over the full set of genes), stats_case1.npz
and stats_case2.npz (every measured gene's symbol, soma_joinid, mean count and detection over the
claim's cells), selection.json and extract.log. With --synthetic the same code runs on
real_pilot_select's fabricated Census and fabricated STRING files.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import real_pilot_4_controls as C  # noqa: E402
import real_pilot_select as RS  # noqa: E402  (Census, SyntheticCensus, present_genes, read_blocks; no engine)

CHUNK = 10000


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_case(data_dir: Path, case: dict):
    path = data_dir / case["data"]["file"]
    got = sha256_file(path)
    if got != case["data"]["sha256"]:
        raise ValueError(f"{path.name}: sha256 {got}, select.json has {case['data']['sha256']}")
    z = np.load(path, allow_pickle=False)
    obs = pd.DataFrame({k[4:]: z[k] for k in z.files if k.startswith("obs_")})
    return np.asarray(z["X"], dtype=np.int32), [str(g) for g in z["genes"]], obs, got


def read_full(census, datasets: pd.DataFrame, joinids: np.ndarray, dataset_ids, chunk: int):
    """Every measured gene's counts for the cells `joinids` (sorted), CSR int32, and the genes' var."""
    ids = sorted(set(map(str, dataset_ids)))
    present = RS.present_genes(census, "mus_musculus", datasets, ids)
    X, obs, var = RS.read_blocks(census, "mus_musculus", np.sort(np.asarray(joinids, dtype=np.int64)), present, chunk)
    if X is None:
        raise ValueError("the counts read are not raw (non-negative integers)")
    return X, obs, var, ids, int(len(present))


def same_columns(X_full, symbols, X_case, genes_case, rows_case, names) -> dict:
    """The claim genes' counts in the new extraction against the case file's (the same cells, in order)."""
    out = {}
    for g in names:
        a = np.asarray(X_full[:, symbols.index(g)].todense()).ravel()
        b = X_case[rows_case, genes_case.index(g)]
        out[g] = bool(np.array_equal(a.astype(np.int64), b.astype(np.int64)))
    return out


def obs_arrays(obs: pd.DataFrame) -> dict:
    """The obs as plain arrays (fixed-width strings, int64 joinids), loadable without pickle."""
    return {f"obs_{k}": (obs[k].to_numpy(np.int64) if k == "soma_joinid" else obs[k].astype(str).to_numpy(dtype=str))
            for k in obs.columns}


def dense_columns(X, cols) -> np.ndarray:
    return np.asarray(X[:, cols].todense(), dtype=np.int32)


def run(args, log) -> dict:
    spec_bytes = (args.data / "select.json").read_bytes()
    spec = json.loads(spec_bytes)
    census = RS.SyntheticCensus() if args.synthetic else RS.Census()
    out = dict(synthetic=bool(args.synthetic), release="SYNTHETIC" if args.synthetic else RS.RELEASE,
               select_sha256=hashlib.sha256(spec_bytes).hexdigest(), python=platform.python_version(),
               numpy=np.__version__, pandas=pd.__version__, settings=dict(seed=C.SEED, k=C.K, nearest_pair=C.NEAREST_PAIR,
                                                                         nearest_gene=C.NEAREST_GENE,
                                                                         string_min_score=C.STRING_MIN_SCORE),
               string_sha256={f: sha256_file(args.string / f) for f in C.STRING_FILES}, cases={})
    if not args.synthetic:
        bad = {f: h for f, h in out["string_sha256"].items() if h != C.STRING_SHA256[f]}
        if bad:
            raise ValueError(f"STRING files differ from the first diagnostics': {bad}")
    log("STRING " + "; ".join(f"{k} {v}" for k, v in out["string_sha256"].items())
        + (" (as the first diagnostics')" if not args.synthetic else " (fabricated)"))
    datasets = census.datasets()
    names = C.string_names(args.string / C.STRING_FILES[0])

    # case 1: the males (claim 1a-males)
    t0 = time.time()
    case = spec["cases"]["case1"]
    Xc, genes_c, obs_c, sha = load_case(args.data, case)
    males = np.where(obs_c["sex"].astype(str).to_numpy() == "male")[0]
    order = np.argsort(obs_c["soma_joinid"].to_numpy(np.int64)[males], kind="mergesort")
    rows = males[order]
    joinids = obs_c["soma_joinid"].to_numpy(np.int64)[rows]
    X, obs_r, var, ids, n_present = read_full(census, datasets, joinids, case["data"]["datasets"], args.chunk)
    if not np.array_equal(obs_r["soma_joinid"].to_numpy(np.int64), joinids):
        raise ValueError("case 1: the cells read are not the case file's males, in order")
    symbols = var["feature_name"].astype(str).tolist()
    mean, det = C.gene_stats(X)
    check = same_columns(X, symbols, Xc, genes_c, rows, list(C.PAIR))
    log(f"## case 1, males: {len(rows)} cells x {n_present} measured genes ({time.time() - t0:.0f} s); the claim genes' "
        f"counts as case1.npz's: {check}")
    if not all(check.values()):
        raise ValueError(f"case 1: the claim genes' counts differ from case1.npz's: {check}")
    sel = C.draw_pairs(symbols, mean, det, list(C.PAIR), names, args.string / C.STRING_FILES[1])
    sel["string_mapped_genes"] = int(sum(s in names for s in symbols))
    log(f"1a-males: {len(sel['controls'])} pairs drawn from {sel['eligible']} eligible ({sel['candidates']} candidates; "
        f"{sel['linked_excluded']} pairs excluded by STRING >= {C.STRING_MIN_SCORE}; {sel['string_mapped_genes']} of "
        f"{len(symbols)} genes mapped to STRING); matched: {sel['matched']}")
    cols = sorted({symbols.index(g) for g in C.PAIR} | {j for c in sel["columns"] for j in c})
    keep = [symbols[j] for j in cols]
    obs_keep = obs_c.iloc[rows][["mouse", "sex", "age", "soma_joinid"]].reset_index(drop=True)
    np.savez_compressed(args.out / "full_case1.npz", X=dense_columns(X, cols), genes=np.asarray(keep),
                        **obs_arrays(obs_keep))
    np.savez_compressed(args.out / "stats_case1.npz", symbol=np.asarray(symbols), soma_joinid=var["soma_joinid"].to_numpy(),
                        mean=mean, detection=det)
    out["cases"]["case1"] = dict(claim="1a-males", case_sha256=sha, cells=int(len(rows)), measured_genes=n_present,
                                 datasets=ids, claim_counts_as_case_file=check, columns_saved=len(cols),
                                 selection={"1a-males": {k: v for k, v in sel.items() if k != "columns"}},
                                 seconds=round(time.time() - t0, 1))
    del X

    # case 2: every cell (claims 2a and 2b)
    t0 = time.time()
    case = spec["cases"]["case2"]
    Xc, genes_c, obs_c, sha = load_case(args.data, case)
    order = np.argsort(obs_c["soma_joinid"].to_numpy(np.int64), kind="mergesort")
    joinids = obs_c["soma_joinid"].to_numpy(np.int64)[order]
    X, obs_r, var, ids, n_present = read_full(census, datasets, joinids, case["data"]["datasets"], args.chunk)
    if not np.array_equal(obs_r["soma_joinid"].to_numpy(np.int64), joinids):
        raise ValueError("case 2: the cells read are not the case file's, in order")
    symbols = var["feature_name"].astype(str).tolist()
    mean, det = C.gene_stats(X)
    check = same_columns(X, symbols, Xc, genes_c, order, ["Xist", *C.Y_GENES])
    log(f"## case 2: {len(order)} cells x {n_present} measured genes ({time.time() - t0:.0f} s); the claim genes' counts "
        f"as case2.npz's: {check}")
    if not all(check.values()):
        raise ValueError(f"case 2: the claim genes' counts differ from case2.npz's: {check}")
    s2a = C.draw_genes(symbols, mean, det, "Xist")
    s2b = C.draw_sets(symbols, mean, det, list(C.Y_GENES))
    log(f"2a: {len(s2a['controls'])} genes drawn ({s2a['candidates']} candidates); matched: {s2a['matched']}",
        f"2b: {len(s2b['controls'])} sets drawn in {s2b['draws']} draws ({s2b['candidates']} candidates); matched: "
        f"{s2b['matched']}")
    cols = sorted({symbols.index(g) for g in ["Xist", *C.Y_GENES]} | {j for c in s2a["columns"] + s2b["columns"] for j in c})
    keep = [symbols[j] for j in cols]
    tot = np.asarray(X.sum(axis=1), dtype=np.int64).ravel()
    obs_keep = obs_c.iloc[order][["mouse", "sex", "soma_joinid"]].reset_index(drop=True)
    np.savez_compressed(args.out / "full_case2.npz", X=dense_columns(X, cols), genes=np.asarray(keep), tot=tot,
                        **obs_arrays(obs_keep))
    np.savez_compressed(args.out / "stats_case2.npz", symbol=np.asarray(symbols), soma_joinid=var["soma_joinid"].to_numpy(),
                        mean=mean, detection=det)
    out["cases"]["case2"] = dict(claims=["2a", "2b"], case_sha256=sha, cells=int(len(order)), measured_genes=n_present,
                                 datasets=ids, claim_counts_as_case_file=check, columns_saved=len(cols),
                                 selection={"2a": {k: v for k, v in s2a.items() if k != "columns"},
                                            "2b": {k: v for k, v in s2b.items() if k != "columns"}},
                                 seconds=round(time.time() - t0, 1))
    census.close()
    for f in ("full_case1.npz", "full_case2.npz", "stats_case1.npz", "stats_case2.npz"):
        out.setdefault("files", {})[f] = sha256_file(args.out / f)
    return out


def check_census(args) -> int:
    """Metadata only: the Census opens at the release, and each case's datasets have measured genes."""
    spec = json.loads((args.data / "select.json").read_text())
    census = RS.SyntheticCensus() if args.synthetic else RS.Census()
    datasets = census.datasets()
    out = dict(release="SYNTHETIC" if args.synthetic else RS.RELEASE, summary=census.summary(), cases={})
    for name in ("case1", "case2"):
        ids = sorted(set(map(str, spec["cases"][name]["data"]["datasets"])))
        present = RS.present_genes(census, "mus_musculus", datasets, ids)
        out["cases"][name] = dict(datasets=ids, measured_genes=int(len(present)))
    census.close()
    (args.out / "census_check.json").write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps(out, indent=1, default=str))
    return 0 if all(c["measured_genes"] > 0 for c in out["cases"].values()) else 1


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--data", type=Path, required=True)
    p.add_argument("--string", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--synthetic", action="store_true", help="the fabricated Census of real_pilot_select (offline)")
    p.add_argument("--chunk", type=int, default=CHUNK)
    p.add_argument("--check-census", action="store_true",
                   help="the dry run's check: open the Census and count each case's measured genes, read no counts")
    args = p.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    if args.check_census:
        return check_census(args)
    lines = []

    def log(*msg):
        for m in msg:
            print(m, flush=True)
            lines.append(str(m))
    tag = "SYNTHETIC " if args.synthetic else ""
    log(f"# {tag}The second diagnostics' extraction (validation/exploratory/real_pilot_diagnostics_2_plan.md, check (ii′))",
        "", f"python {platform.python_version()}; numpy {np.__version__}; pandas {pd.__version__}; Census release "
        f"{'SYNTHETIC' if args.synthetic else RS.RELEASE}", "")
    t0 = time.time()
    status = 0
    try:
        out = run(args, log)
    except Exception as exc:  # reported as it is; the analyses do not start
        out = dict(error=repr(exc), traceback=traceback.format_exc().splitlines())
        log(f"EXCEPTION {exc!r}", *out["traceback"])
        status = 1
    out["seconds"] = round(time.time() - t0, 1)
    log("", f"{tag}extraction: {out['seconds']} s")
    (args.out / "selection.json").write_text(json.dumps(out, indent=1, default=str))
    (args.out / "extract.log").write_text("\n".join(lines) + "\n")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
