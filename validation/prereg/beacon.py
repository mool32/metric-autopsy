"""The public randomness that sets the key of the confirmatory panel (validation/prereg/v1.md, 3.3).

The key is the randomness of one round of drand's *quicknet* chain, run by the League of
Entropy: every 3 seconds a threshold of independent organisations signs the round number with a
BLS signature, and the round's randomness is the sha256 of that signature. Nobody knows a
future round's value in advance, no single member can bias it, and anyone can check a round with
the chain's fixed public key. The run tag names a round at least an hour after the tag is pushed;
the workflow waits for it, fetches it from the public relays, verifies the signature, and uses
the randomness (64 hex characters) as the key (``panel.assign``).

Checking a round by hand: ``curl https://api.drand.sh/<CHAIN_HASH>/public/<ROUND>`` (or any relay
of RELAYS) returns its signature and randomness; ``python beacon.py verify --round R --signature S``
checks the signature against the public key below and prints the randomness, which must equal
the key in the results' manifest. The chain's constants can be compared with
``curl https://api.drand.sh/<CHAIN_HASH>/info`` and the League of Entropy's documentation.

Dependencies: py_ecc (pure-Python BLS12-381; requirements-beacon.txt) for the signature, the
standard library for the rest. This file never imports the engine.

    python validation/prereg/beacon.py round --delay 3600           # the round to name in the run tag
    python validation/prereg/beacon.py wait --round R --out key.json  # wait, fetch, verify, write
    python validation/prereg/beacon.py verify --round R --signature HEX
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
import urllib.request
from datetime import datetime, timezone

# drand quicknet (League of Entropy mainnet): unchained rounds, signatures on G1 (RFC 9380).
CHAIN = dict(
    name="quicknet",
    hash="52db9ba70e0cc0f6eaf7803dd07447a1f5477735fd3f661792ba94600c84e971",
    public_key=("83cf0f2896adee7eb8b5f01fcad3912212c437e0073e911fb90022d3e760183c8c4b450b6a0a6c3ac6a5776a2d1064510d"
                "1fec758c921cc22b0e17e63aaf4bcb5ed66304de9cf809bd274ca73bab4af5a6e9c76a4bc09e76eae8991ef5ece45a"),
    genesis_time=1692803367,
    period=3,
    scheme="bls-unchained-g1-rfc9380",
)
DST_G1 = b"BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_NUL_"
RELAYS = ("https://api.drand.sh", "https://api2.drand.sh", "https://api3.drand.sh",
          "https://drand.cloudflare.com")
MIN_DELAY = 3600     # the named round lies at least this many seconds after the run tag's push
TAG_DELAY = 3900     # the workflow names the first round this far ahead: MIN_DELAY and a margin
LATEST_BEFORE = 600  # before the round is named (nothing lost yet): how long the relays' newest round is asked for
LATEST_AFTER = 240   # after the push: inside the margin TAG_DELAY - MIN_DELAY, so the gap can still hold
FETCH_HOURS = 5.0    # how long the workflow keeps asking the relays after the round's time


def round_time(rnd: int, chain: dict = CHAIN) -> int:
    """Unix time at which round `rnd` is produced."""
    return chain["genesis_time"] + (int(rnd) - 1) * chain["period"]


def first_round_at(t: float, chain: dict = CHAIN) -> int:
    """The first round produced at or after unix time `t`."""
    return max(1, int(math.ceil((t - chain["genesis_time"]) / chain["period"])) + 1)


def round_for_tag(push_time: float, delay: int = MIN_DELAY, chain: dict = CHAIN) -> int:
    """The round the run tag names: the first one at least `delay` seconds after `push_time`.
    The workflow computes it just before the push and adds a margin, then checks the push time."""
    return first_round_at(push_time + delay, chain)


def message(rnd: int) -> bytes:
    """What quicknet signs: sha256 of the round number as 8 big-endian bytes (unchained)."""
    return hashlib.sha256(int(rnd).to_bytes(8, "big")).digest()


def randomness_of(signature_hex: str) -> str:
    """A round's randomness: sha256 of its signature's bytes, as 64 lowercase hex characters."""
    return hashlib.sha256(bytes.fromhex(signature_hex)).hexdigest()


def verify(rnd: int, signature_hex: str, public_key_hex: str = CHAIN["public_key"]) -> bool:
    """BLS check of a quicknet round: e(signature, g2) == e(H(message), public key), with the
    signature (48 bytes) on G1, the public key (96 bytes) on G2 and H the RFC 9380 hash to G1."""
    from py_ecc.bls.hash_to_curve import hash_to_G1
    from py_ecc.bls.point_compression import decompress_G1, decompress_G2
    from py_ecc.optimized_bls12_381 import G2, curve_order, is_inf, multiply, pairing

    sig_b, pk_b = bytes.fromhex(signature_hex), bytes.fromhex(public_key_hex)
    if len(sig_b) != 48 or len(pk_b) != 96:
        return False
    try:
        sig = decompress_G1(int.from_bytes(sig_b, "big"))
        pk = decompress_G2((int.from_bytes(pk_b[:48], "big"), int.from_bytes(pk_b[48:], "big")))
    except Exception:
        return False
    if is_inf(sig) or is_inf(pk) or not is_inf(multiply(sig, curve_order)):
        return False
    h = hash_to_G1(message(rnd), DST_G1, hashlib.sha256)
    return pairing(G2, sig) == pairing(pk, h)


def _get(url: str, timeout: float = 20.0) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "metric-autopsy-panel"})
    with urllib.request.urlopen(req, timeout=timeout) as fh:
        return json.loads(fh.read().decode())


def chain_info(relay: str, chain: dict = CHAIN) -> dict:
    """The chain's parameters as a relay serves them; compared with CHAIN by `check_chain`."""
    return _get(f"{relay}/{chain['hash']}/info")


def check_chain(info: dict, chain: dict = CHAIN) -> list[str]:
    """Differences between a relay's chain info and the constants above (empty if none)."""
    want = dict(public_key=chain["public_key"], period=chain["period"], genesis_time=chain["genesis_time"],
                hash=chain["hash"], schemeID=chain["scheme"])
    return [f"{k}: relay {info.get(k)!r}, expected {v!r}" for k, v in want.items() if info.get(k) != v]


def fetch(rnd: int, relays=RELAYS, chain: dict = CHAIN) -> dict:
    """Round `rnd` from the relays: every relay is asked; the round is accepted when at least
    one answer carries a signature that verifies, and every verified answer must agree.
    Returns the key record (round, randomness, signature, the relays' answers)."""
    answers = []
    for relay in relays:
        try:
            got = _get(f"{relay}/{chain['hash']}/public/{int(rnd)}")
        except Exception as exc:  # a relay that is down or late is recorded, not fatal
            answers.append(dict(relay=relay, error=repr(exc)))
            continue
        sig = str(got.get("signature", ""))
        ok = int(got.get("round", -1)) == int(rnd) and verify(rnd, sig, chain["public_key"])
        answers.append(dict(relay=relay, round=got.get("round"), signature=sig,
                            randomness=got.get("randomness"), verified=ok,
                            randomness_matches=ok and got.get("randomness") == randomness_of(sig)))
    good = [a for a in answers if a.get("verified")]
    if not good:
        raise RuntimeError(f"round {rnd}: no relay returned a verified signature: {answers}")
    if len({a["signature"] for a in good}) != 1:
        raise RuntimeError(f"round {rnd}: verified answers disagree: {answers}")
    sig = good[0]["signature"]
    return dict(chain=chain["name"], chain_hash=chain["hash"], public_key=chain["public_key"],
                round=int(rnd), round_time_utc=_utc(round_time(rnd, chain)), signature=sig,
                randomness=randomness_of(sig), key=randomness_of(sig), answers=answers,
                fetched_utc=_utc(time.time()))


def latest(relays=RELAYS, chain: dict = CHAIN, within: float = 0.0, step: float = 10.0) -> int:
    """The newest round the relays serve whose signature verifies: the beacon's own clock, so the
    gap between the run tag's push and its round need not rest on the runner's clock alone (the
    fourth review). Rounds that exist when it is asked include every round that existed at the
    push before it. While no relay answers with a verified round, the relays are asked again every
    `step` seconds for up to `within` seconds (the fifth review: one try after the push could lose
    a named run to a short outage); a later answer only makes the bound on the gap stricter."""
    deadline = time.time() + within
    while True:
        best = None
        for relay in relays:
            try:
                got = _get(f"{relay}/{chain['hash']}/public/latest")
                rnd = int(got.get("round", -1))
                if rnd > 0 and verify(rnd, str(got.get("signature", "")), chain["public_key"]):
                    best = rnd if best is None else max(best, rnd)
            except Exception:  # a relay that is down is skipped; one that answers must verify
                continue
        if best is not None:
            return best
        if time.time() + step > deadline:
            raise RuntimeError(f"no relay returned a verified latest round within {within:.0f} s")
        time.sleep(step)


def gap_after_push(rnd: int, pushed: float, newest: int, chain: dict = CHAIN) -> dict:
    """How far round `rnd` lies after the run tag's push, by the runner's clock (`pushed`) and by
    the newest round that existed after the push (`newest`, ``latest``): the smaller counts."""
    by_clock = round_time(rnd, chain) - pushed
    by_beacon = (int(rnd) - int(newest)) * chain["period"]
    return dict(round=int(rnd), newest=int(newest), by_clock=by_clock, by_beacon=by_beacon,
                seconds=min(by_clock, by_beacon))


def wait_and_fetch(rnd: int, relays=RELAYS, chain: dict = CHAIN, hours: float = FETCH_HOURS,
                   step: float = 30.0) -> dict:
    """Sleep until round `rnd` is due, then ask the relays until one verified answer arrives or
    `hours` have passed since the round's time."""
    due = round_time(rnd, chain)
    while time.time() < due + 2 * chain["period"]:
        time.sleep(min(60.0, max(1.0, due + 2 * chain["period"] - time.time())))
    last = None
    while time.time() < due + hours * 3600:
        try:
            return fetch(rnd, relays, chain)
        except RuntimeError as exc:
            last = exc
            time.sleep(step)
    raise RuntimeError(f"round {rnd} not obtained within {hours} h of its time: {last}")


def _utc(t: float) -> str:
    return datetime.fromtimestamp(t, tz=timezone.utc).isoformat(timespec="seconds")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("round", help="the first round at least DELAY seconds from now")
    r.add_argument("--delay", type=int, default=MIN_DELAY)
    w = sub.add_parser("wait", help="wait for a round, fetch and verify it, write the key record")
    w.add_argument("--round", type=int, required=True)
    w.add_argument("--out", required=True)
    w.add_argument("--hours", type=float, default=FETCH_HOURS)
    v = sub.add_parser("verify", help="check a round's signature and print its randomness")
    v.add_argument("--round", type=int, required=True)
    v.add_argument("--signature", required=True)
    i = sub.add_parser("info", help="compare the relays' chain info with the constants")
    i.add_argument("--relay", action="append")
    args = p.parse_args(argv)
    if args.cmd == "round":
        rnd = round_for_tag(time.time(), args.delay)
        print(json.dumps(dict(chain=CHAIN["name"], chain_hash=CHAIN["hash"], round=rnd,
                              round_time_utc=_utc(round_time(rnd)))))
    elif args.cmd == "wait":
        rec = wait_and_fetch(args.round, hours=args.hours)
        with open(args.out, "w") as fh:
            json.dump(rec, fh, indent=1)
        print(json.dumps({k: rec[k] for k in ("chain", "round", "round_time_utc", "randomness")}))
    elif args.cmd == "verify":
        ok = verify(args.round, args.signature)
        print(json.dumps(dict(round=args.round, verified=ok,
                              randomness=randomness_of(args.signature) if ok else None)))
        raise SystemExit(0 if ok else 1)
    else:
        bad = 0
        for relay in args.relay or RELAYS:
            try:
                diff = check_chain(chain_info(relay))
            except Exception as exc:
                diff = [f"unreachable: {exc!r}"]
            bad += bool(diff)
            print(json.dumps(dict(relay=relay, matches=not diff, differences=diff)))
        raise SystemExit(1 if bad == len(args.relay or RELAYS) else 0)


if __name__ == "__main__":
    main()
