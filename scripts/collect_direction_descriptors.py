#!/usr/bin/env python3
"""Collect the ISOTROPY program's **official** direction descriptors (milestone R7).

`iso`'s `SHOW DIRECTION VECTOR` prints, for every order-parameter direction of an
irrep, the program's own component string (e.g. ``(a,0,0)``, ``(a,a,0)``).  Those
strings are the ones users copy out of ISOTROPY, while the bundled tables carry
cryspglib's *internal* notation for the cases where the program's string depends
on the irrep (``dim = 2``) or prints per-irrep component lists (``dim >= 4``).

This script walks every machine record ``(parent space group, Miller-Love
label)`` and records both notations plus the record's source hashes, then writes
a JSON that the Rust side freezes into a table.  It is deliberately
**resumable** and **per-record isolated**: each record is one `iso` session,
results are appended to a JSONL cache as they arrive, and one failing record is
reported rather than dropping the whole run.

Usage::

    python3 scripts/collect_direction_descriptors.py            # collect all
    python3 scripts/collect_direction_descriptors.py --limit 20 # smoke test
    python3 scripts/collect_direction_descriptors.py --jobs 8

The cache defaults to ``target/direction-descriptors/cache.jsonl`` and the JSON
to ``scripts/data/direction_descriptors_v1.json``.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import sys
import threading
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import verify_isotropy_oracle as oracle  # noqa: E402

SCHEMA = "direction-descriptors-v1"


def sha256_of(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_hashes() -> dict:
    """The files whose content the collected strings depend on."""
    sources = {
        "iso": os.path.join(oracle.ISO_DIR, "iso"),
        "data_irreps.txt": os.path.join(oracle.ISO_DIR, "data_irreps.txt"),
        "data_isotropy.txt": os.path.join(oracle.ISO_DIR, "data_isotropy.txt"),
    }
    return {name: sha256_of(path) for name, path in sources.items()}


def record_key(sg: int, ml: str) -> str:
    return f"{sg}|{ml}"


def collect_one(sg: int, ml: str) -> dict:
    """One record: the program's direction strings, or why they are missing."""
    started = time.time()
    try:
        official = oracle.run_oracle_direction_vectors(sg, ml)
    except Exception as error:  # per-record isolation: never drop the run
        return {
            "key": record_key(sg, ml),
            "sg": sg,
            "label": ml,
            "status": "failed",
            "error": f"{type(error).__name__}: {error}",
            "seconds": round(time.time() - started, 3),
        }
    return {
        "key": record_key(sg, ml),
        "sg": sg,
        "label": ml,
        "status": "ok" if official else "empty",
        "official": official,
        "seconds": round(time.time() - started, 3),
    }


def load_cache(path: str) -> dict:
    done = {}
    if not os.path.isfile(path):
        return done
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            entry = json.loads(line)
            done[entry["key"]] = entry
    return done


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", default="target/direction-descriptors/cache.jsonl")
    parser.add_argument("--output", default="scripts/data/direction_descriptors_v1.json")
    parser.add_argument("--manifest", default="")
    parser.add_argument("--jobs", type=int, default=0, help="0 = min(8, cores)")
    parser.add_argument("--limit", type=int, default=0, help="collect at most N records")
    parser.add_argument("--retry-failed", action="store_true")
    args = parser.parse_args()

    records = oracle.machine_records()
    keys = sorted(records)
    if args.limit:
        keys = keys[: args.limit]

    os.makedirs(os.path.dirname(os.path.abspath(args.cache)), exist_ok=True)
    done = load_cache(args.cache)
    todo = [
        (sg, ml)
        for (sg, ml) in keys
        if record_key(sg, ml) not in done
        or (args.retry_failed and done[record_key(sg, ml)]["status"] == "failed")
    ]
    jobs = args.jobs or min(8, os.cpu_count() or 1)
    print(
        f"{len(keys)} machine records, {len(done)} cached, {len(todo)} to collect, {jobs} jobs",
        flush=True,
    )

    write_lock = threading.Lock()
    with open(args.cache, "a", encoding="utf-8") as cache:
        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = {
                pool.submit(collect_one, sg, ml): (sg, ml) for (sg, ml) in todo
            }
            finished = 0
            for future in concurrent.futures.as_completed(futures):
                entry = future.result()
                with write_lock:
                    cache.write(json.dumps(entry, sort_keys=True) + "\n")
                    cache.flush()
                    done[entry["key"]] = entry
                finished += 1
                if finished % 25 == 0 or finished == len(todo):
                    print(f"  {finished}/{len(todo)}", flush=True)

    # Assemble the JSON in machine-record order, with the internal notation
    # alongside the official one so the two can be compared downstream.
    out_records = {}
    status_counts = {"ok": 0, "empty": 0, "failed": 0}
    directions = 0
    for (sg, ml) in keys:
        entry = done.get(record_key(sg, ml))
        status_counts[entry["status"]] = status_counts.get(entry["status"], 0) + 1
        official = entry.get("official", {})
        directions += len(official)
        internal = {
            row["label"]: row["direction"] for row in records[(sg, ml)]
        }
        out_records[record_key(sg, ml)] = {
            "official": official,
            "internal": internal,
            "status": entry["status"],
        }
        if entry["status"] == "failed":
            out_records[record_key(sg, ml)]["error"] = entry.get("error", "")

    payload = {
        "schema": SCHEMA,
        "sources": source_hashes(),
        "counts": {
            "records": len(keys),
            "directions": directions,
            **status_counts,
        },
        "records": out_records,
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=1, sort_keys=True)
        handle.write("\n")
    manifest_path = args.manifest or args.output.replace(".json", ".manifest.json")
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "schema": SCHEMA,
                "sources": payload["sources"],
                "counts": payload["counts"],
                "output": os.path.basename(args.output),
                "sha256": sha256_of(args.output),
            },
            handle,
            indent=1,
            sort_keys=True,
        )
        handle.write("\n")
    print(
        f"wrote {args.output}: {payload['counts']} (manifest {manifest_path})",
        flush=True,
    )
    return 0 if status_counts["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
