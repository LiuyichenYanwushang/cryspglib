#!/usr/bin/env python3
"""Task 9 ordinary-table setting census against the pinned ISOTROPY oracle.

For every (parent space group, Miller-Love irrep) in the pinned
``isotropy_subgroup/iso.zip`` this collector runs the archived ``iso`` binary
once, in its own private working directory, with the fixed setting command
``SET I ALL OR 1``, and compares the printed subgroup geometry with the pinned
machine table ``data_isotropy.txt``.

For every machine record it derives the candidate change of basis

    U = W . P_parent . (P_child . B)^-1

(``W`` = stored ``isotropy_basis`` in the parent primitive frame, ``P_*`` =
parent/child primitive bases, ``B`` = basis the program prints).  A record is a
``candidate`` only when ``U`` is integral with ``|det U| = 1`` **and** the
converted stored origin is exactly the printed origin.  This is a candidate
census: it does not resolve deltas, choose settings, prove canonical operation
matching or test engine coverage.

The archive, machine data and binary are pinned by SHA256; binary and data are
extracted into a private temporary directory, never into the repository.  Every
machine record produces exactly one output record, including empty official
tables and failed queries, and the summary partitions are asserted to sum.

Usage::

    python3 scripts/audit_subduction_settings.py --output /tmp/task9.jsonl
    python3 scripts/audit_subduction_settings.py --output /tmp/smoke.json --parent 43 --parent 3

Exit codes: 0 = the requested census was collected completely (candidate
coverage is reported separately and is *not* success); 1 = parse/inconsistency
failure; 2 = usage/archive/IO failure.
"""

import argparse
import collections
import concurrent.futures
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
import zipfile
from fractions import Fraction

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(SCRIPT_DIR)
ISO_DIR = os.path.join(REPO_DIR, "isotropy_subgroup")
ARCHIVE = os.path.join(ISO_DIR, "iso.zip")

if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

import verify_isotropy_oracle as geometry  # noqa: E402
import verify_isotropy_operations as operations  # noqa: E402

PINNED_ZIP_SHA256 = "568667bfc8027095537d642297b319c872d00016b868143c666f90d5931d9f7b"
PINNED_ISO_VERSION = "9.6.1"
EXPECTED_IRREPS = 4777
EXPECTED_RECORDS = 15239
ORACLE_TIMEOUT_SECONDS = 300
PROGRESS_EVERY = 250
MAX_ERRORS_IN_SUMMARY = 200

SETTING_COMMAND = "SET I ALL OR 1"
ORACLE_COMMANDS = (
    "PAGE 1000",
    "SC 250",
    SETTING_COMMAND,
    "VALUE PARENT {parent}",
    "VALUE IRREP {irrep}",
    "SHOW SUBGROUP",
    "SHOW BASIS",
    "SHOW ORIGIN",
    "SHOW SIZE",
    "SHOW DIRECTION",
    "DISPLAY ISOTROPY",
    "QUIT",
)
# Everything the archived binary needs at run time; extracted from the archive
# so that no state outside the pinned SHA256 can influence the census.
ORACLE_MEMBERS = (
    "iso", "const.dat", "data_diperiodic.txt", "data_images.txt", "data_irreps.txt",
    "data_isotropy.txt", "data_little.txt", "data_magnetic.txt", "data_space.txt",
    "data_ssg.txt", "data_ssgmag.txt", "data_wyckoff.txt",
)
BANNER_TAIL = (
    "Harold T. Stokes, Dorian M. Hatch, and Branton J. Campbell",
    "Brigham Young University",
    "Current setting is International (new ed.) with conventional basis vectors.",
)
HEADER_RE = re.compile(r"^\*+Subgroup\s+Size\s+Dir\s+Basis Vectors\s+Origin\s*$")
STAR_RE = re.compile(r"^\*+$")
ROW_RE = re.compile(
    r"^(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+"
    r"\(([^()]*)\),\(([^()]*)\),\(([^()]*)\)\s+\(([^()]*)\)$"
)
TOKEN_RE = re.compile(r"^-?\d+(?:/\d+)?$")
CLASSIFIED = ("candidate", "nonunimodular", "origin_mismatch")
UNRETURNED = "unreturned_empty"


class CollectorError(Exception):
    """Fatal archive/usage problem (exit 2)."""


class OracleParseError(Exception):
    """The pinned oracle printed something outside its pinned table grammar."""


def sha256(payload):
    return hashlib.sha256(payload).hexdigest()


def archive_manifest():
    if not os.path.isfile(ARCHIVE):
        raise CollectorError(
            f"pinned archive missing: {ARCHIVE}\n"
            "  restore the repository copy of the pinned archive before running"
        )
    digest = operations.file_digest(ARCHIVE)
    if digest != PINNED_ZIP_SHA256:
        raise CollectorError(
            f"archive {ARCHIVE} has SHA256 {digest}, expected {PINNED_ZIP_SHA256}; "
            "refusing to run an unpinned oracle"
        )
    return {
        "path": os.path.relpath(ARCHIVE, REPO_DIR),
        "sha256": digest,
        "size": os.path.getsize(ARCHIVE),
    }


def extract_oracle(root):
    """Extract the pinned binary/data into ``root`` and return their hashes."""
    target = os.path.join(root, "oracle")
    os.makedirs(target)
    manifest = {}
    with zipfile.ZipFile(ARCHIVE) as archive:
        names = set(archive.namelist())
        missing = [name for name in ORACLE_MEMBERS if name not in names]
        if missing:
            raise CollectorError(f"pinned archive lacks required members: {missing}")
        for name in ORACLE_MEMBERS:
            payload = archive.read(name)
            with open(os.path.join(target, name), "wb") as handle:
                handle.write(payload)
            manifest[name] = {"size": len(payload), "sha256": sha256(payload)}
    os.chmod(os.path.join(target, "iso"), 0o755)
    return target, manifest


def load_inventory():
    """The 4777 (parent SG, 1-based irrep index, Miller-Love label) queries."""
    lines = geometry.read_file("data_irreps.txt")
    sections = geometry.get_sections(lines)
    labels = geometry.parse_labels(lines, sections, "irrep_label")
    parents = geometry.parse_ints(lines, sections, "irrep_space_group")
    if not labels or len(labels) != len(parents):
        raise CollectorError("data_irreps.txt has inconsistent label/parent arrays")
    if len(labels) != EXPECTED_IRREPS:
        raise CollectorError(
            f"data_irreps.txt has {len(labels)} irreps, expected {EXPECTED_IRREPS}"
        )
    return [
        {"parent": int(parent), "index": position + 1, "irrep": label}
        for position, (parent, label) in enumerate(zip(parents, labels))
    ]


def load_machine_records(inventory):
    """All machine isotropy rows with their true ordinal key, in file order.

    The row payload comes from the existing geometry gate reader
    (``verify_isotropy_oracle.machine_records``); only the ordinal key --
    parent / source irrep index / direction label / child -- is added here.
    """
    grouped = geometry.machine_records()
    lines = geometry.read_file("data_isotropy.txt")
    sections = geometry.get_sections(lines)
    parents = geometry.parse_ints(lines, sections, "isotropy_parent")
    irreps = geometry.parse_ints(lines, sections, "isotropy_irrep")
    if len(parents) != len(irreps) or len(parents) != EXPECTED_RECORDS:
        raise CollectorError(
            f"data_isotropy.txt has {len(parents)} records, expected {EXPECTED_RECORDS}"
        )
    by_index = {entry["index"]: entry for entry in inventory}
    if {(entry["parent"], entry["irrep"]) for entry in inventory} != set(grouped):
        raise CollectorError("data_irreps.txt and data_isotropy.txt disagree on the inventory")
    cursor = {}
    records = []
    for offset in range(len(parents)):
        entry = by_index.get(irreps[offset])
        if entry is None or entry["parent"] != parents[offset]:
            raise CollectorError(f"isotropy record {offset + 1} references an unknown irrep")
        key = (entry["parent"], entry["irrep"])
        rows = grouped[key]
        position = cursor.get(key, 0)
        if position >= len(rows):
            raise CollectorError(f"more isotropy records than gate rows for {key}")
        cursor[key] = position + 1
        row = rows[position]
        records.append(
            {
                "ordinal": offset,
                "parent": entry["parent"],
                "source_irrep_index": irreps[offset],
                "irrep": entry["irrep"],
                "direction_label": row["label"],
                "child": row["subgroup"],
                "size": row["size"],
                "dim": row["dim"],
                "direction": row["direction"],
                "basis": row["basis"],
                "origin": row["origin"],
            }
        )
    if any(cursor.get(key, 0) != len(rows) for key, rows in grouped.items()):
        raise CollectorError("gate reader returned rows no isotropy record claims")
    keys = {
        (r["parent"], r["source_irrep_index"], r["direction_label"], r["child"])
        for r in records
    }
    if len(keys) != len(records):
        raise CollectorError("machine ordinal key (parent, irrep index, label, child) is not unique")
    seen = collections.defaultdict(set)
    for record in records:
        query = (record["parent"], record["source_irrep_index"])
        if record["direction_label"] in seen[query]:
            raise CollectorError(f"{query}: duplicate machine direction label")
        seen[query].add(record["direction_label"])
    return records


def rational_vector(text):
    parts = text.split(",")
    if len(parts) != 3:
        raise OracleParseError(f"expected 3 rational components, got {text!r}")
    values = []
    for part in parts:
        token = part.strip()
        if not TOKEN_RE.match(token):
            raise OracleParseError(f"malformed rational component {part!r}")
        values.append(Fraction(token))
    return values


def parse_oracle_table(stdout):
    """Parse one oracle transcript; anything outside the pinned grammar fails.

    An empty official table prints a star-only marker and no header; a
    non-empty table prints the column header followed by rows.  Malformed rows,
    a second header (pagination) or stray text are hard errors, never skipped.
    """
    lines = stdout.splitlines()
    if len(lines) < 5 or not lines[0].startswith("Isotropy, Version "):
        raise OracleParseError("missing or malformed version banner")
    version = lines[0].split("Version", 1)[1].split(",")[0].strip()
    if version != PINNED_ISO_VERSION:
        raise OracleParseError(f"oracle version {version!r}, expected {PINNED_ISO_VERSION!r}")
    if tuple(lines[1:4]) != BANNER_TAIL:
        raise OracleParseError(f"unexpected banner lines: {lines[1:4]!r}")
    rows = []
    header_seen = False
    star_seen = False
    for raw in lines[4:]:
        line = raw.rstrip()
        if not line:
            continue
        if STAR_RE.match(line):
            star_seen = True
            continue
        if not header_seen and HEADER_RE.match(line):
            header_seen = True
            continue
        if not header_seen:
            raise OracleParseError(f"unexpected line before table header: {line!r}")
        match = ROW_RE.match(line)
        if not match:
            raise OracleParseError(f"malformed table row: {line!r}")
        rows.append(
            {
                "child": int(match.group(1)),
                "symbol": match.group(2),
                "size": int(match.group(3)),
                "label": match.group(4),
                "basis": [
                    rational_vector(match.group(5)),
                    rational_vector(match.group(6)),
                    rational_vector(match.group(7)),
                ],
                "origin": rational_vector(match.group(8)),
                "raw": line,
            }
        )
    if header_seen and not rows:
        raise OracleParseError("table header without rows")
    if not header_seen and not star_seen:
        raise OracleParseError("neither a table nor an empty-table marker")
    return rows


def run_query(binary, data_dir, work_root, query):
    commands = [
        command.format(parent=query["parent"], irrep=query["irrep"])
        for command in ORACLE_COMMANDS
    ]
    env = dict(os.environ, ISODATA=data_dir + os.sep)
    with tempfile.TemporaryDirectory(dir=work_root) as workdir:
        try:
            result = subprocess.run(
                [binary],
                input="\n".join(commands) + "\n",
                capture_output=True,
                text=True,
                cwd=workdir,
                env=env,
                timeout=ORACLE_TIMEOUT_SECONDS,
            )
        except subprocess.TimeoutExpired:
            return {"error": {"kind": "timeout", "detail": f">{ORACLE_TIMEOUT_SECONDS}s"},
                    "rows": [], "stdout_tail": ""}
        except OSError as error:
            return {"error": {"kind": "spawn_failed", "detail": str(error)},
                    "rows": [], "stdout_tail": ""}
    if result.returncode != 0:
        return {
            "error": {"kind": "nonzero_exit",
                      "detail": f"rc={result.returncode} stderr={result.stderr.strip()[:200]}"},
            "rows": [], "stdout_tail": "",
        }
    try:
        rows = parse_oracle_table(result.stdout)
    except OracleParseError as error:
        return {"error": {"kind": "parse_error", "detail": str(error)},
                "rows": [], "stdout_tail": "\n".join(result.stdout.splitlines()[-4:])}
    return {"error": None, "rows": rows, "stdout_tail": ""}


def base_record(row, parent):
    return {
        "type": "record",
        "ordinal": row["ordinal"],
        "parent": parent,
        "source_irrep_index": row["source_irrep_index"],
        "irrep": row["irrep"],
        "direction_label": row["direction_label"],
        "child": row["child"],
        "size": row["size"],
        "dim": row["dim"],
        "direction": row["direction"],
        "parent_centering": geometry.CENTERING_LETTER[parent],
        "child_centering": geometry.CENTERING_LETTER[row["child"]],
        "basis_machine": [[int(value) for value in vector] for vector in row["basis"]],
        "origin_machine": [str(value) for value in row["origin"]],
    }


def blank_evidence():
    """The per-record evidence keys that only a matched oracle row can fill."""
    return {"u": None, "u_det": None, "u_integral": None, "origin_exact": None,
            "origin_delta": None, "origin_converted": None, "origin_printed": None,
            "basis_converted": None, "basis_printed": None, "oracle_row": None}


def derive_record(row, oracle_row, parent):
    """Derive U and the origin comparison for one matched record."""
    child = row["child"]
    p_parent = geometry.PRIMITIVE_BASIS[geometry.CENTERING_LETTER[parent]]
    p_child = geometry.PRIMITIVE_BASIS[geometry.CENTERING_LETTER[child]]
    converted = operations.matmul3(
        [[Fraction(value) for value in vector] for vector in row["basis"]],
        [[Fraction(value) for value in vector] for vector in p_parent],
    )
    printed_primitive = operations.matmul3(
        [[Fraction(value) for value in vector] for vector in p_child], oracle_row["basis"]
    )
    evidence = blank_evidence()
    try:
        u = operations.matmul3(converted, operations.inverse3(printed_primitive))
    except RuntimeError:
        evidence.update(status="singular_basis", detail="printed basis is singular")
        return evidence
    det_u = operations.determinant3([value for vector in u for value in vector])
    u_integral = all(value.denominator == 1 for vector in u for value in vector)
    origin_converted = geometry.to_conventional(row["origin"], p_parent)
    delta = [oracle_row["origin"][i] - origin_converted[i] for i in range(3)]
    origin_exact = all(value == 0 for value in delta)
    if not u_integral:
        status, detail = "nonunimodular", "u_nonintegral"
    elif abs(det_u) != 1:
        status, detail = "nonunimodular", "u_determinant"
    elif not origin_exact:
        status, detail = "origin_mismatch", "u_integral_unimodular_origin_delta"
    else:
        status, detail = "candidate", "origin_exact_u_unimodular"
    evidence.update(
        status=status,
        detail=detail,
        u=[[str(value) for value in vector] for vector in u],
        u_det=str(det_u),
        u_integral=u_integral,
        origin_exact=origin_exact,
        origin_delta=[str(value) for value in delta],
        origin_converted=[str(value) for value in origin_converted],
        origin_printed=[str(value) for value in oracle_row["origin"]],
        basis_converted=[[str(value) for value in vector] for vector in converted],
        basis_printed=[[str(value) for value in vector] for vector in oracle_row["basis"]],
        oracle_row=oracle_row["raw"],
    )
    return evidence


def classify_query(machine_rows, oracle_rows, parent):
    """Return (query_status, per-record evidence, query error or None)."""
    labels = [row["direction_label"] for row in machine_rows]
    blank = blank_evidence()

    def failure(kind, detail):
        return "error", [
            dict(blank, status=kind, detail=detail) for _ in machine_rows
        ], {"kind": kind, "detail": detail}

    if len(set(labels)) != len(labels):
        return failure("duplicate_machine_label", "machine labels are not unique")
    if len({row["label"] for row in oracle_rows}) != len(oracle_rows):
        return failure("duplicate_oracle_label", "oracle labels are not unique")
    if not oracle_rows:
        return "empty", [
            dict(blank, status=UNRETURNED, detail="official table is empty")
            for _ in machine_rows
        ], None
    if len(oracle_rows) != len(machine_rows):
        return failure(
            "row_count_mismatch",
            f"oracle {len(oracle_rows)} rows, machine {len(machine_rows)} rows",
        )
    by_label = {row["label"]: row for row in oracle_rows}
    extra = sorted(set(by_label) - set(labels))
    if extra:
        return failure("unexpected_oracle_row", ",".join(extra))
    results = []
    for row in machine_rows:
        oracle_row = by_label.get(row["direction_label"])
        if oracle_row is None:
            return failure("missing_oracle_row", row["direction_label"])
        if oracle_row["child"] != row["child"]:
            return failure(
                "child_mismatch",
                f"{row['direction_label']}: oracle child {oracle_row['child']} "
                f"!= machine {row['child']}",
            )
        if oracle_row["size"] != row["size"]:
            return failure(
                "size_mismatch",
                f"{row['direction_label']}: oracle size {oracle_row['size']} "
                f"!= machine {row['size']}",
            )
        results.append(derive_record(row, oracle_row, parent))
    return "nonempty", results, None


def collect(queries, records_by_query, binary, data_dir, work_root, workers):
    records = []
    outcomes = {}
    errors = []
    completed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(run_query, binary, data_dir, work_root, query): query
            for query in queries
        }
        for future in concurrent.futures.as_completed(futures):
            query = futures[future]
            outcome = future.result()
            machine_rows = records_by_query[(query["parent"], query["index"])]
            if outcome["error"] is not None:
                status, query_error = "error", dict(
                    outcome["error"], stdout_tail=outcome["stdout_tail"]
                )
                results = [
                    dict(blank_evidence(), status="oracle_error",
                         detail=outcome["error"]["detail"])
                    for _ in machine_rows
                ]
            else:
                status, results, query_error = classify_query(
                    machine_rows, outcome["rows"], query["parent"]
                )
            outcomes[(query["parent"], query["index"])] = status
            for row, evidence in zip(machine_rows, results):
                records.append(dict(base_record(row, query["parent"]), **evidence))
            if query_error is not None:
                errors.append(dict(query_error, parent=query["parent"],
                                   source_irrep_index=query["index"], irrep=query["irrep"]))
            completed += 1
            if completed % PROGRESS_EVERY == 0:
                print(f"collected {completed}/{len(queries)} queries", file=sys.stderr, flush=True)
    records.sort(key=lambda record: record["ordinal"])
    return records, outcomes, errors


def build_summary(queries, records, outcomes, errors, expected_records):
    counts = collections.Counter((record["status"], record["detail"]) for record in records)
    by_status = collections.Counter()
    for (status, _), count in counts.items():
        by_status[status] += count
    classified = sum(by_status[status] for status in CLASSIFIED)
    unreturned = by_status[UNRETURNED]
    failed = len(records) - classified - unreturned
    query_counts = collections.Counter(outcomes.values())
    irrep_of = {(query["parent"], query["index"]): query["irrep"] for query in queries}
    checks = {
        "queries_partition": (
            query_counts["nonempty"] + query_counts["empty"] + query_counts["error"]
            == len(queries)
        ),
        "records_partition": classified + unreturned + failed == len(records),
        "classification_partition": (
            by_status["candidate"] + by_status["nonunimodular"]
            + by_status["origin_mismatch"] == classified
        ),
        "every_machine_record_emitted": len(records) == expected_records,
    }
    if not all(checks.values()):
        errors = list(errors) + [{"kind": "internal_accounting_error", "detail": str(checks)}]
    return {
        "type": "summary",
        "queries_total": len(queries),
        "queries_nonempty": query_counts["nonempty"],
        "queries_empty": query_counts["empty"],
        "queries_error": query_counts["error"],
        "empty_queries": [
            {"parent": parent, "source_irrep_index": index,
             "irrep": irrep_of.get((parent, index), "")}
            for (parent, index), status in sorted(outcomes.items()) if status == "empty"
        ],
        "records_total": len(records),
        "records_classified": classified,
        "records_unreturned_empty": unreturned,
        "records_failed": failed,
        "classification": {
            "candidate": by_status["candidate"],
            "nonunimodular": by_status["nonunimodular"],
            "origin_mismatch": by_status["origin_mismatch"],
        },
        "origin_exact": sum(1 for record in records if record["origin_exact"] is True),
        "u_integral": sum(1 for record in records if record["u_integral"] is True),
        "u_unimodular": sum(
            1 for record in records
            if record["u_integral"] is True and record["u_det"] in ("1", "-1")
        ),
        "detail_counts": {
            f"{status}|{detail}": count for (status, detail), count in sorted(counts.items())
        },
        "error_total": len(errors),
        "errors": errors[:MAX_ERRORS_IN_SUMMARY],
        "errors_truncated": len(errors) > MAX_ERRORS_IN_SUMMARY,
        "collector_completed": not errors and all(checks.values()),
        "engine_coverage_claim": "none; collector success means the census completed",
        "accounting": checks,
    }


def write_output(path, manifest, records, summary):
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    partial = path + ".partial"
    if path.endswith(".jsonl"):
        with open(partial, "w", encoding="utf-8") as handle:
            handle.write(json.dumps(manifest, sort_keys=True) + "\n")
            for record in records:
                handle.write(json.dumps(record, sort_keys=True) + "\n")
            handle.write(json.dumps(summary, sort_keys=True) + "\n")
    elif path.endswith(".json"):
        with open(partial, "w", encoding="utf-8") as handle:
            json.dump(
                {"manifest": manifest, "summary": summary, "records": records},
                handle, sort_keys=True, indent=1,
            )
    else:
        raise CollectorError("--output must end in .jsonl or .json")
    os.replace(partial, path)


def parse_parents(raw):
    if raw is None:
        return list(range(1, 231))
    values = []
    for item in raw:
        for token in str(item).split(","):
            token = token.strip()
            if not token.isdigit() or not 1 <= int(token) <= 230:
                raise CollectorError(f"invalid --parent {token!r}; expected a space group 1..230")
            values.append(int(token))
    return sorted(set(values))


def main(argv=None):
    parser = argparse.ArgumentParser(description="Task 9 oracle setting census")
    parser.add_argument("--output", required=True, help="census path ending in .jsonl or .json")
    parser.add_argument("--parent", action="append", metavar="SG",
                        help="restrict to a parent space group (repeatable/comma separated)")
    parser.add_argument("--workers", type=int, default=4, help="parallel iso processes (default 4)")
    args = parser.parse_args(argv)
    try:
        if args.workers < 1:
            raise CollectorError("--workers must be >= 1")
        if not (args.output.endswith(".jsonl") or args.output.endswith(".json")):
            raise CollectorError("--output must end in .jsonl or .json")
        parents = parse_parents(args.parent)
        archive = archive_manifest()
        inventory = load_inventory()
        all_records = load_machine_records(inventory)
        in_scope = [entry for entry in inventory if entry["parent"] in parents]
        if not in_scope:
            raise CollectorError("no irreps selected")
        records_by_query = collections.defaultdict(list)
        for record in all_records:
            if record["parent"] in parents:
                records_by_query[(record["parent"], record["source_irrep_index"])].append(record)
        missing = [(q["parent"], q["index"]) for q in in_scope
                   if (q["parent"], q["index"]) not in records_by_query]
        if missing:
            raise CollectorError(f"queries without machine records: {missing[:5]}")
        full = len(parents) == 230
        if full and (len(in_scope) != EXPECTED_IRREPS or len(all_records) != EXPECTED_RECORDS):
            raise CollectorError("pinned inventory counts changed")
        with tempfile.TemporaryDirectory(prefix="task9-oracle-") as root:
            data_dir, members = extract_oracle(root)
            work_root = os.path.join(root, "work")
            os.makedirs(work_root)
            records, outcomes, errors = collect(
                in_scope, records_by_query, os.path.join(data_dir, "iso"),
                data_dir, work_root, args.workers,
            )
        summary = build_summary(in_scope, records, outcomes, errors, len(all_records)
                                if full else sum(len(v) for v in records_by_query.values()))
        manifest = {
            "type": "manifest",
            "schema": "task9-oracle-census/1",
            "ordinal_base": 0,
            "tool": "scripts/audit_subduction_settings.py",
            "iso_version": PINNED_ISO_VERSION,
            "setting_command": SETTING_COMMAND,
            "pinned_commands": list(ORACLE_COMMANDS),
            "archive": archive,
            "members": members,
            "workers": args.workers,
            "parents": parents,
            "scope": "full" if full else "subset",
            "timeout_seconds": ORACLE_TIMEOUT_SECONDS,
            "semantics": {
                "collector_success": "every in-scope machine record was collected and classified",
                "engine_coverage": "not claimed by this collector",
                "candidate": "exact origin and integral unimodular U; no delta resolution, "
                             "no canonical operation matching, no complete decomposition",
            },
        }
        write_output(args.output, manifest, records, summary)
    except CollectorError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    except OSError as error:
        print(f"I/O error: {error}", file=sys.stderr)
        return 2
    print(
        f"scope={manifest['scope']} parents={len(parents)} "
        f"queries={summary['queries_total']} "
        f"(nonempty={summary['queries_nonempty']} empty={summary['queries_empty']} "
        f"error={summary['queries_error']}) records={summary['records_total']} "
        f"candidate={summary['classification']['candidate']} "
        f"nonunimodular={summary['classification']['nonunimodular']} "
        f"origin_mismatch={summary['classification']['origin_mismatch']} "
        f"unreturned={summary['records_unreturned_empty']} "
        f"failed={summary['records_failed']} "
        f"collector_completed={summary['collector_completed']}"
    )
    print(f"wrote {args.output} (engine coverage not claimed)")
    return 0 if summary["collector_completed"] else 1


if __name__ == "__main__":
    sys.exit(main())
