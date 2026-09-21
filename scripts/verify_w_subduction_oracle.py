#!/usr/bin/env python3
"""Cross-check the pinned other-wave-vector rows against the live ISOTROPY program.

`data_isotropy.txt` carries two subduction tables for every isotropy record:

* `isotropy_subduce_*` -- the irreps whose subduction onto the record's subgroup
  contains the subgroup identity (anchored per record, and covered by the task 9
  audit's engine comparison),
* `isotropy_w_subduce_*` -- 5,756 additional rows naming irreps **at other wave
  vectors**: the `DT`/`SM` star labels of the eleven cubic parents 196, 202, 203,
  209, 210, 216, 219, 225, 226, 227, 228.  `data_irreps.txt` stores no k vector
  and no character row for those sources, so the audit reported them separately
  instead of computing them.

The official program computes those very numbers.  For a chosen
`(parent, irrep)` the command

    SHOW SUBGROUP
    SHOW SIZE
    SHOW DIRECTION
    SHOW FREQUENCY
    DISPLAY ISOTROPY

prints one row per order-parameter direction,

    <subgroup#> <symbol> <Size> <Dir> <frequency list>

where the frequency list is "which irreps (at any k vector) subduce this
subgroup, with which multiplicity".  The entries whose labels live in the
other-wave-vector table are exactly the pinned `isotropy_w_subduce_*` rows of
that record, so this script turns the 5,756 rows into a live oracle comparison.

The comparison is keyed on `(subgroup number, direction label)`: the pinned
`isotropy_orderparam_label` is the program's own `Dir` column, so no frame
conversion is needed.  Both directions are checked for every `(parent, irrep)`
that carries such rows:

* every oracle entry with an other-wave-vector label must have a pinned row with
  the same frequency (a missing pinned row is a failure),
* every pinned other-wave-vector row must appear in the oracle row.

Usage::

    python3 scripts/verify_w_subduction_oracle.py [--json report.json]
    python3 scripts/verify_w_subduction_oracle.py --limit 5 --verbose

Exit codes: 0 = every compared row agreed, 1 = a row disagreed, 2 = the pinned
archive or the extracted binary is missing.
"""

import argparse
import json
import os
import re
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(SCRIPT_DIR)
ISO_DIR = os.path.join(REPO_DIR, "isotropy_subgroup")

sys.path.insert(0, SCRIPT_DIR)

import generate_irrep_data as generator  # noqa: E402


class PinnedArchive:
    """The pinned tables this gate compares, in record order."""

    def __init__(self):
        irrep_lines = generator.read_file("data_irreps.txt")
        isotropy_lines = generator.read_file("data_isotropy.txt")
        little_lines = generator.read_file("data_little.txt")
        irrep_sec = generator.get_sections(irrep_lines)
        isotropy_sec = generator.get_sections(isotropy_lines)
        little_sec = generator.get_sections(little_lines)

        self.parent = generator.parse_ints(isotropy_lines, isotropy_sec, "isotropy_parent")
        self.subgroup = generator.parse_ints(
            isotropy_lines, isotropy_sec, "isotropy_subgroup"
        )
        self.irrep = generator.parse_ints(isotropy_lines, isotropy_sec, "isotropy_irrep")
        self.direction_label = generator.parse_labels(
            isotropy_lines, isotropy_sec, "isotropy_orderparam_label"
        )
        # Other-wave-vector rows: packed in record order, counted per record.
        self.w_count = generator.parse_ints(
            isotropy_lines, isotropy_sec, "isotropy_w_subduce_count"
        )
        w_irrep = generator.parse_ints(
            isotropy_lines, isotropy_sec, "isotropy_w_subduce_irrep"
        )
        w_frequency = generator.parse_ints(
            isotropy_lines, isotropy_sec, "isotropy_w_subduce_frequency"
        )
        self.w_label = generator.parse_labels(irrep_lines, irrep_sec, "irrep_w_label")

        self.irrep_label = generator.parse_labels(irrep_lines, irrep_sec, "irrep_label")
        # The program only accepts the compact little-irrep spelling; 112 of the
        # stored labels use the legacy Miller-Love spelling (`W1W1` vs `W1WA1`).
        compact = generator.parse_labels(little_lines, little_sec, "little_irr_full_label")
        old_map = generator.parse_ints(little_lines, little_sec, "little_irr_old_map")
        self.compact_label = [
            compact[old_map[index] - 1].strip() for index in range(len(self.irrep_label))
        ]

        self.other_wave_vector = {}
        packed = 0
        for record, count in enumerate(self.w_count):
            rows = []
            for _ in range(count):
                rows.append(
                    (self.w_label[w_irrep[packed] - 1].strip(), w_frequency[packed])
                )
                packed += 1
            if rows:
                self.other_wave_vector[record + 1] = rows

    def records_with_other_wave_vectors(self):
        return sorted(self.other_wave_vector)

    def context(self, ordinal):
        index = ordinal - 1
        return {
            "parent": self.parent[index],
            "subgroup": self.subgroup[index],
            "irrep_label": self.irrep_label[self.irrep[index] - 1].strip(),
            "compact_label": self.compact_label[self.irrep[index] - 1],
            "direction_label": self.direction_label[index].strip(),
        }

    def records_of_irrep(self, sg, compact_label):
        """Every pinned record of one `(parent, compact irrep label)`.

        Returns `(subgroup, direction label) -> (ordinal, w rows)`.  Direction
        labels are unique per record within one irrep's table, so a collision
        would mean the pinned table itself is ambiguous.
        """
        found = {}
        for index in range(len(self.parent)):
            context = self.context(index + 1)
            if context["parent"] != sg or context["compact_label"] != compact_label:
                continue
            key = (context["subgroup"], context["direction_label"])
            if key in found:
                raise RuntimeError(
                    f"SG {sg} {compact_label}: two records share subgroup/direction {key}"
                )
            found[key] = (index + 1, list(self.other_wave_vector.get(index + 1, [])))
        return found


ROW_RE = re.compile(r"^(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s*(.*)$")
FREQUENCY_RE = re.compile(r"^(\d+)\s+(\S+)$")


def parse_frequency_list(text):
    """`"1 GM1, 2 GM2GM3"` -> `[("GM1", 1), ("GM2GM3", 2)]`."""
    rows = []
    for chunk in text.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        match = FREQUENCY_RE.match(chunk)
        if not match:
            raise ValueError(f"unparsable frequency entry {chunk!r}")
        rows.append((match.group(2), int(match.group(1))))
    return rows


def parse_isotropy_rows(stdout):
    """The `DISPLAY ISOTROPY` table keyed by `(subgroup, direction label)`.

    A subgroup/direction pair can appear more than once (one row per domain
    origin); each key therefore keeps a list of frequency lists.  A long
    frequency list is wrapped by the program onto indented continuation lines
    (`SC 1000` avoids the wrap, and the parser joins them anyway so the gate
    cannot silently drop the tail of a row).
    """
    rows = {}
    pending = None
    for line in stdout.splitlines():
        stripped = line.rstrip("*").strip()
        if not stripped or stripped.startswith("Subgroup"):
            continue
        # Table rows start in column 0; an indented line continues the last
        # row's frequency list.  A continuation can itself look like a row
        # (`1 DT2, 1 SM1`), so the indentation decides, not the pattern.
        if line[:1].isspace() and pending is not None:
            pending[1] = f"{pending[1]} {stripped}"
            continue
        match = ROW_RE.match(stripped)
        if not match:
            continue
        key = (int(match.group(1)), match.group(4))
        pending = [key, match.group(5)]
        rows.setdefault(key, []).append(pending)
    for key, entries in rows.items():
        rows[key] = [parse_frequency_list(text) for _, text in entries]
    return rows


def run_oracle(sg, compact_label, timeout=300):
    """Frequency lists of one `(parent, irrep)` from the official program."""
    binary = os.path.join(ISO_DIR, "iso")
    if not os.path.isfile(binary):
        raise FileNotFoundError(
            "the ISOTROPY binary is not extracted; run\n"
            f"  unzip -o {os.path.join(ISO_DIR, 'iso.zip')} -d {ISO_DIR}\n"
            "before using this oracle"
        )
    commands = [
        "PAGE 1000",
        # The frequency list of a low-symmetry row can exceed the default screen
        # width; `SC 1000` keeps every row on one line.
        "SC 1000",
        "SET I ALL OR 1",
        f"VALUE PARENT {sg}",
        f"VALUE IRREP {compact_label}",
        "SHOW SUBGROUP",
        "SHOW SIZE",
        "SHOW DIRECTION",
        "SHOW FREQUENCY",
        "DISPLAY ISOTROPY",
        "QUIT",
    ]
    result = subprocess.run(
        [binary],
        input="\n".join(commands) + "\n",
        capture_output=True,
        text=True,
        cwd=ISO_DIR,
        env=dict(os.environ, ISODATA=ISO_DIR + os.sep),
        timeout=timeout,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"iso exited with {result.returncode} for SG {sg} {compact_label}: "
            f"{result.stderr.strip()[:200]}"
        )
    rows = parse_isotropy_rows(result.stdout)
    if not rows:
        raise RuntimeError(f"iso printed no isotropy rows for SG {sg} {compact_label}")
    return rows


def compare_group(archive, sg, compact_label, oracle_rows, w_labels):
    """Compare the other-wave-vector rows of one `(parent, irrep)` group."""
    pinned = archive.records_of_irrep(sg, compact_label)
    failures = []
    compared = {"oracle_rows": 0, "oracle_w_rows": 0, "pinned_w_rows": 0, "records": 0}
    for key, candidates in sorted(oracle_rows.items()):
        pinned_entry = pinned.get(key)
        oracle_w = []
        for candidate in candidates:
            oracle_w.append(
                sorted((label, count) for label, count in candidate if label in w_labels)
            )
        compared["oracle_rows"] += 1
        compared["oracle_w_rows"] += sum(len(rows) for rows in oracle_w)
        if pinned_entry is None:
            if any(oracle_w):
                failures.append(
                    f"SG {sg} {compact_label} subgroup {key[0]} direction {key[1]}: "
                    f"oracle lists {oracle_w} but the pinned table has no such record"
                )
            continue
        ordinal, pinned_rows = pinned_entry
        compared["records"] += 1
        compared["pinned_w_rows"] += len(pinned_rows)
        expected = sorted(pinned_rows)
        if any(rows == expected for rows in oracle_w):
            continue
        failures.append(
            f"SG {sg} {compact_label} subgroup {key[0]} direction {key[1]} "
            f"(ordinal {ordinal}): pinned "
            + (", ".join(f"{count} {label}" for label, count in expected) or "(none)")
            + " | oracle "
            + " / ".join(
                ", ".join(f"{count} {label}" for label, count in rows) or "(none)"
                for rows in oracle_w
            )
        )
    # The other direction: a pinned row whose oracle row is missing or filtered out.
    for key, (ordinal, pinned_rows) in sorted(pinned.items()):
        if not pinned_rows:
            continue
        if key not in oracle_rows:
            failures.append(
                f"SG {sg} {compact_label} subgroup {key[0]} direction {key[1]} "
                f"(ordinal {ordinal}): pinned {len(pinned_rows)} other-wave-vector "
                "rows have no oracle row"
            )
    return compared, failures


def check(archive, limit=None, verbose=False):
    """Compare every `(parent, irrep)` group that carries other-wave-vector rows."""
    w_labels = {label.strip() for label in archive.w_label}
    ordinals = archive.records_with_other_wave_vectors()
    if limit is not None:
        ordinals = ordinals[:limit]
    groups = {}
    for ordinal in ordinals:
        context = archive.context(ordinal)
        groups.setdefault((context["parent"], context["compact_label"]), None)

    totals = {"groups": 0, "records": 0, "oracle_rows": 0, "oracle_w_rows": 0, "pinned_w_rows": 0}
    labels_seen = {
        label for rows in archive.other_wave_vector.values() for label, _ in rows
    }
    failures = []
    for sg, compact_label in sorted(groups):
        oracle_rows = run_oracle(sg, compact_label)
        compared, group_failures = compare_group(
            archive, sg, compact_label, oracle_rows, w_labels
        )
        totals["groups"] += 1
        for name, value in compared.items():
            totals[name] += value
        failures.extend(group_failures)
        if verbose:
            print(
                f"  SG {sg} {compact_label}: {compared['records']} records, "
                f"{compared['oracle_w_rows']} oracle w rows, "
                f"{len(group_failures)} failures so far",
                flush=True,
            )
    return totals, sorted(labels_seen), failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=None, help="check only the first N records")
    parser.add_argument("--json", dest="json_path", default=None)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    try:
        archive = PinnedArchive()
    except FileNotFoundError as error:
        print(f"pinned input missing: {error}", file=sys.stderr)
        return 2
    if not os.path.isfile(os.path.join(ISO_DIR, "iso")):
        print(
            "the ISOTROPY binary is not extracted; run\n"
            f"  unzip -o {os.path.join(ISO_DIR, 'iso.zip')} -d {ISO_DIR}",
            file=sys.stderr,
        )
        return 2

    totals, labels_seen, failures = check(archive, limit=args.limit, verbose=args.verbose)
    print(
        "groups checked: {groups} | records checked: {records} | oracle rows: "
        "{oracle_rows} | oracle other-wave-vector rows: {oracle_w_rows} | pinned "
        "other-wave-vector rows: {pinned_w_rows}".format(**totals)
    )
    print(f"w labels seen: {len(labels_seen)}")
    for message in failures[:20]:
        print(f"FAIL {message}")
    if len(failures) > 20:
        print(f"... and {len(failures) - 20} more failures")
    print(f"mismatches: {len(failures)}")
    if args.json_path:
        with open(args.json_path, "w", encoding="utf-8") as handle:
            json.dump(
                {
                    "totals": totals,
                    "w_labels_seen": labels_seen,
                    "mismatches": failures,
                },
                handle,
                indent=2,
                sort_keys=True,
            )
            handle.write("\n")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
